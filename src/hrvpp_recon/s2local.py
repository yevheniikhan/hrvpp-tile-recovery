"""Reconstruction without Earth Engine: Sentinel-2 L2A from Microsoft Planetary Computer,
computed in Colab (NumPy on CPU or CuPy on GPU).

Same algorithm as gee_recon.py / the Code Editor script (PPI after the HR-VPP ATBD, 10-day
median composites, winter base level, linear gap filling <= 60 days, two 1-2-1 passes,
season metrics of phenology_np.extract_metrics). Known differences to the Earth Engine run:

* cloud mask: SCL only, classes 0,1,3,7,8,9,10,11 excluded (Cloud Score+ is not available
  outside Earth Engine; adding 7 and 11 brought a 2021 pilot closest to the Earth Engine run);
* no L1C fallback for scenes without L2A (relevant for 2017-2018 only);
* the tile-constant MDVI term (mean of DVI_max - DVI_q95) is computed from 160 m overviews;
* percentiles/medians use NumPy definitions (Earth Engine may differ in ties/interpolation).

Output: one Int16 GeoTIFF per year with the bands and factors of the Earth Engine export
(AMPL, LENGTH, RSLOPE, TPROD, SOSD, EOSD, MAXV, MINV, NOBS; nodata -32768), so the same
correction and MONIFUN conversion code can be applied.
"""
from __future__ import annotations

import calendar
import concurrent.futures as cf
import datetime as dt
import os
import time
from dataclasses import dataclass
from typing import Callable, List, Optional

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.warp import transform_bounds

from . import phenology_np as pnp
from .config import GEE_EXPORT_FACTORS, ReconParams

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
BAD_SCL = (0, 1, 3, 7, 8, 9, 10, 11)   # +7 unclassified, +11 snow: closest to the Cloud Score+ run (pilot 2021)
BANDS_OUT = ["AMPL", "LENGTH", "RSLOPE", "TPROD", "SOSD", "EOSD", "MAXV", "MINV", "NOBS"]
NODATA = -32768

os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
os.environ.setdefault("GDAL_HTTP_MERGE_CONSECUTIVE_RANGES", "YES")
os.environ.setdefault("GDAL_HTTP_MULTIPLEX", "YES")
os.environ.setdefault("GDAL_HTTP_VERSION", "2")
os.environ.setdefault("VSI_CACHE", "TRUE")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "5")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "2")


def xp_module(use_gpu: bool = True):
    if use_gpu:
        try:
            import cupy
            cupy.zeros(1)
            return cupy
        except Exception:
            pass
    return np


# --------------------------------------------------------------------------- scenes
@dataclass
class Scene:
    when: dt.datetime
    sza: float
    offset: float          # 1000 for processing baseline >= 04.00 (harmonisation), else 0
    b04: str
    b08: str
    scl: str
    bbox: tuple            # lon/lat


def search_scenes(bbox_ll, start: str, end: str, max_cloud: float = 95.0) -> List[Scene]:
    import planetary_computer as pc
    from pystac_client import Client
    cat = Client.open(STAC, modifier=pc.sign_inplace)
    items = cat.search(collections=["sentinel-2-l2a"], bbox=list(bbox_ll), datetime=f"{start}/{end}",
                       query={"eo:cloud_cover": {"lte": max_cloud}}).item_collection()
    out = []
    for it in items:
        p = it.properties
        sza = p.get("s2:mean_solar_zenith")
        if sza is None:
            continue
        try:
            base = float(str(p.get("s2:processing_baseline", "0")))
        except ValueError:
            base = 0.0
        when = dt.datetime.fromisoformat(p["datetime"].replace("Z", "+00:00")).replace(tzinfo=None)
        out.append(Scene(when, float(sza), 1000.0 if base >= 4.0 else 0.0, it.assets["B04"].href,
                         it.assets["B08"].href, it.assets["SCL"].href, tuple(it.bbox)))
    out.sort(key=lambda s: s.when)
    return out


def _resign(href: str) -> str:
    import planetary_computer as pc
    return pc.sign(href.split("?")[0])


def _read(href, crs, transform, width, height, overview_level=None, resampling=Resampling.nearest):
    kw = {"overview_level": overview_level} if overview_level is not None else {}
    for attempt in range(4):
        try:
            if attempt == 2:
                kw = {}          # no such overview level: read full resolution and resample
            with rasterio.open(_resign(href), **kw) as src:
                with WarpedVRT(src, crs=crs, transform=transform, width=width, height=height,
                               resampling=resampling, src_nodata=0, nodata=0) as vrt:
                    return vrt.read(1)
        except Exception:
            if attempt == 3:
                raise
            time.sleep(3 * (attempt + 1))


def read_dvi(scene: Scene, crs, transform, width, height, coarse=False):
    """Cloud-masked DVI on the target grid (NaN where masked / no data)."""
    lv = (3, 3, 2) if coarse else (None, None, None)
    rs = Resampling.average if coarse else Resampling.nearest
    b4 = _read(scene.b04, crs, transform, width, height, lv[0], rs).astype("float32")
    b8 = _read(scene.b08, crs, transform, width, height, lv[1], rs).astype("float32")
    scl = _read(scene.scl, crs, transform, width, height, lv[2], Resampling.nearest)
    ok = (b4 > 0) & (b8 > 0) & ~np.isin(scl, BAD_SCL)
    dvi = ((b8 - scene.offset) - (b4 - scene.offset)) / 10000.0
    return np.where(ok, dvi, np.nan).astype("float32" if coarse else "float16")


# --------------------------------------------------------------------------- algorithm
def nanquantile(xp, a, q, axis=0):
    """Quantile ignoring NaN (linear interpolation, NumPy definition); works for CuPy."""
    s = xp.sort(a, axis=axis)                     # NaN sorted to the end
    n = (~xp.isnan(a)).sum(axis=axis)
    pos = (n - 1).clip(0) * q
    lo = xp.floor(pos).astype("int64")
    hi = xp.minimum(lo + 1, (n - 1).clip(0))
    take = lambda i: xp.take_along_axis(s, xp.expand_dims(i, axis), axis=axis).squeeze(axis)
    v = take(lo) + (take(hi) - take(lo)) * (pos - lo)
    return xp.where(n > 0, v, xp.nan)


def ppi_from_dvi(xp, dvi, sza, mdvi, p: ReconParams):
    """dvi: (n_obs, P); sza: (n_obs,); mdvi: (P,)."""
    cosz = xp.cos(xp.asarray(sza, dtype="float32") * np.float32(np.pi / 180))[:, None]
    dc = 0.0336 + 0.0477 / cosz
    qe = dc + (1 - dc) * p.leaf_g / cosz
    k = (1.0 / (4 * qe)) * (1 + mdvi[None, :]) / (1 - mdvi[None, :])
    d = xp.minimum(dvi, mdvi[None, :] - 1e-4)
    d = xp.maximum(d, p.dvi_soil)
    ppi = -k * xp.log((mdvi[None, :] - d) / (mdvi[None, :] - p.dvi_soil))
    return xp.where(xp.isnan(dvi), xp.nan, xp.clip(ppi, 0, 3))


def metrics_from_ppi(xp, ppi, obs_bin, year, p: ReconParams):
    """ppi: (n_obs, P) with NaN; obs_bin: (n_obs,) bin index or -1. Returns dict of (P,) arrays."""
    centres, ws, nbins = pnp.window_bins(year, p.bin_days)
    P = ppi.shape[1]
    jan1 = dt.date(year, 1, 1)
    base = nanquantile(xp, ppi, p.winter_base_pct / 100.0, axis=0)
    series = xp.full((nbins, P), xp.nan, dtype="float32")
    nobs = xp.zeros(P, dtype="float32")
    obs_bin = np.asarray(obs_bin)
    for b in range(nbins):
        rows = np.nonzero(obs_bin == b)[0]
        if len(rows):
            sub = ppi[rows]
            cnt = (~xp.isnan(sub)).sum(axis=0)
            nobs += cnt
            med = nanquantile(xp, sub, 0.5, axis=0)
            series[b] = xp.where(cnt > 0, med, xp.nan)
        month = (jan1 + dt.timedelta(days=float(centres[b]))).month
        if month in p.winter_months:
            series[b] = xp.where(xp.isnan(series[b]), base, series[b])
    # linear interpolation between the nearest valid bins within max_gap_days
    valid = ~xp.isnan(series)
    t = xp.asarray(centres, dtype="float32")
    idx = xp.arange(nbins)[:, None] * xp.ones((1, P), dtype="int64")
    prev = xp.where(valid, idx, -1)
    nxt = xp.where(valid, idx, nbins)
    for b in range(1, nbins):
        prev[b] = xp.maximum(prev[b], prev[b - 1])
    for b in range(nbins - 2, -1, -1):
        nxt[b] = xp.minimum(nxt[b], nxt[b + 1])
    pc_ = prev.clip(0, nbins - 1)
    nc_ = nxt.clip(0, nbins - 1)
    tb, ta = t[pc_], t[nc_]
    vb = xp.take_along_axis(series, pc_, axis=0)
    va = xp.take_along_axis(series, nc_, axis=0)
    tt = t[:, None]
    ok = (prev >= 0) & (nxt < nbins) & ((tt - tb) <= p.max_gap_days) & ((ta - tt) <= p.max_gap_days)
    span = xp.where(ta - tb == 0, 1, ta - tb)
    interp = vb + (va - vb) * (tt - tb) / span
    series = xp.where(valid, series, xp.where(ok, interp, xp.nan))
    complete = (~xp.isnan(series)).all(axis=0)
    V = xp.where(xp.isnan(series), 0, series).T                 # (P, nbins)
    for _ in range(p.smooth_passes):
        left = xp.concatenate([V[:, :1], V[:, :-1]], axis=1)
        right = xp.concatenate([V[:, 1:], V[:, -1:]], axis=1)
        V = (left + 2 * V + right) / 4.0
    m = extract_metrics_xp(xp, V, xp.asarray(centres), year, p)
    keep = complete & (nobs >= p.min_obs) & ~xp.isnan(m["AMPL"])
    out = {k: xp.where(keep, v, xp.nan) for k, v in m.items()}
    out["NOBS"] = nobs
    out["_keep"] = keep
    return out


def extract_metrics_xp(xp, V, t, year, p: ReconParams):
    """phenology_np.extract_metrics with an array module argument (NumPy or CuPy)."""
    P, N = V.shape
    idx = xp.arange(N)
    len_y = 366 if calendar.isleap(year) else 365
    half = p.bin_days / 2.0
    in_year = (t >= 0) & (t < len_y)
    Vm = xp.where(in_year[None, :], V, -1.0)
    maxv = Vm.max(axis=1)
    peak = xp.where(Vm == maxv[:, None], idx[None, :], 999).min(axis=1)
    left = idx[None, :] <= peak[:, None]
    right = idx[None, :] >= peak[:, None]
    lmin = xp.where(left, V, 99.0).min(axis=1)
    rmin = xp.where(right, V, 99.0).min(axis=1)
    lmin_i = xp.where(left & (V == lmin[:, None]), idx[None, :], -1).max(axis=1)
    rmin_i = xp.where(right & (V == rmin[:, None]), idx[None, :], 999).min(axis=1)
    minv = (lmin + rmin) / 2.0
    ampl = maxv - minv
    thr_s = lmin + p.sos_frac * (maxv - lmin)
    thr_e = rmin + p.eos_frac * (maxv - rmin)
    rng_s = (idx[None, :] >= lmin_i[:, None]) & left
    rng_e = (idx[None, :] <= rmin_i[:, None]) & right
    sos_i = xp.where(rng_s & (V >= thr_s[:, None]), idx[None, :], 999).min(axis=1)
    eos_i = xp.where(rng_e & (V >= thr_e[:, None]), idx[None, :], -1).max(axis=1)
    rows = xp.arange(P)
    s0 = xp.clip(xp.maximum(sos_i - 1, lmin_i), 0, N - 1)
    s1 = xp.clip(sos_i, 0, N - 1)
    e0 = xp.clip(eos_i, 0, N - 1)
    e1 = xp.clip(xp.minimum(eos_i + 1, rmin_i), 0, N - 1)
    vs0, vs1, ve0, ve1 = V[rows, s0], V[rows, s1], V[rows, e0], V[rows, e1]
    ds = vs1 - vs0
    fs = xp.clip(xp.where(ds != 0, (thr_s - vs0) / xp.where(ds == 0, 1, ds), 1.0), 0, 1)
    de = ve0 - ve1
    fe = xp.clip(xp.where(de != 0, (ve0 - thr_e) / xp.where(de == 0, 1, de), 0.0), 0, 1)
    t0 = float(t[0])
    sos_t = t0 + p.bin_days * (s0 + fs * (s1 - s0))
    eos_t = t0 + p.bin_days * (e0 + fe * (e1 - e0))
    length = eos_t - sos_t
    lo = xp.maximum(t[None, :] - half, sos_t[:, None])
    hi = xp.minimum(t[None, :] + half, eos_t[:, None])
    w = xp.clip(hi - lo, 0, None)
    tprod = (V * w).sum(axis=1)
    rslope = xp.abs(ve0 - ve1) / p.bin_days
    valid = (ampl >= p.min_ampl) & (length >= p.min_length) & (sos_i < 999) & (eos_i >= 0) \
        & (sos_i <= peak) & (eos_i >= peak)
    len_prev = 366 if calendar.isleap(year - 1) else 365

    def code(tt):
        tf = xp.floor(tt)
        yy = xp.where(tf < 0, (year - 1) % 100, xp.where(tf >= len_y, (year + 1) % 100, year % 100))
        doy = xp.where(tf < 0, tf + len_prev + 1, xp.where(tf >= len_y, tf - len_y + 1, tf + 1))
        return yy * 1000 + doy

    out = {"AMPL": ampl, "LENGTH": length, "RSLOPE": rslope, "TPROD": tprod, "MAXV": maxv,
           "MINV": minv, "SOSD": code(sos_t), "EOSD": code(eos_t)}
    return {k: xp.where(valid, v.astype("float32"), xp.nan) for k, v in out.items()}


# --------------------------------------------------------------------------- driver
def run_region(year: int, crs: str, x0: float, y1: float, width: int, height: int, out_tif: str,
               delta_bounds: Optional[tuple] = None, block_rows: int = 512, workers: int = 32,
               p: Optional[ReconParams] = None, use_gpu: bool = True, gpu_chunk: int = 400_000,
               log: Callable = print) -> str:
    """Reconstruct one year for the grid (crs, upper-left x0/y1, 10 m, width x height)."""
    p = p or ReconParams()
    xp = xp_module(use_gpu)
    res = 10.0
    bounds = (x0, y1 - height * res, x0 + width * res, y1)
    bbox_ll = transform_bounds(crs, "EPSG:4326", *bounds)
    centres, ws, nbins = pnp.window_bins(year, p.bin_days)
    start, end = ws.isoformat(), (dt.date(year + 1, 4, 1) - dt.timedelta(days=1)).isoformat()
    t0 = time.time()
    scenes = search_scenes(bbox_ll, start, end)
    jan1 = dt.datetime(year, 1, 1)
    ws_dt = dt.datetime(ws.year, ws.month, ws.day)
    obs_bin_all = np.array([int((s.when - ws_dt).total_seconds() // (86400 * p.bin_days)) for s in scenes])
    obs_bin_all = np.where((obs_bin_all >= 0) & (obs_bin_all < nbins), obs_bin_all, -1)
    log(f"{year}: {len(scenes)} scenes ({time.time() - t0:.0f} s)")

    # pre-pass: tile-constant MDVI term from 160 m overviews over delta_bounds (default: the region)
    db = delta_bounds or bounds
    cw, ch = int((db[2] - db[0]) // 200), int((db[3] - db[1]) // 200)
    ctr = from_origin(db[0], db[3], 200, 200)
    dbll = transform_bounds(crs, "EPSG:4326", *db)
    sel = [s for s in scenes if _overlap(s.bbox, dbll)]
    with cf.ThreadPoolExecutor(workers) as ex:
        coarse = list(ex.map(lambda s: read_dvi(s, crs, ctr, cw, ch, coarse=True), sel))
    cst = np.stack(coarse)
    delta = float(np.nanmean(np.nanmax(cst, axis=0) - np.nanpercentile(cst, 95, axis=0)))
    log(f"{year}: MDVI delta {delta:.4f} from {len(sel)} scenes at 200 m ({time.time() - t0:.0f} s)")

    prof = dict(driver="GTiff", width=width, height=height, count=len(BANDS_OUT), dtype="int16",
                crs=crs, transform=from_origin(x0, y1, res, res), nodata=NODATA, tiled=True,
                blockxsize=512, blockysize=512, compress="deflate", BIGTIFF="IF_SAFER")
    with rasterio.open(out_tif, "w", **prof) as dst:
        for i, b in enumerate(BANDS_OUT, 1):
            dst.set_band_description(i, b)
        for r0 in range(0, height, block_rows):
            tb0 = time.time()
            h = min(block_rows, height - r0)
            btr = from_origin(x0, y1 - r0 * res, res, res)
            bb = (x0, y1 - (r0 + h) * res, x0 + width * res, y1 - r0 * res)
            bll = transform_bounds(crs, "EPSG:4326", *bb)
            ids = [k for k, s in enumerate(scenes) if obs_bin_all[k] >= 0 and _overlap(s.bbox, bll)]
            with cf.ThreadPoolExecutor(workers) as ex:
                stack = list(ex.map(lambda k: read_dvi(scenes[k], crs, btr, width, h), ids))
            dvi = np.stack(stack).reshape(len(ids), -1).astype("float16")
            # drop scenes without any valid pixel in this block
            keep_obs = ~np.isnan(dvi).all(axis=1)
            dvi, ids = dvi[keep_obs], [k for k, m in zip(ids, keep_obs) if m]
            sza = np.array([scenes[k].sza for k in ids], dtype="float32")
            ob = obs_bin_all[ids]
            tread = time.time() - tb0
            P = width * h
            out = {b: np.full(P, np.nan, dtype="float32") for b in BANDS_OUT}
            for c0 in range(0, P, gpu_chunk):
                c1 = min(P, c0 + gpu_chunk)
                d = xp.asarray(dvi[:, c0:c1], dtype="float32")
                q95 = nanquantile(xp, d, 0.95, axis=0)
                mdvi = xp.clip(q95 + delta, None, p.mdvi_cap)
                mdvi = xp.maximum(mdvi, p.dvi_soil + 0.1)
                mdvi = xp.where(xp.isnan(q95), p.dvi_soil + 0.1, mdvi)
                ppi = ppi_from_dvi(xp, d, sza, mdvi, p)
                m = metrics_from_ppi(xp, ppi, ob, year, p)
                for b in BANDS_OUT:
                    v = m[b]
                    out[b][c0:c1] = v.get() if hasattr(v, "get") else v
                del d, ppi, m
            for i, b in enumerate(BANDS_OUT, 1):
                v = out[b].reshape(h, width)
                f = GEE_EXPORT_FACTORS.get(b)
                dn = np.where(np.isnan(v), NODATA, np.round(v * f) if f else np.round(v))
                if b == "NOBS":
                    dn = np.where(np.isnan(out["AMPL"].reshape(h, width)), NODATA, np.round(v))
                dst.write(np.clip(dn, -32768, 32767).astype("int16"), i,
                          window=rasterio.windows.Window(0, r0, width, h))
            log(f"{year}: rows {r0}-{r0 + h} of {height}: {len(ids)} scenes, read {tread:.0f} s, "
                f"total {time.time() - tb0:.0f} s")
    log(f"{year}: done in {(time.time() - t0) / 60:.1f} min -> {out_tif}")
    return out_tif


def _overlap(a, b):
    return not (a[2] < b[0] or a[0] > b[2] or a[3] < b[1] or a[1] > b[3])

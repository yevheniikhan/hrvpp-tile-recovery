"""Validation of the reconstruction against official HR-VPP tiles.

The missing tile T34UGV overlaps its official neighbours T34UFV (west), T34UGA (north)
and T34UGU (south) by ~10 km. In those strips the reconstructed and the official
pixels sit on the *same* 10 m Sentinel-2 grid, so they can be compared pixel by pixel
without resampling. This is the core accuracy test (docs/VALIDATION.md).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import rasterio
from rasterio.windows import from_bounds

from .config import GEE_EXPORT_FACTORS

GEE_NODATA = -32768


def metrics(ref: np.ndarray, est: np.ndarray) -> Dict[str, float]:
    m = np.isfinite(ref) & np.isfinite(est)
    ref, est = ref[m], est[m]
    n = int(m.sum())
    if n < 10:
        return {"n": n}
    d = est - ref
    r = float(np.corrcoef(ref, est)[0, 1])
    slope, intercept = np.polyfit(ref, est, 1)
    mean_ref = float(np.mean(ref))
    rmse = float(np.sqrt(np.mean(d ** 2)))
    return {
        "n": n, "mean_ref": mean_ref, "mean_est": float(np.mean(est)), "bias": float(np.mean(d)),
        "mae": float(np.mean(np.abs(d))), "rmse": rmse,
        "rrmse_pct": float(100 * rmse / abs(mean_ref)) if mean_ref else np.nan,
        "r": r, "r2": r * r, "ols_slope": float(slope), "ols_intercept": float(intercept),
        "median_abs_err": float(np.median(np.abs(d))),
    }


def overlap_arrays(recon_tif: Path, recon_band: str, official_tif: Path,
                   max_pixels: Optional[int] = 2_000_000, seed: int = 42):
    """Return (official, reconstructed) physical values on the common pixels."""
    with rasterio.open(official_tif) as off, rasterio.open(recon_tif) as rec:
        if off.crs != rec.crs:
            raise ValueError(f"CRS differ: {off.crs} vs {rec.crs}")
        b1, b2 = off.bounds, rec.bounds
        left, bottom = max(b1.left, b2.left), max(b1.bottom, b2.bottom)
        right, top = min(b1.right, b2.right), min(b1.top, b2.top)
        if right <= left or top <= bottom:
            raise ValueError("no overlap")
        w_off = from_bounds(left, bottom, right, top, off.transform).round_offsets().round_lengths()
        w_rec = from_bounds(left, bottom, right, top, rec.transform).round_offsets().round_lengths()
        a = off.read(1, window=w_off).astype("float64")
        if off.nodata is not None:
            a[a == off.nodata] = np.nan
        a = a * (off.scales[0] or 1.0) + (off.offsets[0] or 0.0)
        bi = [i for i, d in enumerate(rec.descriptions, 1) if d == recon_band]
        bi = bi[0] if bi else 1
        b = rec.read(bi, window=w_rec).astype("float64")
        if recon_band in GEE_EXPORT_FACTORS and rec.descriptions[bi - 1] == recon_band:
            b[b == GEE_NODATA] = np.nan
            b = b / GEE_EXPORT_FACTORS[recon_band]
        else:  # already a MONIFUN-encoded product file
            if rec.nodata is not None:
                b[b == rec.nodata] = np.nan
            b = b * (rec.scales[bi - 1] or 1.0) + (rec.offsets[bi - 1] or 0.0)
    h, w = min(a.shape[0], b.shape[0]), min(a.shape[1], b.shape[1])
    a, b = a[:h, :w].ravel(), b[:h, :w].ravel()
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if max_pixels and a.size > max_pixels:
        idx = np.random.default_rng(seed).choice(a.size, max_pixels, replace=False)
        a, b = a[idx], b[idx]
    return a, b


def coverage_stats(official_tif: Path, recon_tif: Path, recon_band: str) -> Dict[str, float]:
    """Share of overlap pixels that have a value in official / reconstructed / both."""
    with rasterio.open(official_tif) as off, rasterio.open(recon_tif) as rec:
        b1, b2 = off.bounds, rec.bounds
        bounds = (max(b1.left, b2.left), max(b1.bottom, b2.bottom), min(b1.right, b2.right), min(b1.top, b2.top))
        a = off.read(1, window=from_bounds(*bounds, off.transform).round_offsets().round_lengths())
        bi = [i for i, d in enumerate(rec.descriptions, 1) if d == recon_band] or [1]
        b = rec.read(bi[0], window=from_bounds(*bounds, rec.transform).round_offsets().round_lengths())
        h, w = min(a.shape[0], b.shape[0]), min(a.shape[1], b.shape[1])
        va = a[:h, :w] != off.nodata
        vb = b[:h, :w] != (GEE_NODATA if rec.descriptions[bi[0] - 1] == recon_band else rec.nodata)
    tot = va.size
    return {"pixels": int(tot), "official_valid_pct": 100 * va.mean(), "recon_valid_pct": 100 * vb.mean(),
            "both_valid_pct": 100 * (va & vb).mean(), "official_only_pct": 100 * (va & ~vb).mean(),
            "recon_only_pct": 100 * (~va & vb).mean()}


def strip_calibration(recon_tif: Path, official_root: Path, year: int, products, neighbours,
                      max_pixels: int = 1_000_000):
    """Fit official ~ k * recon + c per product on all overlap strips of `year`, and report a
    leave-one-strip-out check (fit on two strips, test on the third).

    Returns (fits: {product: (k, c)}, loso: list of dicts, raw: list of dicts).
    """
    store, raw = {}, []
    for nb in neighbours:
        for prod in products:
            hits = sorted((Path(official_root) / prod / str(year)).glob(f"*{nb}*_s1_{prod}.tif"))
            if not hits:
                continue
            a, b = overlap_arrays(recon_tif, prod, hits[0], max_pixels=max_pixels)
            store[(nb, prod)] = (a, b)
            m = metrics(a, b)
            m.update(neighbour=nb, product=prod, year=year)
            raw.append(m)
    fits, loso = {}, []
    for prod in products:
        nbs = [nb for nb in neighbours if (nb, prod) in store]
        if not nbs:
            continue
        a = np.concatenate([store[(n, prod)][0] for n in nbs])
        b = np.concatenate([store[(n, prod)][1] for n in nbs])
        k, c = np.polyfit(b, a, 1)
        fits[prod] = (float(k), float(c))
        for test in nbs:
            tr = [n for n in nbs if n != test]
            if not tr:
                continue
            ka, ca = np.polyfit(np.concatenate([store[(n, prod)][1] for n in tr]),
                                np.concatenate([store[(n, prod)][0] for n in tr]), 1)
            ref, est = store[(test, prod)]
            r0, r1 = metrics(ref, est), metrics(ref, ka * est + ca)
            loso.append({"year": year, "product": prod, "test_strip": test, "k": ka, "c": ca, "r": r0["r"],
                         "bias_raw": r0["bias"], "bias_cal": r1["bias"], "rmse_raw": r0["rmse"],
                         "rmse_cal": r1["rmse"], "rrmse_raw": r0["rrmse_pct"], "rrmse_cal": r1["rrmse_pct"]})
    return fits, loso, raw

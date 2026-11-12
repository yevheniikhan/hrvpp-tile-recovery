"""Earth Engine reconstruction of HR-VPP season-1 parameters for one Sentinel-2 tile.

Pipeline (per target year Y, window 1 Oct Y-1 .. 31 Mar Y+1):

1. Sentinel-2 L2A (COPERNICUS/S2_SR_HARMONIZED); for scenes without an L2A counterpart
   (mostly 2017-2018) the L1C scene (COPERNICUS/S2_HARMONIZED) is used and flagged.
2. Cloud/shadow masking with Cloud Score+ (cs_cdf >= threshold) and, for L2A, SCL.
3. Plant Phenology Index (Jin & Eklundh 2014) exactly as in the HR-VPP ATBD, Sec. 5.2.1:
      PPI = -K * ln((MDVI - DVI) / (MDVI - DVIs)),  DVI = B08 - B04
      K   = 1/(4*QE) * (1 + MDVI)/(1 - MDVI),  QE = dc + (1 - dc) * G / cos(sza)
      dc  = 0.0336 + 0.0477 / cos(sza),  G = 0.5,  MDVI <= 0.8,  PPI clipped to [0, 3]
4. 10-day median composites; winter bins (Oct-Mar) without data filled with the 5th
   percentile base level (ATBD); remaining gaps linearly interpolated (<= max_gap_days).
5. Binomial smoothing in time; season metrics with amplitude thresholds of 25 % (SOS)
   and 15 % (EOS) as in the HR-VPP Product User Manual. The season is the one whose
   peak falls in year Y.

Differences to the operational HR-VPP chain (TIMESAT double-logistic fitting, tile-wise
MDVI/DVIs interpolation, QFLAG2 masking) are documented in docs/METHOD.md and quantified
in the validation (docs/VALIDATION.md).
"""
from __future__ import annotations

import calendar
import datetime as dt
from typing import Dict, List, Optional

import ee

from .config import GEE_EXPORT_FACTORS, MISSING_TILE, Grid, ReconParams, S2Tile

S2_SR = "COPERNICUS/S2_SR_HARMONIZED"
S2_L1C = "COPERNICUS/S2_HARMONIZED"
CS_PLUS = "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED"


# --------------------------------------------------------------------------- geometry
def tile_geometry(tile: S2Tile = MISSING_TILE) -> ee.Geometry:
    xmin, ymin, xmax, ymax = tile.bounds
    return ee.Geometry.Rectangle([xmin, ymin, xmax, ymax], f"EPSG:{tile.epsg}", False)


# --------------------------------------------------------------------------- inputs
def _prep(img: ee.Image, src: int, p: ReconParams) -> ee.Image:
    clear = img.select("cs_cdf").gte(p.cs_threshold)
    if src == 1:  # L2A: also drop SCL no-data, saturated, shadow, cloud, cirrus
        scl = img.select("SCL")
        clear = clear.And(scl.neq(0)).And(scl.neq(1)).And(scl.neq(3)) \
                     .And(scl.neq(8)).And(scl.neq(9)).And(scl.neq(10))
    red = img.select("B4").divide(10000)
    nir = img.select("B8").divide(10000)
    dvi = nir.subtract(red).rename("DVI")
    sza = ee.Number(img.get("MEAN_SOLAR_ZENITH_ANGLE"))
    out = (dvi.updateMask(clear)
           .addBands(ee.Image.constant(src).toByte().rename("SRC"))
           .set("sza", sza))
    return ee.Image(out.copyProperties(img, ["system:time_start", "system:index"]))


def s2_dvi_collection(region: ee.Geometry, start: str, end: str, p: ReconParams) -> ee.ImageCollection:
    """Cloud-masked DVI images, L2A preferred, L1C only where no L2A scene exists."""
    cs = ee.ImageCollection(CS_PLUS)
    sr = (ee.ImageCollection(S2_SR).filterBounds(region).filterDate(start, end)
          .linkCollection(cs, ["cs_cdf"]))
    l1c = (ee.ImageCollection(S2_L1C).filterBounds(region).filterDate(start, end)
           .filter(ee.Filter.inList("system:index", sr.aggregate_array("system:index")).Not())
           .linkCollection(cs, ["cs_cdf"]))
    a = sr.map(lambda i: _prep(i, 1, p))
    b = l1c.map(lambda i: _prep(i, 2, p))
    return a.merge(b).filter(ee.Filter.notNull(["sza"]))


def mdvi_image(dvi_col: ee.ImageCollection, region: ee.Geometry, p: ReconParams) -> ee.Image:
    """MDVI = DVI_q95 + mean_tile(DVI_max - DVI_q95), capped at 0.8 (ATBD 5.2.1)."""
    q95 = dvi_col.select("DVI").reduce(ee.Reducer.percentile([95])).rename("q95")
    dmax = dvi_col.select("DVI").max().rename("dmax")
    delta = ee.Number(dmax.subtract(q95).reduceRegion(
        reducer=ee.Reducer.mean(), geometry=region, scale=200, maxPixels=1e9,
        bestEffort=True, tileScale=4).values().get(0))
    return q95.add(delta).min(p.mdvi_cap).max(p.dvi_soil + 0.1).rename("MDVI")


def add_ppi(dvi_col: ee.ImageCollection, mdvi: ee.Image, p: ReconParams) -> ee.ImageCollection:
    dvis = p.dvi_soil

    def f(img):
        cosz = ee.Number(img.get("sza")).multiply(3.141592653589793 / 180).cos()
        dc = ee.Number(0.0336).add(ee.Number(0.0477).divide(cosz))
        qe = dc.add(ee.Number(1).subtract(dc).multiply(p.leaf_g).divide(cosz))
        k = ee.Image(1).divide(ee.Image.constant(qe).multiply(4)) \
            .multiply(mdvi.add(1)).divide(ee.Image(1).subtract(mdvi))
        dvi = img.select("DVI").min(mdvi.subtract(1e-4)).max(dvis)
        ppi = k.multiply(mdvi.subtract(dvi).divide(mdvi.subtract(dvis)).log()).multiply(-1)
        out = ppi.clamp(0, 3).rename("PPI").updateMask(img.select("DVI").mask()) \
            .addBands(img.select("SRC"))
        return ee.Image(out.copyProperties(img, ["system:time_start"]))

    return dvi_col.map(f)


# --------------------------------------------------------------------------- time series
def window_bins(year: int, bin_days: int = 10):
    ws = dt.date(year - 1, 10, 1)
    we = dt.date(year + 1, 3, 31)
    jan1 = dt.date(year, 1, 1)
    nbins = ((we - ws).days + 1) // bin_days
    centres = [(ws - jan1).days + bin_days / 2.0 + bin_days * i for i in range(nbins)]
    return ws, nbins, centres


def dekadal_series(ppi_col: ee.ImageCollection, year: int, p: ReconParams) -> ee.ImageCollection:
    ws, nbins, centres = window_bins(year, p.bin_days)
    jan1 = dt.datetime(year, 1, 1)
    base = ppi_col.select("PPI").reduce(ee.Reducer.percentile([p.winter_base_pct])).rename("PPI")
    empty = ee.Image.constant(0).toFloat().rename("PPI").updateMask(ee.Image.constant(0))
    imgs = []
    for i, c in enumerate(centres):
        s = dt.datetime(ws.year, ws.month, ws.day) + dt.timedelta(days=i * p.bin_days)
        e = s + dt.timedelta(days=p.bin_days)
        centre = jan1 + dt.timedelta(days=c)
        sub = ppi_col.filterDate(s.strftime("%Y-%m-%d"), e.strftime("%Y-%m-%d")).select("PPI")
        comp = ee.Image(ee.Algorithms.If(sub.size().gt(0), sub.median(), empty)).rename("PPI").toFloat()
        nobs = ee.Image(ee.Algorithms.If(sub.size().gt(0), sub.count(), ee.Image.constant(0))) \
            .rename("N").unmask(0)
        if centre.month in p.winter_months:
            comp = comp.unmask(base)
        millis = int((centre - dt.datetime(1970, 1, 1)).total_seconds() * 1000)
        tband = ee.Image.constant(c).toFloat().rename("T").updateMask(comp.mask())
        imgs.append(comp.addBands(tband).addBands(ee.Image(nobs).toFloat())
                    .set({"system:time_start": millis, "t": c, "i": i}))
    return ee.ImageCollection.fromImages(imgs)


def interpolate_gaps(col: ee.ImageCollection, max_gap_days: int) -> ee.ImageCollection:
    """Linear temporal interpolation of masked PPI bins (join-based, nearest before/after)."""
    millis = max_gap_days * 86400000
    diff = ee.Filter.maxDifference(difference=millis, leftField="system:time_start",
                                   rightField="system:time_start")
    after_f = ee.Filter.And(diff, ee.Filter.lessThanOrEquals(leftField="system:time_start",
                                                             rightField="system:time_start"))
    before_f = ee.Filter.And(diff, ee.Filter.greaterThanOrEquals(leftField="system:time_start",
                                                                 rightField="system:time_start"))
    j1 = ee.Join.saveAll(matchesKey="after", ordering="system:time_start", ascending=False) \
        .apply(primary=col, secondary=col, condition=after_f)
    j2 = ee.Join.saveAll(matchesKey="before", ordering="system:time_start", ascending=True) \
        .apply(primary=j1, secondary=j1, condition=before_f)

    def f(img):
        img = ee.Image(img)
        bm = ee.ImageCollection.fromImages(ee.List(img.get("before"))).select(["PPI", "T"]).mosaic()
        am = ee.ImageCollection.fromImages(ee.List(img.get("after"))).select(["PPI", "T"]).mosaic()
        t = ee.Image.constant(img.get("t")).toFloat()
        span = am.select("T").subtract(bm.select("T"))
        ratio = t.subtract(bm.select("T")).divide(span.where(span.eq(0), 1))
        interp = bm.select("PPI").add(am.select("PPI").subtract(bm.select("PPI")).multiply(ratio))
        out = img.select("PPI").unmask(interp).addBands(img.select("N"))
        return ee.Image(out.copyProperties(img, ["system:time_start", "t", "i"]))

    return ee.ImageCollection(j2.map(f)).sort("system:time_start")


# --------------------------------------------------------------------------- metrics
def _smooth121(a: ee.Image, n: int) -> ee.Image:
    left = a.arraySlice(0, 0, 1).arrayCat(a.arraySlice(0, 0, n - 1), 0)
    right = a.arraySlice(0, 1, n).arrayCat(a.arraySlice(0, n - 1, n), 0)
    return left.add(a.multiply(2)).add(right).divide(4)


def smoothed_array(series: ee.ImageCollection, year: int, p: ReconParams) -> ee.Image:
    """Gap-free, smoothed PPI series as a 1-D array image (band 'PPI', length n)."""
    _, n, _ = window_bins(year, p.bin_days)
    complete = series.select("PPI").count().eq(n)      # masked bins -> 0, tracked separately
    V = series.select("PPI").map(lambda i: i.unmask(0)).toArray().arrayProject([0])
    V = V.updateMask(complete)
    for _ in range(p.smooth_passes):
        V = _smooth121(V, n)
    return V


def season_metrics(series: ee.ImageCollection, year: int, p: ReconParams) -> ee.Image:
    """Array implementation of phenology_np.extract_metrics (keep both in sync!)."""
    _, n, centres = window_bins(year, p.bin_days)
    len_y = 366 if calendar.isleap(year) else 365
    len_prev = 366 if calendar.isleap(year - 1) else 365
    half = p.bin_days / 2.0

    V = smoothed_array(series, year, p)
    nobs = series.select("N").sum().rename("NOBS")

    idx = ee.Image(ee.Array(list(range(n))))
    T = ee.Image(ee.Array([float(c) for c in centres]))
    in_year = T.gte(0).And(T.lt(len_y))

    def amin(a):
        return a.arrayReduce(ee.Reducer.min(), [0]).arrayGet([0])

    def amax(a):
        return a.arrayReduce(ee.Reducer.max(), [0]).arrayGet([0])

    def pick(cond, a, other):          # elementwise cond ? a : other (arrays/scalars)
        return cond.multiply(a).add(cond.Not().multiply(other))

    Vm = pick(in_year, V, -1)
    maxv = amax(Vm)
    peak = amin(pick(Vm.eq(maxv), idx, 999)).toInt()

    left = idx.lte(peak)
    right = idx.gte(peak)
    lmin = amin(pick(left, V, 99))
    rmin = amin(pick(right, V, 99))
    lmin_i = amax(pick(left.And(V.eq(lmin)), idx, -1)).toInt()
    rmin_i = amin(pick(right.And(V.eq(rmin)), idx, 999)).toInt()

    minv = lmin.add(rmin).divide(2)
    ampl = maxv.subtract(minv)
    thr_s = lmin.add(maxv.subtract(lmin).multiply(p.sos_frac))
    thr_e = rmin.add(maxv.subtract(rmin).multiply(p.eos_frac))

    rng_s = idx.gte(lmin_i).And(left)
    rng_e = idx.lte(rmin_i).And(right)
    sos_i = amin(pick(rng_s.And(V.gte(thr_s)), idx, 999)).toInt()
    eos_i = amax(pick(rng_e.And(V.gte(thr_e)), idx, -1)).toInt()

    s0 = sos_i.subtract(1).max(lmin_i).clamp(0, n - 1).toInt()
    s1 = sos_i.clamp(0, n - 1).toInt()
    e0 = eos_i.clamp(0, n - 1).toInt()
    e1 = eos_i.add(1).min(rmin_i).clamp(0, n - 1).toInt()
    vs0, vs1 = V.arrayGet(s0), V.arrayGet(s1)
    ve0, ve1 = V.arrayGet(e0), V.arrayGet(e1)

    ds = vs1.subtract(vs0)
    fs = thr_s.subtract(vs0).divide(ds.where(ds.eq(0), 1)).where(ds.eq(0), 1).clamp(0, 1)
    de = ve0.subtract(ve1)
    fe = ve0.subtract(thr_e).divide(de.where(de.eq(0), 1)).where(de.eq(0), 0).clamp(0, 1)
    t0 = float(centres[0])
    sos_t = s0.add(fs.multiply(s1.subtract(s0))).multiply(p.bin_days).add(t0)
    eos_t = e0.add(fe.multiply(e1.subtract(e0))).multiply(p.bin_days).add(t0)
    length = eos_t.subtract(sos_t)

    lo = T.subtract(half).max(sos_t)
    hi = T.add(half).min(eos_t)
    w = hi.subtract(lo).max(0)
    tprod = V.multiply(w).arrayReduce(ee.Reducer.sum(), [0]).arrayGet([0])
    sprod = tprod.subtract(minv.multiply(length))
    lslope = vs1.subtract(vs0).abs().divide(p.bin_days)
    rslope = ve0.subtract(ve1).abs().divide(p.bin_days)

    def code(tt):
        tf = tt.floor()
        yy = ee.Image(year % 100).where(tf.lt(0), (year - 1) % 100).where(tf.gte(len_y), (year + 1) % 100)
        doy = tf.add(1).where(tf.lt(0), tf.add(len_prev + 1)).where(tf.gte(len_y), tf.subtract(len_y - 1))
        return yy.multiply(1000).add(doy)

    valid = ampl.gte(p.min_ampl).And(length.gte(p.min_length)).And(sos_i.lt(999)) \
        .And(eos_i.gte(0)).And(sos_i.lte(peak)).And(eos_i.gte(peak)).And(nobs.gte(p.min_obs))

    out = ee.Image.cat([
        ampl.rename("AMPL"), length.rename("LENGTH"), rslope.rename("RSLOPE"), tprod.rename("TPROD"),
        sprod.rename("SPROD"), maxv.rename("MAXV"), minv.rename("MINV"), lslope.rename("LSLOPE"),
        code(sos_t).rename("SOSD"), code(eos_t).rename("EOSD"), sos_t.rename("SOS_T"),
        eos_t.rename("EOS_T"),
    ]).updateMask(valid)
    return out.addBands(nobs).set({"year": year})


# --------------------------------------------------------------------------- driver
def reconstruct_year(year: int, tile: S2Tile = MISSING_TILE, p: Optional[ReconParams] = None,
                     region: Optional[ee.Geometry] = None) -> Dict[str, ee.Image]:
    """Return dict with 'metrics' (season-1 parameters), 'series' and 'src' for a year."""
    p = p or ReconParams()
    region = region or tile_geometry(tile)
    ws, _, _ = window_bins(year, p.bin_days)
    start = ws.strftime("%Y-%m-%d")
    end = dt.date(year + 1, 4, 1).strftime("%Y-%m-%d")
    dvi = s2_dvi_collection(region, start, end, p)
    mdvi = mdvi_image(dvi, region, p)
    ppi = add_ppi(dvi, mdvi, p)
    series = interpolate_gaps(dekadal_series(ppi, year, p), p.max_gap_days)
    metrics = season_metrics(series, year, p)
    l1c_share = ppi.select("SRC").map(lambda i: i.eq(2)).mean().rename("L1C_SHARE")
    return {"metrics": metrics.addBands(l1c_share), "series": series, "mdvi": mdvi,
            "smoothed": smoothed_array(series, year, p), "n_scenes": dvi.size(), "params": p}


def export_region(year: int, folder: str, grid: Grid, p: Optional[ReconParams] = None,
                  products: Optional[List[str]] = None, prefix: str = "RECON",
                  mdvi_region: Optional[ee.Geometry] = None, start: bool = True) -> ee.batch.Task:
    """Drive export of the season-1 parameters for any area on a 10 m UTM grid.

    Int16 DN with the factors of config.GEE_EXPORT_FACTORS, nodata -32768, bands
    AMPL, LENGTH, RSLOPE, TPROD, SOSD, EOSD, MAXV, MINV, NOBS (same as the Code Editor script).
    Large areas are split by Earth Engine into row shards; scripts/mosaic_shards.py merges them.
    """
    products = products or ["AMPL", "LENGTH", "RSLOPE", "TPROD", "SOSD", "EOSD", "MAXV", "MINV"]
    region = ee.Geometry.Rectangle(list(grid.bounds), grid.crs, False)
    m = reconstruct_year(year, p=p, region=mdvi_region or region)["metrics"]
    bands = []
    for prod in products:
        b = m.select(prod)
        f = GEE_EXPORT_FACTORS.get(prod)
        bands.append(b.multiply(f).round() if f else b.round())
    img = ee.Image.cat(bands).addBands(m.select("NOBS")).unmask(-32768).toInt16()
    task = ee.batch.Export.image.toDrive(
        image=img, description=f"{prefix}_{grid.name}_{year}", folder=folder,
        fileNamePrefix=f"{prefix}_VPP_{year}_S2_{grid.name}-{int(grid.res):03d}m_s1",
        crs=grid.crs, crsTransform=[grid.res, 0, grid.x0, 0, -grid.res, grid.y1],
        dimensions=f"{grid.width}x{grid.height}", maxPixels=1e10,
        fileDimensions=[((grid.width + 255) // 256) * 256, ((grid.height + 255) // 256) * 256],
        formatOptions={"cloudOptimized": True})
    if start:
        task.start()
    return task


def export_year(year: int, folder: str, tile: S2Tile = MISSING_TILE,
                products: Optional[List[str]] = None, p: Optional[ReconParams] = None,
                scale: float = 10.0, prefix: str = "RECON") -> ee.batch.Task:
    """Start a Drive export of the reconstructed parameters for a whole Sentinel-2 tile."""
    products = products or ["AMPL", "LENGTH", "RSLOPE", "TPROD", "SOSD", "EOSD"]
    m = reconstruct_year(year, tile, p)["metrics"]
    bands = []
    for prod in products:
        b = m.select(prod)
        f = GEE_EXPORT_FACTORS.get(prod)
        bands.append(b.multiply(f).round() if f else b.round())
    img = ee.Image.cat(bands).addBands(m.select("NOBS")).unmask(-32768).toInt16()
    span = tile.size_px * tile.res
    dims = int(round(span / scale))
    task = ee.batch.Export.image.toDrive(
        image=img, description=f"{prefix}_{tile.tile_id}_{year}", folder=folder,
        fileNamePrefix=f"{prefix}_VPP_{year}_S2_{tile.tile_id}-{int(scale):03d}m_s1",
        crs=f"EPSG:{tile.epsg}", crsTransform=[scale, 0, tile.ulx, 0, -scale, tile.uly],
        dimensions=f"{dims}x{dims}", maxPixels=1e10,
        fileDimensions=[((dims + 255) // 256) * 256] * 2, formatOptions={"cloudOptimized": True})
    task.start()
    return task


def sample_points(img: ee.Image, points: ee.FeatureCollection, tile: S2Tile = MISSING_TILE,
                  scale: float = 10.0) -> ee.FeatureCollection:
    """Sample on the exact Sentinel-2 grid of `tile` (no resampling)."""
    return img.sampleRegions(collection=points, scale=scale, projection=ee.Projection(f"EPSG:{tile.epsg}"),
                             geometries=True, tileScale=4)

"""Fast validation and Earth-Engine-vs-NumPy cross-check by point sampling.

No full-tile export is needed: reconstructed values are sampled in Earth Engine at
pixel centres inside the overlap strips, official values are read from the downloaded
GeoTIFFs at exactly the same coordinates.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import ee
import numpy as np
import pandas as pd
import rasterio

from .config import MISSING_TILE, NEIGHBOURS, S2Tile
from .phenology_np import extract_metrics, window_bins
from .validate import metrics


def strip_points(neighbour: S2Tile, tile: S2Tile = MISSING_TILE, n: int = 3000, seed: int = 7,
                 margin_px: int = 5) -> np.ndarray:
    """Random 10 m pixel centres inside the overlap of `tile` and `neighbour` (UTM metres)."""
    a, b = tile.bounds, neighbour.bounds
    xmin, ymin = max(a[0], b[0]) + margin_px * 10, max(a[1], b[1]) + margin_px * 10
    xmax, ymax = min(a[2], b[2]) - margin_px * 10, min(a[3], b[3]) - margin_px * 10
    rng = np.random.default_rng(seed)
    cols = rng.integers(0, int((xmax - xmin) // 10), n)
    rows = rng.integers(0, int((ymax - ymin) // 10), n)
    return np.column_stack([xmin + 10 * cols + 5, ymax - 10 * rows - 5])


def to_fc(xy: np.ndarray, epsg: int) -> ee.FeatureCollection:
    proj = f"EPSG:{epsg}"
    return ee.FeatureCollection([
        ee.Feature(ee.Geometry.Point([float(x), float(y)], proj), {"pid": i}) for i, (x, y) in enumerate(xy)])


def sample_ee(img: ee.Image, xy: np.ndarray, epsg: int, bands: List[str], chunk: int = 1000) -> pd.DataFrame:
    rows = []
    for s in range(0, len(xy), chunk):
        fc = to_fc(xy[s:s + chunk], epsg)
        res = img.select(bands).sampleRegions(collection=fc, scale=10, projection=ee.Projection(f"EPSG:{epsg}"),
                                              geometries=False, tileScale=4).getInfo()
        for f in res["features"]:
            p = f["properties"]
            p["pid"] = int(p["pid"]) + s
            rows.append(p)
    df = pd.DataFrame(rows)
    full = pd.DataFrame({"pid": np.arange(len(xy)), "x": xy[:, 0], "y": xy[:, 1]})
    return full.merge(df, on="pid", how="left")


def sample_official(path: Path, xy: np.ndarray) -> np.ndarray:
    with rasterio.open(path) as src:
        v = np.array([s[0] for s in src.sample([(float(x), float(y)) for x, y in xy])], dtype="float64")
        if src.nodata is not None:
            v[v == src.nodata] = np.nan
        return v * (src.scales[0] or 1.0) + (src.offsets[0] or 0.0)


def find_official(root: Path, product: str, year: int, tile_id: str) -> Optional[Path]:
    hits = sorted((root / product / str(year)).glob(f"*{tile_id}*_s1_{product}*.tif"))
    return hits[0] if hits else None


def fast_validation(metrics_img: ee.Image, year: int, official_root: Path, products: List[str],
                    n: int = 3000, seed: int = 7) -> pd.DataFrame:
    out = []
    for nb_id, nb in NEIGHBOURS.items():
        xy = strip_points(nb, n=n, seed=seed)
        rec = sample_ee(metrics_img, xy, MISSING_TILE.epsg, products)
        for prod in products:
            f = find_official(official_root, prod, year, nb_id)
            if f is None:
                out.append({"year": year, "neighbour": nb_id, "product": prod, "n": 0, "note": "official file missing"})
                continue
            ref = sample_official(f, xy)
            est = rec[prod].to_numpy(dtype="float64") if prod in rec else np.full(len(xy), np.nan)
            m = metrics(ref, est)
            m.update({"year": year, "neighbour": nb_id, "product": prod,
                      "coverage_official_pct": 100 * np.isfinite(ref).mean(),
                      "coverage_recon_pct": 100 * np.isfinite(est).mean()})
            out.append(m)
    cols = ["year", "neighbour", "product", "n", "r", "r2", "bias", "mae", "rmse", "rrmse_pct",
            "mean_ref", "mean_est", "ols_slope", "ols_intercept", "coverage_official_pct", "coverage_recon_pct"]
    df = pd.DataFrame(out)
    return df[[c for c in cols if c in df.columns] + [c for c in df.columns if c not in cols]]


def crosscheck_numpy(recon: Dict, year: int, n: int = 300, seed: int = 11) -> pd.DataFrame:
    """Compare EE season metrics with the NumPy reference applied to the EE smoothed series."""
    p = recon["params"]
    xy = strip_points(S2Tile("CENTRE", 32634, 730000.0, 5470000.0, size_px=4000), MISSING_TILE, n=n, seed=seed)
    img = recon["smoothed"].rename("V").addBands(recon["metrics"])
    bands = ["V", "AMPL", "LENGTH", "TPROD", "RSLOPE", "MAXV", "MINV", "SOS_T", "EOS_T"]
    df = sample_ee(img, xy, MISSING_TILE.epsg, bands, chunk=150).dropna(subset=["V"])
    V = np.array(df["V"].tolist(), dtype="float64")
    t, _, _ = window_bins(year, p.bin_days)
    ref = extract_metrics(V, t, year, sos_frac=p.sos_frac, eos_frac=p.eos_frac, bin_days=p.bin_days,
                          min_ampl=p.min_ampl, min_length=p.min_length)
    rows = []
    for k in ["AMPL", "LENGTH", "TPROD", "RSLOPE", "MAXV", "MINV", "SOS_T", "EOS_T"]:
        a = ref[k]
        b = df[k].to_numpy(dtype="float64")
        both = np.isfinite(a) & np.isfinite(b)
        rows.append({"metric": k, "pixels_with_series": len(df), "valid_both": int(both.sum()),
                     "valid_numpy_only": int((np.isfinite(a) & ~np.isfinite(b)).sum()),
                     "valid_ee_only": int((~np.isfinite(a) & np.isfinite(b)).sum()),
                     "max_abs_diff": float(np.nanmax(np.abs(a[both] - b[both]))) if both.any() else np.nan})
    return pd.DataFrame(rows)

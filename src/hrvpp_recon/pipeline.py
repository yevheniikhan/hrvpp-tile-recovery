"""High-level steps shared by the command-line scripts (scripts/*.py).

raw reconstruction (Earth Engine export or s2local.run_region)
    -> mosaic_shards()            one GeoTIFF per year (Earth Engine splits large exports)
    -> calibrate_and_write()      learned correction on the overlap with official tiles,
                                  MONIFUN layout with the official encoding and band metadata
    -> compare_raw() / compare_products()   agreement between two reconstructions
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
import rasterio
from rasterio.merge import merge
from rasterio.windows import Window

from . import learn, postprocess, validate
from .config import MONIFUN_PRODUCTS

GEE_NODATA = -32768


# --------------------------------------------------------------------------- files
def find_official(root: Path, product: str, year: int, tile_id: str) -> Optional[Path]:
    """Official HR-VPP file in the MONIFUN layout root/<PRODUCT>/<YEAR>/*<TILE>*_s1_<PRODUCT>.tif."""
    hits = sorted((Path(root) / product / str(year)).glob(f"*{tile_id}*_s1_{product}.tif"))
    hits = [h for h in hits if "_RECON_" not in h.name]
    return hits[0] if hits else None


def official_set(root: Path, year: int, tiles: Iterable[str],
                 products: Iterable[str] = MONIFUN_PRODUCTS) -> Dict[str, Dict[str, Path]]:
    """{tile: {product: path}} for the official tiles that exist for this year."""
    out = {}
    for t in tiles:
        d = {p: find_official(root, p, year, t) for p in products}
        d = {p: v for p, v in d.items() if v is not None}
        if d:
            out[t] = d
    return out


def product_path(out_root: Path, product: str, year: int, tile_id: str) -> Path:
    return Path(out_root) / product / str(year) / f"VPP_{year}_S2_{tile_id}-010m_RECON_s1_{product}.tif"


def mosaic_shards(parts: List[Path], out_tif: Path) -> Path:
    """Merge Earth Engine export shards (<name>-<row>-<col>.tif) into one tiled GeoTIFF, keeping band names."""
    parts = sorted(Path(p) for p in parts)
    if not parts:
        raise FileNotFoundError("no shards given")
    out_tif = Path(out_tif)
    out_tif.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(parts[0]) as s0:
        prof, desc = s0.profile.copy(), s0.descriptions
    nod = prof.get("nodata")
    nod = GEE_NODATA if nod is None else nod
    arr, tr = merge([str(p) for p in parts], nodata=nod)
    for k in ("photometric", "interleave", "blockxsize", "blockysize", "tiled", "compress"):
        prof.pop(k, None)
    prof.update(driver="GTiff", height=arr.shape[1], width=arr.shape[2], transform=tr, nodata=nod,
                tiled=True, blockxsize=512, blockysize=512, compress="deflate", BIGTIFF="IF_SAFER")
    tmp = out_tif.with_name("_" + out_tif.name)
    with rasterio.open(tmp, "w", **prof) as dst:
        dst.write(arr)
        for i, d in enumerate(desc, 1):
            if d:
                dst.set_band_description(i, d)
    tmp.replace(out_tif)
    return out_tif


# --------------------------------------------------------------------------- correction
def calibrate_and_write(raw_tif: Path, year: int, official_root: Path, neighbours: List[str],
                        out_root: Path, tile_id: str, products: Iterable[str] = MONIFUN_PRODUCTS,
                        method: str = "gbm_qm", n_per_strip: int = 300_000, report: Optional[Path] = None,
                        template_tile: Optional[str] = None, overwrite: bool = False,
                        extra_tags: Optional[Dict[str, str]] = None) -> Dict:
    """Learn the correction on the overlap of `raw_tif` with the official `neighbours` and write products.

    Returns a report with the leave-one-neighbour-out check of the methods raw / linear / `method`.
    """
    products = list(products)
    offs = official_set(official_root, year, neighbours, products)
    if not offs:
        raise RuntimeError(f"no official neighbour files for {year} under {official_root}")
    df = learn.strip_samples(raw_tif, offs, year, n_per_strip=n_per_strip)
    if df.empty:
        raise RuntimeError("the reconstruction does not overlap any of the official neighbours")
    rep = {"year": year, "raw": str(raw_tif), "neighbours": sorted(offs), "method": method,
           "samples": {k: int(v) for k, v in df["strip"].value_counts().items()}}
    if df["strip"].nunique() > 1:
        res = learn.loso(df, products, methods=sorted({"raw", "linear", method}))
        rep["loso"] = res.round(4).to_dict("records")
    cors = {p: learn.Corrector(method, p).fit(df) for p in products if "OFF_" + p in df}
    tmpl_src = official_set(official_root, year, [template_tile], products).get(template_tile, {}) \
        if template_tile else {}
    written = []
    for p, c in cors.items():
        out = product_path(out_root, p, year, tile_id)
        if out.exists() and not overwrite:
            written.append(str(out))
            continue
        tmpl = tmpl_src.get(p) or next(v[p] for v in offs.values() if p in v)
        tags = {"CALIBRATION": f"learned from official {', '.join(sorted(offs))} {year}, {method}"}
        tags.update(extra_tags or {})
        postprocess.write_product(raw_tif, p, year, out_root, tile_id, learn.phys_fn_for(c, year), tmpl,
                                  extra_tags=tags)
        written.append(str(out))
    rep["written"] = written
    if report:
        Path(report).parent.mkdir(parents=True, exist_ok=True)
        json.dump(rep, open(report, "w"), indent=1, default=float)
    return rep


# --------------------------------------------------------------------------- comparison
def compare_raw(a_tif: Path, b_tif: Path, year: int, products: Iterable[str] = MONIFUN_PRODUCTS,
                step: int = 3400, rows: int = 256) -> Dict[str, Dict[str, float]]:
    """Two raw reconstructions on the same grid, sampled in horizontal bands (ref = a, est = b)."""
    products = list(products)
    acc = {p: ([], []) for p in products}
    with rasterio.open(a_tif) as A, rasterio.open(b_tif) as B:
        if (A.width, A.height, A.transform) != (B.width, B.height, B.transform):
            raise ValueError("the two rasters are not on the same grid")
        for r0 in range(0, A.height, step):
            w = Window(0, r0, A.width, min(rows, A.height - r0))
            fa, fb = learn.read_features(A, year, w), learn.read_features(B, year, w)
            for p in products:
                acc[p][0].append(fa[p].ravel())
                acc[p][1].append(fb[p].ravel())
    return {p: validate.metrics(np.concatenate(acc[p][0]), np.concatenate(acc[p][1])) for p in products}


def compare_products(a_tif: Path, b_tif: Path, product: str, max_pixels: int = 1_000_000) -> Dict[str, float]:
    """Two MONIFUN-encoded product files (or a product file and an official tile) on their common area.

    ref = b, est = a (bias = a - b)."""
    ref, est = validate.overlap_arrays(a_tif, product, b_tif, max_pixels=max_pixels)
    return validate.metrics(ref, est)


def summarise_loso(reports: Iterable[Dict]) -> pd.DataFrame:
    """Mean held-out rRMSE per year, product and method from calibrate_and_write reports."""
    rows = []
    for r in reports:
        for x in r.get("loso", []):
            rows.append({"year": r["year"], **x})
    if not rows:
        return pd.DataFrame()
    d = pd.DataFrame(rows)
    return d.pivot_table(index=["method", "year"], columns="product", values="rrmse_pct").round(1)

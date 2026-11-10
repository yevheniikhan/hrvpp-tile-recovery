"""Convert the Earth Engine export into MONIFUN-ready, HR-VPP-encoded GeoTIFFs.

Output layout (identical to what the MONIFUN extraction scripts expect):

    <root>/Copernicus_Phenology/<PRODUCT>/<YEAR>/VPP_<YEAR>_S2_<TILE>-010m_RECON_s1_<PRODUCT>.tif

The digital encoding (dtype, scale, nodata) is copied from an official neighbouring
file whenever one is available, so reconstructed and official tiles can be mosaicked
or extracted together without any special handling.
"""
from __future__ import annotations

from pathlib import Path
import os
import re
import shutil
import tempfile
from typing import Callable, Dict, Optional, Tuple

import numpy as np
import rasterio
from rasterio.windows import Window

from .config import GEE_EXPORT_FACTORS, PRODUCT_SPECS

GEE_NODATA = -32768


def resolve_encoding(official_tif: Optional[Path], product: str) -> Dict:
    """dtype / scale / offset / nodata of the official product (fallback: config defaults)."""
    spec = PRODUCT_SPECS[product]
    enc = {"dtype": spec.dtype, "scale": spec.default_scale, "offset": 0.0, "nodata": spec.nodata,
           "source": "config-default"}
    if official_tif and Path(official_tif).exists():
        with rasterio.open(official_tif) as src:
            enc.update(dtype=src.dtypes[0], scale=float(src.scales[0] or 1.0),
                       offset=float(src.offsets[0] or 0.0), nodata=src.nodata, source=str(official_tif))
            tags = src.tags(1)
            # some HR-VPP files carry the scale only in tags
            for key in ("scale_factor", "SCALE_FACTOR", "scale"):
                if key in tags and enc["scale"] == 1.0 and product not in ("LENGTH", "SOSD", "EOSD"):
                    enc["scale"] = float(tags[key])
    return enc


def band_index(src: rasterio.io.DatasetReader, name: str) -> int:
    for i, d in enumerate(src.descriptions, start=1):
        if d == name:
            return i
    raise KeyError(f"band {name} not found in {src.name}: {src.descriptions}")


def harmonise_metadata(dst, product: str, year: int, official_tif: Optional[Path]):
    """Copy band description, units and band tags of the official product (year substituted),
    so the reconstructed band looks like the official one in GIS and extraction scripts."""
    desc, units, btags = PRODUCT_SPECS[product].description, PRODUCT_SPECS[product].unit, {}
    if official_tif and Path(official_tif).exists():
        with rasterio.open(official_tif) as src:
            if src.descriptions[0]:
                desc = re.sub(r"(19|20)\d\d", str(year), src.descriptions[0], count=1)
            units = src.units[0] or units
            btags = src.tags(1)
    dst.set_band_description(1, desc)
    if units:
        dst.set_band_unit(1, units)
    if btags:
        dst.update_tags(1, **btags)


def write_product(recon_tif: Path, product: str, year: int, out_root: Path, tile_id: str,
                  phys_fn: Callable, official_tif: Optional[Path] = None, block: int = 2048,
                  extra_tags: Optional[Dict[str, str]] = None) -> Path:
    """Write one product with the official encoding; `phys_fn(src, window)` returns physical values
    (NaN = no season)."""
    enc = resolve_encoding(official_tif, product)
    out_dir = out_root / product / str(year)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"VPP_{year}_S2_{tile_id}-010m_RECON_s1_{product}.tif"
    # write locally, then copy: block-wise GeoTIFF writes on a FUSE-mounted Drive are very slow
    tmp = Path(tempfile.gettempdir()) / f"_{os.getpid()}_{out.name}"  # unique per process
    info = np.iinfo(np.dtype(enc["dtype"]))
    with rasterio.open(recon_tif) as src:
        prof = src.profile.copy()
        prof.update(count=1, dtype=enc["dtype"], nodata=enc["nodata"], compress="deflate",
                    predictor=2, tiled=True, blockxsize=512, blockysize=512, BIGTIFF="IF_SAFER")
        prof.pop("photometric", None)
        with rasterio.open(tmp, "w", **prof) as dst:
            for r in range(0, src.height, block):
                for c in range(0, src.width, block):
                    w = Window(c, r, min(block, src.width - c), min(block, src.height - r))
                    phys = np.asarray(phys_fn(src, w), dtype="float64")
                    nod = ~np.isfinite(phys)
                    dn = np.round((np.where(nod, 0, phys) - enc["offset"]) / enc["scale"])
                    dn = np.clip(dn, info.min, info.max)
                    if enc["nodata"] is not None:
                        # never let a valid value collide with nodata
                        dn[(~nod) & (dn == enc["nodata"])] += 1 if enc["nodata"] <= 0 else -1
                        dn[nod] = enc["nodata"]
                    dst.write(dn.astype(enc["dtype"]), 1, window=w)
            dst.scales = (enc["scale"],)
            dst.offsets = (enc["offset"],)
            harmonise_metadata(dst, product, year, official_tif)
            dst.update_tags(
                PRODUCT=product, YEAR=str(year), TILE=tile_id, SEASON="1",
                SOURCE="Reconstruction from Sentinel-2 (Earth Engine) - NOT an official CLMS product",
                METHOD="HR-VPP-like PPI + 10-day composites + amplitude thresholds 25%/15%",
                ENCODING_FROM=enc["source"], SOFTWARE="hrvpp-tile-recovery", **(extra_tags or {}))
    shutil.move(str(tmp), str(out))
    return out


def gee_to_product(recon_tif: Path, product: str, year: int, out_root: Path, tile_id: str,
                   official_tif: Optional[Path] = None, block: int = 2048,
                   calibration: Optional[Tuple[float, float]] = None) -> Path:
    """Write one product of the reconstructed tile with the official HR-VPP encoding.

    `calibration` = (k, c) applies official ~= k * reconstructed + c to the physical values
    (fitted on the overlap strips of the same year, see validate.strip_calibration).
    """
    factor = GEE_EXPORT_FACTORS.get(product, 1.0)

    def phys_fn(src, w):
        a = src.read(band_index(src, product), window=w).astype("float64")
        phys = np.where(a == GEE_NODATA, np.nan, a / factor)
        return phys if calibration is None else calibration[0] * phys + calibration[1]

    cal = ("none" if calibration is None else
           f"official = {calibration[0]:.6f} * recon + {calibration[1]:.6f} (overlap strips, {year})")
    return write_product(recon_tif, product, year, out_root, tile_id, phys_fn, official_tif, block,
                         {"CALIBRATION": cal})


def read_physical(path: Path, band: int = 1, window: Optional[Window] = None) -> np.ndarray:
    """Read a single band as float with scale/offset applied and nodata -> NaN."""
    with rasterio.open(path) as src:
        a = src.read(band, window=window).astype("float64")
        if src.nodata is not None:
            a[a == src.nodata] = np.nan
        return a * (src.scales[band - 1] or 1.0) + (src.offsets[band - 1] or 0.0)

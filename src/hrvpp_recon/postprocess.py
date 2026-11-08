"""Convert the Earth Engine export into MONIFUN-ready, HR-VPP-encoded GeoTIFFs.

Output layout (identical to what the MONIFUN extraction scripts expect):

    <root>/Copernicus_Phenology/<PRODUCT>/<YEAR>/VPP_<YEAR>_S2_<TILE>-010m_RECON_s1_<PRODUCT>.tif

The digital encoding (dtype, scale, nodata) is copied from an official neighbouring
file whenever one is available, so reconstructed and official tiles can be mosaicked
or extracted together without any special handling.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Tuple

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


def gee_to_product(recon_tif: Path, product: str, year: int, out_root: Path, tile_id: str,
                   official_tif: Optional[Path] = None, block: int = 2048,
                   calibration: Optional[Tuple[float, float]] = None) -> Path:
    """Write one product of the reconstructed tile with the official HR-VPP encoding.

    `calibration` = (k, c) applies official ~= k * reconstructed + c to the physical values
    (fitted on the overlap strips of the same year, see validate.strip_calibration).
    """
    enc = resolve_encoding(official_tif, product)
    factor = GEE_EXPORT_FACTORS.get(product, 1.0)
    out_dir = out_root / product / str(year)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"VPP_{year}_S2_{tile_id}-010m_RECON_s1_{product}.tif"
    info = np.iinfo(np.dtype(enc["dtype"]))
    with rasterio.open(recon_tif) as src:
        bi = band_index(src, product)
        prof = src.profile.copy()
        prof.update(count=1, dtype=enc["dtype"], nodata=enc["nodata"], compress="deflate",
                    predictor=2, tiled=True, blockxsize=512, blockysize=512, BIGTIFF="IF_SAFER")
        prof.pop("photometric", None)
        with rasterio.open(out, "w", **prof) as dst:
            for r in range(0, src.height, block):
                for c in range(0, src.width, block):
                    w = Window(c, r, min(block, src.width - c), min(block, src.height - r))
                    a = src.read(bi, window=w).astype("float64")
                    nod = a == GEE_NODATA
                    phys = a / factor
                    if calibration is not None:
                        phys = calibration[0] * phys + calibration[1]
                    dn = np.round((phys - enc["offset"]) / enc["scale"])
                    dn = np.clip(dn, info.min, info.max)
                    dn[nod] = enc["nodata"]
                    if enc["nodata"] is not None:
                        # never let a valid value collide with nodata
                        dn[(~nod) & (dn == enc["nodata"])] += 1 if enc["nodata"] <= 0 else -1
                    dst.write(dn.astype(enc["dtype"]), 1, window=w)
            dst.scales = (enc["scale"],)
            dst.offsets = (enc["offset"],)
            dst.set_band_description(1, product)
            dst.update_tags(
                PRODUCT=product, YEAR=str(year), TILE=tile_id, SEASON="1",
                SOURCE="Reconstruction from Sentinel-2 (Earth Engine) - NOT an official CLMS product",
                METHOD="HR-VPP-like PPI + 10-day composites + amplitude thresholds 25%/15%",
                ENCODING_FROM=enc["source"], SOFTWARE="hrvpp-tile-recovery",
                CALIBRATION=("none" if calibration is None else
                             f"official = {calibration[0]:.6f} * recon + {calibration[1]:.6f} (overlap strips, {year})"))
    return out


def read_physical(path: Path, band: int = 1, window: Optional[Window] = None) -> np.ndarray:
    """Read a single band as float with scale/offset applied and nodata -> NaN."""
    with rasterio.open(path) as src:
        a = src.read(band, window=window).astype("float64")
        if src.nodata is not None:
            a[a == src.nodata] = np.nan
        return a * (src.scales[band - 1] or 1.0) + (src.offsets[band - 1] or 0.0)

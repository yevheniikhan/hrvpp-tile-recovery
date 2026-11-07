"""Static configuration: tile geometry, product definitions and encodings.

All numbers here are taken from public CLMS / VITO documentation (see docs/METHOD.md)
or derived from the Sentinel-2 MGRS tiling grid.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class S2Tile:
    """A Sentinel-2 MGRS tile on its native UTM grid (10 980 x 10 980 px at 10 m)."""

    tile_id: str            # e.g. "T34UGV"
    epsg: int               # UTM zone EPSG code
    ulx: float              # upper-left X (m)
    uly: float              # upper-left Y (m)
    size_px: int = 10980
    res: float = 10.0

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """(xmin, ymin, xmax, ymax) in the tile CRS."""
        span = self.size_px * self.res
        return (self.ulx, self.uly - span, self.ulx + span, self.uly)

    @property
    def crs_transform(self) -> List[float]:
        """Affine transform in Earth Engine order [xScale, xShear, xT, yShear, yScale, yT]."""
        return [self.res, 0, self.ulx, 0, -self.res, self.uly]


# The tile that is missing from the HR-VPP catalogue for every year 2017-2025
# (verified against the VITO OpenSearch catalogue on 08/11/2026).
MISSING_TILE = S2Tile("T34UGV", 32634, 699960.0, 5500020.0)

# Official neighbours that overlap MISSING_TILE by ~10 km and are used for validation.
NEIGHBOURS: Dict[str, S2Tile] = {
    "T34UFV": S2Tile("T34UFV", 32634, 600000.0, 5500020.0),   # west
    "T34UGA": S2Tile("T34UGA", 32634, 699960.0, 5600040.0),   # north
    "T34UGU": S2Tile("T34UGU", 32634, 699960.0, 5400000.0),   # south
}

# Points inside the ~10 km overlap strips between the gap and each official neighbour.
# A tiny bbox around such a point returns exactly that neighbour from WEkEO HDA.
NEIGHBOUR_PROBES: Dict[str, Tuple[float, float]] = {
    "T34UFV": (23.84, 49.10),
    "T34UGA": (24.50, 49.57),
    "T34UGU": (24.50, 48.65),
}

# Approximate lon/lat footprint of the gap (used for catalogue / HDA searches).
GAP_BBOX_LONLAT = [23.71, 48.59, 25.29, 49.62]

# MONIFUN WP3 uses these four VPP parameters (season 1) for 2017-2023.
MONIFUN_PRODUCTS: List[str] = ["AMPL", "LENGTH", "RSLOPE", "TPROD"]
MONIFUN_YEARS = list(range(2017, 2024))

# WEkEO HDA dataset (UTM tiles, as in MONIFUN download_copernicus_pheno.py)
HDA_DATASET_ID = "EO:EEA:DAT:CLMS_HRVPP_VPP"
HDA_PRODUCT_GROUP = "s1"

# VITO OpenSearch catalogue (anonymous, read-only)
VITO_OPENSEARCH = "https://phenology.vgt.vito.be/products"
VITO_COLLECTION_VPP_UTM = "copernicus_r_utm-wgs84_10_m_hrvpp-vpp_p_2017-now_v01"
VITO_COLLECTION_ST_UTM = "copernicus_r_utm-wgs84_10_m_hrvpp-st_p_2017-now_v01"


@dataclass(frozen=True)
class ProductSpec:
    name: str
    unit: str
    dtype: str
    nodata: float
    default_scale: float     # physical = DN * scale (used only if no official file is available)
    description: str = ""


# Encodings follow the CLMS HR-VPP Product User Manual and the VITO catalogue band
# metadata. NOTE: for TPROD/SPROD the two sources disagree (PUM: x10 -> 0.1;
# VITO catalogue: 0.02). postprocess.resolve_scales() therefore always reads the
# scale from an official neighbouring file when one is available.
PRODUCT_SPECS: Dict[str, ProductSpec] = {
    "AMPL": ProductSpec("AMPL", "PPI", "int16", -32768, 1e-4, "Season amplitude (MAXV - MINV)"),
    "LENGTH": ProductSpec("LENGTH", "day", "int16", 0, 1.0, "Season length (EOSD - SOSD)"),
    "RSLOPE": ProductSpec("RSLOPE", "PPI/day", "int16", -32768, 1e-4, "Slope of green-down period at EOS"),
    "TPROD": ProductSpec("TPROD", "PPI*day", "uint16", 65535, 0.1, "Total productivity, SOS..EOS incl. base"),
    "SPROD": ProductSpec("SPROD", "PPI*day", "uint16", 65535, 0.1, "Seasonal productivity above base"),
    "MAXV": ProductSpec("MAXV", "PPI", "int16", -32768, 1e-4, "Maximum PPI of the season"),
    "MINV": ProductSpec("MINV", "PPI", "int16", -32768, 1e-4, "Mean of left and right minima"),
    "LSLOPE": ProductSpec("LSLOPE", "PPI/day", "int16", -32768, 1e-4, "Slope of green-up period at SOS"),
    "SOSD": ProductSpec("SOSD", "YYDOY", "int32", 0, 1.0, "Start-of-season date"),
    "EOSD": ProductSpec("EOSD", "YYDOY", "int32", 0, 1.0, "End-of-season date"),
}

# Internal integer encoding used for the Earth Engine export (keeps files small and
# loss-free at the precision of the official product). physical = DN / factor
GEE_EXPORT_FACTORS: Dict[str, float] = {
    "AMPL": 1e4, "MAXV": 1e4, "MINV": 1e4, "RSLOPE": 1e4, "LSLOPE": 1e4,
    "LENGTH": 1.0, "TPROD": 10.0, "SPROD": 10.0,
}


@dataclass
class ReconParams:
    """Parameters of the HR-VPP-like reconstruction (defaults follow the HR-VPP ATBD/PUM)."""

    cs_threshold: float = 0.60        # Cloud Score+ cs_cdf; >= threshold is clear
    dvi_soil: float = 0.09            # soil DVI (ATBD alternative constant method)
    mdvi_cap: float = 0.80            # upper limit of MDVI (ATBD)
    leaf_g: float = 0.5               # G for spherical leaf angle distribution
    bin_days: int = 10                # temporal compositing step (HR-VPP ST is 10-daily)
    max_gap_days: int = 60            # linear interpolation window
    smooth_passes: int = 2            # 1-2-1 binomial passes in time
    sos_frac: float = 0.25            # SOS threshold (25 % of amplitude, PUM)
    eos_frac: float = 0.15            # EOS threshold (15 % of amplitude, PUM)
    min_ampl: float = 0.05            # below this no season is reported
    min_length: int = 30              # days
    min_obs: int = 8                  # minimum clear observations in the window
    winter_months: tuple = (10, 11, 12, 1, 2, 3)   # base-level gap filling (ATBD)
    winter_base_pct: int = 5          # 5th percentile base level (ATBD)
    extra: dict = field(default_factory=dict)

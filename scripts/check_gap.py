#!/usr/bin/env python3
"""Step 1. Confirm that a Sentinel-2 tile is missing from the HR-VPP catalogue (anonymous, no account).

    python scripts/check_gap.py --tile T34UGV --years 2017-2025

Queries the VITO OpenSearch catalogue behind WEkEO / the CLMS viewer for every year and lists the tiles that
intersect the footprint. Writes <data>/reports/01_catalogue_availability.csv.
"""
import argparse

from _common import data_root, log, parse_years
from hrvpp_recon import catalogue, config

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", default=config.MISSING_TILE.tile_id)
ap.add_argument("--bbox", default=",".join(map(str, config.GAP_BBOX_LONLAT)), help="lon_min,lat_min,lon_max,lat_max")
ap.add_argument("--years", default="2017-2025")
ap.add_argument("--data", default=None, help="working folder (default $HRVPP_DATA or ./data)")
a = ap.parse_args()

bbox = [float(v) for v in a.bbox.split(",")]
df = catalogue.tile_availability(a.tile, bbox, parse_years(a.years))
print(df.to_string(index=False))
out = data_root(a.data) / "reports" / "01_catalogue_availability.csv"
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)
log(f"{a.tile} missing in all years: {bool((~df['available']).all())}; table: {out}")

#!/usr/bin/env python3
"""Step 2. Download official HR-VPP neighbours from WEkEO HDA (MONIFUN layout), for calibration and validation.

    export WEKEO_USER=...  WEKEO_PASSWORD=...        # free account at https://wekeo.copernicus.eu
    python scripts/download_official.py --tiles T34UFV,T34UGA,T34UGU,T35ULP --years 2017-2024

Files land in <data>/official/<PRODUCT>/<YEAR>/ exactly as with MONIFUN download_copernicus_pheno.py.
Each tile is reached with a tiny bounding box around a point inside it (config.NEIGHBOUR_PROBES, or --probe).
Truncated transfers are detected and repeated (see docs/MONIFUN_ISSUE.md). Credentials are read from the
environment (or Colab Secrets) only; nothing is written to disk.
"""
import argparse

from _common import data_root, log, parse_years
from hrvpp_recon import config, wekeo

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tiles", default="T34UFV,T34UGA,T34UGU,T35ULP")
ap.add_argument("--years", default="2017-2024")
ap.add_argument("--products", default=",".join(config.MONIFUN_PRODUCTS))
ap.add_argument("--probe", action="append", default=[], help="TILE:lon,lat for a tile not in config.NEIGHBOUR_PROBES")
ap.add_argument("--data", default=None)
a = ap.parse_args()

if not wekeo.load_credentials():
    raise SystemExit("Set WEKEO_USER and WEKEO_PASSWORD (environment variables or Colab Secrets).")
probes = dict(config.NEIGHBOUR_PROBES)
for p in a.probe:
    t, ll = p.split(":")
    probes[t] = tuple(float(v) for v in ll.split(","))
root = data_root(a.data) / "official"
for tile in a.tiles.split(","):
    if tile not in probes:
        raise SystemExit(f"no probe point for {tile}: pass --probe {tile}:lon,lat")
    for y in parse_years(a.years):
        for prod in a.products.split(","):
            if list((root / prod / str(y)).glob(f"*{tile}*_s1_{prod}.tif")):
                continue
            r = wekeo.download_one(year=y, product=prod, outdir=root, bbox=wekeo.probe_bbox(*probes[tile]))
            log(tile, y, prod, r.status)
n = wekeo.write_checksums(root, data_root(a.data) / "reports" / "02_official_checksums.csv")
log(n, "official files present under", root)

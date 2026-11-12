#!/usr/bin/env python3
"""Step 3a. Reconstruct HR-VPP season-1 parameters in Google Earth Engine and export them to Google Drive.

    earthengine authenticate                      # once
    export EE_PROJECT=<your Cloud project registered for Earth Engine>
    python scripts/ee_export.py --region T34UGV  --years 2017-2024 --folder HRVPP_export
    python scripts/ee_export.py --region T34UGVE --years 2018-2024 --folder HRVPP_export
    python scripts/ee_export.py --grid 32634,799760,5560020,7000,16980 --name MYAREA --years 2021

Regions are defined in src/hrvpp_recon/config.py (REGIONS); --grid takes EPSG,x0,y1,width,height on a 10 m grid.
One export task per year: Int16 bands AMPL, LENGTH, RSLOPE, TPROD, SOSD, EOSD, MAXV, MINV, NOBS (nodata -32768).
Big areas are split into shards by Earth Engine: merge them with scripts/mosaic_shards.py.
A whole tile costs roughly 100-150 EECU-hours per year (noncommercial quota: watch the monthly limit).
The same algorithm without Python: gee/hrvpp_recon_patch.js in the Code Editor.
"""
import argparse
import os

from _common import log, parse_years
from hrvpp_recon import config

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
g = ap.add_mutually_exclusive_group(required=True)
g.add_argument("--region", choices=sorted(config.REGIONS))
g.add_argument("--grid", help="EPSG,x0,y1,width,height (upper-left corner, 10 m pixels)")
ap.add_argument("--name", default="AREA", help="name used in file names with --grid")
ap.add_argument("--years", required=True)
ap.add_argument("--folder", default="HRVPP_export", help="Google Drive folder of the Earth Engine account")
ap.add_argument("--project", default=os.getenv("EE_PROJECT"), help="Cloud project (default $EE_PROJECT)")
ap.add_argument("--dry-run", action="store_true", help="build the tasks but do not start them")
a = ap.parse_args()

import ee  # noqa: E402

ee.Initialize(project=a.project)
from hrvpp_recon import gee_recon  # noqa: E402

if a.region:
    grid = config.REGIONS[a.region]
else:
    e, x0, y1, w, h = a.grid.split(",")
    grid = config.Grid(a.name, int(e), float(x0), float(y1), int(w), int(h))
for y in parse_years(a.years):
    t = gee_recon.export_region(y, a.folder, grid, start=not a.dry_run)
    log(y, grid.name, f"{grid.width}x{grid.height}", "task", t.id or "(not started)")
log("Follow the tasks at https://code.earthengine.google.com/tasks")

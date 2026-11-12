#!/usr/bin/env python3
"""Step 4. Learned correction on the overlap with official tiles, then MONIFUN files with the official encoding.

    python scripts/calibrate.py --raw data/raw_ee/RECON_VPP_2021_S2_T34UGVE-010m_s1.tif --year 2021 \
        --tile-id T34UGVE --neighbours T34UGA,T35ULP
    python scripts/calibrate.py --raw-dir data/raw_ee --region T34UGV --years 2017-2024

The raw reconstruction (Earth Engine or local) is compared pixel by pixel with the official tiles where they
overlap. A gradient-boosting model (features AMPL, LENGTH, RSLOPE, TPROD, SOS/EOS day, number of observations)
followed by quantile mapping is fitted per product and year and applied to the whole area. With two or more
neighbours the report holds a leave-one-neighbour-out check (train on the others, test on the held-out one).

Output: <data>/Copernicus_Phenology/<PRODUCT>/<YEAR>/VPP_<YEAR>_S2_<TILE>-010m_RECON_s1_<PRODUCT>.tif
(dtype, scale, nodata, band description, unit and PhysRange copied from an official file; tags SOURCE and
CALIBRATION mark the file as a reconstruction) and <data>/reports/calibration_<TILE>_<YEAR>.json.
"""
import argparse
import re
from pathlib import Path

from _common import data_root, log, parse_years
from hrvpp_recon import config, pipeline

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
src = ap.add_mutually_exclusive_group(required=True)
src.add_argument("--raw", help="one raw reconstruction GeoTIFF (needs --year)")
src.add_argument("--raw-dir", help="folder with *_VPP_<YEAR>_S2_<TILE>-010m_s1.tif")
ap.add_argument("--year", type=int)
ap.add_argument("--years", default=",".join(map(str, config.CALIBRATION_YEARS)))
ap.add_argument("--region", choices=sorted(config.REGIONS), help="sets --tile-id and --neighbours")
ap.add_argument("--tile-id")
ap.add_argument("--neighbours", help="official tiles, comma separated")
ap.add_argument("--method", default="gbm_qm", choices=["raw", "linear", "quantile", "gbm", "gbm_qm"])
ap.add_argument("--data", default=None)
ap.add_argument("--official", default=None, help="official files (default <data>/official)")
ap.add_argument("--out", default=None, help="default <data>/Copernicus_Phenology")
ap.add_argument("--overwrite", action="store_true")
a = ap.parse_args()

tile = a.tile_id or a.region
nbs = a.neighbours.split(",") if a.neighbours else config.REGION_NEIGHBOURS.get(a.region or "", [])
if not tile or not nbs:
    raise SystemExit("give --region, or --tile-id and --neighbours")
root = data_root(a.data)
official = Path(a.official) if a.official else root / "official"
out = Path(a.out) if a.out else root / "Copernicus_Phenology"
if a.raw:
    if not a.year:
        raise SystemExit("--raw needs --year")
    jobs = [(a.year, Path(a.raw))]
else:
    jobs = []
    for y in parse_years(a.years):
        hits = [f for f in sorted(Path(a.raw_dir).glob(f"*_VPP_{y}_S2_{tile}-010m_s1*.tif"))
                if re.search(rf"_{tile}-010m_s1(\.tif|$)", f.name)]
        if hits:
            jobs.append((y, hits[0]))
        else:
            log(y, "no raw file for", tile)
reports = []
for y, raw in jobs:
    rep = pipeline.calibrate_and_write(raw, y, official, nbs, out, tile, method=a.method, overwrite=a.overwrite,
                                       report=root / "reports" / f"calibration_{tile}_{y}.json")
    reports.append(rep)
    log(y, tile, "samples", rep["samples"], "->", len(rep["written"]), "files")
s = pipeline.summarise_loso(reports)
if len(s):
    print("\nleave-one-neighbour-out rRMSE % (lower is better)\n" + s.to_string())

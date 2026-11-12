#!/usr/bin/env python3
"""Merge Earth Engine export shards into one GeoTIFF per year (band names are kept).

    python scripts/mosaic_shards.py --in exports/ --out data/raw_ee

Earth Engine writes large exports as RECON_VPP_<Y>_S2_<NAME>-010m_s1-<row>-<col>.tif. Files without shards are
copied as they are. Output: <out>/RECON_VPP_<Y>_S2_<NAME>-010m_s1.tif
"""
import argparse
import re
import shutil
from collections import defaultdict
from pathlib import Path

from _common import log
from hrvpp_recon import pipeline

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--in", dest="src", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

groups = defaultdict(list)
for f in sorted(Path(a.src).glob("*.tif")):
    groups[re.sub(r"-\d{10}-\d{10}$", "", f.stem)].append(f)
for name, parts in groups.items():
    out = Path(a.out) / f"{name}.tif"
    if out.exists():
        log(out.name, "exists, skipped")
        continue
    out.parent.mkdir(parents=True, exist_ok=True)
    if len(parts) == 1:
        shutil.copy(parts[0], out)
    else:
        pipeline.mosaic_shards(parts, out)
    log(f"{name}: {len(parts)} part(s) -> {out}")

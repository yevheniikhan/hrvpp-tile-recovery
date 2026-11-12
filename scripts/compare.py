#!/usr/bin/env python3
"""Compare two reconstructions, or a reconstruction with an official tile, on their common pixels.

    # raw against raw (same grid), e.g. Earth Engine against the local GPU run
    python scripts/compare.py raw data/raw_ee/RECON_VPP_2021_S2_T34UGVE-010m_s1.tif \
        data/raw_local/LOCAL_VPP_2021_S2_T34UGVE-010m_s1.tif --year 2021
    # product files (MONIFUN encoding), e.g. the extension against T34UGV in their overlap, or against official
    python scripts/compare.py product A.tif B.tif --product AMPL

Prints n, Pearson r, bias (second minus first for raw, first minus second for product), RMSE and rRMSE %.
"""
import argparse

import pandas as pd

from _common import log
from hrvpp_recon import config, pipeline

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("kind", choices=["raw", "product"])
ap.add_argument("a")
ap.add_argument("b")
ap.add_argument("--year", type=int, help="needed for raw (SOS/EOS decoding)")
ap.add_argument("--product", default="AMPL")
a = ap.parse_args()

if a.kind == "raw":
    if not a.year:
        raise SystemExit("raw needs --year")
    res = pipeline.compare_raw(a.a, a.b, a.year, config.MONIFUN_PRODUCTS)
else:
    res = {a.product: pipeline.compare_products(a.a, a.b, a.product)}
df = pd.DataFrame(res).T[["n", "r", "bias", "rmse", "rrmse_pct"]]
print(df.round(4).to_string())
log("done")

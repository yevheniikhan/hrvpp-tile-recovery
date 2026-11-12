"""Reconstruction of the east extension T34UGVE without Earth Engine (Colab GPU runtime).

For each year:
1. s2local.run_region for three training strips inside T34UGV (west, north, south 10 km strips,
   which overlap the official tiles T34UFV, T34UGA, T34UGU and T35ULP) and for the extension
   itself (same grid and extent as the Earth Engine export RECON_T34UGVE_<year>).
2. Correction learned from the official tiles on the strips (gradient boosting + quantile
   mapping, as for T34UGV), with a leave-one-strip-out check.
3. MONIFUN layout in Copernicus_Phenology_gpu/<PRODUCT>/<YEAR>/..._T34UGVE-010m_RECON_s1_<PRODUCT>.tif
4. Independent checks: the extension against the official T35ULP in its southern part, and
   against the Earth Engine-based corrected T34UGV in their 10 km overlap.

Usage: python production_local_gpu.py 2021,2018,... ; log to reports/12_gpu_production.log
"""
import os
import json
import shutil
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(os.getenv("HRVPP_ROOT", "/content/drive/MyDrive/HRVPP_T34UGV_recovery"))
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd  # noqa: E402
from hrvpp_recon import learn, postprocess, s2local, validate, wekeo  # noqa: E402

LOCAL = Path("/content/local"); LOCAL.mkdir(exist_ok=True)
RAW = ROOT / "raw_gpu"; RAW.mkdir(exist_ok=True)
OUT = ROOT / "Copernicus_Phenology_gpu"
OFFD = ROOT / "official"; REP = ROOT / "reports"
PRODS4 = ["AMPL", "LENGTH", "RSLOPE", "TPROD"]
CRS = "EPSG:32634"
TILE = (699960, 5390220, 809760, 5500020)                  # T34UGV bounds
REGIONS = {                                                  # name: (x0, y1, width, height, delta bounds)
    "STRIP_W": (699960, 5500020, 1000, 10980, TILE),
    "STRIP_N": (699960, 5500020, 10980, 1000, TILE),
    "STRIP_S": (699960, 5400220, 10980, 1000, TILE),
    "T34UGVE": (799760, 5560020, 7000, 16980, (799760, 5390220, 869760, 5560020)),
}
STRIP_NB = {"STRIP_W": ["T34UFV"], "STRIP_N": ["T34UGA"], "STRIP_S": ["T34UGU", "T35ULP"]}
YEARS = [int(a) for a in sys.argv[1].split(",")]


def log(*a):
    print(time.strftime("%H:%M"), *a, flush=True)


def official(y, nb):
    d = {p: next(iter(sorted((OFFD / p / str(y)).glob(f"*{nb}*_s1_{p}.tif"))), None) for p in PRODS4}
    return {p: v for p, v in d.items() if v is not None and wekeo.is_valid_tif(v)}


for y in YEARS:
    rep = {"year": y}
    try:
        files = {}
        for name, (x0, y1, w, h, db) in REGIONS.items():
            f = LOCAL / f"GPU_VPP_{y}_S2_{name}-010m_s1.tif"
            if not f.exists():
                tmp = LOCAL / ("_" + f.name)
                s2local.run_region(y, CRS, x0, y1, w, h, str(tmp), delta_bounds=db, block_rows=512,
                                   workers=48, log=log)
                tmp.rename(f)
                shutil.copy(f, RAW / f.name)
            files[name] = f
        # training data from the strips
        dfs = []
        for name, nbs in STRIP_NB.items():
            offs = {nb: official(y, nb) for nb in nbs}
            offs = {k: v for k, v in offs.items() if v}
            d = learn.strip_samples(files[name], offs, y, n_per_strip=200_000)
            dfs.append(d)
        df = pd.concat(dfs, ignore_index=True)
        res = learn.loso(df, PRODS4, methods=["raw", "linear", "gbm_qm"], train_cap=300_000)
        res.to_csv(REP / f"12_gpu_loso_{y}.csv", index=False)
        rep["loso"] = (res.groupby(["product", "method"]).agg(r=("r", "mean"), rrmse=("rrmse_pct", "mean"),
                       std_ratio=("std_ratio", "mean")).round(4).reset_index().to_dict("records"))
        cors = {p: learn.Corrector("gbm_qm", p).fit(df) for p in PRODS4}
        tmpl = {p: official(y, "T34UFV").get(p) for p in PRODS4}
        for p in PRODS4:
            postprocess.write_product(files["T34UGVE"], p, y, OUT, "T34UGVE", learn.phys_fn_for(cors[p], y), tmpl[p],
                                      extra_tags={"CALIBRATION": f"learned from official overlap strips {y}, gbm_qm",
                                                  "INPUT": "Sentinel-2 L2A from Microsoft Planetary Computer, SCL mask, computed in Colab"})
        log(y, "T34UGVE written")
        chk = {}
        o = official(y, "T35ULP")
        if o:
            dfe = learn.strip_samples(files["T34UGVE"], {"T35ULP": o}, y, n_per_strip=200_000)
            for p in PRODS4:
                if "OFF_" + p in dfe:
                    ref = dfe["OFF_" + p].to_numpy("float64")
                    chk[f"T35ULP_{p}_raw"] = validate.metrics(ref, dfe[p].to_numpy("float64"))
                    chk[f"T35ULP_{p}_corrected"] = validate.metrics(ref, cors[p].predict(dfe))
        for p in PRODS4:
            a = OUT / p / str(y) / f"VPP_{y}_S2_T34UGVE-010m_RECON_s1_{p}.tif"
            b = ROOT / "Copernicus_Phenology" / p / str(y) / f"VPP_{y}_S2_T34UGV-010m_RECON_s1_{p}.tif"
            if b.exists():
                aa, bb = validate.overlap_arrays(a, p, b, max_pixels=500_000)
                chk[f"vs_EE_T34UGV_{p}"] = validate.metrics(bb, aa)
        rep["checks"] = chk
        rep["state"] = "done"
    except Exception as e:
        rep["state"] = "error: " + str(e)[:300]
        traceback.print_exc()
    json.dump(rep, open(REP / f"12_gpu_production_{y}.json", "w"), indent=1, default=float)
    log(y, rep["state"])
log("ALL DONE")

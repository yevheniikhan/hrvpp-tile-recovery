"""Production run (separate process): learned correction for T34UGV and the east extension T34UGVE.

Phase 1, every year: train the correction per product on all overlap strips of T34UGV (official
T34UFV, T34UGA, T34UGU, T35ULP), record the leave-one-strip-out check, write T34UGV.
Phase 2, every year: wait for the T34UGVE export, write it with the T34UGV model of the same year,
check it against T35ULP (independent) and against corrected T34UGV in the 10 km overlap.
Existing outputs are kept (restartable). Log: reports/11_production.log; results: reports/11_production_<year>.json
"""
import json, sys, time, traceback
from pathlib import Path

ROOT = Path("/content/drive/MyDrive/HRVPP_T34UGV_recovery"); sys.path.insert(0, str(ROOT / "src"))
import ee, numpy as np, pandas as pd
from hrvpp_recon import learn, postprocess, eedrive, wekeo, validate

EE_PROJ = "gen-lang-client-0200925587"
ee.Initialize(project=EE_PROJ)
LOCAL = Path("/content/recon"); OFFD = ROOT / "official"; REP = ROOT / "reports"
OUT = ROOT / "Copernicus_Phenology"
NBS = ["T34UFV", "T34UGA", "T34UGU", "T35ULP"]; PRODS4 = ["AMPL", "LENGTH", "RSLOPE", "TPROD"]
YEARS = [int(a) for a in sys.argv[1].split(",")]
CHOICE = json.loads(sys.argv[2])


def log(*a):
    print(time.strftime("%H:%M"), *a, flush=True)


def official_dict(y):
    out = {}
    for nb in NBS:
        d = {p: next(iter(sorted((OFFD / p / str(y)).glob(f"*{nb}*_s1_{p}.tif"))), None) for p in PRODS4}
        d = {p: v for p, v in d.items() if v is not None and wekeo.is_valid_tif(v)}
        if d:
            out[nb] = d
    return out


def task_state(desc):
    for t in ee.data.getTaskList()[:100]:
        if t.get("description") == desc:
            return t["state"]
    return "MISSING"


def wait_for(desc, hours=8):
    t0 = time.time()
    while time.time() - t0 < hours * 3600:
        s = task_state(desc)
        if s in ("COMPLETED", "FAILED", "CANCELLED", "MISSING"):
            return s
        time.sleep(120)
    return "TIMEOUT"


def fetch(prefix):
    hits = sorted(LOCAL.glob(prefix + "*.tif"))
    if hits:
        return hits[0]
    files = [p for p in eedrive.fetch(prefix, LOCAL, EE_PROJ) if p.suffix == ".tif"]
    return files[0] if files else None


def out_path(p, y, tile):
    return OUT / p / str(y) / f"VPP_{y}_S2_{tile}-010m_RECON_s1_{p}.tif"


def rep_load(y):
    f = REP / f"11_production_{y}.json"
    return json.load(open(f)) if f.exists() else {"year": y}


def rep_save(rep):
    json.dump(rep, open(REP / f"11_production_{rep['year']}.json", "w"), indent=1, default=float)


def models(y, rep):
    offs = official_dict(y)
    rep["official_strips"] = {k: sorted(v) for k, v in offs.items()}
    rec = fetch(f"RECON_VPP_{y}_S2_T34UGV-")
    if rec is None:
        rep["export_T34UGV"] = wait_for(f"RECON_T34UGV_{y}")
        rec = fetch(f"RECON_VPP_{y}_S2_T34UGV-")
    df = learn.strip_samples(rec, offs, y, n_per_strip=200_000)
    cors = {p: learn.Corrector(CHOICE[p], p).fit(df) for p in PRODS4}
    tmpl = {p: next(iter(v[p] for v in offs.values() if p in v)) for p in PRODS4}
    return rec, df, cors, tmpl, offs


CACHE = {}
# ---------------------------------------------------------------- phase 1: T34UGV
for y in YEARS:
    rep = rep_load(y)
    try:
        rec, df, cors, tmpl, offs = models(y, rep)
        CACHE[y] = (cors, tmpl, offs)
        if "loso" not in rep:
            res = learn.loso(df, PRODS4, methods=["raw", "linear", "gbm_qm"], train_cap=300_000)
            res.to_csv(REP / f"11_loso_{y}.csv", index=False)
            rep["choice"] = CHOICE
            rep["loso"] = (res.groupby(["product", "method"])
                           .agg(r=("r", "mean"), abs_bias=("bias", lambda v: v.abs().mean()),
                                rrmse=("rrmse_pct", "mean"), std_ratio=("std_ratio", "mean"))
                           .round(4).reset_index().to_dict("records"))
        for p in PRODS4:
            if not out_path(p, y, "T34UGV").exists():
                postprocess.write_product(rec, p, y, OUT, "T34UGV", learn.phys_fn_for(cors[p], y), tmpl[p],
                    extra_tags={"CALIBRATION": f"learned from official overlap strips {y} ({', '.join(offs)}), method {CHOICE[p]}"})
        rep["T34UGV"] = "done"
        log(y, "T34UGV done")
    except Exception as e:
        rep["T34UGV"] = "error: " + str(e)[:300]
        traceback.print_exc()
    rep_save(rep)

# ---------------------------------------------------------------- phase 2: T34UGVE
for y in YEARS:
    rep = rep_load(y)
    try:
        cors, tmpl, offs = CACHE[y] if y in CACHE else models(y, rep)[2:]
        rep["export_T34UGVE"] = wait_for(f"RECON_T34UGVE_{y}")
        recE = fetch(f"RECON_VPP_{y}_S2_T34UGVE-")
        for p in PRODS4:
            if not out_path(p, y, "T34UGVE").exists():
                postprocess.write_product(recE, p, y, OUT, "T34UGVE", learn.phys_fn_for(cors[p], y), tmpl[p],
                    extra_tags={"CALIBRATION": f"model trained on T34UGV overlap strips {y}, method {CHOICE[p]}"})
        chk = {}
        if "T35ULP" in offs:
            dfe = learn.strip_samples(recE, {"T35ULP": offs["T35ULP"]}, y, n_per_strip=200_000)
            for p in PRODS4:
                if "OFF_" + p in dfe:
                    ref = dfe["OFF_" + p].to_numpy("float64")
                    chk[f"T35ULP_{p}_raw"] = validate.metrics(ref, dfe[p].to_numpy("float64"))
                    chk[f"T35ULP_{p}_corrected"] = validate.metrics(ref, cors[p].predict(dfe))
        for p in PRODS4:
            a, b = validate.overlap_arrays(out_path(p, y, "T34UGVE"), p, out_path(p, y, "T34UGV"), max_pixels=500_000)
            chk[f"UGVE_vs_UGV_{p}"] = validate.metrics(b, a)
        rep["checks_T34UGVE"] = chk
        rep["T34UGVE"] = "done"
        log(y, "T34UGVE done")
    except Exception as e:
        rep["T34UGVE"] = "error: " + str(e)[:300]
        traceback.print_exc()
    rep_save(rep)
log("ALL DONE")

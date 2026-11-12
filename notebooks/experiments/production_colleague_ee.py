"""East extension T34UGVE from the Earth Engine exports computed by a colleague (2018-2024).

The colleague ran gee/hrvpp_recon_patch.js in the Code Editor under his own Earth Engine quota. Earth
Engine split each year into row shards (RECON_VPP_<y>_S2_T34UGVE-010m_s1-<row>-<col>.tif).

For each year:
1. Mosaic of the shards into one raw file on the T34UGVE grid (band names kept).
2. Learned correction (gradient boosting + quantile mapping, as everywhere else), trained on the
   extension's own overlap with the official neighbours T34UGA (10 km strip in the north-west) and
   T35ULP (east and south). The T34UGV raw exports used for the 2017 extension are no longer on
   Drive, so here the extension calibrates itself on its official neighbours.
   Check: leave-one-neighbour-out (train on T34UGA, test on T35ULP and the other way round).
3. MONIFUN layout in Copernicus_Phenology/<PRODUCT>/<YEAR>/VPP_<y>_S2_T34UGVE-010m_RECON_s1_<PRODUCT>.tif
4. Further checks: against the corrected T34UGV in the 10 km overlap, and against the GPU
   (Planetary Computer) version where it exists, raw against raw and corrected against corrected.

Usage: python production_colleague_ee.py 2018,2019,... [shard folder, default /content/colleague]
Report: reports/13_colleague_<year>.json, reports/13_colleague_loso_<year>.csv
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
import numpy as np  # noqa: E402
import rasterio  # noqa: E402
from rasterio.merge import merge  # noqa: E402
from rasterio.windows import Window  # noqa: E402
from hrvpp_recon import learn, postprocess, validate, wekeo  # noqa: E402

YEARS = [int(a) for a in sys.argv[1].split(",")]
SHARDS = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("/content/colleague")
LOCAL = Path("/content/col_local"); LOCAL.mkdir(exist_ok=True)
OUT = ROOT / "Copernicus_Phenology"; OUTG = ROOT / "Copernicus_Phenology_gpu"; RAWG = ROOT / "raw_gpu"
OFFD = ROOT / "official"; REP = ROOT / "reports"
PRODS4 = ["AMPL", "LENGTH", "RSLOPE", "TPROD"]
NBS = ["T34UGA", "T35ULP"]
GRID = (7000, 16980, 799760.0, 5560020.0)          # width, height, x0, y1 (EPSG:32634)


def log(*a):
    print(time.strftime("%H:%M"), *a, flush=True)


def official(y, nb):
    d = {p: next(iter(sorted((OFFD / p / str(y)).glob(f"*{nb}*_s1_{p}.tif"))), None) for p in PRODS4}
    return {p: v for p, v in d.items() if v is not None and wekeo.is_valid_tif(v)}


def out_path(root, p, y, tile):
    return root / p / str(y) / f"VPP_{y}_S2_{tile}-010m_RECON_s1_{p}.tif"


def mosaic(y):
    out = LOCAL / f"RECON_VPP_{y}_S2_T34UGVE-010m_s1.tif"
    parts = sorted(SHARDS.glob(f"RECON_VPP_{y}_S2_T34UGVE-*.tif"))
    if out.exists():
        return out, parts
    if not parts:
        raise FileNotFoundError(f"no shards for {y} in {SHARDS}")
    with rasterio.open(parts[0]) as s0:
        prof, desc = s0.profile.copy(), s0.descriptions
    nod = prof.get("nodata")
    nod = postprocess.GEE_NODATA if nod is None else nod
    arr, tr = merge([str(p) for p in parts], nodata=nod)
    for k in ("photometric", "interleave", "blockxsize", "blockysize", "tiled", "compress"):
        prof.pop(k, None)
    prof.update(driver="GTiff", height=arr.shape[1], width=arr.shape[2], transform=tr, nodata=nod,
                tiled=True, blockxsize=512, blockysize=512, compress="deflate", BIGTIFF="IF_SAFER")
    tmp = out.with_name("_" + out.name)
    with rasterio.open(tmp, "w", **prof) as dst:
        dst.write(arr)
        for i, d in enumerate(desc, 1):
            if d:
                dst.set_band_description(i, d)
    del arr
    tmp.rename(out)
    return out, parts


def raw_compare(a_tif, b_tif, y, step=3400, rows=256):
    """Raw features of two reconstructions on the same grid, sampled in horizontal bands."""
    acc = {p: ([], []) for p in PRODS4}
    with rasterio.open(a_tif) as A, rasterio.open(b_tif) as B:
        if (A.width, A.height, A.transform) != (B.width, B.height, B.transform):
            return {"error": "grids differ"}
        for r0 in range(0, A.height, step):
            w = Window(0, r0, A.width, min(rows, A.height - r0))
            fa, fb = learn.read_features(A, y, w), learn.read_features(B, y, w)
            for p in PRODS4:
                acc[p][0].append(fa[p].ravel()); acc[p][1].append(fb[p].ravel())
    return {p: validate.metrics(np.concatenate(acc[p][0]), np.concatenate(acc[p][1])) for p in PRODS4}


for y in YEARS:
    rep = {"year": y, "source": "Earth Engine export by a colleague, gee/hrvpp_recon_patch.js"}
    try:
        rec, parts = mosaic(y)
        rep["shards"] = [p.name for p in parts]
        with rasterio.open(rec) as s:
            g = (s.width, s.height, round(s.transform.c, 2), round(s.transform.f, 2))
            rep["grid"] = list(g) + [str(s.crs)]
            rep["grid_ok"] = g == GRID and str(s.crs) == "EPSG:32634"
        log(y, "mosaic", rep["grid"], "ok" if rep["grid_ok"] else "GRID DIFFERS")
        offs = {nb: official(y, nb) for nb in NBS}
        offs = {k: v for k, v in offs.items() if v}
        df = learn.strip_samples(rec, offs, y, n_per_strip=300_000)
        rep["samples"] = {k: int(v) for k, v in df["strip"].value_counts().items()}
        res = learn.loso(df, PRODS4, methods=["raw", "linear", "gbm_qm"], train_cap=300_000)
        res.to_csv(REP / f"13_colleague_loso_{y}.csv", index=False)
        rep["loso"] = res[["product", "test_strip", "method", "n", "r", "bias", "rrmse_pct", "std_ratio"]] \
            .round(4).to_dict("records")
        log(y, "loso done")
        cors = {p: learn.Corrector("gbm_qm", p).fit(df) for p in PRODS4}
        tmpl = {p: official(y, "T34UFV").get(p) or next(iter(v[p] for v in offs.values() if p in v)) for p in PRODS4}
        for p in PRODS4:
            if not out_path(OUT, p, y, "T34UGVE").exists():
                postprocess.write_product(
                    rec, p, y, OUT, "T34UGVE", learn.phys_fn_for(cors[p], y), tmpl[p],
                    extra_tags={"CALIBRATION": f"learned from the official neighbours {', '.join(offs)} "
                                               f"on the extension itself, {y}, gbm_qm",
                                "INPUT": "Earth Engine, Sentinel-2 SR harmonised, Cloud Score+ and SCL mask "
                                         "(gee/hrvpp_recon_patch.js, run by a colleague)"})
        log(y, "T34UGVE written")
        chk = {}
        for p in PRODS4:
            e = out_path(OUT, p, y, "T34UGVE")
            u = out_path(OUT, p, y, "T34UGV")
            if u.exists():
                a, b = validate.overlap_arrays(e, p, u, max_pixels=500_000)
                chk[f"vs_T34UGV_{p}"] = validate.metrics(b, a)
            gp = out_path(OUTG, p, y, "T34UGVE")
            if gp.exists():
                a, b = validate.overlap_arrays(e, p, gp, max_pixels=1_000_000)
                chk[f"vs_GPU_corrected_{p}"] = validate.metrics(b, a)
        graw = RAWG / f"GPU_VPP_{y}_S2_T34UGVE-010m_s1.tif"
        if graw.exists():
            for p, m in raw_compare(graw, rec, y).items():
                chk[f"vs_GPU_raw_{p}"] = m
        rep["checks"] = chk
        rep["state"] = "done"
    except Exception as e:
        rep["state"] = "error: " + str(e)[:300]
        traceback.print_exc()
    json.dump(rep, open(REP / f"13_colleague_{y}.json", "w"), indent=1, default=float)
    log(y, rep["state"])
log("ALL DONE")

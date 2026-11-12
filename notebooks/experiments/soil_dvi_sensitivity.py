#@title Iteration 2: soil-DVI sensitivity (2021) on 4x4 km blocks in the three overlap strips (batch exports, no timeouts)
# Run in the Colab notebook after the setup cells (needs ROOT, Earth Engine auth and the official 2021 files).
import os
import time, numpy as np, pandas as pd, rasterio, ee
from pathlib import Path
from hrvpp_recon import gee_recon, validate, eedrive
from hrvpp_recon.config import ReconParams, GEE_EXPORT_FACTORS

Y = 2021  # also run with 2020
EE_PROJ = os.environ["EE_PROJECT"]  # your Cloud project registered for Earth Engine
OFFICIAL = Path("/content/official"); REP = Path(ROOT) / "reports"; SENS = Path("/content/sens")
PRODS = ["AMPL", "LENGTH", "RSLOPE", "TPROD"]
BLOCKS = {"T34UFV": (702000, 5452000), "T34UGA": (740000, 5497000), "T34UGU": (740000, 5397000)}  # UL corner, 400x400 px
VARIANTS = [0.05, 0.07, 0.09, 0.11]

tasks = {}
for dvis in VARIANTS:
    m = gee_recon.reconstruct_year(Y, p=ReconParams(dvi_soil=dvis))["metrics"]
    bands = [m.select(b).multiply(GEE_EXPORT_FACTORS[b]).round() for b in PRODS + ["MINV", "MAXV"]]
    img = ee.Image.cat(bands).unmask(-32768).toInt16()
    for nb, (x0, y1) in BLOCKS.items():
        name = f"SENS_{Y}_DVIS{int(round(dvis * 100)):02d}_{nb}"
        t = ee.batch.Export.image.toDrive(image=img, description=name, folder="HRVPP_T34UGV_sens",
                                          fileNamePrefix=name, crs="EPSG:32634",
                                          crsTransform=[10, 0, x0, 0, -10, y1], dimensions="400x400", maxPixels=1e9)
        t.start(); tasks[name] = t
print(len(tasks), "exports started", flush=True)

t0 = time.time()
while time.time() - t0 < 3 * 3600:
    st = {n: t.status()["state"] for n, t in tasks.items()}
    if all(s in ("COMPLETED", "FAILED", "CANCELLED") for s in st.values()):
        break
    time.sleep(60)
print({s: list(st.values()).count(s) for s in set(st.values())}, f"{(time.time() - t0) / 60:.0f} min", flush=True)

files = eedrive.fetch(f"SENS_{Y}_DVIS", SENS, EE_PROJ)
rows = []
for f in sorted(p for p in files if p.suffix == ".tif"):
    dvis, nb = int(f.stem.split("_")[2][4:]) / 100, f.stem.split("_")[3].split("-")[0]
    with rasterio.open(f) as s:
        minv = s.read(5).astype(float)
    for prod in PRODS:
        off = sorted((OFFICIAL / prod / str(Y)).glob(f"*{nb}*_s1_{prod}.tif"))[0]
        a, b = validate.overlap_arrays(f, prod, off, max_pixels=None)
        mm = validate.metrics(a, b)
        mm.update(dvi_soil=dvis, neighbour=nb, product=prod,
                  minv_zero_pct=float(100 * np.mean(minv[minv != -32768] <= 0)))
        rows.append(mm)
df = pd.DataFrame(rows)
df.to_csv(REP / f"08_sensitivity_dvi_soil_{Y}.csv", index=False)
summary = df.groupby(["product", "dvi_soil"])[["n", "r", "bias", "rrmse_pct", "ols_slope", "minv_zero_pct"]].agg(
    {"n": "sum", "r": "mean", "bias": "mean", "rrmse_pct": "mean", "ols_slope": "mean", "minv_zero_pct": "mean"}).round(3)
print(summary.to_csv(sep="|"))

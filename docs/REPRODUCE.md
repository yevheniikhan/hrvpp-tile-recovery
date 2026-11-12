# Reproducing the reconstruction

Everything needed is public: Sentinel-2 (Earth Engine or Microsoft Planetary Computer), the official HR-VPP
tiles (WEkEO, free account) and this code. No file in the repository contains account names, project IDs or
passwords; you provide your own through environment variables (see `.env.example`).

| variable | used by | how to get it |
|---|---|---|
| `EE_PROJECT` | `scripts/ee_export.py`, notebook | a Google Cloud project registered for Earth Engine (noncommercial use is free, with a monthly compute quota) |
| `WEKEO_USER`, `WEKEO_PASSWORD` | `scripts/download_official.py`, notebook | free account at <https://wekeo.copernicus.eu>; accept the CLMS dataset terms once on the website |
| `HRVPP_DATA` | all scripts | working folder, default `./data` |

In Colab put the WEkEO credentials in **Secrets** (key icon), not in a cell.

## 0. Install and test

```bash
git clone https://github.com/yevheniikhan/hrvpp-tile-recovery && cd hrvpp-tile-recovery
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt                 # core + Earth Engine + WEkEO clients
pip install -r requirements-local.txt           # only for the route without Earth Engine
pytest -q                                       # synthetic data, no network, about 5 s
```

## 1. Confirm the gap (no account, < 1 min)

```bash
python scripts/check_gap.py --tile T34UGV --years 2017-2025
```

`available = False` in every year means the tile is a production gap: re-downloading cannot fill it.

## 2. Official neighbours (WEkEO, 10–40 min)

```bash
python scripts/download_official.py --tiles T34UFV,T34UGA,T34UGU,T35ULP --years 2017-2024
```

About 0.25 GB per file; 4 tiles × 4 products × 8 years ≈ 32 GB. Files land in `data/official/<PRODUCT>/<YEAR>/`
exactly as with the MONIFUN download script; truncated transfers are detected and repeated.

## 3. Raw reconstruction: choose one route

### 3a. Earth Engine from Python

```bash
earthengine authenticate
export EE_PROJECT=<your project>
python scripts/ee_export.py --region T34UGV  --years 2017-2024 --folder HRVPP_export
python scripts/ee_export.py --region T34UGVE --years 2017-2024 --folder HRVPP_export
```

* One task per year, exported to the Google Drive of the Earth Engine account. Follow them in the
  [task list](https://code.earthengine.google.com/tasks).
* Cost: about 130 EECU-hours per year for the whole T34UGV tile (the extension is of similar size). The
  noncommercial tier has a monthly limit; when it is used up, tasks wait until the 1st of the next month.
  Plan large runs accordingly, or split the years between several accounts.
* Large areas are written in row shards. Download the folder and merge:
  `python scripts/mosaic_shards.py --in <downloaded folder> --out data/raw`.

### 3b. Earth Engine Code Editor (no Python)

Open <https://code.earthengine.google.com>, choose your Cloud project, paste
[`gee/hrvpp_recon_patch.js`](../gee/hrvpp_recon_patch.js), check the settings block (years, area in EPSG:32634,
Drive folder) and press **Run**, then **Run** for every task in the Tasks tab. `MODE = 'test'` prints values at six
pixels that can be compared with the Python version (they agree exactly). Share the Drive folder with whoever does
step 4; do not rename the files. This is how the extension for 2018–2024 was produced.

### 3c. Without Earth Engine (Planetary Computer, CPU or GPU)

```bash
python scripts/local_reconstruct.py --region T34UGVE --years 2021 --out data/raw
```

* Sentinel-2 L2A is read directly from Microsoft Planetary Computer (no account). Same algorithm and output bands;
  the cloud mask is SCL-based because Cloud Score+ exists only in Earth Engine.
* On a Colab L4 GPU the extension takes about 1.5 h per year, mostly reading 400–650 scenes per block. If the GPU
  runs out of memory (for example with two processes on one GPU), lower `--gpu-chunk`.
* Colab stops idle runtimes: keep the notebook tab open for long runs, or use a machine of your own.

## 4. Correction and MONIFUN files (≈ 10 min per area and year on a CPU)

```bash
python scripts/calibrate.py --raw-dir data/raw --region T34UGV  --years 2017-2024
python scripts/calibrate.py --raw-dir data/raw --region T34UGVE --years 2018-2024
```

* For each year the raw reconstruction is sampled where it overlaps the official tiles (up to 300 000 pixels per
  tile). A small gradient-boosting model (40 trees, 15 leaves) per product, followed by quantile mapping, is
  fitted and applied to the whole area.
* With two or more neighbours the report holds a leave-one-neighbour-out check of `raw`, `linear` and the chosen
  method (`data/reports/calibration_<TILE>_<YEAR>.json`); the script prints a summary table.
* Products are written to `data/Copernicus_Phenology/` with the official encoding and band metadata. Existing files
  are kept unless you pass `--overwrite`.

## 5. Checks

```bash
# two raw runs on the same grid (e.g. Earth Engine against the local route)
python scripts/compare.py raw data/raw/RECON_VPP_2021_S2_T34UGVE-010m_s1.tif data/raw/LOCAL_VPP_2021_S2_T34UGVE-010m_s1.tif --year 2021
# the extension against T34UGV in their 10 km overlap, or a product against an official tile
python scripts/compare.py product data/Copernicus_Phenology/AMPL/2021/VPP_2021_S2_T34UGVE-010m_RECON_s1_AMPL.tif \
    data/Copernicus_Phenology/AMPL/2021/VPP_2021_S2_T34UGV-010m_RECON_s1_AMPL.tif --product AMPL
```

Expected values for this area are in [`VALIDATION.md`](VALIDATION.md). Large deviations usually mean a wrong
year in `--year`, files of different areas, or a raw file that is not complete (check the shard merge).

## Another area

Any rectangle on a 10 m UTM grid works:

1. add it to `REGIONS` and `REGION_NEIGHBOURS` in `src/hrvpp_recon/config.py`, or pass
   `--grid EPSG,x0,y1,width,height --name NAME` to `ee_export.py` / `local_reconstruct.py`;
2. download the official tiles that overlap it (step 2, `--probe TILE:lon,lat` for tiles not in the config);
3. calibrate with `--tile-id NAME --neighbours TILE1,TILE2`.

The correction needs an overlap with at least one official tile; two or more give an honest held-out check.

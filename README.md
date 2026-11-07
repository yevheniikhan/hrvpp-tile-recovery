# hrvpp-tile-recovery

Detect, reconstruct and validate a missing tile of the Copernicus **HR-VPP** (High Resolution Vegetation Phenology
and Productivity) product for the **MONIFUN WP3** auxiliary-data workflow.

The concrete case is Sentinel-2 tile **T34UGV** in western Ukraine (≈ 23.7–25.3° E, 48.6–49.6° N; Ivano-Frankivsk,
Kalush, Stryi). It is missing from the CLMS catalogue for every year 2017–2025, so the MONIFUN download script
(`download_copernicus_pheno.py`, WEkEO HDA) leaves a hole in the phenology layers `AMPL`, `LENGTH`, `RSLOPE`,
`TPROD` for all NFI plots in that area.

```
 official T34UFV │  T34UGA (official, north)
  (west)     ┌───┼──────────────────────┐
             │ ▒▒│  T34UGV  (missing)    │  ▒ = ~10 km overlap strips used for validation
             │ ▒▒│  reconstructed from   │
             │ ▒▒│  Sentinel-2 in GEE    │
             └───┼──────────────────────┘
                 │  T34UGU (official, south)
```

## What it does

1. **Confirms the gap** in the anonymous VITO OpenSearch catalogue (no credentials).
2. **Downloads** official neighbours via WEkEO HDA with the MONIFUN logic (same folder layout and summary CSV).
3. **Reconstructs** season-1 VPP parameters in Google Earth Engine: Sentinel-2 → Cloud Score+ → PPI with the
   HR-VPP ATBD formula → 10-day composites → interpolation and smoothing → SOS/EOS at 25 %/15 % of amplitude.
4. **Validates** pixel by pixel in the ~10 km overlap strips with T34UFV, T34UGA and T34UGU, and cross-checks the
   Earth Engine code against an independent NumPy implementation.
5. **Writes MONIFUN-ready files** with the official encoding:
   `Copernicus_Phenology/<PRODUCT>/<YEAR>/VPP_<YEAR>_S2_T34UGV-010m_RECON_s1_<PRODUCT>.tif`.

Details: [`docs/METHOD.md`](docs/METHOD.md) · [`docs/VALIDATION.md`](docs/VALIDATION.md) ·
[`docs/RUNBOOK.md`](docs/RUNBOOK.md)

## Quick start (Google Colab)

1. Copy the repository folder to Google Drive as `MyDrive/HRVPP_T34UGV_recovery` (the notebook imports
   `src/` from there).
2. Open `notebooks/hrvpp_tile_recovery.ipynb` in Colab.
3. In Colab Secrets (🔑) add `WEKEO_USER` and `WEKEO_PASSWORD` (free account at <https://wekeo.copernicus.eu>).
4. Set `EE_PROJECT` to a Google Cloud project with the Earth Engine API enabled.
5. Run all cells. Tick `START_FULL_EXPORTS` to produce the full 10 m tile for 2017–2023.

## Repository layout

```
src/hrvpp_recon/
  config.py        tile geometry, encodings, algorithm parameters
  catalogue.py     anonymous availability check (VITO OpenSearch)
  wekeo.py         WEkEO HDA download (MONIFUN-compatible)
  gee_recon.py     Earth Engine reconstruction
  phenology_np.py  NumPy reference implementation of the season metrics
  fastval.py       point-based validation and EE-vs-NumPy cross-check
  postprocess.py   export → MONIFUN layout with official encoding
  validate.py      overlap metrics
notebooks/         Colab notebook (main entry point)
tests/             unit tests (pytest), run without Earth Engine or credentials
docs/              method, validation protocol, runbook
```

## Tests

```bash
pip install -r requirements.txt
pytest -q
```

## Status and caveats

* Reconstructed rasters are **not official CLMS products**. They are flagged in file names (`_RECON_`) and GeoTIFF
  tags, and should be reported as gap-filled.
* The method approximates, but does not replicate, the operational TIMESAT processing. See the limitations in
  `docs/METHOD.md` and the validation results before using the data.

## Licence and data

Code: MIT. Copernicus data: free and open under Regulation (EU) No 377/2014 and Commission Delegated Regulation (EU)
No 1159/2013; attribute "Copernicus Land Monitoring Service, European Environment Agency (EEA)" and state that the
data were modified.

---

### Коротко українською

Інструмент для тайла HR-VPP **T34UGV** (Івано-Франківщина / Львівщина), якого немає в каталозі CLMS за жоден рік
2017–2025. Ноутбук підтверджує відсутність тайла, відтворює параметри AMPL, LENGTH, RSLOPE, TPROD із Sentinel-2
в Google Earth Engine за методикою HR-VPP, перевіряє точність на смугах перекриття з офіційними сусідніми тайлами і
записує файли в структуру папок MONIFUN.

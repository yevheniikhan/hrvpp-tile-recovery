# Runbook

## One-off setup (≈ 5 min)

| what | where |
|---|---|
| WEkEO account | <https://wekeo.copernicus.eu> → register (free) |
| Colab Secrets | Colab left bar → 🔑 → add `WEKEO_USER`, `WEKEO_PASSWORD`, enable "Notebook access" |
| Earth Engine project | <https://code.earthengine.google.com> → note the Cloud project ID, or create one registered for Earth Engine |
| Drive folder | `MyDrive/HRVPP_T34UGV_recovery` containing `src/`, `notebooks/`, `docs/` |

## Run order

| step | run time | output |
|---|---|---|
| 0 Parameters, 0.1 Environment | 1–2 min | Drive mounted, libraries installed |
| 1 Catalogue check | < 1 min | `reports/01_catalogue_availability.csv` |
| 2a HDA search, 2b neighbour download | 10–40 min (≈ 0.25 GB per file, 3 tiles × 4 products × validation years) | `Copernicus_Phenology_official/…`, checksums |
| 3 Reconstruction + quick look | 2–5 min | thumbnails, scene count |
| 4 EE vs NumPy | 1–3 min | `reports/04_crosscheck_ee_numpy.csv` |
| 5 Point validation | 5–15 min | `reports/05_validation_points.csv`, scatter plots |
| 6 Full exports | 1–3 h per year on Earth Engine (run in background) | `RECON_VPP_<YEAR>_…tif` in the Drive folder |
| 6.1 Convert to MONIFUN layout | ≈ 2 min per year | `Copernicus_Phenology/<PRODUCT>/<YEAR>/…_RECON_…tif` |
| 6.2 Full-overlap validation | ≈ 1 min per year | `reports/06_validation_full_overlap.csv` |
| 7 Reproducibility + report | 1 min | `reports/REPORT.md`, `run_log.json` |

## Using the result in MONIFUN

Copy (or point the R scripts to) `Copernicus_Phenology/` next to the official downloads. The reconstructed files use
the same folder structure, CRS, grid, dtype, scale and nodata as the official tiles, so
`_run_extract_xbyx_linux.R` / `_windows.R` pick them up without changes. Keep the `_RECON_` files identifiable in the
plot table (e.g. a flag column) so that WP3 can test model sensitivity with and without gap-filled covariates.

## Troubleshooting

* **`ee.Initialize` fails**: set `EE_PROJECT` to a project registered for Earth Engine.
* **HDA returns 0 results for a neighbour**: check the WEkEO account and terms acceptance for CLMS datasets on the
  WEkEO website; the probe points are in `config.NEIGHBOUR_PROBES`.
* **Export "User memory limit exceeded"**: run the export at 20 m (`export_year(..., scale=20)`) or reduce
  `ReconParams.smooth_passes`; 10 m is recommended for MONIFUN 5×5 windows.

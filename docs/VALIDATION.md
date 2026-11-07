# Validation protocol

## Design

T34UGV overlaps three official HR-VPP tiles by ~10 km. In those strips official and reconstructed pixels lie on the
same 10 m Sentinel-2 grid, so they can be compared directly without resampling.

| neighbour | position | overlap strip (EPSG:32634) |
|---|---|---|
| T34UFV | west | x 699 960–709 800 |
| T34UGA | north | y 5 490 240–5 500 020 |
| T34UGU | south | y 5 390 220–5 400 000 |

Landscape in the strips (Carpathian foothills, Dniester valley, mixed forest and arable land) is representative
of the gap itself.

## Levels

1. **Implementation check.** The Earth Engine array code is compared with the independent NumPy implementation
   on the same smoothed series (300 pixels). Expected: differences at float32 rounding level.
2. **Unit tests.** The NumPy implementation is tested on synthetic double-logistic curves with analytically known
   SOS, EOS, amplitude and integral (`tests/test_phenology_np.py`), and the raster encoding and overlap logic are
   tested on synthetic GeoTIFFs (`tests/test_postprocess_validate.py`).
3. **Accuracy check (point-based).** `N_VALIDATION_POINTS` random pixel centres per strip and year (default 3 000
   × 3 strips) for AMPL, LENGTH, RSLOPE, TPROD.
4. **Accuracy check (full overlap).** After the full export, all pixels of the strips.
5. **Reproducibility.** Repeated Earth Engine sampling with a fixed seed must be identical. SHA-256 checksums of
   the official downloads are stored so a re-download can be compared file by file.

## Metrics

n, Pearson r, bias (reconstructed − official), MAE, RMSE, relative RMSE (% of the official mean), OLS slope and
intercept, and the share of pixels with a valid season in each source.

## Acceptance guidance (proposed, to be agreed with WP3)

The reconstructed tile is a covariate source for a model, not a reference dataset. As a pragmatic threshold, its
error should be clearly smaller than the between-pixel variability of the official product in the same landscape:

| product | proposed acceptance |
|---|---|
| AMPL | r ≥ 0.80, abs(bias) ≤ 0.10 PPI |
| TPROD | r ≥ 0.80, rRMSE ≤ 20 % |
| LENGTH | r ≥ 0.60, abs(bias) ≤ 15 days |
| RSLOPE | r ≥ 0.50 (the noisiest parameter) |

If a product fails, the recommended action for MONIFUN is to keep the product as missing for T34UGV plots
(or flag the plots) rather than use the reconstruction.

## Results

Results are written by the notebook to `reports/` (`05_validation_points.csv`, `06_validation_full_overlap.csv`,
`REPORT.md`, scatter plots). See the run report in the project Drive folder.

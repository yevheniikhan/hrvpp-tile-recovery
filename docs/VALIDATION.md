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

## Results for 2021 (run of 08/11/2026)

Earth Engine noncommercial project; 1 828 Sentinel-2 scenes in the 2021 window; 99.7 % of T34UGV
received a valid season.

**Implementation check.** On 291 random pixels, the Earth Engine metrics and the NumPy reference agree exactly
(maximum absolute difference 0 for AMPL, LENGTH, RSLOPE, MAXV, MINV, SOS, EOS; 5.7e-14 for TPROD).

**Accuracy, full overlap** (exported 10 m tile vs official V105 tiles; 1 000 000 random pixels per strip; the
pixel-centre comparison with 600 points per strip gave the same picture):

| strip | product | r | bias | RMSE | rRMSE % | valid in both % |
|---|---|---|---|---|---|---|
| T34UFV | AMPL | 0.832 | +0.165 | 0.426 | 23.8 | 99.0 |
| T34UFV | LENGTH | 0.706 | +7.5 d | 27.5 d | 19.2 | 99.0 |
| T34UFV | RSLOPE | 0.692 | −0.010 | 0.016 | 63.2 | 99.0 |
| T34UFV | TPROD | 0.919 | −1.9 | 26.2 | 15.6 | 99.0 |
| T34UGA | AMPL | 0.813 | +0.294 | 0.486 | 24.7 | 96.3 |
| T34UGA | LENGTH | 0.743 | +3.9 d | 26.2 d | 19.9 | 96.3 |
| T34UGA | RSLOPE | 0.712 | −0.011 | 0.019 | 60.5 | 96.3 |
| T34UGA | TPROD | 0.935 | +7.5 | 23.5 | 14.1 | 96.3 |
| T34UGU | AMPL | 0.827 | +0.164 | 0.424 | 24.3 | 99.5 |
| T34UGU | LENGTH | 0.682 | +7.3 d | 27.6 d | 18.6 | 99.5 |
| T34UGU | RSLOPE | 0.587 | −0.008 | 0.014 | 66.1 | 99.5 |
| T34UGU | TPROD | 0.921 | −2.2 | 25.4 | 15.0 | 99.5 |

**Leave-one-strip-out linear calibration** (official = k · recon + c fitted on two strips, tested on the third; mean
over the three tests):

| product | r | bias raw → cal | rRMSE raw → cal % | decision |
|---|---|---|---|---|
| AMPL | 0.824 | +0.208 → −0.003 | 24.3 → 21.6 | **calibrate** (all-strip fit k = 0.918, c = −0.040) |
| RSLOPE | 0.664 | −0.0094 → −0.0002 | 63.3 → 52.1 | **calibrate** (k = 1.164, c = 0.0067) |
| LENGTH | 0.710 | +6.2 d → +0.5 d | 19.3 → 17.2 | keep uncalibrated: passes thresholds; k = 0.69 would compress natural variance |
| TPROD | 0.925 | +1.1 → 0.0 | 14.9 → 15.2 | keep uncalibrated: no gain |

**Against the acceptance thresholds:** TPROD and LENGTH pass uncalibrated; AMPL passes after calibration (bias within
±0.10 in every held-out strip); RSLOPE passes the r threshold and is unbiased after calibration but remains the
noisiest parameter (rRMSE ≈ 50 %), as expected from the different slope definition (10-day differences vs TIMESAT
derivative). For every year the calibration is re-fitted on that year's own overlap strips
(`validate.strip_calibration`) and written to the GeoTIFF tag `CALIBRATION`.

**Cause of the AMPL bias.** The first guess was a clipped winter minimum (MINV = 0 where winter DVI falls below the
constant soil DVI). The soil-DVI sensitivity test below rejects it: the bias grows as fewer pixels are clipped.
The bias therefore comes from the peak side (MAXV), which has not been checked yet because the official MAXV
files were not downloaded.

**Reproducibility.** Downloads are integrity-checked (two truncated WEkEO transfers were caught and repeated during
this run); the export file is read back and validated before use.

## Results for 2017–2023 (run of 08–09/11/2026)

Same protocol for every year: official neighbours of that year downloaded from WEkEO (72 files; two more truncated
transfers, 2018 AMPL T34UGU and 2020 RSLOPE T34UFV, were detected and re-downloaded), full-overlap comparison on
1 000 000 random pixels per strip, calibration fitted on that year's own strips. All 28 product-year files
(4 products × 7 years) are written in the MONIFUN layout; 98.9–99.1 % of the tile has a valid season in every year.
Tables: `reports/06_validation_full_overlap_all_years.csv`, `06_loso_calibration_all_years.csv`,
`06_per_strip_all_years.csv`, `calibration_YYYY.json`, `07_years_status.json`.

**Uncalibrated accuracy** (mean of the three strips; bias = reconstructed − official):

| year | AMPL r | AMPL bias | LENGTH r | LENGTH bias d | RSLOPE r | TPROD r | TPROD bias | TPROD rRMSE % |
|---|---|---|---|---|---|---|---|---|
| 2017 | 0.824 | +0.39 | 0.682 | +2.8 | 0.558 | 0.931 | +11.4 | 17.8 |
| 2018 | 0.858 | +0.37 | 0.668 | +7.9 | 0.619 | 0.926 | +10.2 | 17.3 |
| 2019 | 0.860 | +0.36 | 0.714 | +4.1 | 0.576 | 0.939 | +14.6 | 17.2 |
| 2020 | 0.785 | +0.55 | 0.623 | −4.1 | 0.494 | 0.896 | +21.1 | 23.4 |
| 2021 | 0.824 | +0.21 | 0.710 | +6.2 | 0.664 | 0.925 | +1.1 | 14.9 |
| 2022 | 0.844 | +0.25 | 0.700 | +2.2 | 0.557 | 0.929 | +4.2 | 17.1 |
| 2023 | 0.854 | +0.28 | 0.691 | +0.8 | 0.662 | 0.920 | +5.8 | 18.0 |

**Calibration actually applied** (official ≈ k · recon + c, all three strips of the year; stored in the GeoTIFF tag
`CALIBRATION`) and its leave-one-strip-out check (largest held-out abs(bias) after calibration):

| year | AMPL k, c | AMPL max held-out abs(bias) | RSLOPE k, c | RSLOPE max held-out abs(bias) |
|---|---|---|---|---|
| 2017 | 0.800, −0.025 | 0.091 | 0.894, 0.0057 | 0.0027 |
| 2018 | 0.881, −0.155 | 0.011 | 1.001, 0.0040 | 0.0017 |
| 2019 | 0.846, −0.079 | 0.041 | 0.908, 0.0060 | 0.0043 |
| 2020 | 0.781, −0.137 | 0.094 | 0.774, 0.0050 | 0.0030 |
| 2021 | 0.918, −0.040 | 0.098 | 1.163, 0.0067 | 0.0020 |
| 2022 | 0.951, −0.159 | 0.058 | 1.016, 0.0052 | 0.0032 |
| 2023 | 0.989, −0.263 | 0.055 | 1.139, 0.0033 | 0.0049 |

Mean held-out AMPL rRMSE falls from 24–51 % to 20–28 % after calibration; RSLOPE from 60–66 % to 52–64 %.

**Against the acceptance thresholds**

- **AMPL** passes in every year after calibration (held-out abs(bias) ≤ 0.10 in all 21 strip-years). Uncalibrated r
  is 0.79–0.86; 2020 is marginally below 0.80. The coefficients differ between years (k 0.78–0.99, c −0.26 to −0.03),
  which confirms that the calibration must be fitted per year, not once.
- **LENGTH** passes uncalibrated in every year (r 0.62–0.71, abs(bias) ≤ 8 days). Calibration is not applied: the
  held-out check is inconsistent (it worsens 2020 and 2023).
- **RSLOPE** passes r ≥ 0.50 in 6 of 7 years (2020: 0.494) and is unbiased after calibration, but remains the noisiest
  parameter (rRMSE 52–64 %). Use with caution.
- **TPROD** has the highest correlation (0.90–0.94) and passes in 6 of 7 years. **2020 fails** the rRMSE threshold
  uncalibrated (23.4 %, bias +21); a linear calibration would bring it to 17.3 % in the held-out check. In 2017–2019
  the positive bias (+10 to +15) is also visible. A year-specific TPROD calibration is therefore an option for WP3
  to decide; it is not applied in the delivered files.
- **2020** is the weakest year for every product (lower r, larger biases). Users should flag T34UGV plots in 2020
  or treat them with lower weight.

The positive AMPL bias in all years is not explained by the soil DVI (see the sensitivity test below). Why 2020 is
weaker has not been investigated yet (candidates: winter observation density, the 2020 spring drought).

## Iteration 2: sensitivity to the soil DVI (09/11/2026)

The constant soil DVI (0.09) was varied on three 4 × 4 km blocks, one in each overlap strip (≈ 478 000 pixels per
setting, batch exports; `reports/08_sensitivity_dvi_soil_2021.csv`, `..._2020.csv`). Values are means over the
three blocks, uncalibrated.

**2021**

| soil DVI | AMPL r / bias | LENGTH r / bias d | RSLOPE r | TPROD r / bias / rRMSE % | MINV = 0 % |
|---|---|---|---|---|---|
| 0.05 | 0.793 / +0.32 | 0.692 / +12.2 | 0.538 | 0.912 / +21.8 / 19.8 | 11 |
| 0.07 | 0.793 / +0.28 | 0.704 / +8.5 | 0.570 | 0.917 / +13.2 / 16.6 | 17 |
| **0.09** (used) | 0.792 / +0.24 | 0.711 / +4.6 | 0.595 | 0.921 / +4.4 / 14.5 | 33 |
| 0.11 | 0.789 / +0.20 | 0.753 / −0.9 | 0.612 | 0.928 / −5.1 / 13.8 | 52 |
| 0.13 | 0.787 / +0.15 | 0.790 / −6.4 | 0.627 | 0.933 / −14.7 / 15.5 | 68 |
| 0.15 | 0.783 / +0.10 | 0.800 / −11.4 | 0.639 | 0.935 / −24.2 / 19.0 | 80 |
| 0.18 | 0.777 / +0.01 | 0.785 / −18.7 | 0.649 | 0.933 / −38.7 / 25.9 | 92 |

**2020 (weakest year)**

| soil DVI | AMPL r / bias | LENGTH r / bias d | RSLOPE r | TPROD r / bias / rRMSE % |
|---|---|---|---|---|
| **0.09** (used) | 0.685 / +0.62 | 0.725 / −5.0 | 0.444 | 0.881 / +26.7 / 25.2 |
| 0.11 | 0.680 / +0.57 | 0.724 / −9.4 | 0.467 | 0.880 / +17.8 / 21.0 |
| 0.13 | 0.675 / +0.52 | 0.718 / −14.0 | 0.491 | 0.877 / +8.6 / 18.1 |

**Conclusions**

1. Raising the soil DVI shifts LENGTH down by ≈ 4–5.5 days and TPROD down by ≈ 9 units per 0.02 step, in both years.
   It moves the biases but does not remove their year-to-year differences: LENGTH is +4.6 days in 2021 and −5.0 days
   in 2020 at the same setting, so any single value that fixes one year worsens the other.
2. Correlations hardly change for AMPL and TPROD; RSLOPE gains slightly at higher values.
3. **Decision:** keep 0.09 and the delivered files unchanged. A per-year linear calibration on the overlap strips
   removes the same biases without re-running the reconstruction, and is already applied to AMPL and RSLOPE.
   Extending it to TPROD in 2020 (and possibly 2017–2019) is the cheaper fix if WP3 wants it.
4. Next candidates for the remaining AMPL bias: the peak side (MAXV, the tile-constant δ in MDVI, the median
   compositing). Testing them needs the official MAXV and MINV files of the neighbours.

## Iteration 3: correction learned from the official tiles (10/11/2026)

Prompted by the visible seam between the patch and the official tiles. Training data: the overlap strips with
T34UFV, T34UGA, T34UGU and, new, T35ULP (south-east corner, reprojected from UTM 35); 200 000 pixels per strip.
Leave-one-strip-out: each strip is predicted by a model trained on the other three. Means over the four held-out
strips, 2021 (`reports/10_learned_correction_loso_2021.csv`):

| product | method | r | abs(bias) | rRMSE % | std ratio (est / official) |
|---|---|---|---|---|---|
| AMPL | none | 0.817 | 0.226 | 24.8 | 0.86 |
| AMPL | linear (v1) | 0.817 | 0.070 | 21.7 | 0.80 |
| AMPL | quantile mapping | 0.814 | 0.098 | 22.9 | 0.98 |
| AMPL | gradient boosting | 0.859 | 0.047 | 19.1 | 0.85 |
| AMPL | **gradient boosting + quantile mapping** | 0.857 | 0.062 | 19.9 | 0.98 |
| LENGTH | none | 0.734 | 5.0 d | 18.0 | 1.02 |
| LENGTH | linear | 0.734 | 2.6 d | 16.6 | 0.74 |
| LENGTH | gradient boosting | 0.848 | 1.0 d | 12.4 | 0.84 |
| LENGTH | **gradient boosting + quantile mapping** | 0.843 | 1.6 d | 13.0 | 0.98 |
| RSLOPE | none | 0.664 | 0.009 | 64.0 | 0.58 |
| RSLOPE | linear (v1) | 0.664 | 0.002 | 53.5 | 0.68 |
| RSLOPE | gradient boosting | 0.852 | 0.001 | 37.1 | 0.85 |
| RSLOPE | **gradient boosting + quantile mapping** | 0.852 | 0.001 | 38.4 | 0.98 |
| TPROD | none | 0.924 | 4.3 | 14.7 | 0.92 |
| TPROD | linear | 0.924 | 5.8 | 14.8 | 0.92 |
| TPROD | gradient boosting | 0.937 | 4.0 | 13.4 | 0.94 |
| TPROD | **gradient boosting + quantile mapping** | 0.935 | 4.1 | 13.9 | 1.01 |

The multivariate model gains most where a single parameter of the reconstruction is a poor predictor of its official
counterpart: RSLOPE (r 0.66 → 0.85) and LENGTH (0.73 → 0.85); the other parameters, above all the season dates,
carry the information. Plain gradient boosting shrinks the spread (std ratio 0.84–0.94), which makes the patch look
smoother than its neighbours; the added quantile mapping restores the spread (0.98–1.01) at a small cost in RMSE.
**Decision:** gradient boosting + quantile mapping for all four products, fitted per year on that year's strips.

**Model size.** The table above uses 300 trees × 63 leaves. Predicting a whole tile (120 M pixels) with that model
takes about 1.5 h per product on the 2-CPU Colab runtime, so production uses 40 trees × 15 leaves: 9 times faster,
held-out rRMSE (2021) AMPL 19.3 instead of 19.0 %, RSLOPE 38.6 instead of 37.1 %, LENGTH 12.9 instead of 12.4 %.
A third-order polynomial ridge regression was also tested (AMPL 19.8, RSLOPE 40.4, LENGTH 13.1, TPROD 13.8 %) and is
slightly worse at similar speed.

### Production results 2017–2023 (10/11/2026)

Same protocol for every year, four strips (T34UFV, T34UGA, T34UGU, T35ULP), production model (40 trees × 15 leaves)
plus quantile mapping; held-out means (`reports/11_loso_all_years.csv`). rRMSE in %, r in brackets:

| year | AMPL none → corrected | LENGTH none → corrected | RSLOPE none → corrected | TPROD none → corrected |
|---|---|---|---|---|
| 2017 | 36.4 (0.81) → 24.3 (0.85) | 19.9 (0.70) → 15.9 (0.80) | 67.5 (0.53) → 49.5 (0.79) | 17.2 (0.93) → 16.0 (0.93) |
| 2018 | 34.2 (0.85) → 21.0 (0.87) | 20.0 (0.67) → 15.4 (0.80) | 67.0 (0.61) → 49.9 (0.81) | 17.2 (0.93) → 15.4 (0.93) |
| 2019 | 34.2 (0.85) → 19.6 (0.87) | 19.0 (0.74) → 13.4 (0.84) | 60.6 (0.58) → 39.0 (0.83) | 18.0 (0.94) → 13.5 (0.94) |
| 2020 | 51.0 (0.78) → 26.1 (0.81) | 21.6 (0.62) → 15.9 (0.79) | 64.2 (0.49) → 50.2 (0.75) | 23.2 (0.90) → 17.0 (0.90) |
| 2021 | 24.8 (0.82) → 20.2 (0.85) | 18.0 (0.73) → 13.7 (0.83) | 64.0 (0.66) → 40.2 (0.84) | 14.7 (0.92) → 14.8 (0.93) |
| 2022 | 30.0 (0.84) → 22.2 (0.88) | 21.9 (0.67) → 19.3 (0.77) | 65.2 (0.53) → 50.7 (0.75) | 17.8 (0.93) → 16.6 (0.94) |
| 2023 | 29.0 (0.85) → 22.4 (0.87) | 19.9 (0.70) → 15.5 (0.82) | 67.3 (0.67) → 46.5 (0.85) | 18.0 (0.91) → 18.2 (0.92) |

The standard deviation of the corrected values is 0.95–1.05 of the official one in every product-year (uncorrected
0.58–1.16), i.e. the patch now has the spread of its neighbours. The gain is largest for RSLOPE and LENGTH and for
the weakest year 2020; for TPROD the correction mainly removes the 2017–2020 bias and matches the spread (rRMSE is
within 0.3 points of the uncorrected value in 2021 and 2023). All 28 files 2017–2023 for T34UGV are in
`Copernicus_Phenology/` with the official band description, unit and PhysRange tag; the first delivery is kept in
`Copernicus_Phenology_v1_linear/`. 2024 and the extension T34UGVE depend on the Earth Engine exports still running.

## East extension T34UGVE 2018–2024 from a colleague's Earth Engine run (12/11/2026)

The extension (EPSG:32634, x 799760–869760, y 5390220–5560020, 7000 × 16980 px) was computed with
`gee/hrvpp_recon_patch.js` in the Code Editor under a colleague's Earth Engine quota (2018–2025; Earth Engine split
each year into three row shards). `notebooks/experiments/production_colleague_ee.py` mosaics the shards (grid
checked against the expected one for every year), learns the gbm_qm correction on the extension's own overlap with
the official T34UGA (10 km strip in the north-west) and T35ULP, and writes `Copernicus_Phenology/<P>/<Y>/
VPP_<Y>_S2_T34UGVE-010m_RECON_s1_<P>.tif`. The T34UGV raw exports used for the 2017 extension are no longer on
Drive, so the extension is calibrated on its own official neighbours here. 2025 is not processed: there is no
official HR-VPP 2025 to calibrate or check against.

Leave-one-neighbour-out (train on T34UGA, test on T35ULP and the other way round), mean rRMSE %,
uncorrected → corrected (`reports/13_colleague_loso_<year>.csv`):

| year | AMPL | LENGTH | RSLOPE | TPROD |
|---|---|---|---|---|
| 2018 | 44.6 → 29.1 | 28.8 → 20.2 | 61.9 → 49.3 | 25.5 → 20.8 |
| 2019 | 33.6 → 24.7 | 19.9 → 18.3 | 57.7 → 42.3 | 17.5 → 15.9 |
| 2020 | 40.2 → 25.2 | 27.8 → 18.5 | 58.0 → 44.8 | 21.8 → 18.4 |
| 2021 | 27.1 → 20.1 | 23.4 → 16.0 | 59.9 → 33.8 | 18.0 → 16.4 |
| 2022 | 38.7 → 27.4 | 23.8 → 20.8 | 62.5 → 47.8 | 21.1 → 18.5 |
| 2023 | 27.8 → 22.0 | 23.7 → 17.6 | 58.8 → 34.8 | 18.2 → 20.1 |
| 2024 | 29.6 → 23.4 | 22.4 → 17.7 | 59.7 → 33.2 | 19.6 → 20.8 |

For TPROD the correction does not help in 2023 and 2024 (the raw values are already close); it is kept for a
consistent spread across years.

Agreement with the corrected T34UGV in the 10 km overlap (`reports/13_colleague_<year>.json`, 500 000 pixels):
r 0.96–0.99 for AMPL, 0.97–0.98 LENGTH, 0.96–0.98 RSLOPE, 0.98–0.99 TPROD; rRMSE 6–12 %, 8–10 %, 16–21 %, 7–10 %.
The extension's season is 5–8 days longer than in T34UGV (T34UGV is calibrated on four strips, the extension on two),
so the seam is small but present in LENGTH.

Cross-check with the independent GPU pipeline (Sentinel-2 L2A from Planetary Computer, SCL mask, correction trained
on the T34UGV strips) on the same grid:

| year | level | AMPL r / rRMSE | LENGTH r / rRMSE / bias | RSLOPE r / rRMSE | TPROD r / rRMSE |
|---|---|---|---|---|---|
| 2021 | raw | 0.90 / 15.0 % | 0.89 / 15.7 % / +1.8 d | 0.81 / 38.2 % | 0.94 / 15.6 % |
| 2021 | corrected | 0.89 / 17.9 % | 0.93 / 13.0 % / −6.9 d | 0.91 / 28.3 % | 0.92 / 15.4 % |
| 2024 | raw | 0.97 / 6.9 % | 0.96 / 10.2 % / −0.8 d | 0.95 / 19.7 % | 0.99 / 6.2 % |
| 2024 | corrected | 0.95 / 13.0 % | 0.96 / 12.9 % / −8.7 d | 0.96 / 21.2 % | 0.96 / 13.0 % |

(bias = GPU − Earth Engine). The two implementations agree well on the raw metrics, best in 2024; after correction
the remaining difference is mostly the season length, where the GPU version is 7–9 days shorter because its model
comes from different strips. The Earth Engine version is kept as the delivery (better seam with T34UGV, calibrated
on the extension's own neighbours); the GPU files stay in `Copernicus_Phenology_gpu/` as an independent check.

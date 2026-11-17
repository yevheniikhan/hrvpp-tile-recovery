# Validation

This document describes how the accuracy of the reconstructed files was assessed and reports the results. Section 4
refers to the delivered files in `Copernicus_Phenology/`. The appendices document the earlier versions and the tests
that led to the final method.

## 1. Design

Every reconstructed area overlaps official HR-VPP tiles by about 10 km. In the overlap strips the official and the
reconstructed pixels lie on the same 10 m Sentinel-2 grid (T35ULP is reprojected from UTM zone 35 by nearest
neighbour), so the two can be compared pixel by pixel.

| area | neighbour | position | overlap strip (EPSG:32634) |
|---|---|---|---|
| T34UGV | T34UFV | west | x 699 960–709 800 |
| T34UGV | T34UGA | north | y 5 490 240–5 500 020 |
| T34UGV | T34UGU | south | y 5 390 220–5 400 000 |
| T34UGV | T35ULP | south-east corner | UTM 35 tile, reprojected |
| T34UGVE | T34UGA | north-west | 10 km strip |
| T34UGVE | T35ULP | south | 10 km strip |

The strips cover the Carpathian foothills, the Dniester valley, forest and arable land, i.e. the same landscape types
as the gap itself.

**Leave-one-strip-out cross-validation.** The correction model is trained on all strips but one and tested on the
remaining strip, which does not take part in training and thus represents the situation of the gap. The procedure is
repeated for each strip and the metrics are averaged. T34UGV has four strips, T34UGVE two (train on one, test on the
other). Each strip contributes 200 000 random pixels per year.

## 2. Metrics

Number of pixels, Pearson's r, bias (reconstructed − official), MAE, RMSE, relative RMSE (rRMSE, % of the official
mean), ordinary least squares slope and intercept, ratio of standard deviations (reconstructed / official) and the
share of pixels with a valid season in each source.

## 3. Implementation checks

1. **Earth Engine against NumPy.** The Earth Engine array code and the independent NumPy implementation
   (`phenology_np.py`) were applied to the same smoothed series at 291 random pixels. The results are identical
   (maximum absolute difference 0 for AMPL, LENGTH, RSLOPE, MAXV, MINV, SOS and EOS; 5.7 × 10⁻¹⁴ for TPROD).
2. **Unit tests.** The NumPy implementation is tested on synthetic double-logistic curves with analytically known
   start and end of season, amplitude and integral (`tests/test_phenology_np.py`); raster encoding, overlap logic and
   download checks are tested on synthetic GeoTIFFs (`tests/test_postprocess_validate.py`).
3. **Reproducibility.** Repeated Earth Engine sampling with a fixed seed gives identical results. Every downloaded
   official file is checked for completeness (`DOWNLOAD_INTEGRITY.md`), and SHA-256 checksums are stored so that a
   later download can be compared file by file.

## 4. Results for the delivered files

### 4.1 T34UGV, 2017–2024

Held-out rRMSE in %, before → after correction, with Pearson's r in brackets (`reports/11_loso_all_years.csv`;
for 2024 only rRMSE is tabulated):

| year | AMPL | LENGTH | RSLOPE | TPROD |
|---|---|---|---|---|
| 2017 | 36.4 (0.81) → 24.3 (0.85) | 19.9 (0.70) → 15.9 (0.80) | 67.5 (0.53) → 49.5 (0.79) | 17.2 (0.93) → 16.0 (0.93) |
| 2018 | 34.2 (0.85) → 21.0 (0.87) | 20.0 (0.67) → 15.4 (0.80) | 67.0 (0.61) → 49.9 (0.81) | 17.2 (0.93) → 15.4 (0.93) |
| 2019 | 34.2 (0.85) → 19.6 (0.87) | 19.0 (0.74) → 13.4 (0.84) | 60.6 (0.58) → 39.0 (0.83) | 18.0 (0.94) → 13.5 (0.94) |
| 2020 | 51.0 (0.78) → 26.1 (0.81) | 21.6 (0.62) → 15.9 (0.79) | 64.2 (0.49) → 50.2 (0.75) | 23.2 (0.90) → 17.0 (0.90) |
| 2021 | 24.8 (0.82) → 20.2 (0.85) | 18.0 (0.73) → 13.7 (0.83) | 64.0 (0.66) → 40.2 (0.84) | 14.7 (0.92) → 14.8 (0.93) |
| 2022 | 30.0 (0.84) → 22.2 (0.88) | 21.9 (0.67) → 19.3 (0.77) | 65.2 (0.53) → 50.7 (0.75) | 17.8 (0.93) → 16.6 (0.94) |
| 2023 | 29.0 (0.85) → 22.4 (0.87) | 19.9 (0.70) → 15.5 (0.82) | 67.3 (0.67) → 46.5 (0.85) | 18.0 (0.91) → 18.2 (0.92) |
| 2024 | 30.2 → 21.8 | 17.8 → 14.4 | 68.6 → 42.7 | 18.6 → 17.8 |

After correction the ratio of standard deviations is 0.95–1.05 for every parameter and year (0.58–1.16 before), so
the reconstructed tile has the spatial variability of its neighbours. The correction gives the largest gain for
RSLOPE and LENGTH and for 2020. For TPROD it mainly removes the positive bias of 2017–2020; in 2021 and 2023 the rRMSE
is unchanged within 0.3 percentage points. 98.9–99.1 % of the tile has a valid season in every year.

### 4.2 T34UGVE, 2018–2024

The extension was reconstructed with `gee/hrvpp_recon_patch.js` in the Earth Engine Code Editor. Earth Engine wrote
each year in three row shards, which were mosaicked with a grid check for every year
(`notebooks/experiments/production_colleague_ee.py`). The correction is trained on the extension's own overlap with
T34UGA and T35ULP. Held-out rRMSE in %, before → after correction (`reports/13_colleague_loso_<year>.csv`):

| year | AMPL | LENGTH | RSLOPE | TPROD |
|---|---|---|---|---|
| 2018 | 44.6 → 29.1 | 28.8 → 20.2 | 61.9 → 49.3 | 25.5 → 20.8 |
| 2019 | 33.6 → 24.7 | 19.9 → 18.3 | 57.7 → 42.3 | 17.5 → 15.9 |
| 2020 | 40.2 → 25.2 | 27.8 → 18.5 | 58.0 → 44.8 | 21.8 → 18.4 |
| 2021 | 27.1 → 20.1 | 23.4 → 16.0 | 59.9 → 33.8 | 18.0 → 16.4 |
| 2022 | 38.7 → 27.4 | 23.8 → 20.8 | 62.5 → 47.8 | 21.1 → 18.5 |
| 2023 | 27.8 → 22.0 | 23.7 → 17.6 | 58.8 → 34.8 | 18.2 → 20.1 |
| 2024 | 29.6 → 23.4 | 22.4 → 17.7 | 59.7 → 33.2 | 19.6 → 20.8 |

In 2023 and 2024 the correction slightly increases the TPROD error, because the uncorrected values are already close
to the official ones. It is applied in all years to keep the variance consistent across years. The extension for 2017
is corrected with the T34UGV model of the same year and has no separate held-out check.

### 4.3 Joint between T34UGVE and T34UGV

In the 10 km overlap of the two corrected reconstructions (500 000 pixels per year) r is 0.96–0.99 for AMPL,
0.97–0.98 for LENGTH, 0.96–0.98 for RSLOPE and 0.98–0.99 for TPROD; rRMSE is 6–12 %, 8–10 %, 16–21 % and 7–10 %
respectively. The season in T34UGVE is 5–8 days longer than in T34UGV, because the two areas are corrected with
models trained on different strips. The difference is small but can be seen in LENGTH maps.

### 4.4 Comparison of two independent implementations

T34UGVE was also reconstructed for 2021 and 2024 with the GPU implementation (`s2local.py`, Sentinel-2 Level-2A
from Microsoft Planetary Computer, SCL cloud mask, correction trained on the T34UGV strips). Agreement with the Earth
Engine version on the same grid (bias = GPU − Earth Engine):

| year | level | AMPL r / rRMSE | LENGTH r / rRMSE / bias | RSLOPE r / rRMSE | TPROD r / rRMSE |
|---|---|---|---|---|---|
| 2021 | raw | 0.90 / 15.0 % | 0.89 / 15.7 % / +1.8 d | 0.81 / 38.2 % | 0.94 / 15.6 % |
| 2021 | corrected | 0.89 / 17.9 % | 0.93 / 13.0 % / −6.9 d | 0.91 / 28.3 % | 0.92 / 15.4 % |
| 2024 | raw | 0.97 / 6.9 % | 0.96 / 10.2 % / −0.8 d | 0.95 / 19.7 % | 0.99 / 6.2 % |
| 2024 | corrected | 0.95 / 13.0 % | 0.96 / 12.9 % / −8.7 d | 0.96 / 21.2 % | 0.96 / 13.0 % |

The uncorrected parameters agree closely, best in 2024. After correction the main difference is season length, which
is 7–9 days shorter in the GPU version because its model was trained on different strips. The Earth Engine version is
delivered, since it joins T34UGV better and is corrected on the extension's own neighbours. The GPU files are kept in
`Copernicus_Phenology_gpu/` as an independent check.

### 4.5 Recommendations for use

The reconstructed files are intended as covariates for models, not as reference data. The following thresholds were
proposed for acceptance (to be agreed with WP3):

| parameter | proposed threshold |
|---|---|
| AMPL | r ≥ 0.80, abs(bias) ≤ 0.10 PPI |
| TPROD | r ≥ 0.80, rRMSE ≤ 20 % |
| LENGTH | r ≥ 0.60, abs(bias) ≤ 15 days |
| RSLOPE | r ≥ 0.50 |

For T34UGV in 2017–2023 the corrected values meet the r and rRMSE thresholds for all four parameters (section 4.1).
For T34UGVE the TPROD rRMSE is 20.8 % in 2018 and 2024, marginally above the proposed 20 %. RSLOPE remains the least
accurate parameter (rRMSE 33–51 %) and should be used with caution. The year 2020 is the least accurate for all parameters; plots in
the reconstructed areas should be flagged for this year or given a lower weight.

## Appendix A. First delivery: uncorrected values and linear calibration, 2017–2023

The first delivery (08–09/11/2026) used three overlap strips (T34UFV, T34UGA, T34UGU) and a linear calibration of
AMPL and RSLOPE only. These files are kept in `Copernicus_Phenology_v1_linear/`. For each year the official
neighbours were downloaded from WEkEO (72 files), and 1 000 000 random pixels per strip were compared.

**Uncorrected accuracy, 2021, by strip** (Earth Engine noncommercial project; 1 828 Sentinel-2 scenes in the 2021
window; 99.7 % of T34UGV with a valid season):

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

**Leave-one-strip-out linear calibration, 2021** (official = k · reconstructed + c, fitted on two strips and tested on
the third; mean of the three tests):

| product | r | bias raw → cal | rRMSE raw → cal % | decision |
|---|---|---|---|---|
| AMPL | 0.824 | +0.208 → −0.003 | 24.3 → 21.6 | **calibrate** (all-strip fit k = 0.918, c = −0.040) |
| RSLOPE | 0.664 | −0.0094 → −0.0002 | 63.3 → 52.1 | **calibrate** (k = 1.164, c = 0.0067) |
| LENGTH | 0.710 | +6.2 d → +0.5 d | 19.3 → 17.2 | keep uncalibrated: passes thresholds; k = 0.69 would compress natural variance |
| TPROD | 0.925 | +1.1 → 0.0 | 14.9 → 15.2 | keep uncalibrated: no gain |

**Uncorrected accuracy by year** (mean of the three strips):

| year | AMPL r | AMPL bias | LENGTH r | LENGTH bias d | RSLOPE r | TPROD r | TPROD bias | TPROD rRMSE % |
|---|---|---|---|---|---|---|---|---|
| 2017 | 0.824 | +0.39 | 0.682 | +2.8 | 0.558 | 0.931 | +11.4 | 17.8 |
| 2018 | 0.858 | +0.37 | 0.668 | +7.9 | 0.619 | 0.926 | +10.2 | 17.3 |
| 2019 | 0.860 | +0.36 | 0.714 | +4.1 | 0.576 | 0.939 | +14.6 | 17.2 |
| 2020 | 0.785 | +0.55 | 0.623 | −4.1 | 0.494 | 0.896 | +21.1 | 23.4 |
| 2021 | 0.824 | +0.21 | 0.710 | +6.2 | 0.664 | 0.925 | +1.1 | 14.9 |
| 2022 | 0.844 | +0.25 | 0.700 | +2.2 | 0.557 | 0.929 | +4.2 | 17.1 |
| 2023 | 0.854 | +0.28 | 0.691 | +0.8 | 0.662 | 0.920 | +5.8 | 18.0 |

**Linear calibration applied** (all three strips of the year; stored in the GeoTIFF tag `CALIBRATION`) and the
largest held-out abs(bias) after calibration:

| year | AMPL k, c | AMPL max held-out abs(bias) | RSLOPE k, c | RSLOPE max held-out abs(bias) |
|---|---|---|---|---|
| 2017 | 0.800, −0.025 | 0.091 | 0.894, 0.0057 | 0.0027 |
| 2018 | 0.881, −0.155 | 0.011 | 1.001, 0.0040 | 0.0017 |
| 2019 | 0.846, −0.079 | 0.041 | 0.908, 0.0060 | 0.0043 |
| 2020 | 0.781, −0.137 | 0.094 | 0.774, 0.0050 | 0.0030 |
| 2021 | 0.918, −0.040 | 0.098 | 1.163, 0.0067 | 0.0020 |
| 2022 | 0.951, −0.159 | 0.058 | 1.016, 0.0052 | 0.0032 |
| 2023 | 0.989, −0.263 | 0.055 | 1.139, 0.0033 | 0.0049 |

Calibration reduced the mean held-out rRMSE of AMPL from 24–51 % to 20–28 % and of RSLOPE from 60–66 % to 52–64 %.
The coefficients differ between years (AMPL k 0.78–0.99), so the calibration has to be fitted per year. LENGTH and
TPROD were left uncalibrated in this version; TPROD in 2020 exceeded the 20 % threshold (23.4 %).

The positive AMPL bias was first attributed to a clipped winter minimum (MINV = 0 where winter DVI falls below the
constant soil DVI). The sensitivity test in Appendix B does not support this: the bias increases as fewer pixels are
clipped. The bias therefore originates on the peak side (MAXV); this was not tested further because the official MAXV
files were not downloaded.

## Appendix B. Sensitivity to the soil DVI

The constant soil DVI (0.09) was varied on three 4 × 4 km blocks, one in each overlap strip (about 478 000 pixels per
setting; `reports/08_sensitivity_dvi_soil_2021.csv`, `..._2020.csv`). Values are means over the three blocks,
without correction.

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

**2020**

| soil DVI | AMPL r / bias | LENGTH r / bias d | RSLOPE r | TPROD r / bias / rRMSE % |
|---|---|---|---|---|
| **0.09** (used) | 0.685 / +0.62 | 0.725 / −5.0 | 0.444 | 0.881 / +26.7 / 25.2 |
| 0.11 | 0.680 / +0.57 | 0.724 / −9.4 | 0.467 | 0.880 / +17.8 / 21.0 |
| 0.13 | 0.675 / +0.52 | 0.718 / −14.0 | 0.491 | 0.877 / +8.6 / 18.1 |

Raising the soil DVI by 0.02 shortens LENGTH by about 4–5.5 days and lowers TPROD by about 9 units in both years. It
shifts the biases but does not remove their differences between years: at the same setting LENGTH is +4.6 days in
2021 and −5.0 days in 2020, so a value that suits one year degrades the other. Correlations of AMPL and TPROD hardly
change; RSLOPE improves slightly at higher values. The value 0.09 was kept, and the biases are removed by the
per-year correction instead.

## Appendix C. Choice of the correction method

Five options were compared by leave-one-strip-out cross-validation on 2021 with four strips (T34UFV, T34UGA, T34UGU,
T35ULP), 200 000 pixels per strip; means over the four held-out strips (`reports/10_learned_correction_loso_2021.csv`):

| product | method | r | abs(bias) | rRMSE % | std ratio (est / official) |
|---|---|---|---|---|---|
| AMPL | none | 0.817 | 0.226 | 24.8 | 0.86 |
| AMPL | linear (v1) | 0.817 | 0.070 | 21.7 | 0.80 |
| AMPL | quantile mapping | 0.814 | 0.098 | 22.9 | 0.98 |
| AMPL | gradient boosting | 0.859 | 0.047 | 19.1 | 0.85 |
| AMPL | gradient boosting + quantile mapping (selected) | 0.857 | 0.062 | 19.9 | 0.98 |
| LENGTH | none | 0.734 | 5.0 d | 18.0 | 1.02 |
| LENGTH | linear | 0.734 | 2.6 d | 16.6 | 0.74 |
| LENGTH | gradient boosting | 0.848 | 1.0 d | 12.4 | 0.84 |
| LENGTH | gradient boosting + quantile mapping (selected) | 0.843 | 1.6 d | 13.0 | 0.98 |
| RSLOPE | none | 0.664 | 0.009 | 64.0 | 0.58 |
| RSLOPE | linear (v1) | 0.664 | 0.002 | 53.5 | 0.68 |
| RSLOPE | gradient boosting | 0.852 | 0.001 | 37.1 | 0.85 |
| RSLOPE | gradient boosting + quantile mapping (selected) | 0.852 | 0.001 | 38.4 | 0.98 |
| TPROD | none | 0.924 | 4.3 | 14.7 | 0.92 |
| TPROD | linear | 0.924 | 5.8 | 14.8 | 0.92 |
| TPROD | gradient boosting | 0.937 | 4.0 | 13.4 | 0.94 |
| TPROD | gradient boosting + quantile mapping (selected) | 0.935 | 4.1 | 13.9 | 1.01 |

The multivariate model gains most for RSLOPE (r 0.66 → 0.85) and LENGTH (0.73 → 0.85), whose official values are
better predicted by the combination of reconstructed parameters and season dates than by the corresponding single
parameter. Gradient boosting alone reduces the variance (standard deviation ratio 0.84–0.94); the added quantile
mapping restores it (0.98–1.01) at a small cost in RMSE. Gradient boosting with quantile mapping was selected for all
four parameters, fitted per year.

**Model size.** The comparison above uses 300 trees with 63 leaves. Applying this model to a full tile (120 million
pixels) takes about 1.5 hours per parameter on a 2-CPU Colab runtime. The production model has 40 trees with 15
leaves and is nine times faster; its held-out rRMSE for 2021 is 19.3 % instead of 19.0 % for AMPL, 38.6 % instead
of 37.1 % for RSLOPE and 12.9 % instead of 12.4 % for LENGTH. A third-order polynomial ridge regression was also
tested (AMPL 19.8 %, RSLOPE 40.4 %, LENGTH 13.1 %, TPROD 13.8 %) and was slightly less accurate at similar speed.

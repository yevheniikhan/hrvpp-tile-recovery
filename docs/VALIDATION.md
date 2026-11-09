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

Earth Engine project `gen-lang-client-0200925587`; 1 828 Sentinel-2 scenes in the 2021 window; 99.7 % of T34UGV
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

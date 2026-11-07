# Method: reconstruction of HR-VPP season-1 parameters for tile T34UGV

## 1. Problem statement

The MONIFUN WP3 workflow downloads the Copernicus Land Monitoring Service (CLMS) High Resolution
Vegetation Phenology and Productivity (HR-VPP) parameters `AMPL`, `LENGTH`, `RSLOPE` and `TPROD`
(season 1, 2017–2023) through the WEkEO Harmonised Data Access (HDA) API
(dataset `EO:EEA:DAT:CLMS_HRVPP_VPP`, `productGroupId = s1`).

For western Ukraine the delivered mosaic has a rectangular hole. The hole corresponds exactly to the
Sentinel-2 MGRS tile **T34UGV**:

| item | value |
|---|---|
| CRS | WGS 84 / UTM 34N (EPSG:32634) |
| upper-left corner | x = 699 960 m, y = 5 500 020 m |
| size | 10 980 × 10 980 px at 10 m (109.8 km) |
| lon/lat footprint | ≈ 23.71–25.29° E, 48.59–49.62° N |

Evidence (08/11/2026):

* VITO OpenSearch catalogue (`https://phenology.vgt.vito.be/products`), collection
  `copernicus_r_utm-wgs84_10_m_hrvpp-vpp_p_2017-now_v01`: 0 products at 24.6° E, 49.1° N for 2017–2025;
  the same for the PPI Seasonal Trajectories collection. Neighbouring tiles T34UFV, T34UGA, T34UGU and
  T35ULP are present for every year.
* The CLMS WMTS (`phenology.vgt.vito.be/wmts`) shows no data inside the footprint for 2017–2024, in both the
  UTM and the LAEA layers.

The gap is therefore a production gap. Re-downloading cannot fill it.

## 2. Reconstruction algorithm

The algorithm follows the published HR-VPP documents wherever they give a definitive specification and
uses documented simplifications elsewhere.

| Stage | HR-VPP operational chain | This reconstruction |
|---|---|---|
| Input | Sentinel-2 L2A (MAJA/Sen2Cor), 10 m | Sentinel-2 L2A `COPERNICUS/S2_SR_HARMONIZED`; L1C only where no L2A scene exists (mostly 2017–2018), flagged in `L1C_SHARE` |
| Cloud mask | QFLAG2 from SCL with morphological buffers | Cloud Score+ `cs_cdf ≥ 0.60` plus SCL 0, 1, 3, 8, 9, 10 removed |
| Vegetation index | PPI (Jin & Eklundh 2014), ATBD §5.2.1 | Same formula (below) |
| MDVI | DVI q95 + δ, δ averaged over 100 sub-areas and interpolated, cap 0.8 | DVI q95 + tile-mean δ, cap 0.8 |
| Soil DVI | adaptive quantile per 10.98 km chunk | constant 0.09 (the alternative given in the ATBD) |
| Temporal model | TIMESAT 4.1 double-logistic fit, 2 seasons | 10-day median composites, winter base filling, linear interpolation (≤ 60 days), two 1-2-1 binomial passes |
| SOS / EOS | 25 % / 15 % of amplitude | same |
| Season-year | year of the peak | same |
| AMPL | MAXV − MINV, MINV = mean of left and right minima | same |
| TPROD | sum of daily PPI between SOS and EOS | integral of the 10-day series between SOS and EOS with fractional edge bins |
| RSLOPE | rate of change of PPI at EOS | absolute slope of the 10-day series around EOS |

**PPI** (ATBD §5.2.1):

```
DVI  = ρNIR(B08) − ρRed(B04)
PPI  = −K · ln[(MDVI − DVI) / (MDVI − DVIs)]       clipped to [0, 3]
K    = 1 / (4 · QE) · (1 + MDVI) / (1 − MDVI)
QE   = dc + (1 − dc) · G / cos θs,   dc = 0.0336 + 0.0477 / cos θs,   G = 0.5
```

θs is the mean solar zenith angle of the scene (`MEAN_SOLAR_ZENITH_ANGLE`).

**Season window.** For target year Y the series runs from 1 October Y−1 to 31 March Y+1 (54 bins of 10 days).
Only peaks inside year Y are considered. Bins from October to March without a clear observation receive the
5th-percentile base level (ATBD gap filling). Remaining gaps are linearly interpolated.

**Quality screening.** A pixel gets no season (nodata) if AMPL < 0.05, LENGTH < 30 days, or fewer than 8 clear
observations fall in the window.

## 3. Output and encoding

The Earth Engine export is an Int16 multi-band GeoTIFF on the exact T34UGV grid (bands AMPL, LENGTH,
RSLOPE, TPROD, SOSD, EOSD, NOBS; internal factors in `config.GEE_EXPORT_FACTORS`). `postprocess.py` converts each
band to the **dtype, scale, offset and nodata of the official neighbouring file** and writes:

```
Copernicus_Phenology/<PRODUCT>/<YEAR>/VPP_<YEAR>_S2_T34UGV-010m_RECON_s1_<PRODUCT>.tif
```

This is the folder layout expected by the MONIFUN extraction scripts (`_run_extract_xbyx_*.R`). Each file carries
the GeoTIFF tag `SOURCE = Reconstruction from Sentinel-2 (Earth Engine) - NOT an official CLMS product`.

Note on TPROD encoding: the Product User Manual gives a factor of 10 (scale 0.1), while the VITO catalogue band
metadata gives scale 0.02. The code never assumes either: it reads the scale from an official file.

## 4. Known limitations

* No TIMESAT double-logistic fitting. Expect slightly noisier SOS/EOS and larger RSLOPE scatter than the official
  product. AMPL and TPROD are least affected.
* Constant soil DVI and a tile-constant δ in MDVI. This mainly affects sparsely vegetated or bright-soil pixels.
* 2017 (and partly 2018) rely on L1C top-of-atmosphere reflectance where L2A is missing in Earth Engine. PPI from TOA
  is biased low. The validation for these years should be read separately.
* Only season 1 is reconstructed (as needed by MONIFUN).

## 5. References

* CLMS (2023). HR-VPP ATBD — Daily raw vegetation indices, issue 1.1 (PPI formula, §5.2.1).
* CLMS. HR-VPP Product User Manual — Seasonal Trajectories and VPP parameters, issue 2.3 (definitions, thresholds,
  encodings).
* Jin, H., Eklundh, L. (2014). A physically based vegetation index for improved monitoring of plant phenology.
  *Remote Sensing of Environment* 152, 512–525.
* Pasquarella, V. J., Brown, C. F., Czerwinski, W., Rucklidge, W. J. (2023). Comprehensive quality assessment of
  optical satellite imagery using weakly supervised video learning. CVPR Workshops (Cloud Score+).
* MONIFUN data wiki: *Auxiliary Data for WP3 update* (gitlab.com/monifun/monifun_data).

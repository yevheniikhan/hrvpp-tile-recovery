# Download integrity check for WEkEO transfers

**Script:** `download_copernicus_pheno.py` (wiki page *Auxiliary Data for WP3 (updated)*), function `download_one()`.

## Context

The MONIFUN download script works as intended when every transfer finishes. The case below arises only
when the connection to the WEkEO broker is interrupted during a download.

## What happens

`download_one()` decides whether a query is complete by comparing the number of files in
`OUTDIR/PRODUCT/YEAR` with the number of search results. A file that was only partly written (the HTTP
stream from the WEkEO broker was interrupted) is still a file, so it is counted as downloaded, and a later
run does not repeat it because the count already matches.

Such a file keeps a valid GeoTIFF header and looks normal in a file listing; the missing part shows up only
when those rows are read.

## Observed in this project (08–09/11/2026)

While downloading 84 official HR-VPP files (3 neighbouring tiles × 4 products × 7 years) through the same HDA
client, the stream broke 4 times. Example: `VPP_2018_S2_T34UGU-010m_V101_s1_AMPL.tif` stopped at 44.3 MB of
236 MB; `VPP_2020_S2_T34UFV-010m_V101_s1_RSLOPE.tif` at 43.9 MB. All four were
caught by the check below and downloaded again.

## Optional addition (as implemented in `src/hrvpp_recon/wekeo.py`)

1. After each pass, open every `.tif` with rasterio and read its **last row**; delete files that fail.
2. Count only the files that belong to the current query (the product ids of `matches.results`), so leftovers
   from another query cannot make a query look complete.
3. Catch exceptions from `matches.download()` inside the pass loop, so an interrupted stream leads to the
   next pass instead of aborting the product.

```python
def is_valid_tif(path):
    try:
        with rasterio.open(path) as src:
            src.read(1, window=Window(0, src.height - 1, src.width, 1))
        return True
    except Exception:
        return False
```

4. Optionally write SHA-256 checksums of the final files (`write_checksums()`), so a later re-download can be
   compared file by file.

Existing downloads can be confirmed in a few minutes by running `is_valid_tif()` once over the
`Copernicus_Phenology` folder; an empty list of failed files means everything is complete.

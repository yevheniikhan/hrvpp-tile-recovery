"""hrvpp_recon: detect, reconstruct and validate missing Copernicus HR-VPP tiles.

Modules
-------
config        tile geometry, product encodings, algorithm parameters
catalogue     anonymous availability check (VITO OpenSearch)
wekeo         WEkEO HDA download (MONIFUN-compatible)
gee_recon     Sentinel-2 -> PPI -> season parameters in Google Earth Engine
phenology_np  NumPy reference implementation (tests, cross-check of the EE code)
postprocess   EE export -> MONIFUN folder layout with official HR-VPP encoding
validate      pixel-wise comparison with official neighbouring tiles
"""
__version__ = "0.1.0"

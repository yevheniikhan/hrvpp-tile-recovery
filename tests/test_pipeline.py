"""End-to-end test of the shared pipeline on synthetic rasters (no Earth Engine, no network)."""
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from hrvpp_recon import config, pipeline  # noqa: E402
from test_learn import _official, _recon  # noqa: E402

YEAR = 2021


def test_grid_regions():
    g = config.REGIONS["T34UGVE"]
    assert g.bounds == (799760.0, 5390220.0, 869760.0, 5560020.0)
    assert config.Grid.from_tile(config.MISSING_TILE).bounds == config.MISSING_TILE.bounds
    assert set(config.REGION_NEIGHBOURS["T34UGVE"]) <= set(config.REGION_NEIGHBOURS["T34UGV"])


def test_mosaic_shards_keeps_grid_and_band_names(tmp_path):
    parts = []
    for i, r0 in enumerate([0, 6, 12]):
        h = 6 if r0 < 12 else 5
        p = tmp_path / f"RECON_VPP_2018_S2_X-010m_s1-{r0:010d}-0000000000.tif"
        with rasterio.open(p, "w", driver="GTiff", width=7, height=h, count=2, dtype="int16", crs="EPSG:32634",
                           transform=from_origin(799760, 5560020 - r0 * 10, 10, 10)) as t:
            t.write(np.full((2, h, 7), i + 1, "int16"))
            t.set_band_description(1, "AMPL")
            t.set_band_description(2, "NOBS")
        parts.append(p)
    out = pipeline.mosaic_shards(parts, tmp_path / "out" / "m.tif")
    with rasterio.open(out) as s:
        a = s.read()
        assert s.descriptions == ("AMPL", "NOBS")
        assert a.shape == (2, 17, 7) and s.transform.f == 5560020
        assert list(np.unique(a[0][:, 0])) == [1, 2, 3]


def test_calibrate_and_write(tmp_path):
    raw = tmp_path / "RECON_VPP_2021_S2_TEST-010m_s1.tif"
    ampl = _recon(raw)
    root = tmp_path / "official"
    # two official neighbours over the left and right thirds, official = 0.8 * recon - 0.05
    for tile, sl, x0 in (("T34AAA", slice(0, 70), 700000.0), ("T34BBB", slice(130, 200), 701300.0)):
        d = root / "AMPL" / str(YEAR)
        d.mkdir(parents=True, exist_ok=True)
        _official(d / f"VPP_{YEAR}_S2_{tile}-010m_V105_s1_AMPL.tif", 0.8 * ampl[:, sl] - 0.05, x0, 5500020.0)
    out = tmp_path / "Copernicus_Phenology"
    rep = pipeline.calibrate_and_write(raw, YEAR, root, ["T34AAA", "T34BBB", "T34CCC"], out, "TEST",
                                       products=["AMPL"], method="linear", report=tmp_path / "rep.json")
    assert rep["neighbours"] == ["T34AAA", "T34BBB"]
    lo = [r for r in rep["loso"] if r["method"] == "linear"]
    assert len(lo) == 2 and all(r["rrmse_pct"] < 2 for r in lo)
    f = pipeline.product_path(out, "AMPL", YEAR, "TEST")
    assert f.exists() and json.load(open(tmp_path / "rep.json"))["written"] == [str(f)]
    with rasterio.open(f) as s:
        v = s.read(1).astype(float) * s.scales[0]
        assert s.descriptions[0] == "Seasonal amplitude 2021-season1"
        assert "CALIBRATION" in s.tags()
    assert np.nanmax(np.abs(v - (0.8 * ampl - 0.05))) < 0.01
    # comparison helpers
    m = pipeline.compare_raw(raw, raw, YEAR, ["AMPL"], step=50, rows=20)["AMPL"]
    assert m["r"] > 0.999 and abs(m["bias"]) < 1e-9
    assert len(pipeline.summarise_loso([rep])) == 2

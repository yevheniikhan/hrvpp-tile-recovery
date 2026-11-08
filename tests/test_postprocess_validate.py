"""End-to-end test of export conversion and overlap validation on synthetic rasters."""
import os
import sys

import numpy as np
import pytest

rasterio = pytest.importorskip("rasterio")
from rasterio.transform import from_origin  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from hrvpp_recon.postprocess import gee_to_product, read_physical  # noqa: E402
from hrvpp_recon.validate import coverage_stats, metrics, overlap_arrays  # noqa: E402

PX = 300   # small synthetic tiles, 10 m pixels


def truth(x, y):
    return 0.8 + 0.4 * np.sin(x / 700.0) * np.cos(y / 900.0)


def make_official(path, ulx, uly):
    xs = ulx + 5 + 10 * np.arange(PX)
    ys = uly - 5 - 10 * np.arange(PX)
    X, Y = np.meshgrid(xs, ys)
    dn = np.round(truth(X, Y) / 1e-4).astype("int16")
    dn[:10, :10] = -32768
    with rasterio.open(path, "w", driver="GTiff", width=PX, height=PX, count=1, dtype="int16",
                       crs="EPSG:32634", transform=from_origin(ulx, uly, 10, 10), nodata=-32768) as d:
        d.write(dn, 1)
        d.scales = (1e-4,)


def make_recon(path, ulx, uly, noise=0.0):
    xs = ulx + 5 + 10 * np.arange(PX)
    ys = uly - 5 - 10 * np.arange(PX)
    X, Y = np.meshgrid(xs, ys)
    v = truth(X, Y) + noise * np.random.default_rng(1).standard_normal(X.shape)
    ampl = np.round(v * 1e4).astype("int16")
    length = np.full((PX, PX), 160, "int16")
    with rasterio.open(path, "w", driver="GTiff", width=PX, height=PX, count=2, dtype="int16",
                       crs="EPSG:32634", transform=from_origin(ulx, uly, 10, 10)) as d:
        d.write(ampl, 1)
        d.write(length, 2)
        d.set_band_description(1, "AMPL")
        d.set_band_description(2, "LENGTH")


def test_overlap_identical(tmp_path):
    off, rec = tmp_path / "off.tif", tmp_path / "rec.tif"
    make_official(off, 600000, 5500020)
    make_recon(rec, 600000 + 2000, 5500020)          # 1 km overlap
    a, b = overlap_arrays(rec, "AMPL", off, max_pixels=None)
    assert a.size == PX * (PX - 200)
    m = metrics(a, b)
    assert m["rmse"] < 1e-4 and m["r"] > 0.9999
    cov = coverage_stats(off, rec, "AMPL")
    assert cov["both_valid_pct"] > 99


def test_metrics_detect_noise(tmp_path):
    off, rec = tmp_path / "off.tif", tmp_path / "rec.tif"
    make_official(off, 600000, 5500020)
    make_recon(rec, 600000, 5500020, noise=0.05)
    m = metrics(*overlap_arrays(rec, "AMPL", off))
    assert 0.04 < m["rmse"] < 0.06 and abs(m["bias"]) < 0.01


def test_gee_to_product_uses_official_encoding(tmp_path):
    off, rec = tmp_path / "off.tif", tmp_path / "rec.tif"
    make_official(off, 600000, 5500020)
    make_recon(rec, 600000, 5500020)
    out = gee_to_product(rec, "AMPL", 2021, tmp_path / "CP", "T34UGV", official_tif=off, block=128)
    assert out.name == "VPP_2021_S2_T34UGV-010m_RECON_s1_AMPL.tif"
    assert out.parent.as_posix().endswith("CP/AMPL/2021")
    a = read_physical(out)
    b = read_physical(off)
    assert np.nanmax(np.abs(a - b)) < 1.5e-4
    with rasterio.open(out) as d:
        assert d.dtypes[0] == "int16" and d.nodata == -32768 and abs(d.scales[0] - 1e-4) < 1e-12
        assert "NOT an official" in d.tags()["SOURCE"]


def test_truncated_download_is_detected_and_removed(tmp_path):
    from hrvpp_recon import wekeo
    good = tmp_path / "VPP_2021_S2_T34UFV-010m_V105_s1_TPROD.tif"
    make_official(good, 600000, 5500020)
    bad = tmp_path / "VPP_2021_S2_T34UGA-010m_V105_s1_TPROD.tif"
    bad.write_bytes(good.read_bytes()[: good.stat().st_size // 3])
    assert wekeo.is_valid_tif(good) and not wekeo.is_valid_tif(bad)
    # only files of this query are counted; the truncated one is deleted
    assert wekeo._valid_count(tmp_path, {"VPP_2021_S2_T34UGA-010m_V105_s1_TPROD"}) == 0
    assert not bad.exists() and good.exists()


def test_strip_calibration_recovers_known_linear_bias(tmp_path):
    from hrvpp_recon.validate import strip_calibration
    root = tmp_path / "official"
    (root / "AMPL" / "2021").mkdir(parents=True)
    # three "neighbours" overlapping the reconstruction, each with official = (recon - 0.04) / 0.9
    rec = tmp_path / "rec.tif"
    make_recon(rec, 600000, 5500020)
    with rasterio.open(rec) as s:
        a = s.read(1).astype("float64") / 1e4
    for i, nb in enumerate(["TA", "TB", "TC"]):
        off = root / "AMPL" / "2021" / f"VPP_2021_S2_{nb}-010m_V105_s1_AMPL.tif"
        dn = np.round(((a - 0.04) / 0.9) / 1e-4).astype("int16")
        with rasterio.open(off, "w", driver="GTiff", width=PX, height=PX, count=1, dtype="int16",
                           crs="EPSG:32634", transform=from_origin(600000, 5500020, 10, 10), nodata=-32768) as d:
            d.write(dn, 1)
            d.scales = (1e-4,)
    fits, loso, raw = strip_calibration(rec, root, 2021, ["AMPL"], ["TA", "TB", "TC"], max_pixels=None)
    k, c = fits["AMPL"]
    assert abs(k - 1 / 0.9) < 1e-3 and abs(c + 0.04 / 0.9) < 1e-3
    assert all(abs(r["bias_cal"]) < 1e-3 for r in loso)
    out = gee_to_product(rec, "AMPL", 2021, tmp_path / "CP", "T34UGV", official_tif=off, block=128,
                         calibration=(k, c))
    with rasterio.open(out) as d:
        assert "official =" in d.tags()["CALIBRATION"]

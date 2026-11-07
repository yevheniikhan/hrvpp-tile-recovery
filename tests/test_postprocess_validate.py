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

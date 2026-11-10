"""Tests of the learned correction (synthetic rasters, no Earth Engine)."""
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from hrvpp_recon import learn, postprocess  # noqa: E402

YEAR = 2021


def _recon(path, x0=700000.0, y1=5500020.0, n=200, seed=1):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n]
    ampl = 1.5 + 0.9 * np.sin(xx / 15.0) * np.cos(yy / 19.0)  # spatially smooth like real fields
    bands = {
        "AMPL": ampl * 1e4, "LENGTH": rng.uniform(150, 250, (n, n)), "RSLOPE": -rng.uniform(0.005, 0.03, (n, n)) * 1e4,
        "TPROD": (ampl * 120) * 10, "SOSD": 21000 + rng.integers(80, 140, (n, n)),
        "EOSD": 21000 + rng.integers(250, 320, (n, n)), "NOBS": rng.integers(20, 60, (n, n)),
    }
    prof = dict(driver="GTiff", width=n, height=n, count=len(bands), dtype="int16", crs="EPSG:32634",
                transform=from_origin(x0, y1, 10, 10), nodata=-32768)
    with rasterio.open(path, "w", **prof) as dst:
        for i, (k, v) in enumerate(bands.items(), 1):
            dst.write(np.round(v).astype("int16"), i)
            dst.set_band_description(i, k)
    return ampl


def _official(path, values, x0, y1, crs="EPSG:32634"):
    prof = dict(driver="GTiff", width=values.shape[1], height=values.shape[0], count=1, dtype="int16", crs=crs,
                transform=from_origin(x0, y1, 10, 10), nodata=-32768)
    with rasterio.open(path, "w", **prof) as dst:
        dst.write(np.round(values / 1e-4).astype("int16"), 1)
        dst.scales = (1e-4,)
        dst.set_band_description(1, "Seasonal amplitude 2019-season1")
        dst.set_band_unit(1, "-")
        dst.update_tags(1, PhysRange="0 to 3")


def test_strip_samples_same_and_other_crs(tmp_path):
    rec = tmp_path / "rec.tif"
    ampl = _recon(rec)
    # same CRS, overlapping the left 50 columns; official = 0.8 * recon - 0.05
    off1 = tmp_path / "VPP_2021_S2_T34AAA-010m_V105_s1_AMPL.tif"
    _official(off1, 0.8 * ampl[:, :50] - 0.05, 700000.0, 5500020.0)
    # other CRS: reproject a field to zone 35 covering the right half
    off2 = tmp_path / "VPP_2021_S2_T35BBB-010m_V105_s1_AMPL.tif"
    src_t = from_origin(700000.0, 5500020.0, 10, 10)
    from rasterio.warp import calculate_default_transform
    t35, w35, h35 = calculate_default_transform("EPSG:32634", "EPSG:32635", 200, 200,
                                                *rasterio.transform.array_bounds(200, 200, src_t), resolution=10)
    dst = np.full((h35, w35), np.nan)
    reproject(0.8 * ampl - 0.05, dst, src_transform=src_t, src_crs="EPSG:32634", dst_transform=t35,
              dst_crs="EPSG:32635", resampling=Resampling.nearest, dst_nodata=np.nan)
    with rasterio.open(off2, "w", driver="GTiff", width=w35, height=h35, count=1, dtype="int16",
                       crs="EPSG:32635", transform=t35, nodata=-32768) as d:
        d.write(np.where(np.isfinite(dst), np.round(dst / 1e-4), -32768).astype("int16"), 1)
        d.scales = (1e-4,)
    df = learn.strip_samples(rec, {"A": {"AMPL": off1}, "B": {"AMPL": off2}}, YEAR, n_per_strip=10**6)
    assert set(df.strip) == {"A", "B"}
    for s in ("A", "B"):
        d = df[(df.strip == s) & df.OFF_AMPL.notna()]
        assert len(d) > 1000
        k, c = np.polyfit(d.AMPL, d.OFF_AMPL, 1)
        assert abs(k - 0.8) < 0.01 and abs(c + 0.05) < 0.01, (s, k, c)
    assert df.SOS_DAY.between(80, 140).all() and df.EOS_DAY.between(250, 320).all()


def test_correctors_and_loso():
    rng = np.random.default_rng(0)
    n = 30000
    x = rng.uniform(0.5, 2.5, n)
    df = {"AMPL": x, "LENGTH": rng.uniform(150, 250, n), "RSLOPE": rng.normal(-0.01, 0.003, n),
          "TPROD": x * 120, "SOS_DAY": rng.uniform(80, 140, n), "EOS_DAY": rng.uniform(250, 320, n),
          "NOBS": rng.uniform(20, 60, n)}
    import pandas as pd
    df = pd.DataFrame(df)
    df["OFF_AMPL"] = 0.8 * x ** 1.2 - 0.1 + rng.normal(0, 0.05, n)
    df["strip"] = np.repeat(["A", "B", "C"], n // 3)
    res = learn.loso(df, ["AMPL"], train_cap=20000)
    best = res.groupby("method").rmse.mean()
    assert best["gbm"] < best["linear"] < best["raw"]
    assert res[res.method == "quantile"].std_ratio.between(0.9, 1.1).all()
    choice = learn.choose_methods(res)
    assert choice["AMPL"] in ("gbm", "gbm_qm", "quantile")


def test_write_product_copies_official_metadata(tmp_path):
    rec = tmp_path / "rec.tif"
    ampl = _recon(rec)
    off = tmp_path / "off.tif"
    _official(off, ampl, 700000.0, 5500020.0)
    out = postprocess.gee_to_product(rec, "AMPL", 2024, tmp_path / "out", "T34UGV", official_tif=off)
    with rasterio.open(out) as s:
        assert s.descriptions[0] == "Seasonal amplitude 2024-season1"
        assert s.units[0] == "-" and s.tags(1)["PhysRange"] == "0 to 3"
        a = s.read(1).astype(float) * s.scales[0]
    assert np.allclose(a, np.round(ampl * 1e4) / 1e4, atol=1e-4)

"""The CuPy/NumPy port of the season metrics must equal phenology_np.extract_metrics."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from hrvpp_recon import phenology_np as pnp, s2local  # noqa: E402
from hrvpp_recon.config import ReconParams  # noqa: E402


def test_extract_metrics_xp_matches_reference():
    year = 2021
    t, _, n = pnp.window_bins(year)
    rng = np.random.default_rng(0)
    V = np.stack([pnp.double_logistic(t, base=rng.uniform(0, .3), amp=rng.uniform(.3, 2.5),
                                      sos=rng.uniform(80, 150), eos=rng.uniform(230, 300))
                  + rng.normal(0, .02, n) for _ in range(500)])
    V = pnp.smooth_121(V, 2)
    ref = pnp.extract_metrics(V, t, year)
    got = s2local.extract_metrics_xp(np, V, t, year, ReconParams())
    for k in ["AMPL", "LENGTH", "RSLOPE", "TPROD", "MAXV", "MINV", "SOSD", "EOSD"]:
        a = np.where(np.isnan(ref[k]), 0, ref[k]) if k not in ("SOSD", "EOSD") else ref[k]
        b = np.nan_to_num(got[k]) if k not in ("SOSD", "EOSD") else np.nan_to_num(got[k])
        assert np.allclose(a, b, atol=1e-4), k


def test_pipeline_on_synthetic_observations():
    """dvi observations -> metrics: a clean double-logistic season is recovered."""
    year = 2021
    p = ReconParams()
    t, ws, n = pnp.window_bins(year)
    days = np.arange(-92, 455, 5.0)                       # an observation every 5 days
    ppi_true = pnp.double_logistic(days, base=0.1, amp=1.5, sos=120, eos=280)
    obs_bin = ((days + 92) // 10).astype(int)
    obs_bin = np.where(obs_bin < n, obs_bin, -1)
    P = 50
    ppi = np.repeat(ppi_true[:, None], P, axis=1).astype("float32")
    ppi[::7, :] = np.nan                                  # some cloudy observations
    m = s2local.metrics_from_ppi(np, ppi, obs_bin, year, p)
    assert np.all(np.abs(m["AMPL"] - 1.5) < 0.15)
    assert np.all((m["LENGTH"] > 150) & (m["LENGTH"] < 230))
    assert np.all(m["NOBS"] >= p.min_obs)


def test_nanquantile_matches_numpy():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(40, 300)); a[rng.random(a.shape) < .3] = np.nan
    for q in (0.05, 0.5, 0.95):
        assert np.allclose(s2local.nanquantile(np, a, q, axis=0), np.nanquantile(a, q, axis=0), equal_nan=True)

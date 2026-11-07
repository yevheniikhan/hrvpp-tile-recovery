"""Unit tests of the reference season-metric extraction on synthetic curves with known answers."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from hrvpp_recon.phenology_np import (double_logistic, extract_metrics, smooth_121,  # noqa: E402
                                      window_bins, yydoy)


def dense_truth(year, base, amp, sos, eos, k1, k2, sos_frac=0.25, eos_frac=0.15):
    t = np.arange(-92, 457, 0.01)
    v = double_logistic(t, base, amp, sos, eos, k1, k2)
    in_year = (t >= 0) & (t < 365)
    peak = np.argmax(np.where(in_year, v, -1))
    lmin, rmin = v[:peak + 1].min(), v[peak:].min()
    maxv = v[peak]
    thr_s = lmin + sos_frac * (maxv - lmin)
    thr_e = rmin + eos_frac * (maxv - rmin)
    s = t[:peak + 1][np.argmax(v[:peak + 1] >= thr_s)]
    e = t[peak:][np.where(v[peak:] >= thr_e)[0].max()]
    m = (t >= s) & (t <= e)
    tprod = float(np.sum(v[m]) * 0.01)
    return dict(SOS_T=s, EOS_T=e, LENGTH=e - s, MAXV=maxv, AMPL=maxv - (lmin + rmin) / 2, TPROD=tprod)


@pytest.mark.parametrize("params", [
    dict(base=0.10, amp=1.20, sos=115, eos=285, k1=0.08, k2=0.06),   # broadleaved forest-like
    dict(base=0.35, amp=0.60, sos=95, eos=300, k1=0.05, k2=0.04),    # evergreen-like, small amplitude
    dict(base=0.05, amp=1.80, sos=140, eos=230, k1=0.12, k2=0.10),   # crop-like, short season
])
def test_metrics_match_dense_truth(params):
    year = 2021
    t, _, n = window_bins(year)
    assert n == 54
    v = double_logistic(t, **params)
    got = extract_metrics(v, t, year, min_ampl=0.0)
    truth = dense_truth(year, **params)
    assert abs(got["SOS_T"][0] - truth["SOS_T"]) < 2.5
    assert abs(got["EOS_T"][0] - truth["EOS_T"]) < 2.5
    assert abs(got["LENGTH"][0] - truth["LENGTH"]) < 4.0
    assert abs(got["AMPL"][0] - truth["AMPL"]) / truth["AMPL"] < 0.02
    assert abs(got["TPROD"][0] - truth["TPROD"]) / truth["TPROD"] < 0.02


def test_smoothing_preserves_constant_and_mass():
    v = np.full(54, 0.4)
    assert np.allclose(smooth_121(v, 3), 0.4)
    rng = np.random.default_rng(0)
    w = rng.random(54)
    s = smooth_121(w, 1)
    assert abs(s[1:-1].sum() - w[1:-1].sum()) < 0.5


def test_flat_series_has_no_season():
    t, _, _ = window_bins(2020)
    out = extract_metrics(np.full_like(t, 0.2), t, 2020)
    assert np.isnan(out["AMPL"][0]) and out["SOSD"][0] == 0


def test_yydoy_encoding():
    assert yydoy(np.array([0]), 2021)[0] == 21001
    assert yydoy(np.array([364]), 2021)[0] == 21365
    assert yydoy(np.array([-1]), 2021)[0] == 20366        # 31 Dec 2020 (leap year)
    assert yydoy(np.array([365]), 2021)[0] == 22001


def test_vectorised_rows_are_independent():
    year = 2022
    t, _, _ = window_bins(year)
    a = double_logistic(t, 0.1, 1.0, 120, 280, 0.08, 0.06)
    b = double_logistic(t, 0.2, 0.5, 100, 300, 0.05, 0.05)
    both = extract_metrics(np.vstack([a, b]), t, year)
    one = extract_metrics(b, t, year)
    for k in ("SOS_T", "EOS_T", "TPROD", "AMPL"):
        assert np.isclose(both[k][1], one[k][0])

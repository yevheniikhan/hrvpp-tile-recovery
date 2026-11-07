"""NumPy reference implementation of the season-metric extraction.

It mirrors, step by step, the array logic used in Earth Engine (gee_recon.py) so that
(i) the algorithm can be unit-tested on synthetic curves with known answers, and
(ii) the Earth Engine output can be cross-checked pixel by pixel.

Conventions
-----------
* V : smoothed, gap-free PPI time series, shape (N,) or (P, N); bins are `bin_days` apart.
* t : bin centres in days relative to 1 January of the target year (t[0] is negative
      because the analysis window starts on 1 October of the previous year).
* A season belongs to the year in which its peak occurs (HR-VPP rule).
"""
from __future__ import annotations

import calendar
from typing import Dict

import numpy as np


def smooth_121(v: np.ndarray, passes: int = 2) -> np.ndarray:
    """Binomial 1-2-1 smoothing along the last axis with edge replication."""
    out = np.asarray(v, dtype="float64")
    for _ in range(passes):
        left = np.concatenate([out[..., :1], out[..., :-1]], axis=-1)
        right = np.concatenate([out[..., 1:], out[..., -1:]], axis=-1)
        out = (left + 2 * out + right) / 4.0
    return out


def yydoy(t: np.ndarray, year: int) -> np.ndarray:
    """Convert day offsets (relative to 1 Jan `year`) to the HR-VPP YYDOY code."""
    t = np.floor(np.asarray(t, dtype="float64"))
    len_y = 366 if calendar.isleap(year) else 365
    len_prev = 366 if calendar.isleap(year - 1) else 365
    yy = np.where(t < 0, year - 1, np.where(t >= len_y, year + 1, year)) % 100
    doy = np.where(t < 0, t + len_prev + 1, np.where(t >= len_y, t - len_y + 1, t + 1))
    return (yy * 1000 + doy).astype("int64")


def extract_metrics(V: np.ndarray, t: np.ndarray, year: int, *, sos_frac=0.25, eos_frac=0.15,
                    bin_days=10, min_ampl=0.05, min_length=30) -> Dict[str, np.ndarray]:
    """Extract season-1 VPP parameters from smoothed dekadal PPI series.

    Returns a dict of float arrays (NaN = no valid season), keys:
    SOS_T, EOS_T, SOSD, EOSD, LENGTH, MAXV, MINV, AMPL, TPROD, SPROD, LSLOPE, RSLOPE.
    """
    V = np.atleast_2d(np.asarray(V, dtype="float64"))
    t = np.asarray(t, dtype="float64")
    P, N = V.shape
    idx = np.arange(N)
    len_y = 366 if calendar.isleap(year) else 365
    half = bin_days / 2.0

    in_year = (t >= 0) & (t < len_y)
    Vm = np.where(in_year, V, -1.0)
    maxv = Vm.max(axis=1)
    peak = np.where(Vm == maxv[:, None], idx, 999).min(axis=1)          # first max

    left = idx[None, :] <= peak[:, None]
    right = idx[None, :] >= peak[:, None]
    lmin = np.where(left, V, 99.0).min(axis=1)
    rmin = np.where(right, V, 99.0).min(axis=1)
    lmin_i = np.where(left & (V == lmin[:, None]), idx, -1).max(axis=1)   # last min before peak
    rmin_i = np.where(right & (V == rmin[:, None]), idx, 999).min(axis=1)  # first min after peak

    minv = (lmin + rmin) / 2.0
    ampl = maxv - minv
    thr_s = lmin + sos_frac * (maxv - lmin)
    thr_e = rmin + eos_frac * (maxv - rmin)

    rng_s = (idx[None, :] >= lmin_i[:, None]) & left
    rng_e = (idx[None, :] <= rmin_i[:, None]) & right
    sos_i = np.where(rng_s & (V >= thr_s[:, None]), idx, 999).min(axis=1)
    eos_i = np.where(rng_e & (V >= thr_e[:, None]), idx, -1).max(axis=1)

    rows = np.arange(P)
    s0 = np.clip(np.maximum(sos_i - 1, lmin_i), 0, N - 1)
    s1 = np.clip(sos_i, 0, N - 1)
    e0 = np.clip(eos_i, 0, N - 1)
    e1 = np.clip(np.minimum(eos_i + 1, rmin_i), 0, N - 1)
    vs0, vs1 = V[rows, s0], V[rows, s1]
    ve0, ve1 = V[rows, e0], V[rows, e1]

    with np.errstate(divide="ignore", invalid="ignore"):
        fs = np.clip(np.where(vs1 != vs0, (thr_s - vs0) / (vs1 - vs0), 1.0), 0, 1)
        fe = np.clip(np.where(ve0 != ve1, (ve0 - thr_e) / (ve0 - ve1), 0.0), 0, 1)
    t0 = t[0]
    sos_t = t0 + bin_days * (s0 + fs * (s1 - s0))
    eos_t = t0 + bin_days * (e0 + fe * (e1 - e0))
    length = eos_t - sos_t

    lo = np.maximum(t[None, :] - half, sos_t[:, None])
    hi = np.minimum(t[None, :] + half, eos_t[:, None])
    w = np.clip(hi - lo, 0, None)
    tprod = (V * w).sum(axis=1)
    sprod = tprod - minv * length

    lslope = np.abs(vs1 - vs0) / bin_days
    rslope = np.abs(ve0 - ve1) / bin_days

    valid = (ampl >= min_ampl) & (length >= min_length) & (sos_i < 999) & (eos_i >= 0) \
        & (sos_i <= peak) & (eos_i >= peak)

    out = {
        "SOS_T": sos_t, "EOS_T": eos_t, "LENGTH": length, "MAXV": maxv, "MINV": minv,
        "AMPL": ampl, "TPROD": tprod, "SPROD": sprod, "LSLOPE": lslope, "RSLOPE": rslope,
    }
    out = {k: np.where(valid, v, np.nan) for k, v in out.items()}
    out["SOSD"] = np.where(valid, yydoy(np.nan_to_num(sos_t), year), 0)
    out["EOSD"] = np.where(valid, yydoy(np.nan_to_num(eos_t), year), 0)
    out["PEAK_I"] = peak
    return out


def window_bins(year: int, bin_days: int = 10, start=(10, 1), end=(3, 31)):
    """Bin centres (days relative to 1 Jan `year`) for the window 1 Oct (Y-1) .. 31 Mar (Y+1)."""
    import datetime as dt

    ws = dt.date(year - 1, *start)
    we = dt.date(year + 1, *end)
    jan1 = dt.date(year, 1, 1)
    nbins = ((we - ws).days + 1) // bin_days
    first_centre = (ws - jan1).days + bin_days / 2.0
    return first_centre + bin_days * np.arange(nbins), ws, nbins


def double_logistic(t, base=0.1, amp=1.0, sos=120.0, eos=280.0, k1=0.08, k2=0.06):
    """Synthetic PPI curve used in the unit tests."""
    t = np.asarray(t, dtype="float64")
    return base + amp * (1 / (1 + np.exp(-k1 * (t - sos))) - 1 / (1 + np.exp(-k2 * (t - eos))))

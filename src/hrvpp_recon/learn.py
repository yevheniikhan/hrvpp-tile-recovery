"""Learn the correction of the reconstruction from the official tiles in the overlap strips.

The reconstruction is not a copy of the operational TIMESAT processing, so its values differ
systematically from the official product (docs/VALIDATION.md). Where the reconstructed tile
overlaps official tiles (~10 km strips), both are known for millions of pixels. This module
uses those pixels as training data:

* features  = all bands of the reconstruction (AMPL, LENGTH, RSLOPE, TPROD, SOS and EOS day, NOBS)
* target    = the official value of one product
* methods   = linear (k*x + c, the earlier calibration), quantile mapping (matches the full
              distribution, keeps texture), gradient boosting (multivariate, non-linear),
              and gradient boosting followed by quantile mapping of its output
* evaluation = leave-one-strip-out: fit on all strips but one, test on the held-out strip,
              which is the situation of the gap (no official data at all).

Official tiles from another UTM zone (e.g. T35ULP for a zone-34 reconstruction) are
reprojected on the fly to the reconstruction grid with nearest-neighbour resampling.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np
import pandas as pd
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.windows import from_bounds, Window
from rasterio.warp import transform_bounds

from .config import GEE_EXPORT_FACTORS
from .validate import metrics

GEE_NODATA = -32768
FEATURES = ["AMPL", "LENGTH", "RSLOPE", "TPROD", "SOS_DAY", "EOS_DAY", "NOBS"]
METHODS = ["raw", "linear", "quantile", "gbm", "gbm_qm"]


# --------------------------------------------------------------------------- reading
def yydoy_to_day(v: np.ndarray, year: int) -> np.ndarray:
    """YYDOY -> day relative to 1 January of `year` (EOS may fall in the next year)."""
    v = np.asarray(v, dtype="float64")
    yy, doy = np.floor(v / 1000.0), np.mod(v, 1000.0)
    return (yy - (year % 100)) * 365.0 + doy


def read_features(src: rasterio.io.DatasetReader, year: int, window: Optional[Window] = None) -> Dict[str, np.ndarray]:
    """Physical feature arrays from a raw Earth Engine export (NaN where no season)."""
    desc = list(src.descriptions)
    out = {}
    for name, band in [("AMPL", "AMPL"), ("LENGTH", "LENGTH"), ("RSLOPE", "RSLOPE"), ("TPROD", "TPROD"),
                       ("SOS_DAY", "SOSD"), ("EOS_DAY", "EOSD"), ("NOBS", "NOBS")]:
        if band not in desc:
            continue
        a = src.read(desc.index(band) + 1, window=window).astype("float64")
        a[a == GEE_NODATA] = np.nan
        if band in GEE_EXPORT_FACTORS:
            a = a / GEE_EXPORT_FACTORS[band]
        elif band in ("SOSD", "EOSD"):
            a = yydoy_to_day(a, year)
        out[name] = a
    valid = np.isfinite(out["AMPL"])
    for k in out:
        out[k] = np.where(valid, out[k], np.nan)
    return out


def _official_on_grid(official_tif: Path, rec: rasterio.io.DatasetReader):
    """(vrt on the reconstruction grid, raw dataset). Nearest neighbour; identity when the grids match."""
    raw = rasterio.open(official_tif)
    vrt = WarpedVRT(raw, crs=rec.crs, transform=rec.transform, width=rec.width, height=rec.height,
                    resampling=Resampling.nearest, nodata=raw.nodata)
    return vrt, raw


def strip_samples(recon_tif: Path, official: Dict[str, Dict[str, Path]], year: int,
                  n_per_strip: int = 200_000, seed: int = 0) -> pd.DataFrame:
    """Pixels of all overlap strips: features of the reconstruction + official targets.

    `official` = {neighbour: {product: path}}.
    """
    rng = np.random.default_rng(seed)
    frames = []
    with rasterio.open(recon_tif) as rec:
        for nb, prods in official.items():
            with rasterio.open(next(iter(prods.values()))) as raw:
                b = transform_bounds(raw.crs, rec.crs, *raw.bounds)
            left, bottom = max(b[0], rec.bounds.left), max(b[1], rec.bounds.bottom)
            right, top = min(b[2], rec.bounds.right), min(b[3], rec.bounds.top)
            if right <= left or top <= bottom:
                continue
            w = from_bounds(left, bottom, right, top, rec.transform).round_offsets().round_lengths()
            w = w.intersection(Window(0, 0, rec.width, rec.height))
            feats = read_features(rec, year, w)
            df = pd.DataFrame({k: v.ravel() for k, v in feats.items()})
            for prod, path in prods.items():
                vrt, raw = _official_on_grid(path, rec)
                with raw, vrt:
                    a = vrt.read(1, window=w).astype("float64")
                    nod, sc, of = raw.nodata, raw.scales[0] or 1.0, raw.offsets[0] or 0.0
                if nod is not None:
                    a[a == nod] = np.nan
                df["OFF_" + prod] = (a * sc + of).ravel()
            df = df[np.isfinite(df["AMPL"])]
            offc = [c for c in df.columns if c.startswith("OFF_")]
            df = df[df[offc].notna().any(axis=1)]
            if len(df) > n_per_strip:
                df = df.iloc[rng.choice(len(df), n_per_strip, replace=False)]
            df["strip"] = nb
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# --------------------------------------------------------------------------- models
class QuantileMap:
    """Monotone mapping of one variable onto the official distribution."""

    def __init__(self, n_q: int = 201):
        self.q = np.linspace(0, 1, n_q)

    def fit(self, x: np.ndarray, y: np.ndarray):
        m = np.isfinite(x) & np.isfinite(y)
        self.xq = np.quantile(x[m], self.q)
        self.yq = np.quantile(y[m], self.q)
        # strictly increasing x for np.interp
        self.xq = np.maximum.accumulate(self.xq + np.arange(len(self.xq)) * 1e-12)
        return self

    def predict(self, x: np.ndarray) -> np.ndarray:
        out = np.interp(x, self.xq, self.yq)
        # linear extrapolation beyond the fitted range keeps extremes ordered
        lo, hi = x < self.xq[0], x > self.xq[-1]
        if lo.any():
            s = (self.yq[5] - self.yq[0]) / max(self.xq[5] - self.xq[0], 1e-12)
            out[lo] = self.yq[0] + s * (x[lo] - self.xq[0])
        if hi.any():
            s = (self.yq[-1] - self.yq[-6]) / max(self.xq[-1] - self.xq[-6], 1e-12)
            out[hi] = self.yq[-1] + s * (x[hi] - self.xq[-1])
        return out


class Corrector:
    """Fit one correction method for one product; predict on feature dicts/frames."""

    def __init__(self, method: str, product: str, seed: int = 0):
        assert method in METHODS
        self.method, self.product, self.seed = method, product, seed

    def _X(self, d) -> np.ndarray:
        return np.column_stack([np.asarray(d[f], dtype="float64").ravel() for f in FEATURES if f in self.features])

    def fit(self, df: pd.DataFrame):
        y = df["OFF_" + self.product].to_numpy("float64")
        x = df[self.product].to_numpy("float64")
        m = np.isfinite(y) & np.isfinite(x)
        df, x, y = df[m], x[m], y[m]
        self.features = [f for f in FEATURES if f in df.columns]
        if self.method == "linear":
            self.k, self.c = np.polyfit(x, y, 1)
        elif self.method == "quantile":
            self.qm = QuantileMap().fit(x, y)
        elif self.method in ("gbm", "gbm_qm"):
            from sklearn.ensemble import HistGradientBoostingRegressor
            # small model: ~9x faster prediction than 300 trees x 63 leaves for <= 1.5 points rRMSE
            # (2021 leave-one-strip-out; whole tiles are 120 M pixels on a 2-CPU Colab runtime)
            self.gbm = HistGradientBoostingRegressor(max_iter=40, learning_rate=0.25, max_leaf_nodes=15,
                                                     min_samples_leaf=200, l2_regularization=1.0,
                                                     random_state=self.seed)
            self.gbm.fit(self._X(df), y)
            if self.method == "gbm_qm":
                self.qm = QuantileMap().fit(self.gbm.predict(self._X(df)), y)
        return self

    def predict(self, d) -> np.ndarray:
        x = np.asarray(d[self.product], dtype="float64")
        shape = x.shape
        x = x.ravel()
        out = np.full(x.shape, np.nan)
        ok = np.isfinite(x)
        if self.method == "raw":
            out[ok] = x[ok]
        elif self.method == "linear":
            out[ok] = self.k * x[ok] + self.c
        elif self.method == "quantile":
            out[ok] = self.qm.predict(x[ok])
        else:
            X = self._X(d)[ok]
            p = self.gbm.predict(X)
            out[ok] = self.qm.predict(p) if self.method == "gbm_qm" else p
        return out.reshape(shape)


def loso(df: pd.DataFrame, products: List[str], methods: List[str] = METHODS,
         train_cap: int = 300_000, seed: int = 0) -> pd.DataFrame:
    """Leave-one-strip-out comparison of the correction methods."""
    rng = np.random.default_rng(seed)
    rows = []
    strips = sorted(df["strip"].unique())
    for prod in products:
        col = "OFF_" + prod
        if col not in df:
            continue
        d = df[np.isfinite(df[col])]
        for test in strips:
            tr, te = d[d.strip != test], d[d.strip == test]
            if len(te) < 1000 or tr.strip.nunique() < 1:
                continue
            if len(tr) > train_cap:
                tr = tr.iloc[rng.choice(len(tr), train_cap, replace=False)]
            ref = te[col].to_numpy("float64")
            for meth in methods:
                est = Corrector(meth, prod, seed).fit(tr).predict(te)
                m = metrics(ref, est)
                rows.append({"product": prod, "test_strip": test, "method": meth, "n": m.get("n"),
                             "r": m.get("r"), "bias": m.get("bias"), "rmse": m.get("rmse"),
                             "rrmse_pct": m.get("rrmse_pct"),
                             "std_ratio": float(np.nanstd(est) / np.nanstd(ref)),
                             "mean_ref": m.get("mean_ref")})
    return pd.DataFrame(rows)


def choose_methods(res: pd.DataFrame, std_tol: float = 0.15, prefer: str = "rmse") -> Dict[str, str]:
    """Per product: lowest mean held-out RMSE among methods whose std ratio stays within
    1 +/- std_tol (so the corrected tile keeps the texture of the official neighbours)."""
    out = {}
    for prod, g in res.groupby("product"):
        s = g.groupby("method").agg(rmse=("rmse", "mean"), std_ratio=("std_ratio", "mean"),
                                    bias=("bias", lambda v: np.mean(np.abs(v))))
        ok = s[(s.std_ratio - 1).abs() <= std_tol]
        out[prod] = (ok if len(ok) else s)[prefer].idxmin()
    return out


# --------------------------------------------------------------------------- writing
def predict_tile(recon_tif: Path, year: int, correctors: Dict[str, Corrector], block: int = 2048):
    """Generator of (product, window, physical array) over the whole reconstructed tile."""
    with rasterio.open(recon_tif) as src:
        for r in range(0, src.height, block):
            for c in range(0, src.width, block):
                w = Window(c, r, min(block, src.width - c), min(block, src.height - r))
                feats = read_features(src, year, w)
                for prod, cor in correctors.items():
                    yield prod, w, cor.predict(feats)


def phys_fn_for(corrector: Corrector, year: int) -> Callable:
    """Callable(src, window) -> physical array, for postprocess.write_product."""
    def f(src, window):
        return corrector.predict(read_features(src, year, window))
    return f

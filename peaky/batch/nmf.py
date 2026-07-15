"""Error-weighted non-negative matrix factorisation (NMF / PMF-lite) of a batch
time-series into a small number of co-varying COMPONENTS (factors).

Why this on top of `cluster.py`: correlation clustering gives each ion a HARD
assignment to one family. NMF gives each ion FRACTIONAL membership across k
factors plus a quantitative per-factor time series -- the natural unit for
apportionment. X (time x peak) ~= G (time x k) . H (k x peak): each factor is a
mass-spectral profile H[f] with a time series G[:,f].

The one thing that makes it useful rather than a toy: ERROR WEIGHTING. Plain NMF
minimises ||X - GH||^2, so a 500k-cps background ion contributes ~10^4x more to
the loss than a 50-cps trace and dominates every factor. We minimise the
PMF objective sum_ij Wt_ij (X_ij - (GH)_ij)^2 with Wt = 1/sigma^2 and a
proportional+floor error model sigma = sqrt((rel*X)^2 + floor^2), so every ion
contributes ~equally regardless of brightness and below-detection points are
down-weighted. Linear intensity space (factors are additive spectra; log would
break additivity).

Pure numpy/pandas/scipy -- no scikit-learn dependency. Reuses
`timeseries.build_matrix` for the peak x time matrix.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from .timeseries import build_matrix, DEFAULT_TOL_PPM

__version__ = "0.1.0"

# error model: sigma_ij = sqrt((REL*X)^2 + floor^2); floor = FLOOR_FRAC * median(X>0)
REL_ERR = 0.15          # proportional error (fractional) -- equalises bright vs dim ions
FLOOR_FRAC = 0.10       # absolute floor as a fraction of the median positive intensity
MISSING_SIGMA_MULT = 3.0  # not-detected cells get this * max-sigma -> strongly down-weighted
N_ITER = 500
SEED = 0
EPS = 1e-10


# ---------------------------------------------------------------------------
# input: X (samples x bin) + per-cell error sigma
# ---------------------------------------------------------------------------
def build_input(peaks: pd.DataFrame, *, tol_ppm: float = DEFAULT_TOL_PPM,
                reagent_series: pd.Series | None = None, min_presence: float = 0.0,
                mz_col="mz", height_col="height", sample_col="sample_item_id",
                time_col="datetime_utc"):
    """Peaks -> (X, sigma, bin_mz, times, cols). X is a dense samples x bin
    intensity matrix (non-detections -> 0). sigma is the matching PMF error
    matrix. `reagent_series` (per sample_item_id) optionally divides each row ->
    a concentration proxy. `min_presence` drops bins detected in fewer than this
    fraction of samples (standard PMF species screening; also cuts run time)."""
    mat, bin_mz = build_matrix(peaks, tol_ppm=tol_ppm, mz_col=mz_col,
                               height_col=height_col, sample_col=sample_col)
    if mat.empty:
        raise ValueError("build_matrix returned an empty matrix")
    if min_presence > 0:
        keep = mat.notna().mean(axis=0) >= min_presence
        mat = mat.loc[:, keep[keep].index]
        bin_mz = bin_mz.reindex(mat.columns)
    # sample -> datetime (one row per sample)
    tmap = (peaks.dropna(subset=[sample_col, time_col])
                 .drop_duplicates(sample_col).set_index(sample_col)[time_col])
    times = pd.to_datetime(mat.index.map(tmap), utc=True)
    # optional reagent normalisation (concentration proxy)
    if reagent_series is not None:
        rs = reagent_series.reindex(mat.index).replace(0, np.nan)
        mat = mat.div(rs, axis=0)
    observed = mat.notna().to_numpy()
    X = np.array(mat.fillna(0.0).to_numpy(dtype=float), copy=True)
    X[X < 0] = 0.0
    # PMF-lite error model
    med = np.median(X[X > 0]) if np.any(X > 0) else 1.0
    floor = FLOOR_FRAC * med
    sigma = np.sqrt((REL_ERR * X) ** 2 + floor ** 2)
    sigma[~observed] = MISSING_SIGMA_MULT * np.nanmax(sigma)  # down-weight non-detections
    return X, sigma, bin_mz, times, list(mat.columns)


# ---------------------------------------------------------------------------
# weighted NMF (PMF-lite) via multiplicative updates
# ---------------------------------------------------------------------------
def wnmf(X: np.ndarray, sigma: np.ndarray, k: int, *, n_iter: int = N_ITER,
         seed: int = SEED, tol: float = 1e-6):
    """Minimise sum Wt (X - G H)^2, Wt = 1/sigma^2, G,H >= 0.

    Returns (G [n x k], H [k x m], q [final weighted residual], history).
    Standard weighted multiplicative updates (Wang & Zhang); factors ordered by
    descending explained weighted variance and each H row scaled to sum 1 (a
    fractional composition profile) with the scale pushed into G (factor time
    series in the input's intensity/concentration units)."""
    n, m = X.shape
    Wt = 1.0 / np.maximum(sigma, EPS) ** 2
    WX = Wt * X
    rng = np.random.default_rng(seed)
    scale = np.sqrt(max(np.mean(X), EPS) / k)
    G = np.abs(rng.normal(scale=scale, size=(n, k))) + EPS
    H = np.abs(rng.normal(scale=scale, size=(k, m))) + EPS
    hist = []
    prev = np.inf
    for it in range(n_iter):
        WGH = Wt * (G @ H)
        H *= (G.T @ WX) / (G.T @ WGH + EPS)
        WGH = Wt * (G @ H)
        G *= (WX @ H.T) / (WGH @ H.T + EPS)
        if it % 20 == 0 or it == n_iter - 1:
            q = float(np.sum(Wt * (X - G @ H) ** 2))
            hist.append(q)
            if abs(prev - q) / max(prev, EPS) < tol:
                break
            prev = q
    # order factors by explained weighted variance (descending)
    contrib = [float(np.sum(Wt * (np.outer(G[:, f], H[f]) ** 2))) for f in range(k)]
    order = np.argsort(contrib)[::-1]
    G, H = G[:, order], H[order]
    # scale: H rows sum to 1
    rowsum = H.sum(axis=1, keepdims=True)
    H = H / np.maximum(rowsum, EPS)
    G = G * rowsum.T
    q = float(np.sum(Wt * (X - G @ H) ** 2))
    return G, H, q, hist


def q_expected(X: np.ndarray, sigma: np.ndarray, k: int) -> float:
    """Expected Q for a good fit ~ (#finite data points) - k*(n+m) degrees of
    freedom. Q/Qexp near 1 => the model explains the data to within its errors."""
    n, m = X.shape
    n_data = int(np.isfinite(sigma).sum())
    return max(1.0, n_data - k * (n + m))


def rank_scan(X: np.ndarray, sigma: np.ndarray, ks, *, seed: int = SEED,
              n_iter: int = N_ITER) -> pd.DataFrame:
    """Run wnmf for each k and report Q, Qexp, Q/Qexp. The elbow / where Q/Qexp
    stops dropping steeply is the usual factor-number pick."""
    rows = []
    for k in ks:
        _, _, q, _ = wnmf(X, sigma, k, seed=seed, n_iter=n_iter)
        qe = q_expected(X, sigma, k)
        rows.append(dict(k=int(k), Q=q, Qexp=qe, Q_over_Qexp=q / qe))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# factor interpretation
# ---------------------------------------------------------------------------
def factor_top_ions(H: np.ndarray, bin_mz: pd.Series, cols, *, top: int = 12) -> list[pd.DataFrame]:
    """For each factor, its most DISTINCTIVE bins -- ranked by exclusivity-weighted
    loading `H[f]*exclusivity` where exclusivity = H[f,j]/sum_f H[:,j] (the fraction
    of ion j's total loading that lands on factor f). Ranking by raw H instead just
    resurfaces the same always-on bright background ions in every factor; the
    exclusivity weight surfaces what actually separates the factors. Columns:
    mz, weight (H[f]), exclusivity (0..1), score."""
    colsum = H.sum(axis=0) + EPS                      # each ion's total loading
    mz = bin_mz.reindex(cols).values
    out = []
    for f in range(H.shape[0]):
        excl = H[f] / colsum
        score = H[f] * excl
        d = pd.DataFrame({"mz": mz, "weight": H[f], "exclusivity": excl, "score": score}, index=cols)
        out.append(d.sort_values("score", ascending=False).head(top).reset_index(drop=True))
    return out


def factor_timeseries(G: np.ndarray, times: pd.DatetimeIndex) -> pd.DataFrame:
    """Factor time series as a tidy frame [datetime_utc, f0, f1, ...]."""
    df = pd.DataFrame(G, columns=[f"f{f}" for f in range(G.shape[1])])
    df.insert(0, "datetime_utc", pd.DatetimeIndex(times).values)
    return df


def run_nmf(peaks: pd.DataFrame, k: int, *, tol_ppm: float = DEFAULT_TOL_PPM,
            reagent_series: pd.Series | None = None, seed: int = SEED,
            n_iter: int = N_ITER):
    """End-to-end on a peaks frame: build input -> factorise -> package.
    Returns a dict with G, H, bin_mz, times, cols, q, q_over_qexp."""
    X, sigma, bin_mz, times, cols = build_input(
        peaks, tol_ppm=tol_ppm, reagent_series=reagent_series)
    G, H, q, hist = wnmf(X, sigma, k, seed=seed, n_iter=n_iter)
    return dict(G=G, H=H, bin_mz=bin_mz, times=times, cols=cols, q=q,
                q_over_qexp=q / q_expected(X, sigma, k), X=X, sigma=sigma,
                ts=factor_timeseries(G, times),
                top_ions=factor_top_ions(H, bin_mz, cols))

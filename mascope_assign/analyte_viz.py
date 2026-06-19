"""Consistent analyte visualisation: Van Krevelen + raw time-series.

Both plots are derived here from one committed ledger (the assignment) plus an
optional batch time-series peak table, so the negative Br-CIMS and positive
urea-CIMS figures are computed *identically* -- same analyte definition, same
H/C-O/C convention, same "changing" threshold, same RAW intensity choice.

Design choices (deliberate, so the figures are comparable):
- ANALYTE = a committed M0 with an organic neutral, EXCLUDING contamination
  (Si-bearing siloxane/silanediol) and the reagent/artifact roles. Polarity- and
  reagent-agnostic: it reads the neutral formula, not the adduct.
- Van Krevelen on the NEUTRAL: x = O/C, y = H/C (halogens are not in these
  neutrals -- the reagent halogen lives in the adduct).
- A compound's time trace is the RAW summed intensity of its bins under the
  given adducts. RAW (not reagent-normalised): a positive urea-CIMS spectrum
  often excludes the reagent ions from its mass range, and TIC is
  analyte-dominated, so reagent/TIC normalisation introduces closure artifacts.
  `mode='reagent'/'tic'` is available when a real normaliser exists.
- `cv` = std/mean of the (mode) trace; `changing` = cv >= CHANGING_CV. Same
  threshold both instruments.

Pure data functions (no plotting deps) + lazy-matplotlib renderers + a
`widget_payload()` that returns the arrays the interactive chat widget consumes.
"""
from __future__ import annotations

import bisect

import numpy as np
import pandas as pd

from . import chemistry as C
from . import timeseries as TS

__version__ = "0.1.0"

CHANGING_CV = 0.30           # cv_raw at/above which a peak is "changing"
DEFAULT_BIN_MIN = 30         # time-bin width (minutes) for the trace plot
_CONTAM_ELEMENTS = ("Si",)   # siloxane / silanediol = instrument contamination


def _is_contaminant(formula: str) -> bool:
    cnt = C.parse_formula(str(formula))
    return any(cnt.get(e, 0) for e in _CONTAM_ELEMENTS)


def analyte_table(ledger: pd.DataFrame, *, exclude_contaminant: bool = True
                  ) -> pd.DataFrame:
    """The committed organic analytes (one row per M0). Adds H/C, O/C, N count
    and a CHO/CHON/CHOS class. Contamination (Si) and non-M0 roles dropped."""
    m0 = ledger[ledger["role"] == "M0"].copy()
    m0 = m0[m0["neutral_formula"].notna() & (m0["neutral_formula"].astype(str) != "")]
    if exclude_contaminant:
        m0 = m0[~m0["neutral_formula"].map(_is_contaminant)]
    cnt = m0["neutral_formula"].map(lambda f: C.parse_formula(str(f)))
    for col, el in (("nC", "C"), ("nH", "H"), ("nO", "O"), ("nN", "N"), ("nS", "S")):
        m0[col] = cnt.map(lambda c, el=el: c.get(el, 0))
    m0 = m0[m0["nC"] > 0].copy()
    m0["hc"] = m0["nH"] / m0["nC"]
    m0["oc"] = m0["nO"] / m0["nC"]
    m0["klass"] = np.where(m0["nS"] > 0, "CHOS",
                           np.where(m0["nN"] > 0, "CHON", "CHO"))
    # channels (adducts) a neutral was seen in -- the "ion" half of the hover.
    if "adduct" in m0.columns:
        chmap = (m0.groupby("neutral_formula")["adduct"]
                 .apply(lambda s: " / ".join(dict.fromkeys(s.dropna().astype(str)))).to_dict())
        m0["channels"] = m0["neutral_formula"].map(chmap)
    else:
        m0["channels"] = ""
    # ONE row per NEUTRAL: a compound seen in two channels ([M+Br]- and [M-H]-,
    # or [M+H]+ and [M+urea+H]+) is a single analyte. Keep the best tier.
    if "tier" in m0.columns:
        m0["_tr"] = m0["tier"].map({"Identified": 0, "Candidate": 1}).fillna(2)
        m0 = m0.sort_values("_tr").drop_duplicates("neutral_formula", keep="first").drop(columns="_tr")
    else:
        m0 = m0.drop_duplicates("neutral_formula", keep="first")
    return m0.reset_index(drop=True)


# ---------------------------------------------------------------------------
# time series
# ---------------------------------------------------------------------------
def _bins_for_formula(arr, idx, cols_set, formula, adducts):
    """Bin indices whose m/z match the neutral under any adduct (<=8 ppm)."""
    out = []
    for a in adducts:
        if a not in C.ADDUCT_SHIFTS:
            continue
        mz = C.ion_mz(formula, a)
        i = bisect.bisect_left(arr, mz)
        for j in (i - 1, i):
            if 0 <= j < len(arr) and abs(arr[j] - mz) / mz * 1e6 <= 8:
                out.append(idx[j])
    return [c for c in set(out) if c in cols_set]


def time_traces(ts_peaks: pd.DataFrame, formulas, adducts, *, mode: str = "raw",
                bin_minutes: int = DEFAULT_BIN_MIN, reagent_mzs=None
                ) -> tuple[np.ndarray, pd.DataFrame]:
    """Raw (or normalised) intensity trace per formula, binned on wall-clock.

    Returns (hours grid, DataFrame indexed by grid-bin, columns = formulas) of
    the binned MEDIAN intensity. NaN where a formula has no matched bin."""
    ts_peaks = ts_peaks.copy()
    ts_peaks["datetime_utc"] = pd.to_datetime(ts_peaks["datetime_utc"], utc=True)
    tstamp = ts_peaks.groupby("sample_item_id")["datetime_utc"].first()
    mat, bin_mz = TS.build_matrix(ts_peaks)
    mat = mat.reindex(tstamp.sort_values().index)
    t0 = tstamp.min()
    hr = (tstamp.reindex(mat.index) - t0).dt.total_seconds().values / 3600.0
    bm = bin_mz.sort_values(); arr = bm.values; idx = bm.index.values
    cols_set = set(mat.columns)
    if mode == "reagent" and reagent_mzs:
        rt = TS.reagent_total(mat, bin_mz, reagent_mzs)
    elif mode == "tic":
        rt = mat.sum(axis=1)
    else:
        rt = None
    grid = np.arange(0, np.nanmax(hr) + bin_minutes / 60.0, bin_minutes / 60.0)
    digit = np.digitize(hr, grid)
    out = {}
    for f in formulas:
        b = _bins_for_formula(arr, idx, cols_set, f, adducts)
        if not b:
            out[f] = pd.Series(np.nan, index=range(1, len(grid) + 1))
            continue
        s = mat[b].sum(axis=1)
        if rt is not None:
            s = s / rt.replace(0, np.nan)
        out[f] = pd.Series(s.values).groupby(digit).median().reindex(range(1, len(grid) + 1))
    return grid, pd.DataFrame(out)


def attach_dynamics(analytes: pd.DataFrame, ts_peaks: pd.DataFrame, adducts, *,
                    mode: str = "raw", reagent_mzs=None) -> pd.DataFrame:
    """Add median_cps / cv / changing to the analyte table from the time series."""
    df = analytes.copy()
    _, traces = time_traces(ts_peaks, df["neutral_formula"].tolist(), adducts,
                            mode=mode, reagent_mzs=reagent_mzs)
    med, cv = {}, {}
    for f in df["neutral_formula"]:
        tr = traces[f].dropna().values if f in traces else np.array([])
        med[f] = float(np.median(tr)) if len(tr) else 0.0
        cv[f] = float(np.std(tr) / np.mean(tr)) if len(tr) and np.mean(tr) > 0 else np.nan
    df["median_cps"] = df["neutral_formula"].map(med)
    df["cv"] = df["neutral_formula"].map(cv)
    df["changing"] = df["cv"] >= CHANGING_CV
    return df


# ---------------------------------------------------------------------------
# widget payload (the arrays the interactive chat plot consumes)
# ---------------------------------------------------------------------------
def widget_payload(analytes: pd.DataFrame, grid: np.ndarray | None = None,
                   traces: pd.DataFrame | None = None, *, top_ts: int = 28) -> dict:
    """{'vk': [[oc,hc,n,changing,logI,formula,channels,tier],...],
        'ts': {'grid':[h...], 'series':[{f,n,ch,y:[...]}, ...]}}. `channels` is
        the adduct(s) the neutral was detected as -- the 'ion' shown on hover."""
    vk = []
    for r in analytes.itertuples():
        logI = round(float(np.log10(max(getattr(r, "median_cps", 0.0), 1))), 2)
        vk.append([round(r.oc, 3), round(r.hc, 3), int(r.nN > 0),
                   int(bool(getattr(r, "changing", False))), logI, r.neutral_formula,
                   str(getattr(r, "channels", "")), str(getattr(r, "tier", ""))])
    ts = None
    if grid is not None and traces is not None and "changing" in analytes:
        chg = analytes[analytes["changing"]].sort_values("median_cps", ascending=False).head(top_ts)
        series = []
        for r in chg.itertuples():
            y = [None if pd.isna(v) or v <= 0 else round(float(v), 0)
                 for v in traces[r.neutral_formula].values]
            series.append({"f": r.neutral_formula, "n": int(r.nN > 0),
                           "ch": str(getattr(r, "channels", "")), "y": y})
        ts = {"grid": [round(float(x), 1) for x in grid], "series": series}
    return {"vk": vk, "ts": ts}


# ---------------------------------------------------------------------------
# matplotlib renderers (lazy import; standalone PNGs)
# ---------------------------------------------------------------------------
def render_van_krevelen(analytes: pd.DataFrame, path: str, *, title: str = "") -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    sz = (np.maximum(2.0, (np.log10(analytes["median_cps"].clip(lower=1)) - 2.5) * 6))
    flat = analytes[~analytes.get("changing", pd.Series(False, index=analytes.index))]
    chg = analytes[analytes.get("changing", pd.Series(False, index=analytes.index))]
    ax.scatter(flat["oc"], flat["hc"], s=sz.loc[flat.index], c="#B4B2A9", alpha=0.4,
               linewidths=0, label="flat background")
    for kl, col in (("CHO", "#1D9E75"), ("CHON", "#7F77DD"), ("CHOS", "#D85A30")):
        g = chg[chg["klass"] == kl]
        if len(g):
            ax.scatter(g["oc"], g["hc"], s=sz.loc[g.index], c=col, alpha=0.85,
                       linewidths=0.3, edgecolors="white", label=f"changing — {kl}")
    ax.set_xlabel("O/C"); ax.set_ylabel("H/C"); ax.set_xlim(0, 1.3); ax.set_ylim(0.3, 3.0)
    ax.set_title(title); ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)
    return path


def render_timeseries(grid: np.ndarray, traces: pd.DataFrame, analytes: pd.DataFrame,
                      path: str, *, top: int = 28, title: str = "") -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nN = dict(zip(analytes["neutral_formula"], analytes["nN"]))
    chg = analytes[analytes.get("changing", pd.Series(False, index=analytes.index))]
    chg = chg.sort_values("median_cps", ascending=False).head(top)
    fig, ax = plt.subplots(figsize=(8, 5))
    for f in chg["neutral_formula"]:
        col = "#7F77DD" if nN.get(f, 0) else "#1D9E75"
        ax.plot(grid, traces[f].values, color=col, lw=1.1, alpha=0.7)
    ax.set_yscale("log"); ax.set_xlabel("hour of experiment (UTC)")
    ax.set_ylabel("raw intensity (cps)"); ax.set_title(title); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)
    return path

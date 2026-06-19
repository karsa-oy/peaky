"""Time-resolved disposition layer.

A sum spectrum cannot tell a bright stable INLET/instrument contaminant from a
real ambient ANALYTE, nor a degenerate reagent-background cluster from chemistry.
A *time series* of the same sample over a day can: in a halide-CIMS the physically
meaningful quantity is the analyte normalised to the reagent ion (removes
instrument-sensitivity + reagent-flow common-mode drift), and then

  * a FLAT normalised trace (low coefficient of variation, no diel) == inlet /
    instrument background or a constant reagent cluster -- NOT ambient chemistry;
  * a VARIABLE trace that co-varies with a known chemical family == real ambient
    analyte.

This module ingests a batch's per-sample peak table, builds the reagent-normalised
intensity matrix, measures each peak's variability (`cv_norm`) and (optionally)
its correlation to reference family traces, and stamps a `ts_*` disposition onto
the ledger. It then applies CONSERVATIVE auto-actions: demote a flat di-bromide /
background-channel commit (TS-confirmed background) and flag inlet contaminants.
It never changes a formula -- only the tier/role annotation, with commentary.

All pure pandas/numpy; no network. Reference (2026-06-16 time-series unlock).
"""
from __future__ import annotations

import bisect
import re

import numpy as np
import pandas as pd

from . import ledger as L

__version__ = "0.1.0"

DEFAULT_TOL_PPM = 5.0
FLAT_CV = 0.25          # cv_norm below this == flat / background
COVARY_R = 0.70         # correlation above this == co-varies with the family


# ---------------------------------------------------------------------------
# matrix construction
# ---------------------------------------------------------------------------
def build_matrix(peaks: pd.DataFrame, *, tol_ppm: float = DEFAULT_TOL_PPM,
                 mz_col="mz", height_col="height", sample_col="sample_item_id"
                 ) -> tuple[pd.DataFrame, pd.Series]:
    """Gap-cluster peaks into m/z bins (ppm tolerance) and pivot to a
    samples x bin intensity matrix. Returns (matrix, bin_mz)."""
    d = peaks[[sample_col, mz_col, height_col]].dropna().sort_values(mz_col).reset_index(drop=True)
    mz = d[mz_col].to_numpy()
    if len(mz) == 0:
        return pd.DataFrame(), pd.Series(dtype=float)
    gaps = np.diff(mz) / mz[:-1] * 1e6
    binid = np.zeros(len(mz), dtype=np.int64)
    binid[1:] = np.cumsum(gaps > tol_ppm)
    d["_bin"] = binid
    wsum = (d[mz_col] * d[height_col]).groupby(d["_bin"]).sum()
    hsum = d[height_col].groupby(d["_bin"]).sum()
    bin_mz = (wsum / hsum).rename("mz")
    mat = d.pivot_table(index=sample_col, columns="_bin", values=height_col, aggfunc="sum")
    return mat, bin_mz


def reagent_total(mat: pd.DataFrame, bin_mz: pd.Series, reagent_mzs, *, tol_ppm=8.0):
    """Per-sample sum of the reagent bins (the normaliser). reagent_mzs is a list
    of reagent ion m/z (e.g. the Br3- isotopologues)."""
    cols = []
    bm = bin_mz.sort_values()
    arr = bm.to_numpy(); idx = bm.index.to_numpy()
    for r in reagent_mzs:
        i = bisect.bisect_left(arr, r)
        for j in (i - 1, i):
            if 0 <= j < len(arr) and abs(arr[j] - r) / r * 1e6 <= tol_ppm:
                cols.append(idx[j])
    cols = [c for c in set(cols) if c in mat.columns]
    if not cols:
        return None
    return mat[cols].sum(axis=1)


def normalize(mat: pd.DataFrame, reagent_series) -> pd.DataFrame:
    """Divide every bin by the per-sample reagent total (concentration proxy)."""
    if reagent_series is None:
        return mat
    return mat.div(reagent_series.replace(0, np.nan), axis=0)


def bin_metrics(norm: pd.DataFrame, bin_mz: pd.Series) -> pd.DataFrame:
    """Per-bin presence + cv_norm on the (reagent-normalised) matrix."""
    n = len(norm)
    presence = norm.notna().sum() / n if n else norm.notna().sum()
    mean = norm.mean(); std = norm.std()
    cv = (std / mean).replace([np.inf, -np.inf], np.nan)
    out = pd.DataFrame({"mz": bin_mz.reindex(norm.columns), "presence": presence,
                        "median": norm.median(), "cv_norm": cv})
    out.index.name = "_bin"
    return out


def family_trace(norm: pd.DataFrame, bin_ids):
    """z-scored mean log-trace of a set of bins (a reference family trace)."""
    bb = [b for b in bin_ids if b in norm.columns]
    if not bb:
        return None
    lg = np.log10(norm[bb].clip(lower=norm[norm > 0].min().min() or 1e-9))
    z = (lg - lg.mean()) / lg.std()
    return z.mean(axis=1)


def correlate(norm: pd.DataFrame, trace) -> pd.Series:
    if trace is None:
        return pd.Series(np.nan, index=norm.columns)
    lg = np.log10(norm.clip(lower=norm[norm > 0].min().min() or 1e-9))
    with np.errstate(invalid="ignore", divide="ignore"):  # flat bins -> NaN r (fine)
        return lg.apply(lambda c: c.corr(trace))


# ---------------------------------------------------------------------------
# disposition + ledger application
# ---------------------------------------------------------------------------
def _nearest(bin_mz_sorted_vals, bin_mz_sorted_idx, mz, tol_ppm):
    i = bisect.bisect_left(bin_mz_sorted_vals, mz)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(bin_mz_sorted_vals):
            ppm = abs(bin_mz_sorted_vals[j] - mz) / mz * 1e6
            if ppm <= tol_ppm and (best is None or ppm < best[1]):
                best = (bin_mz_sorted_idx[j], ppm)
    return best[0] if best else None


def _disposition(row, cv, r_mono, r_formic):
    """Classify one M0 row from its formula + time-series behavior."""
    ion = str(row.get("ion_formula", "")); adduct = str(row.get("adduct", ""))
    neutral = str(row.get("neutral_formula", ""))
    if "CO3" in adduct:
        return "background:CO3-channel (TS-flat)" if (pd.notna(cv) and cv < FLAT_CV) else "CO3-channel"
    if re.search(r"Br[23]", ion):
        return "background:di-bromide cluster (TS-flat)" if (pd.notna(cv) and cv < FLAT_CV) else "di-bromide"
    if pd.notna(cv) and cv < FLAT_CV:
        if "Si" in neutral or "F" in ion:
            return "background:inlet/instrument contaminant (TS-flat)"
        return "background:flat (TS-flat)"
    if pd.notna(r_mono) and r_mono >= COVARY_R:
        return "ambient:biogenic-SOA (co-varies)"
    if pd.notna(r_formic) and r_formic >= 0.9:
        return "ambient:acid/oxygenate pool (co-varies)"
    if pd.notna(cv) and cv >= 0.45:
        return "ambient:variable"
    return "intermediate"


def apply_timeseries(ledger: pd.DataFrame, peaks: pd.DataFrame, *,
                     reagent_mzs=None, mono_anchor_mzs=None, formic_mz=None,
                     tol_ppm: float = DEFAULT_TOL_PPM, demote=True, log=print) -> dict:
    """Annotate `ledger` (in place) with ts_cv_norm / ts_r_mono / ts_r_formic /
    ts_disposition from the time-series `peaks` table, and (if demote) cap a flat
    di-bromide / CO3-channel Identified commit at Candidate (TS-confirmed
    background). Returns a summary dict. Reagent normaliser + anchors are taken
    from the ledger when not supplied.
    """
    summary = {"annotated": 0, "demoted": 0, "ambient": 0, "background": 0}
    for col in ("ts_cv_norm", "ts_r_mono", "ts_r_formic", "ts_disposition"):
        if col not in ledger.columns:
            ledger[col] = np.nan if col != "ts_disposition" else ""

    mat, bin_mz = build_matrix(peaks, tol_ppm=tol_ppm)
    if mat.empty:
        log("[timeseries] no peaks; skipped"); return summary

    # reagent normaliser: explicit, else the ledger's reagent Br_n rows
    if reagent_mzs is None:
        rr = ledger[(ledger["role"] == L.ROLE_REAGENT)
                    & ledger["ion_formula"].astype(str).str.match(r"Br\d-")]
        reagent_mzs = rr["mz"].dropna().tolist()
    rt = reagent_total(mat, bin_mz, reagent_mzs) if reagent_mzs else None
    norm = normalize(mat, rt)
    met = bin_metrics(norm, bin_mz)

    # reference family traces (optional)
    bmz_s = bin_mz.sort_values(); bvals = bmz_s.to_numpy(); bidx = bmz_s.index.to_numpy()
    def bins_for(mzs):
        out = []
        for m in (mzs or []):
            b = _nearest(bvals, bidx, m, tol_ppm)
            if b is not None:
                out.append(b)
        return out
    if mono_anchor_mzs is None:
        mono_anchor_mzs = ledger.loc[
            ledger["neutral_formula"].astype(str).isin(
                {"C10H16O3", "C10H16O4", "C10H16O5", "C10H16O6"}), "mz"].dropna().tolist()
    mono_tr = family_trace(norm, bins_for(mono_anchor_mzs))
    formic_b = _nearest(bvals, bidx, formic_mz, tol_ppm) if formic_mz else \
        _nearest(bvals, bidx, 124.9243, tol_ppm)
    formic_tr = norm[formic_b].pipe(lambda c: np.log10(c.clip(lower=1e-9))) if formic_b in norm.columns else None
    r_mono = correlate(norm, mono_tr)
    r_formic = correlate(norm, formic_tr)

    # stamp the ledger (M0 rows)
    for i in ledger.index[ledger["role"] == L.ROLE_M0]:
        mz = ledger.at[i, "mz"]
        if pd.isna(mz):
            continue
        b = _nearest(bvals, bidx, float(mz), tol_ppm)
        cv = float(met.at[b, "cv_norm"]) if (b is not None and b in met.index and pd.notna(met.at[b, "cv_norm"])) else np.nan
        rm = float(r_mono.get(b, np.nan)) if b is not None else np.nan
        rf = float(r_formic.get(b, np.nan)) if b is not None else np.nan
        disp = _disposition(ledger.loc[i], cv, rm, rf)
        ledger.at[i, "ts_cv_norm"] = cv
        ledger.at[i, "ts_r_mono"] = rm
        ledger.at[i, "ts_r_formic"] = rf
        ledger.at[i, "ts_disposition"] = disp
        summary["annotated"] += 1
        if disp.startswith("ambient"):
            summary["ambient"] += 1
        elif disp.startswith("background"):
            summary["background"] += 1
            # conservative auto-demote: a flat di-bromide / CO3 background commit
            # must not stay Identified once the time series shows it is background
            if demote and str(ledger.at[i, "tier"]) == "Identified" and (
                    "di-bromide" in disp or "CO3-channel" in disp):
                ledger.at[i, "tier"] = "Candidate"
                ledger.at[i, "tier_reason"] = (str(ledger.at[i, "tier_reason"] or "")
                    + " | time-series: flat background (reagent/inlet), demoted").strip(" |")
                summary["demoted"] += 1
    log(f"[timeseries] {summary}")
    return summary

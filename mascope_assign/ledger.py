"""The peak ledger -- the single mutable state carrier for the whole pipeline.

One ledger row == one PHYSICAL peak (deduplicated by peak_id). Passes fill and
annotate rows; they never drop them. Every committed assignment carries
provenance + commentary, and the commit API enforces the structural invariants
so no pass can corrupt the shared state:

  I1. Each peak has exactly one role: 'unexplained' | 'M0' | 'iso_child'
      | 'reagent'.
  I2. An 'iso_child' row points (parent_peak_id) to a peak that owns an M0
      assignment.
  I3. A peak that is already locked is immutable to later passes.
  I4. No peak is claimed twice (an iso_child cannot also be an M0 owner).
  I5. Every M0 / series / contaminant assignment records pass + method +
      confidence + commentary.

This module is pure pandas; it never talks to Mascope.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import chemistry as C

__version__ = "0.2.0"

ROLE_UNEXPLAINED = "unexplained"
ROLE_M0 = "M0"
ROLE_ISO = "iso_child"
ROLE_REAGENT = "reagent"

# Canonical column set. Identity columns come from the peak source; the rest are
# filled by the assignment passes.
_ASSIGN_COLS: dict[str, object] = {
    "neutral_formula": pd.NA,
    "adduct": pd.NA,
    "ion_formula": pd.NA,
    "ion_score": np.nan,
    "compound_score": np.nan,
    "ppm_error": np.nan,
    "dbe": np.nan,
    "confidence": pd.NA,
    "role": ROLE_UNEXPLAINED,
    "parent_peak_id": pd.NA,
    "iso_label": pd.NA,
    "iso_match_score": np.nan,
    "pass_no": pd.NA,
    "method": pd.NA,
    "anchor_peak_id": pd.NA,
    "series_unit": pd.NA,
    "locked": False,
    "commentary": pd.NA,
    "alternatives": pd.NA,     # JSON string
    "isotopologues": pd.NA,    # JSON string (per-isotopologue Mascope scores)
}

_REQUIRED_IDENTITY = ("peak_id", "mz")


class LedgerError(Exception):
    pass


def new_ledger(peaks: pd.DataFrame) -> pd.DataFrame:
    """Build a fresh ledger from a peaks table (must have peak_id, mz; height/
    area optional). Deduplicates by peak_id keeping the highest-intensity row."""
    for col in _REQUIRED_IDENTITY:
        if col not in peaks.columns:
            raise LedgerError(f"peaks table missing required column {col!r}")
    df = peaks.copy()
    sort_col = "height" if "height" in df.columns else ("area" if "area" in df.columns else None)
    if sort_col is not None:
        df = df.sort_values(sort_col, ascending=False)
    df = df.drop_duplicates(subset="peak_id", keep="first").reset_index(drop=True)
    if "height" not in df.columns:
        df["height"] = np.nan
    if "area" not in df.columns:
        df["area"] = np.nan
    for col, default in _ASSIGN_COLS.items():
        df[col] = default
    return df


def _row_index(ledger: pd.DataFrame, peak_id) -> int:
    idx = ledger.index[ledger["peak_id"] == peak_id]
    if len(idx) == 0:
        raise LedgerError(f"peak_id {peak_id!r} not in ledger")
    if len(idx) > 1:
        raise LedgerError(f"peak_id {peak_id!r} is duplicated in ledger")
    return int(idx[0])


def is_locked(ledger: pd.DataFrame, peak_id) -> bool:
    return bool(ledger.at[_row_index(ledger, peak_id), "locked"])


def role_of(ledger: pd.DataFrame, peak_id) -> str:
    return str(ledger.at[_row_index(ledger, peak_id), "role"])


def unassigned_peaks(ledger: pd.DataFrame) -> pd.DataFrame:
    return ledger[ledger["role"] == ROLE_UNEXPLAINED]


def assigned_peaks(ledger: pd.DataFrame) -> pd.DataFrame:
    return ledger[ledger["role"] == ROLE_M0]


def commit_assignment(
    ledger: pd.DataFrame,
    peak_id,
    *,
    neutral_formula: str,
    adduct: str,
    ion_formula: str | None = None,
    ion_score: float,
    compound_score: float | None = None,
    ppm_error: float | None = None,
    pass_no: int,
    method: str,
    confidence: str,
    commentary: str,
    alternatives: list[dict] | None = None,
    isotopologues: list[dict] | None = None,
    anchor_peak_id=None,
    series_unit: str | None = None,
    overwrite: bool = False,
) -> pd.DataFrame:
    """Commit an M+0 (monoisotopic) owner to a peak. Enforces I1-I5.

    Returns the ledger (mutated in place and returned for chaining)."""
    i = _row_index(ledger, peak_id)
    cur_role = str(ledger.at[i, "role"])

    # I3: locked rows are immutable
    if bool(ledger.at[i, "locked"]) and not overwrite:
        raise LedgerError(f"peak {peak_id!r} is locked; refusing to overwrite")
    # I4: an iso_child of a locked parent cannot be reassigned as M0
    if cur_role == ROLE_ISO:
        parent = ledger.at[i, "parent_peak_id"]
        if parent is not pd.NA and not pd.isna(parent) and is_locked(ledger, parent) and not overwrite:
            raise LedgerError(
                f"peak {peak_id!r} is the isotopologue child of locked parent "
                f"{parent!r}; refusing to claim it as an M0 owner")
    # I5: provenance required
    if not commentary or not method or not confidence:
        raise LedgerError("commit requires non-empty commentary, method, confidence")

    ledger.at[i, "neutral_formula"] = neutral_formula
    ledger.at[i, "adduct"] = adduct
    ledger.at[i, "ion_formula"] = ion_formula if ion_formula is not None else neutral_formula
    ledger.at[i, "ion_score"] = float(ion_score)
    ledger.at[i, "compound_score"] = (np.nan if compound_score is None else float(compound_score))
    ledger.at[i, "ppm_error"] = (np.nan if ppm_error is None else float(ppm_error))
    ledger.at[i, "dbe"] = C.dbe(neutral_formula)
    ledger.at[i, "confidence"] = confidence
    ledger.at[i, "role"] = ROLE_M0
    ledger.at[i, "parent_peak_id"] = pd.NA
    ledger.at[i, "iso_label"] = pd.NA
    ledger.at[i, "iso_match_score"] = np.nan
    ledger.at[i, "pass_no"] = int(pass_no)
    ledger.at[i, "method"] = method
    ledger.at[i, "anchor_peak_id"] = pd.NA if anchor_peak_id is None else anchor_peak_id
    ledger.at[i, "series_unit"] = pd.NA if series_unit is None else series_unit
    ledger.at[i, "commentary"] = commentary
    ledger.at[i, "alternatives"] = json.dumps(alternatives or [])
    ledger.at[i, "isotopologues"] = json.dumps(isotopologues or [])
    return ledger


def attach_isotopologue(
    ledger: pd.DataFrame,
    child_peak_id,
    parent_peak_id,
    *,
    iso_label: str,
    iso_match_score: float | None = None,
    overwrite: bool = False,
) -> pd.DataFrame:
    """Mark child_peak_id as an isotopologue of an M0-owning parent. Enforces
    I2 (parent must own an M0 assignment) and I3/I4 (child must be free)."""
    pi = _row_index(ledger, parent_peak_id)
    if str(ledger.at[pi, "role"]) != ROLE_M0:
        raise LedgerError(f"parent {parent_peak_id!r} does not own an M0 assignment")
    ci = _row_index(ledger, child_peak_id)
    if child_peak_id == parent_peak_id:
        raise LedgerError("a peak cannot be its own isotopologue child")
    cur = str(ledger.at[ci, "role"])
    if cur == ROLE_M0 and not overwrite:
        raise LedgerError(f"peak {child_peak_id!r} owns an M0 assignment; not overwriting with iso role")
    if bool(ledger.at[ci, "locked"]) and not overwrite:
        raise LedgerError(f"peak {child_peak_id!r} is locked")
    ledger.at[ci, "role"] = ROLE_ISO
    ledger.at[ci, "parent_peak_id"] = parent_peak_id
    ledger.at[ci, "iso_label"] = iso_label
    ledger.at[ci, "iso_match_score"] = (np.nan if iso_match_score is None else float(iso_match_score))
    return ledger


_ASSIGNMENT_FIELDS = {
    "neutral_formula": pd.NA, "adduct": pd.NA, "ion_formula": pd.NA,
    "ion_score": np.nan, "compound_score": np.nan, "ppm_error": np.nan,
    "dbe": np.nan, "confidence": pd.NA, "pass_no": np.nan, "method": pd.NA,
    "anchor_peak_id": pd.NA, "series_unit": pd.NA,
    "alternatives": pd.NA, "isotopologues": pd.NA,
}


def clear_assignment(ledger: pd.DataFrame, peak_id, *, reason: str) -> pd.DataFrame:
    """Demote an M0 owner back to unexplained (e.g. failed a post-hoc audit).

    Its isotopologue children are orphaned by I2, so they are cleared back to
    unexplained too. Locked peaks are refused. `reason` is recorded in the
    commentary so the demotion is auditable."""
    i = _row_index(ledger, peak_id)
    if bool(ledger.at[i, "locked"]):
        raise LedgerError(f"peak {peak_id!r} is locked; refusing to clear")
    if str(ledger.at[i, "role"]) != ROLE_M0:
        raise LedgerError(f"peak {peak_id!r} does not own an M0 assignment")
    for ci in ledger.index[ledger["parent_peak_id"] == peak_id]:
        ledger.at[ci, "role"] = ROLE_UNEXPLAINED
        ledger.at[ci, "parent_peak_id"] = pd.NA
        ledger.at[ci, "iso_label"] = pd.NA
        ledger.at[ci, "iso_match_score"] = np.nan
    old = ledger.at[i, "commentary"]
    for col, na in _ASSIGNMENT_FIELDS.items():
        ledger.at[i, col] = na
    ledger.at[i, "role"] = ROLE_UNEXPLAINED
    ledger.at[i, "commentary"] = (f"CLEARED ({reason}). Was: {old}"
                                  if old is not pd.NA and pd.notna(old) else
                                  f"CLEARED ({reason}).")
    return ledger


def displace_to_isotopologue(
    ledger: pd.DataFrame,
    child_peak_id,
    parent_peak_id,
    *,
    iso_label: str,
    iso_match_score: float | None = None,
) -> pd.DataFrame:
    """Convert a peak that owns an M0 assignment into the isotopologue child of
    a stronger parent (M0-vs-iso-child arbitration). The child's own former
    iso children are re-parented to the new parent with combined labels: a 13C
    satellite of a displaced 81Br twin is the parent's 13C+81Br line."""
    ci = _row_index(ledger, child_peak_id)
    if bool(ledger.at[ci, "locked"]):
        raise LedgerError(f"peak {child_peak_id!r} is locked; refusing to displace")
    if str(ledger.at[ci, "role"]) != ROLE_M0:
        raise LedgerError(f"peak {child_peak_id!r} does not own an M0 assignment")
    old = ledger.at[ci, "commentary"]
    for gi in ledger.index[ledger["parent_peak_id"] == child_peak_id]:
        ledger.at[gi, "parent_peak_id"] = parent_peak_id
        glab = str(ledger.at[gi, "iso_label"])
        ledger.at[gi, "iso_label"] = "+".join(sorted(set(
            glab.split("+") + iso_label.split("+"))))
    for col, na in _ASSIGNMENT_FIELDS.items():
        ledger.at[ci, col] = na
    attach_isotopologue(ledger, child_peak_id, parent_peak_id,
                        iso_label=iso_label, iso_match_score=iso_match_score,
                        overwrite=True)
    ledger.at[ci, "commentary"] = (
        f"DISPLACED to {iso_label} isotopologue of {parent_peak_id!r}. Was: {old}"
        if old is not pd.NA and pd.notna(old) else
        f"DISPLACED to {iso_label} isotopologue of {parent_peak_id!r}.")
    return ledger


def mark_reagent(ledger: pd.DataFrame, peak_id, label: str) -> pd.DataFrame:
    i = _row_index(ledger, peak_id)
    if bool(ledger.at[i, "locked"]):
        raise LedgerError(f"peak {peak_id!r} is locked")
    ledger.at[i, "role"] = ROLE_REAGENT
    ledger.at[i, "commentary"] = label
    return ledger


def lock_peaks(ledger: pd.DataFrame, peak_ids) -> pd.DataFrame:
    """Make these peaks immutable to later passes."""
    mask = ledger["peak_id"].isin(list(peak_ids))
    ledger.loc[mask, "locked"] = True
    return ledger


def validate(ledger: pd.DataFrame) -> list[str]:
    """Return a list of invariant violations (empty == healthy)."""
    problems: list[str] = []
    # duplicate peak ids
    dups = ledger["peak_id"][ledger["peak_id"].duplicated()].tolist()
    if dups:
        problems.append(f"duplicate peak_id rows: {dups[:5]}")
    # role domain
    bad_roles = set(ledger["role"]) - {ROLE_UNEXPLAINED, ROLE_M0, ROLE_ISO, ROLE_REAGENT}
    if bad_roles:
        problems.append(f"unknown roles: {bad_roles}")
    # I2: every iso_child points to an M0 owner
    m0_ids = set(ledger.loc[ledger["role"] == ROLE_M0, "peak_id"])
    iso = ledger[ledger["role"] == ROLE_ISO]
    for _, r in iso.iterrows():
        p = r["parent_peak_id"]
        if p is pd.NA or pd.isna(p):
            problems.append(f"iso_child {r['peak_id']!r} has no parent")
        elif p not in m0_ids:
            problems.append(f"iso_child {r['peak_id']!r} parent {p!r} is not an M0 owner")
    # I5: M0 rows carry provenance
    m0 = ledger[ledger["role"] == ROLE_M0]
    miss = m0[m0["commentary"].isna() | m0["method"].isna() | m0["confidence"].isna()]
    if len(miss):
        problems.append(f"{len(miss)} M0 rows missing provenance")
    return problems


def stats(ledger: pd.DataFrame) -> dict:
    """Coverage summary by count and by signal (height)."""
    n = len(ledger)
    h_total = ledger["height"].sum(skipna=True)
    out = {"n_peaks": n, "by_role": {}, "signal_by_role": {},
           "count_frac_by_role": {}}
    for role in (ROLE_M0, ROLE_ISO, ROLE_REAGENT, ROLE_UNEXPLAINED):
        sub = ledger[ledger["role"] == role]
        out["by_role"][role] = int(len(sub))
        out["signal_by_role"][role] = (
            float(sub["height"].sum(skipna=True) / h_total) if h_total else 0.0)
        # signal-% is dominated by the reagent ions in CIMS; the peak-count
        # fraction is the honest residual metric, so report both
        out["count_frac_by_role"][role] = float(len(sub) / n) if n else 0.0
    if "confidence" in ledger.columns:
        out["by_confidence"] = (
            ledger.loc[ledger["role"] == ROLE_M0, "confidence"]
            .value_counts(dropna=True).to_dict())
    return out

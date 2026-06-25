"""Assignment tiers -- the report-level verdict on every M0 claim (ROADMAP 2).

The pipeline commits ONE winning formula per peak, but presenting every commit
as "the identification" overstates what the evidence supports. This module
splits the committed assignments into two report tiers by MECHANICAL rules on
ledger columns (reproducible, no judgment calls at report time):

  Assigned   -- the formula is unique in the calibrated mass window, or it is
                corroborated by independent evidence (Mascope-confirmed
                isotopologues / attached satellites, the same neutral assigned
                in a second ionization channel, or series-anchor support), and
                nothing about the chemistry contradicts the validated sample
                profile.

  Candidate  -- a plausible formula, honestly ambiguous. Reasons: base
                confidence Low/Suspect; an effective-score near-tie; close
                alternatives in the window without isotope/cross-channel
                discrimination; the honest cross-family mass-degeneracy audit
                (degeneracy.py) finding many distinct plausible ions on the
                mass with no corroboration to break the tie (a per-pass
                "unique" claim that is only unique inside its narrow box);
                oxygen count beyond the validated chemistry (the O>=12 "lattice
                monsters" -- mass-fit fantasies sitting on the unexplained
                multi-halogen C/H lattice); or the mixed-BrCl family where the
                isotope pattern pins the halogens but not the backbone.

This module is degeneracy-AWARE: it reads the degeneracy_density /
degeneracy_note columns that degeneracy.apply_degeneracy() stamps, so degeneracy
MUST run before tiers (assign.run orders them that way; report.build_sheets only
re-tiers ledgers that already carry the stamp). On a ledger without those columns
(predates the audit) the degeneracy rule is simply inert.

The third tier of the report -- Below assignability -- is the UNEXPLAINED
residual, characterized peak-by-peak in residual.characterize_residual()
(isotope-partner / has-constraints / isolated).

Candidate DENSITY is the confidence currency: how many distinct formulas
live within CLOSE_MARGIN effective score of the winner (winner included).
Density 1 == unique. The stored alternatives list is capped upstream, so a
density equal to 1 + the cap is a lower bound, rendered as ">=N".

Works on OLD ledgers (CSV round-trips without the eff_score/eff_margin/tied
columns): the effective-score margin is recovered from the mechanical
commentary ("nearest competitor ... trails by X"), falling back to the raw
score gap. New ledgers carry the arbitration data directly.
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd

from . import chemistry as C
from . import ledger as L

__version__ = "0.5.0"  # offset-tolerant calibration (large systematic ppm offsets)

TIER_ASSIGNED = "Assigned"
TIER_CANDIDATE = "Candidate"

TIE_MARGIN = 0.05        # arbitrate()'s own near-tie window
CLOSE_MARGIN = 0.10      # alternative within this eff-score of winner = "close"
O_MAX_IDENTIFIED = 11    # validated chemistry tops out at O7 (monoterpene
                         # rungs); O8-11 is plausible HOM-type oxidation; O>=12
                         # only ever appeared as lattice-monster mass fits

# Mass-error-distribution test (Gao et al. 2024, Anal. Chem. 96:10210): the
# instrument's true mass error is a tight Gaussian (self-calibrated from the
# corroborated CHO/CHON core). A commit whose ppm error sits in the tail of
# that distribution AND carries no independent corroboration is an
# off-calibration mass fit -- the formula is at the edge of tolerance because
# the TRUE carrier is a different species. Demote it to Candidate. Corroborated
# commits (isotopologue / cross-channel / series-anchor) are never touched, so
# every flagship is structurally safe.
Z_TAIL_DEMOTE = 2.6      # |z| beyond which an UNCORROBORATED commit is demoted
CAL_MIN_N = 20           # need this many core peaks to trust the calibration
CAL_SIGMA_FLOOR = 0.15   # ppm; a lucky-tight core must not reject everything

# Background air-ion channels: opportunistic adducts the enumerator tries on top
# of the sample's real reagent ions (carbonate / superoxide / electron
# attachment). They are NOT a halide-CIMS sample's primary reagent, so they give
# the search a free +CO3/+O2 degree of freedom that O-rich mass fits exploit. A
# commit on such a channel is only Assigned-grade if (a) the channel is one of
# the sample's primary detected channels (e.g. a genuine CO3-CIMS run) or (b) the
# assignment is independently corroborated.
BACKGROUND_CHANNELS = ("[M+CO3]-", "[M+HBr+CO3]-", "[M+O2]-", "[M]-.")

# Honest cross-family mass degeneracy (degeneracy.measure_degeneracy, stamped as
# degeneracy_density / degeneracy_note). The per-pass candidate_density only
# counts competitors inside the ONE narrow element box that peak's pass
# enumerated; the degeneracy audit re-counts how many distinct plausible IONS
# fall in the calibrated window across ALL families (CHO/CHON, fluorinated, Si,
# S, halogen ...). A high count means the mass is not identifiable from accurate
# mass alone. So a commit that is degenerate at this honest level AND carries no
# extra-spectral corroboration (committed isotopologue child / stored
# isotopologue / second ionization channel / series anchor) must be capped at
# Candidate -- isotopes or a second channel are exactly the evidence that would
# break the degeneracy, and without them the unique-in-my-narrow-box claim is an
# artifact of the box, not of the spectrum. Corroborated commits are spared (the
# corroboration IS the tie-breaker), so isotopologue-confirmed / cross-channel
# backbones stay Assigned. Demote when the honest density exceeds this (>=3
# plausible ions) or the note carries the MASS-SATURATED flag (density > 8).
DEGEN_DEMOTE_DENSITY = 2

_TRAILS_RE = re.compile(r"trails by ([0-9.]+)")


def base_confidence(conf) -> str:
    """'Good (fluorinated)' -> 'Good'."""
    m = re.match(r"\s*([A-Za-z]+)", str(conf) if conf is not None else "")
    return m.group(1) if m else ""


def _alts(cell) -> list[dict]:
    try:
        v = json.loads(cell) if isinstance(cell, str) else (cell or [])
        return v if isinstance(v, list) else []
    except Exception:
        return []


def _truthy(v) -> bool | None:
    """Robust bool for a ledger 'tied' cell that may have round-tripped CSV."""
    if v is None or (isinstance(v, float) and np.isnan(v)) or v is pd.NA:
        return None
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    s = str(v).strip().lower()
    if s in ("true", "1"):
        return True
    if s in ("false", "0"):
        return False
    return None


def _degeneracy(row) -> tuple[int | None, bool]:
    """(honest degeneracy_density, is_mass_degenerate) from the degeneracy.py
    stamp. Degenerate when the cross-family count exceeds DEGEN_DEMOTE_DENSITY
    (>=3 plausible ions) OR the note carries the unambiguous MASS-SATURATED flag
    (density > 8). Returns (None, False) on ledgers that predate the degeneracy
    audit (no column) so the rule is simply inert there -- never an error."""
    d = row.get("degeneracy_density")
    density = None
    if d is not None and pd.notna(d):
        try:
            density = int(float(d))
        except (TypeError, ValueError):
            density = None
    note = row.get("degeneracy_note")
    note = "" if note is None or pd.isna(note) else str(note)
    flagged = "MASS-SATURATED" in note
    is_degen = (density is not None and density > DEGEN_DEMOTE_DENSITY) or flagged
    return density, is_degen


def _winner_raw(row) -> float | None:
    vals = [row.get("ion_score"), row.get("compound_score")]
    vals = [float(v) for v in vals if v is not None and pd.notna(v)]
    return min(vals) if vals else None


def _alt_raw(a: dict) -> float | None:
    for k in ("raw_score", "ion_score", "eff_score"):
        if a.get(k) is not None:
            return float(a[k])
    return None


# optional '^' after the sign = heavy-isotope-labelled reagent (e.g. +^NO3 for the
# ¹⁵N nitrate cluster); element COUNTS are isotope-independent so drop the marker.
_ADDUCT_TOKENS = re.compile(r"([+-])\^?([A-Za-z0-9]+)")


def _ion_counts(neutral, adduct) -> dict | None:
    """Element counts of the ION for a (neutral, adduct) reading, or None when
    the adduct string is not parseable. '[M+HBr+Br]-' adds H, 2x Br, etc."""
    if not neutral or not adduct:
        return None
    s = str(adduct).strip()
    if not s.startswith("[M"):
        return None
    cnt = dict(C.parse_formula(str(neutral)))
    for sign, tok in _ADDUCT_TOKENS.findall(s.split("]")[0][2:]):
        for el, n in C.parse_formula(tok).items():
            cnt[el] = cnt.get(el, 0) + (n if sign == "+" else -n)
    return {k: v for k, v in cnt.items() if v}


def _drop_decomposition_aliases(row, alts: list[dict]) -> tuple[list[dict], int]:
    """Remove alternatives that are the SAME ION as the winner under a
    different neutral/adduct split (covalent-vs-cluster decomposition).
    No spectral evidence can ever distinguish those readings -- the adduct
    reading is preferred by policy (2026-06-11) -- so they are NOT competing
    candidates and must not count toward ties or density."""
    ion0 = _ion_counts(row.get("neutral_formula"), row.get("adduct"))
    if ion0 is None:
        return alts, 0
    kept = [a for a in alts
            if _ion_counts(a.get("formula"), a.get("adduct")) != ion0]
    return kept, len(alts) - len(kept)


def _margin_density_tie(row, alts: list[dict], n_aliased: int,
                        n_stored: int) -> tuple[float | None, int, bool, bool]:
    """(margin to best real alternative, candidate density, capped?, tied?).

    margin None == no competing candidate in the search window. Aliases are
    already removed from `alts`; when any were removed, the stored tie flag /
    commentary tie may refer to an alias, so the verdict is recomputed."""
    if not alts:
        return None, 1, False, False
    eff = row.get("eff_score")
    if eff is not None and pd.notna(eff) \
            and any(a.get("eff_score") is not None for a in alts):
        # new ledgers: the arbitration's own effective scores, exact
        margins = [float(eff) - float(a["eff_score"]) for a in alts
                   if a.get("eff_score") is not None]
        margin = min(margins)
        n_close = sum(1 for m in margins if m < CLOSE_MARGIN)
        tied = margin < TIE_MARGIN
    else:
        # old ledger: the mechanical commentary holds the true eff margin --
        # but only trust it when no alias was filtered (it may name the alias)
        m = _TRAILS_RE.search(str(row.get("commentary") or ""))
        wr = _winner_raw(row)
        if m and n_aliased == 0:
            margin = float(m.group(1))
            tied = "(TIE)" in str(row.get("commentary") or "") or margin < TIE_MARGIN
            n_close = 1 if margin < CLOSE_MARGIN else 0
            if wr is not None and len(alts) > 1:
                n_close = max(n_close, sum(
                    1 for a in alts
                    if _alt_raw(a) is not None and wr - _alt_raw(a) < CLOSE_MARGIN))
        elif wr is not None:
            raws = [r for r in (_alt_raw(a) for a in alts) if r is not None]
            if not raws:
                return None, 1, False, False
            margin = wr - max(raws)
            n_close = sum(1 for r in raws if wr - r < CLOSE_MARGIN)
            tied = margin < TIE_MARGIN
        else:
            return None, 1 + len(alts), False, False
    # stored alternatives are capped upstream (arbitrate keeps the top few);
    # if every stored one is close, the true density is a lower bound
    capped = n_stored >= 3 and n_close >= len(alts) and n_close > 0
    return margin, 1 + n_close, capped, tied


def density_text(density: int, capped: bool) -> str:
    return f">={density}" if capped else str(density)


def _calibrate(m0: pd.DataFrame, kids_of: pd.Series) -> tuple[float, float] | None:
    """(mu, sigma) of the ppm-error distribution from the corroborated CHO/CHON
    core, or None if the core is too small. Robust (median + scaled MAD). The
    core is deliberately pure-organic (no halogen/Si/S, N<=1), High/Good, and
    isotopologue-backed -- the assignments we are most certain are correct, so
    their mass errors define the instrument's real accuracy for this run."""
    ppms = []
    for _, r in m0.iterrows():
        counts = C.parse_formula(str(r.get("neutral_formula") or ""))
        if any(counts.get(e, 0) for e in ("F", "Cl", "Br", "Si", "S")):
            continue
        if counts.get("N", 0) > 1:
            continue
        if base_confidence(r.get("confidence")) not in ("High", "Good"):
            continue
        if not ((kids_of.get(r["peak_id"], 0) > 0) or bool(_alts(r.get("isotopologues")))):
            continue
        p = r.get("ppm_error")
        if p is None or pd.isna(p) or abs(float(p)) > 15:   # gross-garbage guard
            continue
        ppms.append(float(p))
    if len(ppms) < CAL_MIN_N:
        return None
    s = pd.Series(ppms)
    # Outlier guard is RELATIVE to the robust median, not a 0-centered |ppm|<=2
    # window: an absolute window excludes the entire backbone of a source with a
    # large systematic offset (the uronium instrument sits at ~-2.4 ppm), leaving
    # the tier engine uncalibrated. Centering on the median keeps the same core
    # for a well-calibrated source (where the median is ~0) AND for an offset one.
    center = float(s.median())
    core = s[(s - center).abs() <= 2.0]
    if len(core) < CAL_MIN_N:
        core = s
    mu = float(core.median())
    sigma = max(float(1.4826 * (core - mu).abs().median()), CAL_SIGMA_FLOOR)
    return mu, sigma


def compute_tiers(ledger: pd.DataFrame) -> pd.DataFrame:
    """One row per M0 peak: [peak_id, tier, tier_reason, candidate_density,
    density_capped]. Pure; does not mutate the ledger."""
    m0 = ledger[ledger["role"] == L.ROLE_M0]
    # corroboration sources
    kids_of = ledger.loc[ledger["role"] == L.ROLE_ISO, "parent_peak_id"].value_counts()
    chan_count = m0.groupby("neutral_formula")["adduct"].nunique()
    cal = _calibrate(m0, kids_of)   # (mu, sigma) ppm, or None when uncalibrated
    # primary detected channels = adducts carrying the High pass-1 backbone
    # (the unambiguous real reagent ions). A background channel that shows up
    # here (a genuine CO3-CIMS run) is treated as primary, not demoted.
    _bb = m0[m0["method"].astype(str).str.startswith("cheminfo+grid")
             & m0["confidence"].astype(str).str.startswith("High")]
    primary_channels = set(_bb["adduct"].dropna())

    rows = []
    for _, r in m0.iterrows():
        formula = str(r.get("neutral_formula") or "")
        counts = C.parse_formula(formula)
        base = base_confidence(r.get("confidence"))
        alts_all = _alts(r.get("alternatives"))
        alts, n_aliased = _drop_decomposition_aliases(r, alts_all)
        margin, density, capped, tied = _margin_density_tie(
            r, alts, n_aliased, len(alts_all))
        if n_aliased == 0:
            # the arbitration's stored verdict is authoritative when no alias
            # polluted the runner-up slot
            stored = _truthy(r.get("tied"))
            if stored is not None:
                tied = stored
            elif "(TIE)" in str(r.get("commentary") or ""):
                tied = True
        iso_ev = (kids_of.get(r["peak_id"], 0) > 0) or bool(_alts(r.get("isotopologues")))
        cross_channel = int(chan_count.get(formula, 0)) >= 2
        has_anchor = pd.notna(r.get("anchor_peak_id")) or pd.notna(r.get("series_unit"))
        corroborated = iso_ev or cross_channel or has_anchor
        degen_density, mass_degenerate = _degeneracy(r)

        method = str(r.get("method") or "")
        tier, reason = TIER_ASSIGNED, ""
        if method.startswith("known:"):
            reason = ("known species (pass-0 locked list, mass + own-twin "
                      "self-consistency gated)")
        elif base in ("Low", "Suspect"):
            tier = TIER_CANDIDATE
            reason = (f"{base} confidence: score/mass evidence below the "
                      "identification bar")
        elif counts.get("O", 0) > O_MAX_IDENTIFIED:
            tier = TIER_CANDIDATE
            reason = (f"O{counts['O']} exceeds validated chemistry for this "
                      "matrix (flagships top out at O7); mass sits in the "
                      "unexplained C/H-lattice region -- likely a lattice "
                      "family member wearing a CHO(N) mass fit")
        elif counts.get("Br", 0) >= 1 and counts.get("Cl", 0) >= 1:
            tier = TIER_CANDIDATE
            reason = ("mixed Br/Cl halogenation: the isotope envelope pins the "
                      "halogen count but the backbone candidates stay ambiguous")
        elif tied and not (cross_channel or has_anchor):
            # a spectral eff-score tie cannot be broken by isotopes (they are
            # already in the score) -- only extra-spectral corroboration
            # (second channel, series anchor) rescues a tied winner
            alt0 = alts[0].get("formula", "?") if alts else "?"
            tier = TIER_CANDIDATE
            reason = (f"near-tie: best alternative ({alt0}) within "
                      f"{TIE_MARGIN} effective score")
        elif density > 1 and not corroborated:
            tier = TIER_CANDIDATE
            reason = (f"{density - 1} alternative(s) within {CLOSE_MARGIN} "
                      "effective score and no isotope / cross-channel / "
                      "series corroboration to discriminate")
        elif (str(r.get("adduct")) in BACKGROUND_CHANNELS
              and str(r.get("adduct")) not in primary_channels
              and not corroborated):
            tier = TIER_CANDIDATE
            reason = (f"background air-ion channel {r.get('adduct')} is not a "
                      "primary detected reagent channel for this sample and the "
                      "assignment has no isotope / cross-channel / series "
                      "corroboration (an unsupported +CO3/+O2 mass fit)")
        elif (cal is not None and not corroborated
              and pd.notna(r.get("ppm_error"))
              and abs((float(r["ppm_error"]) - cal[0]) / cal[1]) > Z_TAIL_DEMOTE):
            z = (float(r["ppm_error"]) - cal[0]) / cal[1]
            tier = TIER_CANDIDATE
            reason = (f"mass error {float(r['ppm_error']):+.2f} ppm is {z:+.1f}"
                      f"σ off the calibrated instrument accuracy "
                      f"({cal[0]:+.2f}±{cal[1]:.2f} ppm) with no isotope / "
                      "cross-channel / series corroboration "
                      "(mass-error-distribution test)")
        elif mass_degenerate and not corroborated:
            # the honest cross-family degeneracy audit says many distinct
            # plausible ions share this mass; "unique in the calibrated window"
            # was only true inside the narrow per-pass box. With no isotope /
            # cross-channel / series corroboration there is nothing to break the
            # tie, so the mass is NOT identifiable from accurate mass alone.
            tier = TIER_CANDIDATE
            if degen_density is not None:
                reason = (f"mass-degenerate: {degen_density} plausible formulas "
                          "share this mass within the calibrated window "
                          "(degeneracy audit) and no isotope / cross-channel / "
                          "series corroboration to break the tie — not "
                          "identifiable from accurate mass alone")
            else:
                reason = ("mass-saturated window (degeneracy audit) with no "
                          "isotope / cross-channel / series corroboration to "
                          "break the tie — not identifiable from accurate mass "
                          "alone")
        else:
            parts = []
            if density == 1:
                parts.append("unique formula in the calibrated window"
                             + (f" ({n_aliased} same-ion decomposition "
                                "reading(s) excluded)" if n_aliased else ""))
            else:
                parts.append(f"best of {density_text(density, capped)} candidates "
                             f"(margin {margin:.2f})" if margin is not None else
                             f"best of {density_text(density, capped)} candidates")
            if iso_ev:
                parts.append("isotopologue-confirmed")
            if cross_channel:
                parts.append("seen in a second ionization channel")
            if has_anchor:
                parts.append("series-anchor support")
            reason = "; ".join(parts)
        rows.append({"peak_id": r["peak_id"], "tier": tier, "tier_reason": reason,
                     "candidate_density": density, "density_capped": capped})
    return pd.DataFrame(rows, columns=["peak_id", "tier", "tier_reason",
                                       "candidate_density", "density_capped"])


def apply_tiers(ledger: pd.DataFrame) -> pd.DataFrame:
    """Stamp tier / tier_reason / candidate_density onto the M0 rows of the
    ledger (in place; returns the ledger). Non-M0 rows keep NA."""
    for col in ("tier", "tier_reason", "candidate_density"):
        if col not in ledger.columns:
            ledger[col] = pd.Series(pd.NA, index=ledger.index, dtype="object")
        elif ledger[col].dtype != object:
            # candidate_density holds '>=N' strings; a float column (e.g. an
            # all-NaN CSV round-trip) must widen before the stamp
            ledger[col] = ledger[col].astype("object")
    t = compute_tiers(ledger)
    if not len(t):
        return ledger
    idx = ledger.index[ledger["peak_id"].isin(t["peak_id"])]
    by_pid = t.set_index("peak_id")
    for i in idx:
        pid = ledger.at[i, "peak_id"]
        ledger.at[i, "tier"] = by_pid.at[pid, "tier"]
        ledger.at[i, "tier_reason"] = by_pid.at[pid, "tier_reason"]
        ledger.at[i, "candidate_density"] = density_text(
            int(by_pid.at[pid, "candidate_density"]),
            bool(by_pid.at[pid, "density_capped"]))
    flag_below_assignability(ledger)
    return ledger


def flag_below_assignability(ledger: pd.DataFrame) -> int:
    """Mark M0 commits that are high-oxygen 'monsters' AND mass-saturated as
    BELOW reliable assignability: their base mass fits, but ~dozens of distinct
    plausible ions sit within <=1 ppm (degeneracy.py MASS-SATURATED), so the
    formula is one arbitrary pick of a sub-ppm-degenerate set, not an ID. Stamp a
    `below_assignability` flag so the report lists them as a constrained mass, not
    a confident formula. They are already capped at Candidate by the tier rules;
    this is the explicit do-not-trust-the-formula disposition."""
    if "below_assignability" not in ledger.columns:
        ledger["below_assignability"] = False
    n = 0
    for i in ledger.index[ledger["role"] == L.ROLE_M0]:
        nf = str(ledger.at[i, "neutral_formula"])
        o = C.parse_formula(nf).get("O", 0) if nf and nf != "nan" else 0
        _density, is_degen = _degeneracy(ledger.loc[i])
        if o >= 11 and is_degen:
            ledger.at[i, "below_assignability"] = True
            ledger.at[i, "tier_reason"] = (str(ledger.at[i, "tier_reason"] or "")
                + " | below-assignability (O>=11, mass-saturated)").strip(" |")
            n += 1
    return n

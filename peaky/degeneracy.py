"""Honest mass-degeneracy measurement (single-file / sum-spectrum only).

The pipeline commits one winner per peak and its stored ``candidate_density`` is
relative to whichever NARROW element box that peak's pass enumerated (CHO/CHON,
or a single opened contaminant family). That undercounts cross-family
degeneracy: at heavy m/z a fluorinated, a CHNOS and a Si formula can all sit
within the instrument's mass accuracy and the per-pass density never compares
them. This module re-measures, for every committed M0 peak, how many distinct,
chemically-plausible IONS fall inside the *calibrated* mass window across all
sample channels -- the honest "how many things could this peak be" count.

Policy (user decision 2026-06-13, revised same day): this module is MEASUREMENT
ONLY -- it never sets a tier itself; it stamps ``degeneracy_density`` +
``degeneracy_note`` (the competing tie set) so a reader can see exactly how
identifiable each mass really is. The VERDICT lives in tiers.py, which now reads
these columns: an uncorroborated commit whose honest cross-family density is high
(or MASS-SATURATED) is capped at Candidate there. (The initial decision to stamp
without demoting was revised once the contradiction surfaced -- a row reading
``tier=Assigned, "unique formula"`` while also carrying ``degeneracy_density=27,
MASS-SATURATED`` is self-contradictory. Corroborated commits -- committed
isotopologue child / second channel / series anchor -- are still spared, because
that corroboration is exactly the extra-spectral evidence that breaks the tie;
the honest-but-uncorroborated F-heavy backbones land as Candidate, not a forced
Assigned.) Separation of concerns is preserved: degeneracy MEASURES, tiers
JUDGES, and degeneracy.apply_degeneracy must run before tiers.apply_tiers.

Same-ion decomposition readings (covalent ``Y(Br)[M-H]-`` vs cluster
``Y'.HBr.Br-``) are the SAME ion and are deduplicated, so they never inflate the
count -- exactly as tiers.py treats them.
"""
from __future__ import annotations

import dataclasses
from bisect import bisect_left, bisect_right

import pandas as pd

from . import chemistry as C
from . import contexts as X
from . import ledger as L
from . import tiers as T

__version__ = "0.1.0"

# Bounded enumeration box: wide enough to admit the real cross-family
# competitors the pipeline can open (F/Si/S/Cl/Br contaminants + CHON), tight
# enough that one full-range grid build is affordable. O is capped at
# fluorochemical/HOM levels (the v16 lesson: an open O range grids O-monster
# mass fits onto every heavy peak).
RELAXED_BOX = "C0-20 H0-36 N0-3 O0-12 S0-1 F0-17 Cl0-2 Br0-2 Si0-3"
# adducts a peak could plausibly be read under (the sample's channels + the
# opportunistic background ones the enumerator tries)
DEFAULT_ADDUCTS = ("[M+Br]-", "[M-H]-", "[M+CO3]-", "[M+HBr+Br]-", "[M+HBr+CO3]-")
K_SIGMA = 3.0        # calibrated half-window, in sigma
MAX_ALTS = 6         # competitors listed in the note
GRID_MASS_MAX = 730.0  # peaks top out ~720 Da; no need to enumerate to 900
# Realistic-competitor gate: real atmospheric/contaminant species carry at most
# a few heteroatom TYPES (CHO, CHON, CHOS, a single-halogen organic, a siloxane
# ...). Counting a formula that simultaneously mixes Br+Cl+F+S+Si as a rival
# wildly overstates the degeneracy (the relaxed box admits ~10^7 such monsters).
# So a competitor must use at most this many distinct heteroatom types.
_HET = ("N", "S", "P", "F", "Cl", "Br", "Si", "I")
MAX_HET_TYPES = 3
SATURATION_DENSITY = 8   # beyond this the mass is not identifiable by mass alone


def _het_types(counts: dict) -> int:
    return sum(1 for el in _HET if counts.get(el, 0) > 0)


def _disk_grid(box: str, mass_min: float = 30.0,
               mass_max: float = GRID_MASS_MAX) -> list[tuple[float, str]]:
    """Sorted (mass, formula) grid for the relaxed box. The full-box
    enumeration is ~minutes (the loop has no inner mass pruning), so it is
    cached IN-PROCESS by chemistry._grid_cached -- one build per run, reused by
    every peak. (An on-disk pickle was tried but the relaxed box is ~10M
    formulas / ~280 MB, which overran the temp filesystem; in-process caching is
    the right scope -- the degeneracy pass runs once per pipeline run.)"""
    return C._grid_cached(C.parse_ranges(box), mass_min, mass_max)


def relaxed_profile(base: X.ContextProfile) -> X.ContextProfile:
    """A copy of the context profile with the strict ambient F/Si/S/halogen caps
    raised to the levels the pipeline's contaminant families can open -- so the
    degeneracy count is over the formula space the pipeline COULD have committed,
    not the strict default box. Van Krevelen windows + DBE/Senior/O-cap are
    kept (they are real chemistry, not a heteroatom budget)."""
    return dataclasses.replace(base, max_F=18, max_S=2, max_Si=6,
                               max_Cl=2, max_Br=2, max_N=3)


def _canonical_ion(neutral: str, adduct: str) -> str | None:
    """A canonical ion-composition key, so covalent-vs-cluster aliases collapse
    to one (same physical ion = one candidate)."""
    cnt = T._ion_counts(neutral, adduct)
    if not cnt:
        return None
    return "".join(f"{el}{cnt[el]}" for el in sorted(cnt))


def measure_degeneracy(ledger: pd.DataFrame, *, cal: tuple[float, float] | None,
                       context: str = "ambient-air",
                       adducts=DEFAULT_ADDUCTS, box: str = RELAXED_BOX,
                       k_sigma: float = K_SIGMA, max_alts: int = MAX_ALTS,
                       log=lambda *_: None) -> dict[str, dict]:
    """Return {peak_id: {density, note, alts}} for every committed M0 peak.
    Pure; does not mutate the ledger. ``cal`` is (mu, sigma) ppm from
    tiers._calibrate (None -> uncalibrated, returns {})."""
    if cal is None:
        log("[degeneracy] uncalibrated; skipped")
        return {}
    mu, sigma = cal
    lo_ppm, hi_ppm = mu - k_sigma * sigma, mu + k_sigma * sigma
    grid = _disk_grid(box)
    masses = [g[0] for g in grid]
    profile = relaxed_profile(X.get_context(context))
    shifts = {a: C.ADDUCT_SHIFTS[a] for a in adducts if a in C.ADDUCT_SHIFTS}
    log(f"[degeneracy] grid {len(grid)} formulas; window [{lo_ppm:+.2f},{hi_ppm:+.2f}] ppm")

    m0 = ledger[ledger["role"] == L.ROLE_M0]
    out: dict[str, dict] = {}
    for _, r in m0.iterrows():
        mz = r.get("mz")
        if mz is None or pd.isna(mz):
            continue
        mz = float(mz)
        found: dict[str, tuple[str, str, float]] = {}   # ion -> (neutral, adduct, ppm)
        for adduct, shift in shifts.items():
            # ion mass window so the peak's ppm vs the candidate lands in
            # [lo_ppm, hi_ppm]; convert to a neutral-mass window for the grid.
            ion_lo = mz / (1 + hi_ppm * 1e-6)
            ion_hi = mz / (1 + lo_ppm * 1e-6)
            m_lo, m_hi = ion_lo - shift, ion_hi - shift
            if m_hi < 1:
                continue
            for i in range(bisect_left(masses, m_lo), bisect_right(masses, m_hi)):
                neutral = grid[i][1]
                if _het_types(C.parse_formula(neutral)) > MAX_HET_TYPES:
                    continue
                ion = _canonical_ion(neutral, adduct)
                if ion is None:
                    continue
                keep, _why = X.filter_by_profile(neutral, profile)
                if not keep:
                    continue
                theo_ion = grid[i][0] + shift
                ppm = (mz - theo_ion) / theo_ion * 1e6
                prev = found.get(ion)
                if prev is None or abs(ppm - mu) < abs(prev[2] - mu):
                    found[ion] = (neutral, adduct, ppm)

        # describe relative to the committed reading
        assigned_ion = _canonical_ion(r.get("neutral_formula"), r.get("adduct"))
        density = len(found)
        comp = sorted((v for k, v in found.items() if k != assigned_ion),
                      key=lambda t: abs(t[2] - mu))
        alts = [f"{n} {a} ({p:+.2f} ppm)" for n, a, p in comp[:max_alts]]
        if density <= 1:
            note = f"unique within ±{k_sigma:.0f}σ calibrated mass window"
        elif density > SATURATION_DENSITY:
            note = (f"MASS-SATURATED: {density} plausible formulas (≤{MAX_HET_TYPES} "
                    f"heteroatom types) within ±{k_sigma:.0f}σ calibrated window — "
                    "not identifiable from accurate mass alone; needs isotope "
                    "envelope / time-series corroboration")
        else:
            shown = "; ".join(alts)
            more = "" if len(comp) <= max_alts else f" (+{len(comp) - max_alts} more)"
            note = (f"MASS-DEGENERATE: {density} plausible ions within ±"
                    f"{k_sigma:.0f}σ calibrated window — competitors: {shown}{more}")
        out[str(r["peak_id"])] = {"density": density, "note": note, "alts": alts}
    return out


def apply_degeneracy(ledger: pd.DataFrame, *, cal=None, context="ambient-air",
                     log=lambda *_: None, **kw) -> pd.DataFrame:
    """Stamp degeneracy_density / degeneracy_note onto the M0 rows (in place)."""
    for col in ("degeneracy_density", "degeneracy_note"):
        if col not in ledger.columns:
            ledger[col] = pd.Series(pd.NA, index=ledger.index, dtype="object")
        elif ledger[col].dtype != object:
            ledger[col] = ledger[col].astype("object")
    if cal is None:
        cal = T._calibrate(
            ledger[ledger["role"] == L.ROLE_M0],
            ledger.loc[ledger["role"] == L.ROLE_ISO, "parent_peak_id"].value_counts())
    res = measure_degeneracy(ledger, cal=cal, context=context, log=log, **kw)
    if not res:
        return ledger
    for i in ledger.index[ledger["role"] == L.ROLE_M0]:
        pid = str(ledger.at[i, "peak_id"])
        if pid in res:
            ledger.at[i, "degeneracy_density"] = res[pid]["density"]
            ledger.at[i, "degeneracy_note"] = res[pid]["note"]
    return ledger

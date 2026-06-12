"""Isotope-pattern prescan -> grid constraints.

This is deliberately NOT a scorer. Mascope (match_compounds) is the authoritative
judge of whether an isotopologue pattern fits a formula. The prescan only looks
at the raw peak list to answer two cheap questions that shrink the candidate
search before we ever call the server:

  1. Which heteroatoms show isotope-pair evidence (Br, Cl, S, Si)? -> only put
     those elements in the grid ranges (huge combinatorial saving).
  2. What is the largest carbon number implied by the brightest 13C satellites?
     -> cap C in the grid.

It walks the deduplicated, intensity-sorted peak list looking for satellite
pairs at characteristic delta-m with plausible intensity ratios.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

__version__ = "0.1.0"

# delta-m (Da) between an isotopologue satellite and its monoisotopic parent
D_13C = 1.003355
D_37CL = 1.997050
D_81BR = 1.997795
D_34S = 1.995796
D_29SI = 0.999568
D_18O = 2.004246

# natural-abundance per-atom ratio of the +N satellite to the monoisotopic peak
R_13C_PER_C = 0.0107          # 1.1% per carbon
R_37CL_PER_CL = 0.3196        # 37Cl/35Cl
R_81BR_PER_BR = 0.9728        # 81Br/79Br
R_34S_PER_S = 0.0443          # 34S/32S
R_29SI_PER_SI = 0.0510        # 29Si/28Si


@dataclass
class PrescanResult:
    has_Br: bool = False
    has_Cl: bool = False
    has_S: bool = False
    has_Si: bool = False
    has_multi_Br: bool = False          # Br2 triplet seen
    estimated_max_C: int = 0
    evidence: list = field(default_factory=list)   # human-readable hits

    def as_dict(self) -> dict:
        return {
            "has_Br": self.has_Br, "has_Cl": self.has_Cl, "has_S": self.has_S,
            "has_Si": self.has_Si, "has_multi_Br": self.has_multi_Br,
            "estimated_max_C": self.estimated_max_C, "n_evidence": len(self.evidence),
        }


def _find_partner(mz_sorted, height_by_mz, target_mz, ppm_tol):
    """Return (mz, height) of a peak near target_mz within ppm_tol, else None."""
    tol = target_mz * ppm_tol * 1e-6
    lo, hi = target_mz - tol, target_mz + tol
    import bisect
    i = bisect.bisect_left(mz_sorted, lo)
    best = None
    while i < len(mz_sorted) and mz_sorted[i] <= hi:
        m = mz_sorted[i]
        h = height_by_mz[m]
        if best is None or h > best[1]:
            best = (m, h)
        i += 1
    return best


def prescan(peaks: pd.DataFrame, *, mz_col="mz", height_col="height",
            ppm_tol=8.0, min_height=0.0,
            reagent_mzs: list[float] | None = None,
            reagent_ppm=15.0) -> PrescanResult:
    """Scan a peak table for isotope-pair signatures.

    `reagent_mzs` lets the caller strip known reagent-cluster peaks (e.g. bare
    Br_n clusters) before the Br scan, so they don't masquerade as analyte Br.
    """
    res = PrescanResult()
    df = peaks[[mz_col, height_col]].dropna().copy()
    df = df[df[height_col] >= min_height]
    if reagent_mzs:
        keep = []
        for mz in df[mz_col]:
            is_reagent = any(abs(mz - r) / r * 1e6 <= reagent_ppm for r in reagent_mzs)
            keep.append(not is_reagent)
        df = df[keep]
    df = df.sort_values(height_col, ascending=False)
    mz_sorted = sorted(df[mz_col].tolist())
    height_by_mz = dict(zip(df[mz_col], df[height_col]))

    n_scan = min(len(df), 400)   # brightest peaks carry the isotope information
    for parent_mz, parent_h in zip(df[mz_col].head(n_scan), df[height_col].head(n_scan)):
        if parent_h <= 0:
            continue
        # --- 13C: estimate carbon count from the +1.00336 satellite ratio ---
        p = _find_partner(mz_sorted, height_by_mz, parent_mz + D_13C, ppm_tol)
        if p:
            ratio = p[1] / parent_h
            if 0.003 <= ratio <= 0.9:
                n_c = round(ratio / R_13C_PER_C)
                if n_c > res.estimated_max_C:
                    res.estimated_max_C = int(n_c)
        # --- 81Br: +1.99795. One Br -> M+2/M ~ 1.0; Br2 (1:2:1) -> M+2/M ~ 2.0
        #     with an M+4 at ~1.0*M. Either pattern confirms Br. ---
        p = _find_partner(mz_sorted, height_by_mz, parent_mz + D_81BR, ppm_tol)
        if p:
            ratio = p[1] / parent_h
            if 0.6 <= ratio <= 1.4:        # single Br
                res.has_Br = True
                res.evidence.append(("Br", round(parent_mz, 4), round(ratio, 2)))
                continue
            if 1.5 <= ratio <= 2.6:        # Br2: check M+4 ~ 1.0*M
                p2 = _find_partner(mz_sorted, height_by_mz, parent_mz + 2 * D_81BR, ppm_tol)
                if p2 and 0.6 <= p2[1] / parent_h <= 1.4:
                    res.has_Br = True
                    res.has_multi_Br = True
                    res.evidence.append(("Br2", round(parent_mz, 4), round(ratio, 2)))
                    continue
        # --- 37Cl: +1.99705, ratio ~0.32 (one Cl) ---
        p = _find_partner(mz_sorted, height_by_mz, parent_mz + D_37CL, ppm_tol)
        if p:
            ratio = p[1] / parent_h
            if 0.22 <= ratio <= 0.45:
                res.has_Cl = True
                res.evidence.append(("Cl", round(parent_mz, 4), round(ratio, 2)))
        # --- 34S: +1.9958, ratio ~0.045 per S ---
        p = _find_partner(mz_sorted, height_by_mz, parent_mz + D_34S, ppm_tol)
        if p:
            ratio = p[1] / parent_h
            if 0.025 <= ratio <= 0.09:
                res.has_S = True
                res.evidence.append(("S", round(parent_mz, 4), round(ratio, 2)))
        # --- 29Si: +0.99957, ratio ~0.05 per Si ---
        p = _find_partner(mz_sorted, height_by_mz, parent_mz + D_29SI, ppm_tol)
        if p:
            ratio = p[1] / parent_h
            if 0.035 <= ratio <= 0.08:
                res.has_Si = True
                res.evidence.append(("Si", round(parent_mz, 4), round(ratio, 2)))
    return res


def constrain_ranges(base_ranges: dict[str, tuple[int, int]],
                     pre: PrescanResult,
                     context_caps: dict[str, int]) -> dict[str, tuple[int, int]]:
    """Apply prescan evidence to grid ranges:
      * cap C at estimated_max_C (+ small headroom) when we have an estimate
      * zero out Br/Cl/S/Si that show NO spectral evidence (subject to context)
    `context_caps` gives the per-element max the context allows."""
    r = dict(base_ranges)
    if pre.estimated_max_C and "C" in r:
        cmax = min(r["C"][1], pre.estimated_max_C + 4)
        r["C"] = (r["C"][0], max(cmax, r["C"][0]))
    for el, flag in (("Br", pre.has_Br), ("Cl", pre.has_Cl),
                     ("S", pre.has_S), ("Si", pre.has_Si)):
        cap = context_caps.get(el, 0)
        if not flag or cap <= 0:
            r[el] = (0, 0)
        else:
            r[el] = (0, min(r.get(el, (0, cap))[1] or cap, cap))
    return r

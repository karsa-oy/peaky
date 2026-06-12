"""Pass 4 -- the residual explainer.

Runs after Passes 1-3 + reagent labeling and attacks what is left, using
PATTERN evidence instead of chemical-prior filters:

  Stage A  Isotope-pair resolution. The residual of a halide-CIMS spectrum is
           dominated by ~1.998-Da doublets (Br/Cl isotope envelopes). Each
           doublet's light member becomes an M0 hypothesis whose ion MUST
           contain the implied halogen count; candidates come from the local
           grid; Mascope match_compounds validates the full isotope pattern.
           Explaining the light member automatically claims the heavy partner
           as the 81Br/37Cl child.

  Stage B  Deep series propagation. Bright residual peaks that sit 1-2 exact
           repeat/cluster units (CH2, O, H2O, CO, CO2, C2H2O, C2H4O2, HX) from
           an assigned anchor are proposed as series/cluster extensions and
           validated by Mascope.

Acceptance policy (user directive 2026-06-10):
  * |ppm| <= ppm_strict (1.0): accept on score alone.
  * ppm_strict < |ppm| <= ppm_pattern (4.0): accept ONLY with pattern evidence
    (Mascope-confirmed isotopologue partner, or >=2 supporting series anchors),
    capped at Low.
  * Plausibility filter is DBE-ONLY (integer DBE >= 0 + Senior's rule, enforced
    structurally by the grid). No H/C / O/C ratio gates in this pass --
    pattern evidence replaces the chemical priors.
  * Confidence is capped at Good (never High) -- pattern-driven assignments
    stay below the locked Pass-1 backbone.

All stages take an injectable `score_fn` (defaults to io_mascope.score_candidates)
so the acceptance logic is unit-testable offline.
"""
from __future__ import annotations

import bisect

import numpy as np
import pandas as pd

from . import chemistry as C
from . import io_mascope as IO
from . import ledger as L
from . import series_gka as G
from .passes import PassConfig, arbitrate, confidence_label, z_of, _f

__version__ = "0.2.0"

# isotope spacings
D_PAIR_BR = 1.997795
D_PAIR_CL = 1.997050


# ---------------------------------------------------------------------------
# Stage A helpers
# ---------------------------------------------------------------------------
def find_iso_pairs(ledger: pd.DataFrame, *, ppm_tol: float = 8.0,
                   min_height: float = 0.0) -> pd.DataFrame:
    """Find ~1.998-Da doublets within the UNEXPLAINED residual.

    Returns DataFrame[light_pid, heavy_pid, light_mz, ratio, element, n_halogen]
    where element is 'Br' (ratio ~1 per Br) or 'Cl' (ratio ~0.32 per Cl) and
    n_halogen is the implied halogen count in the ION (1 or 2)."""
    un = ledger[(ledger["role"] == L.ROLE_UNEXPLAINED)
                & (ledger["height"].fillna(0) >= min_height)]
    rows = []
    mzs = un["mz"].to_numpy()
    order = np.argsort(mzs)
    ms = mzs[order]
    hs = un["height"].to_numpy()[order]
    pids = un["peak_id"].to_numpy()[order]
    claimed_heavy: set = set()
    for i in range(len(ms)):
        if hs[i] <= 0:
            continue
        for delta, element in ((D_PAIR_BR, "Br"), (D_PAIR_CL, "Cl")):
            target = ms[i] + delta
            tol = target * ppm_tol * 1e-6
            j = bisect.bisect_left(ms, target - tol)
            best = None
            while j < len(ms) and ms[j] <= target + tol:
                if pids[j] not in claimed_heavy and j != i:
                    if best is None or hs[j] > best[1]:
                        best = (j, hs[j])
                j += 1
            if best is None:
                continue
            ratio = best[1] / hs[i]
            n_hal = 0
            if element == "Br":
                if 0.6 <= ratio <= 1.4:
                    n_hal = 1
                elif 1.5 <= ratio <= 2.6:
                    n_hal = 2
            else:
                if 0.22 <= ratio <= 0.45:
                    n_hal = 1
                elif 0.5 <= ratio <= 0.85:
                    n_hal = 2
            if n_hal:
                rows.append({"light_pid": pids[i], "heavy_pid": pids[best[0]],
                             "light_mz": float(ms[i]), "ratio": float(ratio),
                             "element": element, "n_halogen": n_hal})
                claimed_heavy.add(pids[best[0]])
                break   # one element interpretation per light peak
    return pd.DataFrame(rows)


def _halogen_count_in_ion(neutral: str, adduct: str, element: str) -> int:
    n = C.parse_formula(neutral).get(element, 0)
    if element == "Br" and "Br" in adduct:
        n += 1
    if element == "Cl" and "Cl" in adduct:
        n += 1
    return n


def candidates_for_pair(light_mz: float, element: str, n_halogen: int,
                        adducts: list[str], *, ranges: dict,
                        ppm: float) -> set[str]:
    """Grid candidates for a pair's light member: the ion must contain exactly
    n_halogen atoms of `element`. DBE-only filtering (grid-structural)."""
    out: set[str] = set()
    for adduct in adducts:
        if adduct not in C.ADDUCT_SHIFTS:
            continue
        # halogen needed in the NEUTRAL given this adduct
        need = n_halogen - (1 if element in adduct else 0)
        if need < 0:
            continue
        r = dict(ranges)
        r[element] = (need, need)
        for f in C.candidates_for_peaks([light_mz], r, [adduct], ppm_tolerance=ppm):
            if C.parse_formula(f).get(element, 0) == need:
                out.add(f)
    return out


def _accept(score: float, ppm: float | None, has_pattern: bool,
            cfg: PassConfig) -> tuple[bool, str]:
    """The Pass-4 acceptance rule. Returns (accept, reason).

    When the run is self-calibrated (cfg.cal_mu/cal_sigma fitted on the pass-1
    backbone) the strict/pattern bands are calibrated z-scores; otherwise the
    static residual_ppm_* windows apply."""
    if score is None or score < cfg.tau_suspect:
        return False, "score below floor"
    if ppm is None or pd.isna(ppm):
        return False, "no mass error reported"
    z = z_of(ppm, cfg)
    if z is not None:
        if z <= cfg.cal_z_accept:
            return True, f"z={z:.1f} within calibrated accuracy"
        if z <= cfg.cal_z_pattern and has_pattern:
            return True, f"z={z:.1f} > accept band but pattern evidence corroborates"
        return False, f"z={z:.1f} without pattern evidence"
    a = abs(ppm)
    if a <= cfg.residual_ppm_strict:
        return True, f"|ppm|={a:.2f} within strict tolerance"
    if a <= cfg.residual_ppm_pattern and has_pattern:
        return True, f"|ppm|={a:.2f} > strict but pattern evidence corroborates"
    return False, f"|ppm|={a:.2f} without pattern evidence"


def _cap_conf(conf: str) -> str:
    return conf.replace("High", "Good")


# ---------------------------------------------------------------------------
# Stage A -- isotope-pair resolution
# ---------------------------------------------------------------------------
def stage_a_iso_pairs(client, sample_id: str, ledger: pd.DataFrame, profile,
                      pre, cfg: PassConfig, adducts: list[str], *,
                      score_fn=None, log=print) -> dict:
    score_fn = score_fn or IO.score_candidates
    out = {"committed": 0, "locked": 0, "iso_attached": 0}
    pairs = find_iso_pairs(ledger, min_height=cfg.height_cutoff)
    if len(pairs) == 0:
        log("[pass4.A] no isotope pairs in residual")
        return out
    log(f"[pass4.A] {len(pairs)} isotope pairs in residual "
        f"({(pairs.element == 'Br').sum()} Br-like, {(pairs.element == 'Cl').sum()} Cl-like)")
    cmax = max(12, (pre.estimated_max_C + 4) if getattr(pre, "estimated_max_C", 0) else 40)
    cmax = min(cmax, 40)
    base_ranges = {"C": (1, cmax), "H": (0, 2 * cmax + 4), "O": (0, 30),
                   "N": (0, profile.max_N), "S": (0, min(1, profile.max_S))}
    # enumerate per pair, pooled scoring
    cand_by_formula: dict[str, list] = {}
    for _, p in pairs.iterrows():
        cands = candidates_for_pair(p["light_mz"], p["element"], p["n_halogen"],
                                    adducts, ranges=base_ranges,
                                    ppm=cfg.residual_ppm_pattern)
        for f in cands:
            cand_by_formula.setdefault(f, []).append(p["light_pid"])
    if not cand_by_formula:
        log("[pass4.A] no grid candidates for any pair")
        return out
    log(f"[pass4.A] {len(cand_by_formula)} candidate formulas (DBE-only filter)")
    scored = score_fn(client, sample_id, sorted(cand_by_formula), mechanism_ids=cfg.mechanism_ids)
    if scored is None or len(scored) == 0:
        return out
    arb = arbitrate(scored, cfg)
    win = arb.get("winners", pd.DataFrame())
    kids = arb.get("iso_children", pd.DataFrame())
    pair_by_light = {p["light_pid"]: p for _, p in pairs.iterrows()}
    for _, w in win.iterrows():
        pid = w["peak_id"]
        if pid not in pair_by_light:
            continue   # only commit onto pair light members in this stage
        p = pair_by_light[pid]
        has_pattern = w["n_iso"] >= 1   # Mascope-confirmed satellite(s)
        ok, why = _accept(w["raw_score"], w["ppm_error"], has_pattern, cfg)
        if not ok:
            continue
        try:
            if L.role_of(ledger, pid) != L.ROLE_UNEXPLAINED:
                continue
            conf = _cap_conf(confidence_label(
                w["raw_score"], w["ppm_error"], w["n_iso"], w["tied"], cfg,
                suffix="iso-pair"))
            if conf == "Reject":
                continue
            L.commit_assignment(
                ledger, pid, neutral_formula=w["neutral"], adduct=w["adduct"],
                ion_formula=w["ion_formula"], ion_score=w["ion_score"],
                compound_score=w["compound_score"], ppm_error=w["ppm_error"],
                pass_no=4, method="residual:iso-pair", confidence=conf,
                commentary=(f"Pass 4 (iso-pair): {p['element']} doublet "
                            f"(ratio {p['ratio']:.2f}, n_{p['element']}="
                            f"{p['n_halogen']} in ion) anchors this peak; "
                            f"{w['neutral']} {w['adduct']} scored "
                            f"{w['raw_score']:.2f} by Mascope. Accepted: {why}. "
                            f"DBE-only plausibility (no ratio filters)."),
                alternatives=w["alternatives"])
            out["committed"] += 1
            # attach the heavy partner: prefer Mascope's own attribution,
            # fall back to the observed-pair labeling
            attached = False
            if len(kids):
                mine = kids[(kids["parent_peak_id"] == pid)]
                for _, k in mine.iterrows():
                    try:
                        L.attach_isotopologue(ledger, k["peak_id"], pid,
                                              iso_label=k["iso_label"],
                                              iso_match_score=k["iso_score"])
                        out["iso_attached"] += 1
                        if k["peak_id"] == p["heavy_pid"]:
                            attached = True
                    except L.LedgerError:
                        continue
            if not attached:
                try:
                    L.attach_isotopologue(
                        ledger, p["heavy_pid"], pid,
                        iso_label=("81Br" if p["element"] == "Br" else "37Cl")
                        + "(pair)")
                    out["iso_attached"] += 1
                except L.LedgerError:
                    pass
        except L.LedgerError:
            continue
    log(f"[pass4.A] {out}")
    return out


# ---------------------------------------------------------------------------
# Stage B -- deep series / cluster propagation
# ---------------------------------------------------------------------------
def stage_b_series(client, sample_id: str, ledger: pd.DataFrame, profile,
                   cfg: PassConfig, adducts: list[str], *, reagent: str | None,
                   score_fn=None, log=print) -> dict:
    score_fn = score_fn or IO.score_candidates
    out = {"committed": 0, "locked": 0, "iso_attached": 0}
    anchors = set(ledger.loc[ledger["role"] == L.ROLE_M0, "neutral_formula"].dropna())
    if not anchors:
        return out
    units = tuple(G.ORGANIC_UNITS) + ("C2H4O2",)
    if reagent in ("Br", "Cl"):
        units = units + ("H" + reagent,)
    un = ledger[(ledger["role"] == L.ROLE_UNEXPLAINED)
                & (ledger["height"].fillna(0) >= cfg.height_cutoff)]
    proposals: dict[str, dict] = {}
    for _, prow in un.iterrows():
        for prop in G.propose_for_peak(prow["mz"], anchors, adducts, units=units,
                                       ppm=cfg.residual_ppm_pattern,
                                       max_steps=cfg.residual_max_steps):
            f = prop.neutral_formula
            ok, _why = C.dbe_ok(f)        # structural gates: DBE + oxygen cap
            if ok:
                ok, _why = C.oxygen_ok(f)
            if not ok:
                continue
            cur = proposals.get(f)
            if cur is None or prop.n_supporting_anchors > cur["support"]:
                proposals[f] = {"support": prop.n_supporting_anchors,
                                "anchor": prop.anchor_formula,
                                "unit": prop.unit, "steps": prop.n_steps}
    proposals = {f: v for f, v in proposals.items() if f not in anchors}
    if not proposals:
        log("[pass4.B] no series proposals")
        return out
    log(f"[pass4.B] {len(proposals)} deep-series proposals (<= {cfg.residual_max_steps} steps)")
    scored = score_fn(client, sample_id, sorted(proposals), mechanism_ids=cfg.mechanism_ids)
    if scored is None or len(scored) == 0:
        return out
    arb = arbitrate(scored, cfg)
    for _, w in arb.get("winners", pd.DataFrame()).iterrows():
        f = w["neutral"]
        meta = proposals.get(f)
        if meta is None:
            continue
        has_pattern = (w["n_iso"] >= 1) or (meta["support"] >= 2)
        ok, why = _accept(w["raw_score"], w["ppm_error"], has_pattern, cfg)
        if not ok:
            continue
        pid = w["peak_id"]
        try:
            if L.role_of(ledger, pid) != L.ROLE_UNEXPLAINED:
                continue
            conf = _cap_conf(confidence_label(
                w["raw_score"], w["ppm_error"], w["n_iso"], w["tied"], cfg,
                suffix="deep-series"))
            if conf == "Reject":
                continue
            L.commit_assignment(
                ledger, pid, neutral_formula=f, adduct=w["adduct"],
                ion_formula=w["ion_formula"], ion_score=w["ion_score"],
                compound_score=w["compound_score"], ppm_error=w["ppm_error"],
                pass_no=4, method="residual:series", confidence=conf,
                commentary=(f"Pass 4 (deep series): {meta['steps']:+d}x"
                            f"{meta['unit']} from anchor {meta['anchor']} "
                            f"({meta['support']} supporting anchors); scored "
                            f"{w['raw_score']:.2f}. Accepted: {why}."),
                alternatives=w["alternatives"], series_unit=meta["unit"])
            out["committed"] += 1
        except L.LedgerError:
            continue
    log(f"[pass4.B] {out}")
    return out


def explain_residual(client, sample_id: str, ledger: pd.DataFrame, profile,
                     pre, cfg: PassConfig, adducts: list[str], *,
                     reagent: str | None = None, score_fn=None,
                     log=print) -> dict:
    a = stage_a_iso_pairs(client, sample_id, ledger, profile, pre, cfg, adducts,
                          score_fn=score_fn, log=log)
    b = stage_b_series(client, sample_id, ledger, profile, cfg, adducts,
                       reagent=reagent, score_fn=score_fn, log=log)
    return {k: a[k] + b[k] for k in a}

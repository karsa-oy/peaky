"""Post-assignment residual cleanup (single-file, runs after pass 6).

Three honest reclassifications of the leftover 'unexplained' residual, none of
which invent chemistry the spectrum doesn't support:

  1. flag_ringing_artifacts -- a weak peak sitting within a few mDa of a MUCH
     brighter peak is a peak-shape / ringing / shoulder artifact of that bright
     ion, not a real ion. Role -> 'artifact'.

  2. label_bromide_clusters -- a strongly negative mass-defect peak (defect <
     -0.16 => >=2 reagent Br) that carries a Br isotope partner and has no sane
     covalent reading is a bromide reagent-cluster ion. Role -> 'reagent'.

  3. recover_isotope_gated -- the genuine molecules the score gate dropped:
     enumerate ONLY low-complexity CHO/CHON/CHOS (+<=1 covalent Br/Cl), and
     commit a winner ONLY when the measured halogen isotope pattern CONFIRMS the
     candidate's halogen count (an independent corroboration that justifies the
     relaxed score floor). The halogen-isotope requirement + complexity clamp is
     what keeps the N2S2 / O-monster mass-fit flood out.
"""
from __future__ import annotations

import bisect

import pandas as pd

from . import chemistry as C
from . import contexts as X
from . import degeneracy as D
from . import io_mascope as IO
from . import ledger as L

__version__ = "0.4.0"   # + reclaim_satellites (attach leaked 13C/81Br/37Cl satellites)
# (history) v43-review fixes: ringing brightness floor (H4),
                        # CHO-only isotope-confirmed recovery (H1/H2), covalent
                        # cluster commentary (M4)

BR = 1.99795          # 81Br - 79Br
C13 = 1.003355

# ---- 1. ringing / shoulder artifacts -------------------------------------
# FT ringing / detector-shoulder satellites only appear next to a SATURATING
# ion, and their amplitude scales UP with the parent intensity. So a real
# artifact needs (a) the parent to be genuinely bright (>= MIN_PARENT cps -- a
# dim 6k-cps peak cannot ring), (b) a large brightness factor (>= the floor, so
# the satellite is <~1% of the parent), and (c) sub-resolution proximity. These
# three gates are what separate a true sidelobe from a resolved independent ion
# that merely happens to sit nearby (the 273.05 false flag, parent only 6385 cps
# at 43x and 2.3% -- physics runs backwards, v43 review H4).
MIN_RING_PARENT = 50000.0   # cps; below this a peak is too weak to ring
RING_FACTOR = 100.0         # parent must be >= this x brighter (=> satellite <1%)


def flag_ringing_artifacts(ledger: pd.DataFrame, *, factor: float = RING_FACTOR,
                           dmz: float = 0.012, min_parent: float = MIN_RING_PARENT,
                           log=print) -> dict:
    """Mark an unexplained peak as 'artifact' when a SATURATING peak (>= min_parent
    cps) that is >=`factor`x brighter sits within `dmz`. The tight dmz + the
    brightness floor + the high factor together ensure only sub-resolution
    sidelobes of a saturating ion are flagged, never a resolved independent
    neighbour."""
    pk = ledger.dropna(subset=["mz"]).sort_values("mz")
    mzs = pk["mz"].tolist(); hts = pk["height"].tolist()
    n = 0
    for pid, mz, h in zip(ledger["peak_id"], ledger["mz"], ledger["height"]):
        if ledger.loc[ledger.peak_id == pid, "role"].iloc[0] != L.ROLE_UNEXPLAINED:
            continue
        if pd.isna(mz) or pd.isna(h) or h <= 0:
            continue
        lo = bisect.bisect_left(mzs, mz - dmz)
        hit = None
        for i in range(lo, len(mzs)):
            if mzs[i] > mz + dmz:
                break
            if (abs(mzs[i] - mz) > 1e-4 and hts[i] >= min_parent
                    and hts[i] >= factor * h):
                if hit is None or hts[i] > hit[1]:
                    hit = (mzs[i], hts[i])
        if hit is not None:
            L.mark_artifact(ledger, pid,
                            f"FT ringing/sidelobe of {hit[0]:.4f} "
                            f"({hit[1]:.0f} cps, {hit[1] / h:.0f}x brighter, "
                            f"Δ{(hit[0] - mz) * 1000:+.0f} mDa, satellite "
                            f"{100 * h / hit[1]:.1f}% of parent)")
            n += 1
    log(f"[cleanup] flagged {n} ringing/sidelobe artifacts")
    return {"flagged": n}


# ---- 2. bromide reagent clusters -----------------------------------------
# A covalent di-/tri-bromo organic is formula-DEGENERATE with the bromide-
# cluster ion (same ion formula), so the reagent-adduct reading stays the
# correct CALL -- but the commentary must not assert "no covalent reading" when
# the oracle gives that covalent neutral a high score. So we score the
# mass-feasible covalent bromo-organics before wording the note.
CLUSTER_COVALENT_BOX = "C0-12 H0-22 N0-1 O0-8 S0-1 Cl0-1 Br1-3"
CLUSTER_COVALENT_TAU = 0.70


def label_bromide_clusters(ledger: pd.DataFrame, client=None, sample_id=None, *,
                           score_fn=None, defect_max: float = -0.16,
                           covalent_tau: float = CLUSTER_COVALENT_TAU,
                           log=print) -> dict:
    """Label strongly negative-defect unexplained peaks that carry a Br isotope
    partner (1.998 Da, ratio 0.4-3.0) as reagent bromide-cluster ions.

    When the oracle is available, each peak's commentary first tests whether a
    covalent di-/tri-bromo organic fits (the formula-degenerate reading): a
    high-scoring tie gets the "reagent-adduct reading preferred over degenerate
    di-bromo organic" wording (the pipeline's existing precedent), NOT the false
    "no covalent reading". The role stays reagent either way."""
    pk = ledger.dropna(subset=["mz"]).sort_values("mz")
    mzs = pk["mz"].tolist(); hts = pk["height"].tolist()
    def near(t, tol=0.006):
        i = bisect.bisect_left(mzs, t - tol)
        best = None
        while i < len(mzs) and mzs[i] <= t + tol:
            if best is None or abs(mzs[i] - t) < abs(best[0] - t):
                best = (mzs[i], hts[i])
            i += 1
        return best

    # phase 1: collect the cluster peaks (defect + Br-twin gate), before marking
    targets = []   # (peak_id, mz)
    for pid, mz, h in zip(list(ledger["peak_id"]), list(ledger["mz"]), list(ledger["height"])):
        if ledger.loc[ledger.peak_id == pid, "role"].iloc[0] != L.ROLE_UNEXPLAINED:
            continue
        if pd.isna(mz) or pd.isna(h):
            continue
        defect = float(mz) - round(float(mz))
        if defect >= defect_max:
            continue
        up, dn = near(mz + BR), near(mz - BR)
        twin = (up is not None and 0.4 < up[1] / h < 3.0) or \
               (dn is not None and h > 0 and 0.4 < (h / dn[1] if dn[1] else 9) < 3.0)
        if twin:
            targets.append((pid, float(mz)))

    # phase 2: covalent-fit oracle check (only if a scorer is available)
    best_cov = {}   # peak_id -> (ion_formula, score)
    if targets and (score_fn is not None or client is not None):
        sf = score_fn or IO.score_candidates
        box = C.parse_ranges(CLUSTER_COVALENT_BOX)
        per, allf = {}, set()
        for pid, mz in targets:
            cands = {f for f in C.candidates_for_peaks(
                        [mz], box, ["[M-H]-", "[M+Br]-"], ppm_tolerance=2.0,
                        mass_min=mz - 90, mass_max=mz + 2)
                     if C.parse_formula(f).get("Br", 0) >= 1}
            per[pid] = cands; allf |= cands
        if allf:
            try:
                fr = sf(client, sample_id, sorted(allf), allow_partial=True)
                fr = fr[fr["sample_peak_id"].notna() & (fr["sample_peak_intensity"] > 0)]
                for pid, mz in targets:
                    sub = fr[(fr["sample_peak_mz"] - mz).abs() < 0.006]
                    if not sub.empty:
                        r = sub.loc[sub["ion_score"].idxmax()]
                        if float(r["ion_score"]) >= covalent_tau:
                            best_cov[pid] = (str(r["ion_formula"]), float(r["ion_score"]))
            except Exception as e:   # the oracle must never crash the label step
                log(f"[cleanup] cluster covalent-check skipped: {e}")

    # phase 3: mark, with honest wording
    for pid, mz in targets:
        defect = mz - round(mz)
        if pid in best_cov:
            ion, s = best_cov[pid]
            note = (f"bromide cluster ion (defect {defect:+.3f}); reagent-adduct "
                    f"reading preferred over degenerate di-bromo organic "
                    f"(ion {ion} also fits, score {s:.2f})")
        elif best_cov or (score_fn is not None or client is not None):
            # the oracle ran and found no covalent tie for THIS peak
            note = (f"bromide cluster ion (defect {defect:+.3f}, multi-Br reagent "
                    "region; no covalent reading scores above threshold)")
        else:
            # offline: no oracle was consulted -> make no claim about covalency
            note = (f"bromide cluster ion (defect {defect:+.3f}, multi-Br reagent "
                    "region)")
        # known formula -> assigned: record the ion formula when we have one
        L.mark_reagent(ledger, pid, note,
                       ion_formula=(best_cov[pid][0] if pid in best_cov else None))
    log(f"[cleanup] labelled {len(targets)} bromide reagent-cluster ions "
        f"({len(best_cov)} with a degenerate covalent tie)")
    return {"labelled": len(targets), "covalent_ties": len(best_cov)}


# ---- 3. isotope-confirmed low-complexity recovery ------------------------
# CHO neutrals + (optionally) a covalent halogen the isotope pattern can confirm.
# We deliberately do NOT enumerate N or S in the recovery: the only corroboration
# the recovery has is the halogen isotope envelope, and that confirms the adduct
# Br -- it says nothing about N or S. A borderline-score CHON/CHOS fit is exactly
# the degenerate class the review flagged (C18H14N2O3S: N2 + S, 34S absent, beaten
# by other formulas). So an uncorroborated heteroatom is not recoverable here;
# the confident CHON/CHOS are committed in the main passes, not in cleanup.
RECOVERY_BOX = "C0-20 H0-36 O0-12 Cl0-1 Br0-2"
RECOVERY_ADDUCTS = ("[M-H]-", "[M+Br]-")


def _measure_pattern(near, mz, h):
    """(r2, r4) = measured M+2/M0 and M+4/M0 intensity ratios."""
    p2, p4 = near(mz + BR), near(mz + 2 * BR)
    r2 = (p2[1] / h) if (p2 is not None and h > 0) else 0.0
    r4 = (p4[1] / h) if (p4 is not None and h > 0) else 0.0
    return r2, r4


def _pattern_ok(nbr_ion, ncl_ion, r2, r4):
    """Does the measured (r2, r4) confirm the candidate's halogen count?
    Recovery REQUIRES a halogen handle, so a (0 Br, 0 Cl) ion is rejected here
    (it has no independent corroboration for the relaxed score floor)."""
    if nbr_ion == 1 and ncl_ion == 0:
        return 0.78 <= r2 <= 1.20            # 1 Br: M+2/M0 ~ 0.97
    if nbr_ion == 2 and ncl_ion == 0:
        return 1.55 <= r2 <= 2.35 and 0.55 <= r4 <= 1.35   # 2 Br: 1.95 / 0.95
    if nbr_ion == 0 and ncl_ion == 1:
        return 0.20 <= r2 <= 0.48            # 1 Cl: M+2/M0 ~ 0.32
    if nbr_ion == 1 and ncl_ion == 1:
        return 1.10 <= r2 <= 1.55            # BrCl: ~1.29
    return False


def recover_isotope_gated(client, sample_id, ledger, profile, cfg, *,
                          score_fn=None, score_floor: float = 0.65,
                          z_max: float = 2.5, log=print) -> dict:
    """Recover genuine molecules the score gate dropped, gated on a CONFIRMED
    halogen isotope pattern + a low-complexity clamp."""
    score_fn = score_fn or IO.score_candidates
    mu = getattr(cfg, "cal_mu", None)
    sigma = getattr(cfg, "cal_sigma", None) or 0.5
    pk = ledger.dropna(subset=["mz"]).sort_values("mz")
    mzs = pk["mz"].tolist(); hts = pk["height"].tolist()
    def near(t, tol=0.006):
        i = bisect.bisect_left(mzs, t - tol); best = None
        while i < len(mzs) and mzs[i] <= t + tol:
            if best is None or abs(mzs[i] - t) < abs(best[0] - t):
                best = (mzs[i], hts[i])
            i += 1
        return best

    un = ledger[ledger["role"] == L.ROLE_UNEXPLAINED]
    ranges = C.parse_ranges(RECOVERY_BOX)
    # enumerate per unexplained peak, union for ONE batched score
    per = {}
    allf = set()
    for _, p in un.iterrows():
        mz = float(p["mz"])
        cands = {f for f in C.candidates_for_peaks(
                    [mz], ranges, list(RECOVERY_ADDUCTS), ppm_tolerance=2.0,
                    mass_min=mz - 82, mass_max=mz + 2)
                 if X.filter_by_profile(f, profile)[0]
                 and D._het_types(C.parse_formula(f)) <= 2}
        if cands:
            per[p["peak_id"]] = (mz, float(p["height"]), cands)
            allf |= cands
    if not allf:
        log("[cleanup] recovery: no low-complexity candidates")
        return {"recovered": 0}
    fr = score_fn(client, sample_id, sorted(allf), allow_partial=True,
                  mechanism_ids=getattr(cfg, "mechanism_ids", None))
    fr = fr[fr["sample_peak_id"].notna() & (fr["sample_peak_intensity"] > 0)]

    recovered = 0
    for pid, (mz, h, _c) in per.items():
        if ledger.loc[ledger.peak_id == pid, "role"].iloc[0] != L.ROLE_UNEXPLAINED:
            continue
        sub = fr[(fr["sample_peak_mz"] - mz).abs() < 0.006].copy()
        if sub.empty:
            continue
        if mu is not None:
            sub = sub[((sub["ppm_error"] - mu) / sigma).abs() <= z_max]
        sub = sub[sub["ion_score"] >= score_floor]
        if sub.empty:
            continue
        r2, r4 = _measure_pattern(near, mz, h)
        sub = sub.sort_values("ion_score", ascending=False)
        chosen = None
        for _, r in sub.iterrows():
            ic = C.parse_formula(r["ion_formula"])
            if _pattern_ok(ic.get("Br", 0), ic.get("Cl", 0), r2, r4):
                chosen = r
                break
        if chosen is None:
            continue
        # parse neutral/adduct back out of the candidate set for this peak,
        # preferring the reagent-adduct reading (fewest halogens in the neutral)
        neutral, adduct = _decompose(chosen["ion_formula"], _c)
        if neutral is None:
            continue
        # H1/H2: the recovery's only corroboration is the halogen isotope (it
        # confirms the adduct Br, not N or S). A CHON/CHOS neutral here is an
        # uncorroborated heteroatom fit -> refuse it (C18H14N2O3S: N2+S, 34S
        # absent). Only CHO + isotope-confirmed covalent halogen is recoverable.
        nc = C.parse_formula(neutral)
        if nc.get("N", 0) > 0 or nc.get("S", 0) > 0 or nc.get("P", 0) > 0:
            continue
        ppm = float(chosen["ppm_error"])
        z = (ppm - mu) / sigma if mu is not None else 0.0
        L.commit_assignment(
            ledger, pid, neutral_formula=neutral, adduct=adduct,
            ion_formula=str(chosen["ion_formula"]),
            ion_score=float(chosen["ion_score"]),
            compound_score=float(chosen["ion_score"]),
            ppm_error=ppm, pass_no=7, method="cleanup:iso-recovery",
            confidence="Good (recovered)",
            commentary=(f"Recovery: halogen isotope pattern confirms the ion "
                        f"(M+2/M0={r2:.2f}); score {chosen['ion_score']:.2f}, "
                        f"z={z:+.1f} within calibrated accuracy."))
        recovered += 1
    log(f"[cleanup] recovered {recovered} isotope-confirmed molecules")
    return {"recovered": recovered}


def _decompose(ion_formula, candidate_set):
    """Map a scored ion back to (neutral, adduct), PREFERRING the reading with
    the fewest halogens in the neutral. The reagent halogen's ion isotope cannot
    prove it is covalent (C5H10O3.Br- == C5H11BrO3 [M-H]-, same ion), and the
    pipeline policy is the reagent-adduct reading -- so a tie goes to the
    Br-in-the-adduct interpretation, deterministically."""
    target = str(ion_formula).rstrip("-")
    best = None   # (neutral_halogen_count, neutral, adduct)
    for f in sorted(candidate_set):
        for a in RECOVERY_ADDUCTS:
            if _ion_str(f, a) == target:
                cnt = C.parse_formula(f)
                nh = cnt.get("Br", 0) + cnt.get("Cl", 0)
                if best is None or nh < best[0]:
                    best = (nh, f, a)
    return (best[1], best[2]) if best else (None, None)


def _ion_str(neutral, adduct):
    cnt = dict(C.parse_formula(neutral))
    for sign, tok in __import__("re").findall(r"([+-])\^?([A-Za-z0-9]+)",
                                              adduct.split("]")[0][2:]):
        for el, k in C.parse_formula(tok).items():
            cnt[el] = cnt.get(el, 0) + (k if sign == "+" else -k)
    cnt = {k: v for k, v in cnt.items() if v}
    return C.format_formula(cnt)


def reclaim_satellites(ledger: pd.DataFrame, *, ppm: float = 6.0, log=print) -> dict:
    """Final sweep: attach an UNEXPLAINED peak that is a clean 13C / 81Br / 37Cl
    satellite of an assigned M0 (parent carries that element AND the intensity
    ratio is physically consistent) as an iso_child. Catches the few satellites
    the envelope passes leak. Touches only unexplained rows -- it can never demote
    a real M0. The 13C gate is carbon-count-aware (ratio ~ nC*1.1%), so it cannot
    mis-grab an unrelated neighbour at +1.003 Da."""
    from . import isotopes as ISO
    m0 = ledger[ledger["role"] == L.ROLE_M0].dropna(subset=["mz"]).sort_values("mz")
    if not len(m0):
        return {"reclaimed": 0}
    mz = m0["mz"].to_numpy(); pid = m0["peak_id"].to_numpy()
    ionf = m0["ion_formula"].astype(str).to_numpy(); ph = m0["height"].to_numpy()
    DELT = [(ISO.D_13C, "13C", "C"), (ISO.D_81BR, "81Br", "Br"), (ISO.D_37CL, "37Cl", "Cl")]
    n = 0
    for i in ledger.index[ledger["role"] == L.ROLE_UNEXPLAINED]:
        cmz = ledger.at[i, "mz"]; chh = ledger.at[i, "height"]
        if pd.isna(cmz):
            continue
        cmz = float(cmz); chh = float(chh) if pd.notna(chh) else 0.0
        for d, label, el in DELT:
            t = cmz - d; tol = cmz * ppm * 1e-6
            j = bisect.bisect_left(mz, t - tol); best = None
            while j < len(mz) and mz[j] <= t + tol:
                if best is None or abs(mz[j] - t) < abs(mz[best] - t):
                    best = j
                j += 1
            if best is None:
                continue
            cnt = C.parse_formula(ionf[best]); nel = cnt.get(el, 0)
            if nel < 1 or ph[best] <= 0:
                continue
            ratio = chh / ph[best]
            if el == "C":                      # carbon-count-aware 13C gate
                exp = nel * 0.0107
                ok = 0.3 * exp <= ratio <= 2.5 * exp and ratio < 1.0
            elif el == "Br":                   # ~0.97 per Br adduct/atom
                ok = 0.55 <= ratio <= 1.4 * nel
            else:                              # 37Cl ~0.32 per Cl
                ok = 0.18 <= ratio <= 0.5 * nel
            if not ok:
                continue
            try:
                L.attach_isotopologue(ledger, ledger.at[i, "peak_id"], pid[best],
                                      iso_label=label)
                ledger.at[i, "commentary"] = (
                    f"reclaimed {label} satellite of {ionf[best]} "
                    f"(ratio {ratio:.2f}); envelope-pass leak")
                n += 1
            except L.LedgerError:
                pass
            break
    log(f"[cleanup] reclaimed {n} leaked isotopologue satellites")
    return {"reclaimed": n}


def reclaim_envelope_tails(ledger: pd.DataFrame, *, ppm: float = 6.0, log=print) -> dict:
    """Attach the DEEP multi-halogen isotope envelope (k>=2 of ³⁷Cl / ⁸¹Br) of a
    poly-halogen M0 to leaked unexplained peaks. reclaim_satellites only does the
    single k=1 step, so a chlorinated paraffin's M+4/M+6/... ³⁷Cl tail -- which for
    many Cl is BRIGHTER than M0 -- leaks into 'unexplained' (the dominant satellite
    leak measured on the ¹⁵NO₃⁻ batch). Walk k=2..nX from each M0 and accept an
    unexplained peak whose intensity ratio matches the binomial C(nX,k)(p/q)^k.
    Touches only unexplained rows; the ratio gate (and exact +k·Δ mass) prevent
    grabbing an unrelated neighbour. p/q: ³⁷Cl 0.3199, ⁸¹Br 0.9728.

    KNOWN LIMITATION (do not trust this to clear real tails): on real batches this
    has been observed to be a no-op — the deep-tail leak it targets is, in
    practice, already absorbed upstream (reclaim_satellites k=1 + the isotope-locked
    known-species/CP recovery), so by the time this runs there are typically no
    leaked unexplained tails left for the binomial gate to match. The unit test
    exercises only the synthetic poly-³⁷Cl case. Kept (harmless, additive — touches
    only unexplained rows) pending a live-data re-evaluation of the gate; see
    docs/ROADMAP.md."""
    from math import comb
    from . import isotopes as ISO
    m0 = ledger[ledger["role"] == L.ROLE_M0].dropna(subset=["mz"]).sort_values("mz")
    if not len(m0):
        return {"tails": 0}
    mz = m0["mz"].to_numpy(); pid = m0["peak_id"].to_numpy()
    ionf = m0["ion_formula"].astype(str).to_numpy(); ph = m0["height"].to_numpy()
    STEPS = [(ISO.D_37CL, "37Cl", "Cl", 0.3199), (ISO.D_81BR, "81Br", "Br", 0.9728)]
    n = 0
    for i in ledger.index[ledger["role"] == L.ROLE_UNEXPLAINED]:
        cmz = ledger.at[i, "mz"]
        if pd.isna(cmz):
            continue
        cmz = float(cmz); chh = float(ledger.at[i, "height"]) if pd.notna(ledger.at[i, "height"]) else 0.0
        done = False
        for D, label, el, pq in STEPS:
            for k in range(2, 11):
                t = cmz - k * D; tol = cmz * ppm * 1e-6
                j = bisect.bisect_left(mz, t - tol); best = None
                while j < len(mz) and mz[j] <= t + tol:
                    if best is None or abs(mz[j] - t) < abs(mz[best] - t):
                        best = j
                    j += 1
                if best is None:
                    continue
                nX = C.parse_formula(ionf[best]).get(el, 0)
                if k > nX or ph[best] <= 0:
                    continue
                exp = comb(nX, k) * (pq ** k)            # binomial intensity vs M0
                ratio = chh / ph[best]
                if 0.35 * exp <= ratio <= 2.8 * exp:
                    try:
                        L.attach_isotopologue(ledger, ledger.at[i, "peak_id"], pid[best],
                                              iso_label=f"{k}x{label}")
                        ledger.at[i, "commentary"] = (
                            f"reclaimed {k}x{label} envelope tail of {ionf[best]} "
                            f"(ratio {ratio:.2f} vs binomial {exp:.2f})")
                        n += 1; done = True
                    except L.LedgerError:
                        pass
                    break
            if done:
                break
    log(f"[cleanup] reclaimed {n} deep halogen-envelope tails (k>=2)")
    return {"tails": n}


F_DEMOTE_MIN = 4   # F count above which an unanchored, non-PFCA fit is an F-monster


def _is_pfca(cnt: dict) -> bool:
    """Perfluorocarboxylic acid CnHF(2n-1)O2 (TFA, PFPrA, ... PFOA)."""
    nC = cnt.get("C", 0)
    return (nC >= 2 and cnt.get("H", 0) == 1 and cnt.get("O", 0) == 2
            and cnt.get("F", 0) == 2 * nC - 1)


def demote_unconfirmed_fluorine(ledger: pd.DataFrame, *, f_min: int = F_DEMOTE_MIN,
                                log=print) -> dict:
    """Demote M0 assignments resting on UNCONFIRMED fluorine. ¹⁹F is 100%
    monoisotopic, so a high-F formula has NO isotope twin to corroborate it; when it
    is also not a known PFCA and carries no Cl/Br/S anchor (whose ³⁷Cl/⁸¹Br/³⁴S
    pattern WOULD confirm a heavy atom), a high-F fit is almost always a mass
    coincidence ('F-monster', e.g. C11H6F16). Demote Identified->Candidate and flag
    below_assignability. Real PFCAs and isotope-anchored F species are kept."""
    n = 0
    has_ba = "below_assignability" in ledger.columns
    target = (ledger.index[ledger["role"] == L.ROLE_M0]
              if "role" in ledger.columns else ledger.index)   # merged ledger has no role col
    for i in target:
        cnt = C.parse_formula(str(ledger.at[i, "neutral_formula"] or ""))
        if cnt.get("F", 0) < f_min:
            continue
        if _is_pfca(cnt) or cnt.get("Cl", 0) or cnt.get("Br", 0) or cnt.get("S", 0):
            continue
        if str(ledger.at[i, "tier"]) == "Identified":
            ledger.at[i, "tier"] = "Candidate"
        if has_ba:
            ledger.at[i, "below_assignability"] = True
        if "commentary" in ledger.columns:
            note = (f"unconfirmed fluorine (F{cnt.get('F', 0)}; ¹⁹F monoisotopic, no "
                    "Cl/Br/S anchor, not a PFCA) -- likely mass coincidence")
            prev = str(ledger.at[i, "commentary"] or "")
            ledger.at[i, "commentary"] = (prev + "; " + note) if prev and prev != "nan" else note
        n += 1
    log(f"[cleanup] demoted {n} unconfirmed-fluorine M0 (F>={f_min}, F-monsters)")
    return {"f_demoted": n}


def run_cleanup(client, sample_id, ledger, profile, cfg, *, log=print) -> dict:
    """Orchestrate the cleanup steps (recovery first, so a recovered molecule
    isn't then mislabelled a cluster/artifact; satellite reclaim last)."""
    rec = recover_isotope_gated(client, sample_id, ledger, profile, cfg, log=log)
    clu = label_bromide_clusters(ledger, client, sample_id, log=log)
    art = flag_ringing_artifacts(ledger, log=log)
    sat = reclaim_satellites(ledger, log=log)
    tails = reclaim_envelope_tails(ledger, log=log)   # deep poly-halogen envelope (k>=2)
    # NB: demote_unconfirmed_fluorine is NOT called here -- it must run AFTER
    # tiers.apply_tiers (which recomputes tier and would re-promote the F-monster).
    # assign.run calls it post-tiering.
    return {"recovered": rec["recovered"], "clusters": clu["labelled"],
            "cluster_covalent_ties": clu.get("covalent_ties", 0),
            "artifacts": art["flagged"], "reclaimed_satellites": sat["reclaimed"],
            "envelope_tails": tails["tails"]}


def prefer_amine_over_ammonium(ledger: pd.DataFrame, *, ts_peaks=None, r_min: float = 0.7,
                               log=print) -> dict:
    """Uronium / positive urea-CIMS: a [M+NH4]+ adduct of a CHO neutral X is mass-
    AND isotope-identical to [M+H]+ of the amine X+NH3 (the SAME ion formula), so
    the data cannot distinguish them. In an N-rich urea source the protonated amine
    is the simpler explanation than an ammonium side-adduct, so RE-READ each
    [M+NH4]+ assignment as [M+H]+ of X+NH3 -- UNLESS:

      * X is CORROBORATED. Without a time series (`ts_peaks=None`) corroboration is
        PRESENCE-based: X is also assigned as [M+H]+ or [M+(CH4N2O)H]+. With a batch
        time series it is the stronger CO-VARIATION test: the NH4 trace must
        correlate (log Pearson r >= r_min) with the [M+H]+ OR urea-cluster trace of
        the SAME m/z -- a faint NH4 peak that does not track the parent is NOT
        corroboration; or
      * the amine X+NH3 is valence-impossible (saturated X -> negative DBE); the
        NH4 adduct is then the only valid reading and is FORCED/kept.

    Only neutral_formula + adduct change (the ion, m/z, score, tier, ppm are
    identical). Relabels in place; returns a summary.
    """
    if not {"neutral_formula", "adduct"} <= set(ledger.columns):
        return {"relabeled": 0, "kept_corroborated": 0, "forced_nh4": 0}
    is_m0 = (ledger["role"] == L.ROLE_M0) if "role" in ledger.columns \
        else pd.Series(True, index=ledger.index)

    keeps = _make_corroboration_test(ledger, ts_peaks, r_min)

    relabeled = kept = forced = 0
    for i in ledger.index[is_m0 & (ledger["adduct"] == "[M+NH4]+")]:
        X = str(ledger.at[i, "neutral_formula"] or "")
        if not X or X == "nan":
            continue
        if keeps(X):                                      # corroborated -> keep NH4
            kept += 1
            continue
        cnt = C.parse_formula(X)
        cnt["N"] = cnt.get("N", 0) + 1
        cnt["H"] = cnt.get("H", 0) + 3
        ok, _ = C.dbe_ok(cnt)
        if not ok:                                        # amine impossible -> NH4 forced
            forced += 1
            continue
        ledger.at[i, "neutral_formula"] = C.format_formula(cnt)
        ledger.at[i, "adduct"] = "[M+H]+"
        if "tier_reason" in ledger.columns:
            ledger.at[i, "tier_reason"] = (
                str(ledger.at[i, "tier_reason"] or "")
                + " | NH4-adduct re-read as protonated +NH3 amine (uronium parsimony)").strip(" |")
        relabeled += 1
    mode = f"co-variation r>={r_min}" if ts_peaks is not None else "presence"
    log(f"[uronium-amine] ({mode}) {relabeled} [M+NH4]+ re-read as [M+H]+ amine; "
        f"kept {kept} corroborated NH4 adducts; {forced} forced (no valid amine)")
    return {"relabeled": relabeled, "kept_corroborated": kept, "forced_nh4": forced}


def _make_corroboration_test(ledger, ts_peaks, r_min):
    """Return keeps(neutral) -> bool. Presence-based unless a time series is given,
    then NH4 is kept only when its trace co-varies (log r>=r_min) with the [M+H]+ or
    urea-cluster trace."""
    if ts_peaks is None:
        corrob = set(ledger.loc[ledger["adduct"].isin(["[M+H]+", "[M+(CH4N2O)H]+"]),
                                "neutral_formula"].dropna().astype(str))
        return lambda X: X in corrob

    import numpy as np
    from . import timeseries as TS
    mat, bin_mz = TS.build_matrix(ts_peaks)
    bm = bin_mz.sort_values(); arr = bm.to_numpy(); idx = bm.index.to_numpy()

    def _logtrace(neutral, adduct):
        try:
            mz = C.ion_mz(neutral, adduct)
        except Exception:
            return None
        i = np.searchsorted(arr, mz); best = None
        for j in (i - 1, i):
            if 0 <= j < len(arr):
                ppm = abs(arr[j] - mz) / mz * 1e6
                if ppm <= 8 and (best is None or ppm < best[1]):
                    best = (idx[j], ppm)
        if best is None:
            return None
        return np.log10(mat[best[0]].clip(lower=1).to_numpy())

    def keeps(X):
        nh4 = _logtrace(X, "[M+NH4]+")
        if nh4 is None:
            return False
        for ref in ("[M+H]+", "[M+(CH4N2O)H]+"):
            t = _logtrace(X, ref)
            if t is None:
                continue
            ok = np.isfinite(nh4) & np.isfinite(t)
            if ok.sum() >= 6 and np.corrcoef(nh4[ok], t[ok])[0, 1] >= r_min:
                return True
        return False
    return keeps

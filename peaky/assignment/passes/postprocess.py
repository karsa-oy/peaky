"""passes.postprocess — split from the former passes.py monolith."""

from __future__ import annotations


import pandas as pd

from peaky.chem import chemistry as C
from peaky.chem import isotopes as ISO
from peaky.assignment import ledger as L


from .config import PassConfig
from .core import z_of

__all__ = [
    "_D13C",
    "_DBR",
    "_R13C",
    "_R81BR",
    "_D29SI",
    "_R29SI",
    "SI_M1_MIN_FRAC",
    "_peak_near",
    "_si_m1_consistent",
    "complete_isotope_envelopes",
    "_M1_RATIO",
    "detect_composites",
    "split_composites",
    "demote_carbon_inconsistent",
    "demote_massgate_monsters",
    "audit_isotopes",
    "audit_mass_gate",
]


_D13C = 1.0033548  # 13C - 12C


_DBR = 1.9979535  # 81Br - 79Br


_R13C = 0.0107  # 13C natural abundance per carbon


_R81BR = 0.9728  # 81Br/79Br abundance ratio


_D29SI = 0.999568  # 29Si - 28Si (the Si M+1)


_R29SI = 0.0468  # 29Si natural abundance per Si


SI_M1_MIN_FRAC = 0.6  # observed (M+1)/(M0) must be >= this * predicted Si M+1


def _peak_near(mzs: "pd.Series", target: float, ppm: float = 5.0):
    """Index of the closest ledger peak within ppm of target, else None."""
    tol = target * ppm * 1e-6
    d = (mzs - target).abs()
    i = d.idxmin()
    return i if d.loc[i] <= tol else None


def _si_m1_consistent(
    ledger: pd.DataFrame, m0_mz: float, m0_h: float, n_si: int, n_c: int
) -> bool:
    """Is the 29Si M+1 intensity consistent with the CLAIMED Si count? A Si_k species'
    M+1 is dominated by 29Si (n_si*4.68%) plus 13C (n_c*1.07%); when the observed
    (M+1)/(M0) ratio is far below that, the Si count is over-claimed -- a high-O
    organic with only a 13C M+1 masquerading as a siloxane (the C10H18O11 vs
    C8H26O5Si4 degeneracy at m/z 393). True = OK to commit; False = skip."""
    if n_si <= 0 or m0_h <= 0:
        return True
    pred = n_si * _R29SI + n_c * _R13C
    if pred <= 0:
        return True
    j = _peak_near(ledger["mz"], m0_mz + _D29SI, ppm=15.0)
    obs = (
        float(ledger.at[j, "height"]) / m0_h
        if j is not None and pd.notna(ledger.at[j, "height"])
        else 0.0
    )
    return obs >= SI_M1_MIN_FRAC * pred


def complete_isotope_envelopes(
    ledger: pd.DataFrame,
    cfg: PassConfig,
    *,
    min_rel: float = 0.06,
    ppm: float = 12.0,
    log=print,
) -> dict:
    """Claim the FULL predicted isotope envelope (M+1/M+2/M+4...) of every
    committed M0, so multi-isotope species (Si-rich silanediols, multi-Br/Cl
    compounds) don't leak satellites into the residual.

    Two actions per predicted satellite line:
      * an UNEXPLAINED peak at the right mass + consistent intensity is attached
        as an iso_child (the envelope was incompletely claimed by the server);
      * a WEAK committed M0 (not High/Assigned, not locked) that is really a
        parent's satellite is DISPLACED into the iso_child role -- this is the
        393/395 silanediol bug, where the Si4+Br M+2 at 395 got mis-assigned a
        Cl-F-S formula because its M+4/M+2 ratio (~0.26) mimicked a Cl doublet.

    Processed parent-before-satellite (ascending m/z): a satellite is always
    heavier than its parent, so the true parent claims it first. The pattern is
    formula-specific (a CHO ion predicts only 13C, plus a 13C2 M+2 above ~28 C;
    halogen/Si M+2 lines need the actual heteroatom), and an intensity-
    consistency gate is the discriminator against coincidental neighbours -- a
    real satellite sits at the predicted height, an independent compound does
    not. The match tolerance is tight for M+1/M+2 (to separate 13C from 29Si,
    3.8 mDa apart) and looser for the multi-isotope M+4+ centroid."""
    out = {"attached": 0, "displaced": 0}
    mzs = ledger["mz"]
    order = ledger[ledger["role"] == L.ROLE_M0].sort_values("mz")["peak_id"].tolist()
    for pid in order:
        idx = ledger.index[ledger["peak_id"] == pid]
        if not len(idx):
            continue
        i = idx[0]
        if str(ledger.at[i, "role"]) != L.ROLE_M0:
            continue  # displaced by an earlier parent
        ionf = ledger.at[i, "ion_formula"]
        if ionf is pd.NA or pd.isna(ionf) or not str(ionf).strip():
            continue
        pmz = float(ledger.at[i, "mz"])
        ph = float(ledger.at[i, "height"])
        if not (ph > 0):
            continue
        try:
            # max_shift 12: keep the M+7/M+8 envelope of 4+ heavy-halogen ions
            # (a Br4 M+8 is ~0.9x M0) instead of leaking it into the residual
            pattern = ISO.isotope_pattern(str(ionf), min_rel=min_rel, max_shift=12.0)
        except Exception:
            continue
        for dmass, rel, label in pattern:
            # tolerance is shift-aware: M+1/M+2 must separate 13C (+1.0034) from
            # 29Si (+0.9996) -- 3.8 mDa apart -- so they use a tight window; the
            # multi-isotope M+4+ centroid is approximate and uses the loose one
            line_ppm = 5.0 if dmass < 2.5 else ppm
            j = _peak_near(mzs, pmz + dmass, ppm=line_ppm)
            if j is None:
                continue
            tpid = ledger.at[j, "peak_id"]
            if tpid == pid:
                continue
            th = float(ledger.at[j, "height"])
            if not (th > 0):
                continue
            ratio = th / (ph * rel)
            score = min(ratio, 1.0 / ratio) if ratio > 0 else 0.0
            role_j = str(ledger.at[j, "role"])
            if role_j == L.ROLE_UNEXPLAINED:
                if 0.3 <= ratio <= 3.5:
                    try:
                        L.attach_isotopologue(
                            ledger, tpid, pid, iso_label=label, iso_match_score=score
                        )
                        out["attached"] += 1
                    except L.LedgerError:
                        pass
            elif role_j == L.ROLE_M0:
                if bool(ledger.at[j, "locked"]):
                    continue
                conf_j = str(ledger.at[j, "confidence"])
                sc_j = ledger.at[j, "ion_score"]
                weak_score = pd.isna(sc_j) or float(sc_j) < cfg.tau_high
                # only displace a WEAK victim, on a tight intensity match. The
                # tier column is NA here (tiers run later), so protect on
                # CONFIDENCE + standalone SCORE instead: High-confidence or
                # near-High-scoring fits are real compounds, not satellites.
                if (
                    not conf_j.startswith("High")
                    and weak_score
                    and 0.45 <= ratio <= 2.2
                ):
                    try:
                        L.displace_to_isotopologue(
                            ledger, tpid, pid, iso_label=label, iso_match_score=score
                        )
                        out["displaced"] += 1
                    except L.LedgerError:
                        pass
    if out["attached"] or out["displaced"]:
        log(
            f"[iso-envelope] attached {out['attached']} unclaimed satellites, "
            f"displaced {out['displaced']} mis-assigned satellites onto their "
            f"true parents"
        )
    return out


_M1_RATIO = {
    "C": 0.0107,
    "Si": 0.0508,
    "N": 0.003653,
    "S": 0.007896,
    "O": 0.000381,
    "H": 0.000115,
}


def detect_composites(
    ledger: pd.DataFrame,
    cfg: PassConfig,
    *,
    min_m1_rel: float = 0.06,
    excess_frac: float = 0.25,
    min_excess: float = 400.0,
    ppm: float = 8.0,
    log=print,
) -> dict:
    """Flag committed M0 peaks that are UNRESOLVED COMPOSITES -- their M0 (and
    M+2/M+4) intensity exceeds what their own M+1 satellite implies, because a
    coincident co-eluting compound shares the m/z.

    The discriminator is the even/odd isotope split: the M+1 region (13C, 29Si,
    15N -- all halogen-free) scales ONLY with the assigned compound, so it gives
    the assigned compound's true intensity S = M+1_obs / M+1_predicted. If the
    observed M0 markedly exceeds S, the excess is a co-component, and its
    halogen content is read off the EVEN-shift residual (M+2/M0, M+4/M+2 ~ Br /
    BrCl / Br2). This is the silanediol case the isotope-pattern 'mismatch'
    flagged: C8H26O5Si4 (Si4) at 393 sits on a ~45% BrCl compound -- formula and
    prediction are both correct; the peak is mixed. n=2 (clean) is not flagged.

    Flags only (does not demote): the assigned compound IS present; the note
    records the co-component fraction + halogen guess so the report is honest."""
    out = {"flagged": 0}
    mzs = ledger["mz"]
    if "composite_note" not in ledger.columns:
        ledger["composite_note"] = pd.Series(pd.NA, index=ledger.index, dtype="object")

    def _sum_window(lo, hi):
        m = (mzs >= lo) & (mzs <= hi)
        return float(ledger.loc[m, "height"].sum(skipna=True))

    for i, r in ledger[ledger["role"] == L.ROLE_M0].iterrows():
        ionf = r["ion_formula"]
        if ionf is pd.NA or pd.isna(ionf) or not str(ionf).strip():
            continue
        cnt = C.parse_formula(str(ionf))
        m1_rel = sum(_M1_RATIO.get(el, 0.0) * n for el, n in cnt.items())
        if m1_rel < min_m1_rel:
            continue  # too few C/Si to diagnose a composite
        m0 = float(r["mz"])
        h0 = float(r["height"])
        if not (h0 > 0):
            continue
        # observed M+1 region: 13C(+1.0034) + 29Si(+0.9996) + 15N(+0.997)
        h1 = _sum_window(m0 + 0.9940, m0 + 1.0070)
        if h1 <= 0:
            continue
        s_assigned = h1 / m1_rel  # implied true intensity of the M0 owner
        excess = h0 - s_assigned
        if excess < min_excess or excess / h0 < excess_frac:
            continue
        # characterise the co-component via the even-shift residual. SUM (never
        # overwrite) lines that round to the same integer shift, else the big
        # 81Br M+2 (rel ~1.14) is clobbered by the tiny 13C2 line at +2.007.
        pat: dict[int, float] = {}
        try:
            for d, rel, _ in ISO.isotope_pattern(str(ionf), min_rel=0.01):
                pat[round(d)] = pat.get(round(d), 0.0) + rel
        except Exception:
            pat = {}
        h2 = _sum_window(m0 + 1.992, m0 + 2.004)
        h4 = _sum_window(m0 + 3.990, m0 + 4.002)
        x2 = h2 - s_assigned * pat.get(2, 0.0)  # co-component M+2
        x4 = h4 - s_assigned * pat.get(4, 0.0)  # co-component M+4
        hal = "unknown"
        if x2 > min_excess:
            r2 = x2 / excess
            r4 = x4 / x2 if x2 > 0 else 0.0
            if r2 >= 1.6:
                hal = "Br2"
            elif r2 >= 1.15 and r4 >= 0.22:
                hal = "BrCl"
            elif r2 >= 0.7:
                hal = "Br"
            elif 0.22 <= r2 <= 0.45:
                hal = "Cl"
        ledger.at[i, "composite_note"] = (
            f"composite: ~{100 * excess / h0:.0f}% co-eluting {hal} component "
            f"(~{excess:.0f} cps); M+1 implies {str(r['neutral_formula'])} "
            f"= ~{s_assigned:.0f} of {h0:.0f} cps"
        )
        # structured fields for the de-blending step (split_composites)
        ledger.at[i, "assigned_fraction"] = max(0.0, min(1.0, s_assigned / h0))
        ledger.at[i, "co_height"] = excess
        ledger.at[i, "co_halogen"] = hal
        out["flagged"] += 1
    if out["flagged"]:
        log(
            f"[composite] flagged {out['flagged']} unresolved composite peaks "
            f"(M0 inflated beyond the M+1-implied owner intensity)"
        )
    return out


def split_composites(
    ledger: pd.DataFrame, cfg: PassConfig, *, log=print
) -> pd.DataFrame:
    """De-blend the peaks `detect_composites` flagged: the owner keeps its
    `assigned_fraction` of the measured height, and a SYNTHETIC sub-peak
    ('<host>.2', same m/z) is created carrying the co-eluting compound's share
    (co_height) plus its halogen guess. The sub-peak is a characterised residual
    (role unexplained, synthetic=True, host_peak_id->host) -- a target for a
    later constrained match that NAMES the co-component. Signal is conserved:
    effective = height*assigned_fraction, so host + sub-peak sum to the original
    height (the host's measured `height` is never altered).

    Returns the (possibly grown) ledger -- new synthetic rows are appended, so
    the caller must rebind: `led = split_composites(led, cfg)`."""
    if "co_height" not in ledger.columns:
        return ledger
    syn = ledger["synthetic"].fillna(False).astype(bool)
    hosts = ledger[(ledger["role"] == L.ROLE_M0) & ledger["co_height"].notna() & (~syn)]
    existing = set(ledger["peak_id"])
    new_rows = []
    for _, r in hosts.iterrows():
        co_h = float(r["co_height"])
        sub_id = f"{r['peak_id']}.2"
        if co_h < 1.0 or sub_id in existing:
            continue
        hal = str(r["co_halogen"])
        row = {c: pd.NA for c in ledger.columns}
        row.update(
            {
                "peak_id": sub_id,
                "mz": float(r["mz"]),
                "height": co_h,
                "area": float("nan"),
                "role": L.ROLE_UNEXPLAINED,
                "synthetic": True,
                "host_peak_id": r["peak_id"],
                "assigned_fraction": 1.0,
                "locked": False,
                "co_height": float("nan"),
                "commentary": (
                    f"co-eluting {hal} component split from peak "
                    f"{r['mz']:.4f} (~{co_h:.0f} cps); host owner "
                    f"{str(r['neutral_formula'])} keeps "
                    f"{100 * float(r['assigned_fraction']):.0f}%"
                ),
            }
        )
        new_rows.append(row)
        existing.add(sub_id)
    if not new_rows:
        return ledger
    add = pd.DataFrame(new_rows)[list(ledger.columns)]
    log(
        f"[composite] split {len(new_rows)} composite peaks into fractional "
        f"sub-peaks (owner keeps assigned_fraction; co-component -> '<id>.2')"
    )
    return pd.concat([ledger, add], ignore_index=True)


def demote_carbon_inconsistent(
    ledger: pd.DataFrame, cfg: PassConfig, *, log=print
) -> int:
    """Clear committed M0s whose carbon count is contradicted by their 13C
    satellite -- the 'O15 monster' class. Run BEFORE pass 4 (not just in the
    end-of-run audit) so the freed bright peaks are re-offered the correct
    carbon-clamped interpretation. Without this, pass 1 grabs a lattice peak
    with a low-carbon CHON mass-fit (e.g. C11H10N2O15 on the 4.7k-cps 409.0015,
    whose 13C satellite measures ~C16), pass 4 skips it because it is no longer
    unexplained, and the audit only clears it after every pass has run -- too
    late to re-assign as the di-bromide SOA cluster (C15H22O3 [M+HBr+Br]-)."""
    mzs = ledger["mz"]
    n = 0
    for _, r in ledger[
        (ledger["role"] == L.ROLE_M0) & ~ledger["locked"].astype(bool)
    ].iterrows():
        cnt = C.parse_formula(str(r["ion_formula"]))
        n_c = cnt.get("C", 0)
        if n_c < 8:
            continue
        # the carbon clamp reads the M+1 region as 13C only; a Si-bearing formula
        # has a 29Si M+1 (4.7%/Si) far larger than 13C, so the measured "13C
        # ratio" over-estimates carbon and would wrongly clear a real siloxane.
        # Skip them -- their carbon is corroborated by the Si-isotope envelope.
        if cnt.get("Si", 0) > 0:
            continue
        # measure carbon from a committed 13C child, else an unclaimed satellite
        k = ledger[
            (ledger["role"] == L.ROLE_ISO)
            & (ledger["parent_peak_id"] == r["peak_id"])
            & (ledger["iso_label"].astype(str) == "13C")
        ]
        if len(k):
            h_sat = float(k.iloc[0]["height"])
        else:
            j = _peak_near(mzs, float(r["mz"]) + _D13C)
            if j is None or ledger.at[j, "role"] != L.ROLE_UNEXPLAINED:
                continue
            h_sat = float(ledger.at[j, "height"])
        h0 = float(r["height"])
        if not (h0 > 0 and h_sat > 0):
            continue
        # only clamp on a RELIABLY-measured 13C satellite (>= the peak-detection
        # floor). Below it the ratio is noise: a genuine low-intensity M0 whose
        # weak/peak-picker-lost 13C reads as too-few-carbons would be falsely
        # cleared (the real [M+15NO3]- M0s at ~2k cps whose ~150 cps 13C sits
        # near the floor). The over-claim O-monster always has a BRIGHT 13C.
        if h_sat < cfg.height_cutoff:
            continue
        c_est = (h_sat / h0) / _R13C
        if abs(c_est - n_c) > max(2.5, 0.35 * n_c):
            try:
                L.clear_assignment(
                    ledger,
                    r["peak_id"],
                    reason=f"carbon-clamp (pre-pass-4): 13C ratio measures ~C"
                    f"{c_est:.0f}, formula claims C{n_c}",
                )
                n += 1
            except L.LedgerError:
                continue
    if n:
        log(
            f"[pre-pass4] carbon-clamp demoted {n} C-inconsistent M0 monsters "
            f"-> re-offered to the residual passes"
        )
    return n


def demote_massgate_monsters(
    ledger: pd.DataFrame, cfg: PassConfig, *, log=print
) -> int:
    """Clear pre-calibration M0s whose calibrated mass error is egregious
    (z > cal_z_pattern) BEFORE pass 4, the mass-gate twin of
    demote_carbon_inconsistent. Pass 1 grabs bright halogen-doublet peaks with
    high-O CHON mass-fits that the END mass-gate audit clears (e.g. C11H10N2O16
    on the 424.99 di-bromide peak, z=7.3) -- too late for pass 4 to re-claim.
    Only the clear monsters (z > pattern band) are cleared here; the 2..4-sigma
    tier is left for the end-of-run audit to keep this conservative."""
    if cfg.cal_mu is None:
        return 0
    n = 0
    for _, r in ledger[
        (ledger["role"] == L.ROLE_M0) & ~ledger["locked"].astype(bool)
    ].iterrows():
        z = z_of(r["ppm_error"], cfg)
        if z is not None and z > cfg.cal_z_pattern:
            try:
                L.clear_assignment(
                    ledger,
                    r["peak_id"],
                    reason=f"mass-gate (pre-pass-4): z={z:.1f} > {cfg.cal_z_pattern}",
                )
                n += 1
            except L.LedgerError:
                continue
    if n:
        log(
            f"[pre-pass4] mass-gate demoted {n} z>{cfg.cal_z_pattern} monsters "
            f"-> re-offered to the residual passes"
        )
    return n


def audit_isotopes(ledger: pd.DataFrame, cfg: PassConfig, *, log=print) -> dict:
    """Post-run isotope-physics audit. Validates committed M0s against what
    the isotope pattern REQUIRES, independent of match scores (v16 audit):

    1. Br-doublet repair: two M0s 1.99795 apart at ~1:1 height are one
       single-Br compound, not two formulas. If the lighter ion carries Br,
       the heavier peak becomes its 81Br child; if neither formula carries
       Br, both are wrong (the doublet proves Br) and both are cleared.
    2. 13C sweeper: attach the obvious unclaimed 13C satellite (right place,
       right magnitude) as evidence instead of leaving it unexplained.
    3. 13C carbon clamp: a committed 13C child measures the carbon count;
       a formula whose C is far outside it is wrong (C19 claimed, ~C11 seen).
    4. 13C completeness: a formula predicting a comfortably-visible 13C
       satellite that has NO peak at +1.0034 is wrong.
    """
    out = {
        "doublet_child": 0,
        "doublet_cleared": 0,
        "c13_attached": 0,
        "c13_clamp": 0,
        "c13_missing": 0,
    }
    mzs = ledger["mz"]

    # --- 1. Br-doublet repair over committed M0s ---
    m0 = ledger[
        (ledger["role"] == L.ROLE_M0) & ~ledger["locked"].astype(bool)
    ].sort_values("mz")
    rows = list(m0[["peak_id", "mz", "height", "ion_formula"]].itertuples(index=False))
    for a in range(len(rows)):
        lt = rows[a]
        for b in range(a + 1, len(rows)):
            hv = rows[b]
            d = hv.mz - lt.mz
            if d > _DBR + 0.005:
                break
            if abs(d - _DBR) > 0.004:
                continue
            hr = hv.height / lt.height
            if not (0.6 <= hr <= 1.45):
                continue
            try:
                if (
                    L.role_of(ledger, lt.peak_id) != L.ROLE_M0
                    or L.role_of(ledger, hv.peak_id) != L.ROLE_M0
                ):
                    continue
                n_br = C.parse_formula(str(lt.ion_formula)).get("Br", 0)
                if n_br >= 1:
                    # the lighter formula genuinely carries Br -> a ~1:1 twin
                    # 1.998 above IS its ⁸¹Br isotopologue (valid regardless of
                    # the reagent system).
                    L.clear_assignment(
                        ledger,
                        hv.peak_id,
                        reason=f"isotope audit: 81Br twin of "
                        f"{lt.mz:.4f} (ratio {hr:.2f})",
                    )
                    L.attach_isotopologue(
                        ledger, hv.peak_id, lt.peak_id, iso_label="81Br"
                    )
                    out["doublet_child"] += 1
                elif cfg.reagent_element == "Br":
                    # clear-both ONLY in Br-CIMS. There a ~1:1 1.998 doublet is
                    # strong evidence of an (unassigned) bromine, so two non-Br
                    # formulas are both wrong. With any OTHER reagent (e.g.
                    # ¹⁵N-nitrate) bromine is not in play: unrelated CHON
                    # compounds routinely sit ~1.998 apart at ~1:1, and the
                    # spacing is NOT halogen evidence -- clearing both destroys
                    # real M0s (54 genuine [M+¹⁵NO₃]⁻ M0s on the ¹⁵NO₃⁻ batch).
                    L.clear_assignment(
                        ledger,
                        lt.peak_id,
                        reason=f"isotope audit: Br doublet with {hv.mz:.4f} "
                        f"(ratio {hr:.2f}) but no Br in formula",
                    )
                    L.clear_assignment(
                        ledger,
                        hv.peak_id,
                        reason=f"isotope audit: Br doublet with {lt.mz:.4f} "
                        f"(ratio {hr:.2f}) but no Br in formula",
                    )
                    out["doublet_cleared"] += 2
            except L.LedgerError:
                continue

    # --- 2-4. 13C physics on every surviving M0 ---
    kids = ledger[ledger["role"] == L.ROLE_ISO]
    for _, r in ledger[ledger["role"] == L.ROLE_M0].iterrows():
        if bool(r["locked"]):
            continue
        n_c = C.parse_formula(str(r["ion_formula"])).get("C", 0)
        if n_c < 1:
            continue
        expected = float(r["height"]) * _R13C * n_c
        k = kids[
            (kids["parent_peak_id"] == r["peak_id"])
            & (kids["iso_label"].astype(str) == "13C")
        ]
        if not len(k):
            j = _peak_near(mzs, r["mz"] + _D13C)
            if (
                j is not None
                and ledger.at[j, "role"] == L.ROLE_UNEXPLAINED
                and expected > 0
                and 0.3 <= ledger.at[j, "height"] / expected <= 2.5
            ):
                try:
                    L.attach_isotopologue(
                        ledger, ledger.at[j, "peak_id"], r["peak_id"], iso_label="13C"
                    )
                    out["c13_attached"] += 1
                    k = ledger.loc[[j]]
                except L.LedgerError:
                    pass
        if len(k):
            h_sat = float(k.iloc[0]["height"])
            c_est = (h_sat / float(r["height"])) / _R13C
            # clamp ONLY on a reliably-measured 13C satellite (>= the detection
            # floor). A sub-floor 13C ratio is noise and under-reads carbon,
            # which would falsely clear genuine low-intensity M0s (the ~2k cps
            # [M+15NO3]- compounds whose ~150 cps 13C sits near the floor). The
            # over-claim O-monster case always carries a BRIGHT 13C, so it fires.
            if (
                n_c >= 8
                and h_sat >= cfg.height_cutoff
                and abs(c_est - n_c) > max(2.5, 0.35 * n_c)
            ):
                try:
                    L.clear_assignment(
                        ledger,
                        r["peak_id"],
                        reason=f"isotope audit: 13C ratio measures ~C"
                        f"{c_est:.0f}, formula claims C{n_c}",
                    )
                    out["c13_clamp"] += 1
                except L.LedgerError:
                    pass
        elif (
            expected >= 1.5 * cfg.height_cutoff
            and _peak_near(mzs, r["mz"] + _D13C) is None
        ):
            # twin-satellite fallback: when the peak has a halogen isotope
            # twin, the twin's OWN 13C satellite (13C+81Br / 13C+37Cl) is
            # equally valid carbon evidence. v20 falsely cleared C3H6O3.Br-
            # (10.3k cps): its plain 13C is peak-picker-lost, but the twin's
            # satellite at +1.998+1.0034 exists and is carbon-consistent.
            twins = ledger[
                (ledger["parent_peak_id"] == r["peak_id"])
                & ledger["iso_label"].astype(str).str.contains("Br|Cl", regex=True)
            ]
            if any(
                _peak_near(mzs, float(t["mz"]) + _D13C) is not None
                for _, t in twins.iterrows()
            ):
                continue
            # cross-channel fallback: the SAME neutral independently assigned
            # High/Good on another peak (other adduct) is positive evidence
            # that outweighs one absent satellite -- an absent 13C can be a
            # peak-picker loss, an agreeing second channel cannot. (v21
            # cleared five sub-ppm [M+Br]- partners of Good [M-H]-
            # assignments, e.g. C10H16O6 at 311.013 / 2.4k cps.)
            others = ledger[
                (ledger["role"] == L.ROLE_M0)
                & (ledger["peak_id"] != r["peak_id"])
                & (ledger["neutral_formula"] == r["neutral_formula"])
                & ledger["confidence"].astype(str).str.startswith(("High", "Good"))
            ]
            if len(others):
                continue
            try:
                L.clear_assignment(
                    ledger,
                    r["peak_id"],
                    reason=f"isotope audit: predicted 13C satellite "
                    f"({expected:.0f} cps) absent from spectrum",
                )
                out["c13_missing"] += 1
            except L.LedgerError:
                continue

    n = sum(out.values())
    if n:
        log(
            f"[audit] isotope physics: {out['doublet_child']} doublet twins "
            f"re-attached, {out['doublet_cleared']} no-Br doublet formulas "
            f"cleared, {out['c13_attached']} 13C satellites attached, "
            f"{out['c13_clamp']} carbon-clamp clears, "
            f"{out['c13_missing']} missing-13C clears"
        )
    return out


def audit_mass_gate(ledger: pd.DataFrame, cfg: PassConfig, *, log=print) -> dict:
    """Post-run sweep: apply the calibrated mass gate to commits that predate
    calibration (pass 1 runs before the backbone exists). Clears, never
    rewrites: a >4-sigma mass error means the formula is wrong no matter what
    the isotope pattern says; a 2..4-sigma Low/Suspect with no pattern
    evidence is just the closest of many candidates."""
    out = {"cleared_z": 0, "cleared_z_noiso": 0, "cleared_nan": 0}
    if cfg.cal_mu is None:
        return out
    m0 = ledger[(ledger["role"] == L.ROLE_M0) & ~ledger["locked"].astype(bool)]
    parents_with_kids = set(
        ledger.loc[ledger["role"] == L.ROLE_ISO, "parent_peak_id"].dropna()
    )
    for _, r in m0.iterrows():
        weak = not str(r["confidence"]).startswith(("High", "Good"))
        has_kids = r["peak_id"] in parents_with_kids
        z = z_of(r["ppm_error"], cfg)
        try:
            if z is None:
                if pd.isna(r["ppm_error"]) and weak and not has_kids:
                    L.clear_assignment(
                        ledger, r["peak_id"], reason="mass-gate: no ppm error"
                    )
                    out["cleared_nan"] += 1
            elif z > cfg.cal_z_pattern:
                L.clear_assignment(
                    ledger,
                    r["peak_id"],
                    reason=f"mass-gate: z={z:.1f} > {cfg.cal_z_pattern}",
                )
                out["cleared_z"] += 1
            elif z > cfg.cal_z_accept and weak and not has_kids:
                L.clear_assignment(
                    ledger,
                    r["peak_id"],
                    reason=f"mass-gate: z={z:.1f} without pattern evidence",
                )
                out["cleared_z_noiso"] += 1
        except L.LedgerError:
            continue
    n = sum(out.values())
    if n:
        log(
            f"[audit] mass gate cleared {n} assignments "
            f"(z>{cfg.cal_z_pattern}: {out['cleared_z']}, "
            f"{cfg.cal_z_accept}<z<={cfg.cal_z_pattern} no-evidence: "
            f"{out['cleared_z_noiso']}, no-ppm: {out['cleared_nan']})"
        )
    return out

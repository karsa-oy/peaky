"""The three-pass assignment director.

Consumes the ledger and the Mascope oracle; produces a fully annotated ledger.

  Pass 1  Lock the high-confidence CHO / CHON backbone. Enumerate candidates
          (cheminfo + grid fallback), score with match_compounds, arbitrate per
          peak (complexity-penalised), commit M0 owners, attach the
          isotopologue children Mascope attributed, and LOCK the peaks that
          clear the High bar.
  Pass 2  GKA series expansion from the locked anchors (CHO/CHON + siloxane +
          CF2). Every propagated formula is validated by match_compounds.
  Pass 3  Low-quality recovery: open the context's contaminant families
          (sulfate/organosulfate, nitrate, siloxane, ...), score, commit at a
          lower floor with explicit Low/Suspect commentary.

Arbitration (`arbitrate`) is a PURE function over the flat scored table, so it
is unit-tested offline. The pass drivers wrap it with the live oracle.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from . import chemistry as C
from . import contexts as X
from . import io_mascope as IO
from . import ledger as L
from . import series_gka as G

__version__ = "0.2.0"


@dataclass
class PassConfig:
    ppm: float = 1.0                 # user m/z trust
    search_ppm: float = 5.0          # tolerance used for enumeration/scoring
    height_cutoff: float = 500.0
    limit_per_peak: int = 25
    workers: int = 12
    # confidence thresholds (on the RAW min(ion,compound) score)
    tau_high: float = 0.90
    tau_good: float = 0.80
    tau_low: float = 0.70
    tau_suspect: float = 0.50
    complexity_cap: float = 0.20
    require_iso_for_high: bool = True
    series_ppm: float = 3.0
    series_min_score: float = 0.60
    series_max_iter: int = 3   # iterative GKA: chain confirmed members as anchors
    # Pass 4 (residual explainer) acceptance policy: <=strict ppm on score
    # alone; up to pattern ppm ONLY with pattern evidence (confirmed isotope
    # partner / >=2 series anchors). DBE-only plausibility in pass 4.
    residual_ppm_strict: float = 1.0
    residual_ppm_pattern: float = 4.0
    residual_max_steps: int = 2
    # explicit ionization-mechanism ids for match_compounds. None = server
    # auto-selects the sample's configured channels; set by assign.run to the
    # sample's channels PLUS extras like +CO3- so background air-ion adducts
    # get scored too.
    mechanism_ids: list | None = None
    # enumeration: the local grid is the primary, reliable candidate source.
    # cheminfo is an optional best-effort enrichment (compound names) and is the
    # flaky/slow dependency, so it is OFF by default in the search path.
    use_cheminfo: bool = False
    # isotopologue gating: a heteroatom in the NEUTRAL must be backed by its
    # diagnostic isotope confirmed by Mascope, else the candidate is penalised.
    # Cl/Br satellites are large (always visible if real) -> strong penalty;
    # 34S is small (4.4%) -> softer penalty.
    het_iso_penalty_halogen: float = 0.30
    het_iso_penalty_S: float = 0.12
    # The reagent halogen (e.g. Br in Br-CIMS) is special: its heavy isotope in
    # the ION cannot prove the halogen sits in the NEUTRAL (covalent X(Br)[M-H]-
    # and Y.HBr.Br- / Y[M+Br]- aliases share the ion). Confirmation therefore
    # waives only the gate penalty, never the complexity prior, so the
    # adduct/cluster interpretation wins ties. Set by assign.run.
    reagent_element: str | None = None
    # Self-calibration mass gate (ROADMAP 1): mu/sigma of the ppm error fitted
    # on the pass-1 High/Good CHO-CHON backbone (set by assign.run via
    # calibrate()). A candidate is judged by z = |ppm - mu| / sigma:
    # z <= cal_z_accept on score alone; up to cal_z_pattern only WITH pattern
    # evidence (confirmed isotopologue or series membership); beyond that the
    # best fit within tolerance is just the closest of many -- reject. A match
    # with NO ppm at all carries no mass evidence and is never committed.
    cal_mu: float | None = None
    cal_sigma: float | None = None
    cal_z_accept: float = 2.0
    cal_z_pattern: float = 4.0
    cal_sigma_floor: float = 0.25   # don't let a lucky tight fit reject everything
    cal_min_n: int = 20             # min backbone size to trust a fit
    # Channel priors: the reagent / deprotonation channels are PRIMARY; the
    # background air-ion channels (carbonate, superoxide, electron attachment)
    # are MINOR -- real but rare, and offering them to every peak doubles the
    # alias space. A minor-channel candidate pays a ranking penalty (so a
    # near-tie goes to the primary channel) and a minor-channel WINNER may only
    # commit with corroboration: a Good+ score, series-evidence method, or the
    # same neutral independently assigned via a primary channel.
    minor_channels: tuple = ("[M+CO3]-", "[M+O2]-", "[M]-.")
    minor_channel_penalty: float = 0.12


# ---------------------------------------------------------------------------
# Confidence
# ---------------------------------------------------------------------------
def confidence_label(score: float, ppm: float | None, n_iso: int, tied: bool,
                     cfg: PassConfig, suffix: str = "") -> str:
    a = abs(ppm) if ppm is not None and pd.notna(ppm) else 99.0
    if score >= cfg.tau_high and a <= cfg.ppm * 1.5 and n_iso >= 1 and not tied:
        lab = "High"
    elif score >= cfg.tau_good and a <= cfg.ppm * 2:
        lab = "Good"
    elif score >= cfg.tau_low:
        lab = "Low"
    elif score >= cfg.tau_suspect:
        lab = "Suspect"
    else:
        lab = "Reject"
    return f"{lab} ({suffix})" if suffix and lab not in ("Reject",) else lab


# ---------------------------------------------------------------------------
# Self-calibration mass gate  (ROADMAP 1)
# ---------------------------------------------------------------------------
_BACKBONE_ELEMENTS = {"C", "H", "O", "N"}


def calibrate(ledger: pd.DataFrame, cfg: PassConfig, *, log=print) -> tuple | None:
    """Fit the instrument's real mass accuracy (mu, sigma of the ppm error)
    over the committed High/Good CHO-CHON backbone and store it on cfg.

    Robust fit: median + scaled MAD, so a few bad commits can't widen the
    window they would then pass through. Returns (mu, sigma, n) or None when
    the backbone is too small to trust."""
    m0 = ledger[(ledger["role"] == L.ROLE_M0) & ledger["ppm_error"].notna()]
    rows = m0[m0["confidence"].astype(str).str.startswith(("High", "Good"))]
    keep = rows["neutral_formula"].astype(str).map(
        lambda f: set(C.parse_formula(f)) <= _BACKBONE_ELEMENTS)
    ppm = rows.loc[keep, "ppm_error"].astype(float)
    if len(ppm) < cfg.cal_min_n:
        log(f"[calibrate] backbone too small (n={len(ppm)} < {cfg.cal_min_n}); "
            "mass gate stays off")
        return None
    mu = float(ppm.median())
    sigma = max(float(1.4826 * (ppm - mu).abs().median()), cfg.cal_sigma_floor)
    cfg.cal_mu, cfg.cal_sigma = mu, sigma
    log(f"[calibrate] backbone n={len(ppm)}: ppm mu={mu:+.3f} sigma={sigma:.3f} "
        f"-> accept |z|<={cfg.cal_z_accept} ({mu - cfg.cal_z_accept * sigma:+.2f}"
        f"..{mu + cfg.cal_z_accept * sigma:+.2f} ppm), pattern-evidence up to "
        f"|z|<={cfg.cal_z_pattern}")
    return mu, sigma, len(ppm)


def z_of(ppm, cfg: PassConfig) -> float | None:
    """Calibrated z-score of a ppm error; None when uncalibrated or ppm is NaN."""
    if cfg.cal_mu is None or cfg.cal_sigma is None:
        return None
    if ppm is None or pd.isna(ppm):
        return None
    return abs(float(ppm) - cfg.cal_mu) / cfg.cal_sigma


# ---------------------------------------------------------------------------
# Arbitration  (pure)
# ---------------------------------------------------------------------------
def arbitrate(scored: pd.DataFrame, cfg: PassConfig) -> dict:
    """Decide a single best M0 owner per peak from the flat per-isotopologue
    scored table, with complexity-penalised effective scores.

    Returns:
      {
        'winners':  DataFrame[peak_id, neutral, adduct, ion_formula, ion_score,
                              compound_score, raw_score, eff_score, ppm_error,
                              n_iso, tied, alternatives(list)],
        'iso_children': DataFrame[peak_id, parent_peak_id, iso_label, iso_score],
      }
    """
    if scored is None or len(scored) == 0:
        return {"winners": pd.DataFrame(), "iso_children": pd.DataFrame()}

    # base ions that matched a real peak with a usable ion score
    base = scored[scored["is_base"]
                  & scored["sample_peak_id"].notna()
                  & scored["ion_score"].notna()].copy()
    if len(base) == 0:
        return {"winners": pd.DataFrame(), "iso_children": pd.DataFrame()}

    base["raw_score"] = base[["ion_score", "compound_score"]].min(axis=1, skipna=True)

    # Mascope-confirmed isotopologues per (compound, ion): non-base rows for the
    # same ion that matched a real peak with score > 0.4
    iso = scored[(~scored["is_base"])
                 & scored["sample_peak_id"].notna()
                 & (pd.to_numeric(scored["iso_score"], errors="coerce").fillna(0) > 0.4)]
    iso_count = (iso.groupby(["compound_formula", "ion_formula"])["sample_peak_id"]
                 .nunique().to_dict())

    # diagnostic isotope labels confirmed per (compound, ion), e.g. {'13C','81Br','34S'}
    iso_labels: dict[tuple, set] = {}
    for (cf, ifl), g in iso.groupby(["compound_formula", "ion_formula"]):
        labs: set[str] = set()
        for lab in g["iso_label"].astype(str):
            labs.update(lab.split("+"))
        iso_labels[(cf, ifl)] = labs

    # Evidence-adjusted penalty on the NEUTRAL. The complexity prior (skepticism
    # about heteroatoms) applies, BUT a heteroatom whose diagnostic isotope is
    # Mascope-confirmed has its skepticism WAIVED -- direct evidence overrides
    # the prior. Heteroatoms with a usable isotope diagnostic and no confirmation
    # additionally take the gate penalty. Elements without a light-isotope
    # diagnostic (N, P, Si, I, F) keep the plain complexity prior.
    _DIAG = {"Cl": ("37Cl", cfg.het_iso_penalty_halogen),
             "Br": ("81Br", cfg.het_iso_penalty_halogen),
             "S":  ("34S", cfg.het_iso_penalty_S)}

    def _evidence_penalty(row) -> float:
        cnt = C.parse_formula(row["compound_formula"])
        labs = iso_labels.get((row["compound_formula"], row["ion_formula"]), set())
        pen = 0.0
        for el, n in cnt.items():
            if n <= 0 or el not in C._COMPLEXITY_WEIGHT:
                continue
            prior = min(C._COMPLEXITY_WEIGHT[el] * n * 0.01, cfg.complexity_cap)
            if el in _DIAG:
                diag, gate = _DIAG[el]
                confirmed = any(s.startswith(diag) for s in labs)
                if el == cfg.reagent_element:
                    # ion isotope can't prove NEUTRAL ownership of the reagent
                    # halogen: keep the prior; gate only if not even ion-level
                    # confirmation exists.
                    pen += prior + (0.0 if confirmed else gate)
                elif confirmed:
                    continue            # confirmed -> waive skepticism
                else:
                    pen += prior + gate
            else:
                pen += prior            # no diagnostic -> plain prior
        return min(pen, 0.50)

    base["penalty"] = base.apply(_evidence_penalty, axis=1)
    # Channel prior: minor background channels (CO3-/O2-/electron attachment)
    # pay a ranking penalty so a near-tie goes to the primary channel. The
    # commit-side corroboration gate lives in commit_winners.
    base["adduct_label"] = base.apply(_mech_to_adduct, axis=1)
    base["eff_score"] = (base["raw_score"] - base["penalty"]
                         - base["adduct_label"].isin(cfg.minor_channels)
                         * cfg.minor_channel_penalty)

    winners = []
    for pid, grp in base.groupby("sample_peak_id"):
        grp = grp.sort_values("eff_score", ascending=False)
        top = grp.iloc[0]
        n_iso = int(iso_count.get((top["compound_formula"], top["ion_formula"]), 0))
        runner_eff = float(grp.iloc[1]["eff_score"]) if len(grp) > 1 else None
        tied = runner_eff is not None and (float(top["eff_score"]) - runner_eff) < 0.05
        alts = []
        for _, r in grp.iloc[1:4].iterrows():
            alts.append({"formula": r["compound_formula"], "adduct": r["adduct_label"],
                         "ion_score": _f(r["ion_score"]), "raw_score": _f(r["raw_score"]),
                         "eff_score": _f(r["eff_score"]), "ppm": _f(r["ppm_error"])})
        winners.append({
            "peak_id": pid,
            "neutral": top["compound_formula"],
            "ion_formula": top["ion_formula"],
            "adduct": top["adduct_label"],
            "ion_score": _f(top["ion_score"]),
            "compound_score": _f(top["compound_score"]),
            "raw_score": _f(top["raw_score"]),
            "eff_score": _f(top["eff_score"]),
            "ppm_error": _f(top["ppm_error"]),
            "n_iso": n_iso,
            "tied": bool(tied),
            "alternatives": alts,
        })
    win_df = pd.DataFrame(winners)

    # iso children attributed to the winning (compound, ion) pairs only
    children = []
    win_keys = {(w["neutral"], w["ion_formula"]): w["peak_id"] for w in winners}
    for _, r in iso.iterrows():
        key = (r["compound_formula"], r["ion_formula"])
        if key in win_keys:
            parent = win_keys[key]
            if r["sample_peak_id"] != parent:
                children.append({"peak_id": r["sample_peak_id"],
                                 "parent_peak_id": parent,
                                 "iso_label": r["iso_label"],
                                 "iso_score": _f(r["iso_score"])})
    return {"winners": win_df, "iso_children": pd.DataFrame(children)}


def _f(v):
    return None if v is None or pd.isna(v) else float(v)


# Exact (ion - compound) element difference -> adduct label. Keys are the
# alphabetically sorted nonzero element deltas. A diff with no entry here is a
# NEW mechanism leaking through unlabeled — fall back to [M-H]- but only after
# every registered channel is covered (the old heuristic had no CO3 branch, so
# all 168 [M+CO3]- matches in v13 were silently labeled [M-H]-).
_DIFF_TO_ADDUCT = {
    (("H", -1),): "[M-H]-",
    (("Br", 1),): "[M+Br]-",
    (("Cl", 1),): "[M+Cl]-",
    (("I", 1),): "[M+I]-",
    (("N", 1), ("O", 3)): "[M+NO3]-",
    (("H", 1), ("O", 4), ("S", 1)): "[M+HSO4]-",
    (("C", 1), ("H", 1), ("O", 2)): "[M+CHO2]-",
    (("C", 2), ("H", 3), ("O", 2)): "[M+C2H3O2]-",
    (("C", 1), ("O", 3)): "[M+CO3]-",
    (("O", 2),): "[M+O2]-",
    (): "[M]-.",
    (("H", 1),): "[M+H]+",
    (("Na", 1),): "[M+Na]+",
    (("H", 4), ("N", 1)): "[M+NH4]+",
    (("K", 1),): "[M+K]+",
}


def _mech_to_adduct(row) -> str:
    """Adduct label from the exact ion-vs-compound element difference."""
    ci = C.parse_formula(str(row.get("ion_formula") or ""))
    cc = C.parse_formula(str(row.get("compound_formula") or ""))
    diff = tuple(sorted(
        (el, ci.get(el, 0) - cc.get(el, 0))
        for el in set(ci) | set(cc) if ci.get(el, 0) != cc.get(el, 0)))
    return _DIFF_TO_ADDUCT.get(diff, "[M-H]-")


def commit_winners(ledger: pd.DataFrame, arb: dict, *, pass_no: int, method: str,
                   context: str, cfg: PassConfig, lock: bool,
                   min_raw_score: float, confidence_suffix: str = "",
                   claim_unexplained_only: bool = False,
                   only_peaks: set | None = None) -> dict:
    """Commit arbitration winners + their isotope children into the ledger.
    Skips peaks that are locked or already own a better assignment. Returns a
    summary dict."""
    win = arb.get("winners", pd.DataFrame())
    kids = arb.get("iso_children", pd.DataFrame())
    if len(win):
        # stronger claims commit first, so within-pass conflicts (two winners
        # wanting the same isotope child) resolve toward the better evidence
        win = win.sort_values("eff_score", ascending=False, na_position="last")
    committed, locked_ids = [], []
    rejected = {"nan_ppm": 0, "mass_gate": 0, "minor_channel": 0}
    # series/GKA proposals and known-neutral completions carry structural
    # evidence by construction (chain membership / an already-assigned neutral)
    series_like = ("series" in method) or ("gka" in method) \
        or ("completion" in method)
    for _, w in win.iterrows():
        pid = w["peak_id"]
        if only_peaks is not None and pid not in only_peaks:
            continue   # evidence-scoped pass: only its own member peaks
        if w["raw_score"] is None or w["raw_score"] < min_raw_score:
            continue
        # A match with no mass error at all carries no mass evidence -- the 9
        # score-only bromides of v13 came through here. Never commit.
        if w["ppm_error"] is None or pd.isna(w["ppm_error"]):
            rejected["nan_ppm"] += 1
            continue
        # Calibrated mass gate: judge the ppm error against the instrument's
        # real accuracy, not a fixed window. Pattern evidence (a Mascope-
        # confirmed isotopologue or a series-derived proposal) buys the
        # 2..4 sigma band; nothing buys more.
        z = z_of(w["ppm_error"], cfg)
        if z is not None:
            pattern = (w["n_iso"] or 0) >= 1 or series_like
            if z > cfg.cal_z_pattern or (z > cfg.cal_z_accept and not pattern):
                rejected["mass_gate"] += 1
                continue
        # Minor-channel corroboration gate: a CO3-/O2-/M-. winner must be Good+
        # on raw score, or come from series evidence, or have its neutral
        # independently assigned via a primary channel elsewhere in the ledger.
        if w["adduct"] in cfg.minor_channels and w["raw_score"] < cfg.tau_good \
                and not series_like:
            others = ledger[(ledger["role"] == L.ROLE_M0)
                            & (ledger["neutral_formula"] == w["neutral"])
                            & ~ledger["adduct"].isin(cfg.minor_channels)]
            if len(others) == 0:
                rejected["minor_channel"] += 1
                continue
        # Reagent-halogen decomposition policy: covalent X(Br)[M-H]- and
        # Y.Br- / Y.HBr.Br- produce the IDENTICAL ion, so no spectral evidence
        # can distinguish them. The reagent halogen is read as part of the
        # adduct/cluster whenever the decomposition Y = X - HX is structurally
        # valid (user policy, 2026-06-11).
        w = _prefer_adduct_reading(w, cfg)
        try:
            if L.is_locked(ledger, pid):
                continue
            cur_role = L.role_of(ledger, pid)
        except L.LedgerError:
            continue  # peak not in ledger (sub-threshold) -> skip
        # later passes claim free peaks only -- a contaminant family must not
        # displace an earlier pass's isotope-confirmed assignment
        if claim_unexplained_only and cur_role != L.ROLE_UNEXPLAINED:
            continue
        # don't overwrite an existing higher-or-equal raw score
        if cur_role == L.ROLE_M0:
            existing = ledger.loc[ledger.peak_id == pid, "ion_score"].iloc[0]
            if pd.notna(existing) and existing >= w["raw_score"]:
                continue
        conf = confidence_label(w["raw_score"], w["ppm_error"], w["n_iso"],
                                w["tied"], cfg, suffix=confidence_suffix)
        if conf == "Reject":
            continue
        commentary = _commentary(w, pass_no, method)
        try:
            L.commit_assignment(
                ledger, pid, neutral_formula=w["neutral"], adduct=w["adduct"],
                ion_formula=w["ion_formula"], ion_score=w["ion_score"],
                compound_score=w["compound_score"], ppm_error=w["ppm_error"],
                pass_no=pass_no, method=method, confidence=conf,
                commentary=commentary, alternatives=w["alternatives"],
                isotopologues=_iso_list(kids, pid),
                overwrite=(cur_role == L.ROLE_M0))
            committed.append(pid)
            if lock and conf == "High":
                locked_ids.append(pid)
        except L.LedgerError:
            continue
    # attach isotope children (parents must now own M0)
    n_iso_attached = 0
    n_displaced = 0
    if len(kids):
        for _, k in kids.iterrows():
            try:
                if not (k["parent_peak_id"] in set(committed) or
                        L.role_of(ledger, k["parent_peak_id"]) == L.ROLE_M0):
                    continue
                # M0-vs-iso-child arbitration: when the predicted child peak
                # already owns its own M0, a Mascope-confirmed two-peak
                # explanation (parent + correct-ratio satellite) beats a
                # weaker independent formula on the child -- displace it.
                # (Cases: 428.99 'Br2 amine' was the 81Br twin of 426.99;
                # 297.02 'organosilicon CO3' was the 81Br twin of 295.02.)
                if L.role_of(ledger, k["peak_id"]) == L.ROLE_M0 \
                        and not L.is_locked(ledger, k["peak_id"]):
                    crow = ledger.loc[ledger.peak_id == k["peak_id"]].iloc[0]
                    prow = ledger.loc[ledger.peak_id == k["parent_peak_id"]].iloc[0]
                    child_conf = str(crow["confidence"])
                    # only High children are immune: a Mascope-confirmed
                    # two-peak explanation (parent + correct-ratio satellite)
                    # beats an equal-or-weaker single-peak formula (v16 audit:
                    # the 462.99/464.99 'Good vs Good' doublet never displaced)
                    child_immune = child_conf.startswith("High")
                    parent_stronger = (pd.isna(crow["ion_score"])
                                       or (pd.notna(prow["ion_score"]) and
                                           prow["ion_score"] >= crow["ion_score"]))
                    if not child_immune and parent_stronger:
                        L.displace_to_isotopologue(
                            ledger, k["peak_id"], k["parent_peak_id"],
                            iso_label=k["iso_label"], iso_match_score=k["iso_score"])
                        n_displaced += 1
                        n_iso_attached += 1
                    continue
                L.attach_isotopologue(ledger, k["peak_id"], k["parent_peak_id"],
                                      iso_label=k["iso_label"],
                                      iso_match_score=k["iso_score"])
                n_iso_attached += 1
            except L.LedgerError:
                continue
    if locked_ids:
        L.lock_peaks(ledger, locked_ids)
    return {"committed": len(committed), "locked": len(locked_ids),
            "iso_attached": n_iso_attached, "iso_displaced": n_displaced,
            "rejected": rejected}


# Known instrument-contaminant series, labeled BEFORE pass 1 (reagent-style):
# these compositions are H-rich (H/C up to 4) and Si-bearing, so the CHO/CHON
# pass-1 grid mis-claims their peaks with O-rich fantasies and LOCKS them
# (v24: 244.9668 'C5H10O6 High' and 318.9856 'C7H16O7Si' were really the n=2/3
# dimethylsilanediol oligomers -- both sat on the z>2 watch-list).
def _silanediol_series(n_max: int = 8) -> list[str]:
    """HO-(Si(CH3)2-O)n-H: PDMS hydrolysis products, the classic inlet/tubing
    contamination. Composition C(2n)H(6n+2)O(n+1)Si(n)."""
    return [f"C{2*n}H{6*n+2}O{n+1}Si{n}" for n in range(1, n_max + 1)]


def run_pass0_contaminants(client, sample_id: str, ledger: pd.DataFrame,
                           profile, cfg: PassConfig, adducts: list[str], *,
                           score_fn=None, log=print) -> dict:
    """Pass 0 -- label known instrument-contaminant series before the organic
    passes run. Mascope scores the explicit compositions (its isotope model
    covers 29Si/30Si + the reagent halogen); commits are LOCKED so pass 1
    cannot displace them with grid CHO fits."""
    score_fn = score_fn or IO.score_candidates
    out = {"committed": 0, "locked": 0, "iso_attached": 0}
    formulas = _silanediol_series()
    scored = score_fn(client, sample_id, formulas,
                      mechanism_ids=cfg.mechanism_ids)
    if scored is None or len(scored) == 0:
        log("[pass0] WARNING scoring returned EMPTY for the contaminant list")
        out["scoring_empty"] = True
        return out
    base = scored[scored["is_base"] & scored["sample_peak_id"].notna()
                  & scored["ion_score"].notna()]
    kids = scored[(~scored["is_base"]) & scored["sample_peak_id"].notna()
                  & (pd.to_numeric(scored["iso_score"], errors="coerce")
                     .fillna(0) > 0.4)]
    mzs = ledger["mz"]
    for _, r in base.iterrows():
        ppm = r["ppm_error"]
        if ppm is None or pd.isna(ppm) or abs(float(ppm)) > 2.0:
            continue   # the known-list privilege still requires the MASS
        pid = r["sample_peak_id"]
        try:
            if L.role_of(ledger, pid) != L.ROLE_UNEXPLAINED:
                continue
            # self-twin consistency: a [M+Br]- contaminant claim must own a
            # consistent 81Br twin of its OWN. v25 lesson: silanediol n=1
            # (170.9482) collided with lactic acid's 81Br child (170.9485);
            # the 12k cps peak's twin at 172.946 was 427 cps (ratio 0.04) --
            # the peak belongs to the lactic-acid envelope, not the
            # contaminant. A composite minor component cannot be LOCKED.
            if "Br" in str(r["ion_formula"]):
                i0 = ledger.index[ledger["peak_id"] == pid][0]
                m0 = float(ledger.at[i0, "mz"])
                h0 = float(ledger.at[i0, "height"])
                tw = _peak_near(mzs, m0 + _DBR, ppm=8.0)
                rt = (float(ledger.at[tw, "height"]) / h0) if tw is not None else 0.0
                if not (0.5 <= rt <= 1.7):
                    log(f"[pass0] skip {r['compound_formula']} @{m0:.4f}: "
                        f"own-81Br-twin ratio {rt:.2f} inconsistent "
                        f"(composite or wrong claim)")
                    continue
            fam_kids = kids[kids["compound_formula"] == r["compound_formula"]]
            n_kids = int((fam_kids["sample_peak_id"] != pid).sum())
            conf = ("Good (contaminant)" if float(r["ion_score"]) >= 0.7
                    or n_kids >= 2 else "Low (contaminant)")
            L.commit_assignment(
                ledger, pid, neutral_formula=r["compound_formula"],
                adduct=_mech_to_adduct(r), ion_formula=r["ion_formula"],
                ion_score=float(r["ion_score"]),
                compound_score=_f(r.get("compound_score")),
                ppm_error=float(ppm), pass_no=0,
                method="contaminant:silanediol", confidence=conf,
                commentary=(f"Pass 0 (known contaminant): dimethylsilanediol "
                            f"oligomer {r['compound_formula']} "
                            f"{_mech_to_adduct(r)}, ppm {float(ppm):.2f}; "
                            f"PDMS hydrolysis series (inlet/tubing), "
                            f"GKA-discovered 2026-06-12"))
            out["committed"] += 1
            L.lock_peaks(ledger, [pid])
            out["locked"] += 1
            for _, k in fam_kids.iterrows():
                if k["sample_peak_id"] == pid:
                    continue
                try:
                    L.attach_isotopologue(ledger, k["sample_peak_id"], pid,
                                          iso_label=k["iso_label"],
                                          iso_match_score=_f(k["iso_score"]))
                    out["iso_attached"] += 1
                except L.LedgerError:
                    continue
        except L.LedgerError:
            continue
    log(f"[pass0] {out}")
    return out


def run_pass5_completion(client, sample_id: str, ledger: pd.DataFrame, profile,
                         cfg: PassConfig, adducts: list[str], *,
                         score_fn=None, log=print) -> dict:
    """Pass 5 -- known-neutral completion. Opens NO new formula space: it only
    proposes neutrals the run already believes, onto unexplained peaks that are

      (a) cross-channel partners: another registered adduct of an assigned
          neutral (e.g. the [M+Br]- of a Good [M-H]- compound, TFA's [M-H]-), or
      (b) series-gap members: a CH2-bracketed gap inside an assigned homolog
          ladder (the C2/C5/C6 hydroxy-acid ladder whose missing C3 rung,
          C3H6O3.Br- at 10.3k cps, was the biggest unexplained peak of v20).

    Mascope still scores every proposal (its ppm/isotope attribution is
    authoritative) and the normal commit gates apply; 'completion' in the
    method name grants the pattern-evidence band, since the evidence is an
    independently assigned neutral."""
    score_fn = score_fn or IO.score_candidates
    out = {"committed": 0, "locked": 0, "iso_attached": 0}
    m0 = ledger[ledger["role"] == L.ROLE_M0].dropna(subset=["neutral_formula"])
    anchors = m0[m0["confidence"].astype(str).str.startswith(("High", "Good"))]
    assigned = set(anchors["neutral_formula"])
    if not assigned:
        log("[pass5] no High/Good anchors; skipping")
        return out
    un = ledger[ledger["role"] == L.ROLE_UNEXPLAINED]
    gadducts = [a for a in adducts if a in C.ADDUCT_SHIFTS]

    def un_peak_near(target: float):
        if not len(un):
            return None
        d = (un["mz"] - target).abs()
        i = d.idxmin()
        return un.at[i, "peak_id"] if d.loc[i] <= target * cfg.search_ppm * 1e-6 else None

    targets: dict[str, set] = {}
    n_cross = n_gap = 0
    # (a) cross-channel partners of assigned neutrals
    for nf in assigned:
        for ad in gadducts:
            try:
                pid = un_peak_near(C.ion_mz(nf, ad))
            except Exception:
                continue
            if pid is not None:
                targets.setdefault(nf, set()).add(pid)
                n_cross += 1
    # (b) CH2-bracketed gaps between assigned ladder anchors
    for nf in sorted(assigned):
        for k in (2, 3):
            if G.formula_add(nf, "CH2", k) not in assigned:
                continue
            for j in range(1, k):
                mid = G.formula_add(nf, "CH2", j)
                if not mid or mid in assigned or not C.dbe_ok(mid)[0]:
                    continue
                for ad in gadducts:
                    pid = un_peak_near(C.ion_mz(mid, ad))
                    if pid is not None:
                        targets.setdefault(mid, set()).add(pid)
                        n_gap += 1
    if not targets:
        log("[pass5] no completion targets")
        return out
    only = set().union(*targets.values())
    log(f"[pass5] {len(targets)} known neutrals -> {len(only)} target peaks "
        f"({n_cross} cross-channel, {n_gap} series-gap)")
    scored = score_fn(client, sample_id, sorted(targets),
                      mechanism_ids=cfg.mechanism_ids)
    if scored is None or len(scored) == 0:
        log(f"[pass5] WARNING scoring returned EMPTY for {len(targets)} known "
            f"neutrals -- server likely degraded; completion skipped")
        out["scoring_empty"] = True
        return out
    arb = arbitrate(scored, cfg)
    s = commit_winners(ledger, arb, pass_no=5, method="completion:known-neutral",
                       context=profile.label, cfg=cfg, lock=False,
                       min_raw_score=cfg.tau_suspect,
                       confidence_suffix="completion",
                       claim_unexplained_only=True, only_peaks=only)
    log(f"[pass5] {s}")
    return s


# isotope constants for the post-run physics audit
_D13C = 1.0033548          # 13C - 12C
_DBR = 1.9979535           # 81Br - 79Br
_R13C = 0.0107             # 13C natural abundance per carbon
_R81BR = 0.9728            # 81Br/79Br abundance ratio


def _peak_near(mzs: "pd.Series", target: float, ppm: float = 5.0):
    """Index of the closest ledger peak within ppm of target, else None."""
    tol = target * ppm * 1e-6
    d = (mzs - target).abs()
    i = d.idxmin()
    return i if d.loc[i] <= tol else None


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
    out = {"doublet_child": 0, "doublet_cleared": 0, "c13_attached": 0,
           "c13_clamp": 0, "c13_missing": 0}
    mzs = ledger["mz"]

    # --- 1. Br-doublet repair over committed M0s ---
    m0 = ledger[(ledger["role"] == L.ROLE_M0)
                & ~ledger["locked"].astype(bool)].sort_values("mz")
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
                if L.role_of(ledger, lt.peak_id) != L.ROLE_M0 \
                        or L.role_of(ledger, hv.peak_id) != L.ROLE_M0:
                    continue
                n_br = C.parse_formula(str(lt.ion_formula)).get("Br", 0)
                if n_br >= 1:
                    L.clear_assignment(ledger, hv.peak_id,
                                       reason=f"isotope audit: 81Br twin of "
                                              f"{lt.mz:.4f} (ratio {hr:.2f})")
                    L.attach_isotopologue(ledger, hv.peak_id, lt.peak_id,
                                          iso_label="81Br")
                    out["doublet_child"] += 1
                else:
                    L.clear_assignment(
                        ledger, lt.peak_id,
                        reason=f"isotope audit: Br doublet with {hv.mz:.4f} "
                               f"(ratio {hr:.2f}) but no Br in formula")
                    L.clear_assignment(
                        ledger, hv.peak_id,
                        reason=f"isotope audit: Br doublet with {lt.mz:.4f} "
                               f"(ratio {hr:.2f}) but no Br in formula")
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
        k = kids[(kids["parent_peak_id"] == r["peak_id"])
                 & (kids["iso_label"].astype(str) == "13C")]
        if not len(k):
            j = _peak_near(mzs, r["mz"] + _D13C)
            if j is not None and ledger.at[j, "role"] == L.ROLE_UNEXPLAINED \
                    and expected > 0 \
                    and 0.3 <= ledger.at[j, "height"] / expected <= 2.5:
                try:
                    L.attach_isotopologue(ledger, ledger.at[j, "peak_id"],
                                          r["peak_id"], iso_label="13C")
                    out["c13_attached"] += 1
                    k = ledger.loc[[j]]
                except L.LedgerError:
                    pass
        if len(k):
            c_est = (float(k.iloc[0]["height"]) / float(r["height"])) / _R13C
            if n_c >= 8 and abs(c_est - n_c) > max(2.5, 0.35 * n_c):
                try:
                    L.clear_assignment(
                        ledger, r["peak_id"],
                        reason=f"isotope audit: 13C ratio measures ~C"
                               f"{c_est:.0f}, formula claims C{n_c}")
                    out["c13_clamp"] += 1
                except L.LedgerError:
                    pass
        elif expected >= 1.5 * cfg.height_cutoff \
                and _peak_near(mzs, r["mz"] + _D13C) is None:
            # twin-satellite fallback: when the peak has a halogen isotope
            # twin, the twin's OWN 13C satellite (13C+81Br / 13C+37Cl) is
            # equally valid carbon evidence. v20 falsely cleared C3H6O3.Br-
            # (10.3k cps): its plain 13C is peak-picker-lost, but the twin's
            # satellite at +1.998+1.0034 exists and is carbon-consistent.
            twins = ledger[(ledger["parent_peak_id"] == r["peak_id"])
                           & ledger["iso_label"].astype(str)
                           .str.contains("Br|Cl", regex=True)]
            if any(_peak_near(mzs, float(t["mz"]) + _D13C) is not None
                   for _, t in twins.iterrows()):
                continue
            # cross-channel fallback: the SAME neutral independently assigned
            # High/Good on another peak (other adduct) is positive evidence
            # that outweighs one absent satellite -- an absent 13C can be a
            # peak-picker loss, an agreeing second channel cannot. (v21
            # cleared five sub-ppm [M+Br]- partners of Good [M-H]-
            # assignments, e.g. C10H16O6 at 311.013 / 2.4k cps.)
            others = ledger[(ledger["role"] == L.ROLE_M0)
                            & (ledger["peak_id"] != r["peak_id"])
                            & (ledger["neutral_formula"] == r["neutral_formula"])
                            & ledger["confidence"].astype(str)
                            .str.startswith(("High", "Good"))]
            if len(others):
                continue
            try:
                L.clear_assignment(
                    ledger, r["peak_id"],
                    reason=f"isotope audit: predicted 13C satellite "
                           f"({expected:.0f} cps) absent from spectrum")
                out["c13_missing"] += 1
            except L.LedgerError:
                continue

    n = sum(out.values())
    if n:
        log(f"[audit] isotope physics: {out['doublet_child']} doublet twins "
            f"re-attached, {out['doublet_cleared']} no-Br doublet formulas "
            f"cleared, {out['c13_attached']} 13C satellites attached, "
            f"{out['c13_clamp']} carbon-clamp clears, "
            f"{out['c13_missing']} missing-13C clears")
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
    parents_with_kids = set(ledger.loc[ledger["role"] == L.ROLE_ISO,
                                       "parent_peak_id"].dropna())
    for _, r in m0.iterrows():
        weak = not str(r["confidence"]).startswith(("High", "Good"))
        has_kids = r["peak_id"] in parents_with_kids
        z = z_of(r["ppm_error"], cfg)
        try:
            if z is None:
                if pd.isna(r["ppm_error"]) and weak and not has_kids:
                    L.clear_assignment(ledger, r["peak_id"],
                                       reason="mass-gate: no ppm error")
                    out["cleared_nan"] += 1
            elif z > cfg.cal_z_pattern:
                L.clear_assignment(ledger, r["peak_id"],
                                   reason=f"mass-gate: z={z:.1f} > {cfg.cal_z_pattern}")
                out["cleared_z"] += 1
            elif z > cfg.cal_z_accept and weak and not has_kids:
                L.clear_assignment(
                    ledger, r["peak_id"],
                    reason=f"mass-gate: z={z:.1f} without pattern evidence")
                out["cleared_z_noiso"] += 1
        except L.LedgerError:
            continue
    n = sum(out.values())
    if n:
        log(f"[audit] mass gate cleared {n} assignments "
            f"(z>{cfg.cal_z_pattern}: {out['cleared_z']}, "
            f"{cfg.cal_z_accept}<z<={cfg.cal_z_pattern} no-evidence: "
            f"{out['cleared_z_noiso']}, no-ppm: {out['cleared_nan']})")
    return out


def _prefer_adduct_reading(w, cfg: PassConfig):
    """Relabel a winner whose NEUTRAL carries the reagent halogen so the Br
    sits in the adduct/cluster, not the neutral (user reagent rule). With
    Y = X - HBr (structurally valid):
      X(Br) [M-H]-       -> Y [M+Br]-          (deprotonation == Br adduct)
      X(Br) [M+Br]-      -> Y [M+HBr+Br]-       (covalent == HBr cluster)
      X(Br) [M+<chan>]-  -> Y [M+HBr+<chan>]-   (covalent+air-ion == HBr cluster
                                                 on the same background channel,
                                                 e.g. the 426.976 CO3 case)
    The ion formula, score and ppm are unchanged -- only the decomposition is.
    A relabel onto a cluster adduct only fires when that adduct's exact mass is
    registered, so we never invent an unmodelled channel."""
    el = cfg.reagent_element
    if not el:
        return w
    x = w["neutral"]
    if C.parse_formula(x).get(el, 0) < 1:
        return w
    hx = "H" + el
    if hx not in G.REPEAT_UNITS:
        return w
    y = G.formula_add(x, hx, -1)
    if not y or not C.dbe_ok(y)[0] or not C.oxygen_ok(y)[0] \
            or C.parse_formula(y).get("C", 0) < 1:
        return w
    adduct = w["adduct"]
    if adduct == "[M-H]-":
        new_adduct = f"[M+{el}]-"
    elif adduct == f"[M+{el}]-":
        new_adduct = f"[M+{hx}+{el}]-"
    elif adduct.startswith("[M+") and adduct.endswith("]-") and hx not in adduct:
        # background channel (CO3/O2/...): insert the HBr cluster unit
        new_adduct = f"[M+{hx}+{adduct[3:-2]}]-"
    else:
        return w
    if new_adduct not in C.ADDUCT_SHIFTS:
        return w   # unmodelled cluster channel -> keep the covalent reading
    w = w.copy()
    w["adduct"] = new_adduct
    w["neutral"] = y
    w["_relabel_note"] = (
        f" Ion identical to covalent {x} {adduct}; reagent-adduct reading "
        f"preferred ({el} assigned to the adduct/cluster, not the neutral).")
    return w


def _iso_list(kids: pd.DataFrame, parent_pid) -> list[dict]:
    if kids is None or len(kids) == 0:
        return []
    sub = kids[kids["parent_peak_id"] == parent_pid]
    return [{"label": r["iso_label"], "score": r["iso_score"], "peak_id": r["peak_id"]}
            for _, r in sub.iterrows()]


def _commentary(w, pass_no, method) -> str:
    base = (f"Pass {pass_no} ({method}): {w['neutral']} {w['adduct']}, "
            f"ion score {w['ion_score']:.2f}, ppm "
            f"{w['ppm_error']:.2f}" if w['ppm_error'] is not None
            else f"Pass {pass_no} ({method}): {w['neutral']} {w['adduct']}, "
                 f"ion score {w['ion_score']:.2f}")
    if w["n_iso"]:
        base += f"; {w['n_iso']} isotopologue(s) confirmed by Mascope"
    if w["alternatives"]:
        a = w["alternatives"][0]
        margin = (w["eff_score"] - a["eff_score"]) if (w["eff_score"] is not None
                  and a.get("eff_score") is not None) else None
        if margin is not None:
            base += (f". Nearest competitor {a['formula']} trails by {margin:.2f}"
                     + (" (TIE)" if w["tied"] else ""))
    note = w.get("_relabel_note") if hasattr(w, "get") else None
    if note:
        base += note
    return base


# ---------------------------------------------------------------------------
# Range building
# ---------------------------------------------------------------------------
def build_ranges(profile: X.ContextProfile, pre, *, include_N: bool,
                 extra_elements: dict[str, tuple[int, int]] | None = None,
                 o_max: int = 30, c_max: int = 40) -> dict[str, tuple[int, int]]:
    """Build a NEUTRAL-formula grid box.

    Pass 1/2 are CHO(N) only: heteroatoms are NOT auto-added from the (noisy)
    prescan -- they enter the neutral exclusively via `extra_elements` (Pass 3
    contaminant families). This is what prevents the [M+Br]- alias from being
    mis-read as a brominated neutral: in a Br-CIMS sample the Br lives in the
    ADDUCT, not the neutral. The prescan only caps C here.
    """
    cmax = c_max
    if pre is not None and getattr(pre, "estimated_max_C", 0):
        cmax = min(c_max, max(12, pre.estimated_max_C + 4))
    r = {"C": (0, cmax), "H": (0, cmax * 2 + 4), "O": (0, o_max),
         "N": (0, profile.max_N if include_N else 0),
         "S": (0, 0), "P": (0, 0), "Si": (0, 0),
         "F": (0, 0), "Cl": (0, 0), "Br": (0, 0), "I": (0, 0)}
    if extra_elements:
        for el, (lo, hi) in extra_elements.items():
            cap = getattr(profile, f"max_{el}", hi) or hi
            r[el] = (lo, min(hi, cap))
    return r


def ranges_to_string(r: dict[str, tuple[int, int]]) -> str:
    return " ".join(f"{el}{lo}-{hi}" for el, (lo, hi) in r.items() if hi > 0)


def _resolve_hx_clusters(client, sample_id: str, ledger: pd.DataFrame, profile,
                         cfg: PassConfig, reagent: str, hx: str, *, log=print) -> dict:
    """Explain unassigned peaks as anchor.HX clusters (Y.HBr.Br- etc.).

    For each anchor Y, the cluster composition X = Y+HX under [M+X]- is scored
    by Mascope (identical ion to the covalent alias), but committed with
    neutral = Y and adduct = '[M+HX+X]-' so the target list reports the real
    analyte, with commentary naming the cluster interpretation."""
    out = {"committed": 0, "locked": 0, "iso_attached": 0,
           "claimed_formulas": set()}
    anchors = ledger.loc[ledger["role"] == L.ROLE_M0, ["peak_id", "neutral_formula"]]
    anchor_by_formula = dict(zip(anchors["neutral_formula"], anchors["peak_id"]))
    if not anchor_by_formula:
        return out
    adduct = f"[M+{reagent}]-"
    if adduct not in C.ADDUCT_SHIFTS:
        return out
    # cluster bases: anchors plus their +/-CH2 homologs (GKA-validated bridge --
    # e.g. glutaric sits between anchored succinic and adipic acids).
    ys: dict[str, tuple[object, str]] = {
        y: (apid, "anchor") for y, apid in anchor_by_formula.items()}
    for y, apid in anchor_by_formula.items():
        for s in (+1, -1):
            y2 = G.formula_add(y, "CH2", s)
            if y2 and y2 not in ys:
                keep, _ = X.filter_by_context(y2, profile.label)
                if keep:
                    ys[y2] = (apid, f"homolog of anchor {y} ({s:+d}CH2)")
    tgt = _target_peaks(ledger, cfg)
    # propose: peak mz near ANY isotopologue line of Y+HX under [M+X]- is fine;
    # at minimum the base line within series ppm
    proposals: dict[str, tuple[str, object, str]] = {}  # X -> (Y, anchor_pid, note)
    tmz = sorted(tgt["mz"].tolist())
    import bisect as _bs
    for y, (apid, note) in ys.items():
        x = G.formula_add(y, hx, +1)
        if not x:
            continue
        theo = C.ion_mz(x, adduct)
        # accept proposal if a target peak sits at the base line OR the +2
        # heavy-isotope line (Br2 envelopes often have the base line weak)
        for line in (theo, theo + 1.99795):
            tol = line * cfg.series_ppm * 1e-6
            j = _bs.bisect_left(tmz, line - tol)
            if j < len(tmz) and tmz[j] <= line + tol:
                proposals[x] = (y, apid, note)
                break
    if not proposals:
        return out
    scored = IO.score_candidates(client, sample_id, sorted(proposals),
                                 mechanism_ids=cfg.mechanism_ids)
    if len(scored) == 0:
        return out
    n_committed = 0
    n_iso = 0
    want_x = scored["compound_formula"].isin(proposals)
    for (x, ion_f), grp in scored[want_x].groupby(["compound_formula", "ion_formula"]):
        # only the CLUSTER ion form (reagent count = neutral's + 1 from adduct)
        if C.parse_formula(ion_f).get(reagent, 0) != \
                C.parse_formula(x).get(reagent, 0) + 1:
            continue
        y, apid, note = proposals[x]
        brow = grp[grp["is_base"]].iloc[0] if grp["is_base"].any() else None
        if brow is None:
            continue
        ion_score = brow["ion_score"]
        if pd.isna(ion_score) or float(ion_score) < cfg.series_min_score:
            continue
        attributed_iso = grp[(~grp["is_base"]) & grp["sample_peak_id"].notna()
                             & (pd.to_numeric(grp["iso_score"], errors="coerce")
                                .fillna(0) > 0.4)]
        # target peak: base-line attribution preferred; for split Br2 envelopes
        # fall back to the ledger peak at the base line's theoretical m/z,
        # requiring at least one Mascope-attributed heavy isotopologue.
        pid = brow["sample_peak_id"]
        ppm_err = _f(brow["ppm_error"])
        envelope_note = ""
        if pid is None or pd.isna(pid):
            if len(attributed_iso) == 0:
                continue
            theo0 = float(brow["theo_mz"])
            tol = theo0 * cfg.series_ppm * 1e-6
            cand = tgt[(tgt["mz"] - theo0).abs() <= tol]
            if len(cand) == 0:
                continue
            pid = cand.sort_values("height", ascending=False)["peak_id"].iloc[0]
            ppm_err = float((cand["mz"].iloc[0] - theo0) / theo0 * 1e6)
            envelope_note = (" Base isotopologue attribution recovered from the "
                             "heavy-isotope line (split halogen envelope).")
        try:
            if L.is_locked(ledger, pid) or L.role_of(ledger, pid) != L.ROLE_UNEXPLAINED:
                continue
            score = float(ion_score)
            conf = confidence_label(score, ppm_err, len(attributed_iso), False,
                                    cfg, suffix=f"{hx}-cluster")
            if conf == "Reject":
                continue
            L.commit_assignment(
                ledger, pid, neutral_formula=y, adduct=f"[M+{hx}+{reagent}]-",
                ion_formula=ion_f, ion_score=score,
                compound_score=_f(brow["compound_score"]), ppm_error=ppm_err,
                pass_no=3, method=f"cluster:{hx}", confidence=conf,
                commentary=(f"Pass 3 (cluster): {hx} cluster of {y} ({note}, "
                            f"ref peak {apid}); ion {ion_f} scored {score:.2f} "
                            f"by Mascope. Composition identical to covalent {x};"
                            f" cluster reading preferred.{envelope_note}"),
                anchor_peak_id=apid, series_unit=hx)
            out["claimed_formulas"].add(x)
            n_committed += 1
            for _, k in attributed_iso.iterrows():
                kp = k["sample_peak_id"]
                if kp == pid:
                    continue
                try:
                    L.attach_isotopologue(ledger, kp, pid,
                                          iso_label=k["iso_label"],
                                          iso_match_score=_f(k["iso_score"]))
                    n_iso += 1
                except L.LedgerError:
                    continue
        except L.LedgerError:
            continue
    out["committed"] = n_committed
    out["iso_attached"] = n_iso
    log(f"[pass3:cluster-{hx}] {{'committed': {n_committed}, 'iso_attached': {n_iso}}}")
    return out


# ---------------------------------------------------------------------------
# Candidate enumeration helpers
# ---------------------------------------------------------------------------
def _target_peaks(ledger: pd.DataFrame, cfg: PassConfig) -> pd.DataFrame:
    un = L.unassigned_peaks(ledger)
    return un[un["height"].fillna(0) >= cfg.height_cutoff]


def _family_ok(formula: str, ranges: dict[str, tuple[int, int]]) -> bool:
    """Structural gates (integer DBE + Senior + oxygen cap) plus the family's
    element ceilings. Used for Pass-3 contaminant families, where the family's
    ranges -- not the context caps -- are the elemental authority."""
    ok, _ = C.dbe_ok(formula)
    if not ok:
        return False
    ok, _ = C.oxygen_ok(formula)
    if not ok:
        return False
    cnt = C.parse_formula(formula)
    for el, n in cnt.items():
        lo, hi = ranges.get(el, (0, 0))
        if n > hi:
            return False
    return True


def _context_filter(formulas, context: str) -> list[str]:
    out = []
    for f in formulas:
        keep, _ = X.filter_by_context(f, context)
        if keep:
            out.append(f)
    return out


def _enumerate(client, mzs, mech_ids, ranges: dict, cfg: PassConfig,
               adducts: list[str], *, use_grid: bool = True) -> set[str]:
    """Candidate NEUTRAL formulas for these m/z. The local grid is primary
    (complete for CHO/CHON in-range, never fails); cheminfo is best-effort and
    only consulted when cfg.use_cheminfo is set."""
    formulas: set[str] = set()
    if use_grid:
        gadducts = [a for a in adducts if a in C.ADDUCT_SHIFTS]
        formulas.update(C.candidates_for_peaks(
            list(mzs), ranges, gadducts, ppm_tolerance=cfg.search_ppm))
    if cfg.use_cheminfo and mech_ids:
        rng_str = ranges_to_string(ranges)
        if rng_str:
            bulk = IO.query_candidates_bulk(
                client, list(mzs), mech_ids, formula_ranges=rng_str,
                ppm=cfg.search_ppm, limit=cfg.limit_per_peak, workers=cfg.workers)
            for cands in bulk.values():
                formulas.update(cands)
    return formulas


def _mech_ids_for(client, adducts: list[str]) -> list[str]:
    names = [IO.ADDUCT_TO_MECH[a] for a in adducts if a in IO.ADDUCT_TO_MECH]
    return list(IO.resolve_mechanism_ids(client, names).values())


# ---------------------------------------------------------------------------
# Pass 1 -- lock the CHO / CHON backbone
# ---------------------------------------------------------------------------
def run_pass1(client, sample_id: str, ledger: pd.DataFrame, profile, pre,
              cfg: PassConfig, adducts: list[str], *, log=print) -> dict:
    tgt = _target_peaks(ledger, cfg)
    mzs = tgt["mz"].tolist()
    mech_ids = _mech_ids_for(client, adducts)
    log(f"[pass1] {len(mzs)} target peaks; adducts={adducts}")
    # single CHO+CHON enumeration; arbitration's complexity penalty handles the
    # CHO-before-CHON preference, so no need for two separate sub-passes.
    ranges = build_ranges(profile, pre, include_N=True)
    formulas = _enumerate(client, mzs, mech_ids, ranges, cfg, adducts)
    formulas = set(_context_filter(formulas, profile.label))
    log(f"[pass1] {len(formulas)} context-plausible CHO/CHON candidate formulas")
    scored = IO.score_candidates(client, sample_id, sorted(formulas),
                                 mechanism_ids=cfg.mechanism_ids)
    log(f"[pass1] scored rows={len(scored)}")
    arb = arbitrate(scored, cfg)
    summary = commit_winners(ledger, arb, pass_no=1, method="cheminfo+grid",
                             context=profile.label, cfg=cfg, lock=True,
                             min_raw_score=cfg.tau_low)
    log(f"[pass1] {summary}")
    return summary


# ---------------------------------------------------------------------------
# Pass 2 -- GKA series expansion from locked anchors
# ---------------------------------------------------------------------------
def run_pass2(client, sample_id: str, ledger: pd.DataFrame, profile,
              cfg: PassConfig, adducts: list[str], *, log=print) -> dict:
    units = tuple(G.ORGANIC_UNITS)
    if "siloxane" in profile.pass3_families:
        units = units + ("C2H6OSi",)
    if "fluorinated" in profile.pass3_families:
        units = units + ("CF2",)
    total = {"committed": 0, "locked": 0, "iso_attached": 0}
    tried: set[str] = set()
    # Iterative GKA: each round, confirmed members (incl. last round's) act as
    # anchors, so homologous series are walked outward step by step.
    for it in range(cfg.series_max_iter):
        anchors = set(ledger.loc[ledger["role"] == L.ROLE_M0, "neutral_formula"].dropna())
        if not anchors:
            break
        tgt = _target_peaks(ledger, cfg)
        proposals: set[str] = set()
        for mz in tgt["mz"]:
            for p in G.propose_for_peak(mz, anchors, adducts, units=units,
                                        ppm=cfg.series_ppm, max_steps=1):
                proposals.add(p.neutral_formula)
        proposals = (set(_context_filter(proposals, profile.label))
                     - anchors - tried)
        if not proposals:
            log(f"[pass2.{it}] no new proposals; stopping")
            break
        tried |= proposals
        scored = IO.score_candidates(client, sample_id, sorted(proposals),
                                     mechanism_ids=cfg.mechanism_ids)
        arb = arbitrate(scored, cfg)
        s = commit_winners(ledger, arb, pass_no=2, method="gka-series",
                           context=profile.label, cfg=cfg, lock=False,
                           min_raw_score=cfg.series_min_score,
                           confidence_suffix="series",
                           claim_unexplained_only=True)
        for k in total:
            total[k] += s[k]
        log(f"[pass2.{it}] {len(anchors)} anchors -> {len(proposals)} proposals -> {s}")
        if s["committed"] == 0:
            break
    log(f"[pass2] total {total}")
    return total


# ---------------------------------------------------------------------------
# Pass 3 -- contaminant / low-quality recovery
# ---------------------------------------------------------------------------
def run_pass3(client, sample_id: str, ledger: pd.DataFrame, profile, pre,
              cfg: PassConfig, adducts: list[str], *, log=print) -> dict:
    tgt = _target_peaks(ledger, cfg)
    if len(tgt) == 0:
        log("[pass3] nothing unassigned; skipping")
        return {"committed": 0, "locked": 0, "iso_attached": 0}
    mzs = tgt["mz"].tolist()
    total = {"committed": 0, "locked": 0, "iso_attached": 0}
    from . import reagents as _RG
    reagent = _RG.reagent_for_adducts(adducts)

    # --- HX-cluster resolution (halide CIMS) -------------------------------
    # A peak at anchor+HX under [M+X]- is the analyte's HX cluster
    # (Y.HX.X-), NOT a new covalent organohalogen. Resolve these against the
    # anchors FIRST so the organohalogen family below never sees them.
    cluster_claimed: set[str] = set()
    if reagent in ("Br", "Cl"):
        hx = "H" + reagent
        s = _resolve_hx_clusters(client, sample_id, ledger, profile, cfg,
                                 reagent, hx, log=log)
        for k in total:
            total[k] += s[k]
        cluster_claimed = s.get("claimed_formulas", set())

    # halide-CIMS: also try covalent organohalogens. The arbitration keeps the
    # complexity prior on the reagent element (its ion isotope can't prove
    # neutral ownership), so these only win with a real score margin.
    families = list(profile.pass3_families)
    if reagent == "Br" and "bromo_organic" not in families:
        families.append("bromo_organic")
    if reagent == "Cl" and "chloro_organic" not in families:
        families.append("chloro_organic")

    # --- automatic GKA series detection (the machine 'rotating plot') -------
    # Repeat-unit structure in the residual opens the matching contaminant
    # family even when the context has it off (e.g. CF2 links -> fluorinated).
    from . import series_detect as SD
    evidence = SD.detect_series(ledger, ppm=cfg.search_ppm,
                                min_height=cfg.height_cutoff)
    log("[pass3] series evidence: " + ", ".join(
        f"{r.unit}:{r.n_links}x{r.enrichment}" + ("*" if r.significant else "")
        for r in evidence.itertuples()))
    fam_members: dict[str, set] = {}
    for r in evidence.itertuples():
        if r.significant and r.action:
            fam_members.setdefault(r.action, set()).update(
                SD.unit_members(ledger, r.mass, ppm=cfg.search_ppm,
                                min_height=cfg.height_cutoff))
    for fam in SD.families_from_evidence(evidence):
        if fam not in families:
            families.append(fam)
            log(f"[pass3] GKA evidence opened family: {fam} "
                f"({len(fam_members.get(fam, []))} chain-member targets)")
    anchors_now = set(ledger.loc[ledger["role"] == L.ROLE_M0,
                                 "neutral_formula"].dropna())
    for fam_key in families:
        fam = X.CONTAMINANT_FAMILIES.get(fam_key)
        if not fam:
            continue
        try:
            ranges = build_ranges(profile, pre, include_N=True,
                                  extra_elements=fam["add"])
            # family-specific adducts unioned with the sample's reagent adducts
            fam_adducts = list(dict.fromkeys(
                [a for a in fam["adducts"] if a in C.ADDUCT_SHIFTS] + adducts))
            mech_ids = _mech_ids_for(client, fam_adducts)
            if fam_key in fam_members:
                # CHAIN-BASED generation for evidence-opened families: grid the
                # chain HEADS only, then propagate arithmetically along each
                # detected chain. Cheaper than gridding every member, and it
                # imposes series consistency -- a member's formula must be its
                # neighbour's formula +/- the unit.
                link_ppm = max(2.0, 2 * cfg.ppm)
                formulas = set()
                for r in evidence.itertuples():
                    if r.action != fam_key or not r.significant \
                            or r.unit not in G.REPEAT_UNITS:
                        continue
                    for chain in SD.unit_chains(ledger, r.mass, ppm=link_ppm,
                                                min_height=cfg.height_cutoff,
                                                min_len=2):
                        heads = C.candidates_for_peaks(
                            [chain[0][1]], ranges,
                            [a for a in fam_adducts if a in C.ADDUCT_SHIFTS],
                            ppm_tolerance=cfg.search_ppm)
                        for f0 in heads:
                            if not _family_ok(f0, ranges):
                                continue
                            f = f0
                            formulas.add(f)
                            for _k in range(1, len(chain)):
                                f = G.formula_add(f, r.unit, +1)
                                if f and _family_ok(f, ranges):
                                    formulas.add(f)
                                else:
                                    break
            else:
                formulas = _enumerate(client, mzs, mech_ids, ranges, cfg,
                                      fam_adducts)
            # Structural-only filtering is EARNED BY EVIDENCE: GKA-opened
            # families bypass the context caps (ambient max_F=0 would veto the
            # very candidates the detected CF2 chains justify), with the chain
            # membership + arbitration priors as the guard. Profile-default
            # families keep the full context filter -- without chain evidence,
            # dropping the ratio priors lets mass-fit junk flood in (v9/v10
            # lesson: 142 'amines').
            if fam_key in fam_members:
                formulas = {f for f in formulas if _family_ok(f, ranges)}
            else:
                formulas = set(_context_filter(formulas, profile.label))
            if fam_key in ("bromo_organic", "chloro_organic"):
                # drop covalent-X aliases of anchor.HX clusters: if stripping
                # one HX from X yields an existing anchor, the cluster reading
                # owns that composition.
                el = "Br" if fam_key == "bromo_organic" else "Cl"
                formulas = {f for f in formulas
                            if G.formula_add(f, "H" + el, -1) not in anchors_now}
                formulas -= cluster_claimed
            if not formulas:
                continue
            scored = IO.score_candidates(client, sample_id, sorted(formulas),
                                         mechanism_ids=cfg.mechanism_ids)
            arb = arbitrate(scored, cfg)
            s = commit_winners(ledger, arb, pass_no=3, method=f"contaminant:{fam_key}",
                               context=profile.label, cfg=cfg, lock=False,
                               min_raw_score=cfg.tau_suspect,
                               confidence_suffix=fam_key,
                               claim_unexplained_only=True,
                               only_peaks=fam_members.get(fam_key))
        except Exception as e:
            log(f"[pass3:{fam_key}] FAILED: {type(e).__name__}: {e}")
            continue
        for k in total:
            total[k] += s[k]
        log(f"[pass3:{fam_key}] {s}")
    log(f"[pass3] total {total}")
    total["series_evidence"] = evidence.to_dict("records")
    return total


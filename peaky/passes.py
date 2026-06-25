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

import re
from dataclasses import dataclass

import pandas as pd

from . import chemistry as C
from . import contexts as X
from . import io_mascope as IO
from . import isotopes as ISO
from . import ledger as L
from . import series_gka as G

__version__ = "0.8.0"  # offset-tolerant gates + prior_offset + calibration-aware arbitration + carbon-clamp Si-skip


@dataclass
class PassConfig:
    ppm: float = 1.0                 # user m/z trust
    # Grid-enumeration tolerance. The measured instrument accuracy is
    # sigma~0.35 ppm (self-calibration), so 5 ppm was ~14 sigma and enumerated
    # a large candidate cloud the calibrated z-gate then rejected -- pure wasted
    # scoring (and bigger match_compounds requests that time out on a flaky
    # server). 3 ppm is still ~8 sigma: safely past local calibration drift,
    # but ~1.7x fewer candidates. Enumeration only (a formula never gridded
    # can never be scored); match_compounds keeps its 5 ppm window so it still
    # attributes real 29Si/81Br satellites, and the z-gate owns ppm rejection.
    search_ppm: float = 3.0          # grid enumeration tolerance
    height_cutoff: float = 100.0
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
    # rough mass offset (ppm) seeded from the sample's own matches BEFORE the
    # pass-1 self-calibration, so the pre-calibration pass-0 known-species gate is
    # not blind to a large systematic instrument offset (set by assign.run).
    prior_offset: float = 0.0
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
    # Reference-list selection prior: a candidate neutral on an ACTIVE reference
    # peaklist (a published product of the sample's chemistry, or a known
    # contaminant) is far more likely real than a mass-coincidence monster of
    # similar score. Add a small TIE-BREAK bonus to its eff_score -- enough to win
    # a near-tie (gap < the 0.05 tie window), never enough to override a clearly
    # better isotope-scored fit. Empty set / 0.0 -> no-op. Set by assign.run from
    # the run's context-active reference lists (reflists.active_lists).
    reflist_formulas: frozenset = frozenset()
    reflist_prior: float = 0.04


# ---------------------------------------------------------------------------
# Confidence
# ---------------------------------------------------------------------------
def confidence_label(score: float, ppm: float | None, n_iso: int, tied: bool,
                     cfg: PassConfig, suffix: str = "") -> str:
    # ppm proximity is judged against the CALIBRATED mass center when known, so a
    # uniform instrument offset (e.g. the -2.4 ppm of the uronium source) does not
    # read as "off-mass" and cap every commit at Low. Pre-calibration commits
    # (pass 1, cal_mu still None) use 0 and are re-graded by relabel_confidence
    # once calibrate() has fitted the center.
    center = cfg.cal_mu if cfg.cal_mu is not None else 0.0
    a = abs(ppm - center) if ppm is not None and pd.notna(ppm) else 99.0
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

# arbitration penalty per sigma that a candidate's ppm sits BEYOND the calibrated
# accept band (offset-aware; 0 when uncalibrated). Modest, so it only overturns a
# genuinely off-trend mass-coincidence, not an on-trend near-tie.
CAL_ARB_WEIGHT = 0.04


def calibrate(ledger: pd.DataFrame, cfg: PassConfig, *, log=print) -> tuple | None:
    """Fit the instrument's real mass accuracy (mu, sigma of the ppm error)
    over the committed High/Good CHO-CHON backbone and store it on cfg.

    Robust fit: median + scaled MAD, so a few bad commits can't widen the
    window they would then pass through. Returns (mu, sigma, n) or None when
    the backbone is too small to trust."""
    m0 = ledger[(ledger["role"] == L.ROLE_M0) & ledger["ppm_error"].notna()].copy()
    # Select the backbone by SCORE, not the confidence LABEL: the label embeds a
    # |ppm|<=2 sub-gate centered on 0, so at a large systematic mass offset (e.g.
    # the -2.4 ppm uronium source) it would exclude exactly the genuine backbone
    # and bias mu toward 0. Score is offset-independent; the median+MAD then finds
    # the true center robustly. (For a well-calibrated source the two selections
    # coincide.)
    score = pd.to_numeric(m0.get("ion_score"), errors="coerce").fillna(0.0)
    chon = m0["neutral_formula"].astype(str).map(
        lambda f: bool(C.parse_formula(f))
        and set(C.parse_formula(f)) <= _BACKBONE_ELEMENTS)
    ppm = m0.loc[(score >= cfg.tau_good) & chon, "ppm_error"].astype(float)
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


def _conf_suffix(conf) -> str:
    """The parenthetical method tag of a confidence label ('Good (series)' ->
    'series'); '' when there is none."""
    m = re.search(r"\(([^)]+)\)", str(conf or ""))
    return m.group(1) if m else ""


def relabel_confidence(ledger: pd.DataFrame, cfg: PassConfig, *, log=print) -> int:
    """Re-grade committed-M0 confidence labels against the calibrated mass center.

    Pass 1 commits BEFORE calibrate() runs, so at a large systematic mass offset
    its labels -- judged against 0 ppm -- read the whole high-score backbone as
    'Low', which then caps the report tier at Candidate and (because the tier
    engine recalibrates off High/Good rows) starves the mass-error gate. Re-judged
    against cal_mu the backbone recovers its true High/Good grade, while an
    off-trend mass monster that looked Good near 0 is correctly demoted. The
    method suffix (series / siloxane / ...) is preserved. No-op when uncalibrated.
    """
    if cfg.cal_mu is None:
        return 0
    # only the UNLOCKED pass-1 backbone: locked commits (pass-0 known species,
    # the siloxane ladder, any pass-1 High) carry a DELIBERATE grade -- e.g. a
    # known composite contaminant (silanediol n=4, ~35% Mascope score) is locked
    # by the known-species privilege and must not be re-graded to Reject here.
    m0 = ledger[(ledger["role"] == L.ROLE_M0) & ~ledger["locked"].astype(bool)]
    kids = (ledger[ledger["role"] == L.ROLE_ISO]
            .groupby("parent_peak_id")["peak_id"].nunique())
    n = 0
    for i, r in m0.iterrows():
        vals = [float(v) for v in (r.get("ion_score"), r.get("compound_score"))
                if v is not None and pd.notna(v)]
        score = min(vals) if vals else 0.0
        n_iso = int(kids.get(r["peak_id"], 0))
        tied = bool(r.get("tied")) if pd.notna(r.get("tied")) else False
        new = confidence_label(score, r.get("ppm_error"), n_iso, tied, cfg,
                               _conf_suffix(r.get("confidence")))
        if new != str(r.get("confidence")):
            ledger.at[i, "confidence"] = new
            n += 1
    if n:
        log(f"[relabel] re-graded {n} confidence labels to the calibrated center "
            f"(mu={cfg.cal_mu:+.2f} ppm, sigma={cfg.cal_sigma:.2f})")
    return n


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
    # Calibration-aware off-trend penalty: once the mass calibration is known, a
    # candidate whose ppm sits far from the calibrated center is unlikely to be
    # the true formula even if its raw (OFFSET-BLIND) oracle score is high (the
    # server scores |ppm| vs theoretical, not vs the instrument's real center).
    # Without this, at a large systematic offset a mass-coincidence nearer 0 ppm
    # out-scores the true on-trend formula, WINS arbitration, then is z-rejected
    # at commit -- leaving the peak unexplained instead of taking the on-trend
    # formula (the mass-degenerate high-Si PDMS failure at the uronium -2.45 ppm
    # offset). Pre-calibration (pass 1) cal_mu is None -> no penalty.
    def _cal_offtrend(ppm):
        z = z_of(ppm, cfg)
        return 0.0 if z is None else max(0.0, z - cfg.cal_z_accept) * CAL_ARB_WEIGHT
    base["cal_penalty"] = base["ppm_error"].map(_cal_offtrend)
    base["eff_score"] = (base["raw_score"] - base["penalty"]
                         - base["adduct_label"].isin(cfg.minor_channels)
                         * cfg.minor_channel_penalty
                         - base["cal_penalty"])
    # reference-list selection prior: a neutral on an active reference peaklist gets
    # a small tie-break bonus so a known literature/contaminant formula wins a
    # near-tie over a mass-coincidence monster (a soft prior, not an override).
    if cfg.reflist_formulas and cfg.reflist_prior:
        base["eff_score"] += (base["compound_formula"].isin(cfg.reflist_formulas)
                              * cfg.reflist_prior)

    winners = []
    for pid, grp in base.groupby("sample_peak_id"):
        grp = grp.sort_values("eff_score", ascending=False)
        top = grp.iloc[0]
        n_iso = int(iso_count.get((top["compound_formula"], top["ion_formula"]), 0))
        runner_eff = float(grp.iloc[1]["eff_score"]) if len(grp) > 1 else None
        tied = runner_eff is not None and (float(top["eff_score"]) - runner_eff) < 0.05
        # keep up to 6 alternatives: candidate DENSITY is the report's
        # confidence currency (tiers.py), and a 3-deep list saturates too early
        alts = []
        for _, r in grp.iloc[1:7].iterrows():
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
            "eff_margin": (None if runner_eff is None
                           else float(top["eff_score"]) - runner_eff),
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
    # protonated-urea (uronium) adduct: ion = neutral + CH4N2O + H. Without this
    # entry the urea-channel assignments fall through to the "[M-H]-" default and
    # are mislabeled (positive-mode urea-CIMS: ~1/3 of the backbone is this
    # channel). diff sorts alphabetically C,H,N,O.
    (("C", 1), ("H", 5), ("N", 2), ("O", 1)): "[M+(CH4N2O)H]+",
}


def _mech_to_adduct(row) -> str:
    """Adduct label from the exact ion-vs-compound element difference."""
    ion = str(row.get("ion_formula") or "")
    ci = C.parse_formula(ion)
    cc = C.parse_formula(str(row.get("compound_formula") or ""))
    diff = tuple(sorted(
        (el, ci.get(el, 0) - cc.get(el, 0))
        for el in set(ci) | set(cc) if ci.get(el, 0) != cc.get(el, 0)))
    add = _DIFF_TO_ADDUCT.get(diff, "[M-H]-")
    # The element diff cannot see isotopic labelling: a ¹⁵N-nitrate reagent cluster
    # has the SAME (N+1, O+3) diff as a ¹⁴N one, but the server writes the heavy N
    # as '^N' in the ion formula. Without this, the ¹⁵N adduct is labelled [M+NO3]-
    # and the +61.99 (¹⁴N) shift puts ion_mz / jitter ~1 Da off.
    if add == "[M+NO3]-" and "^N" in ion:
        return "[M+^NO3]-"
    return add


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
                eff_score=w.get("eff_score"), eff_margin=w.get("eff_margin"),
                tied=w.get("tied"),
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


# Pass-0 known-species registry: explicit compositions that the general
# organic grid CANNOT reach -- either instrument contaminants, or real small
# molecules excluded by the integer-DBE / C>=1 organic priors. Each is scored
# by Mascope and LOCKED before pass 1, but still gated on mass + 81Br twin.
#   family -> (formula -> human label)
def _known_species(polarity: str = "negative") -> dict:
    # The known-species privilege is reagent/polarity-specific. The lists below
    # are NEGATIVE-mode (Br/halide-CIMS): small atmospheric acids/radicals seen
    # as Br- adducts, [M-H]- nitroaromatics, and the silanediol [M+Br]-/[M-H]-
    # contaminant series. In POSITIVE mode (urea-CIMS) none of these apply -- the
    # N-base / oxygenated-VOC analytes are reachable by the organic grid, and the
    # silanediol series would be scored under the wrong (anion) ion form -- so
    # pass 0 is a no-op and the grid + pass-3 families carry the sample.
    if str(polarity) == "positive":
        # Positive (urea-CIMS): the N-base / oxygenated-VOC analytes are reachable
        # by the organic grid, so most of pass 0 is a no-op. EXCEPTION:
        # ORGANOPHOSPHATE esters / phosphine oxides -- ubiquitous lab & indoor
        # contaminants that ionise well as [M+H]+ / [M+(urea)H]+ but are INVISIBLE
        # to the CHNOS grid (P is off by default; opening the P grid floods mass-
        # degeneracy with P2/P3, N-rich monsters). Supply them as explicit known
        # formulas. P is monoisotopic -> no isotope twin to confirm, so the commit
        # is gated on CROSS-CHANNEL corroboration (>=2 ion channels) in pass 0.
        organophosphate = {
            "C6H15O4P":  "triethyl phosphate (TEP)",
            "C9H21O4P":  "tripropyl phosphate (TPrP)",
            "C12H27O4P": "tri-n-butyl phosphate (TBP / TiBP)",
            "C18H15O4P": "triphenyl phosphate (TPhP)",
            "C18H15OP":  "triphenylphosphine oxide (TPPO)",
            "C18H39O7P": "tris(2-butoxyethyl) phosphate (TBEP)",
            "C24H51O4P": "tris(2-ethylhexyl) phosphate (TEHP)",
            "C21H21O4P": "tricresyl phosphate (TMPP / TCrP)",
        }
        return {"organophosphate": organophosphate}
    atmos = {
        # small atmospheric acids / radicals detected as Br- adducts -- the
        # PRIMARY analytes of a Br-CIMS, all invisible to the organic grid:
        # HO2 is a radical (half-integer DBE); the rest are C0 inorganics.
        # HNO3/HNO2 confirmed AMBIENT analytes by the user (2026-06-12),
        # so they were removed from the reagent-cluster library.
        "HO2":  "hydroperoxyl radical",
        "HNO3": "nitric acid",
        "HNO2": "nitrous acid",
        "HNO4": "peroxynitric acid",
    }
    # Atmospheric nitroaromatics (brown-carbon tracers from NOx + aromatic VOC /
    # biomass burning), detected as [M-H]-. These are H-POOR / high-DBE, so the
    # ambient Van Krevelen floor + DBE/C ceiling block them from the organic grid
    # -- they must be supplied as known tracers (the v45->v46 fix; dinitrophenol
    # was independently confirmed present by an Orbitool assignment). Only those
    # actually present + |ppm|<=2 commit, so listing absent ones is harmless.
    nitroaromatic = {
        "C6H4N2O5": "2,4-dinitrophenol",
        "C6H5NO3":  "nitrophenol",
        "C6H5NO4":  "nitrocatechol",
        "C7H6N2O5": "dinitrocresol",
    }
    contam = {f: "dimethylsilanediol oligomer (PDMS hydrolysis, inlet/tubing)"
              for f in _silanediol_series()}
    # Perfluorocarboxylic acids CnHF(2n-1)O2 -- ubiquitous environmental / lab
    # contaminants and a classic CIMS signal (TFA, PFPrA, PFBA, ... PFOA). They are
    # F-rich so the organic grid (max_F=0, or pass-4's clamped-F path) misses the
    # clean low-F ones (TFA was the brightest peak in the ¹⁵NO₃⁻ batch yet went
    # unexplained). The PFCA formula is HIGHLY specific (high negative mass defect),
    # so supply the series as known formulas; only those present + on-cal commit.
    # In a narrow high-m/z window TFA's [M-H]- (112.99) is out of range and it is
    # seen ONLY as the reagent adduct, so do NOT require the deprotonation channel.
    perfluoroacid = {f"C{n}HF{2 * n - 1}O2": f"perfluoro-C{n} acid (PFCA)"
                     for n in range(2, 13)}
    # Chlorinated paraffins (SCCP/MCCP/LCCP): saturated CnH(2n+2-x)Clx, the
    # persistent-organic-pollutant family the user's screenshot showed (C10H17Cl5 ..
    # C14H22Cl8). Cl 3-15 is FAR above the organic grid's max_Cl (<=2), so they are
    # never enumerated -> they land in 'unexplained'. Supply the (tight, 2-parameter)
    # family as known formulas; they commit ONLY with a confirmed ³⁷Cl envelope, so
    # this is safe despite the wide n/x range. Listing absent ones is harmless.
    chlorinated_paraffin = {}
    for _n in range(10, 31):
        for _x in range(3, 16):
            _h = 2 * _n + 2 - _x
            if _h >= 1:
                chlorinated_paraffin[f"C{_n}H{_h}Cl{_x}"] = f"chlorinated paraffin C{_n}Cl{_x}"
    return {"atmospheric": atmos, "nitroaromatic": nitroaromatic,
            "perfluoroacid": perfluoroacid, "chlorinated_paraffin": chlorinated_paraffin,
            "contaminant:silanediol": contam}


# ³⁷Cl - ³⁵Cl mass spacing (for ledger-side ³⁷Cl-envelope confirmation).
_D37CL = 1.9970499
# Pass-0 known-species families whose identity is ISOTOPE-confirmable (a Cl/Br/S
# envelope), so a real isotope-laddered peak may be committed even when the
# server's aggregate compound_score is too low to anchor it (the ¹⁵N-phantom /
# wide-envelope score depression documented for ¹⁵N-labelled poly-Cl). The
# monoisotopic F/P families (perfluoroacid, organophosphate) are deliberately NOT
# recoverable -- they have no isotope twin to corroborate a mass coincidence.
_RECOVERABLE_KNOWN_FAMS = {"chlorinated_paraffin"}


def run_pass0_known(client, sample_id: str, ledger: pd.DataFrame,
                    profile, cfg: PassConfig, adducts: list[str], *,
                    score_fn=None, log=print) -> dict:
    """Pass 0 -- assign explicit KNOWN species (contaminant series + small
    atmospheric acids/radicals) before the organic passes run. Mascope scores
    the compositions (its isotope model covers 29Si/30Si + the reagent
    halogen); commits are LOCKED so pass 1 cannot displace them with grid CHO
    fits. Each still passes the mass gate (|ppm|<=2) and the 81Br-twin
    consistency check, so a composite collision is refused, not locked."""
    score_fn = score_fn or IO.score_candidates
    out = {"committed": 0, "locked": 0, "iso_attached": 0}
    registry = _known_species(getattr(profile, "polarity", "negative"))
    label_of = {f: (fam, lbl) for fam, d in registry.items()
                for f, lbl in d.items()}
    formulas = sorted(label_of)
    if not formulas:                       # positive mode: pass 0 is a no-op
        log(f"[pass0] no known-species list for polarity="
            f"{getattr(profile, 'polarity', 'negative')!r}; skipping")
        return out
    scored = score_fn(client, sample_id, formulas,
                      mechanism_ids=cfg.mechanism_ids)
    if scored is None or len(scored) == 0:
        log("[pass0] WARNING scoring returned EMPTY for the known-species list")
        out["scoring_empty"] = True
        return out
    base = scored[scored["is_base"] & scored["sample_peak_id"].notna()
                  & scored["ion_score"].notna()]
    # cross-channel corroboration count (on-cal): how many distinct ion channels
    # each known formula matches. Monoisotopic-P organophosphates require >=2.
    if "mechanism_id" in base.columns:
        _onc = base[(pd.to_numeric(base["ppm_error"], errors="coerce")
                     - cfg.prior_offset).abs() <= 2.0]
        ope_channels = _onc.groupby("compound_formula")["mechanism_id"].nunique().to_dict()
    else:
        ope_channels = {}
    kids = scored[(~scored["is_base"]) & scored["sample_peak_id"].notna()
                  & (pd.to_numeric(scored["iso_score"], errors="coerce")
                     .fillna(0) > 0.4)]
    mzs = ledger["mz"]
    for _, r in base.iterrows():
        ppm = r["ppm_error"]
        # the known-list privilege still requires the MASS, but judged against the
        # rough offset (cfg.prior_offset) so a uniformly-shifted instrument (e.g.
        # -1.9 ppm) doesn't drop on-trend contaminants and hand the peak to an
        # off-trend CHO mass-fit (the silanediol-vs-C5H10O6 collision at -1.9 ppm).
        if ppm is None or pd.isna(ppm) or abs(float(ppm) - cfg.prior_offset) > 2.0:
            continue
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
            fam, lbl = label_of[r["compound_formula"]]
            # organophosphates are monoisotopic in P -> require >=2 ion channels
            # (e.g. [M+H]+ AND [M+(urea)H]+) before locking, since there is no
            # isotope twin to confirm a single-channel mass coincidence.
            if fam == "organophosphate" and ope_channels.get(r["compound_formula"], 0) < 2:
                log(f"[pass0] skip {r['compound_formula']} @{float(r['sample_peak_mz']):.4f}: "
                    f"single ion channel (monoisotopic P needs >=2 to corroborate)")
                continue
            tag = ("atmospheric" if fam == "atmospheric"
                   else "nitroaromatic" if fam == "nitroaromatic"
                   else "organophosphate" if fam == "organophosphate"
                   else "perfluoroacid" if fam == "perfluoroacid"
                   else "chlorinated-paraffin" if fam == "chlorinated_paraffin"
                   else "contaminant")
            fam_kids = kids[kids["compound_formula"] == r["compound_formula"]]
            n_kids = int((fam_kids["sample_peak_id"] != pid).sum())
            # chlorinated paraffins (Cl is off the organic grid at Cl>2): commit ONLY
            # when the ³⁷Cl envelope is confirmed (>=2 matched ³⁷Cl satellites), so a
            # CnHmClx mass coincidence is rejected. Cl IS isotope-confirmable (unlike
            # monoisotopic F/P), and the server's aggregate score is artificially low
            # for ¹⁵N-labelled poly-Cl (¹⁴N phantoms + wide envelope), so we lock on
            # the isotope evidence instead of the depressed compound_score.
            if fam == "chlorinated_paraffin" and n_kids < 2:
                log(f"[pass0] skip {r['compound_formula']} @{float(r['sample_peak_mz']):.4f}: "
                    f"³⁷Cl envelope not confirmed (n_kids={n_kids})")
                continue
            # silanediol / any Si-rich known species: the 29Si M+1 must MATCH the Si
            # count, not merely exist. A high-O organic is mass-degenerate with a Si_k
            # oligomer; if the M+1 is too small the Si is over-claimed -> skip, leaving
            # the peak for the organic grid (the C10H18O11 vs C8H26O5Si4 case @393).
            _c0 = C.parse_formula(r["compound_formula"])
            if _c0.get("Si", 0) > 0:
                _ix = ledger.index[ledger["peak_id"] == pid]
                _m0h = (float(ledger.at[_ix[0], "height"])
                        if len(_ix) and pd.notna(ledger.at[_ix[0], "height"]) else 0.0)
                if not _si_m1_consistent(ledger, float(r["sample_peak_mz"]), _m0h,
                                         _c0.get("Si", 0), _c0.get("C", 0)):
                    log(f"[pass0] skip {r['compound_formula']} @{float(r['sample_peak_mz']):.4f}: "
                        f"29Si M+1 too small for Si{_c0.get('Si', 0)} (over-claimed; likely "
                        "a high-O organic) -- left for the grid")
                    out["si_underclaimed"] = out.get("si_underclaimed", 0) + 1
                    continue
            conf = (f"Good ({tag})" if float(r["ion_score"]) >= 0.7
                    or n_kids >= 2 else f"Low ({tag})")
            L.commit_assignment(
                ledger, pid, neutral_formula=r["compound_formula"],
                adduct=_mech_to_adduct(r), ion_formula=r["ion_formula"],
                ion_score=float(r["ion_score"]),
                compound_score=_f(r.get("compound_score")),
                ppm_error=float(ppm), pass_no=0,
                method=f"known:{fam}", confidence=conf,
                commentary=(f"Pass 0 (known {tag}): {r['compound_formula']} "
                            f"{_mech_to_adduct(r)} = {lbl}, ppm "
                            f"{float(ppm):.2f}, ion score {float(r['ion_score']):.2f}"
                            + ("; excluded from the organic grid (radical / C0 "
                               "inorganic)" if fam == "atmospheric"
                               else "; H-poor nitroaromatic blocked by the ambient "
                               "VK floor/DBE ceiling -- assigned as a known BrC tracer"
                               if fam == "nitroaromatic"
                               else "; organophosphate contaminant (P off the grid); "
                               f"corroborated across {ope_channels.get(r['compound_formula'], 0)} "
                               "ion channels (monoisotopic P, no isotope twin)"
                               if fam == "organophosphate"
                               else "; perfluorocarboxylic acid (F off the grid); "
                               "known PFCA series formula, exact-mass committed"
                               if fam == "perfluoroacid"
                               else "; chlorinated paraffin (Cl off the grid); ³⁷Cl "
                               f"envelope confirmed ({n_kids} satellites), isotope-locked"
                               if fam == "chlorinated_paraffin" else "")))
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
    # isotope-confirmed RECOVERY of known species the server scored too low to
    # anchor (e.g. ¹⁵N-labelled chlorinated paraffins, whose aggregate score
    # collapses so the base ion comes back UNANCHORED and the main loop above --
    # which iterates only server-anchored bases -- never sees it).
    rec = _recover_isotope_locked_known(ledger, scored, label_of, cfg, log=log)
    for k, v in rec.items():
        out[k] = out.get(k, 0) + v
    log(f"[pass0] {out}")
    return out


def _recover_isotope_locked_known(ledger: pd.DataFrame, scored: pd.DataFrame,
                                  label_of: dict, cfg: PassConfig, *,
                                  anchor_tol: float = 2.0, sat_ppm: float = 7.0,
                                  min_sats: int = 2, height_floor: float = 20.0,
                                  log=print) -> dict:
    """Recover an isotope-confirmable KNOWN species (Cl/Br/S family) whose server
    compound_score was too low to anchor a real peak.

    For ¹⁵N-labelled poly-Cl the server's aggregate match_score collapses (the
    ¹⁴N phantom lines + the wide Cl envelope drag it under possible_match_
    threshold), so ``match_compounds`` returns the base ion UNANCHORED
    (sample_peak_id / ppm NaN). ``run_pass0_known``'s main loop only iterates
    server-anchored bases, so those congeners never reach the ³⁷Cl confirmation
    and are left unexplained -- exactly the "too low score on the server" miss.

    Re-anchor against the LEDGER by exact mass (the server's theoretical M0 mz)
    and commit ONLY when BOTH hold:

      * a real, still-unexplained ledger peak sits within ``anchor_tol`` ppm of
        the theoretical M0 mass (offset-aware, like the main loop), AND
      * >= ``min_sats`` ³⁷Cl satellites (M0 + k·Δ³⁷Cl) are present in the ledger.

    This is the SAME isotope evidence the committed congeners passed; only the
    server's depressed aggregate score is bypassed. It cannot fabricate: no real
    peak, or no ³⁷Cl envelope, means no commit. Returns counter deltas
    (committed / locked / iso_attached / recovered) for the caller to fold in."""
    out = {"committed": 0, "locked": 0, "iso_attached": 0, "recovered": 0}
    if scored is None or len(scored) == 0 or "is_base" not in scored.columns:
        return out
    fam_of = {f: fam for f, (fam, _lbl) in label_of.items()}
    is_base = scored["is_base"].fillna(False)
    recoverable = scored["compound_formula"].map(
        lambda f: fam_of.get(f) in _RECOVERABLE_KNOWN_FAMS)
    bases = scored[is_base & recoverable]
    if not len(bases):
        return out
    mzs = ledger["mz"]
    for _, r in bases.iterrows():
        cf = r["compound_formula"]
        theo = r.get("theo_mz")
        if theo is None or pd.isna(theo):
            continue
        theo = float(theo)
        # offset-aware exact-mass anchor onto a still-unexplained ledger peak
        i = _peak_near(mzs, theo * (1 + cfg.prior_offset * 1e-6), ppm=anchor_tol)
        if i is None:
            continue
        pid = ledger.at[i, "peak_id"]
        try:
            if L.role_of(ledger, pid) != L.ROLE_UNEXPLAINED:
                continue
        except L.LedgerError:
            continue
        bh = ledger.at[i, "height"]
        if pd.isna(bh) or float(bh) < height_floor:
            continue
        led_ppm = (float(ledger.at[i, "mz"]) - theo) / theo * 1e6
        # ³⁷Cl envelope confirmation against the LEDGER (NOT the server iso_score,
        # which is itself depressed). Server satellite theo_mz preferred, with the
        # M0 + k·Δ³⁷Cl ladder as a fallback; count DISTINCT, still-unexplained
        # ledger peaks so a server/ladder mz that resolve to one peak count once.
        sib = scored[(scored["compound_formula"] == cf)
                     & (scored["ion_formula"] == r["ion_formula"])
                     & (~scored["is_base"].fillna(False))]
        sat_theos = {round(float(tm), 3) for tm in sib["theo_mz"].dropna()
                     if float(tm) > theo + 0.5}
        sat_theos.update(round(theo + k * _D37CL, 3) for k in (1, 2, 3, 4))
        sat_pids: list = []
        for tm in sorted(sat_theos):
            j = _peak_near(mzs, tm, ppm=sat_ppm)
            if j is None:
                continue
            jh = ledger.at[j, "height"]
            if pd.isna(jh) or float(jh) < height_floor * 0.5:
                continue
            jpid = ledger.at[j, "peak_id"]
            if jpid == pid or jpid in sat_pids:
                continue
            try:
                if L.role_of(ledger, jpid) == L.ROLE_UNEXPLAINED:
                    sat_pids.append(jpid)
            except L.LedgerError:
                continue
        if len(sat_pids) < min_sats:
            continue
        fam, lbl = label_of[cf]
        sc = _f(r.get("ion_score"))
        adduct = _mech_to_adduct(r)
        L.commit_assignment(
            ledger, pid, neutral_formula=cf, adduct=adduct,
            ion_formula=r["ion_formula"], ion_score=(0.0 if sc is None else float(sc)),
            compound_score=_f(r.get("compound_score")), ppm_error=float(led_ppm),
            pass_no=0, method=f"known:{fam}",
            confidence="Good (chlorinated-paraffin, recovered)",
            commentary=(f"Pass 0 (known chlorinated-paraffin, RECOVERED): {cf} "
                        f"{adduct} = {lbl}, ppm {led_ppm:.2f}; the server "
                        f"compound_score ({_f(r.get('compound_score'))}) was too low "
                        "to anchor a peak (¹⁵N-phantom / wide-envelope depression), "
                        "so it was exact-mass anchored to the ledger and "
                        f"isotope-locked on its ³⁷Cl envelope ({len(sat_pids)} "
                        "satellites)."))
        out["committed"] += 1
        out["recovered"] += 1
        L.lock_peaks(ledger, [pid])
        out["locked"] += 1
        for spid in sat_pids:
            try:
                L.attach_isotopologue(ledger, spid, pid, iso_label="37Cl")
                out["iso_attached"] += 1
            except L.LedgerError:
                continue
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
_D29SI = 0.999568          # 29Si - 28Si (the Si M+1)
_R29SI = 0.0468            # 29Si natural abundance per Si
SI_M1_MIN_FRAC = 0.6       # observed (M+1)/(M0) must be >= this * predicted Si M+1


def _peak_near(mzs: "pd.Series", target: float, ppm: float = 5.0):
    """Index of the closest ledger peak within ppm of target, else None."""
    tol = target * ppm * 1e-6
    d = (mzs - target).abs()
    i = d.idxmin()
    return i if d.loc[i] <= tol else None


def _si_m1_consistent(ledger: pd.DataFrame, m0_mz: float, m0_h: float,
                      n_si: int, n_c: int) -> bool:
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
    obs = (float(ledger.at[j, "height"]) / m0_h
           if j is not None and pd.notna(ledger.at[j, "height"]) else 0.0)
    return obs >= SI_M1_MIN_FRAC * pred


def complete_isotope_envelopes(ledger: pd.DataFrame, cfg: PassConfig, *,
                               min_rel: float = 0.06, ppm: float = 12.0,
                               log=print) -> dict:
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
    order = (ledger[ledger["role"] == L.ROLE_M0]
             .sort_values("mz")["peak_id"].tolist())
    for pid in order:
        idx = ledger.index[ledger["peak_id"] == pid]
        if not len(idx):
            continue
        i = idx[0]
        if str(ledger.at[i, "role"]) != L.ROLE_M0:
            continue                       # displaced by an earlier parent
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
                        L.attach_isotopologue(ledger, tpid, pid, iso_label=label,
                                              iso_match_score=score)
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
                if (not conf_j.startswith("High") and weak_score
                        and 0.45 <= ratio <= 2.2):
                    try:
                        L.displace_to_isotopologue(ledger, tpid, pid,
                                                   iso_label=label,
                                                   iso_match_score=score)
                        out["displaced"] += 1
                    except L.LedgerError:
                        pass
    if out["attached"] or out["displaced"]:
        log(f"[iso-envelope] attached {out['attached']} unclaimed satellites, "
            f"displaced {out['displaced']} mis-assigned satellites onto their "
            f"true parents")
    return out


# per-atom +1 (M+1) satellite abundance ratios -- all halogen-FREE, so the
# observed M+1 is dominated by the assigned compound even under a coincident
# halogen interferent (whose M+1 is weak); this is the composite discriminator
_M1_RATIO = {"C": 0.0107, "Si": 0.0508, "N": 0.003653, "S": 0.007896,
             "O": 0.000381, "H": 0.000115}


def detect_composites(ledger: pd.DataFrame, cfg: PassConfig, *,
                      min_m1_rel: float = 0.06, excess_frac: float = 0.25,
                      min_excess: float = 400.0, ppm: float = 8.0,
                      log=print) -> dict:
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
            continue                      # too few C/Si to diagnose a composite
        m0 = float(r["mz"]); h0 = float(r["height"])
        if not (h0 > 0):
            continue
        # observed M+1 region: 13C(+1.0034) + 29Si(+0.9996) + 15N(+0.997)
        h1 = _sum_window(m0 + 0.9940, m0 + 1.0070)
        if h1 <= 0:
            continue
        s_assigned = h1 / m1_rel          # implied true intensity of the M0 owner
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
        x2 = h2 - s_assigned * pat.get(2, 0.0)     # co-component M+2
        x4 = h4 - s_assigned * pat.get(4, 0.0)     # co-component M+4
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
            f"= ~{s_assigned:.0f} of {h0:.0f} cps")
        # structured fields for the de-blending step (split_composites)
        ledger.at[i, "assigned_fraction"] = max(0.0, min(1.0, s_assigned / h0))
        ledger.at[i, "co_height"] = excess
        ledger.at[i, "co_halogen"] = hal
        out["flagged"] += 1
    if out["flagged"]:
        log(f"[composite] flagged {out['flagged']} unresolved composite peaks "
            f"(M0 inflated beyond the M+1-implied owner intensity)")
    return out


def split_composites(ledger: pd.DataFrame, cfg: PassConfig, *,
                     log=print) -> pd.DataFrame:
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
    hosts = ledger[(ledger["role"] == L.ROLE_M0)
                   & ledger["co_height"].notna() & (~syn)]
    existing = set(ledger["peak_id"])
    new_rows = []
    for _, r in hosts.iterrows():
        co_h = float(r["co_height"])
        sub_id = f"{r['peak_id']}.2"
        if co_h < 1.0 or sub_id in existing:
            continue
        hal = str(r["co_halogen"])
        row = {c: pd.NA for c in ledger.columns}
        row.update({
            "peak_id": sub_id, "mz": float(r["mz"]), "height": co_h,
            "area": float("nan"), "role": L.ROLE_UNEXPLAINED, "synthetic": True,
            "host_peak_id": r["peak_id"], "assigned_fraction": 1.0,
            "locked": False, "co_height": float("nan"),
            "commentary": (f"co-eluting {hal} component split from peak "
                           f"{r['mz']:.4f} (~{co_h:.0f} cps); host owner "
                           f"{str(r['neutral_formula'])} keeps "
                           f"{100 * float(r['assigned_fraction']):.0f}%"),
        })
        new_rows.append(row)
        existing.add(sub_id)
    if not new_rows:
        return ledger
    add = pd.DataFrame(new_rows)[list(ledger.columns)]
    log(f"[composite] split {len(new_rows)} composite peaks into fractional "
        f"sub-peaks (owner keeps assigned_fraction; co-component -> '<id>.2')")
    return pd.concat([ledger, add], ignore_index=True)


def demote_carbon_inconsistent(ledger: pd.DataFrame, cfg: PassConfig, *,
                               log=print) -> int:
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
    for _, r in ledger[(ledger["role"] == L.ROLE_M0)
                       & ~ledger["locked"].astype(bool)].iterrows():
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
        k = ledger[(ledger["role"] == L.ROLE_ISO)
                   & (ledger["parent_peak_id"] == r["peak_id"])
                   & (ledger["iso_label"].astype(str) == "13C")]
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
                    ledger, r["peak_id"],
                    reason=f"carbon-clamp (pre-pass-4): 13C ratio measures ~C"
                           f"{c_est:.0f}, formula claims C{n_c}")
                n += 1
            except L.LedgerError:
                continue
    if n:
        log(f"[pre-pass4] carbon-clamp demoted {n} C-inconsistent M0 monsters "
            f"-> re-offered to the residual passes")
    return n


def demote_massgate_monsters(ledger: pd.DataFrame, cfg: PassConfig, *,
                             log=print) -> int:
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
    for _, r in ledger[(ledger["role"] == L.ROLE_M0)
                       & ~ledger["locked"].astype(bool)].iterrows():
        z = z_of(r["ppm_error"], cfg)
        if z is not None and z > cfg.cal_z_pattern:
            try:
                L.clear_assignment(
                    ledger, r["peak_id"],
                    reason=f"mass-gate (pre-pass-4): z={z:.1f} > {cfg.cal_z_pattern}")
                n += 1
            except L.LedgerError:
                continue
    if n:
        log(f"[pre-pass4] mass-gate demoted {n} z>{cfg.cal_z_pattern} monsters "
            f"-> re-offered to the residual passes")
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
                    # the lighter formula genuinely carries Br -> a ~1:1 twin
                    # 1.998 above IS its ⁸¹Br isotopologue (valid regardless of
                    # the reagent system).
                    L.clear_assignment(ledger, hv.peak_id,
                                       reason=f"isotope audit: 81Br twin of "
                                              f"{lt.mz:.4f} (ratio {hr:.2f})")
                    L.attach_isotopologue(ledger, hv.peak_id, lt.peak_id,
                                          iso_label="81Br")
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
            h_sat = float(k.iloc[0]["height"])
            c_est = (h_sat / float(r["height"])) / _R13C
            # clamp ONLY on a reliably-measured 13C satellite (>= the detection
            # floor). A sub-floor 13C ratio is noise and under-reads carbon,
            # which would falsely clear genuine low-intensity M0s (the ~2k cps
            # [M+15NO3]- compounds whose ~150 cps 13C sits near the floor). The
            # over-claim O-monster case always carries a BRIGHT 13C, so it fires.
            if n_c >= 8 and h_sat >= cfg.height_cutoff \
                    and abs(c_est - n_c) > max(2.5, 0.35 * n_c):
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
                 o_max: int | None = None, c_max: int | None = None
                 ) -> dict[str, tuple[int, int]]:
    """Build a NEUTRAL-formula grid box.

    Pass 1/2 are CHO(N) only: heteroatoms are NOT auto-added from the (noisy)
    prescan -- they enter the neutral exclusively via `extra_elements` (Pass 3
    contaminant families). This is what prevents the [M+Br]- alias from being
    mis-read as a brominated neutral: in a Br-CIMS sample the Br lives in the
    ADDUCT, not the neutral. The prescan only caps C here.

    The box width defaults to the context's grid_c_max / grid_o_max (40 / 30 for
    the ambient Br-CIMS profiles; wider for a heavier positive source like
    urea-CIMS). An explicit c_max / o_max argument overrides the profile.
    """
    if o_max is None:
        o_max = getattr(profile, "grid_o_max", 30)
    if c_max is None:
        c_max = getattr(profile, "grid_c_max", 40)
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
    if {"siloxane", "pdms"} & set(profile.pass3_families):
        units = units + ("C2H6OSi",)   # the PDMS dimethylsiloxane rung (+74.019)
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
        action = r.action
        if r.significant and isinstance(action, str) and action:
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

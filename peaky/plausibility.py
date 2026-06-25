"""Chemical-plausibility QC for assigned formulas.

High-resolution accurate mass alone lets the fitter hit almost any target once it
is allowed many heteroatoms (every extra N/O/S/halogen adds free parameters), so a
small fraction of assignments are mass-coincidence "monsters" rather than real
molecules in an organic-aerosol matrix: extreme heteroatom counts, implausibly
carbon-rich (very low H/C) skeletons, or a covalent halogen in a positive-mode
spectrum whose reagent provides no halogen.

`scan` flags these for SCRUTINY — it does not delete or re-assign anything. It is
deliberately conservative and only inspects CANDIDATE-tier neutrals: an Assigned
assignment cleared the server's isotope-pattern score, which is independent
corroboration we don't second-guess from element ratios. A neutral seen as
Assigned in ANY channel is therefore never flagged.

Pure formula arithmetic; deterministic. Thresholds are intentionally loose so the
flagged set is small and defensible (the clear coincidences), not a dragnet.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import chemistry as C
from . import ledger as L

__version__ = "0.2.0"   # + Stage 3 demote/relabel hardening (O-monster, carbon-
                        # cluster, adduct-less fragment, series coherence)

# thresholds (loose on purpose — flag the clear coincidences only)
N_HIGH_OC = 3       # N>=3 combined with...
OC_HIGH = 1.0       # ...O/C >= this  -> high-heteroatom coincidence
N_VERY_HIGH = 4     # N>=4 combined with...
O_HIGH = 8          # ...O>= this
HC_FLOOR = 0.35     # H/C below this -> implausibly carbon-rich (F-FREE formulas only)
F_HIGH = 4          # F>=this -> heavily fluorinated; F is monoisotopic, so the fit
                    # has NO isotope twin to confirm it (a fluorine mass coincidence).
                    # Flagged with the RIGHT reason: low H/C here is F displacing H
                    # (fluorine-rich), NOT a carbon-rich skeleton.

# --- the SHARED plausibility oracle (Stage 3) ------------------------------
# The demote/relabel steps (cleanup-level + this module's scrutiny scan) all read
# the SAME predicates so a flagged formula and a demoted formula are never out of
# step. Two gates were hardened (live-data calibrated, see assign.py wiring):
#
#   * O-MONSTER: an extreme O/C (> OC_MONSTER) is the oxygen-lattice mass-fit
#     signature. It is NOT a niso gate -- a 13C satellite confirms the CARBON
#     count, not the O count, so an O-monster carrying a real 13C twin is still an
#     O-monster. Real HOMs top out at O/C ~1.14, so OC_MONSTER=1.3 spares every
#     genuine oxidation product. (The DEMOTE additionally requires mass-saturation
#     from the degeneracy audit; the plain reason-string oracle reports the ratio.)
#
#   * CARBON-CLUSTER: DBE/C >= DBE_PER_C_MONSTER (equivalently H <= N+2) on an
#     F-free C>=2 skeleton is a bare-carbon mass coincidence (e.g. C5H2, C24H2 --
#     DBE/C >= 1). H-poor-but-DBE/C<1 skeletons like C27H8 (0.89) are NOT this gate;
#     they are the separate H/C<0.35 carbon-rich demote in cleanup.py.
#     A HALF-INTEGER DBE is EXEMPT: radicals carry half-integer DBE, whereas the
#     carbon-cluster monsters are all integer-DBE. The >= 1.0 cutoff (NOT the
#     earlier 0.75 proposal) is deliberate -- pyridine (0.80), coumarin (0.78),
#     umbelliferone (0.78), furfural (0.80) and phthalic anhydride (0.88) are real
#     aromatics that sit below 1.0 and MUST be spared.
OC_MONSTER = 1.3            # O/C strictly above this -> oxygen-lattice monster
DBE_PER_C_MONSTER = 1.0     # DBE/C at or above this (F-free, C>=2) -> carbon cluster


def _oc(cnt: dict) -> float:
    """O/C ratio (0 when carbon-free; carbon-free is handled elsewhere)."""
    nc = cnt.get("C", 0)
    return cnt.get("O", 0) / nc if nc else 0.0


def is_oxygen_monster(cnt: dict) -> bool:
    """O/C strictly above OC_MONSTER (the oxygen-lattice mass-fit ratio). Pure
    arithmetic -- the DEMOTE additionally gates on degeneracy mass-saturation."""
    return cnt.get("C", 0) > 0 and _oc(cnt) > OC_MONSTER


def is_carbon_cluster(cnt: dict) -> bool:
    """F-free C>=2 skeleton whose DBE/C >= DBE_PER_C_MONSTER, EXCLUDING radicals
    (half-integer DBE are exempt). H <= N+2 is the equivalent integer test, but we
    compute the real DBE so the half-integer radical exemption is exact."""
    nc = cnt.get("C", 0)
    if nc < 2 or cnt.get("F", 0) > 0:
        return False
    d = C.dbe(cnt)
    if abs(d - round(d)) > 1e-9:       # half-integer DBE -> radical, EXEMPT
        return False
    return d / nc >= DBE_PER_C_MONSTER


def implausible(neutral_formula: str, *, tier: str | None = None,
                polarity: str | None = None) -> str | None:
    """Return a short reason string if `neutral_formula` looks like a mass-coincidence
    fit rather than a real molecule, else None. Only Candidate-tier is scrutinised
    (pass tier=None to scrutinise regardless). `polarity` ('+'/'-') enables the
    wrong-mode-halogen check."""
    if tier is not None and str(tier) != "Candidate":
        return None
    c = C.parse_formula(str(neutral_formula))
    nc = c.get("C", 0)
    if nc == 0:
        return None                      # carbon-free handled elsewhere (reagent/inorganic)
    h, n, o = c.get("H", 0), c.get("N", 0), c.get("O", 0)
    f = c.get("F", 0)
    br, cl = c.get("Br", 0), c.get("Cl", 0)
    hc, oc = h / nc, o / nc
    # Terse labels (the full meaning is spelled out in the scrutiny-page legend);
    # keeping them short stops the table overflowing the page width.
    if is_oxygen_monster(c):    # O/C beyond the HOM ceiling -> oxygen-lattice monster
        return f"O/C {oc:.1f} (oxygen-lattice monster)"
    if n >= N_HIGH_OC and oc >= OC_HIGH:
        return f"N{n}, O/C {oc:.1f} (heteroatom coincidence)"
    if n >= N_VERY_HIGH and o >= O_HIGH:
        return f"N{n}O{o} (heteroatom coincidence)"
    if f >= F_HIGH:           # heavily fluorinated: 19F is 100% monoisotopic
        # NB any 13C/81Br satellites the row carries confirm the CARBON count / the
        # adduct halogen, NOT the fluorine -- 19F has no heavier stable isotope, so
        # the F COUNT is never isotope-confirmable (do NOT say "no isotope twin").
        return f"F{f}: 19F monoisotopic, fluorine count not isotope-confirmable"
    if is_carbon_cluster(c):  # DBE/C>=1.0, F-free, integer-DBE (radicals exempt)
        return f"DBE/C {C.dbe(c) / nc:.2f} (carbon cluster, H<=N+2)"
    if f == 0 and hc < HC_FLOOR:     # genuine carbon-rich skeleton (F not displacing H)
        return f"H/C {hc:.2f} (carbon-rich)"
    if polarity == "+" and (br > 0 or cl > 0):
        return "halogen in neutral, +mode"
    return None


def scan(merged, *, polarity: str | None = None) -> list[dict]:
    """Flag Candidate-only neutrals that look implausible. Returns one dict per
    distinct neutral: {neutral_formula, reason, ion_score, tier}. A neutral that is
    Assigned in any ion channel is excluded (it is corroborated)."""
    if merged is None or "neutral_formula" not in getattr(merged, "columns", []):
        return []
    g = merged.dropna(subset=["neutral_formula"]).copy()
    if not len(g):
        return []
    g["neutral_formula"] = g["neutral_formula"].astype(str)
    has_tier = "tier" in g.columns
    out = []
    for f, sub in g.groupby("neutral_formula"):
        best = ("Assigned" if has_tier and (sub["tier"] == "Assigned").any()
                else "Candidate")
        reason = implausible(f, tier=best, polarity=polarity)
        if reason:
            sc = sub["ion_score"].max() if "ion_score" in sub.columns else None
            out.append({"neutral_formula": f, "reason": reason, "tier": best,
                        "ion_score": (float(sc) if sc is not None and sc == sc else None)})
    out.sort(key=lambda d: (d["reason"], d["neutral_formula"]))
    return out


# ===========================================================================
# Stage 3: demote / relabel-ONLY hardening (never deletes a row)
#
# Every function below either DEMOTES an Assigned M0 to Candidate (+ stamps
# below_assignability) or RELABELS an M0's role to fragment. None of them can
# clear/delete a row, so the worst-case failure is an over-cautious Candidate,
# never a lost peak. Each appends one dict per touched peak to `audit` (when a
# list is passed) so assign/assign_batch can write tables/plausibility_audit_*.
# ===========================================================================

def _is_saturated(note) -> bool:
    """A degeneracy_note that flags the mass as saturated/degenerate -- the
    second leg of the O-monster demote (the ratio alone is not enough; the mass
    must also be one arbitrary pick of a degenerate set)."""
    s = "" if note is None or (isinstance(note, float) and pd.isna(note)) else str(note)
    low = s.lower()
    return "satur" in low or "degener" in low


def _iso_count(s) -> int:
    """Number of server-confirmed isotopologues recorded on a row (0 if none)."""
    if isinstance(s, str) and s.strip().startswith("["):
        try:
            return len(json.loads(s))
        except Exception:
            return 0
    return 0


def _m0_index(ledger):
    """M0 rows when the ledger has a role column; otherwise every row (the merged
    ledger is already M0-only and carries no role column)."""
    return (ledger.index[ledger["role"] == L.ROLE_M0]
            if "role" in ledger.columns else ledger.index)


def _append_note(ledger, i, note):
    if "commentary" in ledger.columns:
        cur = ledger.at[i, "commentary"]
        prev = "" if cur is None or (isinstance(cur, float) and pd.isna(cur)) or cur is pd.NA else str(cur)
        ledger.at[i, "commentary"] = (prev + "; " + note) if prev and prev != "nan" else note


def _demote_row(ledger, i, *, reason, audit, evidence, degeneracy_note, n_iso):
    """Demote one M0 -> Candidate + below_assignability, append the note, and log
    one audit record. Demote-only: tier moves Assigned->Candidate, nothing is
    cleared."""
    before = str(ledger.at[i, "tier"]) if "tier" in ledger.columns else ""
    if "tier" in ledger.columns and before == "Assigned":
        ledger.at[i, "tier"] = "Candidate"
    if "below_assignability" in ledger.columns:
        ledger.at[i, "below_assignability"] = True
    _append_note(ledger, i, reason)
    if audit is not None:
        audit.append({
            "mz": ledger.at[i, "mz"] if "mz" in ledger.columns else None,
            "neutral_formula": ledger.at[i, "neutral_formula"],
            "before_tier": before, "after_tier_or_role": "Candidate",
            "reason": reason, "evidence": evidence,
            "degeneracy_note": ("" if degeneracy_note is None
                                or (isinstance(degeneracy_note, float) and pd.isna(degeneracy_note))
                                else str(degeneracy_note)),
            "n_iso": n_iso})


def demote_oxygen_monsters(ledger: pd.DataFrame, *, audit=None, log=print) -> dict:
    """Demote M0 assignments that are oxygen-lattice 'monsters': O/C > OC_MONSTER
    AND mass-saturated (the degeneracy audit flags ~dozens of plausible ions on
    the mass). NOT niso-gated -- a 13C satellite confirms the carbon count, not the
    oxygen count, so it would wrongly exempt a real O-monster. Real HOMs (O/C<=1.14)
    are spared by the ratio cut; non-saturated high-O fits are spared by the second
    leg. Assigned->Candidate + below_assignability. Demote-only."""
    n = 0
    has_note = "degeneracy_note" in ledger.columns
    for i in _m0_index(ledger):
        cnt = C.parse_formula(str(ledger.at[i, "neutral_formula"] or ""))
        if not is_oxygen_monster(cnt):
            continue
        note = ledger.at[i, "degeneracy_note"] if has_note else None
        if not _is_saturated(note):     # ratio alone is not enough -- needs saturation
            continue
        ni = _iso_count(ledger.at[i, "isotopologues"]) if "isotopologues" in ledger.columns else 0
        reason = (f"oxygen-lattice monster (O/C {_oc(cnt):.2f} > {OC_MONSTER}, "
                  "mass-saturated) -- one arbitrary pick of a sub-ppm-degenerate set")
        _demote_row(ledger, i, reason=reason, audit=audit,
                    evidence=f"O/C={_oc(cnt):.2f}", degeneracy_note=note, n_iso=ni)
        n += 1
    log(f"[plausibility] demoted {n} oxygen-lattice monsters (O/C>{OC_MONSTER}, mass-saturated)")
    return {"o_demoted": n}


def demote_carbon_clusters(ledger: pd.DataFrame, *, audit=None, log=print) -> dict:
    """Demote M0 assignments resting on a bare-carbon skeleton: DBE/C >=
    DBE_PER_C_MONSTER (H<=N+2), F-free, C>=2, with the HALF-INTEGER-DBE radical
    EXEMPTION (radicals carry half-integer DBE; carbon-cluster monsters are
    integer-DBE). This is distinct from the H/C<0.35 carbon-rich demote in
    cleanup.py (kept unchanged): the two together cover the carbon-coincidence
    family without catching real aromatics (pyridine/coumarin/furfural sit below
    DBE/C 1.0). Assigned->Candidate + below_assignability. Demote-only."""
    n = 0
    for i in _m0_index(ledger):
        cnt = C.parse_formula(str(ledger.at[i, "neutral_formula"] or ""))
        if not is_carbon_cluster(cnt):
            continue
        nc = cnt.get("C", 0)
        dpc = C.dbe(cnt) / nc
        ni = _iso_count(ledger.at[i, "isotopologues"]) if "isotopologues" in ledger.columns else 0
        reason = (f"carbon cluster (DBE/C {dpc:.2f} >= {DBE_PER_C_MONSTER}, H<=N+2, "
                  "F-free) -- bare-carbon mass coincidence, not a molecule")
        note = ledger.at[i, "degeneracy_note"] if "degeneracy_note" in ledger.columns else None
        _demote_row(ledger, i, reason=reason, audit=audit,
                    evidence=f"DBE/C={dpc:.2f}", degeneracy_note=note, n_iso=ni)
        n += 1
    log(f"[plausibility] demoted {n} carbon clusters (DBE/C>={DBE_PER_C_MONSTER}, "
        "F-free, integer-DBE)")
    return {"c_cluster_demoted": n}


def demote_implausible(ledger: pd.DataFrame, *, audit=None, log=print) -> dict:
    """The two shared-oracle demotes that fire on a single-file or merged ledger
    without a time series: O-monster + carbon-cluster. Both are demote-only and
    feed the same audit list. (The adduct-less-fragment relabel and the
    series-coherence dissolve need the spectrum / time series, so they run
    separately at the batch level.)"""
    o = demote_oxygen_monsters(ledger, audit=audit, log=log)
    c = demote_carbon_clusters(ledger, audit=audit, log=log)
    return {**o, **c}


# --- adduct-less in-source fragment -> role=fragment (batch/merged level) ----
# A real molecule ionises through its adduct channels (in a urea source: +NH4,
# +urea·H, +Na on top of the bare +H). An in-source FRAGMENT is generated in the
# inlet from a heavier parent, so it shows up almost entirely as the bare protonated
# ion with NO adduct partners. The relabel needs the FULL triangulation:
#   (a) adduct ratio (Σ adduct channels)/[M+H]+ < FRAG_ADDUCT_RATIO with a big bare
#       [M+H]+, AND
#   (b) a mass-consistent CO-VARYING parent: a heavier assigned neutral at a facile
#       neutral loss whose time-series trace beats the incidental co-rider ceiling.
# (a) alone only earns a scrutiny COMMENTARY flag (never a relabel) -- a low-adduct
# species could simply be a poor adduct-former. Composition gate: hydrocarbon / low-O
# (fragments are the de-oxygenated/-dehydrated cores, not fresh O-rich oxidation).
FRAG_ADDUCT_RATIO = 0.05         # Σ(adduct)/[M+H]+ below this -> adduct-less
FRAG_BARE_MIN = 0.0              # the bare [M+H]+ must be present (height>this)
FRAG_PARENT_R_MIN = 0.5          # parent TS correlation must beat this co-rider floor
FRAG_OC_MAX = 0.5                # composition gate: fragment cores are low-O
# facile neutral losses an in-source fragment forms from its parent (formula, label)
_FRAG_LOSSES = [
    ("H2O", "H2O"),
    ("CO", "CO"),
    ("CO2", "CO2"),
    ("CH2O3", "CO+H2O"),   # CO + H2O combined (formic-equivalent)
]
_POS_ADDUCTS = ("[M+NH4]+", "[M+Na]+", "[M+K]+", "[M+(CH4N2O)H]+")


def _facile_loss(parent_cnt: dict, child_cnt: dict):
    """If child = parent - (a facile neutral loss), return the loss label, else None.
    Element-exact: every element of the loss must be removable from the parent."""
    child = {k: v for k, v in child_cnt.items() if v}
    for loss_formula, label in _FRAG_LOSSES:
        cand = dict(parent_cnt)
        ok = True
        for el, k in C.parse_formula(loss_formula).items():
            cand[el] = cand.get(el, 0) - k
            if cand[el] < 0:
                ok = False
                break
        if ok and {k: v for k, v in cand.items() if v} == child:
            return label
    return None


def relabel_adduct_less_fragments(merged: pd.DataFrame, *, ts_peaks=None,
                                  polarity: str = "+", audit=None,
                                  r_min: float = FRAG_PARENT_R_MIN, log=print) -> dict:
    """Relabel an adduct-less protonated M0 as an in-source FRAGMENT of a heavier
    co-varying parent -- ONLY on the full triangulation (adduct-ratio AND a
    mass-consistent co-varying parent). When only the adduct-ratio leg holds, add a
    scrutiny commentary flag but DO NOT relabel. Positive mode only (the adduct
    grammar is the +NH4/+urea/+Na family). Relabel-only: role M0->fragment, the
    neutral/mz/score are unchanged."""
    if polarity != "+" or not {"neutral_formula", "adduct"} <= set(merged.columns):
        return {"relabeled": 0, "flagged": 0}
    if "commentary" not in merged.columns:    # both legs annotate commentary
        merged["commentary"] = pd.NA
    # adduct inventory per neutral: bare [M+H]+ height vs the Σ adduct-channel height
    h = (merged["height"] if "height" in merged.columns
         else merged.get("ion_score", pd.Series(1.0, index=merged.index)))
    bare, addsum, present = {}, {}, {}
    for i in merged.index:
        nf = str(merged.at[i, "neutral_formula"] or "")
        if not nf or nf == "nan":
            continue
        ad = str(merged.at[i, "adduct"] or "")
        hv = float(h.at[i]) if pd.notna(h.at[i]) else 0.0
        present.setdefault(nf, {})[ad] = i
        if ad == "[M+H]+":
            bare[nf] = max(bare.get(nf, 0.0), hv)
        elif ad in _POS_ADDUCTS:
            addsum[nf] = addsum.get(nf, 0.0) + hv

    keeps_parent = _make_parent_test(merged, ts_peaks, r_min)
    relabeled = flagged = 0
    for nf, chans in present.items():
        if "[M+H]+" not in chans:
            continue
        b = bare.get(nf, 0.0)
        if b <= FRAG_BARE_MIN:
            continue
        ratio = addsum.get(nf, 0.0) / b if b else 9.99
        if ratio >= FRAG_ADDUCT_RATIO:
            continue                                  # forms adducts -> a real molecule
        cnt = C.parse_formula(nf)
        if cnt.get("C", 0) < 1 or _oc(cnt) > FRAG_OC_MAX:
            continue                                  # composition gate: low-O cores only
        i = chans["[M+H]+"]
        # leg (b): a heavier assigned co-varying parent at a facile loss
        parent = _find_covarying_parent(merged, nf, cnt, keeps_parent)
        if parent is None:
            # only leg (a): scrutiny flag, NO relabel
            _append_note(merged, i,
                         f"adduct-less (Σadduct/[M+H]+ {ratio:.2f} < {FRAG_ADDUCT_RATIO}) "
                         "-- scrutinise as a possible in-source fragment")
            flagged += 1
            continue
        pf, loss = parent
        note = f"in-source fragment of {pf} ({loss})"
        try:
            # a real per-file ledger (peak_id + role + locked) goes through the
            # invariant-enforcing API; the merged align() frame (no locked column)
            # is relabelled by index.
            if {"peak_id", "role", "locked"} <= set(merged.columns):
                L.mark_fragment(merged, merged.at[i, "peak_id"], note)
            else:
                _mark_fragment_byidx(merged, i, note)
        except Exception as e:
            log(f"[plausibility] fragment relabel skipped for {nf}: {e}")
            continue
        if audit is not None:
            audit.append({
                "mz": merged.at[i, "mz"] if "mz" in merged.columns else None,
                "neutral_formula": nf, "before_tier": str(merged.at[i, "tier"])
                if "tier" in merged.columns else "",
                "after_tier_or_role": "fragment", "reason": note,
                "evidence": f"Σadduct/[M+H]+={ratio:.2f}", "degeneracy_note": "",
                "n_iso": _iso_count(merged.at[i, "isotopologues"])
                if "isotopologues" in merged.columns else 0})
        relabeled += 1
    log(f"[plausibility] adduct-less fragments: {relabeled} relabeled (full "
        f"triangulation), {flagged} flagged (adduct-ratio only)")
    return {"relabeled": relabeled, "flagged": flagged}


def _mark_fragment_byidx(merged, i, note):
    """Role->fragment on a merged frame that has no peak_id / role column: set the
    role column (creating it if absent, defaulting M0) so report.py's role filter
    and the VK/cluster drivers see the fragment. commentary is ensured by the
    caller."""
    if "role" not in merged.columns:
        merged["role"] = L.ROLE_M0
    merged.at[i, "role"] = L.ROLE_FRAGMENT
    merged.at[i, "commentary"] = note


def _find_covarying_parent(merged, child_nf, child_cnt, keeps_parent):
    """Return (parent_formula, loss_label) for a heavier ASSIGNED neutral that is a
    facile-loss parent of `child_nf` AND co-varies above the co-rider ceiling, else
    None. Deterministic: scans candidate parents in sorted order, takes the first
    co-varying facile-loss match."""
    for pf in sorted({str(x) for x in merged["neutral_formula"].dropna()}):
        if pf == child_nf:
            continue
        pcnt = C.parse_formula(pf)
        if C.neutral_mass(pcnt) <= C.neutral_mass(child_cnt):
            continue
        loss = _facile_loss(pcnt, child_cnt)
        if loss is None:
            continue
        if keeps_parent(child_nf, pf):
            return pf, loss
    return None


def _make_parent_test(merged, ts_peaks, r_min):
    """Return covaries(child_nf, parent_nf) -> bool. The FULL triangulation's leg (b)
    needs CO-VARIATION evidence, so WITHOUT a time series no parent ever qualifies
    (the relabel is suppressed and the species only earns the scrutiny flag) --
    deliberately conservative: a structural facile-loss link alone is not proof of
    in-source fragmentation. WITH a batch time series it is the log1p-trace
    correlation >= r_min between the child [M+H]+ and the parent's brightest
    channel, which must beat the incidental co-rider ceiling."""
    if ts_peaks is None or not len(getattr(ts_peaks, "index", [])):
        return lambda child, parent: False   # no TS -> no co-variation -> no relabel

    from . import timeseries as TS
    mat, bin_mz = TS.build_matrix(ts_peaks)
    if not len(mat):
        return lambda child, parent: False
    bm = bin_mz.sort_values()
    arr = bm.to_numpy(); idx = bm.index.to_numpy()

    def _logtrace(nf, adduct):
        try:
            mz = C.ion_mz(nf, adduct)
        except Exception:
            return None
        j = np.searchsorted(arr, mz); best = None
        for k in (j - 1, j):
            if 0 <= k < len(arr):
                ppm = abs(arr[k] - mz) / mz * 1e6
                if ppm <= 8 and (best is None or ppm < best[1]):
                    best = (idx[k], ppm)
        if best is None or best[0] not in mat.columns:
            return None
        return np.log1p(mat[best[0]].fillna(0.0).clip(lower=0).to_numpy())

    def covaries(child, parent):
        ct = _logtrace(child, "[M+H]+")
        if ct is None:
            return False
        for ad in ("[M+H]+", "[M+NH4]+", "[M+(CH4N2O)H]+", "[M+Na]+"):
            pt = _logtrace(parent, ad)
            if pt is None:
                continue
            ok = np.isfinite(ct) & np.isfinite(pt)
            if ok.sum() >= 6 and np.std(ct[ok]) > 0 and np.std(pt[ok]) > 0 \
                    and np.corrcoef(ct[ok], pt[ok])[0, 1] >= r_min:
                return True
        return False
    return covaries


# --- series coherence: dissolve mutually-uncorrelated detected series ---------
# A real homolog / dehydrogenation series co-elutes (its members share a source),
# so its members' time traces are mutually correlated. A SPURIOUS series (the
# C5H2/C5H4 -H2 mass ladder the residual explainer chains through noise) has members
# whose traces are internally uncorrelated (r ~ 0-0.25). When the median pairwise
# log1p-correlation of a series falls below SERIES_R_MIN, the series is incoherent:
# un-commit its members (demote to Candidate + below_assignability so they are never
# silently lost). A co-varying real series is never touched.
SERIES_R_MIN = 0.5


def dissolve_incoherent_series(merged: pd.DataFrame, *, ts_peaks=None, audit=None,
                               r_min: float = SERIES_R_MIN, log=print) -> dict:
    """Demote the members of a detected series whose time traces are mutually
    UNCORRELATED (median pairwise log1p-correlation < r_min). Needs a series
    membership column (series_unit) AND a time series; inert without either.
    Demote-only (Candidate + below_assignability), never deletes a member."""
    if ts_peaks is None or "series_unit" not in merged.columns \
            or "neutral_formula" not in merged.columns:
        return {"series_dissolved": 0, "members_demoted": 0}
    from . import timeseries as TS
    mat, bin_mz = TS.build_matrix(ts_peaks)
    if not len(mat):
        return {"series_dissolved": 0, "members_demoted": 0}
    bm = bin_mz.sort_values()
    arr = bm.to_numpy(); idx = bm.index.to_numpy()

    def _logtrace(nf, ad):
        try:
            mz = C.ion_mz(str(nf), str(ad))
        except Exception:
            return None
        j = np.searchsorted(arr, mz); best = None
        for k in (j - 1, j):
            if 0 <= k < len(arr):
                ppm = abs(arr[k] - mz) / mz * 1e6
                if ppm <= 8 and (best is None or ppm < best[1]):
                    best = (idx[k], ppm)
        if best is None or best[0] not in mat.columns:
            return None
        return np.log1p(mat[best[0]].fillna(0.0).clip(lower=0).to_numpy())

    has_role = "role" in merged.columns
    series_n = members_n = 0
    for unit, g in merged.dropna(subset=["series_unit"]).groupby("series_unit"):
        if has_role:
            g = g[g["role"] == L.ROLE_M0]
        if len(g) < 3:                       # need >=3 members for a coherence read
            continue
        traces = []
        for _, r in g.iterrows():
            t = _logtrace(r["neutral_formula"], r.get("adduct", "[M+H]+"))
            if t is not None and np.std(t) > 0:
                traces.append(t)
        if len(traces) < 3:
            continue
        rs = []
        for a in range(len(traces)):
            for b in range(a + 1, len(traces)):
                ok = np.isfinite(traces[a]) & np.isfinite(traces[b])
                if ok.sum() >= 6:
                    rs.append(np.corrcoef(traces[a][ok], traces[b][ok])[0, 1])
        if not rs or float(np.median(rs)) >= r_min:
            continue                          # coherent (or unmeasurable) -> spare it
        med = float(np.median(rs))
        for i in g.index:
            reason = (f"incoherent series {unit} (median pairwise r {med:.2f} < "
                      f"{r_min}) -- members do not co-vary, likely a spurious mass ladder")
            ni = _iso_count(merged.at[i, "isotopologues"]) if "isotopologues" in merged.columns else 0
            _demote_row(merged, i, reason=reason, audit=audit,
                        evidence=f"series_r={med:.2f}", degeneracy_note=None, n_iso=ni)
            members_n += 1
        series_n += 1
    log(f"[plausibility] dissolved {series_n} incoherent series ({members_n} members demoted)")
    return {"series_dissolved": series_n, "members_demoted": members_n}


_AUDIT_COLS = ["mz", "neutral_formula", "before_tier", "after_tier_or_role",
               "reason", "evidence", "degeneracy_note", "n_iso"]


def write_audit(audit: list, path: str) -> int:
    """Write the plausibility audit (one row per touched peak) to `path`, sorted
    deterministically by mz then formula. Always writes the header (an empty audit
    still produces a 1-line CSV so the artifact set is stable). Returns the row
    count."""
    df = pd.DataFrame(audit, columns=_AUDIT_COLS)
    if len(df):
        df = df.sort_values(["mz", "neutral_formula"], na_position="last").reset_index(drop=True)
    df.to_csv(path, index=False)
    return len(df)

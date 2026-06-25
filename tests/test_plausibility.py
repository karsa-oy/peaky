"""Offline tests for plausibility.py (chemical-plausibility QC of assignments).
Run: python3 tests/test_plausibility.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky import plausibility as PL  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


# --- implausible() rules ---
check("N-monster high O/C flagged (Candidate)",
      PL.implausible("C9H12N4O12", tier="Candidate") is not None)
check("N5O20 monster flagged", PL.implausible("C15H37N5O20", tier="Candidate") is not None)
check("very low H/C flagged", PL.implausible("C35H4", tier="Candidate") is not None)
check("positive-mode halogen flagged",
      PL.implausible("C10H18Br2O12", tier="Candidate", polarity="+") is not None)
check("same halogen NOT flagged in negative mode",
      PL.implausible("C10H18Br2O12", tier="Candidate", polarity="-") is None
      or "halogen" not in (PL.implausible("C10H18Br2O12", tier="Candidate", polarity="-") or ""))

# real molecules are NOT flagged
check("ordinary CHO not flagged", PL.implausible("C10H16O2", tier="Candidate") is None)
check("monoterpene-oxidation product not flagged", PL.implausible("C10H16O5", tier="Candidate") is None)
check("small N1 amine not flagged", PL.implausible("C10H19NO2", tier="Candidate") is None)
check("modest N3 low-O species not flagged (O/C<1)",
      PL.implausible("C12H21N3O3", tier="Candidate") is None)
check("HOM dimer not flagged", PL.implausible("C20H30O14", tier="Candidate") is None)

# Assigned tier is never second-guessed
check("Assigned N-monster NOT flagged (corroborated by isotope score)",
      PL.implausible("C9H12N4O12", tier="Assigned") is None)

# --- scan() de-duplicates and respects best-tier-per-neutral ---
merged = pd.DataFrame([
    dict(neutral_formula="C9H12N4O12", adduct="[M+Na]+", tier="Candidate", ion_score=0.94),
    dict(neutral_formula="C10H16O2", adduct="[M+H]+", tier="Assigned", ion_score=0.9),
    # same monster appears Assigned in one channel -> must NOT be flagged
    dict(neutral_formula="C8H13N3O10", adduct="[M+H]+", tier="Candidate", ion_score=0.6),
    dict(neutral_formula="C8H13N3O10", adduct="[M+NH4]+", tier="Assigned", ion_score=0.8),
    dict(neutral_formula="C10H18Br2O12", adduct="[M+H]+", tier="Candidate", ion_score=0.65),
])
flagged = PL.scan(merged, polarity="+")
names = {d["neutral_formula"] for d in flagged}
check("scan flags the Candidate-only N-monster", "C9H12N4O12" in names, names)
check("scan flags the positive-mode dibromo", "C10H18Br2O12" in names, names)
check("scan excludes a neutral that is Assigned in any channel",
      "C8H13N3O10" not in names, names)
check("scan excludes ordinary CHO", "C10H16O2" not in names, names)
check("scan carries reason + score", all("reason" in d and "ion_score" in d for d in flagged))

# negative mode: the dibromo is NOT flagged for halogen (covalent organohalogen is plausible)
flagged_neg = PL.scan(merged, polarity="-")
check("negative mode does not flag the dibromo on halogen grounds",
      "C10H18Br2O12" not in {d["neutral_formula"] for d in flagged_neg},
      {d["neutral_formula"] for d in flagged_neg})

# empty / missing columns degrade gracefully
check("scan empty frame -> []", PL.scan(pd.DataFrame()) == [])


# ===========================================================================
# Stage 3: HARDENED demote/relabel gates (the shared oracle + the demotes)
# ===========================================================================
import numpy as np                      # noqa: E402
from peaky import chemistry as C        # noqa: E402
from peaky import ledger as L           # noqa: E402

# --- shared oracle: is_oxygen_monster / is_carbon_cluster ---
def cf(s): return C.parse_formula(s)

check("oracle: O-monster O/C 1.6 flagged (C5H4O8)", PL.is_oxygen_monster(cf("C5H4O8")))
check("oracle: O-monster lattice fit (C3H5ClO17)", PL.is_oxygen_monster(cf("C3H5ClO17")))
check("oracle: real HOM O/C 0.7 NOT an O-monster (C10H16O7)", not PL.is_oxygen_monster(cf("C10H16O7")))
check("oracle: HOM dimer O/C 0.7 NOT an O-monster (C20H30O14)", not PL.is_oxygen_monster(cf("C20H30O14")))
check("oracle: O/C exactly 1.3 NOT a monster (strict >)", not PL.is_oxygen_monster(cf("C10H10O13")))

check("oracle: carbon cluster DBE/C 1.0 flagged (C24H2)", PL.is_carbon_cluster(cf("C24H2")))
# real aromatics sit below DBE/C 1.0 and MUST be spared (the 0.75 cutoff caught them)
check("oracle: pyridine C5H5N (DBE/C 0.80) NOT a carbon cluster", not PL.is_carbon_cluster(cf("C5H5N")))
check("oracle: coumarin C9H6O2 (DBE/C 0.78) NOT a carbon cluster", not PL.is_carbon_cluster(cf("C9H6O2")))
check("oracle: furfural C5H4O2 (DBE/C 0.80) NOT a carbon cluster", not PL.is_carbon_cluster(cf("C5H4O2")))
check("oracle: umbelliferone C9H6O3 NOT a carbon cluster", not PL.is_carbon_cluster(cf("C9H6O3")))
check("oracle: phthalic anhydride C8H4O3 (0.88) NOT a carbon cluster", not PL.is_carbon_cluster(cf("C8H4O3")))
# HALF-INTEGER DBE (radical) is EXEMPT even at high DBE/C
check("oracle: half-integer-DBE radical C10H3 (DBE 9.5) EXEMPT",
      abs(C.dbe(cf("C10H3")) - round(C.dbe(cf("C10H3")))) > 1e-9 and not PL.is_carbon_cluster(cf("C10H3")))
check("oracle: F-rich low-H/C NOT a carbon cluster (F-free rule)", not PL.is_carbon_cluster(cf("C11H6F16")))
check("oracle: single carbon C1 not a cluster (C>=2)", not PL.is_carbon_cluster(cf("CH2")))

# implausible() shares the same oracle (one source of truth)
check("implausible: O-monster reason via shared oracle",
      "oxygen-lattice monster" in (PL.implausible("C5H4O8", tier="Candidate") or ""))
check("implausible: carbon-cluster reason via shared oracle",
      "carbon cluster" in (PL.implausible("C24H2", tier="Candidate") or ""))
check("implausible: pyridine NOT flagged by either gate", PL.implausible("C5H5N", tier="Candidate") is None)
check("implausible: coumarin NOT flagged", PL.implausible("C9H6O2", tier="Candidate") is None)
check("implausible: furfural NOT flagged", PL.implausible("C5H4O2", tier="Candidate") is None)

# --- demote_oxygen_monsters: O/C>1.3 AND mass-saturated (NOT niso-gated) ---
ledo = pd.DataFrame([
    dict(role="M0", mz=300.0, neutral_formula="C5H4O8", tier="Assigned", commentary="",
         below_assignability=False, degeneracy_note="MASS-SATURATED: 27 plausible formulas", isotopologues="[]"),
    # O-monster carrying a real 13C twin -> STILL demoted (niso must NOT exempt it)
    dict(role="M0", mz=305.0, neutral_formula="C6H4O9", tier="Assigned", commentary="",
         below_assignability=False, degeneracy_note="MASS-SATURATED: 19 plausible formulas",
         isotopologues='[{"label": "13C", "score": 0.9}]'),
    # high O/C but NOT saturated -> spared (the second leg)
    dict(role="M0", mz=310.0, neutral_formula="C4H4O7", tier="Assigned", commentary="",
         below_assignability=False, degeneracy_note="unique within 3 sigma window", isotopologues="[]"),
    # real HOM O/C 0.7, even if saturated -> spared (the ratio leg)
    dict(role="M0", mz=320.0, neutral_formula="C10H16O7", tier="Assigned", commentary="",
         below_assignability=False, degeneracy_note="MASS-SATURATED: 14 plausible formulas", isotopologues="[]"),
])
audit_o = []
outo = PL.demote_oxygen_monsters(ledo, audit=audit_o, log=lambda *a: None)
check("O-monster demote: 2 saturated O-monsters demoted (incl. one with a 13C twin)",
      outo == {"o_demoted": 2}, outo)
check("O-monster demote: C5H4O8 -> Candidate + below_assignability",
      ledo.loc[0, "tier"] == "Candidate" and bool(ledo.loc[0, "below_assignability"]))
check("O-monster demote: NOT niso-gated (C6H4O9 with 13C still demoted)",
      ledo.loc[1, "tier"] == "Candidate" and bool(ledo.loc[1, "below_assignability"]))
check("O-monster demote: high-O but NOT saturated spared (C4H4O7)",
      ledo.loc[2, "tier"] == "Assigned" and not bool(ledo.loc[2, "below_assignability"]))
check("O-monster demote: real HOM C10H16O7 spared (ratio leg)",
      ledo.loc[3, "tier"] == "Assigned" and not bool(ledo.loc[3, "below_assignability"]))
check("O-monster demote: audit one row per touched peak", len(audit_o) == 2, audit_o)
check("O-monster demote: audit carries O/C evidence + degeneracy note",
      all("O/C" in a["evidence"] and "SATUR" in a["degeneracy_note"].upper() for a in audit_o))

# --- demote_carbon_clusters: DBE/C>=1.0, F-free, radical-exempt ---
ledc = pd.DataFrame([
    dict(role="M0", mz=290.0, neutral_formula="C24H2", tier="Assigned", commentary="",
         below_assignability=False, isotopologues="[]"),                    # DBE/C 1.0 -> demote
    dict(role="M0", mz=300.0, neutral_formula="C5H5N",  tier="Assigned", commentary="",
         below_assignability=False, isotopologues="[]"),                    # pyridine -> spare
    dict(role="M0", mz=310.0, neutral_formula="C9H6O2", tier="Assigned", commentary="",
         below_assignability=False, isotopologues="[]"),                    # coumarin -> spare
    dict(role="M0", mz=320.0, neutral_formula="C5H4O2", tier="Assigned", commentary="",
         below_assignability=False, isotopologues="[]"),                    # furfural -> spare
    dict(role="M0", mz=330.0, neutral_formula="C10H3",  tier="Assigned", commentary="",
         below_assignability=False, isotopologues="[]"),                    # half-int DBE radical -> EXEMPT
    dict(role="M0", mz=340.0, neutral_formula="C10H16O4", tier="Assigned", commentary="",
         below_assignability=False, isotopologues="[]"),                    # ordinary SOA -> spare
])
audit_c = []
outc = PL.demote_carbon_clusters(ledc, audit=audit_c, log=lambda *a: None)
check("carbon-cluster demote: only the DBE/C>=1.0 integer cluster (C24H2)", outc == {"c_cluster_demoted": 1}, outc)
check("carbon-cluster demote: C24H2 -> Candidate + below_assignability",
      ledc.loc[0, "tier"] == "Candidate" and bool(ledc.loc[0, "below_assignability"]))
check("carbon-cluster demote: pyridine C5H5N spared", ledc.loc[1, "tier"] == "Assigned")
check("carbon-cluster demote: coumarin C9H6O2 spared", ledc.loc[2, "tier"] == "Assigned")
check("carbon-cluster demote: furfural C5H4O2 spared", ledc.loc[3, "tier"] == "Assigned")
check("carbon-cluster demote: half-integer-DBE radical C10H3 EXEMPT",
      ledc.loc[4, "tier"] == "Assigned" and not bool(ledc.loc[4, "below_assignability"]))
check("carbon-cluster demote: ordinary SOA C10H16O4 spared", ledc.loc[5, "tier"] == "Assigned")

# demote_implausible runs both, and NEVER deletes a row (row count is preserved)
ledboth = pd.concat([ledo, ledc], ignore_index=True)
ledboth["tier"] = "Assigned"; ledboth["below_assignability"] = False
n_before = len(ledboth)
PL.demote_implausible(ledboth, audit=[], log=lambda *a: None)
check("demote_implausible: demote-only, never deletes a row", len(ledboth) == n_before)
check("demote_implausible: no row left without a tier value",
      ledboth["tier"].isin(["Assigned", "Candidate"]).all())

# --- adduct-less fragment relabel (full triangulation needs TS) ---
_parent, _child = "C6H12O3", "C5H12O2"     # child = parent - CO (a facile loss)
check("fragment: facile-loss CO detected",
      PL._facile_loss(C.parse_formula(_parent), C.parse_formula(_child)) == "CO")
_mzc, _mzp = C.ion_mz(_child, "[M+H]+"), C.ion_mz(_parent, "[M+H]+")
_t = np.linspace(1, 100, 30)
_rows = []
for _s, _sc in enumerate(_t):
    _rows.append(dict(sample_item_id=f"s{_s}", mz=_mzc, height=_sc * 10))
    _rows.append(dict(sample_item_id=f"s{_s}", mz=_mzp, height=_sc * 8 + 0.01 * _s))
_ts = pd.DataFrame(_rows)
mfrag = pd.DataFrame([
    dict(peak_id="pc", mz=_mzc, neutral_formula=_child, adduct="[M+H]+", tier="Assigned", ion_score=0.9, height=5000.0, role="M0"),
    dict(peak_id="pp", mz=_mzp, neutral_formula=_parent, adduct="[M+H]+", tier="Assigned", ion_score=0.9, height=8000.0, role="M0"),
    dict(peak_id="pn", mz=C.ion_mz(_parent, "[M+NH4]+"), neutral_formula=_parent, adduct="[M+NH4]+", tier="Assigned", ion_score=0.8, height=2000.0, role="M0"),
])
afr = []
of = PL.relabel_adduct_less_fragments(mfrag, ts_peaks=_ts, polarity="+", audit=afr, log=lambda *a: None)
check("fragment: adduct-less child + co-varying parent -> relabeled (full triangulation)",
      of["relabeled"] == 1, of)
check("fragment: child role -> fragment", mfrag.loc[0, "role"] == L.ROLE_FRAGMENT)
check("fragment: commentary names the parent + loss",
      mfrag.loc[0, "commentary"] == f"in-source fragment of {_parent} (CO)", mfrag.loc[0, "commentary"])
check("fragment: parent (forms adducts) untouched", mfrag.loc[1, "role"] == "M0")
check("fragment: audit records the relabel as a role change",
      len(afr) == 1 and afr[0]["after_tier_or_role"] == "fragment")

# only leg (a): adduct-less but NO time series -> scrutiny flag, NO relabel
mfrag2 = pd.DataFrame([
    dict(mz=_mzc, neutral_formula=_child, adduct="[M+H]+", tier="Assigned", ion_score=0.9, height=5000.0),
    dict(mz=_mzp, neutral_formula=_parent, adduct="[M+H]+", tier="Assigned", ion_score=0.9, height=8000.0),
    dict(mz=C.ion_mz(_parent, "[M+NH4]+"), neutral_formula=_parent, adduct="[M+NH4]+", tier="Assigned", ion_score=0.8, height=2000.0),
])
of2 = PL.relabel_adduct_less_fragments(mfrag2, ts_peaks=None, polarity="+", audit=[], log=lambda *a: None)
check("fragment: leg-(a)-only (no TS) -> flag, never relabel", of2 == {"relabeled": 0, "flagged": 1}, of2)
check("fragment: scrutiny flag written to commentary",
      "scrutinise" in str(mfrag2.loc[0, "commentary"]))
check("fragment: negative mode -> no-op",
      PL.relabel_adduct_less_fragments(mfrag2.copy(), ts_peaks=None, polarity="-", audit=[], log=lambda *a: None)
      == {"relabeled": 0, "flagged": 0})

# --- series coherence: dissolve mutually-uncorrelated series only ---
_forms = ["C5H2", "C5H4", "C5H6"]
_rng = np.random.RandomState(0)
def _series_ts(corr):
    base = _rng.rand(40) * 100
    rows = []
    for fm in _forms:
        mz = C.ion_mz(fm, "[M+H]+")
        trace = (base * (1 + 0.01 * _rng.rand(40)) if corr else _rng.rand(40) * 100) + 1
        for s in range(40):
            rows.append(dict(sample_item_id=f"s{s}", mz=mz, height=trace[s]))
    return pd.DataFrame(rows)
def _series_merged():
    return pd.DataFrame([dict(mz=C.ion_mz(fm, "[M+H]+"), neutral_formula=fm, adduct="[M+H]+",
                              tier="Assigned", ion_score=0.8, role="M0", series_unit="-H2 ladder",
                              below_assignability=False, commentary="", isotopologues="[]") for fm in _forms])
ms_inc = _series_merged()
audit_s = []
os_inc = PL.dissolve_incoherent_series(ms_inc, ts_peaks=_series_ts(False), audit=audit_s, log=lambda *a: None)
check("series: mutually-uncorrelated series dissolved", os_inc == {"series_dissolved": 1, "members_demoted": 3}, os_inc)
check("series: all members -> Candidate + below_assignability",
      (ms_inc["tier"] == "Candidate").all() and ms_inc["below_assignability"].all())
check("series: audit one row per demoted member", len(audit_s) == 3)
ms_coh = _series_merged()
os_coh = PL.dissolve_incoherent_series(ms_coh, ts_peaks=_series_ts(True), audit=[], log=lambda *a: None)
check("series: co-varying real series spared", os_coh["series_dissolved"] == 0 and (ms_coh["tier"] == "Assigned").all())
check("series: no series_unit column -> inert",
      PL.dissolve_incoherent_series(ms_inc.drop(columns=["series_unit"]), ts_peaks=_series_ts(False),
                                    audit=[], log=lambda *a: None) == {"series_dissolved": 0, "members_demoted": 0})
check("series: no time series -> inert",
      PL.dissolve_incoherent_series(_series_merged(), ts_peaks=None, audit=[], log=lambda *a: None)
      == {"series_dissolved": 0, "members_demoted": 0})

# --- write_audit: deterministic + always a header ---
import tempfile, os as _os    # noqa: E402
with tempfile.TemporaryDirectory() as _d:
    _n0 = PL.write_audit([], _os.path.join(_d, "empty.csv"))
    _hdr = open(_os.path.join(_d, "empty.csv")).readline().strip()
    check("audit: empty -> 0 rows but a header line", _n0 == 0 and _hdr.startswith("mz,neutral_formula"))
    _rows_a = [dict(mz=300.0, neutral_formula="C5H4O8", before_tier="Assigned",
                    after_tier_or_role="Candidate", reason="x", evidence="O/C=1.6",
                    degeneracy_note="SAT", n_iso=0),
               dict(mz=200.0, neutral_formula="C24H2", before_tier="Assigned",
                    after_tier_or_role="Candidate", reason="y", evidence="DBE/C=1.0",
                    degeneracy_note="", n_iso=0)]
    _n1 = PL.write_audit(_rows_a, _os.path.join(_d, "a.csv"))
    _df = pd.read_csv(_os.path.join(_d, "a.csv"))
    check("audit: nonempty rows written + sorted by mz", _n1 == 2 and list(_df["mz"]) == [200.0, 300.0])


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

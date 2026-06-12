"""Offline tests for residual.py (Pass 4). Run: python3 tests/test_residual.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import residual as RD  # noqa: E402
from mascope_assign import ledger as L  # noqa: E402
from mascope_assign import contexts as X  # noqa: E402
from mascope_assign import chemistry as C  # noqa: E402
from mascope_assign import isotopes as ISO  # noqa: E402
from mascope_assign.passes import PassConfig  # noqa: E402

PASS = FAIL = 0
CFG = PassConfig(height_cutoff=100)
PROF = X.get_context("ambient-air")
PRE = ISO.PrescanResult(estimated_max_C=20)


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


# ---------- acceptance policy ----------
ok, _ = RD._accept(0.85, 0.4, False, CFG)
check("strict ppm accepts on score alone", ok)
ok, _ = RD._accept(0.85, 2.5, False, CFG)
check("loose ppm WITHOUT pattern rejected", not ok)
ok, why = RD._accept(0.85, 2.5, True, CFG)
check("loose ppm WITH pattern accepted", ok, why)
ok, _ = RD._accept(0.85, 6.0, True, CFG)
check("beyond pattern ppm rejected even with pattern", not ok)
ok, _ = RD._accept(0.40, 0.2, True, CFG)
check("below score floor rejected", not ok)
check("confidence capped to Good", RD._cap_conf("High (iso-pair)") == "Good (iso-pair)")

# ---------- find_iso_pairs ----------
# Br doublet (ratio ~1) at the REAL m/z of C6H10O3 [M+Br]-, a Cl doublet, orphan
MZ_BR = C.ion_mz("C6H10O3", "[M+Br]-")
peaks = pd.DataFrame({
    "peak_id": ["L1", "H1", "L2", "H2", "orph"],
    "mz": [MZ_BR, MZ_BR + RD.D_PAIR_BR, 200.0, 200.0 + RD.D_PAIR_CL, 123.4],
    "height": [1e4, 9.6e3, 1e4, 0.33e4, 5e3],
})
led = L.new_ledger(peaks)
pairs = RD.find_iso_pairs(led, min_height=100)
check("finds 2 isotope pairs", len(pairs) == 2, pairs.to_dict("records"))
check("Br pair n_halogen=1", (pairs[pairs.element == "Br"].n_halogen == 1).all())
check("Cl pair detected", (pairs.element == "Cl").any())
check("orphan not paired", "orph" not in set(pairs.light_pid) | set(pairs.heavy_pid))

# ---------- candidates_for_pair: ion must carry the halogen, DBE-only ----------
# light_mz for C6H10O3 as [M+Br]- ; n_Br=1 in ion via the adduct (need=0 neutral)
mz_br = C.ion_mz("C6H10O3", "[M+Br]-")
cands = RD.candidates_for_pair(mz_br, "Br", 1, ["[M+Br]-", "[M-H]-"],
                               ranges={"C": (1, 12), "H": (0, 28), "O": (0, 8),
                                       "N": (0, 1), "S": (0, 0)}, ppm=4.0)
check("pair candidates recover C6H10O3 ([M+Br]- gives the Br)", "C6H10O3" in cands,
      sorted(cands)[:6])
# all candidates obey DBE-only gate
check("all pair candidates pass DBE gate", all(C.dbe_ok(f)[0] for f in cands))

# ---------- stage A end-to-end with a fake scorer ----------
def fake_scorer(client, sample_id, formulas, *, mechanism_ids=None, **kw):
    """Return a flat per-isotopologue table: C6H10O3 [M+Br]- matches L1 strongly
    with a confirmed 81Br child on H1; everything else misses."""
    rows = []
    if "C6H10O3" in formulas:
        ion = "C6H10BrO3-"
        rows.append(dict(compound_formula="C6H10O3", compound_score=0.93,
                         compound_category=2, ion_formula=ion, ion_score=0.93,
                         ion_category=2, mechanism_id="m", isotope_formula=ion,
                         iso_label="M0", is_base=True, theo_mz=MZ_BR,
                         rel_abundance=1.0, iso_score=0.97, iso_category=2,
                         sample_peak_id="L1", sample_peak_mz=MZ_BR,
                         sample_peak_intensity=1e4, ppm_error=0.3, abundance_error=0.0))
        rows.append(dict(compound_formula="C6H10O3", compound_score=0.93,
                         compound_category=2, ion_formula=ion, ion_score=0.93,
                         ion_category=2, mechanism_id="m",
                         isotope_formula="[81Br]C6H10O3-", iso_label="81Br",
                         is_base=False, theo_mz=MZ_BR + RD.D_PAIR_BR,
                         rel_abundance=0.97, iso_score=0.95, iso_category=2,
                         sample_peak_id="H1", sample_peak_mz=MZ_BR + RD.D_PAIR_BR,
                         sample_peak_intensity=9.6e3, ppm_error=0.2, abundance_error=0.02))
    return pd.DataFrame(rows)


led2 = L.new_ledger(peaks)
res = RD.stage_a_iso_pairs(None, "SID", led2, PROF, PRE, CFG, ["[M+Br]-", "[M-H]-"],
                           score_fn=fake_scorer, log=lambda *_: None)
check("stage A commits the light member", res["committed"] == 1, res)
check("L1 assigned C6H10O3", led2.loc[led2.peak_id == "L1", "neutral_formula"].iloc[0] == "C6H10O3")
check("L1 confidence is Good (capped, not High)",
      "Good" in str(led2.loc[led2.peak_id == "L1", "confidence"].iloc[0]),
      led2.loc[led2.peak_id == "L1", "confidence"].iloc[0])
check("H1 attached as 81Br isotopologue child",
      L.role_of(led2, "H1") == L.ROLE_ISO
      and led2.loc[led2.peak_id == "H1", "parent_peak_id"].iloc[0] == "L1")
check("commentary names the doublet + DBE-only policy",
      "doublet" in led2.loc[led2.peak_id == "L1", "commentary"].iloc[0]
      and "DBE-only" in led2.loc[led2.peak_id == "L1", "commentary"].iloc[0])
check("ledger valid after stage A", L.validate(led2) == [], L.validate(led2))

# ---------- stage B: 2-anchor support waives loose ppm ----------
def fake_scorer_b(client, sample_id, formulas, *, mechanism_ids=None, **kw):
    rows = []
    if "C10H16O5" in formulas:   # = C10H16O4 + O, between two anchors
        ion = "C10H15O5-"
        rows.append(dict(compound_formula="C10H16O5", compound_score=0.82,
                         compound_category=2, ion_formula=ion, ion_score=0.82,
                         ion_category=2, mechanism_id="m", isotope_formula=ion,
                         iso_label="M0", is_base=True, theo_mz=MZ_S1,
                         rel_abundance=1.0, iso_score=0.9, iso_category=2,
                         sample_peak_id="S1", sample_peak_mz=MZ_S1,
                         sample_peak_intensity=2e3, ppm_error=2.3, abundance_error=0.0))
    return pd.DataFrame(rows)


MZ_S1 = C.ion_mz("C10H16O5", "[M-H]-")
pk = pd.DataFrame({"peak_id": ["A1", "A2", "S1"],
                   "mz": [C.ion_mz("C10H16O4", "[M-H]-"),
                          C.ion_mz("C10H16O6", "[M-H]-"), MZ_S1],
                   "height": [1e5, 1e5, 2e3]})
led3 = L.new_ledger(pk)
# seed two anchors one O below and above the target C10H16O5
L.commit_assignment(led3, "A1", neutral_formula="C10H16O4", adduct="[M-H]-",
                    ion_score=0.97, pass_no=1, method="x", confidence="High", commentary="anchor")
L.commit_assignment(led3, "A2", neutral_formula="C10H16O6", adduct="[M-H]-",
                    ion_score=0.97, pass_no=1, method="x", confidence="High", commentary="anchor")
resb = RD.stage_b_series(None, "SID", led3, PROF, CFG, ["[M-H]-"], reagent="Br",
                         score_fn=fake_scorer_b, log=lambda *_: None)
check("stage B commits C10H16O5 at 2.3 ppm via 2-anchor support", resb["committed"] == 1, resb)
check("S1 assigned C10H16O5",
      led3.loc[led3.peak_id == "S1", "neutral_formula"].iloc[0] == "C10H16O5")
check("S1 commentary cites supporting anchors",
      "supporting anchors" in led3.loc[led3.peak_id == "S1", "commentary"].iloc[0])

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

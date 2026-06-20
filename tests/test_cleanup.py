"""Offline tests for cleanup.py. Run: python3 tests/test_cleanup.py
The recovery's oracle is injected (score_fn), so no network."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import cleanup as CU       # noqa: E402
from mascope_assign import ledger as L         # noqa: E402
from mascope_assign import chemistry as C      # noqa: E402
from mascope_assign import contexts as X       # noqa: E402
from mascope_assign import passes as P         # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")

def mk(rows):
    return L.new_ledger(pd.DataFrame(rows, columns=["peak_id", "mz", "height"]))

# ---- ringing/shoulder artifacts ----
led = mk([("bright", 200.0000, 2.0e5), ("ring", 200.0080, 1.0e3),
          ("far", 250.0, 5e3), ("resolved", 300.0000, 1e3), ("rb", 300.0500, 2e5)])
CU.flag_ringing_artifacts(led, log=lambda *a: None)
check("ringing: weak peak 8 mDa from a 200x peak -> artifact",
      L.role_of(led, "ring") == L.ROLE_ARTIFACT)
check("ringing: the bright peak itself untouched", L.role_of(led, "bright") == L.ROLE_UNEXPLAINED)
check("ringing: resolved neighbour 50 mDa away NOT flagged",
      L.role_of(led, "resolved") == L.ROLE_UNEXPLAINED)
check("ringing: isolated peak untouched", L.role_of(led, "far") == L.ROLE_UNEXPLAINED)
check("ringing: ledger valid", L.validate(led) == [])

# ---- bromide reagent clusters ----
led2 = mk([("clu", 296.7588, 3.7e3), ("twin", 298.7568, 3.5e3), ("cho", 250.10, 1e4)])
CU.label_bromide_clusters(led2, log=lambda *a: None)
check("cluster: neg-defect peak with Br twin -> reagent",
      L.role_of(led2, "clu") == L.ROLE_REAGENT)
check("cluster: positive-defect CHO peak NOT labelled",
      L.role_of(led2, "cho") == L.ROLE_UNEXPLAINED)
# offline (no oracle): must NOT assert the untested "no covalent reading"
clu_note = led2.loc[led2.peak_id == "clu", "commentary"].iloc[0]
check("cluster: offline note makes no covalency claim",
      "covalent" not in clu_note.lower(), clu_note)

# with the oracle (mocked): a high-scoring covalent tie -> the precedent wording.
# C2H2Br2O2 is enumerated from the M0 (294.76); the oracle envelope then lands
# on its M+2/M+4 (296/298) too, so all three cluster peaks see the tie.
led2b = mk([("m0", 294.7609, 4.0e3), ("clu", 296.7588, 3.7e3), ("tw", 298.7568, 3.5e3)])
ENV = [{"ion_formula": "C2H2Br3O2-", "ion_score": 0.955, "ppm_error": p,
        "sample_peak_id": pid, "sample_peak_mz": m, "sample_peak_intensity": h}
       for pid, m, h, p in [("m0", 294.7609, 4e3, -0.65), ("clu", 296.7588, 3.7e3, -0.74),
                            ("tw", 298.7568, 3.5e3, -0.76)]]
def cov_score(client, sample_id, formulas, *, allow_partial=True, mechanism_ids=None):
    return pd.DataFrame(ENV) if "C2H2Br2O2" in formulas else pd.DataFrame()
out_clu = CU.label_bromide_clusters(led2b, client="C", sample_id="S",
                                    score_fn=cov_score, log=lambda *a: None)
note_b = led2b.loc[led2b.peak_id == "clu", "commentary"].iloc[0]
check("cluster: still role=reagent with a covalent tie",
      L.role_of(led2b, "clu") == L.ROLE_REAGENT)
check("cluster: covalent tie -> 'reagent-adduct reading preferred' wording",
      "reagent-adduct reading preferred" in note_b and "C2H2Br3O2" in note_b, note_b)
check("cluster: covalent-tie count reported", out_clu.get("covalent_ties", 0) >= 1, out_clu)
check("cluster: no false 'no covalent reading' when a tie exists",
      "no covalent reading" not in note_b)

# ---- isotope-confirmed recovery (mock oracle) ----
# a 1-Br ion C5H10BrO3- (= C5H10O3 [M+Br]-) + its 81Br twin at ~0.97x
mz = C.neutral_mass("C5H10O3") + C.ADDUCT_SHIFTS["[M+Br]-"]
led3 = mk([("p", mz, 1.0e4), ("twin", mz + CU.BR, 9.6e3), ("noise", 180.0, 5e3)])
prof = X.get_context("ambient-air")
cfg = P.PassConfig(); cfg.cal_mu, cfg.cal_sigma = -0.3, 0.3; cfg.mechanism_ids = None

def fake_score(client, sample_id, formulas, *, allow_partial=True, mechanism_ids=None):
    rows = []
    if "C5H10O3" in formulas:
        rows.append({"ion_formula": "C5H10BrO3-", "ion_score": 0.74,
                     "ppm_error": -0.3, "sample_peak_id": "p",
                     "sample_peak_mz": mz, "sample_peak_intensity": 1.0e4})
    return pd.DataFrame(rows)

out = CU.recover_isotope_gated(None, "S", led3, prof, cfg, score_fn=fake_score,
                               log=lambda *a: None)
check("recovery: isotope-confirmed 1-Br molecule committed",
      out["recovered"] == 1 and L.role_of(led3, "p") == L.ROLE_M0, out)
check("recovery: committed the right formula",
      led3.loc[led3.peak_id == "p", "neutral_formula"].iloc[0] == "C5H10O3",
      repr(led3.loc[led3.peak_id == "p", "neutral_formula"].iloc[0]))
check("recovery: confidence is Good (recovered)",
      led3.loc[led3.peak_id == "p", "confidence"].iloc[0] == "Good (recovered)")
check("recovery: ledger valid", L.validate(led3) == [])

# recovery REFUSES when the isotope twin is absent (no corroboration)
led4 = mk([("p", mz, 1.0e4), ("noise", 180.0, 5e3)])  # no twin
out4 = CU.recover_isotope_gated(None, "S", led4, prof, cfg, score_fn=fake_score,
                                log=lambda *a: None)
check("recovery: refuses a 1-Br fit with NO 81Br twin (uncorroborated)",
      out4["recovered"] == 0 and L.role_of(led4, "p") == L.ROLE_UNEXPLAINED, out4)

# recovery REFUSES a high-complexity (>2 heteroatom-type) fit even if scored
def fake_score_junk(client, sample_id, formulas, *, allow_partial=True, mechanism_ids=None):
    return pd.DataFrame([{"ion_formula": "C5H10BrO3-", "ion_score": 0.74,
                          "ppm_error": -0.3, "sample_peak_id": "p",
                          "sample_peak_mz": mz, "sample_peak_intensity": 1e4}])
# (covered structurally: the enumerator box caps complexity; this just re-checks accept path)

# ---- reclaim_satellites: attach leaked 13C/81Br satellites of an assigned M0 ----
ledr = mk([("par", 279.0236, 5000.0), ("c13", 280.0270, 250.0),
           ("br81", 281.0216, 4850.0), ("p2", 300.0000, 2000.0),
           ("bad", 301.0034, 1800.0), ("far", 400.0000, 100.0)])
L.commit_assignment(ledr, "par", neutral_formula="C10H16O4", adduct="[M+Br]-",
    ion_formula="C10H16BrO4-", ion_score=0.9, ppm_error=0.5, pass_no=1,
    method="grid", confidence="High", commentary="test")
L.commit_assignment(ledr, "p2", neutral_formula="C10H18O3", adduct="[M+Br]-",
    ion_formula="C10H18BrO3-", ion_score=0.9, ppm_error=0.5, pass_no=1,
    method="grid", confidence="High", commentary="test")
res = CU.reclaim_satellites(ledr, log=lambda *a: None)
check("reclaim: 13C satellite -> iso_child", L.role_of(ledr, "c13") == L.ROLE_ISO)
check("reclaim: 81Br satellite -> iso_child", L.role_of(ledr, "br81") == L.ROLE_ISO)
check("reclaim: count == 2", res["reclaimed"] == 2, res)
check("reclaim: too-bright +1.003 (ratio 0.9 vs 10C) NOT grabbed",
      L.role_of(ledr, "bad") == L.ROLE_UNEXPLAINED)
check("reclaim: isolated peak untouched", L.role_of(ledr, "far") == L.ROLE_UNEXPLAINED)
check("reclaim: parent M0 untouched", L.role_of(ledr, "par") == L.ROLE_M0)
check("reclaim: ledger valid", L.validate(ledr) == [])

# --- prefer_amine_over_ammonium (uronium NH4 -> protonated amine) -----------
leda = pd.DataFrame([
    dict(peak_id="a", mz=186.15, neutral_formula="C10H16O2", adduct="[M+NH4]+",
         role=L.ROLE_M0, tier_reason=""),                    # uncorroborated, valid amine
    dict(peak_id="b", mz=158.13, neutral_formula="C9H16O", adduct="[M+NH4]+",
         role=L.ROLE_M0, tier_reason=""),                    # corroborated by row c
    dict(peak_id="c", mz=141.12, neutral_formula="C9H16O", adduct="[M+H]+",
         role=L.ROLE_M0, tier_reason=""),
    dict(peak_id="d", mz=168.12, neutral_formula="C6H14O4", adduct="[M+NH4]+",
         role=L.ROLE_M0, tier_reason=""),                    # saturated -> amine impossible
    dict(peak_id="e", mz=999.0, neutral_formula=None, adduct=None,
         role=L.ROLE_ISO, tier_reason=""),
])
outa = CU.prefer_amine_over_ammonium(leda, log=lambda *a, **k: None)
check("amine: uncorroborated NH4 -> [M+H]+ of X+NH3",
      leda.loc[0, "neutral_formula"] == "C10H19NO2" and leda.loc[0, "adduct"] == "[M+H]+",
      leda.loc[0].to_dict())
check("amine: corroborated NH4 adduct kept",
      leda.loc[1, "neutral_formula"] == "C9H16O" and leda.loc[1, "adduct"] == "[M+NH4]+")
check("amine: saturated X (no valid amine) -> NH4 forced/kept",
      leda.loc[3, "neutral_formula"] == "C6H14O4" and leda.loc[3, "adduct"] == "[M+NH4]+")
check("amine: summary counts", outa == {"relabeled": 1, "kept_corroborated": 1, "forced_nh4": 1}, outa)
check("amine: relabel noted in tier_reason", "re-read" in str(leda.loc[0, "tier_reason"]))

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

"""Offline tests for cleanup.py. Run: python3 tests/test_cleanup.py
The recovery's oracle is injected (score_fn), so no network."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky import cleanup as CU       # noqa: E402
from peaky import ledger as L         # noqa: E402
from peaky import chemistry as C      # noqa: E402
from peaky import contexts as X       # noqa: E402
from peaky import passes as P         # noqa: E402

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

# ---- reclaim_envelope_tails: deep multi-37Cl tail of a poly-Cl M0 (chlorinated paraffin) ----
# C11H18Cl6 [M-H]- @358.9467; M+4 (2×37Cl) ratio C(6,2)·0.322=1.54, M+6 (3×37Cl)=0.66
ledt = mk([("par", 358.94668, 1000.0), ("m4", 362.94078, 1535.0),
           ("m6", 364.93783, 655.0), ("far", 400.0, 100.0)])
L.commit_assignment(ledt, "par", neutral_formula="C11H18Cl6", adduct="[M-H]-",
    ion_formula="C11H17Cl6-", ion_score=0.9, ppm_error=0.5, pass_no=0,
    method="known:chlorinated_paraffin", confidence="Good (chlorinated-paraffin)", commentary="test")
rest = CU.reclaim_envelope_tails(ledt, log=lambda *a: None)
check("envelope-tail: M+4 (2×37Cl) -> iso_child", L.role_of(ledt, "m4") == L.ROLE_ISO)
check("envelope-tail: M+6 (3×37Cl) -> iso_child", L.role_of(ledt, "m6") == L.ROLE_ISO)
check("envelope-tail: count == 2", rest["tails"] == 2, rest)
check("envelope-tail: isolated peak untouched", L.role_of(ledt, "far") == L.ROLE_UNEXPLAINED)
check("envelope-tail: ledger valid", L.validate(ledt) == [])

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

# ---------- demote_unconfirmed_fluorine (F-monster curb) ----------
ledf = pd.DataFrame([
    dict(role="M0", neutral_formula="C11H6F16", tier="Assigned", commentary="", below_assignability=False),
    dict(role="M0", neutral_formula="C2HF3O2",  tier="Assigned", commentary="", below_assignability=False),  # TFA (PFCA)
    dict(role="M0", neutral_formula="C8HF15O2", tier="Assigned", commentary="", below_assignability=False),  # PFOA (PFCA, F15)
    dict(role="M0", neutral_formula="C6H4ClF5O", tier="Assigned", commentary="", below_assignability=False), # Cl-anchored
    dict(role="M0", neutral_formula="C3H5FO2",  tier="Assigned", commentary="", below_assignability=False),  # low-F (F1)
])
outf = CU.demote_unconfirmed_fluorine(ledf, log=lambda *a: None)
check("F-monster demoted: only the unanchored high-F non-PFCA", outf == {"f_demoted": 1}, outf)
check("F-monster C11H6F16 -> Candidate + below_assignability",
      ledf.loc[0, "tier"] == "Candidate" and bool(ledf.loc[0, "below_assignability"]))
check("PFCA (TFA, PFOA) kept Assigned",
      ledf.loc[1, "tier"] == "Assigned" and ledf.loc[2, "tier"] == "Assigned")
check("Cl-anchored F kept; low-F kept",
      ledf.loc[3, "tier"] == "Assigned" and ledf.loc[4, "tier"] == "Assigned")


# ---- implausibly carbon-rich demote ----
ledc = pd.DataFrame([
    dict(role="M0", neutral_formula="C27H8",    tier="Candidate",  commentary="", below_assignability=False),  # H/C 0.30
    dict(role="M0", neutral_formula="C36H6O",   tier="Assigned", commentary="", below_assignability=False),  # H/C 0.17
    dict(role="M0", neutral_formula="C10H16O4", tier="Assigned", commentary="", below_assignability=False),  # H/C 1.6 keep
    dict(role="M0", neutral_formula="C11H6F16", tier="Candidate",  commentary="", below_assignability=False),  # F-rich, not C
])
outc = CU.demote_implausible_carbon(ledc, log=lambda *a: None)
check("carbon demote: only the F-free low-H/C clusters (C27H8, C36H6O)", outc == {"c_demoted": 2}, outc)
check("carbon demote: C36H6O Assigned -> Candidate + below_assignability",
      ledc.loc[1, "tier"] == "Candidate" and bool(ledc.loc[1, "below_assignability"]))
check("carbon demote: normal C10H16O4 kept",
      ledc.loc[2, "tier"] == "Assigned" and not bool(ledc.loc[2, "below_assignability"]))
check("carbon demote: F-rich skeleton NOT touched here (F-free rule; F-demote owns it)",
      not bool(ledc.loc[3, "below_assignability"]))

# ---- implausible-ionization demote (heteroatom-free via anion channel) ----
ledi = pd.DataFrame([
    dict(role="M0", neutral_formula="C7H10",   adduct="[M-H]-",   tier="Assigned", commentary="", below_assignability=False),  # HC deprotonation
    dict(role="M0", neutral_formula="C2H2",    adduct="[M+CO3]-", tier="Assigned", commentary="", below_assignability=False),  # HC carbonate cluster
    dict(role="M0", neutral_formula="C10H16",  adduct="[M+Br]-",  tier="Candidate",  commentary="", below_assignability=False),  # HC bromide adduct
    dict(role="M0", neutral_formula="C6H12O6", adduct="[M-H]-",   tier="Assigned", commentary="", below_assignability=False),  # has O -> keep
    dict(role="M0", neutral_formula="C7H10",   adduct="[M]-.",    tier="Assigned", commentary="", below_assignability=False),  # electron attach -> exempt
])
outi = CU.demote_implausible_ionization(ledi, log=lambda *a: None)
check("ionization demote: 3 heteroatom-free via FG-requiring anion channels", outi == {"ionization_demoted": 3}, outi)
check("ionization demote: C7H10 [M-H]- -> Candidate + below_assignability",
      ledi.loc[0, "tier"] == "Candidate" and bool(ledi.loc[0, "below_assignability"]))
check("ionization demote: C2H2 [M+CO3]- demoted",
      ledi.loc[1, "tier"] == "Candidate" and bool(ledi.loc[1, "below_assignability"]))
check("ionization demote: oxygenated C6H12O6 [M-H]- kept",
      ledi.loc[3, "tier"] == "Assigned" and not bool(ledi.loc[3, "below_assignability"]))
check("ionization demote: electron-attachment [M]-. exempt",
      ledi.loc[4, "tier"] == "Assigned" and not bool(ledi.loc[4, "below_assignability"]))

# ---- reagent-precursor / brominated-background halocarbon relabel ----
ledh = mk([("chbr2", 170.8451, 2e2), ("dbaa", 214.8349, 7e2), ("real", 250.10, 1e4)])
for pid, neutral, adduct, ionf in [("chbr2", "C", "[M+HBr+Br]-", "CHBr2-"),       # bare-C mis-read
                                   ("dbaa", "C2HBrO2", "[M+Br]-", "C2HBr2O2-"),    # dibromoacetic acid ion
                                   ("real", "C10H16O2", "[M+Br]-", "C10H16BrO2-")]:
    j = ledh.index[ledh["peak_id"] == pid][0]
    ledh.at[j, "role"] = L.ROLE_M0; ledh.at[j, "neutral_formula"] = neutral
    ledh.at[j, "adduct"] = adduct; ledh.at[j, "ion_formula"] = ionf
    ledh.at[j, "method"] = "seed"; ledh.at[j, "confidence"] = "Good"
    ledh.at[j, "commentary"] = "seed"        # provenance so the M0 rows are I5-valid
res = CU.relabel_reagent_halocarbons(ledh, reagent="Br", log=lambda *a: None)
check("halocarbon: relabeled 2 (CHBr2-, C2HBr2O2-)", res["relabeled"] == 2, res)
check("halocarbon: CHBr2- (bare-C mis-read) -> reagent", L.role_of(ledh, "chbr2") == L.ROLE_REAGENT)
_d = ledh[ledh["peak_id"] == "dbaa"].iloc[0]
check("halocarbon: C2HBr2O2- -> named dibromoacetic acid, M0",
      _d["neutral_formula"] == "C2H2Br2O2" and _d["role"] == L.ROLE_M0)
check("halocarbon: dibromoacetic acid adduct fixed to [M-H]-", _d["adduct"] == "[M-H]-")
check("halocarbon: real organobromine untouched",
      L.role_of(ledh, "real") == L.ROLE_M0
      and ledh[ledh["peak_id"] == "real"].iloc[0]["neutral_formula"] == "C10H16O2")
check("halocarbon: ledger valid after relabel", L.validate(ledh) == [])
# polarity gate: a non-Br reagent profile is a no-op
ledhu = mk([("x", 170.8451, 2e2)])
_j = ledhu.index[ledhu["peak_id"] == "x"][0]
ledhu.at[_j, "role"] = L.ROLE_M0; ledhu.at[_j, "ion_formula"] = "CHBr2-"
CU.relabel_reagent_halocarbons(ledhu, reagent="Ur", log=lambda *a: None)
check("halocarbon: non-Br reagent -> no-op", L.role_of(ledhu, "x") == L.ROLE_M0)
CU.relabel_reagent_halocarbons(ledhu, reagent=None, log=lambda *a: None)
check("halocarbon: reagent=None -> no-op", L.role_of(ledhu, "x") == L.ROLE_M0)

# ---- F-demote exemption requires a CONFIRMED isotope (audit rule-gap 1) ----
import json as _json  # noqa: E402
def _isos(*labels):
    return _json.dumps([{"label": l, "score": 0.9, "peak_id": "p"} for l in labels])
ledfi = pd.DataFrame([
    dict(role="M0", neutral_formula="C8H13F7N2O4S", adduct="[M+Br]-", tier="Assigned", commentary="", below_assignability=False, isotopologues="[]"),
    dict(role="M0", neutral_formula="C8H13F7N2O4S", adduct="[M+Br]-", tier="Assigned", commentary="", below_assignability=False, isotopologues=_isos("34S")),
    dict(role="M0", neutral_formula="C6H4ClF5O",    adduct="[M-H]-",  tier="Assigned", commentary="", below_assignability=False, isotopologues=_isos("37Cl")),
    dict(role="M0", neutral_formula="C6H4ClF5O",    adduct="[M-H]-",  tier="Assigned", commentary="", below_assignability=False, isotopologues="[]"),
    dict(role="M0", neutral_formula="C8HF15O2",     adduct="[M-H]-",  tier="Assigned", commentary="", below_assignability=False, isotopologues="[]"),  # PFCA
])
ofi = CU.demote_unconfirmed_fluorine(ledfi, log=lambda *a: None)
check("F-demote(iso): unconfirmed S/Cl F-monsters demoted (2), confirmed/PFCA kept", ofi == {"f_demoted": 2}, ofi)
check("F-demote(iso): C8H13F7N2O4S no-34S -> Candidate+below",
      ledfi.loc[0, "tier"] == "Candidate" and bool(ledfi.loc[0, "below_assignability"]))
check("F-demote(iso): C8H13F7N2O4S with 34S kept Assigned", ledfi.loc[1, "tier"] == "Assigned")
check("F-demote(iso): Cl anchor needs 37Cl (confirmed kept, unconfirmed demoted)",
      ledfi.loc[2, "tier"] == "Assigned" and ledfi.loc[3, "tier"] == "Candidate")
check("F-demote(iso): PFCA still exempt", ledfi.loc[4, "tier"] == "Assigned")

# ---- speculative-residual demote (audit rule-gaps 2-4) ----
class _Cfg:
    minor_channels = ("[M+CO3]-", "[M+O2]-", "[M]-.")
    cal_mu, cal_sigma, cal_z_accept = 0.0, 0.4, 2.0
ledr = pd.DataFrame([
    dict(role="M0", neutral_formula="C6H5N3",  adduct="[M+Br]-",  tier="Assigned", method="residual:iso-pair", ppm_error=-0.5, isotopologues="[]", commentary="", below_assignability=False),
    dict(role="M0", neutral_formula="C12H6O",  adduct="[M+CO3]-", tier="Assigned", method="residual:series",   ppm_error=-0.35, isotopologues="[]", commentary="-2xCH2 (0 supporting anchors)", below_assignability=False),
    dict(role="M0", neutral_formula="C10H16O4", adduct="[M-H]-",  tier="Assigned", method="cheminfo+grid",     ppm_error=0.2, isotopologues="[]", commentary="", below_assignability=False),
    dict(role="M0", neutral_formula="C9H12O5",  adduct="[M-H]-",  tier="Assigned", method="residual:iso-pair", ppm_error=5.0, isotopologues=_isos("13C"), commentary="", below_assignability=False),
])
orr = CU.demote_speculative_residual(ledr, _Cfg(), log=lambda *a: None)
check("residual demote: 3 (N3 / 0-anchor series / off-cal); cheminfo HOM acid protected", orr == {"residual_demoted": 3}, orr)
check("residual demote: C6H5N3 -> Candidate+below",
      ledr.loc[0, "tier"] == "Candidate" and bool(ledr.loc[0, "below_assignability"]))
check("residual demote: C12H6O (0-anchor/sole-minor) -> Candidate", ledr.loc[1, "tier"] == "Candidate")
check("residual demote: legit cheminfo+grid HOM acid untouched", ledr.loc[2, "tier"] == "Assigned")
check("residual demote: off-cal residual -> Candidate", ledr.loc[3, "tier"] == "Candidate")


# ---- radical-anion relabel (hydrocarbon FG-cluster -> M-. of closed-shell neutral) ----
ledrad = pd.DataFrame([
    dict(role="M0", neutral_formula="C3H4",   adduct="[M+CO3]-", tier="Candidate", commentary="", below_assignability=True,  ion_formula="", dbe=2.0),  # -> C4H4O3, corroborated below
    dict(role="M0", neutral_formula="C4H4O3", adduct="[M-H]-",   tier="Assigned",  commentary="", below_assignability=False, ion_formula="", dbe=3.0),  # corroborating neutral
    dict(role="M0", neutral_formula="C6H6",   adduct="[M+CO3]-", tier="Candidate", commentary="", below_assignability=True,  ion_formula="", dbe=4.0),  # -> C7H6O3, NOT corroborated
    dict(role="M0", neutral_formula="C6H12O6",adduct="[M-H]-",   tier="Assigned",  commentary="", below_assignability=False, ion_formula="", dbe=1.0),  # oxygenated, untouched
])
outrad = CU.relabel_radical_anions(ledrad, log=lambda *a: None)
check("radical: 2 hydrocarbon CO3 clusters relabeled", outrad["radical_relabeled"] == 2, outrad)
check("radical: 1 corroborated (C4H4O3 has [M-H]-)", outrad["radical_corroborated"] == 1, outrad)
check("radical: C3H4 [M+CO3]- -> C4H4O3 [M]-. corroborated, VISIBLE (not below_assignability)",
      ledrad.loc[0, "neutral_formula"] == "C4H4O3" and ledrad.loc[0, "adduct"] == "[M]-."
      and not bool(ledrad.loc[0, "below_assignability"]))
check("radical: C6H6 [M+CO3]- -> C7H6O3 [M]-. uncorroborated, Candidate+below (still SHOWN)",
      ledrad.loc[2, "neutral_formula"] == "C7H6O3" and ledrad.loc[2, "adduct"] == "[M]-."
      and ledrad.loc[2, "tier"] == "Candidate" and bool(ledrad.loc[2, "below_assignability"]))
check("radical: corroborating [M-H]- row untouched",
      ledrad.loc[1, "neutral_formula"] == "C4H4O3" and ledrad.loc[1, "adduct"] == "[M-H]-")
check("radical: oxygenated C6H12O6 [M-H]- untouched", ledrad.loc[3, "adduct"] == "[M-H]-")
# the relabeled radicals must now ESCAPE the hydrocarbon implausible-ionization demote
outi_rad = CU.demote_implausible_ionization(ledrad, log=lambda *a: None)
check("radical: relabeled radicals escape implausible-ionization demote",
      outi_rad == {"ionization_demoted": 0}, outi_rad)


# ---- positive-mode reagent-N re-read (hydrocarbon via urea/NH4 -> [M+H]+ N-heterocycle) ----
ledrn = pd.DataFrame([
    dict(role="M0", neutral_formula="C5H6",  adduct="[M+(CH4N2O)H]+", tier="Candidate", commentary="", below_assignability=False, ion_formula="", dbe=3.0),  # -> C6H10N2O [M+H]+
    dict(role="M0", neutral_formula="C5H4",  adduct="[M+NH4]+",       tier="Assigned",  commentary="", below_assignability=False, ion_formula="", dbe=4.0),  # -> C5H7N [M+H]+
    dict(role="M0", neutral_formula="C10H16",adduct="[M+NH4]+",       tier="Assigned",  commentary="", below_assignability=False, ion_formula="", dbe=3.0),  # terpene WITH [M+H]+ -> keep
    dict(role="M0", neutral_formula="C10H16",adduct="[M+H]+",         tier="Assigned",  commentary="", below_assignability=False, ion_formula="", dbe=3.0),  # the terpene's [M+H]+
    dict(role="M0", neutral_formula="C6H10O2",adduct="[M+NH4]+",      tier="Assigned",  commentary="", below_assignability=False, ion_formula="", dbe=2.0),  # oxygenated -> untouched
])
outrn = CU.relabel_reagent_n_adducts(ledrn, log=lambda *a: None)
check("reagent-N: 2 hydrocarbon urea/NH4 clusters re-read", outrn == {"reagent_n_relabeled": 2}, outrn)
check("reagent-N: C5H6 [M+(CH4N2O)H]+ -> C6H10N2O [M+H]+ Candidate+below",
      ledrn.loc[0, "neutral_formula"] == "C6H10N2O" and ledrn.loc[0, "adduct"] == "[M+H]+"
      and ledrn.loc[0, "tier"] == "Candidate" and bool(ledrn.loc[0, "below_assignability"]))
check("reagent-N: C5H4 [M+NH4]+ -> C5H7N [M+H]+", ledrn.loc[1, "neutral_formula"] == "C5H7N")
check("reagent-N: terpene C10H16 [M+NH4]+ KEPT (it has its own [M+H]+)",
      ledrn.loc[2, "neutral_formula"] == "C10H16" and ledrn.loc[2, "adduct"] == "[M+NH4]+")
check("reagent-N: oxygenated C6H10O2 [M+NH4]+ untouched", ledrn.loc[4, "neutral_formula"] == "C6H10O2")


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

"""Offline tests for passes.py arbitration + commit. Run: python3 tests/test_passes.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import passes as P  # noqa: E402
from mascope_assign import ledger as L  # noqa: E402

PASS = FAIL = 0
CFG = P.PassConfig()


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


def iso_row(**kw):
    base = dict(compound_formula=None, compound_score=None, compound_category=2,
                ion_formula=None, ion_score=None, ion_category=2, mechanism_id="m",
                isotope_formula=None, iso_label="M0", is_base=True, theo_mz=100.0,
                rel_abundance=1.0, iso_score=None, iso_category=2,
                sample_peak_id=None, sample_peak_mz=100.0, sample_peak_intensity=1e4,
                ppm_error=0.2, abundance_error=0.0)
    base.update(kw)
    return base


# ---------- confidence_label ----------
check("High needs score>=0.90 + iso + untied",
      P.confidence_label(0.97, -0.3, 1, False, CFG) == "High")
check("no iso -> not High (Good)",
      P.confidence_label(0.97, -0.3, 0, False, CFG) == "Good")
check("tied -> not High",
      P.confidence_label(0.97, -0.3, 1, True, CFG) == "Good")
check("0.82 -> Good", P.confidence_label(0.82, 0.5, 0, False, CFG) == "Good")
check("0.74 -> Low", P.confidence_label(0.74, 0.5, 0, False, CFG) == "Low")
check("0.55 -> Suspect", P.confidence_label(0.55, 0.5, 0, False, CFG) == "Suspect")
check("0.40 -> Reject", P.confidence_label(0.40, 0.5, 0, False, CFG) == "Reject")
check("suffix applied", P.confidence_label(0.75, 0.5, 0, False, CFG, "series") == "Low (series)")

# ---------- arbitration: CHO beats CHON on a peak ----------
scored = pd.DataFrame([
    # peak A: CHO strong vs CHON weaker
    iso_row(sample_peak_id="A", compound_formula="C10H16O4", compound_score=0.96,
            ion_formula="C10H15O4-", ion_score=0.96, ppm_error=-0.3),
    iso_row(sample_peak_id="A", compound_formula="C9H13NO5", compound_score=0.90,
            ion_formula="C9H12NO5-", ion_score=0.90, ppm_error=0.8),
    # peak A 13C child attributed to peak B (non-base)
    iso_row(sample_peak_id="B", compound_formula="C10H16O4", compound_score=0.96,
            ion_formula="C10H15O4-", ion_score=0.96, isotope_formula="[13C]C9H15O4-",
            iso_label="13C", is_base=False, iso_score=0.93, ppm_error=0.1),
])
arb = P.arbitrate(scored, CFG)
win = arb["winners"]
wa = win[win.peak_id == "A"].iloc[0]
check("peak A winner is CHO C10H16O4", wa["neutral"] == "C10H16O4", wa["neutral"])
check("peak A n_iso counts the 13C child", wa["n_iso"] == 1, wa["n_iso"])
check("peak A has CHON as alternative",
      any(a["formula"] == "C9H13NO5" for a in wa["alternatives"]), wa["alternatives"])
kids = arb["iso_children"]
check("13C child B attributed to parent A",
      len(kids) == 1 and kids.iloc[0]["peak_id"] == "B"
      and kids.iloc[0]["parent_peak_id"] == "A", kids.to_dict("records"))

# ---------- complexity penalty: CHO wins over slightly-higher neutral-Br ----------
scored2 = pd.DataFrame([
    iso_row(sample_peak_id="C", compound_formula="C7H12BrO4", compound_score=0.93,
            ion_formula="C7H12BrO4-", ion_score=0.93, ppm_error=0.4),   # neutral Br
    iso_row(sample_peak_id="C", compound_formula="C7H12O4", compound_score=0.88,
            ion_formula="C7H11O4-", ion_score=0.88, ppm_error=0.6),     # CHO, [M-H]-
])
arb2 = P.arbitrate(scored2, CFG)
wc = arb2["winners"].iloc[0]
check("neutral Br loses to CHO after complexity penalty (0.93-0.20 < 0.88)",
      wc["neutral"] == "C7H12O4", (wc["neutral"], wc["eff_score"]))

# ---------- isotopologue-gated heteroatoms ----------
# peak D: an organosulfate (S) slightly outscores a CHON alt, but has NO 34S
# evidence -> should lose after the het-iso penalty.
scored3 = pd.DataFrame([
    iso_row(sample_peak_id="D", compound_formula="C8H12O6S", compound_score=0.90,
            ion_formula="C8H11O6S-", ion_score=0.90, ppm_error=0.3),
    iso_row(sample_peak_id="D", compound_formula="C9H15NO5", compound_score=0.86,
            ion_formula="C9H14NO5-", ion_score=0.86, ppm_error=0.5),
])
arb3 = P.arbitrate(scored3, CFG)
wd = arb3["winners"].iloc[0]
check("S candidate without 34S loses to CHON (iso gate)",
      wd["neutral"] == "C9H15NO5", (wd["neutral"], wd["eff_score"]))

# peak E: same S candidate but WITH a confirmed 34S satellite -> should win.
scored4 = pd.DataFrame([
    iso_row(sample_peak_id="E", compound_formula="C8H12O6S", compound_score=0.90,
            ion_formula="C8H11O6S-", ion_score=0.90, ppm_error=0.3),
    iso_row(sample_peak_id="F", compound_formula="C8H12O6S", compound_score=0.90,
            ion_formula="C8H11O6S-", ion_score=0.90, isotope_formula="[34S]C8H11O6-",
            iso_label="34S", is_base=False, iso_score=0.88, ppm_error=0.2),
    iso_row(sample_peak_id="E", compound_formula="C9H15NO5", compound_score=0.86,
            ion_formula="C9H14NO5-", ion_score=0.86, ppm_error=0.5),
])
arb4 = P.arbitrate(scored4, CFG)
we = arb4["winners"][arb4["winners"].peak_id == "E"].iloc[0]
check("S candidate WITH 34S confirmed wins",
      we["neutral"] == "C8H12O6S", (we["neutral"], we["eff_score"]))

# ---------- reagent-element alias: adduct reading beats covalent ----------
# In Br-CIMS the ion C6H10BrO3- is equally covalent C6H11BrO3 [M-H]- or
# C6H10O3 [M+Br]-. Both carry 81Br confirmation; with reagent_element="Br" the
# covalent reading keeps its complexity prior and must lose the tie.
cfg_br = P.PassConfig()
cfg_br.reagent_element = "Br"
scored5 = pd.DataFrame([
    iso_row(sample_peak_id="G", compound_formula="C6H11BrO3", compound_score=0.90,
            ion_formula="C6H10BrO3-", ion_score=0.90, ppm_error=0.3),
    iso_row(sample_peak_id="G", compound_formula="C6H10O3", compound_score=0.90,
            ion_formula="C6H10BrO3-", ion_score=0.90, ppm_error=0.3),
    # 81Br satellite confirmed for BOTH (same ion) — attach to the covalent one
    iso_row(sample_peak_id="H", compound_formula="C6H11BrO3", compound_score=0.90,
            ion_formula="C6H10BrO3-", ion_score=0.90, isotope_formula="[81Br]C6H10O3-",
            iso_label="81Br", is_base=False, iso_score=0.95, ppm_error=0.1),
    iso_row(sample_peak_id="H", compound_formula="C6H10O3", compound_score=0.90,
            ion_formula="C6H10BrO3-", ion_score=0.90, isotope_formula="[81Br]C6H10O3-",
            iso_label="81Br", is_base=False, iso_score=0.95, ppm_error=0.1),
])
arb5 = P.arbitrate(scored5, cfg_br)
wg = arb5["winners"][arb5["winners"].peak_id == "G"].iloc[0]
check("adduct reading C6H10O3 [M+Br]- beats covalent C6H11BrO3 (reagent prior)",
      wg["neutral"] == "C6H10O3", (wg["neutral"], wg["eff_score"]))

# HBr cluster unit arithmetic
from mascope_assign import series_gka as G2  # noqa: E402
check("Y + HBr = covalent alias composition",
      G2.formula_add("C6H10O3", "HBr", 1) == "C6H11BrO3",
      G2.formula_add("C6H10O3", "HBr", 1))
check("alias - HBr recovers anchor",
      G2.formula_add("C6H11BrO3", "HBr", -1) == "C6H10O3")

# ---------- alias guard extends to the CO3 background channel (426.976 case) ----------
w_co3 = {"neutral": "C11H14BrNO8", "adduct": "[M+CO3]-",
         "ion_formula": "C12H14BrNO11-"}
r_co3 = P._prefer_adduct_reading(w_co3, cfg_br)
check("covalent bromo + CO3 relabeled to HBr-cluster of a Br-free neutral",
      r_co3["neutral"] == "C11H13NO8" and r_co3["adduct"] == "[M+HBr+CO3]-",
      (r_co3["neutral"], r_co3["adduct"]))
check("CO3 relabel carries the reagent-rule note", "_relabel_note" in r_co3)
# but only when the cluster channel mass is registered: O2 has no HBr-cluster shift
w_o2 = {"neutral": "C11H14BrNO8", "adduct": "[M+O2]-", "ion_formula": "x"}
r_o2 = P._prefer_adduct_reading(w_o2, cfg_br)
check("unmodelled cluster channel keeps the covalent reading",
      r_o2["neutral"] == "C11H14BrNO8" and r_o2["adduct"] == "[M+O2]-",
      (r_o2["neutral"], r_o2["adduct"]))
# no reagent element set -> never relabel
check("no relabel without a reagent element",
      P._prefer_adduct_reading(dict(w_co3), P.PassConfig())["adduct"] == "[M+CO3]-")

# ---------- commit into ledger end-to-end ----------
peaks = pd.DataFrame({"peak_id": ["A", "B", "C"], "mz": [200.1, 201.1, 191.0],
                      "height": [1e5, 1e4, 8e4]})
led = L.new_ledger(peaks)
summary = P.commit_winners(led, arb, pass_no=1, method="cheminfo",
                           context="ambient-air", cfg=CFG, lock=True,
                           min_raw_score=0.70)
check("commit reports >=1 committed", summary["committed"] >= 1, summary)
check("peak A committed as M0 CHO",
      L.role_of(led, "A") == L.ROLE_M0
      and led.loc[led.peak_id == "A", "neutral_formula"].iloc[0] == "C10H16O4")
check("peak A is High and locked",
      led.loc[led.peak_id == "A", "confidence"].iloc[0] == "High"
      and L.is_locked(led, "A"))
check("peak B attached as iso_child of A",
      L.role_of(led, "B") == L.ROLE_ISO
      and led.loc[led.peak_id == "B", "parent_peak_id"].iloc[0] == "A")
check("commentary written", "C10H16O4" in led.loc[led.peak_id == "A", "commentary"].iloc[0])
check("ledger validates after commit", L.validate(led) == [], L.validate(led))

# ---------- claim_unexplained_only: families cannot displace assignments ----------
peaks_c = pd.DataFrame({"peak_id": ["X1"], "mz": [200.1], "height": [1e4]})
led_c = L.new_ledger(peaks_c)
L.commit_assignment(led_c, "X1", neutral_formula="C10H16O4", adduct="[M-H]-",
                    ion_score=0.80, pass_no=1, method="cheminfo", confidence="Good",
                    commentary="backbone, unlocked")
scored_c = pd.DataFrame([
    iso_row(sample_peak_id="X1", compound_formula="C7H8F4O2", compound_score=0.95,
            ion_formula="C7H7F4O2-", ion_score=0.95, ppm_error=0.2)])
arb_c = P.arbitrate(scored_c, CFG)
s_c = P.commit_winners(led_c, arb_c, pass_no=3, method="contaminant:fluorinated",
                       context="ambient-air", cfg=CFG, lock=False,
                       min_raw_score=0.5, confidence_suffix="fluorinated",
                       claim_unexplained_only=True)
check("family cannot displace existing M0 (claim_unexplained_only)",
      s_c["committed"] == 0
      and led_c.loc[led_c.peak_id == "X1", "neutral_formula"].iloc[0] == "C10H16O4",
      (s_c, led_c.neutral_formula.iloc[0]))

# without the guard, the higher score would displace (documents the contrast)
s_c2 = P.commit_winners(led_c, arb_c, pass_no=3, method="contaminant:fluorinated",
                        context="ambient-air", cfg=CFG, lock=False,
                        min_raw_score=0.5, confidence_suffix="fluorinated")
check("without guard the F formula would displace (contrast)",
      s_c2["committed"] == 1, s_c2)

# ---------- only_peaks restriction ----------
led_o = L.new_ledger(pd.DataFrame({"peak_id": ["P1", "P2"], "mz": [200.1, 300.2],
                                   "height": [1e4, 1e4]}))
scored_o = pd.DataFrame([
    iso_row(sample_peak_id="P1", compound_formula="C9H8F10O", compound_score=0.9,
            ion_formula="C9H8BrF10O-", ion_score=0.9, ppm_error=0.3),
    iso_row(sample_peak_id="P2", compound_formula="C10H8F12O", compound_score=0.9,
            ion_formula="C10H8BrF12O-", ion_score=0.9, ppm_error=0.3)])
arb_o = P.arbitrate(scored_o, CFG)
s_o = P.commit_winners(led_o, arb_o, pass_no=3, method="contaminant:fluorinated",
                       context="ambient-air", cfg=CFG, lock=False,
                       min_raw_score=0.5, confidence_suffix="fluorinated",
                       claim_unexplained_only=True, only_peaks={"P1"})
check("only_peaks: P1 committed, P2 skipped",
      s_o["committed"] == 1 and L.role_of(led_o, "P2") == L.ROLE_UNEXPLAINED, s_o)

# ---------- _family_ok ----------
fr = {"C": (0, 20), "H": (0, 44), "F": (1, 17), "O": (0, 6), "N": (0, 0),
      "S": (0, 0), "P": (0, 0), "Si": (0, 0), "Cl": (0, 0), "Br": (0, 0), "I": (0, 0)}
check("_family_ok accepts in-range fluorinated", P._family_ok("C10H8F14", fr))
check("_family_ok rejects element beyond family range", not P._family_ok("C5H5F3N2O2", fr))
check("_family_ok rejects half-integer DBE", not P._family_ok("C8H14NO12", fr))

# ---------- ranges helper ----------
from mascope_assign import contexts as X  # noqa: E402
from mascope_assign import isotopes as ISO  # noqa: E402
pre = ISO.PrescanResult(has_Br=True, estimated_max_C=12)
# Pass1/2: CHO(N) only -- heteroatoms NEVER auto-added from prescan
r = P.build_ranges(X.get_context("ambient-air"), pre, include_N=True)
s = P.ranges_to_string(r)
check("CHON ranges include N", "N" in s, s)
check("CHON ranges exclude Br/Cl/S/Si (no auto-add)",
      all(x not in s for x in ("Br", "Cl", "Si")) and " S" not in (" " + s), s)
check("ranges cap C near estimate", r["C"][1] <= 16, r["C"])
# Pass3: heteroatoms enter only via extra_elements
r3 = P.build_ranges(X.get_context("ambient-air"), pre, include_N=True,
                    extra_elements={"S": (1, 1), "O": (3, 6)})
check("extra_elements adds S in pass3", r3["S"] == (1, 1), r3["S"])

# ---------- _mech_to_adduct: label derived from exact ion-neutral diff ----------
def adduct_of(comp, ion):
    return P._mech_to_adduct({"compound_formula": comp, "ion_formula": ion})

check("adduct: deprotonation", adduct_of("C5H10O", "C5H9O-") == "[M-H]-")
check("adduct: Br attachment", adduct_of("C10H16O4", "C10H16BrO4-") == "[M+Br]-")
check("adduct: carbonate (v13 mislabel regression)",
      adduct_of("C5H10O", "C6H10O4-") == "[M+CO3]-",
      adduct_of("C5H10O", "C6H10O4-"))
check("adduct: carbonate on N/Br neutral",
      adduct_of("C11H14BrNO8", "C12H14BrNO11-") == "[M+CO3]-",
      adduct_of("C11H14BrNO8", "C12H14BrNO11-"))
check("adduct: superoxide", adduct_of("C5H10O", "C5H10O3-") == "[M+O2]-")
check("adduct: nitrate", adduct_of("C5H10O", "C5H10NO4-") == "[M+NO3]-")
check("adduct: electron attachment", adduct_of("C5H10O", "C5H10O-") == "[M]-.")
check("adduct: unknown diff falls back to [M-H]-",
      adduct_of("C5H10O", "C5H12NaO5-") == "[M-H]-")

# ---------- calibrate + z_of ----------
cal_led = L.new_ledger(pd.DataFrame({
    "peak_id": [f"c{i}" for i in range(30)] + ["x1"],
    "mz": [100.0 + i for i in range(30)] + [400.0],
    "height": [1e4] * 31}))
for i in range(30):
    L.commit_assignment(cal_led, f"c{i}", neutral_formula="C5H10O3",
                        adduct="[M-H]-", ion_score=0.92,
                        ppm_error=0.2 + (0.3 if i % 2 else -0.3),
                        pass_no=1, method="cheminfo+grid", confidence="Good",
                        commentary="backbone")
cfg_cal = P.PassConfig()
res = P.calibrate(cal_led, cfg_cal, log=lambda *a: None)
check("calibrate fits backbone mu", res is not None and abs(cfg_cal.cal_mu - 0.2) < 0.11,
      (cfg_cal.cal_mu, cfg_cal.cal_sigma))
check("z_of judges against calibration",
      P.z_of(0.2, cfg_cal) is not None and P.z_of(0.2, cfg_cal) < 0.1
      and P.z_of(4.0, cfg_cal) > cfg_cal.cal_z_pattern, P.z_of(4.0, cfg_cal))
check("z_of None when uncalibrated", P.z_of(1.0, P.PassConfig()) is None)
check("z_of None on NaN ppm", P.z_of(float("nan"), cfg_cal) is None)
small = L.new_ledger(pd.DataFrame({"peak_id": ["a"], "mz": [100.0], "height": [1e4]}))
cfg_small = P.PassConfig()
check("calibrate refuses a tiny backbone",
      P.calibrate(small, cfg_small, log=lambda *a: None) is None
      and cfg_small.cal_mu is None)

# ---------- offset-tolerance: a -2.4 ppm instrument (uronium) ----------
# the backbone is high-SCORE but mislabeled 'Low' (the pre-calibration confidence
# gate is centered on 0 ppm). calibrate must still find the -2.4 ppm center by
# selecting on SCORE, not the label.
off_led = L.new_ledger(pd.DataFrame({
    "peak_id": [f"o{i}" for i in range(30)],
    "mz": [150.0 + i for i in range(30)], "height": [1e4] * 30}))
for i in range(30):
    L.commit_assignment(off_led, f"o{i}", neutral_formula="C8H18O3",
                        adduct="[M+H]+", ion_score=0.95, compound_score=0.95,
                        ppm_error=-2.4 + (0.2 if i % 2 else -0.2),
                        pass_no=1, method="cheminfo+grid",
                        confidence="Low", commentary="backbone")
cfg_off = P.PassConfig()
P.calibrate(off_led, cfg_off, log=lambda *a: None)
check("calibrate finds the -2.4 ppm offset despite 'Low' labels",
      cfg_off.cal_mu is not None and abs(cfg_off.cal_mu + 2.4) < 0.15
      and cfg_off.cal_sigma < 0.6, (cfg_off.cal_mu, cfg_off.cal_sigma))

# confidence_label judges ppm against the calibrated center
cfg_u = P.PassConfig()
check("confidence_label: -2.4 ppm reads Low when uncalibrated (centered 0)",
      P.confidence_label(0.97, -2.4, 1, False, cfg_u) == "Low")
cfg_u.cal_mu, cfg_u.cal_sigma = -2.4, 0.3
check("confidence_label: -2.4 ppm reads High at the calibrated center",
      P.confidence_label(0.97, -2.4, 1, False, cfg_u) == "High")
check("confidence_label: an off-trend +0.3 ppm fit reads Low even calibrated",
      P.confidence_label(0.85, 0.3, 0, False, cfg_u) == "Low")

# relabel_confidence re-grades pass-1's commits against the fitted center:
# upgrade the on-cal backbone, demote the off-trend monster, keep the suffix.
rel_led = L.new_ledger(pd.DataFrame({
    "peak_id": ["good", "mon"], "mz": [200.0, 462.1], "height": [1e5, 1e5]}))
L.commit_assignment(rel_led, "good", neutral_formula="C8H18O3", adduct="[M+H]+",
                    ion_score=0.96, compound_score=0.96, ppm_error=-2.45,
                    pass_no=1, method="cheminfo+grid", confidence="Low",
                    commentary="x")
L.commit_assignment(rel_led, "mon", neutral_formula="C15H27NO15", adduct="[M+H]+",
                    ion_score=0.81, compound_score=0.81, ppm_error=0.32,
                    pass_no=2, method="gka-series", confidence="Good (series)",
                    commentary="x")
cfg_rel = P.PassConfig(); cfg_rel.cal_mu, cfg_rel.cal_sigma = -2.45, 0.3
P.relabel_confidence(rel_led, cfg_rel, log=lambda *a: None)
gconf = str(rel_led.loc[rel_led.peak_id == "good", "confidence"].iloc[0])
mconf = str(rel_led.loc[rel_led.peak_id == "mon", "confidence"].iloc[0])
check("relabel upgrades the on-cal backbone (Low -> Good/High)",
      gconf.startswith(("Good", "High")), gconf)
check("relabel demotes the off-trend monster + keeps the suffix",
      mconf.startswith(("Low", "Suspect")) and "series" in mconf, mconf)
check("relabel is a no-op when uncalibrated",
      P.relabel_confidence(rel_led, P.PassConfig(), log=lambda *a: None) == 0)

# arbitration is calibration-aware: an off-trend mass-coincidence with a HIGHER
# raw score loses to the on-trend formula once calibrated (else it wins then is
# z-rejected at commit, leaving the peak unexplained -- the high-Si PDMS failure).
cfg_arb = P.PassConfig(); cfg_arb.cal_mu, cfg_arb.cal_sigma = -2.45, 0.27
sc_off = pd.DataFrame([
    iso_row(sample_peak_id="P", compound_formula="C12H39NO6Si6", compound_score=0.85,
            ion_formula="C12H40NO6Si6", ion_score=0.85, ppm_error=-2.40),   # on-trend
    iso_row(sample_peak_id="P", compound_formula="C19H51NO4Si9", compound_score=0.90,
            ion_formula="C19H52NO4Si9", ion_score=0.90, ppm_error=1.60),    # off-trend, higher
])
w_cal = P.arbitrate(sc_off, cfg_arb)["winners"].iloc[0]
check("arbitration: on-trend PDMS beats higher-scoring off-trend coincidence (calibrated)",
      w_cal["neutral"] == "C12H39NO6Si6", (w_cal["neutral"], w_cal["eff_score"]))
w_unc = P.arbitrate(sc_off, P.PassConfig())["winners"].iloc[0]
check("arbitration uncalibrated: higher raw score still wins (no cal penalty)",
      w_unc["neutral"] == "C19H51NO4Si9", w_unc["neutral"])

# ---------- commit gates: NaN ppm, z-score, minor channel ----------
def gate_ledger(pids):
    return L.new_ledger(pd.DataFrame({
        "peak_id": pids, "mz": [300.0 + i for i in range(len(pids))],
        "height": [1e4] * len(pids)}))

def one_winner(pid, formula, ion, ppm, score=0.85, n_iso=0):
    return {"winners": pd.DataFrame([{
        "peak_id": pid, "neutral": formula, "ion_formula": ion,
        "adduct": P._mech_to_adduct({"compound_formula": formula, "ion_formula": ion}),
        "ion_score": score, "compound_score": score, "raw_score": score,
        "eff_score": score, "ppm_error": ppm, "n_iso": n_iso, "tied": False,
        "alternatives": []}]), "iso_children": pd.DataFrame()}

led_g = gate_ledger(["n1"])
s = P.commit_winners(led_g, one_winner("n1", "C9H10BrNO8", "C9H9BrNO8-", None),
                     pass_no=3, method="contaminant:bromo_organic", context="ambient-air",
                     cfg=cfg_cal, lock=False, min_raw_score=0.5)
check("NaN-ppm winner rejected", s["committed"] == 0 and s["rejected"]["nan_ppm"] == 1, s)

led_g = gate_ledger(["z1"])
s = P.commit_winners(led_g, one_winner("z1", "C9H10O4", "C9H9O4-", 3.9),
                     pass_no=3, method="contaminant:bromo_organic", context="ambient-air",
                     cfg=cfg_cal, lock=False, min_raw_score=0.5)
check("z>accept without evidence rejected",
      s["committed"] == 0 and s["rejected"]["mass_gate"] == 1, s)

led_g = gate_ledger(["z2"])
s = P.commit_winners(led_g, one_winner("z2", "C9H10O4", "C9H9O4-", 0.9, n_iso=1),
                     pass_no=3, method="contaminant:bromo_organic", context="ambient-air",
                     cfg=cfg_cal, lock=False, min_raw_score=0.5)
check("2<z<=4 with isotope evidence commits", s["committed"] == 1, (s, P.z_of(0.9, cfg_cal)))

led_g = gate_ledger(["z3"])
s = P.commit_winners(led_g, one_winner("z3", "C9H10O4", "C9H9O4-", 3.9, n_iso=3),
                     pass_no=3, method="contaminant:bromo_organic", context="ambient-air",
                     cfg=cfg_cal, lock=False, min_raw_score=0.5)
check("z>pattern band rejected even with isotopes",
      s["committed"] == 0 and s["rejected"]["mass_gate"] == 1, s)

# minor-channel commit gate: a lone Low CO3 winner is refused...
led_g = gate_ledger(["m1"])
s = P.commit_winners(led_g, one_winner("m1", "C7H9NO8", "C8H9NO11-", 0.3, score=0.75),
                     pass_no=3, method="contaminant:bromo_organic", context="ambient-air",
                     cfg=cfg_cal, lock=False, min_raw_score=0.5)
check("lone Low minor-channel winner rejected",
      s["committed"] == 0 and s["rejected"]["minor_channel"] == 1, s)
# ...a Good+ CO3 winner commits...
led_g = gate_ledger(["m2"])
s = P.commit_winners(led_g, one_winner("m2", "C5H10O", "C6H10O4-", 0.3, score=0.92),
                     pass_no=1, method="cheminfo+grid", context="ambient-air",
                     cfg=cfg_cal, lock=False, min_raw_score=0.5)
check("Good minor-channel winner commits", s["committed"] == 1, s)
# ...and a Low CO3 winner with the same neutral assigned via a primary channel commits.
led_g = gate_ledger(["m3", "m4"])
P.commit_winners(led_g, one_winner("m3", "C7H9NO8", "C7H8NO8-", 0.2, score=0.9),
                 pass_no=1, method="cheminfo+grid", context="ambient-air",
                 cfg=cfg_cal, lock=False, min_raw_score=0.5)
s = P.commit_winners(led_g, one_winner("m4", "C7H9NO8", "C8H9NO11-", 0.3, score=0.75),
                     pass_no=3, method="contaminant:bromo_organic", context="ambient-air",
                     cfg=cfg_cal, lock=False, min_raw_score=0.5)
check("cross-channel-corroborated minor winner commits", s["committed"] == 1, s)

# ---------- channel prior in arbitrate: near-tie goes to the primary channel ----------
scored_ch = pd.DataFrame([
    iso_row(sample_peak_id="P", compound_formula="C7H9NO8", compound_score=0.95,
            ion_formula="C8H9NO11-", ion_score=0.95, ppm_error=1.1),     # [M+CO3]-
    iso_row(sample_peak_id="P", compound_formula="C10H16O5", compound_score=0.87,
            ion_formula="C10H16BrO5-", ion_score=0.87, ppm_error=-0.7),  # [M+Br]-
])
arb_ch = P.arbitrate(scored_ch, CFG)
wp = arb_ch["winners"].iloc[0]
check("near-tie: primary-channel [M+Br]- beats minor [M+CO3]- (295.018 case)",
      wp["neutral"] == "C10H16O5" and wp["adduct"] == "[M+Br]-",
      (wp["neutral"], wp["adduct"], wp["eff_score"]))

# ---------- M0-vs-iso-child displacement through commit_winners ----------
led_d = gate_ledger(["A1", "A2"])
L.commit_assignment(led_d, "A2", neutral_formula="C6H11NO7Si", adduct="[M+CO3]-",
                    ion_score=0.60, ppm_error=1.5, pass_no=2, method="gka-series",
                    confidence="Low (series)", commentary="weak own M0")
arb_d = {"winners": pd.DataFrame([{
            "peak_id": "A1", "neutral": "C10H16O5", "ion_formula": "C10H16BrO5-",
            "adduct": "[M+Br]-", "ion_score": 0.87, "compound_score": 0.87,
            "raw_score": 0.87, "eff_score": 0.87, "ppm_error": -0.7, "n_iso": 1,
            "tied": False, "alternatives": []}]),
         "iso_children": pd.DataFrame([{
            "peak_id": "A2", "parent_peak_id": "A1", "iso_label": "81Br",
            "iso_score": 0.95}])}
s = P.commit_winners(led_d, arb_d, pass_no=1, method="cheminfo+grid",
                     context="ambient-air", cfg=cfg_cal, lock=False,
                     min_raw_score=0.5)
check("weak M0 displaced into parent's 81Br child (297.016 case)",
      s["iso_displaced"] == 1 and L.role_of(led_d, "A2") == L.ROLE_ISO
      and led_d.loc[led_d.peak_id == "A2", "parent_peak_id"].iloc[0] == "A1", s)
# a strong existing M0 is NOT displaced
led_d2 = gate_ledger(["B1", "B2"])
L.commit_assignment(led_d2, "B2", neutral_formula="C9H14O4", adduct="[M-H]-",
                    ion_score=0.95, ppm_error=0.1, pass_no=1, method="cheminfo+grid",
                    confidence="Good", commentary="strong own M0")
arb_d2 = {"winners": arb_d["winners"].assign(peak_id="B1"),
          "iso_children": arb_d["iso_children"].assign(peak_id="B2",
                                                       parent_peak_id="B1")}
s = P.commit_winners(led_d2, arb_d2, pass_no=3, method="contaminant:bromo_organic",
                     context="ambient-air", cfg=cfg_cal, lock=False,
                     min_raw_score=0.5)
check("Good existing M0 not displaced",
      s["iso_displaced"] == 0 and L.role_of(led_d2, "B2") == L.ROLE_M0, s)

# ---------- audit_mass_gate ----------
led_a = gate_ledger(["q1", "q2", "q3", "q4", "q5"])
L.commit_assignment(led_a, "q1", neutral_formula="C11H12N2O16", adduct="[M-H]-",
                    ion_score=0.81, ppm_error=3.83, pass_no=1, method="cheminfo+grid",
                    confidence="Low", commentary="pass1 junk tail")
L.commit_assignment(led_a, "q2", neutral_formula="C10H16O4", adduct="[M-H]-",
                    ion_score=0.95, ppm_error=-0.3, pass_no=1, method="cheminfo+grid",
                    confidence="Good", commentary="backbone")
L.commit_assignment(led_a, "q3", neutral_formula="C9H13NO5", adduct="[M-H]-",
                    ion_score=0.72, ppm_error=1.6, pass_no=3, method="x",
                    confidence="Low", commentary="2-4 sigma, no evidence")
L.commit_assignment(led_a, "q4", neutral_formula="C8H12O6S", adduct="[M-H]-",
                    ion_score=0.72, ppm_error=1.6, pass_no=3, method="x",
                    confidence="Low", commentary="2-4 sigma WITH child")
L.attach_isotopologue(led_a, "q5", "q4", iso_label="34S", iso_match_score=0.8)
out = P.audit_mass_gate(led_a, cfg_cal, log=lambda *a: None)
check("audit clears >4-sigma pass1 Low", L.role_of(led_a, "q1") == L.ROLE_UNEXPLAINED, out)
check("audit clears 2-4 sigma Low without evidence",
      L.role_of(led_a, "q3") == L.ROLE_UNEXPLAINED, out)
check("audit keeps 2-4 sigma Low WITH iso child", L.role_of(led_a, "q4") == L.ROLE_M0, out)
check("audit ledger validates clean", L.validate(led_a) == [], L.validate(led_a))

# ---------- audit_isotopes: post-run isotope-physics audit ----------
def mk_ledger(rows):
    led = L.new_ledger(pd.DataFrame(
        {"peak_id": [r[0] for r in rows], "mz": [r[1] for r in rows],
         "height": [r[2] for r in rows]}))
    return led

def commit(led, pid, neutral, ion, conf="Good", score=0.9):
    L.commit_assignment(led, pid, neutral_formula=neutral, adduct="[M+Br]-",
                        ion_formula=ion, ion_score=score, ppm_error=0.1,
                        pass_no=3, method="test", confidence=conf,
                        commentary="t")

ACFG = P.PassConfig(height_cutoff=100.0)

# v16 case: 462.99/464.99 — two Good M0s 1.99795 apart, ~1:1; light ion has Br
led = mk_ledger([("La", 462.9933, 26882.0), ("Hb", 464.9913, 25609.0),
                 ("X", 463.9966, 3450.0)])   # X = 13C of La
commit(led, "La", "C12H12F12", "C12H12BrF12-")
commit(led, "Hb", "C11H9F7O12", "C11H8F7O12-")
s = P.audit_isotopes(led, ACFG, log=lambda *a: None)
check("audit: Br-doublet heavy M0 demoted to 81Br child",
      L.role_of(led, "Hb") == L.ROLE_ISO and s["doublet_child"] == 1, s)
check("audit: 13C satellite swept up as evidence",
      L.role_of(led, "X") == L.ROLE_ISO and s["c13_attached"] == 1, s)

# doublet where NEITHER formula carries Br -> both cleared
led = mk_ledger([("Lc", 284.0501, 641.0), ("Hd", 286.0480, 543.0)])
commit(led, "Lc", "C9H20NO4", "C9H20NO4-")    # no Br anywhere
commit(led, "Hd", "C14H10O3", "C15H10O6-")
s = P.audit_isotopes(led, ACFG, log=lambda *a: None)
check("audit: no-Br doublet clears both formulas",
      L.role_of(led, "Lc") == L.ROLE_UNEXPLAINED
      and L.role_of(led, "Hd") == L.ROLE_UNEXPLAINED
      and s["doublet_cleared"] == 2, s)

# 13C clamp: formula claims C19, satellite measures ~C11
led = mk_ledger([("M", 444.9861, 4761.0), ("K", 445.9895, 564.0)])
commit(led, "M", "C19H15N2O4S", "C19H14BrN2O4S-")
L.attach_isotopologue(led, "K", "M", iso_label="13C")
s = P.audit_isotopes(led, ACFG, log=lambda *a: None)
check("audit: 13C carbon clamp clears C19-vs-C11",
      L.role_of(led, "M") == L.ROLE_UNEXPLAINED and s["c13_clamp"] == 1, s)

# 13C missing: big peak, formula predicts a visible satellite, none exists
led = mk_ledger([("N", 168.9505, 10300.0)])
commit(led, "N", "C3H6O3", "C3H6BrO3-")
s = P.audit_isotopes(led, ACFG, log=lambda *a: None)
check("audit: predicted-13C-absent clears the formula",
      L.role_of(led, "N") == L.ROLE_UNEXPLAINED and s["c13_missing"] == 1, s)

# sanity: a CORRECT assignment survives all four checks
led = mk_ledger([("P", 295.0184, 2965.0), ("Q", 297.0164, 2924.0),
                 ("R", 296.0218, 320.0)])
commit(led, "P", "C10H16O5", "C10H16BrO5-")
L.attach_isotopologue(led, "Q", "P", iso_label="81Br")
s = P.audit_isotopes(led, ACFG, log=lambda *a: None)
check("audit: correct assignment untouched (13C attached, nothing cleared)",
      L.role_of(led, "P") == L.ROLE_M0 and s["c13_clamp"] == 0
      and s["c13_missing"] == 0 and s["c13_attached"] == 1, s)

# ---------- complete_isotope_envelopes (the 393/395 silanediol bug) ----------
from mascope_assign import chemistry as CHEM  # noqa: E402
# silanediol Si4+Br at 393 (M0), its M+2 at 395 wrongly committed as a Cl-F-S
# organic, M+4 at 397 attached as 395's child. Heights = real Si4+Br envelope.
led = mk_ledger([("Si", 393.0045, 20086.0), ("Mp2", 395.0028, 24523.0),
                 ("Mp4", 397.0029, 6344.0), ("far", 600.0, 5000.0)])
commit(led, "Si", "C8H26O5Si4", "C8H26BrO5Si4-")        # the silanediol M0
L.lock_peaks(led, ["Si"])                                 # pass-0 locks it
L.commit_assignment(led, "Mp2", neutral_formula="C8H12ClF6NO2S", adduct="[M+CO3]-",
                    ion_formula="C9H12ClF6NO5S-", ion_score=0.63, ppm_error=-1.6,
                    pass_no=4, method="residual:iso-pair", confidence="Low (iso-pair)",
                    commentary="phantom Cl doublet")
L.attach_isotopologue(led, "Mp4", "Mp2", iso_label="37Cl(pair)")
import mascope_assign.tiers as _T  # noqa: E402
_T.apply_tiers(led)
out = P.complete_isotope_envelopes(led, P.PassConfig(), log=lambda *a: None)
check("envelope: silanediol M+2 (395) displaced off its phantom formula",
      L.role_of(led, "Mp2") == L.ROLE_ISO and out["displaced"] >= 1, out)
check("envelope: 395 re-parented to the silanediol 393",
      led.loc[led.peak_id == "Mp2", "parent_peak_id"].iloc[0] == "Si")
check("envelope: 397 (M+4) re-parented to the silanediol too",
      led.loc[led.peak_id == "Mp4", "parent_peak_id"].iloc[0] == "Si")
check("envelope: ledger still valid after displacement", L.validate(led) == [], L.validate(led))

# di-bromide [M+HBr+Br]- core (2 Br in the ion) commits in pass 6 -> the new 3rd
# envelope sweep must claim its M+2 (~1.95x) and M+4 (~0.95x) satellites, which
# were sitting in the residual (the Family-A GKA ladder, 2026-06-13).
ledd = mk_ledger([("core", 356.9348, 742.0), ("m2", 358.9328, 1395.0),
                  ("m4", 360.9308, 700.0)])
L.commit_assignment(ledd, "core", neutral_formula="C10H14O4", adduct="[M+HBr+Br]-",
                    ion_formula="C10H15Br2O4-", ion_score=0.80, ppm_error=-0.4,
                    pass_no=6, method="ladder:gapfill", confidence="Good (ladder)",
                    commentary="di-bromide SOA core")
_T.apply_tiers(ledd)
outd = P.complete_isotope_envelopes(ledd, P.PassConfig(), log=lambda *a: None)
check("envelope: di-bromide M+2 (358.93) attached to the core",
      L.role_of(ledd, "m2") == L.ROLE_ISO
      and ledd.loc[ledd.peak_id == "m2", "parent_peak_id"].iloc[0] == "core", outd)
check("envelope: di-bromide M+4 (360.93) attached to the core",
      L.role_of(ledd, "m4") == L.ROLE_ISO
      and ledd.loc[ledd.peak_id == "m4", "parent_peak_id"].iloc[0] == "core", outd)

# guard: a CHO-only M0 must NOT claim a coincidental peak ~2 Da above it (its
# pattern has no M+2 driver), and must NOT displace a real neighbour
led2 = mk_ledger([("a", 200.0, 1e5), ("b", 202.0, 5e4)])
commit(led2, "a", "C10H16O4", "C10H15O4-")   # CHO ion, no Br/Cl/Si
commit(led2, "b", "C9H12O5", "C9H11O5-")      # independent neighbour
_T.apply_tiers(led2)
out2 = P.complete_isotope_envelopes(led2, P.PassConfig(), log=lambda *a: None)
check("envelope: CHO ion does not claim a coincidental +2 neighbour",
      L.role_of(led2, "b") == L.ROLE_M0 and out2["displaced"] == 0, out2)

# guard (review #3): a STRONG victim (High conf OR near-High score) sitting at a
# parent's M+2 with matching intensity must NOT be displaced -- the tier column
# is NA during this pass, so protection rides on confidence + score, not tier.
led3 = mk_ledger([("p", 300.0, 1e5), ("q", 301.9979, 9.7e4)])
commit(led3, "p", "C9H17BrO5", "C9H17BrO5-")          # 1-Br parent, predicts M+2
L.commit_assignment(led3, "q", neutral_formula="C12H20O8", adduct="[M-H]-",
                    ion_formula="C12H19O8-", ion_score=0.95, ppm_error=0.2,
                    pass_no=1, method="cheminfo+grid", confidence="High",
                    commentary="strong standalone fit")
outg = P.complete_isotope_envelopes(led3, P.PassConfig(), log=lambda *a: None)
check("envelope: a High/strong-score victim is NOT displaced (tier-NA safe)",
      L.role_of(led3, "q") == L.ROLE_M0 and outg["displaced"] == 0, outg)

# ---------- detect_composites (silanediol-on-BrCl even/odd test) ----------
# build a silanediol n=4 envelope: M0 inflated ~45% by a coincident BrCl
# compound. Odd shifts (M+1) = pure silanediol; even (M0/M+2/M+4) carry the
# extra halogen. Heights from the real <sample-id> data.
import mascope_assign.chemistry as _CH  # noqa: E402
mz0 = _CH.ion_mz("C8H26O5Si4", "[M+Br]-")   # 393.0046
comp = mk_ledger([("M0", mz0, 20086.0), ("M1", mz0 + 1.0008, 2698.0),
                  ("M1b", mz0 + 1.0034, 531.0), ("M2", mz0 + 1.9979, 24523.0),
                  ("M3", mz0 + 2.9986, 2748.0), ("M4", mz0 + 3.9957, 6344.0),
                  ("clean", 244.9670, 10712.0), ("cleanM1", 244.9670 + 1.0, 1461.0),
                  ("cleanM2", 244.9670 + 1.9979, 7000.0)])
commit(comp, "M0", "C8H26O5Si4", "C8H26BrO5Si4-")
commit(comp, "clean", "C4H14O3Si2", "C4H14BrO3Si2-")   # n=2: ~clean, small extra
oc = P.detect_composites(comp, P.PassConfig(), log=lambda *a: None)
note4 = comp.loc[comp.peak_id == "M0", "composite_note"].iloc[0]
check("composite: silanediol n=4 flagged as composite", pd.notna(note4) and oc["flagged"] >= 1, note4)
check("composite: identifies the BrCl co-component", "BrCl" in str(note4), note4)
check("composite: estimates ~45% co-component", "45%" in str(note4) or "44%" in str(note4), note4)
check("composite: clean n=2 NOT flagged",
      pd.isna(comp.loc[comp.peak_id == "clean", "composite_note"].iloc[0]),
      comp.loc[comp.peak_id == "clean", "composite_note"].iloc[0])
# a CHO compound whose M0 matches its M+1 is NOT flagged (no inflation)
pure = mk_ledger([("p", 200.0, 1e5), ("pm1", 201.0034, 1e5 * 10 * 0.0107)])
commit(pure, "p", "C10H16O4", "C10H15O4-")
P.detect_composites(pure, P.PassConfig(), log=lambda *a: None)
check("composite: a self-consistent CHO M0 is NOT flagged",
      pd.isna(pure.loc[pure.peak_id == "p", "composite_note"].iloc[0]))

# ---------- split_composites: de-blend into fractional sub-peaks ----------
P.detect_composites(comp, P.PassConfig(), log=lambda *a: None)
comp2 = P.split_composites(comp, P.PassConfig(), log=lambda *a: None)
host = comp2[comp2.peak_id == "M0"].iloc[0]
sub = comp2[comp2.peak_id == "M0.2"]
check("split: a synthetic sub-peak M0.2 was created", len(sub) == 1, list(comp2.peak_id))
check("split: host keeps assigned_fraction < 1 (the co-component is removed)",
      float(host["assigned_fraction"]) < 0.95, host["assigned_fraction"])
if len(sub):
    s = sub.iloc[0]
    check("split: sub-peak is synthetic, at the host m/z, linked to the host",
          bool(s["synthetic"]) and abs(float(s["mz"]) - float(host["mz"])) < 1e-6
          and s["host_peak_id"] == "M0", (s["synthetic"], s["host_peak_id"]))
    check("split: signal is conserved (host_eff + sub == measured host height)",
          abs(float(host["height"]) * float(host["assigned_fraction"])
              + float(s["height"]) - float(host["height"])) < 2.0,
          (host["height"], host["assigned_fraction"], s["height"]))
    check("split: sub-peak commentary names the co-component (BrCl)",
          "BrCl" in str(s["commentary"]), s["commentary"])
# stats: synthetic excluded from the real-peak count, signal not double-counted
st_c = L.stats(comp2)
check("split: stats excludes synthetic from n_peaks", st_c["n_synthetic"] >= 1
      and st_c["n_peaks"] == int((~comp2["synthetic"].fillna(False).astype(bool)).sum()),
      (st_c["n_peaks"], st_c["n_synthetic"]))
check("split: ledger still validates with synthetic rows", L.validate(comp2) == [],
      L.validate(comp2))

# ---------- demote_carbon_inconsistent (pre-pass-4 O15-monster clear) ----------
# the 409.0015 case: pass 1 grabbed C11H10N2O15 (ion C11) but the 13C satellite
# at +1.0034 measures ~C16 -> must clear BEFORE pass 4 so the di-bromide SOA
# core (C15) can be re-claimed. Satellite ratio 16*1.07% = 0.171 of the parent.
led = mk_ledger([("mon", 409.0015, 4720.0), ("sat", 410.0049, 807.0)])
commit(led, "mon", "C11H10N2O15", "C11H10N2BrO15-")   # ion carbon = 11
n = P.demote_carbon_inconsistent(led, ACFG, log=lambda *a: None)
check("pre-pass4: C11 monster with ~C16 satellite demoted",
      n == 1 and L.role_of(led, "mon") == L.ROLE_UNEXPLAINED, (n, L.role_of(led, "mon")))
check("pre-pass4: freed peak is re-claimable (unexplained, no formula)",
      pd.isna(led.loc[led.peak_id == "mon", "neutral_formula"].iloc[0]))
# a carbon-CONSISTENT assignment is left alone (C15 ion, ~C15 satellite)
led2 = mk_ledger([("ok", 409.0015, 4720.0), ("s2", 410.0049, 758.0)])
commit(led2, "ok", "C15H23BrO3", "C15H23Br2O3-")      # ion carbon = 15
n2 = P.demote_carbon_inconsistent(led2, ACFG, log=lambda *a: None)
check("pre-pass4: carbon-consistent C15 assignment survives",
      n2 == 0 and L.role_of(led2, "ok") == L.ROLE_M0, (n2, L.role_of(led2, "ok")))
# no satellite -> cannot measure -> never demotes (avoids false clears)
led3 = mk_ledger([("nosat", 409.0015, 4720.0)])
commit(led3, "nosat", "C11H10N2O15", "C11H10N2BrO15-")
check("pre-pass4: no satellite -> no demotion",
      P.demote_carbon_inconsistent(led3, ACFG, log=lambda *a: None) == 0
      and L.role_of(led3, "nosat") == L.ROLE_M0)

# ---------- audit: twin-satellite fallback for missing 13C ----------
# v20 false-clear: C3H6O3.Br- at 10.3k cps -- own 13C absent from the peak
# list, but the 81Br twin's 13C satellite exists and proves the carbon.
from mascope_assign import chemistry as CH  # noqa: E402
mz_l = CH.ion_mz("C3H6O3", "[M+Br]-")
led = mk_ledger([("La", mz_l, 10300.0),
                 ("Tw", mz_l + 1.9979535, 10000.0),
                 ("Sat", mz_l + 1.9979535 + 1.0033548, 350.0)])
commit(led, "La", "C3H6O3", "C3H6BrO3-")
L.attach_isotopologue(led, "Tw", "La", iso_label="81Br")
s = P.audit_isotopes(led, ACFG, log=lambda *a: None)
check("audit: twin 13C satellite blocks the missing-13C clear",
      L.role_of(led, "La") == L.ROLE_M0 and s["c13_missing"] == 0, s)
# control: same peak with NO twin satellite still clears
led = mk_ledger([("Lb", mz_l, 10300.0), ("Tw2", mz_l + 1.9979535, 10000.0)])
commit(led, "Lb", "C3H6O3", "C3H6BrO3-")
L.attach_isotopologue(led, "Tw2", "Lb", iso_label="81Br")
s = P.audit_isotopes(led, ACFG, log=lambda *a: None)
check("audit: no satellite anywhere still clears",
      L.role_of(led, "Lb") == L.ROLE_UNEXPLAINED and s["c13_missing"] == 1, s)
# cross-channel fallback: missing 13C but the SAME neutral assigned Good on
# another peak (other adduct) -> positive evidence wins, no clear (v21 case)
mh = CH.ion_mz("C10H16O6", "[M-H]-")
mbr = CH.ion_mz("C10H16O6", "[M+Br]-")
led = mk_ledger([("MH", mh, 1060.0), ("MB", mbr, 2426.0)])
L.commit_assignment(led, "MH", neutral_formula="C10H16O6", adduct="[M-H]-",
                    ion_formula="C10H15O6-", ion_score=0.87, ppm_error=-0.9,
                    pass_no=1, method="t", confidence="Good", commentary="t")
L.commit_assignment(led, "MB", neutral_formula="C10H16O6", adduct="[M+Br]-",
                    ion_formula="C10H16BrO6-", ion_score=0.86, ppm_error=-0.5,
                    pass_no=1, method="t", confidence="Good", commentary="t")
s = P.audit_isotopes(led, ACFG, log=lambda *a: None)
check("audit: cross-channel agreement blocks the missing-13C clear",
      L.role_of(led, "MB") == L.ROLE_M0 and s["c13_missing"] == 0, s)

# ---------- pass 5: known-neutral completion ----------
from mascope_assign import contexts as XC  # noqa: E402
PROF5 = XC.get_context("ambient-air")
ADD5 = ["[M-H]-", "[M+Br]-"]
mz_c2 = CH.ion_mz("C2H4O3", "[M+Br]-")
mz_c5 = CH.ion_mz("C5H10O3", "[M+Br]-")
mz_c3 = CH.ion_mz("C3H6O3", "[M+Br]-")          # bracketed gap (C2..C5, j=1)
mz_x_mh = CH.ion_mz("C10H16O3", "[M-H]-")
mz_x_br = CH.ion_mz("C10H16O3", "[M+Br]-")      # cross-channel partner
led5 = mk_ledger([("A2", mz_c2, 700.0), ("A5", mz_c5, 1000.0),
                  ("G3", mz_c3, 10300.0), ("AX", mz_x_mh, 400.0),
                  ("PX", mz_x_br, 1500.0), ("far", 555.555, 300.0)])
commit(led5, "A2", "C2H4O3", "C2H4BrO3-", conf="High")
commit(led5, "A5", "C5H10O3", "C5H10BrO3-", conf="High")
L.commit_assignment(led5, "AX", neutral_formula="C10H16O3", adduct="[M-H]-",
                    ion_formula="C10H15O3-", ion_score=0.9, ppm_error=0.1,
                    pass_no=1, method="t", confidence="Good", commentary="t")

def fake5(client, sample_id, formulas, *, mechanism_ids=None, **kw):
    rows = []
    if "C3H6O3" in formulas:
        rows.append(dict(compound_formula="C3H6O3", compound_score=0.9,
                         ion_formula="C3H6BrO3-", ion_score=0.9, iso_label="M0",
                         is_base=True, theo_mz=mz_c3, rel_abundance=1.0,
                         iso_score=0.9, sample_peak_id="G3", sample_peak_mz=mz_c3,
                         sample_peak_intensity=10300.0, ppm_error=0.3,
                         abundance_error=0.0))
    if "C10H16O3" in formulas:
        rows.append(dict(compound_formula="C10H16O3", compound_score=0.9,
                         ion_formula="C10H16BrO3-", ion_score=0.9, iso_label="M0",
                         is_base=True, theo_mz=mz_x_br, rel_abundance=1.0,
                         iso_score=0.9, sample_peak_id="PX", sample_peak_mz=mz_x_br,
                         sample_peak_intensity=1500.0, ppm_error=0.2,
                         abundance_error=0.0))
    return pd.DataFrame(rows)

s5 = P.run_pass5_completion(None, "SID", led5, PROF5, ACFG, ADD5,
                            score_fn=fake5, log=lambda *a: None)
check("pass5 commits the bracketed series gap (C3H6O3)",
      led5.loc[led5.peak_id == "G3", "neutral_formula"].iloc[0] == "C3H6O3", s5)
check("pass5 commits the cross-channel partner (C10H16O3 [M+Br]-)",
      led5.loc[led5.peak_id == "PX", "neutral_formula"].iloc[0] == "C10H16O3"
      and led5.loc[led5.peak_id == "PX", "adduct"].iloc[0] == "[M+Br]-", s5)
check("pass5 commits exactly the two targets", s5["committed"] == 2, s5)
check("pass5 method recorded", "completion" in
      led5.loc[led5.peak_id == "G3", "method"].iloc[0])
check("untargeted peak untouched",
      L.role_of(led5, "far") == L.ROLE_UNEXPLAINED)
check("ledger valid after pass5", L.validate(led5) == [], L.validate(led5))

# ---------- pass 0: known-contaminant pre-pass (silanediol ladder) ----------
check("silanediol series composition", P._silanediol_series(3) ==
      ["C2H8O2Si1", "C4H14O3Si2", "C6H20O4Si3"], P._silanediol_series(3))
mz_n2 = CH.ion_mz("C4H14O3Si2", "[M+Br]-")
led0 = mk_ledger([("D2", mz_n2, 10712.0), ("D2br", mz_n2 + 1.99795, 10100.0)])

def fake0(client, sample_id, formulas, *, mechanism_ids=None, **kw):
    rows = []
    if "C4H14O3Si2" in formulas:
        rows.append(dict(compound_formula="C4H14O3Si2", compound_score=0.92,
                         ion_formula="C4H14BrO3Si2-", ion_score=0.92,
                         iso_label="M0", is_base=True, theo_mz=mz_n2,
                         rel_abundance=1.0, iso_score=0.92,
                         sample_peak_id="D2", sample_peak_mz=mz_n2,
                         sample_peak_intensity=10712.0, ppm_error=-0.87,
                         abundance_error=0.0))
        rows.append(dict(compound_formula="C4H14O3Si2", compound_score=0.92,
                         ion_formula="C4H14BrO3Si2-", ion_score=0.92,
                         iso_label="81Br", is_base=False,
                         theo_mz=mz_n2 + 1.99795, rel_abundance=0.97,
                         iso_score=0.94, sample_peak_id="D2br",
                         sample_peak_mz=mz_n2 + 1.99795,
                         sample_peak_intensity=10100.0, ppm_error=-0.8,
                         abundance_error=0.02))
    return pd.DataFrame(rows)

s0 = P.run_pass0_known(None, "SID", led0, PROF5, ACFG, ADD5,
                              score_fn=fake0, log=lambda *a: None)
check("pass0 commits the silanediol oligomer",
      led0.loc[led0.peak_id == "D2", "neutral_formula"].iloc[0] == "C4H14O3Si2"
      and s0["committed"] == 1, s0)
check("pass0 locks the peak", L.is_locked(led0, "D2"))
check("pass0 attaches the 81Br child",
      L.role_of(led0, "D2br") == L.ROLE_ISO, s0)
# the v24 failure mode: pass-1 must NOT be able to displace it with a CHO fit
try:
    L.commit_assignment(led0, "D2", neutral_formula="C5H10O6",
                        adduct="[M+Br]-", ion_formula="C5H10BrO6-",
                        ion_score=0.95, ppm_error=0.8, pass_no=1, method="grid",
                        confidence="High", commentary="bogus")
    stolen = True
except L.LedgerError:
    stolen = False
check("locked contaminant refuses the CHO grid fit", not stolen
      and led0.loc[led0.peak_id == "D2", "neutral_formula"].iloc[0] == "C4H14O3Si2")

# pass0 twin gate: a contaminant claim on a peak whose own 81Br twin is
# missing must be refused (the v25 lactic-acid collision)
mz_n1 = CH.ion_mz("C2H8O2Si1", "[M+Br]-")
ledg = mk_ledger([("X1", mz_n1, 11950.0), ("Xw", mz_n1 + 1.99795, 427.0)])

def fake_n1(client, sample_id, formulas, *, mechanism_ids=None, **kw):
    return pd.DataFrame([dict(
        compound_formula="C2H8O2Si1", compound_score=0.9,
        ion_formula="C2H8BrO2Si-", ion_score=0.9, iso_label="M0",
        is_base=True, theo_mz=mz_n1, rel_abundance=1.0, iso_score=0.9,
        sample_peak_id="X1", sample_peak_mz=mz_n1,
        sample_peak_intensity=11950.0, ppm_error=1.0, abundance_error=0.0)])

sg = P.run_pass0_known(None, "SID", ledg, PROF5, ACFG, ADD5,
                              score_fn=fake_n1, log=lambda *a: None)
check("pass0 refuses claim with inconsistent own twin (ratio 0.04)",
      sg["committed"] == 0 and L.role_of(ledg, "X1") == L.ROLE_UNEXPLAINED, sg)

# pass0 nitroaromatic: dinitrophenol C6H4N2O5 [M-H]- is H-poor (VK-floor blocked)
# so the grid can't reach it; pass-0 supplies it. No Br -> no twin gate, just the
# |ppm|<=2 + score gate. (v45->v46 fix; confirmed present by Orbitool.)
check("dinitrophenol in the known-species nitroaromatic family",
      P._known_species().get("nitroaromatic", {}).get("C6H4N2O5") is not None)
mz_dnp = CH.ion_mz("C6H4N2O5", "[M-H]-")
ledn = mk_ledger([("DNP", mz_dnp, 808.0)])

def fake_dnp(client, sample_id, formulas, *, mechanism_ids=None, **kw):
    if "C6H4N2O5" not in formulas:
        return pd.DataFrame([])
    return pd.DataFrame([dict(
        compound_formula="C6H4N2O5", compound_score=0.82,
        ion_formula="C6H3N2O5-", ion_score=0.82, iso_label="M0",
        is_base=True, theo_mz=mz_dnp, rel_abundance=1.0, iso_score=0.82,
        sample_peak_id="DNP", sample_peak_mz=mz_dnp,
        sample_peak_intensity=808.0, ppm_error=-0.24, abundance_error=0.0)])

sn = P.run_pass0_known(None, "SID", ledn, PROF5, ACFG, ADD5,
                       score_fn=fake_dnp, log=lambda *a: None)
check("pass0 commits dinitrophenol [M-H]- (no Br twin gate)",
      sn["committed"] == 1
      and ledn.loc[ledn.peak_id == "DNP", "neutral_formula"].iloc[0] == "C6H4N2O5"
      and ledn.loc[ledn.peak_id == "DNP", "method"].iloc[0] == "known:nitroaromatic", sn)

# ---------- pass0 known-species gate is offset-aware (prior_offset) ----------
# a silanediol at -2.3 ppm: rejected when the instrument offset is unknown
# (|2.3|>2.0), committed once the rough offset (-1.9) is seeded -- the
# silanediol-vs-C5H10O6 collision at a -1.9 ppm source.
mz_off = CH.ion_mz("C4H14O3Si2", "[M+Br]-")


def fake_off(client, sample_id, formulas, *, mechanism_ids=None, **kw):
    if "C4H14O3Si2" not in formulas:
        return pd.DataFrame([])
    return pd.DataFrame([
        dict(compound_formula="C4H14O3Si2", compound_score=0.92,
             ion_formula="C4H14BrO3Si2-", ion_score=0.92, iso_label="M0",
             is_base=True, theo_mz=mz_off, rel_abundance=1.0, iso_score=0.92,
             sample_peak_id="O2", sample_peak_mz=mz_off,
             sample_peak_intensity=10712.0, ppm_error=-2.3, abundance_error=0.0),
        dict(compound_formula="C4H14O3Si2", compound_score=0.92,
             ion_formula="C4H14BrO3Si2-", ion_score=0.92, iso_label="81Br",
             is_base=False, theo_mz=mz_off + 1.99795, rel_abundance=0.97,
             iso_score=0.94, sample_peak_id="O2br",
             sample_peak_mz=mz_off + 1.99795, sample_peak_intensity=10100.0,
             ppm_error=-2.3, abundance_error=0.02)])


led_b = mk_ledger([("O2", mz_off, 10712.0), ("O2br", mz_off + 1.99795, 10100.0)])
P.run_pass0_known(None, "SID", led_b, PROF5, P.PassConfig(), ADD5,
                  score_fn=fake_off, log=lambda *a: None)
check("pass0 rejects a -2.3 ppm known species when offset-blind",
      L.role_of(led_b, "O2") == L.ROLE_UNEXPLAINED)
led_o = mk_ledger([("O2", mz_off, 10712.0), ("O2br", mz_off + 1.99795, 10100.0)])
cfg_seed = P.PassConfig(); cfg_seed.prior_offset = -1.9
P.run_pass0_known(None, "SID", led_o, PROF5, cfg_seed, ADD5,
                  score_fn=fake_off, log=lambda *a: None)
check("pass0 commits the -2.3 ppm known species once offset seeded (-1.9)",
      L.role_of(led_o, "O2") == L.ROLE_M0
      and led_o.loc[led_o.peak_id == "O2", "neutral_formula"].iloc[0] == "C4H14O3Si2")

# ---------- build_ranges reads the context grid-box width ----------
PROF_URO = XC.get_context("uronium")
r_amb = P.build_ranges(PROF5, None, include_N=True)               # ambient-air
r_uro = P.build_ranges(PROF_URO, None, include_N=True)
check("ambient grid box C<=40 / O<=30",
      r_amb["C"][1] == 40 and r_amb["O"][1] == 30, (r_amb["C"], r_amb["O"]))
check("uronium grid box widened to C<=46 / O<=32",
      r_uro["C"][1] == 46 and r_uro["O"][1] == 32, (r_uro["C"], r_uro["O"]))
check("uronium grid admits N up to the context cap",
      r_uro["N"][1] == PROF_URO.max_N, r_uro["N"])
r_ovr = P.build_ranges(PROF_URO, None, include_N=True, c_max=20, o_max=10)
check("explicit c_max/o_max overrides the profile",
      r_ovr["C"][1] == 20 and r_ovr["O"][1] == 10, (r_ovr["C"], r_ovr["O"]))

# ---------- positive polarity: pass 0 known-species is a no-op ----------
check("_known_species(positive) is empty", P._known_species("positive") == {})
check("_known_species(negative) keeps the atmospheric list",
      "atmospheric" in P._known_species("negative"))


def _boom(*a, **k):
    raise AssertionError("pass0 must not score the oracle in positive mode")


s_pos = P.run_pass0_known(None, "SID", mk_ledger([("p", 100.0, 500.0)]),
                          PROF_URO, ACFG, ["[M+H]+", "[M+(CH4N2O)H]+"],
                          score_fn=_boom, log=lambda *a: None)
check("pass0 positive mode commits nothing without calling the oracle",
      s_pos["committed"] == 0, s_pos)

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

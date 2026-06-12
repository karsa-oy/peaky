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

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

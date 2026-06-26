"""Offline tests for passes.rearbitrate_offcal_degenerate -- the off-calibration
degenerate-winner swap that applies the tier engine's calibration-sigma +
corroboration gate AT WINNER-SELECTION. Run: python3 tests/test_rearbitrate.py

Reproduces the orange-uronium case the fix targets: m/z 464.143 committed the
mass-degenerate C36H17N [M+H]+ (DBE 29 azabenzo-PAH, ~-4.8 sigma off the
calibrated mass center, no isotope/cross-channel corroboration) over the on-cal
C22H23N3O7 [M+Na]+ that sat in its alternatives list. The winner-selection bug:
pass 1 commits before calibration, so the off-cal gate never sees the monster.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky import ledger as L  # noqa: E402
from peaky import passes  # noqa: E402

PassConfig = passes.PassConfig

PASS = FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


def build_ledger():
    """A ledger with a 20-peak isotopologue-backed CHO/CHON calibration backbone
    (centered near 0 ppm) plus four probe peaks. tiers._calibrate needs >=20 core
    peaks, so the backbone makes the off-cal gate live."""
    rows = []
    # --- backbone: parent + 13C child, on-cal, High, CHO ---
    ppm_cycle = [-0.10, 0.0, 0.10, 0.05, -0.05]
    for i in range(20):
        pid, cid = f"bb{i}", f"bb{i}c"
        rows.append((pid, 150.0 + i, 1e5))
        rows.append((cid, 151.0 + i, 1e3))
    # --- probe peaks ---
    rows += [
        ("P_swap", 464.143, 5e3),       # off-cal monster, on-cal plausible alt -> SWAP
        ("P_badalt", 470.10, 5e3),      # off-cal monster, only on-cal alt is implausible -> KEEP
        ("P_corrob", 480.20, 5e3),      # off-cal monster but isotopologue-backed -> KEEP
        ("P_corrob_c", 481.20, 1.6e3),  # the 13C child of P_corrob
        ("P_oncal", 490.30, 5e3),       # on-cal winner -> KEEP
    ]
    peaks = pd.DataFrame(rows, columns=["peak_id", "mz", "height"])
    led = L.new_ledger(peaks)

    for i in range(20):
        pid, cid = f"bb{i}", f"bb{i}c"
        nf = f"C{8 + i}H{14 + i}O4"
        L.commit_assignment(led, pid, neutral_formula=nf, adduct="[M-H]-",
                            ion_formula=nf, ion_score=0.95, compound_score=0.95,
                            eff_score=0.93, eff_margin=0.2, tied=False,
                            ppm_error=ppm_cycle[i % 5], pass_no=1,
                            method="cheminfo+grid", confidence="High",
                            commentary="backbone",
                            isotopologues=[{"label": "13C", "peak_id": cid}])
        L.attach_isotopologue(led, cid, pid, iso_label="13C", iso_match_score=0.9)

    # P_swap: the case-2 reproduction. Winner C36H17N (DBE 29, DBE/C 0.81) off-cal
    # at -1.0 ppm (~-6.7 sigma); alt C22H23N3O7 (DBE 13) on-cal at +0.15 ppm.
    L.commit_assignment(led, "P_swap", neutral_formula="C36H17N", adduct="[M+H]+",
                        ion_formula="C36H18N+", ion_score=0.88, compound_score=0.88,
                        eff_score=0.85, eff_margin=0.24, tied=False, ppm_error=-1.0,
                        pass_no=1, method="cheminfo+grid", confidence="Good",
                        commentary="Pass 1",
                        alternatives=[{"formula": "C22H23N3O7", "adduct": "[M+Na]+",
                                       "ion_score": 0.70, "raw_score": 0.70,
                                       "eff_score": 0.61, "ppm": 0.15}])
    # P_badalt: off-cal C30H14 monster (DBE 24); its only on-cal alternative is a
    # carbon-cluster (DBE/C >= 1.0 -> implausible) -> must NOT swap.
    L.commit_assignment(led, "P_badalt", neutral_formula="C30H14", adduct="[M+H]+",
                        ion_formula="C30H15+", ion_score=0.88, compound_score=0.88,
                        eff_score=0.85, eff_margin=0.2, tied=False, ppm_error=-1.0,
                        pass_no=1, method="cheminfo+grid", confidence="Good",
                        commentary="Pass 1",
                        alternatives=[{"formula": "C24H2", "adduct": "[M+H]+",
                                       "ion_score": 0.80, "raw_score": 0.80,
                                       "eff_score": 0.70, "ppm": 0.10}])
    # P_corrob: same off-cal monster shape, but isotopologue-backed -> corroborated
    # -> must NOT swap (corroboration is the evidence that breaks a degeneracy).
    L.commit_assignment(led, "P_corrob", neutral_formula="C36H17N", adduct="[M+H]+",
                        ion_formula="C36H18N+", ion_score=0.88, compound_score=0.88,
                        eff_score=0.85, eff_margin=0.2, tied=False, ppm_error=-1.0,
                        pass_no=1, method="cheminfo+grid", confidence="Good",
                        commentary="Pass 1",
                        isotopologues=[{"label": "13C", "peak_id": "P_corrob_c"}],
                        alternatives=[{"formula": "C22H23N3O7", "adduct": "[M+Na]+",
                                       "ion_score": 0.70, "raw_score": 0.70,
                                       "eff_score": 0.61, "ppm": 0.15}])
    L.attach_isotopologue(led, "P_corrob_c", "P_corrob", iso_label="13C",
                          iso_match_score=0.9)
    # P_oncal: aromatic winner but ON calibration (+0.1 ppm) -> the committed
    # reading stands; never displaced even with a plausible alternative present.
    L.commit_assignment(led, "P_oncal", neutral_formula="C36H17N", adduct="[M+H]+",
                        ion_formula="C36H18N+", ion_score=0.88, compound_score=0.88,
                        eff_score=0.85, eff_margin=0.2, tied=False, ppm_error=0.1,
                        pass_no=1, method="cheminfo+grid", confidence="Good",
                        commentary="Pass 1",
                        alternatives=[{"formula": "C22H23N3O7", "adduct": "[M+Na]+",
                                       "ion_score": 0.70, "raw_score": 0.70,
                                       "eff_score": 0.61, "ppm": 0.15}])
    return led


def _nf(led, pid):
    return str(led.loc[led["peak_id"] == pid, "neutral_formula"].iloc[0])


def _adduct(led, pid):
    return str(led.loc[led["peak_id"] == pid, "adduct"].iloc[0])


led = build_ledger()
cfg = PassConfig()
res = passes.rearbitrate_offcal_degenerate(led, cfg, log=lambda *_: None)

check("exactly one off-cal monster swapped", res["swapped"] == 1, res)
check("P_swap -> on-cal C22H23N3O7 (the corroborated/plausible member)",
      _nf(led, "P_swap") == "C22H23N3O7", _nf(led, "P_swap"))
check("P_swap adduct moves to the alternative's [M+Na]+",
      _adduct(led, "P_swap") == "[M+Na]+", _adduct(led, "P_swap"))
check("P_swap ppm re-stamped to the on-cal alternative (+0.15)",
      abs(float(led.loc[led["peak_id"] == "P_swap", "ppm_error"].iloc[0]) - 0.15) < 1e-6)
check("P_swap method records the re-arbitration provenance",
      str(led.loc[led["peak_id"] == "P_swap", "method"].iloc[0]).startswith("rearb<-"))
check("P_swap commentary names the displaced off-cal monster",
      "C36H17N" in str(led.loc[led["peak_id"] == "P_swap", "commentary"].iloc[0]))
check("P_badalt KEPT (only on-cal alternative is an implausible carbon cluster)",
      _nf(led, "P_badalt") == "C30H14", _nf(led, "P_badalt"))
check("P_corrob KEPT (off-cal but isotopologue-corroborated)",
      _nf(led, "P_corrob") == "C36H17N", _nf(led, "P_corrob"))
check("P_oncal KEPT (winner is on calibration)",
      _nf(led, "P_oncal") == "C36H17N", _nf(led, "P_oncal"))

# uncalibrated -> no-op (a small ledger with no backbone fits no calibration)
small = L.new_ledger(pd.DataFrame([("x", 464.143, 5e3)], columns=["peak_id", "mz", "height"]))
L.commit_assignment(small, "x", neutral_formula="C36H17N", adduct="[M+H]+",
                    ion_formula="C36H18N+", ion_score=0.88, compound_score=0.88,
                    ppm_error=-1.0, pass_no=1, method="cheminfo+grid",
                    confidence="Good", commentary="p1",
                    alternatives=[{"formula": "C22H23N3O7", "adduct": "[M+Na]+",
                                   "ion_score": 0.70, "raw_score": 0.70, "ppm": 0.15}])
res0 = passes.rearbitrate_offcal_degenerate(small, cfg, log=lambda *_: None)
check("uncalibrated ledger -> no swap (gate needs the tier calibration)",
      res0["swapped"] == 0 and _nf(small, "x") == "C36H17N", res0)


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

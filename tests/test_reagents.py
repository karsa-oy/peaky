"""Offline tests for reagents.py. Run: python3 tests/test_reagents.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import reagents as RG  # noqa: E402
from mascope_assign import ledger as L  # noqa: E402

PASS = FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


def near(lib, mz, ppm=10):
    return [lbl for lbl, m, _f in lib if abs(m - mz) / mz * 1e6 <= ppm]


def near_f(lib, mz, ppm=10):
    return [f for _lbl, m, f in lib if abs(m - mz) / mz * 1e6 <= ppm]


# --- library contains Br-, [Br3]- and isotopologues at the right masses ---
lib = RG.build_library("Br")
check("Br- present ~78.9189", bool(near(lib, 78.9189)), near(lib, 78.9189))
# tribromide [Br3]- monoisotopic 3*78.9183 + e = 236.7555
check("[Br3]- present ~236.7555", bool(near(lib, 236.7555)), near(lib, 236.7555))
# isotopologue 79Br2 81Br at ~238.7535
check("[Br3]- 79,79,81 isotopologue ~238.7535", bool(near(lib, 238.7535)), near(lib, 238.7535))
# Br . H2O cluster ~ 78.9189 + 18.0106 = 96.929 (HNO3/HNO2 were removed from
# the cluster library 2026-06-12 -- they are ambient analytes assigned in pass 0)
check("[Br+H2O]- present ~96.929", bool(near(lib, 96.929)), near(lib, 96.929))
check("HNO3 NOT in reagent library (now an analyte)", not near(lib, 141.914))
# di-bromide radical anion Br2-. = 2*78.9183 + e = 157.8372 (user registered it
# on the server 2026-06-12); the labeler must catch bare even-n clusters too
check("[Br2]-. present ~157.8372", bool(near(lib, 157.8372)), near(lib, 157.8372))
# ambient ORGANIC ACIDS were removed from the cluster library: [Br+HCOOH]- =
# [CH2O2+Br]- = the analyte channel (formic acid's 124.92/126.92 giants), so it
# must NOT be a reagent label anymore
check("[Br+HCOOH]- (124.924) NOT in reagent library (it is the [M+Br]- analyte)",
      not near(lib, 124.924), near(lib, 124.924))
check("[Br+pinic]- (267.006) NOT in reagent library", not near(lib, 267.006))
# HBr cluster on the di-bromide core stays reagent (pure halogen, no analyte):
# [Br2+HBr]- = 157.8372 + 79.926 = 237.763
check("[Br+HBr]- (HBr2- ~160.843) still reagent", bool(near(lib, 160.843)))
# [Br+HF]- = BrHF- = 98.9251 -- HF background halogen cluster (v47 time-series ID;
# the only clean fit among the variable-unassigned residual). Both isotopologues.
check("[Br+HF]- (BrHF- ~98.9251) present", bool(near(lib, 98.9251)), near(lib, 98.9251))
check("[Br+HF]- 81Br twin (~100.9231) present", bool(near(lib, 100.9231)))
check("[Br+HF]- carries ion_formula HBrF-", "HBrF-" in near_f(lib, 98.9251))

# --- reagent_for_adducts ---
check("Br reagent from [M+Br]-", RG.reagent_for_adducts(["[M-H]-", "[M+Br]-"]) == "Br")
check("I reagent from [M+I]-", RG.reagent_for_adducts(["[M+I]-"]) == "I")
check("None when no halide reagent", RG.reagent_for_adducts(["[M-H]-", "[M+NO3]-"]) is None)

# --- labeler marks the bright Br3 cluster peaks ---
peaks = pd.DataFrame({
    "peak_id": ["b1", "b3", "b3b", "org"],
    "mz": [78.9189, 236.7555, 238.7535, 257.0181],
    "height": [2e5, 1e5, 9e4, 8e4],
})
led = L.new_ledger(peaks)
n = RG.label_reagents(led, "Br", ppm=15)
check("labels >=3 reagent peaks", n >= 3, n)
check("Br3 peak labeled reagent",
      L.role_of(led, "b3") == L.ROLE_REAGENT, L.role_of(led, "b3"))
check("organic peak NOT labeled reagent",
      L.role_of(led, "org") == L.ROLE_UNEXPLAINED, L.role_of(led, "org"))
check("reagent commentary written",
      "reagent ion" in str(led.loc[led.peak_id == "b3", "commentary"].iloc[0]))
# known formula -> assigned: the reagent row must carry its ion_formula
check("Br3 reagent row records ion_formula Br3-",
      str(led.loc[led.peak_id == "b3", "ion_formula"].iloc[0]) == "Br3-",
      led.loc[led.peak_id == "b3", "ion_formula"].iloc[0])

# --- BOTH BrO- isotopologues present (the 81Br twin at 96.91 was being dropped) ---
check("79BrO- present ~94.9138", bool(near(lib, 94.9138)), near(lib, 94.9138))
check("81BrO- present ~96.9118 (the missed twin)", bool(near(lib, 96.9118)), near(lib, 96.9118))
check("BrO- ion formula recorded as BrO-", "BrO-" in near_f(lib, 96.9118), near_f(lib, 96.9118))
# and the labeler assigns it out of the residual
led2 = L.new_ledger(pd.DataFrame({"peak_id": ["bo"], "mz": [96.9117], "height": [2350.0]}))
RG.label_reagents(led2, "Br", ppm=15)
check("81BrO- peak (96.91) now labeled reagent, not unexplained",
      L.role_of(led2, "bo") == L.ROLE_REAGENT, L.role_of(led2, "bo"))
check("81BrO- peak carries ion_formula BrO-",
      str(led2.loc[led2.peak_id == "bo", "ion_formula"].iloc[0]) == "BrO-")

# --- does not touch assigned peaks ---
led2 = L.new_ledger(peaks)
L.commit_assignment(led2, "b3", neutral_formula="C5H8O2", adduct="[M-H]-",
                    ion_score=0.9, pass_no=1, method="x", confidence="High",
                    commentary="real assignment")
RG.label_reagents(led2, "Br", ppm=15)
check("assigned peak not overwritten by reagent labeler",
      L.role_of(led2, "b3") == L.ROLE_M0)

# --- POSITIVE molecular reagent: the urea (uronium) cluster library ----------
ulib = RG.build_library("urea")
# [urea_n + H]+ at 61.0396 / 121.0720 / 181.1044 / 241.1368
check("[urea+H]+ present ~61.0396", bool(near(ulib, 61.0396)), near(ulib, 61.0396))
check("[urea2+H]+ present ~121.0720", bool(near(ulib, 121.0720)), near(ulib, 121.0720))
check("[urea3+H]+ present ~181.1044", bool(near(ulib, 181.1044)), near(ulib, 181.1044))
check("[urea4+H]+ present ~241.1368", bool(near(ulib, 241.1368)), near(ulib, 241.1368))
# spacing is exactly one urea (60.0324)
umasses = sorted(m for _l, m, _f in ulib)
check("urea cluster spacing ~60.0324",
      abs((umasses[1] - umasses[0]) - 60.0324) < 1e-3, umasses[1] - umasses[0])
# ion formulae are CATIONS with the known elemental composition
check("[urea+H]+ ion_formula CH5N2O+", "CH5N2O+" in near_f(ulib, 61.0396), near_f(ulib, 61.0396))
check("[urea2+H]+ ion_formula C2H9N4O2+", "C2H9N4O2+" in near_f(ulib, 121.0720), near_f(ulib, 121.0720))

# reagent_for_adducts maps the urea adduct -> 'urea', halogens unchanged
check("reagent_for_adducts urea", RG.reagent_for_adducts(["[M+(CH4N2O)H]+", "[M+H]+"]) == "urea")
check("reagent_for_adducts Br unchanged", RG.reagent_for_adducts(["[M+Br]-", "[M-H]-"]) == "Br")
check("reagent_for_adducts bare positive -> None",
      RG.reagent_for_adducts(["[M+H]+", "[M+Na]+"]) is None,
      RG.reagent_for_adducts(["[M+H]+", "[M+Na]+"]))

# the labeler pulls urea clusters out of the positive residual + records formula
uled = L.new_ledger(pd.DataFrame({
    "peak_id": ["u1", "u2", "u3", "org"],
    "mz": [61.0396, 121.0720, 181.1044, 158.1536],   # org = a real [M+H]+ analyte
    "height": [5e4, 9e4, 4e4, 2.5e5],
}))
nu = RG.label_reagents(uled, "urea", ppm=15)
check("labels >=3 urea reagent peaks", nu >= 3, nu)
check("urea2 peak labeled reagent", L.role_of(uled, "u2") == L.ROLE_REAGENT, L.role_of(uled, "u2"))
check("urea2 reagent row records ion_formula C2H9N4O2+",
      str(uled.loc[uled.peak_id == "u2", "ion_formula"].iloc[0]) == "C2H9N4O2+",
      uled.loc[uled.peak_id == "u2", "ion_formula"].iloc[0])
check("positive analyte peak NOT labeled reagent",
      L.role_of(uled, "org") == L.ROLE_UNEXPLAINED, L.role_of(uled, "org"))

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

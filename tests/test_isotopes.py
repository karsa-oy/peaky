"""Offline tests for isotopes.py prescan. Run: python3 tests/test_isotopes.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import isotopes as ISO  # noqa: E402

PASS = FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


def peaks(rows):
    return pd.DataFrame(rows, columns=["mz", "height"])


# --- brominated species: parent + 81Br satellite at ~1:1 ---
br = peaks([(200.0, 1.0e5), (200.0 + ISO.D_81BR, 0.98e5),  # Br pair
            (201.00336, 0.10e5)])                          # a 13C on the parent (~9 C)
r = ISO.prescan(br)
check("detects Br (1:1 pair)", r.has_Br, r.as_dict())
check("Br not flagged multi without triplet", not r.has_multi_Br)
check("13C carbon estimate ~9", 7 <= r.estimated_max_C <= 11, r.estimated_max_C)

# --- Br2 triplet -> multi-Br ---
br2 = peaks([(300.0, 1.0e5), (300.0 + ISO.D_81BR, 1.9e5), (300.0 + 2 * ISO.D_81BR, 0.95e5)])
r2 = ISO.prescan(br2)
check("detects multi-Br triplet", r2.has_Br and r2.has_multi_Br, r2.as_dict())

# --- chlorinated: +1.997 at ~0.32 ---
cl = peaks([(150.0, 1.0e5), (150.0 + ISO.D_37CL, 0.32e5)])
rc = ISO.prescan(cl)
check("detects Cl (~0.32 ratio)", rc.has_Cl and not rc.has_Br, rc.as_dict())

# --- sulfur: +1.9958 at ~0.045 ---
s = peaks([(180.0, 1.0e5), (180.0 + ISO.D_34S, 0.045e5)])
rs = ISO.prescan(s)
check("detects S (~0.045 ratio)", rs.has_S, rs.as_dict())

# --- silicon: +0.99957 at ~0.05 ---
si = peaks([(220.0, 1.0e5), (220.0 + ISO.D_29SI, 0.05e5)])
rsi = ISO.prescan(si)
check("detects Si (~0.05 ratio)", rsi.has_Si, rsi.as_dict())

# --- clean CHO: no heteroatoms flagged (only a 13C satellite) ---
cho = peaks([(185.0, 1.0e5), (185.0 + ISO.D_13C, 0.10e5), (171.0, 5e3)])
rcho = ISO.prescan(cho)
check("clean CHO: no Br/Cl/S/Si", not (rcho.has_Br or rcho.has_Cl or rcho.has_S or rcho.has_Si), rcho.as_dict())

# --- reagent stripping: a bare Br pair declared reagent is ignored ---
rr = ISO.prescan(br, reagent_mzs=[200.0, 200.0 + ISO.D_81BR])
check("reagent peaks stripped -> Br not detected", not rr.has_Br, rr.as_dict())

# --- constrain_ranges: zero out unseen elements, cap C ---
base = {"C": (0, 40), "H": (0, 80), "O": (0, 30), "N": (0, 5),
        "S": (0, 1), "Br": (0, 2), "Cl": (0, 2), "Si": (0, 1)}
caps = {"S": 1, "Br": 2, "Cl": 2, "Si": 1}
constrained = ISO.constrain_ranges(base, r, caps)   # r had Br + C~9
check("Br kept (evidence)", constrained["Br"] == (0, 2), constrained["Br"])
check("Cl zeroed (no evidence)", constrained["Cl"] == (0, 0), constrained["Cl"])
check("S zeroed (no evidence)", constrained["S"] == (0, 0), constrained["S"])
check("Si zeroed (no evidence)", constrained["Si"] == (0, 0), constrained["Si"])
check("C capped near estimate+headroom", constrained["C"][1] <= 15, constrained["C"])

# --- context with Br cap 0 forces Br off even with evidence ---
caps_nobr = {"S": 1, "Br": 0, "Cl": 0, "Si": 0}
c2 = ISO.constrain_ranges(base, r, caps_nobr)
check("context Br cap 0 -> Br zeroed despite evidence", c2["Br"] == (0, 0), c2["Br"])

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

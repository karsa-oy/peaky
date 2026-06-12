"""Offline tests for series_gka.py. Run: python3 tests/test_series_gka.py"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import series_gka as G  # noqa: E402
from mascope_assign import chemistry as C  # noqa: E402

PASS = FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


# --- unit masses ---
check("CH2 mass ~14.0157", abs(G.unit_mass("CH2") - 14.01565) < 1e-3, G.unit_mass("CH2"))
check("CF2 mass ~49.9968", abs(G.unit_mass("CF2") - 49.99681) < 1e-3, G.unit_mass("CF2"))
check("siloxane C2H6OSi mass ~74.0188",
      abs(G.unit_mass("C2H6OSi") - 74.01879) < 1e-3, G.unit_mass("C2H6OSi"))

# --- formula arithmetic ---
check("C9H14O4 + CH2 = C10H16O4", G.formula_add("C9H14O4", "CH2", 1) == "C10H16O4",
      G.formula_add("C9H14O4", "CH2", 1))
check("C10H16O4 - CH2 = C9H14O4", G.formula_add("C10H16O4", "CH2", -1) == "C9H14O4")
check("C9H14O4 + O = C9H14O5", G.formula_add("C9H14O4", "O", 1) == "C9H14O5")
check("negative element -> None", G.formula_add("CH4", "CO2", -1) is None)
check("add CF2 introduces F", G.formula_add("C6H4O2", "CF2", 1) == "C7H4F2O2",
      G.formula_add("C6H4O2", "CF2", 1))

# --- GKA math: KMD of a CH2 series is ~constant ---
mz1 = C.ion_mz("C9H14O4", "[M-H]-")
mz2 = C.ion_mz("C10H16O4", "[M-H]-")   # +CH2
d1 = G.gkd(mz1, "CH2", 14)
d2 = G.gkd(mz2, "CH2", 14)
check("CH2 series has near-equal Kendrick defect", abs(d1 - d2) < 1e-3, (d1, d2))

# --- propagation: recover C10H16O4 from anchor C9H14O4 (+CH2) ---
target = C.ion_mz("C10H16O4", "[M-H]-")
props = G.propose_for_peak(target, {"C9H14O4", "C11H18O4"}, ["[M-H]-"],
                           units=G.ORGANIC_UNITS, ppm=3.0, max_steps=1)
check("propagation yields proposals", len(props) > 0, len(props))
best = props[0]
check("top proposal is C10H16O4 via CH2",
      best.neutral_formula == "C10H16O4" and best.unit == "CH2",
      (best.neutral_formula, best.unit, best.n_steps))
check("two-anchor support recognised (C9 below, C11 above)",
      best.n_supporting_anchors >= 2, best.n_supporting_anchors)
check("ppm error tiny", abs(best.ppm_error) < 1.0, best.ppm_error)

# --- propagation respects ppm gate (no spurious matches for a random mass) ---
props_none = G.propose_for_peak(123.45678, {"C9H14O4"}, ["[M-H]-"], ppm=1.0)
check("no proposal for unrelated m/z", len(props_none) == 0, len(props_none))

# --- siloxane series propagation ---
# D4 (C8H24O4Si4) -> D5 (C10H30O5Si5) is +C2H6OSi
d4 = "C8H24O4Si4"
d5_mz = C.ion_mz("C10H30O5Si5", "[M+H]+")
sp = G.propose_for_peak(d5_mz, {d4}, ["[M+H]+"], units=("C2H6OSi",), ppm=3.0)
check("siloxane D4->D5 propagation", sp and sp[0].neutral_formula == "C10H30O5Si5",
      [p.neutral_formula for p in sp[:3]])

# --- homolog series detection ---
series = {f: C.neutral_mass(f) for f in
          ["C6H12O2", "C7H14O2", "C8H16O2", "C9H18O2", "C12H10"]}
chains = G.find_homolog_series(series, "CH2", min_len=3)
check("detects CH2 chain of length 4",
      any(len(ch) == 4 for ch in chains), [len(c) for c in chains])

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

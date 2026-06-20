"""Offline test for gka_figure.py — series detection, KMD math, and a render
smoke test from a synthetic ledger. Run: python3 tests/test_gka_figure.py"""
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import gka_figure as GF  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


# a clean CH2 ladder + a separate O-addition ladder, plus a lone peak
CH2_CHAIN = ["C5H10O2", "C6H12O2", "C7H14O2", "C8H16O2", "C9H18O2"]
O_CHAIN = ["C10H16O2", "C10H16O3", "C10H16O4", "C10H16O5"]
LEDGER = pd.DataFrame(
    [dict(role="M0", neutral_formula=f) for f in CH2_CHAIN + O_CHAIN]
    + [dict(role="M0", neutral_formula="C6H18O3Si3"),       # siloxane singleton
       dict(role="iso_child", neutral_formula="C5H10O2")]   # must be ignored (not M0)
)

fmass = GF._neutral_masses(LEDGER)
check("_neutral_masses: drops non-M0 rows / dedups",
      len(fmass) == len(set(CH2_CHAIN + O_CHAIN)) + 1, len(fmass))

# --- series detection ------------------------------------------------------
ch2 = [s for s in GF.detect_series(fmass, units=["CH2"], min_len=3)]
check("detect_series: finds the CH2 ladder", len(ch2) == 1 and ch2[0].length == 5,
      [(s.unit, s.length) for s in ch2])
check("detect_series: members ascend in mass",
      ch2 and ch2[0].masses == sorted(ch2[0].masses), ch2[0].masses if ch2 else None)

oser = [s for s in GF.detect_series(fmass, units=["O"], min_len=3)]
check("detect_series: finds the O ladder", len(oser) == 1 and oser[0].length == 4,
      [(s.unit, s.length) for s in oser])

check("detect_series: min_len filters short chains",
      GF.detect_series(fmass, units=["CH2"], min_len=6) == [],
      GF.detect_series(fmass, units=["CH2"], min_len=6))

# --- KMD math --------------------------------------------------------------
masses = np.array([fmass[f] for f in CH2_CHAIN])
k = GF.kmd(masses, "CH2")
check("kmd: a CH2 ladder is one horizontal row (constant KMD)",
      float(k.max() - k.min()) < 1e-6, float(k.max() - k.min()))
# O additions move the CH2-KMD (not a CH2 homology) -> not constant
ko = GF.kmd(np.array([fmass[f] for f in O_CHAIN]), "CH2")
check("kmd: an O ladder is NOT flat under the CH2 base",
      float(ko.max() - ko.min()) > 1e-3, float(ko.max() - ko.min()))

# --- family summary --------------------------------------------------------
fam = GF.family_summary(GF.detect_series(fmass, min_len=3))
byname = {f["family"]: f for f in fam}
check("family_summary: alkyl family present", byname["alkyl"]["n_series"] >= 1, byname["alkyl"])
check("family_summary: oxidation family present", byname["oxidation"]["n_series"] >= 1,
      byname["oxidation"])
check("family_summary: every FAMILY represented",
      set(byname) == {lab for lab, *_ in GF.FAMILIES}, set(byname))

# --- render smoke test -----------------------------------------------------
with tempfile.TemporaryDirectory() as d:
    png = f"{d}/gka.png"
    out = GF.render_gka(LEDGER, png, title="synthetic")
    check("render_gka: PNG created", os.path.exists(out) and os.path.getsize(out) > 5000,
          os.path.getsize(out) if os.path.exists(out) else "missing")
    check("render_gka: it is a PNG", open(out, "rb").read(8) == b"\x89PNG\r\n\x1a\n")

    # a merged-style ledger (no 'role' column) must also render
    merged = pd.DataFrame({"neutral_formula": CH2_CHAIN + O_CHAIN, "mz": range(9)})
    out2 = GF.render_gka(merged, f"{d}/gka2.png", title="merged")
    check("render_gka: works on a role-less merged ledger", os.path.exists(out2))

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

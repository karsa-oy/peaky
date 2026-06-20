"""Offline tests for composition.py (signal-weighting + amine-shadow accounting).
Run: python3 tests/test_composition.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import composition as CO  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


# backbone classification
check("backbone CHO", CO.backbone("C10H16O2") == "CHO")
check("backbone CHON", CO.backbone("C10H19NO2") == "CHON")
check("backbone CHOS (S wins over N)", CO.backbone("C5H11NO2S") == "CHOS")
check("backbone Si folds into CHO", CO.backbone("C6H18O3Si3") == "CHO")

# minus_nh3 shadow twin
check("minus_nh3 removes N1H3", CO.minus_nh3("C10H19NO2") == "C10H16O2")
check("minus_nh3 drops N when it hits zero", CO.minus_nh3("C8H17N") == "C8H14")
check("minus_nh3 None when no N", CO.minus_nh3("C10H16O2") is None)
check("minus_nh3 None when <3 H", CO.minus_nh3("CHN") is None)

# a tiny ledger: an amine WITH its CHO twin (shadow), an amine WITHOUT, and CHO/CHOS
merged = pd.DataFrame([
    dict(neutral_formula="C10H16O2", adduct="[M+H]+"),     # CHO, also the twin of the amine below
    dict(neutral_formula="C10H19NO2", adduct="[M+H]+"),    # CHON, shadow of C10H16O2 (= +NH3)
    dict(neutral_formula="C7H9N", adduct="[M+H]+"),        # CHON, no twin (C7H6 absent)
    dict(neutral_formula="C5H11NO2S", adduct="[M+H]+"),    # CHOS
])

cnt = CO.count_by_backbone(merged)
check("count_by_backbone", cnt == {"CHO": 1, "CHON": 2, "CHOS": 1}, cnt)

st = CO.amine_shadow_stats(merged)
check("shadow: 4 neutrals", st["n_neutrals"] == 4, st)
check("shadow: 3 N-bearing (incl. the CHOS amine)", st["n_amine"] == 3, st)
check("shadow: 1 shadowed (the one with a present twin)", st["n_shadowed"] == 1, st)
check("shadow: collapsed drops the duplicate amine", st["collapsed_neutrals"] == 3, st)

asg, coll, ndrop = CO.collapsed_composition(merged)
check("collapsed: as-assigned matches count", asg == {"CHO": 1, "CHON": 2, "CHOS": 1}, asg)
check("collapsed: ammonium-as-CHO drops 1 CHON", coll == {"CHO": 1, "CHON": 1, "CHOS": 1}, coll)
check("collapsed: n_collapsed == 1", ndrop == 1, ndrop)

# signal weighting: make the lone CHO bright and the amines dim
sig = {"C10H16O2": 100000.0, "C10H19NO2": 500.0, "C7H9N": 200.0, "C5H11NO2S": 1000.0}
frac, abso = CO.signal_by_backbone(merged, sig)
check("signal: CHO dominates by signal despite 1/4 of the count",
      frac["CHO"] > 0.95, frac)
check("signal: CHON tiny by signal", frac["CHON"] < 0.02, frac)
check("signal: fractions sum ~1", abs(sum(frac.values()) - 1.0) < 1e-9, frac)

top = CO.top_species_by_signal(merged, sig, n=2)
check("top_species: brightest first", top[0]["neutral_formula"] == "C10H16O2", top)
check("top_species: carries class + frac", top[0]["klass"] == "CHO" and top[0]["frac"] > 0.9, top)

# oligomer flag
olig = pd.DataFrame([dict(neutral_formula="C20H30O14"), dict(neutral_formula="C10H16O2"),
                     dict(neutral_formula="C19H30O12"), dict(neutral_formula="C6H18O3Si3"),
                     dict(neutral_formula="C44H29NO7")])
flag = CO.oligomer_flag(olig)
check("oligomer: high-C high-O flagged, sorted by C desc",
      flag == ["C20H30O14", "C19H30O12"], flag)
check("oligomer: excludes siloxane and small species", "C6H18O3Si3" not in flag and "C10H16O2" not in flag, flag)
check("oligomer: c_max excludes absurd high-C coincidences (C44)", "C44H29NO7" not in flag, flag)

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

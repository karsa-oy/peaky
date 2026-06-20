"""Offline tests for plausibility.py (chemical-plausibility QC of assignments).
Run: python3 tests/test_plausibility.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import plausibility as PL  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


# --- implausible() rules ---
check("N-monster high O/C flagged (Candidate)",
      PL.implausible("C9H12N4O12", tier="Candidate") is not None)
check("N5O20 monster flagged", PL.implausible("C15H37N5O20", tier="Candidate") is not None)
check("very low H/C flagged", PL.implausible("C35H4", tier="Candidate") is not None)
check("positive-mode halogen flagged",
      PL.implausible("C10H18Br2O12", tier="Candidate", polarity="+") is not None)
check("same halogen NOT flagged in negative mode",
      PL.implausible("C10H18Br2O12", tier="Candidate", polarity="-") is None
      or "halogen" not in (PL.implausible("C10H18Br2O12", tier="Candidate", polarity="-") or ""))

# real molecules are NOT flagged
check("ordinary CHO not flagged", PL.implausible("C10H16O2", tier="Candidate") is None)
check("monoterpene-oxidation product not flagged", PL.implausible("C10H16O5", tier="Candidate") is None)
check("small N1 amine not flagged", PL.implausible("C10H19NO2", tier="Candidate") is None)
check("modest N3 low-O species not flagged (O/C<1)",
      PL.implausible("C12H21N3O3", tier="Candidate") is None)
check("HOM dimer not flagged", PL.implausible("C20H30O14", tier="Candidate") is None)

# Identified tier is never second-guessed
check("Identified N-monster NOT flagged (corroborated by isotope score)",
      PL.implausible("C9H12N4O12", tier="Identified") is None)

# --- scan() de-duplicates and respects best-tier-per-neutral ---
merged = pd.DataFrame([
    dict(neutral_formula="C9H12N4O12", adduct="[M+Na]+", tier="Candidate", ion_score=0.94),
    dict(neutral_formula="C10H16O2", adduct="[M+H]+", tier="Identified", ion_score=0.9),
    # same monster appears Identified in one channel -> must NOT be flagged
    dict(neutral_formula="C8H13N3O10", adduct="[M+H]+", tier="Candidate", ion_score=0.6),
    dict(neutral_formula="C8H13N3O10", adduct="[M+NH4]+", tier="Identified", ion_score=0.8),
    dict(neutral_formula="C10H18Br2O12", adduct="[M+H]+", tier="Candidate", ion_score=0.65),
])
flagged = PL.scan(merged, polarity="+")
names = {d["neutral_formula"] for d in flagged}
check("scan flags the Candidate-only N-monster", "C9H12N4O12" in names, names)
check("scan flags the positive-mode dibromo", "C10H18Br2O12" in names, names)
check("scan excludes a neutral that is Identified in any channel",
      "C8H13N3O10" not in names, names)
check("scan excludes ordinary CHO", "C10H16O2" not in names, names)
check("scan carries reason + score", all("reason" in d and "ion_score" in d for d in flagged))

# negative mode: the dibromo is NOT flagged for halogen (covalent organohalogen is plausible)
flagged_neg = PL.scan(merged, polarity="-")
check("negative mode does not flag the dibromo on halogen grounds",
      "C10H18Br2O12" not in {d["neutral_formula"] for d in flagged_neg},
      {d["neutral_formula"] for d in flagged_neg})

# empty / missing columns degrade gracefully
check("scan empty frame -> []", PL.scan(pd.DataFrame()) == [])

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

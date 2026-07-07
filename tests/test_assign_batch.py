"""Offline tests for assign_batch.py pure align/merge. Run: python3 tests/test_assign_batch.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky import assign_batch as AB  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


def m0(rows):
    return pd.DataFrame(rows, columns=["mz", "neutral_formula", "adduct", "tier", "ion_score"])


# fileA + fileB share a peak (C10H16O2 @ ~217.12, mass jitter ~2 ppm), each has
# one unique peak; fileB also assigns the shared peak Candidate (A is Assigned).
A = m0([(217.1200, "C10H16O2", "[M+H]+", "Assigned", 0.91),
        (158.1539, "C9H19NO", "[M+H]+", "Assigned", 0.85)])      # A-only
B = m0([(217.1205, "C10H16O2", "[M+H]+", "Candidate", 0.72),
        (300.1000, "C15H17NO5", "[M+H]+", "Assigned", 0.80)])    # B-only

merged, jitter = AB.align({"A": A, "B": B}, tol_ppm=6.0)

check("3 distinct clusters (1 shared + 2 unique)", len(merged) == 3, len(merged))
shared = merged[merged["neutral_formula"] == "C10H16O2"].iloc[0]
check("shared peak seen in 2 files", shared["n_files"] == 2, shared.to_dict())
check("best tier wins (Assigned over Candidate)", shared["tier"] == "Assigned", shared["tier"])
check("formula_agree True for the shared peak", bool(shared["formula_agree"]))
check("raw mz jitter ~2.3 ppm measured",
      2.0 <= shared["mz_jitter_ppm_raw"] <= 2.6, shared["mz_jitter_ppm_raw"])
check("A-only peak flagged single-file",
      merged[merged["neutral_formula"] == "C9H19NO"].iloc[0]["n_files"] == 1)
check("jitter long-form has 4 rows (2 shared + 2 unique)", len(jitter) == 4, len(jitter))

# --- offset-awareness: same peak, files at different calibrations ------------
# fileC at +3 ppm, fileD at -3 ppm -> raw mz differ by ~6 ppm; with offsets the
# corrected positions coincide and they cluster as ONE peak.
mz0 = 250.0
C = m0([(mz0 * (1 + 3e-6), "C12H20O5", "[M+H]+", "Assigned", 0.88)])
D = m0([(mz0 * (1 - 3e-6), "C12H20O5", "[M+H]+", "Assigned", 0.87)])
mC, _ = AB.align({"C": C, "D": D}, tol_ppm=4.0)                    # raw spread 6ppm > 4
check("without offsets: 6ppm raw split into 2 clusters at tol 4", len(mC) == 2, len(mC))
mO, _ = AB.align({"C": C, "D": D}, tol_ppm=4.0, offsets={"C": 3.0, "D": -3.0})
check("offset-aware: corrected positions coincide -> 1 cluster", len(mO) == 1, len(mO))
check("offset-aware: raw jitter ~6 ppm but cal-adjusted ~0",
      len(mO) == 1 and mO.iloc[0]["mz_jitter_ppm_raw"] >= 5.5
      and mO.iloc[0]["mz_jitter_ppm_caldj"] < 0.5,
      mO.iloc[0][["mz_jitter_ppm_raw", "mz_jitter_ppm_caldj"]].to_dict() if len(mO) else "empty")

# --- formula disagreement is detected ---------------------------------------
E = m0([(400.0000, "C20H25NO7", "[M+H]+", "Assigned", 0.7)])
F = m0([(400.0010, "C16H29NO10", "[M+H]+", "Candidate", 0.6)])     # different formula, same m/z
mEF, _ = AB.align({"E": E, "F": F}, tol_ppm=6.0)
check("formula disagreement flagged (formula_agree False)",
      len(mEF) == 1 and not bool(mEF.iloc[0]["formula_agree"]), mEF.to_dict("records"))

# --- cross-file consensus: corroborated formula beats single-file outlier ----
# The Ur+ m/z 424.218 bug: a mass-degenerate competitor reads on-cal
# (Assigned) with a marginally higher local score in ONE file's calibration,
# while five files agree on the real reflist HOM. The old "best (tier, ion_score)
# row" let the single-file outlier win; the consensus vote must pick the formula
# Assigned across the most files.
def _file(mz, nf, tier, ion):
    return m0([(mz, nf, "[M+NH4]+", tier, ion)])

consensus = {
    "f1": _file(424.2177, "C18H30O10", "Assigned", 0.944),
    "f2": _file(424.2179, "C18H30O10", "Assigned", 0.922),
    "f3": _file(424.2177, "C18H30O10", "Assigned", 0.985),
    "f4": _file(424.2178, "C18H30O10", "Assigned", 0.961),
    "f5": _file(424.2177, "C18H30O10", "Assigned", 0.987),
    "f6": m0([(424.2165, "C17H31N5O6", "[M+Na]+", "Candidate", 0.939)]),
    "f7": m0([(424.2167, "C17H31N5O6", "[M+Na]+", "Assigned", 0.996)]),  # outlier
}
mc, _ = AB.align(consensus, tol_ppm=6.0)
row = mc.iloc[(mc["mz"] - 424.2174).abs().argmin()]
check("consensus winner is the 5-file Assigned formula, not the 1-file outlier",
      row["neutral_formula"] == "C18H30O10", row.to_dict())
check("consensus winner keeps Assigned tier", row["tier"] == "Assigned", row["tier"])
check("consensus cluster spans all 7 files", row["n_files"] == 7, row["n_files"])

# a single-file bright Assigned formula with no competitor is still chosen (no
# spurious override when there is nothing to out-vote).
solo, _ = AB.align({"a": m0([(500.0, "C20H30O8", "[M+H]+", "Assigned", 0.9)])}, tol_ppm=6.0)
check("solo Assigned formula chosen (vote is a no-op without a competitor)",
      len(solo) == 1 and solo.iloc[0]["neutral_formula"] == "C20H30O8")

# --- empty input ------------------------------------------------------------
me, je = AB.align({})
check("empty -> empty merged + jitter with schema",
      len(me) == 0 and "n_files" in me.columns and "cluster" in je.columns)

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

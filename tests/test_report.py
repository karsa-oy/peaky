"""Offline tests for report.py. Run: python3 tests/test_report.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import ledger as L  # noqa: E402
from mascope_assign import report as R  # noqa: E402

PASS = FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


# build a small finished ledger
peaks = pd.DataFrame({"peak_id": ["A", "B", "C", "D"],
                      "mz": [200.1, 201.1, 191.0, 999.0],
                      "height": [1e5, 1e4, 8e4, 50.0]})
led = L.new_ledger(peaks)
L.commit_assignment(led, "A", neutral_formula="C10H16O4", adduct="[M-H]-",
                    ion_formula="C10H15O4-", ion_score=0.97, compound_score=0.96,
                    ppm_error=-0.3, pass_no=1, method="cheminfo", confidence="High",
                    commentary="Pass 1: C10H16O4 [M-H]-, ion score 0.97",
                    alternatives=[{"formula": "C9H13NO5", "ion_score": 0.90, "ppm": 0.8,
                                   "eff_score": 0.87}],
                    isotopologues=[{"label": "13C", "score": 0.93, "peak_id": "B"}])
L.attach_isotopologue(led, "B", "A", iso_label="13C", iso_match_score=0.93)
L.commit_assignment(led, "C", neutral_formula="C7H12O4", adduct="[M-H]-",
                    ion_formula="C7H11O4-", ion_score=0.83, compound_score=0.83,
                    ppm_error=0.5, pass_no=2, method="gka-series", confidence="Good (series)",
                    commentary="Pass 2 series from C8H14O4 -CH2")

sheets = R.build_sheets(led, "ambient-air")
check("has all sheets", {"Assignments", "By class", "Unique formulas", "Target list",
                         "Isotopologues", "Peak ownership", "Unassigned", "Reagent ions"}
      <= set(sheets), set(sheets))
check("assignments has 2 rows", len(sheets["Assignments"]) == 2, len(sheets["Assignments"]))
check("assignment carries commentary",
      sheets["Assignments"]["commentary"].str.contains("C10H16O4").any())
check("alternatives rendered to text",
      sheets["Assignments"]["alternatives_text"].str.contains("C9H13NO5").any(),
      sheets["Assignments"]["alternatives_text"].tolist())
check("isotopologues_text rendered",
      sheets["Assignments"]["isotopologues_text"].str.contains("13C").any())
check("isotopologues sheet has child B",
      "B" in set(sheets["Isotopologues"]["peak_id"]))
check("ownership covers all 4 peaks", len(sheets["Peak ownership"]) == 4)
check("unassigned has D", "D" in set(sheets["Unassigned"]["peak_id"]))
check("compound_class assigned", "C10 monomer" in set(sheets["Assignments"]["compound_class"]))

# summary stats
ss = R.summary_stats(led)
check("summary has n_peaks", (ss["metric"] == "n_peaks").any())

# excel write (needs openpyxl)
try:
    import openpyxl  # noqa: F401
    out = Path("/tmp/_report_test.xlsx")
    R.write_excel(led, out, "ambient-air")
    check("excel written", out.exists() and out.stat().st_size > 0)
    out.unlink(missing_ok=True)
except ImportError:
    print("  (openpyxl not installed; skipping excel write)")

# markdown
result = {"ledger": led, "stats": L.stats(led), "sample_id": "TEST",
          "context": "ambient-air", "prescan": {"has_Br": False}, "problems": []}
mdp = R.write_markdown(result, "/tmp/_report_test.md")
md = mdp.read_text()
check("markdown has top assignments", "C10H16O4" in md)
check("markdown reports signal explained", "Signal explained" in md)
mdp.unlink(missing_ok=True)

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

"""Offline tests for report.py. Run: python3 tests/test_report.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky import ledger as L  # noqa: E402
from peaky import report as R  # noqa: E402

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
peaks = pd.DataFrame({"peak_id": ["A", "B", "C", "D", "E"],
                      "mz": [200.1, 201.1, 191.0, 999.0, 78.9189],
                      "height": [1e5, 1e4, 8e4, 1e4, 1e4]})
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
L.mark_reagent(led, "E", "reagent ion: [Br]-")

sheets = R.build_sheets(led, "ambient-air", sample_id="TEST")
check("has all sheets", {"Summary", "Read me", "Assigned", "Candidates",
                         "Unassigned", "By class", "Unique formulas",
                         "Isotopologues", "Peak ownership", "Target list",
                         "Reagent ions"} <= set(sheets), set(sheets))
check("summary is the first sheet", list(sheets)[0] == "Summary", list(sheets))
ident, cand = sheets["Assigned"], sheets["Candidates"]
check("A and C are both Assigned (iso-confirmed / unique)",
      len(ident) == 2 and not len(cand[cand["rank"] == 1]),
      (ident.get("neutral_formula"), cand.get("formula")))
check("identified carries commentary",
      ident["commentary"].str.contains("C10H16O4").any())
check("identified carries evidence text",
      ident["evidence"].str.contains("isotopologue").any(), ident["evidence"].tolist())
check("alternatives rendered to text",
      ident["alternatives_text"].str.contains("C9H13NO5").any(),
      ident["alternatives_text"].tolist())
check("isotopologues_text rendered",
      ident["isotopologues_text"].str.contains("13C").any())
check("isotopologues sheet has child B",
      "B" in set(sheets["Isotopologues"]["peak_id"]))
check("isotopologues joined to parent formula",
      (sheets["Isotopologues"]["parent_formula"] == "C10H16O4").any())
check("ownership covers all 5 peaks", len(sheets["Peak ownership"]) == 5)
check("ownership carries tier", "tier" in sheets["Peak ownership"].columns)
check("unassigned has D", "D" in set(sheets["Unassigned"]["peak_id"]))
check("unassigned evidence interpreted",
      sheets["Unassigned"]["interpretation"].notna().all())
check("compound_class assigned", "C10 monomer" in set(ident["compound_class"]))

# a tied Good with no corroboration must land in Candidates, expanded per formula
L.commit_assignment(led, "D", neutral_formula="C6H10O4", adduct="[M-H]-",
                    ion_formula="C6H9O4-", ion_score=0.85, compound_score=0.85,
                    eff_score=0.84, eff_margin=0.02, tied=True,
                    ppm_error=0.3, pass_no=1, method="cheminfo+grid",
                    confidence="Good", commentary="Pass 1 near-tie",
                    alternatives=[{"formula": "C2H6N2O6", "ion_score": 0.84,
                                   "raw_score": 0.84, "eff_score": 0.82, "ppm": 0.4}])
sheets2 = R.build_sheets(led, "ambient-air")
cand2 = sheets2["Candidates"]
check("tied peak lands in Candidates", (cand2["peak_id"] == "D").any())
check("candidate expanded one row per formula",
      len(cand2[cand2["peak_id"] == "D"]) == 2, len(cand2))
check("rank-1 row is the committed winner",
      cand2[(cand2["peak_id"] == "D") & (cand2["rank"] == 1)]["formula"].iloc[0] == "C6H10O4")
check("rank-2 row is the alternative",
      cand2[(cand2["peak_id"] == "D") & (cand2["rank"] == 2)]["formula"].iloc[0] == "C2H6N2O6")
check("why_candidate explains the tie",
      cand2[(cand2["peak_id"] == "D") & (cand2["rank"] == 1)]["why_candidate"]
      .str.contains("near-tie").all())

# summary stats
ss = R.summary_stats(led, context="ambient-air", sample_id="TEST")
check("summary has peak count", (ss["metric"] == "peaks total").any())
check("summary has tier rows", (ss["section"] == "Tiers").any())
check("summary has sample id", (ss["value"] == "TEST").any())

# fragment role: relabelled M0 -> fragment flows to its own sheet + ownership +
# coverage, but is kept OUT of the analyte (Assigned/Candidates/Unique) sheets.
fled = L.new_ledger(pd.DataFrame(
    {"peak_id": ["P", "Q"], "mz": [180.10, 162.09], "height": [5e4, 2e4]}))
L.commit_assignment(fled, "P", neutral_formula="C6H12O3", adduct="[M+H]+",
                    ion_formula="C6H13O3+", ion_score=0.95, compound_score=0.95,
                    ppm_error=0.1, pass_no=1, method="cheminfo+grid",
                    confidence="High", commentary="parent")
qi = fled.index[fled["peak_id"] == "Q"][0]
fled.at[qi, "role"] = L.ROLE_FRAGMENT
fled.at[qi, "neutral_formula"] = "C6H10O2"
fled.at[qi, "adduct"] = "[M+H]+"
fled.at[qi, "commentary"] = "in-source fragment of C6H12O3 (H2O)"
fsheets = R.build_sheets(fled, "ambient-air", sample_id="FRAG")
check("Fragment ions sheet present", "Fragment ions" in fsheets, list(fsheets))
fr = fsheets["Fragment ions"]
check("fragment row carries parent from commentary",
      len(fr) == 1 and fr["parent"].iloc[0] == "C6H12O3", fr.to_dict("records"))
check("fragment excluded from Assigned (not an analyte)",
      "C6H10O2" not in set(fsheets["Assigned"]["neutral_formula"]))
check("fragment excluded from Unique formulas",
      "C6H10O2" not in set(fsheets["Unique formulas"].get("neutral_formula", pd.Series(dtype=object))))
check("fragment kept in Peak ownership", "Q" in set(fsheets["Peak ownership"]["peak_id"]))
fss = R.summary_stats(fled, sample_id="FRAG")
check("coverage lists fragment ions row",
      ((fss["section"] == "Coverage") & (fss["metric"] == "fragment ions")).any())
# a ledger with NO fragments must not grow a fragment coverage row
nss = R.summary_stats(led, sample_id="TEST")
check("no fragment coverage row when there are no fragments",
      not ((nss["section"] == "Coverage") & (nss["metric"] == "fragment ions")).any())

# excel write (needs openpyxl)
try:
    import openpyxl  # noqa: F401
    out = Path("/tmp/_report_test.xlsx")
    R.write_excel(led, out, "ambient-air")
    check("excel written", out.exists() and out.stat().st_size > 0)
    out.unlink(missing_ok=True)
except ImportError:
    print("  (openpyxl not installed; skipping excel write)")

# markdown (on a tier-stamped ledger, as assign.run leaves it)
from peaky import tiers as T  # noqa: E402
T.apply_tiers(led)
result = {"ledger": led, "stats": L.stats(led), "sample_id": "TEST",
          "context": "ambient-air", "prescan": {"has_Br": False}, "problems": []}
mdp = R.write_markdown(result, "/tmp/_report_test.md")
md = mdp.read_text()
check("markdown has top assignments", "C10H16O4" in md)
check("markdown reports signal explained", "Signal explained" in md)
# every peak is now explained (D was committed in the Candidates test above)
check("markdown signal explained includes reagent ions", "Signal explained: 100.0%" in md, md)
check("markdown reports tiers", "Tiers:" in md, md)
mdp.unlink(missing_ok=True)

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

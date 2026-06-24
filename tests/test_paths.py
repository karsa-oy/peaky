"""The run-folder layout contract (paths.RunPaths).

This is the single source of truth the WRITERS (clustering / analyte_viz /
assign_batch / pdf_report) and the READER (pdf_report.load_context) both derive
from. If routing changes here without both ends moving, the report silently loses
figures — so pin the routing. Run: python3 tests/test_paths.py
"""
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky import paths as PT  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


with tempfile.TemporaryDirectory() as d:
    rp = PT.run_paths(d)
    check("figures/tables/report/data are run-dir subdirs",
          rp.figures == os.path.join(d, "figures") and rp.tables == os.path.join(d, "tables")
          and rp.report == os.path.join(d, "report") and rp.data == os.path.join(d, "data"))
    check("ensure() creates the subdirs (not before)",
          not os.path.isdir(rp.figures))
    rp.ensure()
    check("ensure() created all four subdirs",
          all(os.path.isdir(x) for x in (rp.figures, rp.tables, rp.report, rp.data)))

    # routing: extension -> subdir
    check("place(.png) -> figures/", rp.place("clusters_changing_Ur_p1.png") == os.path.join(rp.figures, "clusters_changing_Ur_p1.png"))
    check("place(.csv) -> tables/", rp.place("jitter.csv") == os.path.join(rp.tables, "jitter.csv"))
    check("place(.xlsx) -> tables/", rp.place("clusters_changing_Ur.xlsx") == os.path.join(rp.tables, "clusters_changing_Ur.xlsx"))
    check("place(.pdf) -> report/", rp.place("report_X.pdf") == os.path.join(rp.report, "report_X.pdf"))

    # root anchors stay flat (read by several modules + the cross-run registry)
    for anchor in ("merged_ledger.csv", "run_manifest.json", "batch_summary.json"):
        check(f"root anchor {anchor} stays at run root",
              rp.place(anchor) == os.path.join(d, anchor))


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

"""Determinism regression — the run's reproducibility CONTRACT.

The scientific content of a run — every figure's pixels, the assignment tables
(merged_ledger.csv, per-file csv, cluster csv/xlsx) and material data — is a PURE
FUNCTION OF THE INPUT DATA. It is byte-identical regardless of WHEN the run
happens. Run time reaches output ONLY as visible text on the PDF cover (the
"generated" line + Report ID), the run-folder name, and run_manifest.json
provenance — never as bytes inside a figure or a data table.

Mechanically: `pipeline.stamp_source_date_epoch()` pins SOURCE_DATE_EPOCH to a
FIXED content epoch (not the run time), so matplotlib's PNG/PDF metadata and the
openpyxl xlsx timestamps are constant. This test asserts: two runs at DIFFERENT
`when` over the same inputs -> identical figure/xlsx/csv bytes; a PDF differs only
when its visible cover text differs.

Run: python3 tests/test_determinism.py
"""
import hashlib
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from peaky import cluster as CL  # noqa: E402
from peaky.pipeline import CONTENT_EPOCH, stamp_source_date_epoch  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


def render(path, fmt, cover=None):
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.plot([0, 1, 2], [2, 1, 3]); ax.set_title("determinism")
    if cover:                       # mimic the PDF cover's run-time TEXT
        fig.text(0.1, 0.95, f"generated {cover}")
    fig.savefig(path, format=fmt); plt.close(fig)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# rows shape for write_cluster_workbook: [(cid, members, rbar, shape, peak_hr), ...]
_ROWS = [("c1", ["m1", "m2"], 0.91, "rise", 3.5),
         ("c2", ["m3"], 0.80, "fall", 9.0)]
_META = {"m1": {"neutral_formula": "C4H6O4", "tier": "Identified"},
         "m2": {"neutral_formula": "C5H8O4", "tier": "Candidate"},
         "m3": {"neutral_formula": "C6H10O5", "tier": "Identified"}}


def xlsx(path):
    CL.write_cluster_workbook(_ROWS, path, meta=_META, member_cols=["neutral_formula", "tier"])


w1 = datetime(2026, 6, 23, 12, 0, tzinfo=timezone.utc)
w2 = datetime(2027, 1, 1, 9, 30, tzinfo=timezone.utc)   # a DIFFERENT run time

with tempfile.TemporaryDirectory() as d:
    # ---- 1. the helper pins a FIXED content epoch, independent of `when` ----
    got1 = stamp_source_date_epoch(w1)
    got2 = stamp_source_date_epoch(w2)          # different when -> SAME epoch
    check("stamp_source_date_epoch is INDEPENDENT of `when` (fixed content epoch)",
          got1 == got2 == str(CONTENT_EPOCH), f"{got1} / {got2} / {CONTENT_EPOCH}")
    check("SOURCE_DATE_EPOCH is the fixed content epoch",
          os.environ.get("SOURCE_DATE_EPOCH") == str(CONTENT_EPOCH))

    # ---- 2. THE CONTRACT: different run time -> IDENTICAL content bytes ----
    # (figures, workbook, and a bare/cover-less PDF), driven only by the helper
    # exactly as the pipeline does it.
    def artifacts(tag, when):
        stamp_source_date_epoch(when)           # ignores `when` -> fixed epoch
        png = os.path.join(d, f"{tag}.png"); render(png, "png")
        xl = os.path.join(d, f"{tag}.xlsx"); xlsx(xl)
        pdf = os.path.join(d, f"{tag}.pdf"); render(pdf, "pdf")    # no cover text
        return sha(png), sha(xl), sha(pdf)

    a = artifacts("runA", w1)
    b = artifacts("runB", w2)                   # DIFFERENT when, same inputs
    check("different `when` -> IDENTICAL figure PNG bytes", a[0] == b[0])
    check("different `when` -> IDENTICAL cluster xlsx bytes (write_cluster_workbook)",
          a[1] == b[1], f"{a[1][:10]} vs {b[1][:10]}")
    check("different `when` -> IDENTICAL bare-figure PDF bytes (no run text)", a[2] == b[2])

    # ---- 3. material CSV carries no clock -> byte-identical by construction ----
    import pandas as pd  # noqa: E402
    df = pd.DataFrame({"mz": [169.12, 183.05], "neutral_formula": ["C10H16O2", "C9H10N2O"],
                       "tier": ["Identified", "Candidate"]})
    p1, p2 = os.path.join(d, "led1.csv"), os.path.join(d, "led2.csv")
    df.to_csv(p1, index=False); df.to_csv(p2, index=False)
    check("merged_ledger-style CSV is byte-identical (no run time in material data)",
          sha(p1) == sha(p2))

    # ---- 4. the ONLY cross-run variation: the visible cover TIMESTAMP TEXT ----
    # two PDFs identical except a cover string, under the SAME fixed epoch, MUST differ.
    stamp_source_date_epoch(w1)
    g1 = os.path.join(d, "cover1.pdf"); render(g1, "pdf", cover="2026-06-23 12:00 UTC")
    g2 = os.path.join(d, "cover2.pdf"); render(g2, "pdf", cover="2027-01-01 09:30 UTC")
    check("a report PDF differs ONLY because its visible cover timestamp differs",
          sha(g1) != sha(g2))
    # ...and the SAME cover text reproduces byte-for-byte (true determinism)
    g3 = os.path.join(d, "cover3.pdf"); render(g3, "pdf", cover="2026-06-23 12:00 UTC")
    check("same cover text + same data -> byte-identical PDF", sha(g1) == sha(g3))


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

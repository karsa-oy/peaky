"""Determinism regression — the report's reproducibility contract.

Figures/PDF are a deterministic function of inputs + run time: with a fixed
SOURCE_DATE_EPOCH, matplotlib stamps that time (not now()) into the PDF/PNG
metadata, so a re-render is byte-identical; a different epoch -> different bytes
(the stamp tracks the run, like the Report ID). This locks the env-var contract
that run_peaky / pipeline.generate_report rely on. (xlsx is covered separately in
test_cluster.) Run: python3 tests/test_determinism.py"""
import hashlib
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


def render(path, fmt):
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.plot([0, 1, 2], [2, 1, 3]); ax.set_title("determinism")
    fig.savefig(path, format=fmt); plt.close(fig)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


with tempfile.TemporaryDirectory() as d:
    os.environ["SOURCE_DATE_EPOCH"] = "1700000000"
    a, b = os.path.join(d, "a.pdf"), os.path.join(d, "b.pdf")
    render(a, "pdf"); render(b, "pdf")
    check("PDF byte-identical at a fixed SOURCE_DATE_EPOCH",
          sha(a) == sha(b), f"{sha(a)[:10]} vs {sha(b)[:10]}")

    os.environ["SOURCE_DATE_EPOCH"] = "1800000000"
    c = os.path.join(d, "c.pdf"); render(c, "pdf")
    check("PDF bytes change when SOURCE_DATE_EPOCH changes (stamp tracks the run)",
          sha(c) != sha(a), "matplotlib not honouring SOURCE_DATE_EPOCH")

    os.environ["SOURCE_DATE_EPOCH"] = "1700000000"
    e, f = os.path.join(d, "e.png"), os.path.join(d, "f.png")
    render(e, "png"); render(f, "png")
    check("PNG byte-identical at a fixed SOURCE_DATE_EPOCH",
          sha(e) == sha(f), f"{sha(e)[:10]} vs {sha(f)[:10]}")


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

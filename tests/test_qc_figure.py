"""Offline test for qc_figure.py — category splitting, ppm-point selection, mass-
defect math, and a render smoke test from a synthetic FULL per-sample ledger.
Run: python3 tests/test_qc_figure.py"""
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky import qc_figure as QC  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


# A FULL per-sample ledger: M0 owners (Assigned + Candidate), their iso_child
# satellites (which carry NO tier of their own -> inherit the parent's), reagent,
# and the unexplained residual. peak_id links iso_child -> parent.
LEDGER = pd.DataFrame([
    dict(peak_id=1, mz=169.1223, role="M0", neutral_formula="C10H16O2",
         adduct="[M+H]+", tier="Assigned", ppm_error=0.4, parent_peak_id=None),
    dict(peak_id=2, mz=170.1256, role="iso_child", neutral_formula=None,
         adduct=None, tier=None, ppm_error=None, parent_peak_id=1),
    dict(peak_id=3, mz=263.0288, role="M0", neutral_formula="C10H16O3",
         adduct="[M+Br]-", tier="Candidate", ppm_error=-1.2, parent_peak_id=None),
    dict(peak_id=4, mz=265.0268, role="iso_child", neutral_formula=None,
         adduct=None, tier=None, ppm_error=None, parent_peak_id=3),
    dict(peak_id=5, mz=124.9244, role="M0", neutral_formula="CH2O2",
         adduct="[M+Br]-", tier="Assigned", ppm_error=0.05, parent_peak_id=None),
    dict(peak_id=6, mz=250.0, role="reagent", neutral_formula=None,
         adduct=None, tier=None, ppm_error=None, parent_peak_id=None),
    dict(peak_id=7, mz=301.7, role="unexplained", neutral_formula=None,
         adduct=None, tier=None, ppm_error=None, parent_peak_id=None),
    dict(peak_id=8, mz=412.3, role="unexplained", neutral_formula=None,
         adduct=None, tier=None, ppm_error=None, parent_peak_id=None),
])

# --- mass-defect math ------------------------------------------------------
md = QC.mass_defect([100.0, 100.4, 99.6])
check("mass_defect: mz - round(mz)",
      abs(md[0]) < 1e-9 and abs(md[1] - 0.4) < 1e-9 and abs(md[2] + 0.4) < 1e-9, md)

# --- category split --------------------------------------------------------
cats = QC.split_categories(LEDGER)
check("split: Assigned·M0 = the two Assigned owners",
      sorted(cats["assigned_parent"]["peak_id"]) == [1, 5], cats["assigned_parent"]["peak_id"].tolist())
check("split: Candidate·M0 = the one Candidate owner",
      cats["cand_parent"]["peak_id"].tolist() == [3], cats["cand_parent"]["peak_id"].tolist())
check("split: iso-child inherits its parent's tier (Assigned)",
      cats["assigned_iso"]["peak_id"].tolist() == [2], cats["assigned_iso"]["peak_id"].tolist())
check("split: iso-child inherits its parent's tier (Candidate)",
      cats["cand_iso"]["peak_id"].tolist() == [4], cats["cand_iso"]["peak_id"].tolist())
check("split: unexplained cloud = the two residual peaks",
      sorted(cats["unexplained"]["peak_id"]) == [7, 8], cats["unexplained"]["peak_id"].tolist())
check("split: reagent role is NOT a category (excluded)",
      all(6 not in g["peak_id"].tolist() for g in cats.values()), cats.keys())
check("split: every category frame carries the mass-defect column",
      all("md" in g.columns for g in cats.values() if len(g)))

# --- legacy 'Identified' tier spelling still maps onto Assigned ------------
legacy = LEDGER.copy()
legacy.loc[legacy["tier"] == "Assigned", "tier"] = "Identified"
cleg = QC.split_categories(legacy)
check("split: legacy 'Identified' maps onto the Assigned category",
      sorted(cleg["assigned_parent"]["peak_id"]) == [1, 5]
      and cleg["assigned_iso"]["peak_id"].tolist() == [2],
      cleg["assigned_parent"]["peak_id"].tolist())

# --- ppm points: Assigned + Candidate M0 with a finite ppm_error ----------
pts = QC.ppm_points(LEDGER)
check("ppm_points: only the three M0 rows with a ppm_error",
      len(pts) == 3 and set(pts["tier"]) == {"Assigned", "Candidate"}, pts.to_dict("records"))
check("ppm_points: iso-child / reagent / unexplained excluded",
      pts["ppm_error"].notna().all() and "mz" in pts.columns, pts.columns.tolist())

# --- render smoke test -----------------------------------------------------
import matplotlib                                                    # noqa: E402
matplotlib.use("Agg")

with tempfile.TemporaryDirectory() as d:
    png = f"{d}/qc.png"
    out = QC.render_qc(LEDGER, png, title="synthetic")
    check("render_qc: PNG created", os.path.exists(out) and os.path.getsize(out) > 5000,
          os.path.getsize(out) if os.path.exists(out) else "missing")
    check("render_qc: it is a PNG", open(out, "rb").read(8) == b"\x89PNG\r\n\x1a\n")

    # an unexplained-only ledger (no calibrated M0) still renders panel (a)
    un = pd.DataFrame([dict(peak_id=i, mz=100.0 + i * 7.3, role="unexplained",
                            neutral_formula=None, tier=None, ppm_error=None,
                            parent_peak_id=None) for i in range(5)])
    out2 = QC.render_qc(un, f"{d}/qc2.png", title="residual only")
    check("render_qc: renders with no calibrated M0 rows", os.path.exists(out2))

    # determinism: same ledger -> byte-identical PNG (SOURCE_DATE_EPOCH pinned)
    os.environ["SOURCE_DATE_EPOCH"] = "315532800"
    a = f"{d}/det_a.png"; b = f"{d}/det_b.png"
    QC.render_qc(LEDGER, a, title="x"); QC.render_qc(LEDGER, b, title="x")
    check("render_qc: deterministic (byte-identical re-render)",
          open(a, "rb").read() == open(b, "rb").read())

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

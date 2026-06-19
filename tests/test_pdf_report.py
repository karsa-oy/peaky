"""Offline smoke test for pdf_report.py — builds a PDF from a tiny synthetic run.
Run: python3 tests/test_pdf_report.py"""
import os
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import pdf_report as R  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


with tempfile.TemporaryDirectory() as d:
    os.makedirs(f"{d}/per_file", exist_ok=True)
    # minimal merged ledger (M0 rows) spanning a few classes
    pd.DataFrame([
        dict(mz=169.1223, neutral_formula="C10H16O2", adduct="[M+H]+", tier="Identified", ion_score=0.9),
        dict(mz=183.0, neutral_formula="C9H10N2O", adduct="[M+H]+", tier="Candidate", ion_score=0.6),
        dict(mz=223.06, neutral_formula="C6H18O3Si3", adduct="[M+H]+", tier="Candidate", ion_score=0.5),
        dict(mz=247.0, neutral_formula="C3H2F6O", adduct="[M+Br]-", tier="Identified", ion_score=0.8),
    ]).to_csv(f"{d}/merged_ledger.csv", index=False)
    # one per-file ledger with roles (drives the role breakdown)
    pd.DataFrame([
        dict(mz=169.1223, role="M0", neutral_formula="C10H16O2", height=10000),
        dict(mz=170.1256, role="iso_child", neutral_formula=None, height=1100),
        dict(mz=250.0, role="reagent", neutral_formula=None, height=500),
        dict(mz=300.0, role="unexplained", neutral_formula=None, height=80),
    ]).to_csv(f"{d}/per_file/s1_ledger.csv", index=False)
    pd.DataFrame({"neutral_formula": ["C10H16O2", "C10H14O2", "C9H12O2"],
                  "cluster": [1, 1, 1], "cv": [0.5, 0.6, 0.4],
                  "median_cps": [9000, 4000, 2000]}).to_csv(f"{d}/clusters_changing_Ur.csv", index=False)

    ctx = R.load_context(d, tag="Ur", label="Ur⁺ CIMS")
    check("load_context: merged loaded", ctx["n_m0"] == 4, ctx.get("n_m0"))
    check("load_context: tiers counted", ctx["tiers"].get("Identified") == 2, ctx["tiers"])
    check("load_context: composition is CHO/CHON/CHOS backbone",
          set(ctx["composition"]) <= {"CHO", "CHON", "CHOS"}, ctx["composition"])
    check("load_context: heteroatom side-counts Si + F",
          ctx["hetero"]["Si-bearing"] == 1 and ctx["hetero"]["F-bearing"] == 1,
          ctx.get("hetero"))
    check("load_context: per-adduct channel counts present", bool(ctx.get("adduct_counts")),
          ctx.get("adduct_counts"))
    check("load_context: role breakdown present", "unexplained" in ctx.get("role_count", {}),
          ctx.get("role_count"))

    out = R.build(d, tag="Ur", label="Ur⁺ CIMS", generated="2026-01-01")
    check("build: PDF created", os.path.exists(out), out)
    check("build: PDF non-trivial size", os.path.getsize(out) > 5000, os.path.getsize(out))
    check("build: PDF magic header", open(out, "rb").read(4) == b"%PDF")

    # a failing section must not kill the build (resilience)
    def boom(ctx, pdf):
        raise ValueError("intentional")
    out2 = R.build(d, tag="Ur", label="Ur⁺ CIMS", out_pdf=f"{d}/r2.pdf",
                   sections=[R.cover, boom, R.methods])
    check("build: resilient to a failing section", os.path.exists(out2) and os.path.getsize(out2) > 3000)

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

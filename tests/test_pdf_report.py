"""Offline smoke test for pdf_report.py — builds a PDF from a tiny synthetic run.
Run: python3 tests/test_pdf_report.py"""
import os
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky import pdf_report as R  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


with tempfile.TemporaryDirectory() as d:
    os.makedirs(f"{d}/per_file", exist_ok=True)
    os.makedirs(f"{d}/tables", exist_ok=True)   # cluster tables live under tables/ (paths.RunPaths)
    # minimal merged ledger (M0 rows) spanning a few classes; C10H19NO2 is the
    # NH3-shifted shadow of C10H16O2 (the ammonium/amine degeneracy), and one row
    # has formula_agree=False (drives the single-source disagreement count).
    pd.DataFrame([
        dict(mz=169.1223, neutral_formula="C10H16O2", adduct="[M+H]+", tier="Assigned", ion_score=0.9, n_files=6, formula_agree=True),
        dict(mz=200.0, neutral_formula="C10H19NO2", adduct="[M+NH4]+", tier="Candidate", ion_score=0.7, n_files=3, formula_agree=False),
        dict(mz=183.0, neutral_formula="C9H10N2O", adduct="[M+H]+", tier="Candidate", ion_score=0.6, n_files=2, formula_agree=True),
        dict(mz=223.06, neutral_formula="C6H18O3Si3", adduct="[M+H]+", tier="Candidate", ion_score=0.5, n_files=1, formula_agree=True),
        dict(mz=247.0, neutral_formula="C3H2F6O", adduct="[M+Br]-", tier="Assigned", ion_score=0.8, n_files=4, formula_agree=True),
        dict(mz=400.0, neutral_formula="C9H12N4O12", adduct="[M+Na]+", tier="Candidate", ion_score=0.94, n_files=1, formula_agree=True),  # N-monster, flagged
    ]).to_csv(f"{d}/merged_ledger.csv", index=False)
    # one per-file ledger with roles (drives the role breakdown); the M0 row
    # carries a tier + ppm_error + peak_id so the mass-defect/mass-error QC figure
    # (qc_massdefect) has a calibrated point and a parented iso-child.
    pd.DataFrame([
        dict(peak_id=1, mz=169.1223, role="M0", neutral_formula="C10H16O2", height=10000,
             tier="Assigned", ppm_error=0.3, parent_peak_id=None),
        dict(peak_id=2, mz=170.1256, role="iso_child", neutral_formula=None, height=1100,
             tier=None, ppm_error=None, parent_peak_id=1),
        dict(peak_id=3, mz=250.0, role="reagent", neutral_formula=None, height=500,
             tier=None, ppm_error=None, parent_peak_id=None),
        dict(peak_id=4, mz=300.0, role="unexplained", neutral_formula=None, height=80,
             tier=None, ppm_error=None, parent_peak_id=None),
    ]).to_csv(f"{d}/per_file/s1_ledger.csv", index=False)
    pd.DataFrame({"neutral_formula": ["C10H16O2", "C10H14O2", "C9H12O2"],
                  "cluster": [1, 1, 1], "cv": [0.5, 0.6, 0.4],
                  "median_cps": [9000, 4000, 2000]}).to_csv(f"{d}/tables/clusters_changing_Ur.csv", index=False)

    ctx = R.load_context(d, tag="Ur", label="Ur⁺ CIMS")
    check("load_context: merged loaded", ctx["n_m0"] == 6, ctx.get("n_m0"))
    check("load_context: tiers counted", ctx["tiers"].get("Assigned") == 2, ctx["tiers"])
    check("load_context: composition is CHO/CHON/CHOS backbone",
          set(ctx["composition"]) <= {"CHO", "CHON", "CHOS"}, ctx["composition"])
    check("load_context: heteroatom side-counts Si + F",
          ctx["hetero"]["Si-bearing"] == 1 and ctx["hetero"]["F-bearing"] == 1,
          ctx.get("hetero"))
    check("load_context: per-adduct channel counts present", bool(ctx.get("adduct_counts")),
          ctx.get("adduct_counts"))
    check("load_context: role breakdown present", "unexplained" in ctx.get("role_count", {}),
          ctx.get("role_count"))
    check("load_context: brightest full per-file ledger retained for the QC figure",
          ctx.get("bright_ledger") is not None and "role" in ctx["bright_ledger"].columns,
          None if ctx.get("bright_ledger") is None else list(ctx["bright_ledger"].columns)[:3])
    check("load_context: role-signal split (analyte/reagent/unexplained)",
          set(ctx.get("role_signal_frac", {})) == {"analyte", "reagent", "unexplained"},
          ctx.get("role_signal_frac"))
    check("load_context: signal-weighted composition present",
          "CHO" in ctx.get("sig_comp_frac", {}), ctx.get("sig_comp_frac"))
    check("load_context: amine-shadow detected (C10H19NO2 = C10H16O2+NH3)",
          ctx.get("shadow", {}).get("n_shadowed") == 1, ctx.get("shadow"))
    check("load_context: two-way collapsed composition present",
          ctx.get("n_collapsed") == 1, ctx.get("n_collapsed"))
    check("load_context: single-source disagreements from formula_agree (=1)",
          ctx.get("n_disagree") == 1 and ctx.get("n_multifile") == 4,
          (ctx.get("n_disagree"), ctx.get("n_multifile")))
    check("load_context: top species by signal carries the bright CHO",
          ctx.get("top_species") and ctx["top_species"][0]["neutral_formula"] == "C10H16O2",
          ctx.get("top_species"))
    check("load_context: polarity detected positive (NH4/[M+H]+ present)",
          ctx.get("positive") is True, ctx.get("positive"))
    check("load_context: plausibility flags the Candidate N-monster only",
          [f["neutral_formula"] for f in ctx.get("flagged", [])] == ["C9H12N4O12"],
          ctx.get("flagged"))
    # findings + scrutiny sections build (degrade gracefully without an event TS)
    out_f = R.build(d, tag="Ur", label="Ur⁺ CIMS", out_pdf=f"{d}/rf.pdf",
                    sections=[R.findings, R.scrutiny])
    check("build: findings+scrutiny sections standalone OK",
          os.path.exists(out_f) and os.path.getsize(out_f) > 1500)

    # the mass-defect / mass-error QC section builds standalone + writes its figure
    out_q = R.build(d, tag="Ur", label="Ur⁺ CIMS", out_pdf=f"{d}/rq.pdf",
                    sections=[R.qc_massdefect])
    check("build: qc_massdefect section standalone OK",
          os.path.exists(out_q) and os.path.getsize(out_q) > 1500)
    check("build: qc_massdefect wrote the figure under figures/",
          os.path.exists(f"{d}/figures/massdefect_masserror_Ur.png"))

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

    # run versioning: run_id + a date+time 'generated' on the cover (title page)
    RID = "Orange-peeling-Ur-CIMS_2026-06-20_143512"
    ctx_r = R.load_context(d, tag="Ur", label="Ur⁺ CIMS", run_id=RID)
    check("load_context carries run_id for the cover", ctx_r.get("run_id") == RID)
    out3 = R.build(d, tag="Ur", label="Ur⁺ CIMS", generated="2026-06-20 14:35",
                   run_id=RID, out_pdf=f"{d}/r3.pdf")
    check("build with run_id + timestamped generated -> PDF", os.path.exists(out3))
    # default PDF filename carries the Report ID (self-identifying outside its folder)
    out4 = R.build(d, tag="Ur", label="Ur⁺ CIMS", run_id=RID)
    check("default PDF filename includes the Report ID", os.path.basename(out4) == f"report_{RID}.pdf",
          os.path.basename(out4))
    try:
        import fitz  # PyMuPDF — verify the cover text if available
        cover = fitz.open(out3)[0].get_text()
        check("cover shows the Report ID (with time)", RID in cover, cover[:400])
        check("cover 'generated' carries date AND time", "2026-06-20 14:35" in cover, cover[:400])
    except ImportError:
        pass

    # --- compress_pdf: optional size-reduced companion ---
    check("compress_pdf is a no-op when the input is already small (returns None)",
          R.compress_pdf(out4, min_mb=100.0) is None)
    try:
        import io  # noqa: F401
        import fitz
        from PIL import Image  # noqa: F401
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        big = os.path.join(d, "big.pdf")
        with PdfPages(big) as pp:
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.imshow(np.random.rand(1600, 1600, 3))     # large embedded raster
            pp.savefig(fig, dpi=200); plt.close(fig)
        small = R.compress_pdf(big, min_mb=0.0, max_px=200, quality=40)
        check("compress_pdf shrinks an image-heavy PDF",
              bool(small) and os.path.getsize(small) < os.path.getsize(big),
              f"{os.path.getsize(big)} -> {os.path.getsize(small) if small else None}")
        check("compressed PDF is valid + keeps the page count",
              bool(small) and open(small, "rb").read(4) == b"%PDF"
              and fitz.open(small).page_count == fitz.open(big).page_count)
    except ImportError:
        pass

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

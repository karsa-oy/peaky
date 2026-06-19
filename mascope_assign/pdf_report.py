"""Standard PDF report for a representative-batch assignment run.

Assembles the assignment findings into one PDF per batch: cover + headline,
coverage stats (assigned vs unassigned, by count AND signal), composition / full
Van Krevelen, analyte families, correlated-cluster TS figures, and methods.

ITERABLE BY DESIGN: the report is an ordered list of SECTIONS, each a function
`section(ctx, pdf)` that draws one or more pages. To change the report, edit /
reorder / add a section function and list it in SECTIONS — nothing else couples.
`ctx` is a dict loaded once by `load_context()` from the run's on-disk artifacts
(merged_ledger.csv, per_file/*, the figures, *_summary.json), so a section just
reads what it needs and degrades gracefully when an artifact is missing.

Pure matplotlib (PdfPages) — no new dependencies, no web stack. Build with
`build(out_dir, tag=..., label=..., ts_path=...)`.
"""
from __future__ import annotations

import glob
import json
import os

import numpy as np
import pandas as pd

__version__ = "0.1.0"

A4 = (8.27, 11.69)                       # portrait inches
INK = "#222222"
GREY = "#777777"


# ---------------------------------------------------------------------------
# context: load everything the report needs, once
# ---------------------------------------------------------------------------
def _skill_version() -> str:
    """Identify the skill build: assign version + git short SHA of this repo."""
    import subprocess
    d = os.path.dirname(os.path.abspath(__file__))
    sha = "?"
    try:
        sha = subprocess.check_output(["git", "-C", d, "rev-parse", "--short", "HEAD"],
                                      text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        pass
    try:
        from . import assign
        av = assign.__version__
    except Exception:
        av = "?"
    return f"mascope_assign (assign v{av}) · git {sha}"


def load_context(out_dir: str, *, tag: str, label: str, ts_path: str | None = None,
                 generated: str = "", batch_name: str | None = None) -> dict:
    from . import analyte_viz as V
    from . import chemistry as C
    out_dir = os.path.expanduser(out_dir)
    ctx: dict = {"out_dir": out_dir, "tag": tag, "label": label, "fig": {},
                 "generated": generated, "version": _skill_version(),
                 "batch_name": batch_name}

    merged = pd.read_csv(f"{out_dir}/merged_ledger.csv")
    ctx["merged"] = merged
    ctx["n_m0"] = len(merged)
    ctx["tiers"] = merged["tier"].value_counts().to_dict()
    u = merged.drop_duplicates("neutral_formula").copy()
    ctx["n_neutrals"] = len(u)
    # composition by CHO/CHON/CHOS backbone (Si/F/halogen folded in)
    ctx["composition"] = u["neutral_formula"].map(V.backbone_class).value_counts().to_dict()
    # heteroatom side-counts (additions to the backbone)
    cnt = u["neutral_formula"].map(lambda f: C.parse_formula(str(f)))
    ctx["hetero"] = {
        "Si-bearing (siloxane)": int(cnt.map(lambda c: c.get("Si", 0) > 0).sum()),
        "F-bearing": int(cnt.map(lambda c: c.get("F", 0) > 0).sum()),
        "Cl/Br in neutral": int(cnt.map(lambda c: c.get("Cl", 0) + c.get("Br", 0) > 0).sum()),
    }

    # per-file pooled role breakdown (count + signal) + the explained m/z set
    rows = []
    for f in sorted(glob.glob(f"{out_dir}/per_file/*_ledger.csv")):
        d = pd.read_csv(f)
        d["h"] = pd.to_numeric(d.get("height"), errors="coerce").fillna(0)
        rows.append(d)
    ctx["n_files"] = len(rows)
    if rows:
        a = pd.concat(rows, ignore_index=True)
        ctx["role_count"] = a["role"].value_counts().to_dict()
        ctx["role_signal"] = {k: float(v) for k, v in a.groupby("role")["h"].sum().items()}
        ctx["expl_mz"] = np.sort(a.loc[a["role"].astype(str) != "unexplained",
                                       "mz"].dropna().to_numpy())

    # whole-spectrum coverage (TS bins explained vs unexplained, count + signal)
    if ts_path and os.path.exists(os.path.expanduser(ts_path)):
        from . import timeseries as TS
        ts = pd.read_parquet(os.path.expanduser(ts_path))
        if not ctx.get("batch_name") and "sample_batch_name" in ts.columns:
            names = ts["sample_batch_name"].dropna().unique()
            if len(names):
                ctx["batch_name"] = str(names[0])
        mat, bin_mz = TS.build_matrix(ts)
        binsig = mat.sum(axis=0)
        expl = ctx.get("expl_mz", np.array([]))

        def matched(mz, tol=8.0):
            i = np.searchsorted(expl, mz)
            return any(0 <= j < len(expl) and abs(expl[j] - mz) / mz * 1e6 <= tol
                       for j in (i - 1, i))
        ex = np.array([matched(float(bin_mz[b])) for b in mat.columns])
        ctx["ts"] = {"nbins": int(len(ex)), "expl_count": int(ex.sum()),
                     "expl_signal": float(binsig.values[ex].sum()),
                     "tot_signal": float(binsig.sum())}

    for key, fn in [("jitter", "jitter_summary.json"), ("batch", "batch_summary.json")]:
        p = f"{out_dir}/{fn}"
        if os.path.exists(p):
            ctx[key] = json.load(open(p))
    vk = f"{out_dir}/van_krevelen_full_{tag}.png"
    if os.path.exists(vk):
        ctx["fig"]["vk"] = vk           # single (scatter)
    # cluster figures are PAGED (clusters_<set>_<tag>_p<i>.png) — collect ALL pages
    for key, stem in [("changing", f"clusters_changing_{tag}"),
                      ("flat", f"clusters_flat_{tag}"),
                      ("unassigned", f"clusters_unassigned_{tag}")]:
        paged = sorted(glob.glob(f"{out_dir}/{stem}_p*.png"),
                       key=lambda s: int(s.rsplit("_p", 1)[1].split(".")[0]))
        if paged:
            ctx["fig"][key] = paged
        elif os.path.exists(f"{out_dir}/{stem}.png"):
            ctx["fig"][key] = [f"{out_dir}/{stem}.png"]
    cc = f"{out_dir}/clusters_changing_{tag}.csv"
    if os.path.exists(cc):
        ctx["changing_csv"] = pd.read_csv(cc)
    return ctx


# ---------------------------------------------------------------------------
# page helpers
# ---------------------------------------------------------------------------
def _close(pdf, fig):
    import matplotlib.pyplot as plt
    pdf.savefig(fig)
    plt.close(fig)


def _text_lines(fig, lines, *, x=0.08, y0=0.90, dy=0.026, size=10):
    """Render (style, text) lines top-down. style: 'h' head / 'b' body / 'm' mono."""
    y = y0
    for style, txt in lines:
        if style == "gap":
            y -= dy * (txt or 1)
            continue
        kw = dict(fontsize=size, color=INK, va="top", family="sans-serif")
        if style == "h":
            kw.update(fontsize=size + 3, weight="bold")
        elif style == "m":
            kw.update(family="monospace", fontsize=size - 0.5)
        elif style == "dim":
            kw.update(color=GREY, fontsize=size - 1)
        fig.text(x, y, txt, **kw)
        y -= dy
    return y


def _image_page(pdf, png, title, *, landscape=False, dpi=200, native=False, src_dpi=170):
    """Embed a PNG as one page. native=True makes the PAGE the image's own size
    (page = pixels/src_dpi inches) so a tall cluster figure embeds 1:1 and stays
    fully legible (the page is tall; scroll). Otherwise fit to A4 (landscape opt)."""
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    img = mpimg.imread(png)
    if native:
        ih, iw = img.shape[0] / src_dpi, img.shape[1] / src_dpi
        Th = 0.5                                    # title strip (inches)
        fig = plt.figure(figsize=(iw, ih + Th), dpi=src_dpi)
        fig.text(0.01, 1 - 0.20 / (ih + Th), title, fontsize=12, weight="bold", color=INK, va="top")
        ax = fig.add_axes([0.005, 0.002, 0.99, ih / (ih + Th)])
        ax.imshow(img); ax.axis("off")
        pdf.savefig(fig); plt.close(fig)
        return
    figsize = (A4[1], A4[0]) if landscape else A4
    fig = plt.figure(figsize=figsize, dpi=dpi)
    fig.text(0.04, 0.975, title, fontsize=13, weight="bold", color=INK)
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.93])
    ax.imshow(img); ax.axis("off")
    pdf.savefig(fig, dpi=dpi)
    plt.close(fig)


def _pct(part, whole):
    return 100.0 * part / whole if whole else 0.0


# ---------------------------------------------------------------------------
# sections  (each draws page(s); add/reorder freely)
# ---------------------------------------------------------------------------
def cover(ctx, pdf):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=A4)
    batch = ctx.get("batch_name") or ctx["label"]
    fig.text(0.08, 0.93, "Peak Assignment Report", fontsize=20, weight="bold", color=INK)
    fig.text(0.08, 0.895, batch, fontsize=14, color=INK)
    fig.text(0.08, 0.872, f"{ctx['label']} · representative-sample pipeline", fontsize=11,
             color=GREY)
    meta = ctx.get("version", "")
    if ctx.get("generated"):
        meta = f"{meta}  ·  generated {ctx['generated']}" if meta else f"generated {ctx['generated']}"
    if meta:
        fig.text(0.08, 0.852, meta, fontsize=8.5, color=GREY)

    tiers = ctx["tiers"]; idn = tiers.get("Identified", 0); cn = tiers.get("Candidate", 0)
    ts = ctx.get("ts"); sig = (_pct(ts["expl_signal"], ts["tot_signal"]) if ts else None)
    ex_c = (_pct(ts["expl_count"], ts["nbins"]) if ts else None)
    head = [
        ("h", "Headline"),
        ("gap", 0.3),
        ("b", f"Unique analytes assigned (M0):   {ctx['n_m0']}   "
              f"({idn} Identified / {cn} Candidate)"),
        ("b", f"Distinct neutral compounds:       {ctx['n_neutrals']}"),
    ]
    if ts:
        head += [
            ("b", f"Spectral coverage:                {ex_c:.0f}% of m/z bins, "
                  f"{sig:.0f}% of signal explained"),
            ("b", f"Unexplained:                      {ts['nbins']-ts['expl_count']} bins "
                  f"({100-ex_c:.0f}% count, {100-sig:.0f}% signal)"),
        ]
    b = ctx.get("batch", {})
    sids = b.get("sample_ids", [])
    head += [
        ("gap", 1),
        ("h", "Representative samples assigned"),
        ("gap", 0.3),
        ("b", f"{ctx['n_files']} files: 5 evenly time-spaced + the max-TIC sample, "
              "merged by m/z."),
    ]
    for s in sids[:8]:
        head.append(("m", f"   {s}"))
    j = ctx.get("jitter", {})
    if j:
        head += [
            ("gap", 1),
            ("h", "File-to-file reproducibility"),
            ("gap", 0.3),
            ("b", f"Per-file calibration spread:  {j.get('offset_spread_ppm','?')} ppm"),
            ("b", f"Mass jitter (median / p95):   {j.get('mz_jitter_raw_median','?')} / "
                  f"{j.get('mz_jitter_raw_p95','?')} ppm  (≈ genuine peak noise)"),
            ("b", f"Formula disagreements:        {j.get('formula_disagreements','?')}  ·  "
                  f"tier-unstable: {j.get('tier_unstable','?')}"),
        ]
    _text_lines(fig, head, y0=0.80, dy=0.030)
    _close(pdf, fig)


def coverage(ctx, pdf):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.965, "Coverage — assigned vs unassigned", fontsize=15,
             weight="bold", color=INK)
    fig.text(0.08, 0.94, "Read both by count and by signal: the unexplained pool is "
             "large by count but small by signal.", fontsize=9, color=GREY)

    rc = ctx.get("role_count", {}); rs = ctx.get("role_signal", {})
    roles = ["M0", "iso_child", "reagent", "artifact", "unexplained"]
    cols = {"M0": "#1D9E75", "iso_child": "#9AD1BE", "reagent": "#378ADD",
            "artifact": "#BBBBBB", "unexplained": "#D85A30"}
    Nc = sum(rc.values()) or 1; Ns = sum(rs.values()) or 1

    # (a) role breakdown — count vs signal stacked bars
    ax = fig.add_axes([0.10, 0.66, 0.82, 0.22])
    left_c = left_s = 0
    for r in roles:
        c = _pct(rc.get(r, 0), Nc); s = _pct(rs.get(r, 0), Ns)
        ax.barh(1, c, left=left_c, color=cols[r], edgecolor="white")
        ax.barh(0, s, left=left_s, color=cols[r], edgecolor="white",
                label=f"{r} ({rc.get(r,0)})")
        if c > 5:
            ax.text(left_c + c / 2, 1, f"{c:.0f}%", ha="center", va="center", fontsize=7.5, color="white")
        if s > 5:
            ax.text(left_s + s / 2, 0, f"{s:.0f}%", ha="center", va="center", fontsize=7.5, color="white")
        left_c += c; left_s += s
    ax.set_yticks([0, 1]); ax.set_yticklabels(["by signal", "by count"]); ax.set_xlim(0, 100)
    ax.set_xlabel("% of peak-rows / signal (pooled over the representative files)")
    ax.set_title("Role breakdown", loc="left", fontsize=11)
    ax.legend(fontsize=7.5, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.22),
              frameon=False, columnspacing=1.0, handletextpad=0.4)

    # (b) M0 tier split
    ax2 = fig.add_axes([0.10, 0.33, 0.36, 0.15])
    t = ctx["tiers"]; idn = t.get("Identified", 0); cn = t.get("Candidate", 0)
    ax2.bar(["Identified", "Candidate"], [idn, cn], color=["#1D9E75", "#E0A93B"])
    for i, v in enumerate([idn, cn]):
        ax2.text(i, v, str(v), ha="center", va="bottom", fontsize=9)
    ax2.set_title(f"Confidence of the {ctx['n_m0']} assignments", loc="left", fontsize=11)
    ax2.set_ylabel("M0 count")

    # (c) spectral coverage donut-ish bar
    ts = ctx.get("ts")
    ax3 = fig.add_axes([0.56, 0.33, 0.36, 0.15])
    if ts:
        ec = _pct(ts["expl_count"], ts["nbins"]); es = _pct(ts["expl_signal"], ts["tot_signal"])
        ax3.barh(["by signal", "by count"], [es, ec], color="#1D9E75", label="explained")
        ax3.barh(["by signal", "by count"], [100 - es, 100 - ec], left=[es, ec],
                 color="#D85A30", label="unexplained")
        for yi, val in zip([1, 0], [ec, es]):
            ax3.text(val, yi, f" {val:.0f}%", va="center", fontsize=8)
        ax3.set_xlim(0, 100); ax3.set_title("Spectral coverage (TS bins)", loc="left", fontsize=11)
        ax3.legend(fontsize=7.5, loc="lower right", frameon=False)
    else:
        ax3.axis("off"); ax3.text(0, 0.5, "(TS not loaded)", color=GREY)

    # text takeaways
    lines = [("h", "What it means"), ("gap", 0.3)]
    if ts:
        lines.append(("b", f"• {_pct(ts['expl_signal'], ts['tot_signal']):.0f}% of the ion "
                      f"signal is explained; the unexplained {ts['nbins']-ts['expl_count']} "
                      "bins are mostly small, near-noise peaks."))
    lines += [
        ("b", f"• {ctx['tiers'].get('Identified',0)} of {ctx['n_m0']} assignments are "
              f"high-confidence (Identified); the rest Candidate."),
        ("b", "• 'explained' = M0 compounds + their isotope satellites + reagent + "
              "artifact peaks (not just M0)."),
    ]
    _text_lines(fig, lines, y0=0.22, dy=0.028)
    _close(pdf, fig)


def composition(ctx, pdf):
    if "vk" in ctx["fig"]:
        _image_page(pdf, ctx["fig"]["vk"],
                    f"{ctx['label']} — Van Krevelen (all assigned peaks, by composition)")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.93, "Composition of the assigned peaks", fontsize=15, weight="bold", color=INK)
    comp = ctx.get("composition", {})
    het = ctx.get("hetero", {})
    lines = [("h", f"Distinct neutral compounds by backbone ({ctx['n_neutrals']} total)"),
             ("gap", 0.3)]
    for kl in ("CHO", "CHON", "CHOS"):
        if kl in comp:
            lines.append(("m", f"   {kl:6s} {comp[kl]}"))
    lines += [("gap", 1),
              ("h", "Heteroatom additions (within the backbone classes above)"),
              ("gap", 0.3)]
    for k, v in het.items():
        lines.append(("m", f"   {k:24s} {v}"))
    lines += [("gap", 1),
              ("dim", "Si/F/halogen are folded into the CHO/CHON/CHOS backbone, not split out"),
              ("dim", "(a siloxane with no N is CHO; a fluorinated species with N is CHON)."),
              ("dim", "Si = PDMS/silicone inlet bleed; F/halogen are reagent/contaminant ladders.")]
    _text_lines(fig, lines, y0=0.86, dy=0.030)
    _close(pdf, fig)


def families(ctx, pdf):
    for p in ctx["fig"].get("changing", []):
        _image_page(pdf, p, f"{ctx['label']} — co-varying analyte families (changing peaks)",
                    native=True)
    cc = ctx.get("changing_csv")
    if cc is None or not len(cc):
        return
    import matplotlib.pyplot as plt
    from . import chemistry as C
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.93, "Analyte families (temporal clusters)", fontsize=15, weight="bold", color=INK)
    sizes = cc.groupby("cluster").size().sort_values(ascending=False)
    lines = [("h", "Largest co-varying clusters"), ("gap", 0.3),
             ("m", "  cluster   n   median O/C   example members"),
             ("gap", 0.2)]
    for cid in sizes.index[:12]:
        g = cc[cc.cluster == cid]
        ocs = []
        for f in g["neutral_formula"]:
            cnt = C.parse_formula(str(f)); nc = cnt.get("C", 0)
            if nc:
                ocs.append(cnt.get("O", 0) / nc)
        oc = np.median(ocs) if ocs else float("nan")
        ex = ", ".join(g.sort_values("median_cps", ascending=False)["neutral_formula"].head(4))
        lines.append(("m", f"  {int(cid):>5}   {len(g):>2}     {oc:>5.2f}      {ex}"))
    _text_lines(fig, lines, y0=0.86, dy=0.026, size=9)
    _close(pdf, fig)


def clusters(ctx, pdf):
    for p in ctx["fig"].get("flat", []):
        _image_page(pdf, p, f"{ctx['label']} — flat background (RAW intensity bands)", native=True)
    for p in ctx["fig"].get("unassigned", []):
        _image_page(pdf, p, f"{ctx['label']} — unexplained-peak clusters (RAW intensity)", native=True)


def methods(ctx, pdf):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.93, "Methods & caveats", fontsize=15, weight="bold", color=INK)
    lines = [
        ("h", "Pipeline"), ("gap", 0.3),
        ("b", "• Representative-sample rule: assign 5 evenly time-spaced samples + the"),
        ("b", "  max-TIC sample per batch, then merge by m/z — a single averaged file"),
        ("b", "  misses analytes present only part of the run."),
        ("b", "• Each file: multi-pass formula assignment, server isotope-scored matching"),
        ("b", "  (match_compounds), isotope-envelope completion, calibrated tiering."),
        ("b", "• Clusters: TIC/reagent-normalised log-correlation, complete-linkage at r>0.6"),
        ("b", "  on the full-batch time series (signed distance keeps anti-phase apart)."),
        ("gap", 1),
        ("h", "Caveats"), ("gap", 0.3),
        ("b", "• Report coverage both by count and by signal — they differ a lot."),
        ("b", "• 'Unexplained' is dominated by small near-noise peaks (low signal)."),
        ("b", "• Confidence (Identified/Candidate) applies to M0 compounds only."),
        ("b", "• Mass-degenerate peaks can flip formula across files (see jitter)."),
    ]
    _text_lines(fig, lines, y0=0.86, dy=0.028)
    if ctx.get("generated"):
        fig.text(0.08, 0.06, f"mascope_assign · generated {ctx['generated']}", fontsize=8, color=GREY)
    _close(pdf, fig)


SECTIONS = [cover, coverage, composition, families, clusters, methods]


def build(out_dir: str, *, tag: str, label: str, ts_path: str | None = None,
          out_pdf: str | None = None, generated: str = "", batch_name: str | None = None,
          sections=SECTIONS) -> str:
    """Build the PDF report for one batch run. `out_dir` holds the run artifacts.
    `batch_name` titles the report (else taken from the TS, else the reagent label)."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_pdf import PdfPages
    ctx = load_context(out_dir, tag=tag, label=label, ts_path=ts_path,
                       generated=generated, batch_name=batch_name)
    out_pdf = out_pdf or os.path.join(os.path.expanduser(out_dir), f"report_{tag}.pdf")
    with PdfPages(out_pdf) as pdf:
        for section in sections:
            try:
                section(ctx, pdf)
            except Exception as e:                       # one bad section must not kill the report
                import matplotlib.pyplot as plt
                fig = plt.figure(figsize=A4)
                fig.text(0.08, 0.9, f"[section '{section.__name__}' failed: {e}]",
                         fontsize=10, color="#B00020", wrap=True)
                _close(pdf, fig)
    return out_pdf

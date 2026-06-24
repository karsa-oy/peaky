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

from . import paths as PT

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
    return f"peaky (assign v{av}) · git {sha}"


def load_context(out_dir: str, *, tag: str, label: str, ts_path: str | None = None,
                 generated: str = "", batch_name: str | None = None,
                 run_id: str | None = None) -> dict:
    from . import analyte_viz as V
    from . import chemistry as C
    out_dir = os.path.expanduser(out_dir)
    RP = PT.run_paths(out_dir)
    FIG, TAB = RP.figures, RP.tables       # MUST mirror the writers (clustering/analyte_viz)
    ctx: dict = {"out_dir": out_dir, "fig_dir": FIG, "tag": tag, "label": label,
                 "fig": {}, "generated": generated, "version": _skill_version(),
                 "batch_name": batch_name, "run_id": run_id}

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
        "Si-bearing": int(cnt.map(lambda c: c.get("Si", 0) > 0).sum()),
        "F-bearing": int(cnt.map(lambda c: c.get("F", 0) > 0).sum()),
        "Cl/Br in neutral": int(cnt.map(lambda c: c.get("Cl", 0) + c.get("Br", 0) > 0).sum()),
    }
    # match score per tier + adduct-channel breakdown (assignment-quality page)
    if "ion_score" in merged.columns:
        gs = merged.groupby("tier")["ion_score"].mean()
        ctx["score_by_tier"] = {t: float(gs[t]) for t in gs.index}
    ctx["adduct_counts"] = merged["adduct"].value_counts().to_dict()   # actual channels
    # the representative sample NAMES (timestamps), not just ids
    ss = f"{TAB}/selected_samples.csv"
    if os.path.exists(ss):
        s = pd.read_csv(ss)
        name = s["sample_item_name"] if "sample_item_name" in s.columns else s["sample_item_id"]
        ctx["samples"] = list(zip(name.astype(str), s.get("role", pd.Series([""] * len(s))).astype(str)))

    # per-file pooled role breakdown (count + signal) + the explained m/z set
    rows, rep_ids = [], []
    for f in sorted(glob.glob(f"{out_dir}/per_file/*_ledger.csv")):
        rep_ids.append(os.path.basename(f).replace("_ledger.csv", ""))
        d = pd.read_csv(f)
        d["h"] = pd.to_numeric(d.get("height"), errors="coerce").fillna(0)
        rows.append(d)
    ctx["n_files"] = len(rows)
    ctx["rep_ids"] = rep_ids
    if rows:
        a = pd.concat(rows, ignore_index=True)
        ctx["role_count"] = a["role"].value_counts().to_dict()
        ctx["role_signal"] = {k: float(v) for k, v in a.groupby("role")["h"].sum().items()}
        # signal share split into analyte (M0+iso) / reagent / unexplained — for an
        # honest "explained" headline: a Br- spectrum is mostly the reagent ion, so
        # "98% explained" must not read as "98% analyte characterised".
        rs = ctx["role_signal"]; rtot = sum(rs.values()) or 1.0
        ctx["role_signal_frac"] = {
            "analyte": (rs.get("M0", 0) + rs.get("iso_child", 0)) / rtot,
            "reagent": rs.get("reagent", 0) / rtot,
            "unexplained": (rs.get("unexplained", 0) + rs.get("artifact", 0)) / rtot}
        # per-neutral summed M0 signal (drives signal-weighted composition + findings)
        if "neutral_formula" in a.columns:
            m0a = a[a["role"] == "M0"]
            ns = m0a.groupby(m0a["neutral_formula"].astype(str))["h"].sum()
            ctx["neutral_signal"] = {k: float(v) for k, v in ns.items()
                                     if k and k != "nan"}
        ctx["expl_mz"] = np.sort(a.loc[a["role"].astype(str) != "unexplained",
                                       "mz"].dropna().to_numpy())
        # mass error (ppm) by category for the accuracy plot
        if "ppm_error" in a.columns:
            m0 = a[a["role"] == "M0"]
            pbc = {}
            for t in ("Identified", "Candidate"):
                v = m0[m0.get("tier") == t]["ppm_error"].dropna().tolist()
                if v:
                    pbc[t] = v
            iso = _iso_ppm(rows)
            if iso:
                pbc["isotopologue"] = iso
            ctx["ppm_by_cat"] = pbc
        if "adduct" in a.columns:                     # signal share per ion channel
            asig = a[a["role"] == "M0"].groupby("adduct")["h"].sum()
            tot = float(asig.sum()) or 1.0
            ctx["adduct_signal"] = {k: float(v) / tot for k, v in asig.items()}

    # signal-weighted composition + ammonium/amine degeneracy (composition page)
    from . import composition as CMP
    nsig = ctx.get("neutral_signal", {})
    ctx["sig_comp_frac"], ctx["sig_comp_abs"] = CMP.signal_by_backbone(merged, nsig)
    ctx["shadow"] = CMP.amine_shadow_stats(merged)
    ctx["comp_asg"], ctx["comp_collapsed"], ctx["n_collapsed"] = CMP.collapsed_composition(merged)
    ctx["top_species"] = CMP.top_species_by_signal(merged, nsig, n=8)
    ctx["oligomers"] = CMP.oligomer_flag(merged)
    # polarity (gates positive-only messaging: the amine re-read, the shadow note)
    # + chemical-plausibility QC of the assignments
    ctx["positive"] = any(str(k).rstrip().endswith("+") for k in ctx.get("adduct_counts", {}))
    from . import plausibility as PL
    ctx["flagged"] = PL.scan(merged, polarity=("+" if ctx["positive"] else "-"))
    # single source of truth for "formula disagreements": the merged ledger's own
    # formula_agree column (the artifact the report publishes), with a denominator
    # — not the separate jitter pass, which clusters differently and prints a
    # slightly different number (the 65-vs-68 mismatch).
    if "formula_agree" in merged.columns:
        ctx["n_disagree"] = int((~merged["formula_agree"].astype(bool)).sum())
    if "n_files" in merged.columns:
        ctx["n_multifile"] = int((merged["n_files"] >= 2).sum())

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
        # event overview: per-sample total signal vs wall-clock over the FULL batch
        # (the reader needs to SEE the burst; the cluster pages only show a
        # normalised 0-1 'hour' axis). Mark where the representative files sit.
        try:
            tt = pd.to_datetime(ts["datetime_utc"], utc=True)
            tstamp = tt.groupby(ts["sample_item_id"]).first()
            order = tstamp.sort_values().index
            t0 = tstamp.min()
            hrs = (tstamp.reindex(order) - t0).dt.total_seconds().to_numpy() / 3600.0
            tot = mat.reindex(order).sum(axis=1).to_numpy()
            rep_h = [float((tstamp[r] - t0).total_seconds() / 3600.0)
                     for r in ctx.get("rep_ids", []) if r in tstamp.index]
            ctx["event"] = {"hours": hrs, "total": tot, "rep_hours": rep_h}
        except Exception:
            pass

    for key, fn in [("jitter", "jitter_summary.json"), ("batch", "batch_summary.json")]:
        p = f"{out_dir}/{fn}"
        if os.path.exists(p):
            ctx[key] = json.load(open(p))
    vk = f"{FIG}/van_krevelen_full_{tag}.png"
    if os.path.exists(vk):
        ctx["fig"]["vk"] = vk           # single (scatter)
    # cluster figures are PAGED (clusters_<set>_<tag>_p<i>.png) — collect ALL pages
    for key, stem in [("changing", f"clusters_changing_{tag}"),
                      ("changers", f"clusters_changers_{tag}"),
                      ("flat", f"clusters_flat_{tag}"),
                      ("unassigned", f"clusters_unassigned_{tag}")]:
        paged = sorted(glob.glob(f"{FIG}/{stem}_p*.png"),
                       key=lambda s: int(s.rsplit("_p", 1)[1].split(".")[0]))
        if paged:
            ctx["fig"][key] = paged
        elif os.path.exists(f"{FIG}/{stem}.png"):
            ctx["fig"][key] = [f"{FIG}/{stem}.png"]
    cc = f"{TAB}/clusters_changing_{tag}.csv"
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
        Th = 0.5 if title else 0.0                  # title strip (inches); the figures
        fig = plt.figure(figsize=(iw, ih + Th), dpi=src_dpi)   # carry their own title
        if title:
            fig.text(0.01, 1 - 0.20 / (ih + Th), title, fontsize=12, weight="bold",
                     color=INK, va="top")
        ax = fig.add_axes([0.005, 0.002, 0.99, ih / (ih + Th)])
        ax.imshow(img); ax.axis("off")
        pdf.savefig(fig); plt.close(fig)
        return
    figsize = (A4[1], A4[0]) if landscape else A4
    fig = plt.figure(figsize=figsize, dpi=dpi)
    if title:
        fig.text(0.04, 0.975, title, fontsize=13, weight="bold", color=INK)
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.93] if title else [0.02, 0.02, 0.96, 0.96])
    ax.imshow(img); ax.axis("off")
    pdf.savefig(fig, dpi=dpi)
    plt.close(fig)


def _pct(part, whole):
    return 100.0 * part / whole if whole else 0.0


# readable descriptor for each ion channel (so the breakdown names the actual adduct)
_ADDUCT_DESC = {
    "[M+H]+": "protonated", "[M-H]-": "deprotonated",
    "[M+Br]-": "Br- cluster", "[M+HBr+Br]-": "di-bromide cluster", "[M+Br2]-": "Br2 cluster",
    "[M+CO3]-": "CO3- cluster", "[M+HBr+CO3]-": "HBr.CO3- cluster", "[M+HSO4]-": "HSO4- cluster",
    "[M+(CH4N2O)H]+": "urea cluster", "[M+Na]+": "Na+ adduct", "[M+NH4]+": "NH4+ adduct",
    "[M+Cl]-": "Cl- cluster", "[M+I]-": "I- cluster",
    "[M+NO3]-": "NO3- cluster", "[M+^NO3]-": "¹⁵NO3- cluster",
}


def _adduct_label(a) -> str:
    a = str(a)
    d = _ADDUCT_DESC.get(a)
    return f"{a}  ({d})" if d else a


_ISO_C13 = 1.0033548   # 13C - 12C; used to recover isotopologue mass error


def _iso_ppm(per_file_frames):
    """Isotopologue mass error (ppm): per file, match each M0's predicted 13C M+1
    (parent ion m/z + 1.00335) to the nearest iso_child peak. iso_child rows carry
    no ppm_error of their own, so recompute it here."""
    import numpy as np
    out = []
    for d in per_file_frames:
        if "role" not in d.columns or "mz" not in d.columns:
            continue
        iso = np.sort(d[d["role"] == "iso_child"]["mz"].dropna().to_numpy())
        if not len(iso):
            continue
        for mz in d[d["role"] == "M0"]["mz"].dropna():
            tgt = float(mz) + _ISO_C13
            j = np.searchsorted(iso, tgt)
            for k in (j - 1, j):
                if 0 <= k < len(iso):
                    p = (iso[k] - tgt) / tgt * 1e6
                    if abs(p) <= 6:
                        out.append(p); break
    return out


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
    if ctx.get("run_id"):
        fig.text(0.08, 0.834, f"Report ID:  {ctx['run_id']}", fontsize=8.5, color=GREY)

    tiers = ctx["tiers"]; idn = tiers.get("Identified", 0); cn = tiers.get("Candidate", 0)
    ts = ctx.get("ts"); sig = (_pct(ts["expl_signal"], ts["tot_signal"]) if ts else None)
    ex_c = (_pct(ts["expl_count"], ts["nbins"]) if ts else None)
    head = [
        ("h", "Summary"),
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
        rf = ctx.get("role_signal_frac", {})
        if rf.get("reagent", 0) >= 0.05:        # reagent-dominated (e.g. Br-): be honest
            head += [("dim", "   of detected signal:  analyte (M0+iso) "
                             f"{rf['analyte']*100:.0f}%,  reagent ion {rf['reagent']*100:.0f}%,  "
                             f"unexplained {rf['unexplained']*100:.0f}%")]
    head += [
        ("gap", 1),
        ("h", "Representative samples assigned"),
        ("gap", 0.3),
        ("b", f"{ctx['n_files']} files: 5 evenly time-spaced + the max-TIC sample, "
              "merged by m/z."),
    ]
    for name, role in ctx.get("samples", [])[:8]:
        head.append(("m", f"   {name}   [{role}]"))
    j = ctx.get("jitter", {})
    if j:
        nm = ctx.get("n_multifile")
        nd = ctx.get("n_disagree", j.get("formula_disagreements", "?"))
        tu = j.get("tier_unstable", "?")
        dis = (f"{nd} of {nm} multi-file ({_pct(nd, nm):.0f}%)"
               if nm and isinstance(nd, int) else f"{nd}")
        tus = (f"{tu} of {nm} multi-file ({_pct(tu, nm):.0f}%)"
               if nm and isinstance(tu, int) else f"{tu}")
        head += [
            ("gap", 1),
            ("h", "File-to-file reproducibility"),
            ("gap", 0.3),
            ("b", f"Per-file calibration spread:  {j.get('offset_spread_ppm','?')} ppm"),
            ("b", f"Mass jitter (median / p95):   {j.get('mz_jitter_raw_median','?')} / "
                  f"{j.get('mz_jitter_raw_p95','?')} ppm  (≈ genuine peak noise)"),
            ("b", f"Formula disagreements:        {dis}  (same m/z, different formula)"),
            ("b", f"Tier-unstable assignments:    {tus}  (flip Identified <-> Candidate)"),
            ("dim", "On a disagreement/flip the merge keeps the highest tier, then the "
                    "highest match score."),
        ]
    _text_lines(fig, head, y0=0.80, dy=0.029)
    _close(pdf, fig)


def findings(ctx, pdf):
    """Plain-language findings page (right after the cover): the event time-trace
    plus data-driven takeaways — top species by signal, signal-weighted
    composition, and the oligomer/HOM fingerprint. Answers 'what changed during
    the run?' before the QC detail. Everything is derived from the data, so it
    works for any batch."""
    import matplotlib.pyplot as plt
    ev = ctx.get("event"); top = ctx.get("top_species", [])
    scf = ctx.get("sig_comp_frac", {}); olig = ctx.get("oligomers", [])
    if not (ev or top):
        return
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.955, "Findings", fontsize=16, weight="bold", color=INK)
    rise_txt = None
    if ev and len(ev.get("total", [])):
        ax = fig.add_axes([0.10, 0.60, 0.82, 0.27])
        h = np.asarray(ev["hours"], float); tt = np.asarray(ev["total"], float)
        ax.plot(h, tt, color="#1D9E75", lw=1.5)
        ax.fill_between(h, tt, color="#1D9E75", alpha=0.12)
        for rh in ev.get("rep_hours", []):
            ax.axvline(rh, color="#888", lw=0.7, ls=":")
        ax.set_xlabel("hour of experiment (UTC)", fontsize=9)
        ax.set_ylabel("total signal (cps)", fontsize=9)
        ax.set_title("Event overview — total signal vs time (dotted = assigned samples)",
                     loc="left", fontsize=10.5)
        ax.grid(alpha=0.3)
        ok = np.isfinite(tt)
        if ok.sum() >= 3:
            tail = tt[ok][max(1, int(ok.sum() * 2 / 3)):]
            base = float(np.median(tail)) if len(tail) else float(np.median(tt[ok]))
            peak = float(np.nanmax(tt)); pk_h = float(h[int(np.nanargmax(tt))])
            if base > 0:
                rise_txt = (f"Total signal peaks at hour {pk_h:.1f} — {peak/base:.1f}x the "
                            f"late-run baseline — then decays (a transient event).")
    lines = []
    if rise_txt:
        lines += [("b", "• " + rise_txt)]
    if scf:
        cc = ctx.get("comp_asg", {}); nn = ctx.get("n_neutrals", 1)
        cho = scf.get("CHO", 0) * 100; chon = scf.get("CHON", 0) * 100
        lines += [("b", f"• By signal the assigned chemistry is {cho:.0f}% CHO / {chon:.0f}% CHON"
                        f" — vs {_pct(cc.get('CHON', 0), nn):.0f}% CHON by compound count")]
        if ctx.get("positive"):     # the CHON count inflation is the amine re-read (positive only)
            lines += [("b", "  (the count is inflated by mass-degenerate ammonium/amine re-reads; "
                            "see Composition).")]
        else:
            lines += [("b", "  (a few bright CHO species carry most of the signal).")]
    if top:
        lines += [("gap", 0.6), ("h", "Top species by signal"), ("gap", 0.25),
                  ("m", "   share   class   neutral")]
        for r in top[:8]:
            lines.append(("m", f"   {r['frac']*100:>4.1f}%   {r['klass']:5s}   {r['neutral_formula']}"))
    if olig:
        nsig = ctx.get("neutral_signal", {})
        olig = sorted(olig, key=lambda f: nsig.get(f, 0.0), reverse=True)[:12]
        lines += [("gap", 0.6), ("h", "Accretion / oligomer products (high C & O, by signal)"),
                  ("gap", 0.25)]
        for k in range(0, len(olig), 6):           # wrap ~6 formulas per line (no edge clip)
            lines.append(("m", "   " + ", ".join(olig[k:k + 6])))
        lines += [("dim", "high-carbon high-oxygen neutrals — candidate HOM dimers / oligomers,"),
                  ("dim", "often the most event-specific signal.")]
    _text_lines(fig, lines, y0=0.52, dy=0.027, size=9.5)
    _close(pdf, fig)


def coverage(ctx, pdf):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.965, "Assignment quality", fontsize=15, weight="bold", color=INK)

    # (a) mean match score, Identified vs Candidate
    sbt = ctx.get("score_by_tier", {})
    tiers = [t for t in ("Identified", "Candidate") if t in sbt]
    ax = fig.add_axes([0.11, 0.74, 0.33, 0.17])
    ax.bar(tiers, [sbt[t] for t in tiers], color=["#1D9E75", "#E0A93B"][:len(tiers)], width=0.55)
    for i, t in enumerate(tiers):
        ax.text(i, sbt[t], f"{sbt[t]:.2f}", ha="center", va="bottom", fontsize=10)
    ax.set_ylim(0, 1.05); ax.set_ylabel("mean match score", fontsize=9)
    ax.set_title("Match score by confidence", loc="left", fontsize=11)

    # (b) mass error (ppm) by category: Identified / Candidate / isotopologues
    pbc = ctx.get("ppm_by_cat", {})
    ax2 = fig.add_axes([0.58, 0.74, 0.34, 0.17])
    cats = [c for c in ("Identified", "Candidate", "isotopologue") if c in pbc]
    if cats:
        bp = ax2.boxplot([pbc[c] for c in cats], vert=True, showfliers=False, widths=0.5,
                         patch_artist=True, medianprops=dict(color="#111"))
        for patch, c in zip(bp["boxes"], ["#1D9E75", "#E0A93B", "#9AD1BE"]):
            patch.set_facecolor(c)
        ax2.set_xticklabels([f"{c}\n(n={len(pbc[c])})" for c in cats], fontsize=8)
        ax2.axhline(0, color="0.7", lw=0.7)
    ax2.set_ylabel("mass error (ppm)", fontsize=9)
    ax2.set_title("Mass accuracy", loc="left", fontsize=11)

    # (c) signal & peak share by ROLE — analyte (M0+iso) / reagent ion / unexplained.
    # Splitting out the reagent is the honest "explained" picture: a Br- spectrum is
    # mostly the reagent ion, so a single "98% explained" number is misleading.
    rc = ctx.get("role_count", {}); rs = ctx.get("role_signal", {})
    def _roles(d):
        tot = sum(d.values()) or 1.0
        return (100 * (d.get("M0", 0) + d.get("iso_child", 0)) / tot,
                100 * d.get("reagent", 0) / tot,
                100 * (d.get("unexplained", 0) + d.get("artifact", 0)) / tot)
    ax3 = fig.add_axes([0.11, 0.47, 0.33, 0.17])
    if rc or rs:
        for y, d in zip([1, 0], [rc, rs]):
            an, rg, ot = _roles(d)
            ax3.barh(y, an, color="#1D9E75")
            ax3.barh(y, rg, left=an, color="#378ADD")
            ax3.barh(y, ot, left=an + rg, color="#D85A30")
            for v, x in ((an, an / 2), (rg, an + rg / 2), (ot, an + rg + ot / 2)):
                if v >= 7:
                    ax3.text(x, y, f"{v:.0f}", va="center", ha="center", fontsize=7, color="white")
        ax3.set_yticks([1, 0]); ax3.set_yticklabels(["by count", "by signal"], fontsize=9)
        ax3.set_xlim(0, 100)
        ax3.set_xlabel("% — green analyte · blue reagent · red unexpl.", fontsize=7.3)
    ax3.set_title("Signal & peaks by role", loc="left", fontsize=11)

    # (d) assignments by ACTUAL ion channel (which adduct each peak was assigned on)
    adc = ctx.get("adduct_counts", {}); asig = ctx.get("adduct_signal", {})
    items = sorted(adc.items(), key=lambda kv: kv[1], reverse=True)
    ax4 = fig.add_axes([0.58, 0.47, 0.34, 0.17])
    y = list(range(len(items)))[::-1]
    ax4.barh(y, [v for _, v in items], color="#378ADD")
    ax4.set_yticks(y); ax4.set_yticklabels([k for k, _ in items], fontsize=8, family="monospace")
    xmax = max([v for _, v in items], default=1)
    for yi, (k, v) in zip(y, items):
        st = f" · {asig[k] * 100:.0f}% sig" if k in asig else ""
        ax4.text(v + xmax * 0.02, yi, f"{v} · {_ADDUCT_DESC.get(k, '')}{st}", va="center", fontsize=6.8)
    ax4.set_xlim(0, xmax * 1.9); ax4.set_xlabel("M0 count (a compound can appear in several channels)", fontsize=7.5)
    ax4.set_title("Reagent / ion channels assigned on", loc="left", fontsize=11)

    idn = ctx["tiers"].get("Identified", 0); cn = ctx["tiers"].get("Candidate", 0)
    lines = [("h", "Reading this page"), ("gap", 0.3),
             ("b", f"• {idn} Identified vs {cn} Candidate. Match score = the server "
                   "isotope-scored compound match (0-1)."),
             ("b", "• Mass accuracy: ppm error of the matched peaks (boxes = IQR, line = median; "
                   "near 0 = well calibrated)."),
             ("b", "• Ion channel names the adduct each compound was assigned on. It counts a "
                   "compound once per adduct, so bright species recur across channels — read the"),
             ("b", "  signal % on the bars (e.g. protonation usually dominates the signal)."),
             ("b", "• Unexplained peaks are a third of the m/z bins but only a few % of the signal "
                   "(dim, near-noise)."),
             ("b", "• 'By role' splits explained signal into analyte (M0 + isotopes) vs the reagent "
                   "ion vs unexplained —"),
             ("b", "  so a reagent-dominated negative-mode spectrum isn't read as fully characterised.")]
    # the NH4->amine caveat only applies to positive urea-CIMS (where NH4 adducts
    # exist); never print it on a negative-mode (e.g. Br-) report.
    pos = any(("NH4" in str(k)) or ("CH4N2O" in str(k)) for k in ctx.get("adduct_counts", {}))
    if pos:
        lines += [("dim", "[M+NH4]+ is mass/isotope-identical to [M+H]+ of the +NH3 amine -- kept as NH4 only"),
                  ("dim", "when its trace co-varies (r>=0.7) with the protonated/urea parent, else re-read as the"),
                  ("dim", "amine (a parsimony prior, not a measurement). Peak roles -> Methods page.")]
    else:
        lines += [("dim", "Peak roles are defined on the Methods page.")]
    _text_lines(fig, lines, y0=0.36, dy=0.029)
    _close(pdf, fig)


def composition(ctx, pdf):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=A4)            # text page FIRST, then the VK figure
    fig.text(0.08, 0.95, "Composition of the assigned peaks", fontsize=15, weight="bold", color=INK)
    comp = ctx.get("composition", {})
    het = ctx.get("hetero", {})
    scf = ctx.get("sig_comp_frac", {})
    lines = [("h", f"Distinct neutral compounds by backbone ({ctx['n_neutrals']} total)"),
             ("gap", 0.3),
             ("dim", "by COUNT — each compound once, regardless of how bright it is")]
    for kl in ("CHO", "CHON", "CHOS"):
        if kl in comp:
            lines.append(("m", f"   {kl:6s} {comp[kl]:>4}   ({_pct(comp[kl], ctx['n_neutrals']):.0f}%)"))
    if scf:
        lines += [("gap", 0.8),
                  ("h", "Same compounds, weighted by signal"),
                  ("gap", 0.3),
                  ("dim", "where the chemistry actually is — a few bright species carry most signal")]
        for kl in ("CHO", "CHON", "CHOS"):
            if kl in scf:
                lines.append(("m", f"   {kl:6s} {scf[kl]*100:>3.0f}% of assigned signal"))
    sh = ctx.get("shadow", {}); coll = ctx.get("comp_collapsed", {})
    if sh.get("n_shadowed") and ctx.get("positive"):    # the re-read is positive urea-CIMS only
        lines += [("gap", 0.8),
                  ("dim", f"Note: {sh['n_shadowed']} CHON neutrals share an exact NH3-shifted CHO twin that is"),
                  ("dim", "also assigned — mass-degenerate [M+NH4]+/[M+H]+ pairs from the amine re-read,"),
                  ("dim", f"counted twice. Read as their CHO parent: {coll.get('CHO','?')} CHO / "
                          f"{coll.get('CHON','?')} CHON ({sh['collapsed_neutrals']} distinct).")]
    lines += [("gap", 0.8),
              ("h", "Heteroatom additions (within the backbone classes above)"),
              ("gap", 0.3)]
    for k, v in het.items():
        lines.append(("m", f"   {k:24s} {v}"))
    present = "/".join(k for k in ("CHO", "CHON", "CHOS") if k in comp)
    lines += [("gap", 0.8),
              ("dim", f"Si/F/halogen are folded into the {present} backbone, not split out"),
              ("dim", "(a siloxane with no N is CHO; a fluorinated species with N is CHON)."),
              ("dim", "Si = PDMS/silicone inlet bleed; F/halogen are reagent/contaminant ladders.")]
    _text_lines(fig, lines, y0=0.89, dy=0.028)
    _close(pdf, fig)
    if "vk" in ctx["fig"]:
        _image_page(pdf, ctx["fig"]["vk"], "")          # figure carries its own title


def scrutiny(ctx, pdf):
    """Assignments flagged for chemical-plausibility scrutiny (Candidate-tier only):
    high-heteroatom mass coincidences, implausibly carbon-rich skeletons, or a
    halogen in a positive-mode neutral. Flagged, NOT removed — listed so a reader
    can discount them. Renders only when something is flagged."""
    fl = ctx.get("flagged", [])
    if not fl:
        return
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.95, "Assignments flagged for scrutiny", fontsize=15, weight="bold", color=INK)
    lines = [("dim", f"{len(fl)} Candidate-tier neutral(s) whose formula looks more like a mass"),
             ("dim", "coincidence than a real molecule (many heteroatoms, very low H/C, or a"),
             ("dim", "wrong-mode halogen). Flagged for review, NOT removed; Identified (isotope-"),
             ("dim", "scored) assignments are never flagged."),
             ("gap", 0.8),
             ("m", "   score   neutral                  why"),
             ("gap", 0.2)]
    for d in fl[:40]:
        sc = f"{d['ion_score']:.2f}" if d.get("ion_score") is not None else "  -  "
        lines.append(("m", f"   {sc}   {d['neutral_formula']:22s}  {d['reason']}"))
    if len(fl) > 40:
        lines.append(("dim", f"… and {len(fl) - 40} more."))
    _text_lines(fig, lines, y0=0.88, dy=0.025, size=8.5)
    _close(pdf, fig)


def gka(ctx, pdf):
    """GKA homologous-series findings — a small-multiple grid of Kendrick
    mass-defect plots, one per repeat-unit family (alkyl / oxidation / alkoxylate
    / siloxane / fluorinated), each rotated to flatten its own series into
    horizontal ladders. The print counterpart of the interactive rotating-GKA
    widget. Rendered on demand from the merged ledger (neutral_formula only)."""
    from . import gka_figure as GK
    merged = ctx.get("merged")
    if merged is None or "neutral_formula" not in getattr(merged, "columns", []) \
            or not len(merged):
        return
    fig_dir = ctx.get("fig_dir") or ctx["out_dir"]
    os.makedirs(fig_dir, exist_ok=True)
    png = os.path.join(fig_dir, f"gka_{ctx['tag']}.png")
    GK.render_gka(merged, png, title=ctx.get("batch_name") or ctx["label"])
    ctx["fig"]["gka"] = png
    _image_page(pdf, png, "")           # figure carries its own title


def families(ctx, pdf):
    for p in ctx["fig"].get("changing", []):
        _image_page(pdf, p, "")                         # A4-portrait page, own title
    cc = ctx.get("changing_csv")
    if cc is None or not len(cc):
        return
    import matplotlib.pyplot as plt
    from . import chemistry as C
    fig = plt.figure(figsize=A4)
    fig.text(0.08, 0.93, "Analyte families (temporal clusters)", fontsize=15, weight="bold", color=INK)
    sizes = cc.groupby("cluster").size().sort_values(ascending=False)
    big = [c for c in sizes.index if sizes[c] >= 3]      # only the PLOTTED clusters
    singletons = int((sizes < 3).sum())
    lines = [("h", f"Co-varying clusters ({len(big)} with >=3 members)"), ("gap", 0.3),
             ("dim", "members = ion channels (a neutral's [M+H]+ / cluster / adduct ions are clustered "
                     "separately — they often don't co-vary; see the per-cluster workbook for the breakdown)"),
             ("gap", 0.3),
             ("m", "  cluster   n   median O/C   top neutral formulas (by intensity)"),
             ("gap", 0.2)]
    for cid in big:
        g = cc[cc.cluster == cid]
        ocs = [cnt.get("O", 0) / cnt.get("C", 1) for cnt in
               (C.parse_formula(str(f)) for f in g["neutral_formula"]) if cnt.get("C", 0)]
        oc = np.median(ocs) if ocs else float("nan")
        top = ", ".join(g.sort_values("median_cps", ascending=False)["neutral_formula"]
                        .drop_duplicates().head(4))
        lines.append(("m", f"  {int(cid):>5}   {len(g):>2}     {oc:>5.2f}      {top}"))
    if singletons:
        lines += [("gap", 0.5),
                  ("dim", f"+ {singletons} singleton / <3-member peaks — shown together in the "
                          "last 'remaining peaks' panel of the changing-cluster figure.")]
    _text_lines(fig, lines, y0=0.86, dy=0.026, size=9)
    _close(pdf, fig)


def changers(ctx, pdf):
    """Large standalone changes: single channels that change a lot (>= fold) with no
    family — pulled out so they're not buried in the flat panel (user-requested).
    A4-portrait pages (own title), embedded like the other cluster figures."""
    for p in ctx["fig"].get("changers", []):
        _image_page(pdf, p, "")             # fit-to-A4 (the PNG is already A4 portrait)


def clusters(ctx, pdf):
    for p in ctx["fig"].get("flat", []) + ctx["fig"].get("unassigned", []):
        _image_page(pdf, p, "")                         # A4-portrait pages, own titles


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
        ("b", "• Clusters: log-correlation (raw or reagent-normalised) of the full-batch"),
        ("b", "  time series, complete-linkage at r>0.6 (signed distance keeps anti-phase apart)."),
        ("gap", 0.5),
        ("dim", f"Parameters: m/z merge tol {ctx.get('batch', {}).get('tol_ppm', 6.0)} ppm · "
                "cluster r>0.6 · amine co-variation r>=0.7."),
        ("gap", 1),
        ("h", "Peak roles"), ("gap", 0.3),
        ("b", "• M0 = an assigned compound's monoisotopic peak (the identification)."),
        ("b", "• iso_child = its isotope satellite (13C / 81Br / 37Cl ...)."),
        ("b", "• reagent = reagent-ion / reagent-cluster peaks (e.g. Br3-)."),
        ("b", "• artifact = instrument peaks (FT ringing / sidelobes of a bright peak)."),
        ("b", "• unexplained = no confident formula."),
        ("gap", 1),
        ("h", "Caveats"), ("gap", 0.3),
        ("b", "• Report coverage both by count and by signal — they differ a lot."),
        ("b", "• 'Unexplained' is dominated by small near-noise peaks (low signal)."),
        ("b", "• Confidence (Identified/Candidate) applies to M0 compounds only."),
        ("b", "• Mass-degenerate peaks can flip formula across files (see jitter)."),
    ]
    rf = ctx.get("role_signal_frac", {})
    if rf.get("reagent", 0) >= 0.05:
        lines += [("b", f"• The reagent ion carries ~{rf['reagent']*100:.0f}% of the signal; "
                        "'explained signal'"),
                  ("b", "  is split into analyte vs reagent on the quality page.")]
    if any(("NH4" in str(k)) or ("CH4N2O" in str(k)) for k in ctx.get("adduct_counts", {})):
        lines += [("b", "• Positive urea-CIMS: [M+NH4]+ adducts are mass/isotope-identical to the"),
                  ("b", "  protonated +NH3 amine; uncorroborated ones are re-read as the amine"),
                  ("b", "  (co-variation r>=0.7), which raises the CHON count — a parsimony prior,"),
                  ("b", "  not a measurement (see Composition for the degeneracy it leaves).")]
    _text_lines(fig, lines, y0=0.88, dy=0.026)
    if ctx.get("generated"):
        fig.text(0.08, 0.06, f"peaky · generated {ctx['generated']}", fontsize=8, color=GREY)
    _close(pdf, fig)


SECTIONS = [cover, findings, coverage, composition, scrutiny, gka, families, changers, clusters, methods]


def build(out_dir: str, *, tag: str, label: str, ts_path: str | None = None,
          out_pdf: str | None = None, generated: str = "", batch_name: str | None = None,
          run_id: str | None = None, sections=SECTIONS) -> str:
    """Build the PDF report for one batch run. `out_dir` holds the run artifacts.
    `batch_name` titles the report (else taken from the TS, else the reagent label).
    `run_id` (the timestamped run folder name) is stamped on the cover as the Report ID."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_pdf import PdfPages
    ctx = load_context(out_dir, tag=tag, label=label, ts_path=ts_path,
                       generated=generated, batch_name=batch_name, run_id=run_id)
    if out_pdf is None:
        # name the PDF with the Report ID when we have one, so the file is
        # self-identifying even when moved out of its run folder.
        fname = f"report_{run_id}.pdf" if run_id else f"report_{tag}.pdf"
        out_pdf = os.path.join(PT.run_paths(out_dir).ensure().report, fname)
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


def compress_pdf(in_pdf: str, out_pdf: str | None = None, *, max_px: int = 850,
                 quality: int = 58, min_mb: float = 2.0, log=print) -> str | None:
    """Write a size-reduced COMPANION of a report PDF, for emailing/sharing.

    The bulk of a report is the embedded figure rasters; this downsamples each to
    `max_px` on its long edge and re-encodes it as JPEG (`quality`), while leaving
    the page TEXT vector (so formulas/labels stay crisp). The original `in_pdf` is
    never modified — a sibling `<name>_compressed.pdf` is written and its path
    returned. Returns None (no-op) when: PyMuPDF/Pillow are not installed (optional
    deps — `pip install 'mascope-assign[compress]'`), the input is already under
    `min_mb`, or anything goes wrong. Deliberately kept out of `build()` so the
    primary report stays byte-for-byte deterministic (see test_determinism)."""
    if out_pdf is None:
        out_pdf = (in_pdf[:-4] if in_pdf.lower().endswith(".pdf") else in_pdf) + "_compressed.pdf"
    try:
        if os.path.getsize(in_pdf) < min_mb * 1e6:
            return None
    except OSError:
        return None
    try:
        import io
        import fitz                      # PyMuPDF (optional)
        from PIL import Image            # Pillow (optional)
    except Exception:
        log("[report] compress skipped (install 'mascope-assign[compress]' for PyMuPDF+Pillow)")
        return None
    try:
        doc = fitz.open(in_pdf)
        seen: set = set()
        for page in doc:
            for img in page.get_images(full=True):
                xref = img[0]
                if xref in seen:
                    continue
                seen.add(xref)
                try:
                    raw = doc.extract_image(xref)
                    im = Image.open(io.BytesIO(raw["image"]))
                    w, h = im.size
                    if max(w, h) > max_px:
                        s = max_px / max(w, h)
                        im = im.resize((max(1, int(w * s)), max(1, int(h * s))), Image.LANCZOS)
                    buf = io.BytesIO()
                    im.convert("RGB").save(buf, "JPEG", quality=quality, optimize=True)
                    page.replace_image(xref, stream=buf.getvalue())
                except Exception:
                    pass                  # one bad image must not abort compression
        doc.save(out_pdf, garbage=4, deflate=True, clean=True)
        doc.close()
        return out_pdf
    except Exception as e:                # never let compression break a run
        log(f"[report] compress failed ({e}); keeping the full report only")
        return None

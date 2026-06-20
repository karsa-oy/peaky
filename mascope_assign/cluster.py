"""Correlation clustering of a batch time-series into co-varying families.

Consolidates the proven scratch clustering (cluster_analytes / flat_raw /
cluster_unassigned) into one parameterised module. Three peak SETS run through
the SAME engine:

  * changing    assigned analytes (M0, organic, Si-excluded) with cv >= CHANGING
  * flat        assigned background (cv < CHANGING, Si-INCLUSIVE) — clustered on
                RAW intensity (these resolve into intensity bands, not chemistry)
  * unassigned  TS bins not matched to any assigned M0 (the 'unexplained' set)

Engine: per-item time traces (analyte_viz.time_traces / a TS bin matrix) ->
normalise (TIC or reagent) for `changing`/`unassigned`, RAW for `flat` ->
log10 -> Pearson correlation -> distance (1 - r, SIGNED so anti-phase stays
apart) -> COMPLETE-linkage hierarchical clustering, cut at r > (1 - DIST_T).

COMPLETE linkage is essential — average linkage chains every trace into one blob
(seen on Br). Pure numpy/scipy/pandas; lazy matplotlib renderer. Reagent params
(adducts, normaliser) come from a profiles.ReagentProfile.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

__version__ = "0.1.0"

CHANGING = 0.30          # cv at/above which a trace is "changing"
FLAT_CV = CHANGING       # cv BELOW which a trace is "flat / non-varying" and is
                         # bunched out of correlation clustering (the gate knob —
                         # lower it to keep weak-but-coherent families as clusters)
DIST_T = 0.40            # 1-r cut: r > 0.6 merges
MIN_POINTS = 8           # finite trace points needed to correlate
MIN_MEMBERS = 3          # smallest reported cluster
LINK = "complete"
PALETTE = ["#1D9E75", "#7F77DD", "#D85A30", "#378ADD", "#BA7517",
           "#D4537E", "#639922", "#888780", "#0F6E56", "#534AB7"]


# ---------------------------------------------------------------------------
# pure clustering core
# ---------------------------------------------------------------------------
def smooth(y, w=5):
    y = np.asarray(y, float)
    if np.all(np.isnan(y)):
        return y
    f = np.where(np.isnan(y), np.nanmedian(y), y)
    return np.convolve(f, np.ones(w) / w, mode="same") if len(f) >= w else f


def shape_of(mean_z, gap=0.5):
    s, e = np.nanmean(mean_z[:6]), np.nanmean(mean_z[-6:])
    return "rise" if e - s > gap else ("fall" if s - e > gap else "peak")


def correlate(traces_for_corr: pd.DataFrame, cols, *, min_points=MIN_POINTS):
    """log10 -> Pearson corr matrix over `cols` of the (already normalised or raw)
    trace frame. Returns (Lg z-input, corr DataFrame)."""
    cols = [c for c in cols if c in traces_for_corr.columns]
    if not cols:
        return pd.DataFrame(), pd.DataFrame()
    vals = traces_for_corr[cols].values
    floor = np.nanmin(vals[vals > 0]) if np.any(vals > 0) else 1e-9
    Lg = np.log10(traces_for_corr[cols].clip(lower=floor))
    return Lg, Lg.corr(min_periods=min_points)


def cluster(cm: pd.DataFrame, *, dist_t=DIST_T, link=LINK, min_members=MIN_MEMBERS):
    """Complete-linkage clustering of the correlation matrix. Returns
    (labels Series indexed by column, ordered list of clusters >= min_members).
    Degrades gracefully when there are too few items to cluster."""
    cols = list(cm.columns)
    if len(cols) < 2:
        return pd.Series({c: 1 for c in cols}, dtype=int), []
    dist = (1 - cm.fillna(0)).values.copy()
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2
    Z = linkage(squareform(dist, checks=False), method=link)
    lab = pd.Series(fcluster(Z, t=dist_t, criterion="distance"), index=cols)
    sizes = lab.value_counts()
    big = [int(c) for c in sizes.index if sizes[c] >= min_members]
    return lab, big


def trace_cv(traces: pd.DataFrame, col) -> float:
    """Coefficient of variation (std/mean) of a trace over its finite positive
    points. 0.0 if the trace is empty/flat-at-zero. This is the 'does it vary'
    metric — the same one the driver uses to split changing vs flat analytes."""
    if col not in traces:
        return 0.0
    tr = pd.to_numeric(traces[col], errors="coerce").to_numpy()
    tr = tr[np.isfinite(tr) & (tr > 0)]
    m = tr.mean() if len(tr) else 0.0
    return float(tr.std() / m) if m > 0 else 0.0


def split_varying(traces: pd.DataFrame, cols, *, cv_min=FLAT_CV):
    """Partition `cols` into (varying, flat) by trace CV. A flat trace (cv < cv_min)
    has no reliable SHAPE — its pairwise correlation with anything is dominated by
    noise, so flat traces spuriously shatter into many tiny clusters. Pull them out
    BEFORE correlation clustering and bunch them into one 'flat / non-varying' group
    instead. Returns (varying, flat) preserving input order, deduped."""
    varying, flat = [], []
    for c in dict.fromkeys(cols):
        (varying if trace_cv(traces, c) >= cv_min else flat).append(c)
    return varying, flat


def threshold_scan(cm: pd.DataFrame, ts=(0.3, 0.4, 0.5, 0.6, 0.7), link=LINK):
    """Diagnostic: cluster sizes (>=3) vs the 1-r cut, to pick DIST_T."""
    if len(cm.columns) < 2:
        return {}
    dist = (1 - cm.fillna(0)).values.copy()
    np.fill_diagonal(dist, 0.0); dist = (dist + dist.T) / 2
    Z = linkage(squareform(dist, checks=False), method=link)
    out = {}
    for t in ts:
        s = pd.Series(fcluster(Z, t=t, criterion="distance")).value_counts()
        out[t] = sorted([int(x) for x in s if x >= 3], reverse=True)
    return out


def cluster_rows(cols, lab, big, cm, traces_raw, grid, *, median_cps=None, maxc=None):
    """Build per-cluster summary rows (sorted by size) for rendering:
    (cid, members[brightest-first], r_bar, shape, peak_hour). maxc=None -> ALL
    clusters (the default — dissect all signal); pass an int to cap."""
    Lg = np.log10(traces_raw[cols].clip(lower=1e-9))
    Z = (Lg - Lg.mean()) / Lg.std()
    if median_cps is None:
        median_cps = {c: float(np.nanmedian(traces_raw[c].values)) for c in cols}
    rows = []
    for cid in (big if maxc is None else big[:maxc]):
        mem = sorted([c for c in cols if lab[c] == cid],
                     key=lambda c: -float(median_cps.get(c, 0)))
        mz = smooth(Z[mem].mean(axis=1).values)
        ph = float(grid[int(np.nanargmax(mz))])
        sub = cm.loc[mem, mem].values
        iu = np.triu_indices(len(sub), 1)
        rbar = float(np.nanmean(sub[iu])) if iu[0].size else float("nan")
        rows.append((cid, mem, rbar, shape_of(mz), ph))
    return rows, Z


# ---------------------------------------------------------------------------
# renderer (stacked panels: member traces + black median + item labels under)
# ---------------------------------------------------------------------------
def _wrap(toks, width=118):
    lines, cur = [], ""
    for t in toks:
        add = t if not cur else ", " + t
        if len(cur) + len(add) > width:
            lines.append(cur); cur = t
        else:
            cur += add
    if cur:
        lines.append(cur)
    return lines


def render_paged(rows, grid, traces_z, traces_raw, item_label, out_prefix, *,
                 per_page=10, **kw):
    """Render `rows` clusters across multiple pages of <= per_page each, so ALL
    clusters are shown while each panel stays legible. Writes <out_prefix>_p<i>.png
    and returns the list of paths. kw passes through to render_panels (mode/ylim)."""
    if not rows:
        return []
    title = kw.pop("title", "")
    pages = [rows[i:i + per_page] for i in range(0, len(rows), per_page)]
    paths = []
    for i, chunk in enumerate(pages, 1):
        t = f"{title}  (page {i}/{len(pages)})" if len(pages) > 1 else title
        out = f"{out_prefix}_p{i}.png"
        render_panels(chunk, grid, traces_z, traces_raw, item_label, out, title=t, **kw)
        paths.append(out)
    return paths


def render_panels(rows, grid, traces_z, traces_raw, item_label, out, *,
                  mode="raw", title="", event_span=None, ylim=None, labels=True):
    """Stacked per-cluster panels. mode='z' (z-scored shape, mean) or 'raw'
    (log cps, median). `item_label` maps a column -> its printed label
    (formula or m/z). `event_span` shades an (h0,h1) window. labels=False drops
    the per-member label block (for a dense 'remaining peaks' overview panel)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not rows:
        return None
    Wf, Lm, Rm, Tm, Bm = 11.0, 0.95, 0.5, 0.5, 0.25
    plot_h, lh, gap = 1.9, 0.205, 0.5
    wrapped = ([_wrap([item_label(m) for m in mem], width=88) for _, mem, _, _, _ in rows]
               if labels else [[] for _ in rows])
    heights = [(plot_h, 0.14 + len(w) * lh) for w in wrapped]
    Hf = Tm + Bm + sum(ph + th + gap for ph, th in heights)
    fig = plt.figure(figsize=(Wf, Hf))
    desc = "z-scored shape; black = mean" if mode == "z" else "raw cps, log y; black = median"
    fig.suptitle(f"{title}  ({desc})", fontsize=11.5, y=1 - 0.18 / Hf)
    lf, wf = Lm / Wf, (Wf - Lm - Rm) / Wf
    cur = Hf - Tm
    for k, ((cid, mem, rbar, sh, ph), w, (h0, th)) in enumerate(zip(rows, wrapped, heights)):
        ax = fig.add_axes([lf, (cur - h0) / Hf, wf, h0 / Hf])
        if event_span:
            ax.axvspan(*event_span, color="#D85A30", alpha=0.06)
        col = PALETTE[k % len(PALETTE)]
        a = 0.6 if len(mem) <= 12 else (0.4 if len(mem) <= 30 else 0.25)
        M = []
        for m in mem:
            y = (traces_z[m].values if mode == "z" else traces_raw[m].values).astype(float)
            yy = y if mode == "z" else np.where(y > 0, y, np.nan)
            ax.plot(grid, yy, color=col, lw=0.8, alpha=a)
            M.append(yy)
        with np.errstate(all="ignore"):
            med = smooth(np.nanmean(np.array(M), axis=0)) if mode == "z" \
                else np.nanmedian(np.array(M), axis=0)
        ax.plot(grid, med, color="#111", lw=2.4, alpha=0.9, zorder=5)
        if mode == "z":
            ax.set_ylim(-3.3, 2.7); ax.set_ylabel("z")
        else:
            ax.set_yscale("log")
            if ylim:
                ax.set_ylim(*ylim)
            ax.set_ylabel("cps")
        ax.set_xlim(0, float(grid[-1])); ax.grid(alpha=0.18, which="both")
        ax.tick_params(labelsize=10)
        rbar_s = f"r̄={rbar:.2f} · " if rbar == rbar else ""        # skip for the remaining panel
        head = (f"{cid} · n={len(mem)}" if isinstance(cid, str)
                else f"cluster {cid} · n={len(mem)}")
        ax.set_title(f"{head} · {rbar_s}{sh} (peak~h{ph:.1f})", fontsize=11.5, loc="left")
        ax.set_xlabel("hour of experiment (UTC)", fontsize=10) if k == len(rows) - 1 else ax.set_xticklabels([])
        if labels:
            tx = fig.add_axes([lf, (cur - h0 - th) / Hf, wf, th / Hf]); tx.axis("off")
            tx.text(0, 1, "\n".join(w), va="top", ha="left", fontsize=9.5,
                    family="monospace", color="#333", transform=tx.transAxes)
        cur -= (h0 + th + gap)
    fig.savefig(out, dpi=170, bbox_inches="tight"); plt.close(fig)
    return out


def render_a4(rows, grid, traces_z, traces_raw, item_label, out_prefix, *,
              mode="raw", ylim=None, title="", subtitle="", labels=True, base_color=0):
    """Render cluster panels onto A4 PORTRAIT pages: each panel spans the A4 text
    width, panels are packed top-down by height and a new page starts when the page
    is full (a panel never straddles a page). A clear GAP separates each trace from
    its formula/m-z list. Returns the list of A4 page PNG paths."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not rows:
        return []
    PAGE_W, PAGE_H = 8.27, 11.69
    L, R, BOTM = 0.62, 0.30, 0.40
    TOPM = 1.0 if subtitle else 0.70          # extra room for a subtitle/legend key
    UW, UH = PAGE_W - L - R, PAGE_H - TOPM - BOTM
    TRACE_H, GAP_TL, LINE_H, PANEL_GAP, WRAP = 1.30, 0.40, 0.165, 0.55, 96
    XLAB_H = 0.26          # extra clearance for the last panel's "hour..." x-axis label

    pan = []
    for row in rows:
        w = _wrap([item_label(m) for m in row[1]], width=WRAP) if labels else []
        lab_h = (0.05 + len(w) * LINE_H) if w else 0.0
        pan.append((row, w, TRACE_H + (GAP_TL + lab_h if w else 0.0), lab_h))

    pages, cur, curh = [], [], 0.0
    for p in pan:
        if cur and curh + PANEL_GAP + p[2] > UH:
            pages.append(cur); cur, curh = [], 0.0
        curh += (PANEL_GAP if cur else 0.0) + p[2]
        cur.append(p)
    if cur:
        pages.append(cur)

    paths, idx = [], base_color
    for pi, page in enumerate(pages, 1):
        fig = plt.figure(figsize=(PAGE_W, PAGE_H))
        ty = (PAGE_H - 0.38) / PAGE_H
        fig.text(L / PAGE_W, ty, title, fontsize=11, weight="bold", color="#222")
        if subtitle:
            fig.text(L / PAGE_W, ty - 0.165 / PAGE_H, subtitle, fontsize=8.5, color="#666")
        if len(pages) > 1:                            # page number RIGHT-aligned so it never clips
            fig.text(1 - R / PAGE_W, ty, f"page {pi}/{len(pages)}", fontsize=10,
                     color="#777", ha="right")
        ytop = PAGE_H - TOPM
        for k, (row, w, h, lab_h) in enumerate(page):
            cid, mem, rbar, sh, ph = row
            ax = fig.add_axes([L / PAGE_W, (ytop - TRACE_H) / PAGE_H, UW / PAGE_W, TRACE_H / PAGE_H])
            col = PALETTE[idx % len(PALETTE)]; idx += 1
            al = 0.6 if len(mem) <= 12 else (0.4 if len(mem) <= 30 else 0.25)
            M = []
            for m in mem:
                y = (traces_z[m].values if mode == "z" else traces_raw[m].values).astype(float)
                yy = y if mode == "z" else np.where(y > 0, y, np.nan)
                ax.plot(grid, yy, color=col, lw=0.8, alpha=al); M.append(yy)
            with np.errstate(all="ignore"):
                med = smooth(np.nanmean(np.array(M), axis=0)) if mode == "z" \
                    else np.nanmedian(np.array(M), axis=0)
            ax.plot(grid, med, color="#111", lw=2.2, alpha=0.9, zorder=5)
            if mode == "z":
                ax.set_ylim(-3.3, 2.7); ax.set_ylabel("z", fontsize=9)
            else:
                ax.set_yscale("log"); ax.set_ylabel("cps", fontsize=9)
                if ylim:
                    ax.set_ylim(*ylim)
            ax.set_xlim(0, float(grid[-1])); ax.grid(alpha=0.18, which="both"); ax.tick_params(labelsize=9)
            rbar_s = f"r̄={rbar:.2f} · " if rbar == rbar else ""
            head = f"{cid} · n={len(mem)}" if isinstance(cid, str) else f"cluster {cid} · n={len(mem)}"
            ax.set_title(f"{head} · {rbar_s}{sh} (peak~h{ph:.1f})", fontsize=10.5, loc="left")
            last = (k == len(page) - 1)
            if last:
                ax.set_xlabel("hour of experiment (UTC)", fontsize=9)
            if w:
                # the last panel carries the x-axis label in the gap -> push its
                # formula list down by XLAB_H so the two never collide
                ty = ytop - TRACE_H - GAP_TL - (XLAB_H if last else 0.0)
                tx = fig.add_axes([L / PAGE_W, (ty - lab_h) / PAGE_H, UW / PAGE_W, lab_h / PAGE_H])
                tx.axis("off")
                tx.text(0, 1, "\n".join(w), va="top", ha="left", fontsize=8.8,
                        family="monospace", color="#333", transform=tx.transAxes)
            ytop -= (h + PANEL_GAP)
        out = f"{out_prefix}_p{pi}.png"
        fig.savefig(out, dpi=200); plt.close(fig)
        paths.append(out)
    return paths


def remaining_row(cols, lab, big, traces_raw, grid):
    """A single overview 'cluster' of every peak NOT in a >=MIN_MEMBERS cluster
    (singletons + <3-member groups) so ALL signal is plotted. Returns (row, members)
    or (None, [])."""
    bigset = set(int(b) for b in big)
    rem = [c for c in cols if int(lab.get(c, -1)) not in bigset]
    if not rem:
        return None, []
    Lg = np.log10(traces_raw[rem].clip(lower=1e-9))
    Z = (Lg - Lg.mean()) / Lg.std()
    mz = smooth(Z.mean(axis=1).values)
    ph = float(grid[int(np.nanargmax(mz))])
    return (f"remaining {len(rem)} peaks", rem, float("nan"), shape_of(mz), ph), rem


def _sheet_name(cid, used: set) -> str:
    """A valid, unique Excel sheet name (<=31 chars, no []:*?/\\)."""
    base = (f"c{int(cid)}" if not isinstance(cid, str)
            else "".join(ch for ch in cid if ch not in "[]:*?/\\")[:28]) or "cluster"
    name = base[:31]
    k = 1
    while name in used:
        suf = f"_{k}"; name = base[:31 - len(suf)] + suf; k += 1
    used.add(name)
    return name


def write_cluster_workbook(rows, out_xlsx, *, meta=None, item_label=None,
                           member_cols=None, title="clusters"):
    """Per-cluster Excel workbook: a 'summary' sheet (one row per cluster: id, n,
    shape, peak hour, mean r) + ONE sheet PER CLUSTER listing its members. `rows`
    is the cluster_rows output [(cid, members, rbar, shape, peak_hr), ...]; `meta`
    maps a member key -> a dict of columns (e.g. neutral_formula / channel / m_z /
    match_score / tier). `member_cols` orders/selects those columns. Returns the
    path (or None if no rows)."""
    if not rows:
        return None
    summary = pd.DataFrame([{
        "cluster": (cid if isinstance(cid, str) else int(cid)), "n": len(mem),
        "shape": sh, "peak_hour": round(ph, 2),
        "mean_r": (round(rbar, 3) if rbar == rbar else None),
    } for cid, mem, rbar, sh, ph in rows])
    used: set = set()
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
        summary.to_excel(xw, sheet_name="summary", index=False)
        for cid, mem, rbar, sh, ph in rows:
            recs = []
            for m in mem:
                d = {"member": (item_label(m) if item_label else str(m))}
                if meta:
                    d.update(meta.get(m, {}))
                recs.append(d)
            df = pd.DataFrame(recs)
            if member_cols:
                cols = ["member"] + [c for c in member_cols if c in df.columns]
                df = df[[c for c in cols if c in df.columns]]
            df.to_excel(xw, sheet_name=_sheet_name(cid, used), index=False)
    return out_xlsx


def render_flat_panel(cols, traces_raw, grid, out, item_label, *,
                      label="flat / non-varying", ylim=None, list_max=72, title=""):
    """One A4-portrait page: a SINGLE overview panel of every 'flat / non-varying'
    trace overlaid (raw cps, log y) + their median. These were pulled OUT of
    correlation clustering (flat traces have no reliable shape, so they only
    bloat the cluster count). If few enough (<= list_max) the members are listed
    below; otherwise just the count (the full membership is in the set's CSV).
    Returns the path, or None for an empty set."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if cols is None or not len(cols):
        return None
    listing = len(cols) <= list_max
    PAGE_W, PAGE_H = 8.27, 11.69
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    fig.text(0.075, 0.955, f"{title}", fontsize=12, weight="bold", color="#222")
    fig.text(0.075, 0.935, f"n={len(cols)} · cv < {FLAT_CV:g} — pulled from correlation "
             "clustering (flat traces have no reliable shape, so they only bloat the cluster count)",
             fontsize=8.5, color="#666")
    # taller panel when there's no member list to show below it
    ax = (fig.add_axes([0.10, 0.50, 0.86, 0.36]) if listing
          else fig.add_axes([0.10, 0.34, 0.86, 0.52]))
    import warnings
    M = []
    for c in cols:
        y = traces_raw[c].values.astype(float)
        yy = np.where(y > 0, y, np.nan)
        ax.plot(grid, yy, color="#888780", lw=0.6, alpha=0.30)
        M.append(yy)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN time bins are fine
        med = np.nanmedian(np.array(M), axis=0)
    ax.plot(grid, med, color="#111", lw=2.2, alpha=0.9, zorder=5, label="median")
    ax.set_yscale("log"); ax.set_xlim(0, float(grid[-1]))
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_ylabel("cps", fontsize=9); ax.set_xlabel("hour of experiment (UTC)", fontsize=9)
    ax.grid(alpha=0.18, which="both"); ax.tick_params(labelsize=9)
    ax.legend(fontsize=9, loc="upper right")
    if listing:
        w = _wrap([item_label(c) for c in cols], width=96)
        tx = fig.add_axes([0.075, 0.06, 0.89, 0.40]); tx.axis("off")
        tx.text(0, 1, "\n".join(w), va="top", ha="left", fontsize=8.8,
                family="monospace", color="#333")
    else:
        fig.text(0.075, 0.30, f"{len(cols)} peaks — full membership in the set CSV.",
                 fontsize=9, color="#666")
    fig.savefig(out, dpi=200); plt.close(fig)
    return out


def render_clusters(rows, grid, traces_z, traces_raw, item_label, out_prefix, *,
                    remaining=None, per_chunk=36, title="", **kw):
    """A4-portrait cluster pages PLUS the remaining peaks (singletons / <3-member)
    split into LABELLED overview panels of <=per_chunk each, so every peak is both
    plotted and listed. `remaining` is the member list (or None). Returns paths."""
    all_rows = list(rows)
    if remaining:
        rem = list(remaining)
        for i in range(0, len(rem), per_chunk):
            chunk = rem[i:i + per_chunk]
            Lg = np.log10(traces_raw[chunk].clip(lower=1e-9))
            Z = (Lg - Lg.mean()) / Lg.std()
            mz = smooth(Z.mean(axis=1).values)
            ph = float(grid[int(np.nanargmax(mz))])
            lbl = f"remaining {i + 1}-{i + len(chunk)} of {len(rem)} (singletons / <3-member)"
            all_rows.append((lbl, chunk, float("nan"), shape_of(mz), ph))
    return render_a4(all_rows, grid, traces_z, traces_raw, item_label, out_prefix,
                     title=title, **kw)

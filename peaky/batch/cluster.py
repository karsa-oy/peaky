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

import textwrap

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

__version__ = "0.1.0"

CHANGING = 0.30          # cv at/above which a trace is "changing"
FLAT_CV = CHANGING       # cv BELOW which a trace counts as flat (the sustained-
                         # variation half of the gate; lower it to keep weak families)
SMOOTH_W = 3             # smoothing window for the transient-burst detector
PEAK_RANGE = 1.7         # smoothed max/median at/above which a trace has a coherent
                         # transient burst -> clustered even when its cv < FLAT_CV
                         # (a brief spike barely moves cv; this catches it)
FLAT_CLUSTER_RANGE = 1.4  # a cluster whose MEMBER-MEAN smoothed max/median is below this
                          # is a flat family (members co-vary but the family doesn't move)
                          # -> demoted to the flat-background overview
BIG_CHANGE_FOLD = 3.0    # a single channel whose smoothed max/median is >= this is a "big
                         # standalone change" (~>=5-10x raw) — surfaced even with no family
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


def panel_median(M):
    """Median of raw-cps member traces for a cluster panel. The per-trace NaN holes
    are 'below-detection' LOWS (the peak dropped under the bin floor, e.g. during a
    zero-air event), NOT missing-at-random -- so fill them to the cluster's detection
    floor (lowest detected value) and take a PLAIN median. np.nanmedian would DROP the
    holes and median only the surviving bright traces -> survivorship bias -> the bold
    line wrongly stays flat/rises through a zeroing event instead of dipping with it."""
    Mr = np.asarray(M, dtype=float)
    pos = Mr[np.isfinite(Mr) & (Mr > 0)]
    floor = float(np.min(pos)) if pos.size else 1.0
    return np.median(np.where(np.isfinite(Mr), Mr, floor), axis=0)


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
    points. 0.0 if the trace is empty/flat-at-zero. Catches SUSTAINED variation
    but is blind to a brief spike (1-2 bins barely move std/mean)."""
    if col not in traces:
        return 0.0
    tr = pd.to_numeric(traces[col], errors="coerce").to_numpy()
    tr = tr[np.isfinite(tr) & (tr > 0)]
    m = tr.mean() if len(tr) else 0.0
    return float(tr.std() / m) if m > 0 else 0.0


def _smoothed_range(y, smooth_w=SMOOTH_W) -> float:
    """Smoothed max / median of a 1-D trace — how much it MOVES. ~1 = flat; a real
    rise/fall/burst scores higher. Smoothing first rejects single-bin noise spikes."""
    ys = smooth(np.asarray(y, float), smooth_w)
    pos = ys[np.isfinite(ys) & (ys > 0)]
    if not len(pos):
        return 0.0
    med = np.median(pos)
    return float(np.max(pos) / med) if med > 0 else 0.0


def trace_dynamic_range(traces: pd.DataFrame, col, *, smooth_w=SMOOTH_W) -> float:
    """Smoothed max / median of one trace — sensitive to a coherent TRANSIENT burst
    the global CV misses (a synchronized spike scores high; jitter stays ~1)."""
    if col not in traces:
        return 0.0
    return _smoothed_range(pd.to_numeric(traces[col], errors="coerce").to_numpy(), smooth_w)


def cluster_flatness(members, traces: pd.DataFrame, *, smooth_w=SMOOTH_W) -> float:
    """Smoothed max/median of a cluster's MEMBER-MEAN raw trace — how much the FAMILY
    as a whole moves. ~1 means the members co-vary but the family is flat (correlated
    background: it passed the correlation cut yet has no real dynamics)."""
    import warnings
    cols = [m for m in members if m in traces.columns]
    if not cols:
        return 0.0
    M = np.vstack([np.where(traces[m].to_numpy(float) > 0, traces[m].to_numpy(float), np.nan)
                   for m in cols])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN time bins are fine
        mean = np.nanmean(M, axis=0)
    return _smoothed_range(mean, smooth_w)


def split_flat_clusters(rows, traces: pd.DataFrame, *, range_min=FLAT_CLUSTER_RANGE,
                        smooth_w=SMOOTH_W):
    """Partition cluster `rows` (cluster_rows output) into (dynamic, flat). A cluster
    is FLAT if its member-mean trace has no real excursion (cluster_flatness <
    range_min) — its members correlate but the family itself doesn't move, so it
    belongs in the flat-background overview, not shown as a 'co-varying family'."""
    dyn, flat = [], []
    for row in rows:
        (dyn if cluster_flatness(row[1], traces, smooth_w=smooth_w) >= range_min
         else flat).append(row)
    return dyn, flat


def trace_varies(traces: pd.DataFrame, col, *, cv_min=FLAT_CV, range_min=PEAK_RANGE,
                 smooth_w=SMOOTH_W) -> bool:
    """A trace gets correlation-clustered (vs bunched as flat) if it has SUSTAINED
    variation (cv >= cv_min) OR a coherent TRANSIENT burst (smoothed max/median >=
    range_min). CV alone left real co-varying burst families in the flat bucket — a
    brief synchronized spike barely changes std/mean. The clustering's r>0.6 /
    >=3-member rule then filters any noise this lets back in."""
    return (trace_cv(traces, col) >= cv_min
            or trace_dynamic_range(traces, col, smooth_w=smooth_w) >= range_min)


def split_varying(traces: pd.DataFrame, cols, *, cv_min=FLAT_CV, range_min=PEAK_RANGE,
                  smooth_w=SMOOTH_W):
    """Partition `cols` into (varying, flat). A flat trace has no reliable SHAPE
    (neither sustained variation nor a transient burst) — its pairwise correlation
    is noise, so flat traces spuriously shatter into tiny clusters. Pull them out
    BEFORE clustering and bunch them. Returns (varying, flat), order-preserving."""
    varying, flat = [], []
    for c in dict.fromkeys(cols):
        (varying if trace_varies(traces, c, cv_min=cv_min, range_min=range_min,
                                 smooth_w=smooth_w) else flat).append(c)
    return varying, flat


MERGE_R = 0.85           # merge clusters whose MEAN traces correlate at/above this
MERGE_LINK = "complete"  # COMPLETE linkage on the centroids — average/single chain
                         # distinct shapes into one blob (the very thing the primary
                         # clustering uses complete linkage to avoid)


def merge_similar(traces: pd.DataFrame, lab: pd.Series, big, *, merge_r=MERGE_R,
                  link=MERGE_LINK, min_members=MIN_MEMBERS):
    """Second-level merge: collapse clusters whose MEAN traces are near-identical
    (centroid correlation >= merge_r) into one. Clustering a decay-dominated batch
    on shape over-splits the dominant shape into many near-duplicate clusters (and
    leaves similar small clusters); this folds those together while keeping
    genuinely distinct families (e.g. off-decay bursts) apart.

    `traces` = the (log / normalised) per-item trace frame used for correlation
    (index = time bins, columns = items). Returns (new_lab, new_big), merged-cluster
    ids in a namespace disjoint from the un-merged singleton ids."""
    big = list(big)
    if len(big) < 2:
        return lab, big
    cents = {c: traces[[x for x in lab.index if int(lab[x]) == c]].mean(axis=1) for c in big}
    cc = pd.DataFrame(cents).corr()
    dist = (1 - cc.fillna(0)).values.copy()
    np.fill_diagonal(dist, 0.0); dist = (dist + dist.T) / 2
    Z = linkage(squareform(dist, checks=False), method=link)
    grp = pd.Series(fcluster(Z, t=1.0 - merge_r, criterion="distance"), index=list(cc.columns))
    offset = int(pd.to_numeric(lab, errors="coerce").max()) + 1   # disjoint id space
    new_lab = lab.copy()
    for x in lab.index:
        c = int(lab[x])
        if c in grp.index:
            new_lab[x] = offset + int(grp[c])
    sizes = new_lab.value_counts()
    new_big = sorted((int(g) for g in sizes.index if sizes[g] >= min_members),
                     key=lambda g: -int(sizes[g]))
    return new_lab, new_big


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
                    else panel_median(M)
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
                           member_cols=None, title="clusters", when=None):
    """Per-cluster Excel workbook: a 'summary' sheet (one row per cluster: id, n,
    shape, peak hour, mean r) + ONE sheet PER CLUSTER listing its members. `rows`
    is the cluster_rows output [(cid, members, rbar, shape, peak_hr), ...]; `meta`
    maps a member key -> a dict of columns (e.g. neutral_formula / channel / m_z /
    match_score / tier). `member_cols` orders/selects those columns. `when` is the
    timestamp embedded in the file; it defaults (via _resolve_when) to the FIXED
    content epoch, so the workbook is byte-identical for identical data, whenever
    it's written. Returns the path (or None if no rows)."""
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
    _make_xlsx_deterministic(out_xlsx, when=when)
    return out_xlsx


def _resolve_when(when):
    """The timestamp to embed in a workbook. Precedence: explicit `when` (a datetime
    or epoch seconds) -> SOURCE_DATE_EPOCH env (which the run driver pins to the FIXED
    content epoch, not the run time) -> the same fixed content epoch as a last resort.
    Because the embedded stamp is a constant, the workbook is a PURE FUNCTION OF ITS
    DATA: the same data yields byte-identical bytes whenever it's run. (Run time lives
    only on the PDF cover + Report ID, never in the data files.)"""
    import os
    from datetime import datetime, timezone
    # mirror pipeline.CONTENT_EPOCH (1980-01-01Z); duplicated as a bare int to avoid a
    # cluster -> pipeline import cycle (pipeline imports cluster transitively).
    _CONTENT_EPOCH = 315532800
    if when is None:
        sde = os.environ.get("SOURCE_DATE_EPOCH")
        when = int(sde) if sde else None
    if when is None:
        return datetime.fromtimestamp(_CONTENT_EPOCH, tz=timezone.utc)
    if isinstance(when, (int, float)):
        return datetime.fromtimestamp(int(when), tz=timezone.utc)
    return when


def _make_xlsx_deterministic(path, *, when=None) -> None:
    """Rewrite an .xlsx so its bytes are a deterministic function of ITS DATA:
    replace openpyxl's datetime.now() stamps (docProps/core.xml AND every zip
    member's embedded mtime) with `when` (see _resolve_when), which defaults to the
    FIXED content epoch. Same data -> identical bytes, whenever it's written. Zip
    mtimes must be >= 1980 (the format floor); the content epoch is 1980-01-01."""
    import os
    import re
    import zipfile
    w = _resolve_when(when)
    dt = (max(w.year, 1980), w.month, w.day, w.hour, w.minute, w.second)
    iso = w.strftime("%Y-%m-%dT%H:%M:%SZ").encode()
    with zipfile.ZipFile(path) as zin:
        items = [(n, zin.read(n)) for n in sorted(zin.namelist())]
    tmp = f"{path}.tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for name, data in items:
            if name == "docProps/core.xml":
                data = re.sub(rb"(<dcterms:(?:created|modified)[^>]*>)[^<]*(</dcterms:)",
                              lambda m: m.group(1) + iso + m.group(2), data)
            zi = zipfile.ZipInfo(name, date_time=dt)
            zi.compress_type = zipfile.ZIP_DEFLATED
            zi.external_attr = 0o600 << 16
            zout.writestr(zi, data)
    os.replace(tmp, path)


def big_changers(traces: pd.DataFrame, cols, grid, *, fold_min=BIG_CHANGE_FOLD, smooth_w=2):
    """Channels that change A LOT on their own — peak / baseline (smoothed max over
    the 10th-percentile) >= fold_min — regardless of whether they co-vary with
    anything. Using the low-percentile baseline (not the median) means a monotonic
    decay/rise counts as much as a spike: any 'huge change then long tail' is caught,
    not just a peak above a flat median. These have no family so they'd otherwise sit
    unnoticed in the flat panel. Returns [(col, fold, peak_hour), ...], largest first."""
    out = []
    for c in cols:
        if c not in traces.columns:
            continue
        ys = smooth(pd.to_numeric(traces[c], errors="coerce").to_numpy(), smooth_w)
        pos = ys[np.isfinite(ys) & (ys > 0)]
        if len(pos) < MIN_POINTS:
            continue
        base = np.percentile(pos, 10)         # robust low baseline (not the median)
        fold = float(np.max(pos) / base) if base > 0 else 0.0
        if fold >= fold_min:
            ph = float(grid[int(np.nanargmax(np.where(np.isfinite(ys), ys, -np.inf)))])
            out.append((c, fold, ph))
    out.sort(key=lambda t: -t[1])
    return out


def render_changers(items, traces_raw, grid, out_prefix, item_label, *, cap=48,
                    title="", dpi=200):
    """A4-PORTRAIT small-multiples of the big standalone changers, so this section
    keeps the report's A4 page format (it used to embed a native-sized strip that
    broke it). One mini-plot per channel (raw cps, log y) titled `formula+adduct
    Nx · peak hour`; panels are packed top-down at a FIXED legible size (1 column for
    <=2 channels, else 2) and paginated when an A4 page fills. The log y-axis snaps to
    whole decades with only MAJOR ticks labelled (a near-flat trace otherwise crams
    colliding 2..9x10^n minor labels). Writes <out_prefix>_p<i>.png; returns the
    list of page paths (matching render_a4)."""
    import math

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    items = list(items)[:cap]
    if not items:
        return []
    W, H = 8.27, 11.69                                    # A4 portrait (inches)
    L, RM, TOPM, BOTM, GUT = 0.62, 0.40, 0.80, 0.45, 0.55
    ncol = 1 if len(items) <= 2 else 2
    panel_w = (W - L - RM - GUT * (ncol - 1)) / ncol
    PANEL_H = 2.15 if ncol == 1 else 1.70                 # trace height (in)
    TITLE_H, XLAB_H, VGAP = 0.30, 0.34, 0.24
    pitch = TITLE_H + PANEL_H + XLAB_H + VGAP
    rows_pp = max(1, int((H - TOPM - BOTM) // pitch))
    per_page = rows_pp * ncol
    pages = [items[i:i + per_page] for i in range(0, len(items), per_page)]
    paths = []
    for pi, page in enumerate(pages, 1):
        fig = plt.figure(figsize=(W, H))
        fig.text(L / W, (H - 0.45) / H, title, fontsize=11, weight="bold", color="#222")
        if len(pages) > 1:
            fig.text(1 - RM / W, (H - 0.45) / H, f"page {pi}/{len(pages)}",
                     fontsize=10, color="#777", ha="right")
        for k, (c, fold, ph) in enumerate(page):
            r, cc = divmod(k, ncol)
            ax_l = (L + cc * (panel_w + GUT)) / W
            ax_b = (H - TOPM - r * pitch - TITLE_H - PANEL_H) / H
            ax = fig.add_axes([ax_l, ax_b, panel_w / W, PANEL_H / H])
            yp = np.where(traces_raw[c].to_numpy(float) > 0, traces_raw[c].to_numpy(float), np.nan)
            ax.plot(grid, yp, color="#1D9E75", lw=1.2, marker="o", ms=2.5)
            ax.set_yscale("log"); ax.set_xlim(0, float(grid[-1]))
            finite = yp[np.isfinite(yp)]
            if finite.size:                              # snap y-limits to whole decades
                lo, hi = float(np.nanmin(finite)), float(np.nanmax(finite))
                c0 = 10.0 ** math.floor(math.log10(lo)) if lo > 0 else 1.0
                c1 = 10.0 ** math.ceil(math.log10(hi)) if hi > 0 else c0 * 10
                ax.set_ylim(c0, max(c1, c0 * 10))
            ax.yaxis.set_major_locator(mticker.LogLocator(base=10, numticks=6))
            ax.yaxis.set_minor_formatter(mticker.NullFormatter())   # kill 2..9x labels
            ax.grid(alpha=0.2, which="both"); ax.tick_params(labelsize=7.5)
            ax.set_xlabel("hour (UTC)", fontsize=7.5); ax.set_ylabel("cps", fontsize=7.5)
            ax.set_title(f"{item_label(c)}   {fold:.0f}× · h{ph:.1f}", fontsize=8.5, loc="left")
        out = f"{out_prefix}_p{pi}.png"
        fig.savefig(out, dpi=dpi); plt.close(fig)
        paths.append(out)
    return paths


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
    # wrap the caption to the page width (va='top' so it grows DOWN, toward the
    # plot, never up into the bold title) — a long one-liner used to run off the
    # right edge and get clipped.
    _sub = (f"n={len(cols)} — channels with no real family dynamics: uncorrelated, or in a "
            "co-varying cluster whose mean is flat, + Si contamination; bunched so they "
            "don't bloat the cluster count")
    fig.text(0.075, 0.935, textwrap.fill(_sub, width=112), fontsize=8.5, color="#666",
             va="top")
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
        med = panel_median(M)
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

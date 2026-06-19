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


def cluster_rows(cols, lab, big, cm, traces_raw, grid, *, median_cps=None, maxc=10):
    """Build per-cluster summary rows (sorted by size) for rendering:
    (cid, members[brightest-first], r_bar, shape, peak_hour)."""
    Lg = np.log10(traces_raw[cols].clip(lower=1e-9))
    Z = (Lg - Lg.mean()) / Lg.std()
    if median_cps is None:
        median_cps = {c: float(np.nanmedian(traces_raw[c].values)) for c in cols}
    rows = []
    for cid in big[:maxc]:
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


def render_panels(rows, grid, traces_z, traces_raw, item_label, out, *,
                  mode="raw", title="", event_span=None, ylim=None):
    """Stacked per-cluster panels. mode='z' (z-scored shape, mean) or 'raw'
    (log cps, median). `item_label` maps a column -> its printed label
    (formula or m/z). `event_span` shades an (h0,h1) window."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not rows:
        return None
    Wf, Lm, Rm, Tm, Bm = 11.0, 0.95, 0.5, 0.5, 0.25
    plot_h, lh, gap = 1.7, 0.158, 0.5
    wrapped = [_wrap([item_label(m) for m in mem]) for _, mem, _, _, _ in rows]
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
        ax.set_title(f"cluster {cid} · n={len(mem)} · r̄={rbar:.2f} · {sh} (peak~h{ph:.1f})",
                     fontsize=9, loc="left")
        ax.set_xlabel("hour of experiment (UTC)") if k == len(rows) - 1 else ax.set_xticklabels([])
        tx = fig.add_axes([lf, (cur - h0 - th) / Hf, wf, th / Hf]); tx.axis("off")
        tx.text(0, 1, "\n".join(w), va="top", ha="left", fontsize=7.6,
                family="monospace", color="#333", transform=tx.transAxes)
        cur -= (h0 + th + gap)
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    return out

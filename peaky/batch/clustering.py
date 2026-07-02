"""Correlation-cluster figures for a batch — the in-package home of what used to
be the out-of-repo scratch driver `run_clusters.py`.

`cluster_batch(out_dir, ts, profile)` produces, from a representative-merge
`merged_ledger.csv` (+ the per-file ledgers it left) and the full-batch peak time
series, three figure sets + their CSVs:

  * CHANGING  — assigned analytes clustered PER ION CHANNEL on RAW log-correlation,
                near-identical shapes merged, flat-mean clusters demoted, with the
                big standalone changers pulled into their own small-multiples page;
  * FLAT      — the uncorrelated remainder + Si contamination, bunched into one panel;
  * UNASSIGNED — TS bins not matched to any EXPLAINED peak, gated on cv then clustered;

plus a per-cluster workbook and a channel-agreement QC table. Everything is
parameterised by the ReagentProfile (adducts / normaliser / reagent-ion regex /
label) and the output dir — nothing reagent- or batch-specific is hardcoded.

Pure local compute (no network); reads/writes only under `out_dir`.
"""
from __future__ import annotations

import glob
import json
import os

import numpy as np
import pandas as pd

from peaky.reporting import analyte_viz as V
from peaky.batch import cluster as CL
from peaky import paths as PT
from peaky.batch import timeseries as TS

__version__ = "0.1.0"

FLOOR_DEFAULT = 200.0   # min median cps for a channel to enter clustering (legacy gate)


def _longest_detected_run(tr) -> int:
    """Longest run of consecutive DETECTED (nonzero) bins in a raw trace.
    NaN / 0 = not detected. Distinguishes a real multi-bin episode (rise+decay)
    from a sporadic single-bin spike at the detection limit. Intensity-agnostic:
    unlike a median-cps floor it does not penalise sharp transients whose signal
    is zero most of the run (e.g. prompt accretion dimers)."""
    y = np.nan_to_num(np.asarray(tr, dtype=float), nan=0.0)
    best = cur = 0
    for v in y:
        cur = cur + 1 if v > 0 else 0
        if cur > best:
            best = cur
    return best


def cluster_batch(out_dir, ts, profile, *, merged=None, tag=None, label=None,
                  floor: float = FLOOR_DEFAULT, bin_minutes: int | None = None,
                  gate: str = "median", min_run: int = 3,
                  log=print) -> dict:
    """Build the cluster/flat/unassigned figure sets for one batch.

    out_dir : run folder; holds `merged_ledger.csv` + `per_file/*_ledger.csv`
              (written by assign_batch) and receives the figures + CSVs.
    ts      : the full-batch per-sample peak time series (DataFrame or parquet path).
    profile : ReagentProfile (adducts / normaliser / reagent_ion_re / label).
    merged  : merged ledger (DataFrame or path); default `out_dir/merged_ledger.csv`.
    tag     : filename token (default profile.name); label: title text (profile.label).
    """
    OUT = os.path.expanduser(out_dir)
    P = PT.run_paths(OUT).ensure()
    FIG, TAB = P.figures, P.tables       # .png -> figures/, .csv/.xlsx -> tables/
    tag = tag or profile.name
    label = label or profile.label
    ADDUCTS, NORM = profile.adducts, profile.normaliser
    FLOOR = floor

    # clear stale cluster pages so a shorter run (fewer pages) can't leave orphans
    # that the report would still glob in
    for _old in glob.glob(f"{FIG}/clusters_*_{tag}_p*.png"):
        os.remove(_old)

    if merged is None:
        merged = os.path.join(OUT, "merged_ledger.csv")
    merged = pd.read_csv(merged) if isinstance(merged, str) else merged.copy()
    merged["role"] = "M0"                                  # analyte_table expects a role col

    if isinstance(ts, str):
        ts = pd.read_parquet(os.path.expanduser(ts))
    ts = ts.copy()
    ts["datetime_utc"] = pd.to_datetime(ts["datetime_utc"], utc=True)
    # Native per-sample time resolution by default (BIN_MIN=None): the samples are
    # the common time axis, so traces + correlation use one point per sample at its
    # real time — no re-gridding onto a uniform lattice (which aliased into spurious
    # empty bins). Pass bin_minutes=int to time-bin instead.
    BIN_MIN = bin_minutes
    span_min = (ts["datetime_utc"].max() - ts["datetime_utc"].min()).total_seconds() / 60.0
    log(f"=== {tag} ({label}) : merged M0={len(merged)}, TS samples={ts['sample_item_id'].nunique()}, "
        f"span={span_min:.0f}min -> res={'native/sample' if BIN_MIN is None else str(BIN_MIN)+'min'} ===")

    def cv_of(traces, cols):
        out = {}
        for f in cols:
            tr = traces[f].dropna().values if f in traces else np.array([])
            out[f] = float(np.std(tr) / np.mean(tr)) if len(tr) and np.mean(tr) > 0 else np.nan
        return out

    def reagent_mzs():
        # the merged ledger is M0-only (no ion_formula); read the reagent ions from a
        # per-file ledger (keeps role/ion_formula) so Br gets its physical /Br3- norm.
        if not profile.reagent_ion_re:
            return None
        for f in sorted(glob.glob(f"{OUT}/per_file/*_ledger.csv")):
            d = pd.read_csv(f)
            if "ion_formula" not in d.columns:
                continue
            m = d[d["ion_formula"].astype(str).str.match(profile.reagent_ion_re)]["mz"].dropna()
            if len(m):
                return sorted(set(round(float(x), 4) for x in m))
        return None

    # ---- assigned analytes, clustered PER ION CHANNEL (formula+adduct) ----------
    # Cluster each ion channel SEPARATELY (not the per-neutral SUM): channels of one
    # neutral often diverge in time, so summing blends divergent shapes.
    chan = merged[merged["role"] == "M0"].dropna(subset=["neutral_formula"]).copy()
    chan = chan[chan["neutral_formula"].astype(str) != ""]
    chan["key"] = chan["neutral_formula"].astype(str) + "|" + chan["adduct"].astype(str)
    chan = chan.drop_duplicates("key")
    ion_mz = dict(zip(chan["key"], chan["mz"]))
    is_si = dict(zip(chan["key"], chan["neutral_formula"].map(V._is_contaminant)))

    grid, traces_raw = V.ion_traces(ts, ion_mz, mode="raw", bin_minutes=BIN_MIN)
    med = {k: (float(np.nanmedian(traces_raw[k])) if k in traces_raw and traces_raw[k].notna().any() else 0.0)
           for k in ion_mz}
    cv = cv_of(traces_raw, list(ion_mz))

    def _enter(k):
        """Entry gate to clustering.
        gate='median'  : legacy median-cps floor (>=FLOOR) — blind to transients.
        gate='episode' : detected (nonzero) in >=min_run consecutive bins — rescues
                         sharp low-abundance episodes (e.g. accretion dimers) the
                         median floor drops, while still rejecting sporadic spikes."""
        if k not in traces_raw:
            return False
        if gate == "episode":
            return _longest_detected_run(traces_raw[k].values) >= min_run
        return med.get(k, 0) >= FLOOR

    # per-ion legend (formula+adduct + match-score) + per-cluster workbook metadata
    keylab = {r.key: f"{V.ion_label(r.neutral_formula, r.adduct)} ({float(r.ion_score):.2f})"
              for r in chan.itertuples()}
    ion_lbl = lambda k: keylab.get(k, k)
    meta = {r.key: {"neutral_formula": r.neutral_formula,
                    "channel": V.ADDUCT_SUFFIX.get(str(r.adduct), str(r.adduct)),
                    "adduct": r.adduct, "m_z": round(float(r.mz), 4),
                    "match_score": round(float(r.ion_score), 3), "tier": r.tier,
                    "median_cps": round(med.get(r.key, 0), 0),
                    "cv": round(cv.get(r.key, float("nan")), 3)}
            for r in chan.itertuples()}
    MCOLS = ["neutral_formula", "channel", "m_z", "match_score", "tier", "median_cps", "cv"]

    rmz = reagent_mzs()
    norm_mode = NORM if (NORM != "reagent" or rmz) else "tic"

    def goodk(keys):
        return [k for k in dict.fromkeys(keys)
                if k in traces_raw and np.isfinite(traces_raw[k].values.astype(float)).sum() >= CL.MIN_POINTS]

    # Cluster ALL bright organic ion channels on SHAPE (RAW corr) — NO per-trace cv
    # gate (it can't see coherence). Then MERGE near-identical-shape clusters. The
    # non-clustering remainder + Si contamination is the genuinely-flat bucket.
    clust_cols = goodk([k for k in ion_mz if _enter(k) and not is_si[k]])
    log(f"CLUSTERING {len(clust_cols)} organic ion-channels on RAW shape (entry gate='{gate}'"
        f"{f', min_run={min_run}' if gate=='episode' else f', floor={FLOOR:.0f}cps'}; no cv gate)")
    Lg, cm = CL.correlate(traces_raw, clust_cols)
    lab, big = CL.cluster(cm) if len(clust_cols) >= CL.MIN_MEMBERS else (pd.Series(dtype=int), [])
    log(f"  raw: {len(big)} clusters>=3")
    lab, big = CL.merge_similar(Lg, lab, big, merge_r=CL.MERGE_R)   # collapse near-duplicate shapes
    rows, _ = CL.cluster_rows(clust_cols, lab, big, cm, traces_raw, grid, median_cps=med)
    # DEMOTE flat clusters: members co-vary but the family MEAN doesn't move.
    rows, flat_rows = CL.split_flat_clusters(rows, traces_raw)
    dyn_ids = {r[0] for r in rows}
    flat_cluster_members = [m for r in flat_rows for m in r[1]]
    remainder = [k for k in clust_cols if int(lab.get(k, -1)) not in dyn_ids]   # not in a DYNAMIC family
    # BIG STANDALONE CHANGERS: single channels that change a lot on their own.
    changers = CL.big_changers(traces_raw, remainder, grid, fold_min=CL.BIG_CHANGE_FOLD)
    changer_set = {c for c, _, _ in changers}
    remainder = [k for k in remainder if k not in changer_set]
    CL.render_changers(changers, traces_raw, grid, f"{FIG}/clusters_changers_{tag}", ion_lbl,
                       title=f"{label} · Large standalone changes (≥{CL.BIG_CHANGE_FOLD:g}× fold, no family) — {len(changers)} channels")
    pd.DataFrame({"ion": [ion_lbl(c) for c, _, _ in changers],
                  "neutral_formula": [c.split("|")[0] for c, _, _ in changers],
                  "channel": [meta[c]["channel"] for c, _, _ in changers],
                  "fold": [round(f, 1) for _, f, _ in changers],
                  "peak_hour": [round(ph, 2) for _, _, ph in changers],
                  "median_cps": [round(med[c], 0) for c, _, _ in changers]}).to_csv(f"{TAB}/clusters_changers_{tag}.csv", index=False)
    log(f"CHANGING: {len(rows)} dynamic families covering {sum(len(r[1]) for r in rows)}; "
        f"{len(flat_rows)} flat clusters ({len(flat_cluster_members)} ch) demoted; "
        f"{len(changers)} big standalone changers; "
        f"{len(remainder)} flat")
    for cid, mem, rbar, sh, ph in rows:
        log(f"  c{cid}: n={len(mem)} r̄={rbar:.2f} {sh} h{ph:.1f} | {', '.join(ion_lbl(m) for m in mem[:5])}")
    posc = traces_raw[clust_cols].values if clust_cols else np.array([])
    posc = posc[np.isfinite(posc) & (posc > 0)]
    # top = true max (+20% log headroom), NOT a 99.5 pct cap, so the brightest
    # traces are never clipped at the high end; bottom stays a 1-pct/50-cps floor.
    ylimc = (max(50, np.percentile(posc, 1)), float(np.nanmax(posc)) * 1.2) if len(posc) else None
    Zc = (Lg - Lg.mean()) / Lg.std() if len(clust_cols) else Lg
    CL.render_clusters(rows, grid, Zc, traces_raw, ion_lbl, f"{FIG}/clusters_changing_{tag}",
                       remaining=None, mode="raw", ylim=ylimc,
                       title=f"{label} · Co-varying ion channels — {sum(len(r[1]) for r in rows)} channels in {len(rows)} families",
                       subtitle="legend:  formula+adduct  (match-score) — shape clusters; near-identical merged, flat families demoted")
    # per-cluster workbook: one tab per (dynamic) cluster
    CL.write_cluster_workbook(list(rows), f"{TAB}/clusters_changing_{tag}.xlsx",
                              meta=meta, item_label=ion_lbl, member_cols=MCOLS)
    clustered = [c for c in clust_cols if int(lab.get(c, -1)) in dyn_ids]   # CSV = dynamic families only
    pd.DataFrame({"ion": [ion_lbl(c) for c in clustered],
                  "neutral_formula": [c.split("|")[0] for c in clustered],
                  "channel": [meta[c]["channel"] for c in clustered],
                  "match_score": [meta[c]["match_score"] for c in clustered],
                  "tier": [meta[c]["tier"] for c in clustered],
                  "cluster": [int(lab[c]) for c in clustered],
                  "cv": [round(cv[c], 3) for c in clustered],
                  "median_cps": [round(med[c], 0) for c in clustered]}).to_csv(f"{TAB}/clusters_changing_{tag}.csv", index=False)

    # FLAT = the uncorrelated remainder + Si contamination -> bunched overview panel
    flat_cols = goodk(list(dict.fromkeys(remainder + [k for k in ion_mz if is_si[k] and _enter(k)])))
    pos = traces_raw[flat_cols].values if flat_cols else np.array([])
    pos = pos[np.isfinite(pos) & (pos > 0)]
    ylim = (max(50, np.percentile(pos, 1)), float(np.nanmax(pos)) * 1.2) if len(pos) else None
    log(f"FLAT: {len(flat_cols)} uncorrelated/contaminant ion channels (bunched)")
    CL.render_flat_panel(flat_cols, traces_raw, grid, f"{FIG}/clusters_flat_{tag}_p1.png", ion_lbl,
                         label=f"{label} · flat background (uncorrelated + flat families + Si)", ylim=ylim,
                         title=f"{label} · Flat background — {len(flat_cols)} ion channels (bunched)")
    pd.DataFrame({"ion": [ion_lbl(c) for c in flat_cols],
                  "neutral_formula": [c.split("|")[0] for c in flat_cols],
                  "channel": [meta[c]["channel"] for c in flat_cols], "cluster": 0,
                  "cv": [round(cv[c], 3) for c in flat_cols],
                  "median_cps": [round(med[c], 0) for c in flat_cols]}).to_csv(f"{TAB}/clusters_flat_{tag}.csv", index=False)

    # CHANNEL-AGREEMENT QC: do a neutral's ion channels actually track in time?
    ca = V.channel_agreement(ts, merged[["neutral_formula", "adduct", "mz"]], bin_minutes=BIN_MIN)
    ca.to_csv(f"{TAB}/channel_agreement_{tag}.csv", index=False)
    if len(ca):
        log(f"CHANNEL AGREEMENT: {len(ca)} multi-channel neutrals — "
            f"{ca['verdict'].value_counts().to_dict()}")

    # ---- UNASSIGNED: TS bins not matched to any EXPLAINED peak ------------------
    # match against ALL explained peaks (M0 + iso_child + reagent + artifact) from
    # the per-file ledgers, NOT just merged M0 (which would flag every satellite).
    _expl = []
    for _pf in glob.glob(f"{OUT}/per_file/*_ledger.csv"):
        _d = pd.read_csv(_pf)
        _expl += _d.loc[_d["role"].astype(str) != "unexplained", "mz"].dropna().tolist()
    mat, bin_mz = TS.build_matrix(ts)
    tstamp = ts.groupby("sample_item_id")["datetime_utc"].first()
    assigned = np.sort(np.array(_expl)) if _expl else np.sort(merged["mz"].dropna().to_numpy())

    def is_assigned(mz, tol=8.0):
        i = np.searchsorted(assigned, mz)
        return any(0 <= j < len(assigned) and abs(assigned[j] - mz) / mz * 1e6 <= tol for j in (i - 1, i))

    median_h = mat.median()
    un_bins = [b for b in mat.columns if not is_assigned(float(bin_mz[b]))
               and median_h[b] >= 50.0 and mat[b].notna().sum() >= CL.MIN_POINTS]
    log(f"UNASSIGNED bins: {len(un_bins)} (of {mat.shape[1]} TS bins)")

    def to_grid(m, bin_minutes=BIN_MIN):
        m = m.reindex(tstamp.sort_values().index)
        hr = (tstamp.reindex(m.index) - tstamp.min()).dt.total_seconds().values / 3600.0
        if bin_minutes is None:                            # NATIVE — one row per sample
            return hr, m.reset_index(drop=True)
        g = np.arange(0, np.nanmax(hr) + bin_minutes / 60.0, bin_minutes / 60.0)
        return g, m.groupby(np.digitize(hr, g)).median().reindex(range(1, len(g) + 1))

    rt = mat.sum(axis=1) if NORM == "tic" else TS.reagent_total(mat, bin_mz, rmz or [])
    ggrid, un_raw = to_grid(mat[un_bins])
    _, un_norm = to_grid(mat[un_bins].div((rt if rt is not None else 1).replace(0, np.nan), axis=0)
                         if rt is not None else mat[un_bins])
    # GATE: the unassigned set has no inherent cv gate; split flat bins out and
    # cluster only the varying ones, bunch the flat tail.
    un_vary, un_flat = CL.split_varying(un_raw, un_bins)
    log(f"UNASSIGNED gate: {len(un_vary)} varying / {len(un_flat)} flat (bunched)")
    mzlab = lambda b: f"{bin_mz[b]:.4f}"
    posu = un_raw[un_bins].values; posu = posu[np.isfinite(posu) & (posu > 0)]
    ylimu = (max(50, np.percentile(posu, 1)), float(np.nanmax(posu)) * 1.2) if len(posu) else None
    Lgu, cmu = CL.correlate(un_norm, un_vary)
    labu, bigu = CL.cluster(cmu)
    rowsu, Zu = CL.cluster_rows(un_vary, labu, bigu, cmu, un_raw, ggrid,
                                median_cps={b: float(median_h[b]) for b in un_vary})
    log(f"UNASSIGNED: {len(bigu)} clusters>=3 covering {sum(len(r[1]) for r in rowsu)}")
    for cid, mem, rbar, sh, ph in rowsu:
        log(f"  c{cid}: n={len(mem)} r̄={rbar:.2f} {sh} h{ph:.1f} m/z {bin_mz[mem].min():.1f}-{bin_mz[mem].max():.1f}")
    rem_u = CL.remaining_row(un_vary, labu, bigu, un_raw, ggrid)[1] if len(un_vary) else None
    paths_u = CL.render_clusters(rowsu, ggrid, Zu, un_raw, mzlab, f"{FIG}/clusters_unassigned_{tag}",
                                 remaining=rem_u, mode="raw", ylim=ylimu,
                                 title=f"{label} · Unexplained peaks — {len(un_vary)} varying bins in {len(bigu)} clusters")
    if un_flat:
        CL.render_flat_panel(un_flat, un_raw, ggrid, f"{FIG}/clusters_unassigned_{tag}_p{len(paths_u)+1}.png",
                             mzlab, label=f"{label} · flat / non-varying unexplained", ylim=ylimu,
                             title=f"{label} · Unexplained — {len(un_flat)} flat / non-varying bins (bunched)")
    labu_idx = set(labu.index)
    pd.DataFrame({"mz": [round(float(bin_mz[b]), 4) for b in un_bins],
                  "cluster": [int(labu[b]) if b in labu_idx else 0 for b in un_bins],
                  "median_cps": [round(float(median_h[b]), 0) for b in un_bins]}).to_csv(f"{TAB}/clusters_unassigned_{tag}.csv", index=False)

    # GATE / FUNNEL summary -> the single source of truth the PDF report reads to
    # DOCUMENT the thresholds these figures apply and the unexplained funnel counts
    # (so a reader knows why only N of the unexplained peaks are drawn). Thresholds
    # come from the real constants in cluster.py / this module, not magic numbers.
    n_unassigned_any = int(sum(1 for b in mat.columns
                               if not is_assigned(float(bin_mz[b]))))
    summary = {
        "tag": tag,
        "gates": {
            "match_tol_ppm": 8.0,                       # is_assigned() tolerance
            "unassigned_median_cps_floor": 50.0,        # un_bins brightness floor
            "assigned_clustering_floor_cps": FLOOR,     # assigned-channel floor (median gate)
            "entry_gate": gate,                         # 'median' | 'episode'
            "min_consecutive_bins": min_run,            # episode-gate run length
            "min_trace_points": CL.MIN_POINTS,          # persistence gate
            "varying_cv_min": CL.FLAT_CV,               # sustained-change half
            "varying_burst_range": CL.PEAK_RANGE,       # transient-burst half
            "cluster_corr_r": round(1 - CL.DIST_T, 2),  # complete-linkage cut
            "merge_corr_r": CL.MERGE_R,                 # near-duplicate merge
            "min_cluster_members": CL.MIN_MEMBERS,
            "big_change_fold": CL.BIG_CHANGE_FOLD,
        },
        "unassigned": {
            "n_ts_bins": int(mat.shape[1]),
            "n_unassigned_any": n_unassigned_any,
            "n_after_brightness_persistence": len(un_bins),
            "n_varying_plotted": len(un_vary),
            "n_flat_bunched": len(un_flat),
            "n_clusters": len(bigu),
        },
        "assigned": {
            "n_dynamic_families": len(rows),
            "n_channels_in_families": int(sum(len(r[1]) for r in rows)),
            "n_flat_clusters_demoted": len(flat_rows),
            "n_big_changers": len(changers),
            "n_flat_background": len(flat_cols),
        },
    }
    with open(f"{OUT}/clusters_summary.json", "w") as _fh:
        json.dump(summary, _fh, indent=2)

    log(f"DONE — figures in {OUT}")
    return {"changing": rows, "flat_clusters": flat_rows, "changers": changers,
            "unassigned": rowsu, "channel_agreement": ca, "bin_minutes": BIN_MIN,
            "summary": summary, "out_dir": OUT, "tag": tag}

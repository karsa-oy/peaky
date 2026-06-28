# Peaky — Time-series clustering (co-varying families)

This document explains **how a batch time series is turned into co-varying ion
families** — the `clusters_changing` / `clusters_changers` / `clusters_flat`
figures and tables. It is a module deep-dive companion to
[`ARCHITECTURE.md`](ARCHITECTURE.md) (the whole pipeline),
[`ASSIGNMENT.md`](ASSIGNMENT.md) (how the formulas these channels carry are
assigned), and [`OUTPUTS.md`](OUTPUTS.md) (where each artifact lands).

**Code:** `peaky/batch/clustering.py` (`cluster_batch`, orchestration) and
`peaky/batch/cluster.py` (the pure engine + all constants). Trace building lives
in `peaky/reporting/analyte_viz.py` (`ion_traces`).

> Keep this in sync with the code. Every threshold below is a named constant in
> `cluster.py`; if you change one there, change it here.

---

## 1. What this stage does

After assignment, each batch sample has a peak list with committed formulas. This
stage asks a purely **temporal** question: *which assigned ion channels rise and
fall together over the run?* Channels that move as a group ("co-vary") are
reported as **families**; channels that don't move are pushed to a flat-background
overview. No chemistry is used here — only the shapes of the traces in time. It is
deliberately **event-agnostic** (it knows nothing about when an experiment's event
happened), so the same logic works for any batch.

```
merged ledger (M0 rows)  +  per-sample peak time series
        │
        ▼   one trace per ion channel (formula|adduct), native sample cadence
   per-channel raw traces ──► gate ──► log-shape correlation ──► cluster
                                              │
                         merge near-duplicates│
                                              ▼
                            flat-family demotion ──► dynamic families  (clusters_changing)
                                              ├──────► big standalone movers (clusters_changers)
                                              └──────► everything else      (clusters_flat)
```

---

## 2. Inputs: from ledger to per-channel traces

The unit of clustering is the **ion channel** = `neutral_formula | adduct`.
Channels of one neutral are kept **separate** — the `[M+H]⁺`, `[M+NH₄]⁺` and
`[M+(CH₄N₂O)H]⁺` channels of a compound often diverge in time, so summing them
would blend divergent shapes.

For each channel `ion_traces(ts, …, mode="raw", bin_minutes=None)` produces a
**raw trace**: one point per sample, placed at the sample's real timestamp — **no
re-gridding** onto a uniform lattice (re-gridding aliased into spurious empty
bins). Where a peak fell below a sample's detection floor the value is **NaN**
("below detection", not missing-at-random).

Two per-channel statistics are computed up front:

- `med` = `np.nanmedian` of the trace — **median of the *detected* points only**.
- `cv` = std / mean over finite positive points.

---

## 3. The transformation, stage by stage

All thresholds are the named constants from `cluster.py` (see §4).

1. **Persistence filter** (`goodk`). A channel must have **≥ `MIN_POINTS` (8)**
   finite points to be correlatable at all.

2. **Entry gate** — who is allowed into clustering:
   `med ≥ FLOOR (200 cps)` **and** the channel is not a Si-contaminant.
   There is **no cv gate** at entry — the assigned set is clustered on *shape*,
   and cv is only reported, not used to admit/reject.
   *(Optional: `cluster_batch(gate="episode", min_run=3)` swaps the median floor
   for "detected in ≥ `min_run` consecutive bins". Default is `gate="median"`,
   which reproduces historical output exactly.)*

3. **Correlate** (`cluster.correlate`). Clip each trace to the smallest positive
   value, take `log10`, then a **Pearson correlation matrix** (`min_periods=8`)
   over the **raw** log-traces. Correlation is on shape; the assigned path does
   **not** normalize to TIC/reagent (that is only the unassigned path, §7).

4. **Cluster** (`cluster.cluster`). Distance = **1 − r, signed** (so anti-phase
   channels stay far apart, not folded together), symmetrized, zero diagonal.
   **Complete-linkage** hierarchical clustering (average/single linkage chains
   every trace into one blob — observed on Br⁻), cut at distance
   `DIST_T (0.40)` ⇒ members share **r > 0.60**. Keep groups with
   **≥ `MIN_MEMBERS` (3)**.

5. **Merge near-duplicates** (`merge_similar`). A decay-dominated batch
   over-splits the dominant shape into many near-identical clusters. Take each
   cluster's **centroid** (mean log-trace), correlate centroids, and
   complete-linkage merge any whose centroids correlate at **≥ `MERGE_R` (0.85)**
   — folding duplicates together while keeping genuinely distinct shapes apart.

6. **Flat-family demotion** (`split_flat_clusters`). For each surviving family
   compute `cluster_flatness` = smoothed (window `SMOOTH_W` = 3)
   **max ÷ median of the member-mean** trace. If **< `FLAT_CLUSTER_RANGE` (1.4)**
   the family is demoted to flat background ("members correlate but the family as
   a whole doesn't move"). Otherwise it is a **dynamic co-varying family**.
   ⚠ This is a pure **magnitude** test — it measures *how much* the family moves,
   **not when** (see §8).

7. **Shape label** (`shape_of`). On the z-scored family mean, compare
   `mean(first 6)` vs `mean(last 6)` with gap 0.5 → `rise` / `fall` / `peak`,
   plus the `peak_hour`.

8. **Big standalone changers** (`big_changers`). Channels in the *remainder*
   (didn't join a family) whose smoothed max/median is **≥ `BIG_CHANGE_FOLD`
   (3.0)** are surfaced individually — a large lone change with no co-movers.

9. **Flat background.** Everything left — the uncorrelated remainder, the demoted
   flat families, and bright Si contaminants — is bunched into the flat overview.

---

## 4. Constants reference

All in `peaky/batch/cluster.py` (entry floor in `clustering.py`).

| constant | value | role |
| --- | --- | --- |
| `MIN_POINTS` | 8 | finite trace points required to correlate (persistence) |
| `FLOOR_DEFAULT` | 200 cps | entry brightness floor — `nanmedian` of the channel |
| `DIST_T` | 0.40 | clustering cut: `1 − r`, so members share **r > 0.60** |
| `MIN_MEMBERS` | 3 | smallest reported family |
| `MERGE_R` | 0.85 | centroid correlation to merge duplicate-shape families |
| `FLAT_CLUSTER_RANGE` | 1.4 | family-mean max/median below this → demote to flat |
| `BIG_CHANGE_FOLD` | 3.0 | lone-channel smoothed max/median to be a "big changer" |
| `SMOOTH_W` | 3 | smoothing window for max/median (rejects 1-bin spikes) |
| `CHANGING` / `FLAT_CV` | 0.30 | cv "varying" threshold (used only on the unassigned path) |
| `PEAK_RANGE` | 1.7 | transient-burst max/median (used only on the unassigned path) |
| `LINK` | `complete` | linkage method (both primary and merge passes) |

---

## 5. Metrics, defined

- **`cv` (`trace_cv`)** — std/mean over finite>0 points. Catches *sustained*
  variation; nearly blind to a 1–2 bin spike.
- **`_smoothed_range` / `trace_dynamic_range`** — smoothed `max / median` of a
  trace; ~1 = flat, a coherent burst scores high. Smoothing first kills single-bin
  noise.
- **`cluster_flatness`** — `_smoothed_range` of the **member-mean** trace; ~1 means
  the members co-vary but the family is flat (correlated background).
- **`panel_median`** — for the rendered bold line: NaN holes are below-detection
  LOWS, so they are filled to the family's detection floor before a *plain* median
  (`nanmedian` would drop the holes → survivorship bias → the line wrongly stays
  flat through a zeroing event).
- **`shape_of`** — start-vs-end of the z-scored mean → rise/fall/peak.

---

## 6. Outputs

| artifact | content |
| --- | --- |
| `tables/clusters_changing_<tag>.csv` + `.xlsx` | dynamic families (post-merge, post-demotion); one row per member channel with `cluster` id, `cv`, `median_cps` |
| `figures/clusters_changing_<tag>_p*.png` | family panels (raw cps, log y) |
| `tables/clusters_changers_<tag>.csv` | big standalone movers (`fold`, `peak_hour`) |
| `tables/clusters_flat_<tag>.csv` | bunched flat background (`cluster=0`) |
| `clusters_summary.json` | the gate values + funnel counts the PDF report documents |

---

## 7. The unassigned-peak path (variant)

The same engine clusters TS bins that match **no** assigned species ("unexplained").
Differences from the assigned path: brightness floor **50 cps** (median); a
**pre-cluster `split_varying` gate** — a bin is clustered only if `cv ≥ 0.30`
(`FLAT_CV`) **or** smoothed max/median **≥ `PEAK_RANGE` 1.7** (a brief synchronized
spike barely moves cv, so the burst term catches it); and correlation is on
**reagent- or TIC-normalized** traces, not raw. Flat bins are bunched.

---

## 8. Properties, invariants & gotchas

- **Brightness gate uses `nanmedian`** (median of *detected* samples), so it is
  blind to detection *sparsity*: a channel detected in only 1–2 bright bins can
  pass `med ≥ 200`. The `episode` gate option exists to key on temporal coherence
  instead.
- **Dynamic-vs-flat is magnitude-only, not timing.** A family that drifts gently
  at the *end* of a run (not at any event) can still clear
  `FLAT_CLUSTER_RANGE` and be labelled "rise". This is intentional — the stage is
  universal and event-agnostic; do **not** special-case a known event here.
- **Complete linkage + signed distance** are load-bearing: average/single linkage
  collapses everything into one family; unsigned distance folds anti-correlated
  channels together.
- **Channel notation uses Unicode glyphs**: the deprotonation channel is `−H⁻`
  with a **Unicode minus (U+2212)**, not ASCII `-`; adducts carry superscript
  charge and subscript counts (`+NH₄⁺`). When joining cluster tables to other
  tables, canonicalize (strip `⁺⁻`, map subscripts→digits, `−`→`-`,
  `(CH4N2O)H`→`Ur`) or matches silently fail.
- **Native cadence by default** (`bin_minutes=None`): traces use one point per
  sample at its real time. Pass `bin_minutes=int` to time-bin instead.

---

## 9. Code map

| function (`cluster.py` unless noted) | role |
| --- | --- |
| `clustering.cluster_batch` | orchestrates the whole stage; writes tables/figures/summary |
| `analyte_viz.ion_traces` | builds the per-channel raw traces |
| `correlate` | log10 → Pearson matrix |
| `cluster` | complete-linkage cut → labelled clusters ≥ 3 |
| `merge_similar` | centroid merge of duplicate-shape families |
| `cluster_rows` | assemble (id, members, r̄, shape, peak_hour) |
| `split_flat_clusters` / `cluster_flatness` | demote flat families |
| `big_changers` / `trace_dynamic_range` | lone large movers |
| `split_varying` / `trace_varies` | unassigned-path varying gate |
| `shape_of`, `panel_median`, `smooth` | labelling + rendering helpers |

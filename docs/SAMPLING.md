# Peaky — Sampling (which real samples get assigned)

This document explains **how a many-sample batch is reduced to a small subset to
assign** — the two selection strategies (time-representative and
brightest-coverage) and the exact arithmetic each uses. It is a module deep-dive
companion to [`ARCHITECTURE.md`](ARCHITECTURE.md) (the whole pipeline, the *Whole
batch* flow), [`MERGE.md`](MERGE.md) (which combines the per-sample ledgers this
selects), and [`TIMESERIES.md`](TIMESERIES.md) (whose `build_matrix` the
brightest-coverage strategy reuses).

**Code:** `peaky/batch/sampling.py`. Pure pandas/numpy, no network. The selected
`sample_item_id`s feed `assign.run` one at a time — `match_compounds` is
per-sample, so a synthetic union spectrum can't be scored.

> Keep this in sync with the code. Every threshold below is a named constant or a
> literal in `sampling.py`; if you change one there, change it here.

---

## 1. What this stage does

`assign.run` scores **one** sample, so an analyte present only during an event
window is never a candidate when the chosen file predates the spike. To get a peak
list representative of the *whole* batch, assign a small fixed subset and merge it.
Two strategies pick *which* real samples — both feed the identical assign → merge
chain; only the selection changes, never the scoring or the merge.

- **`representative`** (default, *THE RULE*) — `N_TIME` (5) evenly **time**-spaced
  samples (endpoints always included) **+ the single max-TIC sample**.
- **`brightest`** — bin all batch peaks by m/z and assign each significant bin's
  **brightest** sample (greedy set-cover). A *coverage* play: on reagent-CIMS the
  max-TIC pick is dominated by the reagent ion and is brightest for only a small
  fraction of analyte peaks.

```
batch peaks / per-sample table
        │  sample_table → one row/sample: id, time, tic, n_peaks
        ▼
 ┌──────────────── representative (THE RULE) ───────────────┐
 │ n_time targets = linspace(t_first, t_last, 5)            │
 │   nearest DISTINCT sample to each target  → time-grid    │
 │ + argmax(tic) sample                       → max-TIC     │
 └──────────────────────────────────────────────────────────┘
 ┌──────────────── brightest-coverage ──────────────────────┐
 │ build_matrix → samples × m/z bins                        │
 │ significant bins: per-bin max height ≥ height_floor 1000 │
 │ each bin's argmax sample = its winner (sets are DISJOINT) │
 │ greedy: take winners by bins-won desc until covered ≥0.85 │
 │   bounded k_min(6) ≤ n ≤ k_max(10); pad with richest TIC │
 │ + the two time endpoints (include_time_grid)             │
 └──────────────────────────────────────────────────────────┘
        ▼
 selected rows: role + (bins_won) → sample_item_id list → assign each
```

---

## 2. Inputs

- A **batch table**: either a per-**peak** table (has `height` → tic = Σ heights),
  or a per-**sample** table such as `samples.list` (already has `tic`, one row per
  sample, no peak load needed). Only `sample_item_id` is strictly required.
- For brightest-coverage: the per-**peak** table (height per peak) is mandatory —
  the matrix needs heights.

---

## 3. The transformation, stage by stage

1. **`sample_table`.** Collapse to one row per sample: `sample_item_id`,
   `datetime_utc`, `sample_item_name`, `tic`, `n_peaks`. Per-sample input is used
   directly; per-peak input is grouped (`tic = sum(height)`, `n_peaks = size`).
   Sorted by time when a clock is present.

### Representative (THE RULE)

2. **Time-grid pick.** With a clock and `n > n_time`: targets =
   `np.linspace(t_first, t_last, n_time)` over **ns-since-epoch**; for each target
   take the **nearest DISTINCT** sample (`argsort(|t − target|)`, skip already
   used). Selecting in **time** (not row index) is the point — an irregular run
   (dense early, a lone late file) still places a pick on the sparse region, and
   the endpoints are always included. Too few samples / no clock → take all.

3. **Max-TIC.** Add `argmax(tic)` (the richest single spectrum). Role becomes
   `time-grid+max-TIC` if it coincides with a grid pick, else `max-TIC`.

### Brightest-coverage

4. **Build the m/z matrix** (`TS.build_matrix`) → samples × m/z bins (see
   [`TIMESERIES.md`](TIMESERIES.md)).

5. **Significant bins.** Per-bin max height across samples; a bin is *significant*
   if its max ≥ **`height_floor` (1000 cps)** (reagent-relative — lower for a
   quieter dataset).

6. **Winners.** Each significant bin's `idxmax` is the sample where it is
   brightest; `value_counts` ⇒ `sample → bins_won` (descending). Because each bin
   has exactly **one** brightest sample, the winner sets are **disjoint**, so the
   optimal cover is the trivial greedy: take winners in `bins_won` order.

7. **Greedy cover.** Accumulate winners until **`coverage_target` (0.85)** of
   significant bins is covered, bounded **`k_min` (N_TIME+1 = 6) ≤ n ≤ `k_max`
   (10)** (`k_min = min(k_min, k_max)` — `k_max` is the hard cap). If still under
   `k_min`, **pad** with the richest (max-TIC) remaining samples so a too-high
   floor never under-selects.

8. **Time endpoints.** `include_time_grid` unions the first + last sample (cheap
   insurance the run start/end are represented even if they win no bins); role
   `coverage+time-grid` or `time-grid`, `bins_won = 0`.

---

## 4. Constants reference

All in `peaky/batch/sampling.py`.

| constant | value | role |
| --- | --- | --- |
| `N_TIME` | 5 | evenly-time-spaced samples in THE RULE |
| `select_representative_samples` `include_max_tic` | True | add the max-TIC sample |
| `select_brightest_coverage_samples` `coverage_target` | 0.85 | fraction of significant bins to cover before stopping |
| `k_max` | 10 | hard cap on assigned samples |
| `k_min` | `N_TIME + 1` = 6 | floor (clamped to `k_max`); pad target |
| `height_floor` | 1000.0 cps | a bin is "significant" if its max height ≥ this |
| `include_time_grid` | True | union the two time endpoints into the winner set |

---

## 5. Metrics, defined

- **`tic`** — total ion current per sample: Σ peak heights (per-peak input) or the
  server's own `tic` column (per-sample input). The "richest spectrum" ranking.
- **`bins_won`** — significant m/z bins a sample is the *brightest* for; the
  greedy cover's ranking key and the audit column in `selected_samples.csv`.
- **coverage fraction** — `Σ bins_won(selected) / max(#significant bins, 1)`;
  compared to `coverage_target` to decide when to stop.
- **`n_peaks`** — peak count per sample (informational, on the table).

---

## 6. Outputs

| artifact | content |
| --- | --- |
| `select_representative_samples` | selected `sample_table` rows, time-ordered, `role ∈ {time-grid, max-TIC, time-grid+max-TIC}` |
| `select_brightest_coverage_samples` | selected rows + `role ∈ {coverage-winner, time-grid, coverage+time-grid}` and `bins_won` |
| `*_sample_ids` | just the `sample_item_id`s (time order) for `assign_batch.run` |
| `tables/selected_samples.csv` | the chosen subset written by `assign_batch.run` (role + bins_won) |

---

## 7. Properties, invariants & gotchas

- **Selection in time, not index.** Row-index thinning over-samples a dense early
  region and skips a sparse late one; `linspace` over real timestamps fixes both.
- **Max-TIC ≠ analyte coverage on CIMS.** The reagent ion dominates TIC, so the
  max-TIC sample is brightest for only ~13 % of analyte bins (measured, orange Br⁻
  batch) — that is the whole reason brightest-coverage exists. It is a coverage
  play (more analyte peaks), not a speed play (~`k_max` assigns either way).
- **Disjoint winner sets** make set-cover trivial — there is no NP-hard subset
  search, just a sort by `bins_won`.
- **Padding guarantees ≥ k_min.** A high floor or a quiet dataset can leave few
  winners; the richest-TIC pad keeps the subset from collapsing.
- **Both strategies feed the same downstream.** Only *which* real samples get
  assigned changes — the per-sample scoring and the by-m/z merge
  ([`MERGE.md`](MERGE.md)) are identical.
- **Few-samples shortcut.** `n == 0` → empty; `n ≤ k_min` (brightest) or
  `n ≤ n_time` (representative) → take all.

---

## 8. Code map

| function | role |
| --- | --- |
| `sample_table` | batch table → one row/sample (id, time, tic, n_peaks) |
| `select_representative_samples` | THE RULE: time-grid + max-TIC subset (with roles) |
| `select_brightest_coverage_samples` | greedy brightest-bin set-cover subset (with `bins_won`) |
| `select_representative_sample_ids` / `select_brightest_coverage_sample_ids` | id-list convenience wrappers |
| `timeseries.build_matrix` (reused) | the samples × m/z-bin matrix the coverage strategy bins on |

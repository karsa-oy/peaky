# Peaky — Scoring (in-process IsoSpec + `score_pattern`)

This document explains **how a candidate neutral formula gets a match score** —
the in-process backend that replaced the network `match_compounds` round-trip.
It is a module deep-dive companion to [`ARCHITECTURE.md`](ARCHITECTURE.md) (the
whole pipeline), [`DATA_IO.md`](DATA_IO.md) (which dispatches here from
`score_candidates`), and [`ASSIGNMENT.md`](ASSIGNMENT.md) (which consumes the
scores). The cardinal rule still holds: **Mascope is the scorer.** This module
runs Mascope's *own* released maths (`mascope_tools.composition`) locally; Peaky
never invents a mass or isotope score.

**Code:** `peaky/io/local_scoring.py` (`score_candidates_local`, the engine);
the dispatch + the `PEAKY_LOCAL_SCORING` switch live in
`peaky/io/io_mascope.py` (`score_candidates` / `_score_candidates_local`).

> Keep this in sync with the code. Every threshold below is a named constant in
> `local_scoring.py`; if you change one there, change it here.

---

## 1. What this stage does

Given a sample's peaks and a list of candidate **neutral** formulas, score each
`(neutral × adduct)` ion against the spectrum and return one row per predicted
isotopologue — the **same schema** `io_mascope.flatten_match_tree` produces, so
the passes downstream don't know or care which backend ran.

The motivation is structural: the backend `match_compounds` is a *deep-annotation*
primitive — per `(formula × adduct)` it computes the full theoretical envelope and
returns the whole matched-**and-unmatched** tree, `O(candidates × adducts ×
envelope)` work and tens of thousands of rows for an `O(matches)` signal. That
drove the timeouts and OOM. The local path computes the identical envelope with
**IsoSpec** (`predict_isotopes`), scores it with the identical
**`score_pattern`** (`0.6·mass + 0.2·pattern + 0.2·intensity`), and emits **only
matched isotopologues** — no network, no 30k-row trees.

```
sample peaks (mz, height, peak_id)        candidate neutral formulas
        │  dedup on peak_id, sort by mz            │
        └──────────────┬───────────────────────────┘
                       ▼   for each (neutral × adduct)
        combine_formula_and_ionization → ion formula
                       ▼
        predict_isotopes (IsoSpec) → (pred_mz, pred_int, labels)
                       ▼   pred_rel = pred_int / pred_int[0]
        match M0 in ±ppm window  ── not found ──► drop candidate
                       │ found
                       ▼   match each isotope; keep if intensity-error ≤ tol
        score_pattern(obs_mz, obs_mz_err, obs_int, obs_int_err, pred_rel)
                       ▼   one score per ion, copied onto every iso row
        flat per-isotopologue rows  (matched: peak id + ppm; else None)
```

---

## 2. Inputs

- **`peaks`** — the sample's raw peaks. Only `mz` / `height` / `peak_id` are read;
  rows are **deduped on `peak_id`** (raw server peaks repeat per match) and sorted
  by m/z so the window search can use `np.searchsorted`.
- **`formulas`** — candidate neutral formulas (the grid + cheminfo union).
- **channels** — either peaky `adducts` (labels like `[M+Br]-`) or already-resolved
  mascope **`mechanisms`** strings (`+Br-`); the dispatcher passes the latter, which
  skips `adduct_to_mech`. Each is parsed once via `utils.parse_ionization`.

---

## 3. The transformation, stage by stage

All thresholds are the named constants from `local_scoring.py` (see §4).

1. **Prepare peaks.** Keep `[mz, height, peak_id]`, drop null m/z,
   `drop_duplicates(peak_id)`, sort by m/z → parallel `mzs` / `ints` / `pids`
   arrays.

2. **Build the ion** (`utils.combine_formula_and_ionization`). neutral + parsed
   ionization → the ion formula string (e.g. `C6H12BrO6-`). `adduct_to_mech`
   (when adducts, not mechanisms, are passed) collapses multi-part adducts by
   concatenating the added pieces: `[M+HBr+Br]-` → `+HBrBr-` (= +HBr₂).

3. **Predict the envelope** (`predict_isotopes`, IsoSpec). → `pred_mz`,
   `pred_int`, `labels` for the charged ion. Normalize to the base:
   **`pred_rel = pred_int / pred_int[0]`**. Empty envelope → skip.

4. **Match the monoisotopic base (i = 0).** Window half-width
   `d = pred_mz · ppm · 1e-6` with **`ppm` = `MATCH_MZ_TOLERANCE_PPM` (5.0)**;
   `searchsorted` for `[mz−d, mz+d]`, take the **closest** peak. **If M0 is not
   detected, the candidate is dropped entirely** (`base_int is None` → not a
   candidate at all) — there is no scoring an isotope envelope with no anchor.

5. **Match each heavier isotope (i > 0).** Same window + closest peak, but a peak
   is **accepted only if its intensity matches**: relative observed
   `rel_obs = ints[k]/base_int`, intensity error
   `ierr = |pred_rel[i] − rel_obs| / pred_rel[i]`, kept iff
   **`ierr ≤ INTENSITY_TOLERANCE` (0.4)** (= mascope_tools
   `ISOTOPE_MATCHING_INTENSITY_TOLERANCE`). A peak in the mass window with the
   wrong abundance is **not** attributed — it stays unmatched.

6. **Score the ion** (`score_pattern`). One score per ion from the matched
   pattern: `0.6·(mass term) + 0.2·(pattern term) + 0.2·(intensity term)`. The
   same `compound_score` / `ion_score` is copied onto every isotopologue row of
   that ion.

7. **Categorize** (`_category`). `score ≥ PROBABLE_THRESHOLD (0.8)` → `probable`;
   `≥ POSSIBLE_THRESHOLD (0.4)` → `possible`; else `unlikely`. These mirror the
   network scorer's `probable/possible_match_threshold`.

8. **Emit rows.** One row per predicted isotopologue. Matched rows carry the
   attributed `sample_peak_id`, `sample_peak_mz/intensity`, and a real
   `ppm_error`; unmatched isotopologues carry `None` for peak id, score, and ppm
   (so the envelope is described but the gaps are honest).

---

## 4. Constants reference

All in `peaky/io/local_scoring.py`.

| constant | value | role |
| --- | --- | --- |
| `PROBABLE_THRESHOLD` | 0.8 | score ≥ → `probable` (matches `DEFAULT_MATCH_PARAMS`) |
| `POSSIBLE_THRESHOLD` | 0.4 | score ≥ → `possible`; below → `unlikely` |
| `MATCH_MZ_TOLERANCE_PPM` | 5.0 | half-window for matching a predicted line to a peak |
| `INTENSITY_TOLERANCE` | 0.4 | max relative abundance error to attribute an isotope (= `ISOTOPE_MATCHING_INTENSITY_TOLERANCE`) |
| `score_pattern` weights | 0.6 / 0.2 / 0.2 | mass / pattern / intensity terms of the ion score (in `mascope_tools`) |

---

## 5. Metrics, defined

- **ion score** — `score_pattern(obs_mz, obs_mz_err, obs_int, obs_int_err,
  pred_rel)`: a single 0–1 number, **0.6 mass + 0.2 pattern + 0.2 intensity**.
  The whole arbitration downstream ranks on this.
- **`pred_rel`** — predicted isotope intensities normalized to the base
  (`pred_int / pred_int[0]`); the reference the observed pattern is judged against.
- **`ierr` (intensity error)** — `|pred_rel − rel_obs| / pred_rel`; the gate that
  decides whether a mass-window peak is *really* this isotope (≤ 0.4) or a
  coincidence.
- **`ppm_error`** — `(obs_mz − pred_mz)/pred_mz · 1e6`, populated only on matched
  rows.

---

## 6. Outputs

One DataFrame, **identical columns to `flatten_match_tree`** (so it is a true
drop-in):

| column group | content |
| --- | --- |
| compound | `compound_formula`, `compound_score`, `compound_category` |
| ion | `ion_formula`, `ion_score`, `ion_category`, `mechanism_id` (`im.mascope_notation`) |
| isotopologue | `isotope_formula`, `iso_label` (`M0` / `13C` / …), `is_base` (`i == 0`), `theo_mz`, `rel_abundance` |
| match | `iso_score`/`iso_category` (None if unmatched), `sample_peak_id`, `sample_peak_mz`, `sample_peak_intensity`, `ppm_error`, `abundance_error` (None for the base) |

`io_mascope._score_candidates_local` stamps `frame.attrs` with
`match_batches = 0`, `match_batch_failures = []`, `match_formulas = len(formulas)`
so callers see a uniform shape across both backends.

---

## 7. Properties, invariants & gotchas

- **Same maths, same schema — only the *source* moves.** The scores are Mascope's
  authored code (`mascope_tools`, public PyPI, same authors), executed locally and
  pinned to a library version. This keeps the "Mascope is the scorer" invariant.
- **No M0, no candidate.** An envelope whose monoisotopic line isn't detected is
  dropped before scoring — there is nothing to anchor relative abundances to.
- **Intensity-gated isotope attribution.** Being inside the ±5 ppm window is *not*
  enough; the peak's abundance must match the prediction within 40 %, or it stays
  unmatched. This is what stops a dense spectrum from "confirming" every isotope
  by mass coincidence.
- **Only matched isotopologues are scored data.** Unmatched envelope lines are
  still emitted (so the predicted pattern is visible) but with `None` scores/ppm —
  never a fabricated zero.
- **Deterministic.** No network and no RNG, so a re-run over the same inputs is
  byte-identical — this is what lets the default scorer satisfy
  `test_determinism.py` (the `match_compounds` fallback re-introduces server-side
  variation). See [`ARCHITECTURE.md §7`](ARCHITECTURE.md#7-reproducibility--provenance).
- **Mixed +/- adducts are unsupported** by `adduct_to_mech` (it raises); a
  multi-part adduct must be all-add or all-subtract.
- **The switch is opt-out, not opt-in.** Local scoring is the default;
  `PEAKY_LOCAL_SCORING=0/false/no/off` falls back to the network path (see
  [`DATA_IO.md §3`](DATA_IO.md#3-the-transformation-stage-by-stage) and
  [`MASCOPE_TOOLS_INTEGRATION.md`](MASCOPE_TOOLS_INTEGRATION.md)).

---

## 8. Code map

| function (`local_scoring.py` unless noted) | role |
| --- | --- |
| `score_candidates_local` | the engine: peaks + formulas → flat per-isotopologue table |
| `adduct_to_mech` | peaky adduct label → mascope mechanism string (multi-part collapse) |
| `_category` | score → `probable` / `possible` / `unlikely` |
| `utils.parse_ionization` / `combine_formula_and_ionization` (mascope_tools) | parse channel; build the ion formula |
| `predict_isotopes` (mascope_tools, IsoSpec) | theoretical isotope envelope of the ion |
| `score_pattern` (mascope_tools) | the 0.6/0.2/0.2 mass/pattern/intensity ion score |
| `io_mascope.score_candidates` | backend dispatcher (local default ↔ `match_compounds`) |
| `io_mascope._score_candidates_local` / `_local_scoring_enabled` | local bridge (peaks from cache, mechanism names) + `PEAKY_LOCAL_SCORING` switch |

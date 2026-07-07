# Peaky — Data I/O (Mascope → peak / time-series DataFrames)

This document explains **how raw Mascope data becomes the DataFrames the rest of
the pipeline consumes** — the peak table, the per-sample time series, the
candidate-formula lists, and the flat per-isotopologue score table. It is a
module deep-dive companion to [`ARCHITECTURE.md`](ARCHITECTURE.md) (the whole
pipeline), [`SCORING.md`](SCORING.md) (the maths `score_candidates` dispatches
to), and [`OUTPUTS.md`](OUTPUTS.md) (where each artifact lands).

**Code:** `peaky/io/io_mascope.py` — the **only** module that talks to Mascope.
Everything here wraps the `mascope-sdk` `MascopeClient`; nothing else in the
package opens a socket.

> Keep this in sync with the code. Every threshold below is a named constant or a
> literal in `io_mascope.py`; if you change one there, change it here.

---

## 1. What this stage does

This is the boundary between Mascope (the system of record) and Peaky (the
analysis layer). It does five jobs and **no chemistry**: connect to a server,
pull peaks (single sample) or a whole-batch time series, infer which adducts the
sample was ionized through, estimate a coarse mass offset, enumerate candidate
neutral formulas for an m/z, and turn Mascope's score tree into a flat table. The
one genuine *number* transform here is `flatten_match_tree` (tree → per-row ppm);
the heavy scoring maths live in [`SCORING.md`](SCORING.md).

```
~/.mascope/.env  (MASCOPE_URL + MASCOPE_ACCESS_TOKEN)
      │  connect()  — search precedence, never the MCP's stale in-memory token
      ▼
  MascopeClient
      │
 ┌────┴──────────────────────────────┐
 │ single sample                      │ whole batch
 │ fetch_peaks(sample_id)             │ fetch_batch_samples (one row/sample, no peaks)
 │   get_peaks(matches=True)          │ fetch_batch_peaks → load_peaks(dataset,
 │   cache → ~/.mascope-assign-cache  │   escape_batch(batch))  [legacy per-sample fallback]
 ▼                                    ▼
 raw peak table ──► detect_adducts ──► adduct list   ([M-H]⁻ fallback)
      │         ──► estimate_offset ──► seed ppm offset (median, |ppm|≤10)
      ▼
 query_candidates(mz, mechanism_ids, …) ──► neutral-formula candidates (cheminfo)
      │                                       (flaky endpoint → degrade to [])
      ▼
 score_candidates ──► flatten_match_tree ──► flat per-isotopologue table
   (local default / match_compounds opt-in)   (one row per compound·ion·isotopologue)
```

---

## 2. Inputs

- A **server + token** from a `.env` (or process env): `MASCOPE_URL`,
  `MASCOPE_ACCESS_TOKEN`, optional `MASCOPE_WORKSPACE`.
- A **`sample_id`** (single-sample assignment) or a **`(dataset, batch)`** pair
  (batch time series). `dataset` maps to a *workspace name* on a legacy server.
- For scoring: a list of candidate **neutral formulas** and resolved
  **`mechanism_ids`** (from `resolve_mechanism_ids`).

---

## 3. The transformation, stage by stage

1. **Connect** (`connect`). Resolve a `.env` by **precedence**: explicit
   `--env` → `$MASCOPE_ENV` → the `ENV_SEARCH` list (repo-root `.env` →
   `./.env` → `~/.mascope/.env` → `~/mascope-mcp/.env` →
   `~/.claude/skills/mascope-sdk/.env`) → a cwd walk-up (`find_dotenv`) →
   `CANONICAL_ENV`. The file is read from **disk every time** — the long-running
   MCP server holds a stale in-memory token and 401s, so the live file wins.
   Multi-workspace servers require `--workspace`; a legacy build (no
   `/api/datasets`) is tolerated by `_patch_datasets_list_for_legacy_servers`
   (the constructor's `datasets.list()` health-check 404 is degraded to `None`).

2. **Fetch peaks — single sample** (`fetch_peaks`). `get_peaks(matches=True)`
   pulls the raw peak table *with Mascope's own matches flattened in* (so it is
   multi-row-per-peak — dedup is the ledger's job). Cached to
   `~/.mascope-assign-cache/<sample_id>/peaks.parquet` (CSV fallback if parquet
   write fails); `use_cache=True` reuses it.

3. **Fetch time series — whole batch** (`fetch_batch_peaks`). Uses the SDK
   `load_peaks(dataset=, batches=escape_batch(batch), confirm_above=None)` — one
   row per (sample × peak), enough for the TS/cluster layer. `confirm_above=None`
   never prompts (batches exceed 100 samples). On a legacy server (no
   `/api/datasets`) it falls back to `_legacy_load_batch_peaks`, which resolves
   the `sample_batch_id` from the raw `sample/batches` endpoint then runs the
   SDK's own per-sample fetch loop (`matches=False` there — the TS only needs
   mz/height/datetime/peak-id). `fetch_batch_samples` returns one row per sample
   (id, name, `datetime_utc`, `tic`, polarity) **without** loading peaks — enough
   for representative-sample selection.

4. **Detect adducts** (`detect_adducts`). Read the distinct
   `ionization_mechanism` values off the sample's own matches and reverse-map
   them through `MECH_TO_ADDUCT` (the inverse of `ADDUCT_TO_MECH`). This is what
   makes a Br-CIMS sample get `[M+Br]⁻` offered instead of forcing Br into the
   neutral. **Fallback `["[M-H]-"]`** if nothing is recognized.

5. **Estimate offset** (`estimate_offset`). A coarse **median ppm** mass offset
   from the sample's *own* server matches — base ions only (heavy-isotope rows,
   detected by a `[` in `target_isotope_formula`, are skipped). Each match's ppm
   is `(mz − ion_mz)/ion_mz·1e6`; a **gross-outlier guard** drops `|ppm| > 10`;
   needs **≥ `min_n` (8)** survivors or returns `None`. This only *seeds* the
   pre-calibration gates (pass 0's `|ppm| ≤ 2` known-species gate would be blind
   to a −1.9 ppm instrument); pass-1 self-calibration is the authoritative fit.

6. **Enumerate candidates** (`query_candidates` / `query_candidates_bulk`).
   `cheminfo.query_by_mz` returns candidate **neutral** formulas for one m/z at
   **`ppm` 5.0**, **`limit` 25**, deduped. cheminfo is flaky (timeouts / 500s);
   any error **degrades to `[]`** — the local grid covers the same CHO/CHON space,
   so a failure never kills the run. The bulk variant fans out over **`workers`
   12** threads.

7. **Score candidates** (`score_candidates`). Dispatcher. Default is the
   **in-process** backend (`_score_candidates_local` → `local_scoring`, see
   [`SCORING.md`](SCORING.md)); `PEAKY_LOCAL_SCORING=0` (or `false/no/off`) opts
   back to the network `match_compounds`. The network path batches formulas at
   **`MATCH_BATCH` 200** (it times out above ~500), scores chunks concurrently on
   **`MATCH_WORKERS` 5** threads, coerces `mz_tolerance` to **int** ppm, and (by
   default) **raises if any batch fails** — a partial candidate universe is worse
   than a failed pass. Mechanism ids are reverse-mapped to mascope mechanism
   names by `_mechanism_names`, which fixes the deprotonation sign (`-H+` →
   `-H-`).

8. **Flatten the score tree** (`flatten_match_tree`, **pure**). One row per
   (compound · ion · isotopologue), columns listed in §6. `parse_isotope_label`
   turns `[13C]C3H5O2-` → `('13C', is_base=False)` and a base ion → `('M0',
   True)`. **ppm is computed only for a genuinely matched isotope** (a real
   attributed peak with `sample_peak_id`, positive intensity, and a real
   `sample_peak_mz`); unmatched/forced nodes carry `sample_peak_mz == theo` and
   zero intensity, so their `ppm_error` is left `None`. For a **¹⁵N-labelled
   nitrate** reagent (`^N` in the ion formula) the base is re-anchored
   (`_reanchor_labelled_reagent`) off the phantom all-light M0 onto the real
   single-¹⁵N line (**`delta` 0.997035**, tolerance 0.01).

---

## 4. Constants reference

All in `peaky/io/io_mascope.py`.

| constant | value | role |
| --- | --- | --- |
| `ENV_SEARCH` | repo `.env` → `./.env` → `~/.mascope/.env` → `~/mascope-mcp/.env` → `~/.claude/skills/mascope-sdk/.env` | credential-file search list (after `--env` / `$MASCOPE_ENV`) |
| `CACHE_ROOT` | `~/.mascope-assign-cache` | per-sample peak cache (`<sid>/peaks.parquet`) |
| `MATCH_BATCH` | 200 | formulas per `match_compounds` call (times out above ~500) |
| `MATCH_WORKERS` | 5 | concurrent `match_compounds` batches (I/O-bound, server-safe) |
| `DEFAULT_MATCH_PARAMS.mz_tolerance` | 5 (**int** ppm) | mass window for the network scorer |
| `DEFAULT_MATCH_PARAMS.isotope_ratio_tolerance` | 0.2 | per-isotope abundance tolerance |
| `DEFAULT_MATCH_PARAMS.min_isotope_abundance` | 0.15 | smallest predicted isotope considered |
| `DEFAULT_MATCH_PARAMS.min_isotope_correlation` | 0.7 | isotope-pattern correlation floor |
| `DEFAULT_MATCH_PARAMS.probable_match_threshold` | 0.8 | score ≥ → `probable` |
| `DEFAULT_MATCH_PARAMS.possible_match_threshold` | 0.4 | score ≥ → `possible` |
| `query_candidates` `ppm` / `limit` | 5.0 / 25 | cheminfo enumeration window + cap per m/z |
| `query_candidates_bulk` `workers` | 12 | cheminfo fan-out threads |
| `estimate_offset` `min_n` | 8 | min matches to trust an offset (else `None`) |
| `estimate_offset` outlier guard | `\|ppm\| ≤ 10` | drops gross mismatches before the median |
| `_reanchor_labelled_reagent` `delta` | 0.997035 | ¹⁴N→¹⁵N mass step (tol 0.01) |

---

## 5. Metrics, defined

- **`ppm_error`** — `(sample_peak_mz − theo_mz)/theo_mz · 1e6`, **only for a
  genuinely matched isotope**; `None` for forced/unmatched nodes (their
  `sample_peak_mz` is the theoretical mz, which would falsely read 0 ppm).
- **`rel_abundance`** — the predicted isotope's abundance relative to the base.
- **`abundance_error`** — server `match_abundance_error` for the isotope row.
- **offset (`estimate_offset`)** — the *median* of the per-match ppm over base
  ions, after the `|ppm| ≤ 10` guard; a single coarse number seeding pre-cal gates.

---

## 6. Outputs

| artifact | content |
| --- | --- |
| raw peak table (`fetch_peaks`) | full multi-row-per-peak frame (peaks + Mascope's own matches); cached `~/.mascope-assign-cache/<sid>/peaks.parquet` |
| batch time series (`fetch_batch_peaks`) | one row per (sample × peak): mz / height / `datetime_utc` / sample ids; optional `save_path` parquet |
| per-sample table (`fetch_batch_samples`) | one row per sample: `sample_item_id` / name / `datetime_utc` / `tic` / polarity |
| adduct list (`detect_adducts`) | the adduct labels the sample was ionized through (or `["[M-H]-"]`) |
| offset (`estimate_offset`) | a single median-ppm float, or `None` |
| candidate lists (`query_candidates*`) | neutral formulas per m/z |
| flat score table (`flatten_match_tree`) | one row per compound·ion·isotopologue: `compound_formula/score/category`, `ion_formula/score/category`, `mechanism_id`, `isotope_formula`, `iso_label`, `is_base`, `theo_mz`, `rel_abundance`, `iso_score/category`, `sample_peak_id/mz/intensity`, `ppm_error`, `abundance_error` |

---

## 7. Properties, invariants & gotchas

- **`flatten_match_tree` is pure** (no network) and unit-tested against a captured
  fixture — the contract the offline test suite locks.
- **Batch names are matched as a case-insensitive REGEX** (`str.contains`). A
  literal name with metacharacters — `Sample run (Ur+ CIMS)`, a `^Nitrate`
  prefix — silently matches nothing unless run through `escape_batch` (`re.escape`).
- **`-H+` is not a cation.** The server names deprotonation `-H+` (the *removed*
  species' sign), but it yields an anion. `_mechanism_names` normalizes the
  trailing sign to the mechanism's polarity (`-H+` → `-H-`) before handing it to
  the scorer, or the entire `[M-H]⁻` channel is silently dropped.
- **ppm is meaningful only for matched isotopes.** Forced/phantom nodes get
  `ppm_error = None`, never a misleading 0.
- **¹⁵N nitrate phantom base.** The server models the reagent N as
  natural-abundance, tagging the signal-less all-light line as M0 and the real
  100%-¹⁵N line as a `15N` isotopologue; `_reanchor_labelled_reagent` moves
  `is_base` onto the real line so the passes (which commit only `is_base`) see it.
- **Legacy servers.** No `datasets` concept — `dataset` is treated as a workspace
  name, `sample/batches` ignores its dataset filter, and the datasets health-check
  is patched to degrade gracefully. All `list_*` / fetch helpers carry a fallback.
- **cheminfo failures are non-fatal** — they degrade to `[]`; the local grid
  carries that m/z, so coverage is preserved.
- **`mz_tolerance` is integer ppm** on the network path (coerced with `round`).
- **Caching is content-blind**: `fetch_peaks` keys only on `sample_id`. A
  re-fetch after server-side re-processing needs `use_cache=False`.

---

## 8. Code map

| function | role |
| --- | --- |
| `connect` | build a `MascopeClient` from the resolved `.env` (precedence above) |
| `_find_env` / `_patch_datasets_list_for_legacy_servers` | credential search; legacy health-check tolerance |
| `list_workspaces` / `list_datasets` / `list_batches` | discovery (with legacy fallbacks) |
| `resolve_batch_id` / `escape_batch` | legacy batch-name → id; regex-safe name escaping |
| `fetch_peaks` | single-sample raw peaks (+ matches), cached parquet |
| `fetch_batch_peaks` / `_legacy_load_batch_peaks` | whole-batch peak time series (modern + legacy) |
| `fetch_batch_samples` | per-sample table (no peaks) for sample selection |
| `resolve_mechanism_ids` / `detect_adducts` | mechanism name↔id; infer adduct system |
| `estimate_offset` | coarse median-ppm offset to seed pre-cal gates |
| `query_candidates` / `query_candidates_bulk` | cheminfo neutral-formula enumeration |
| `score_candidates` | scoring dispatcher (local default ↔ `match_compounds`) |
| `_score_candidates_local` / `_local_scoring_enabled` | local backend bridge + the `PEAKY_LOCAL_SCORING` switch |
| `_mechanism_names` | id → mascope mechanism name, deprotonation-sign fix |
| `flatten_match_tree` | **pure** score-tree → flat per-isotopologue table |
| `parse_isotope_label` / `_reanchor_labelled_reagent` | isotope-label parse; ¹⁵N base re-anchor |

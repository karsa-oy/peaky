# Local scoring via `mascope_tools` (WIP)

Status: **spike validated, foundation built, not yet wired as default.** Goal: replace
the network `match_compounds` hot loop with in-process `mascope_tools.composition`,
which is the same Mascope scoring maths run locally (IsoSpec + `score_pattern`).

## Why

`match_compounds` is a deep-annotation endpoint: per `(formula x adduct)` it computes
the full isotope envelope and returns the whole compound->ion->isotopologue tree
(matched *and* unmatched). Measured on the Bromide run, **pass 1 alone returned
36k-141k isotopologue rows per sample** to assign ~600-1300 peaks. That `O(candidates
x adducts x envelope)` payload is what forced the 300s timeout and OOM'd the heavier
Uronium run (146k peaks). peaky only needs an `O(matches)` screen (M0 verdict + a
bounded isotope score).

## What was found / done

- `mascope_tools.composition` (public PyPI `mascope_tools`, same authors) provides the
  local equivalents: `find_compositions` (enumerate), `predict_isotopes` (IsoSpec
  envelope), `match_isotopic_pattern` / `score_pattern` (0.6*mass + 0.2*pattern +
  0.2*intensity, 0-1), plus Senior/valence/element-ratio heuristics.
- **Spike (`assign_compositions`)**: full 873-peak sample in **17.9s vs ~5-6min**
  backend. Its built-in arbitration is naive (picks pure-carbon `C14`), so peaky must
  keep its own arbitration and use the library only for enumerate + score.
- **Score parity is "similar"** (sufficient per spec): local vs backend on real ions
  e.g. `Br2-` 0.975 vs 1.0, `Br3-` 0.938 vs 0.994.
- **15N nitrate adduct MASS fixed** in `mascope_tools` (was missing): `parse_ionization`
  could not mass `+^NO3-` (pyteomics can't mass the `^N` symbol). Added
  `utils.CARET_ISOTOPES = {"^N": "N[15]"}` + `utils.composition_mass()`; `+^NO3-` now
  masses to 62.985401, 15N-14N nitrate delta = 0.997035 (exact). *File:
  `mascope/libraries/tools/src/mascope_tools/composition/utils.py` — uncommitted in
  that repo.* **This fixes only the MASS — 15N isotope-pattern SCORING is still broken;
  see "Known limitation: 15N" below.**
- **`peaky/local_scoring.py` built**: `score_candidates_local(peaks, formulas, adducts)`
  returns peaky's exact `flatten_match_tree` schema, looping `predict_isotopes` +
  `score_pattern` per candidate (no per-peak cap). Validated on `6B76`: schema matches,
  assignments are chemically sensible (`C7H9O5-` 0.999, `C9H13O2-` 0.999, tiny ppm),
  782 peaks matched an M0. Adduct label -> mech conversion handles Br/Ur/CO3/uronium/15N.

## Parity evaluation (10 samples: 6 Bromide + 4 Uronium, 6140 ion-pairs)

Harness: `scripts/eval_local_scoring.py <run_dir>...`. Compares local scoring to the
backend assignment in each per-file ledger (the `role=="M0"` rows). Two tests: (A)
score the backend's exact assigned ions locally and compare scores; (B) local
score + argmax over the same candidate pool vs the backend's assigned neutral.

| | assignment agreement | coverage | score MAE | local scoring time |
|---|---|---|---|---|
| Bromide | **0.91** | 0.92 | 0.070 | 0.1-0.7 s/sample |
| Uronium | **0.98** | 1.00 | 0.087 | 0.2-1.6 s/sample |

(vs ~5-6 min/sample on the backend.) Key findings:

- **Adduct universe must match the passes**, not just the base profile: peaky's passes
  add `[M+NH4]+` (37% of Uronium M0s!) and `[M+CO3]-`/`[M+HBr+Br]-`/`[M+HBr+CO3]-`
  (Br). With those included, coverage went 0.66->1.00 (Ur) and 0.84->0.92 (Br).
  `adduct_to_mech` now collapses multi-part adducts.
- **Score gap is a clean systematic offset**, not noise: local is **+0.067 (Br) /
  +0.074 (Ur)** higher than the backend (median ~= mean, std ~0.07). Offsets are
  rank-preserving, which is why argmax agreement is high.
- **Disagreements are benign**: Br 2/298 (local picks a *more* plausible formula where
  the backend chose a low-scoring exotic one); Ur 31/1149 (near-ties, N-containing
  alternative winning by ~0.005).

**Tuning conclusion:** parity is "close enough" — argmax agreement 0.91/0.98 across both
reagents. The match score itself is indicative only; the local score is ~0.07
optimistic vs the backend. So when switching to local, **shift peaky's tier thresholds
by ~+0.07** (probable 0.8->~0.87, possible 0.4->~0.47) rather than recalibrating the
(honest) library score. No score-weight tuning needed.

## Remaining to make it the default (next steps)

1. **Wire dispatch** in `io_mascope.score_candidates`: behind `PEAKY_LOCAL_SCORING`
   env (or config), call `local_scoring.score_candidates_local(<sample peaks>, chunk,
   adducts)` instead of `client.matching.match_compounds(...)` + `flatten_match_tree`.
   The sample peaks are already available (`fetch_peaks`).
2. **Keep candidate generation bounded** — use peaky's existing `query_candidates`
   (cheminfo, limit=25/peak, context-filtered ~294/pass) OR `find_compositions` with a
   `max_result_rows` cap. Uncapped enumeration produced 11,805 candidates/sample (slow);
   peaky's bounded set scores in ~1s.
3. **Threshold recalibration** — local `score_pattern` is close but not identical to the
   backend; re-check `probable=0.8 / possible=0.4` against the pass behaviour. Parity is
   "similar" so expect minor tuning, not a rewrite.
4. **Validate full pipeline** — run all 6 passes with local scoring on the Bromide
   subset, diff assignments vs the committed backend run
   (`peaky-output/2025-08-11-Bromide-...T120446Z/`).
5. **15N-labelled nitrate scoring** — currently broken (mass-only); see "Known
   limitation: 15N" below for the three tested gaps and the fix needed.
6. **Dependency** — add `mascope-tools` to peaky `pyproject.toml` (public PyPI); during
   dev use a `[tool.uv.sources]` editable override on the local checkout (like mascope-sdk).
7. **Optional memory win** — emit only matched isotopologues (+ M0) instead of the full
   predicted envelope to shrink the per-sample frame further.

## Known limitation: 15N-labelled nitrate scoring (NOT handled)

The 15N-nitrate reagent (`NO3_15N`, adduct `[M+^NO3]-`) is **NOT correctly scored
locally** — only the adduct mass is fixed. This is NOT a regression (it is a new path;
the backend still scores 15N via its own re-anchor), but it must be closed before local
scoring can replace the backend for 15N work. Three layered, tested gaps:

1. **Ion construction fails.** `combine_formula_and_ionization(neutral, +^NO3-)` ->
   `PyteomicsError: Invalid formula: O3^N` — pyteomics rejects the `^N` symbol, so the
   ion is never formed.
2. **IsoSpec can't express a heavy isotope via the formula string.**
   `IsoThreshold(formula="...N[15]...")` -> `ValueError: garbage inside "[15]"`. No
   isotope notation is accepted through the formula string at all.
3. **The 98/2 distribution can't be modelled that way anyway.** A named isotope is 100%
   pure; the reagent is **98% 15N / 2% 14N**. Partial labelling is only expressible via
   IsoSpec's low-level custom-abundance API (`IsoSpecPy.Iso` with explicit
   `isotopeMasses` / `isotopeProbabilities`).

**Fix (in `mascope_tools`, not peaky):** `predict_isotopes` needs a labelled-atom path
that (a) detects the `^N` carried by the adduct, (b) builds the distribution via the
custom-abundance `Iso` API — natural atoms at natural abundance + the labelled N as
`{15N: 0.98, 14N: 0.02}` — yielding the 15N base (98%) AND the 14N satellite (2%, ~0.997
Da below). The purity (0.98) should be a parameter (a reagent property), not hardcoded.
Note: the backend does not model 98/2 either (it treats the reagent N as natural
abundance, then peaky `_reanchor_labelled_reagent` re-anchors), so a correct local path
would be MORE accurate than the server.

## Quick repro

```python
from mascope_tools.composition import CompositionSearchConfig
from mascope_tools.composition.finder import find_compositions
from peaky.local_scoring import score_candidates_local
# peaks: DataFrame with mz, height, peak_id
flat = score_candidates_local(peaks, neutral_formulas, ["[M+Br]-", "[M-H]-"])
```

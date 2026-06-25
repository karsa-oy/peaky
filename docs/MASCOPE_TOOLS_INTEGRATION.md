# Local scoring via `mascope_tools` (WIP)

Status: **wired (opt-in `PEAKY_LOCAL_SCORING=1`) and full-pipeline-validated on Bromide;
one integration bug found & fixed; threshold recalibration is the remaining open call.**
Goal: replace the network `match_compounds` hot loop with in-process
`mascope_tools.composition`, the same Mascope scoring maths run locally (IsoSpec +
`score_pattern`).

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
- **15N-labelled nitrate fully handled** (`mascope_tools` >= 2026.06.25): `+^NO3-`
  masses correctly (62.985401) AND `predict_isotopes` produces the labelled envelope
  (15N base + ~2% 14N satellite) via a self-contained custom-element convolution, with
  per-reagent purity (`NO3_15N.purity = 0.98`). Validated by a molmass parity test.
  Merged to mascope `develop` and released on PyPI.
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

## Full-pipeline validation (Bromide, 6 representative samples)

Ran the SAME batch end-to-end three ways on today's code and diffed final per-file M0
assignments peak-by-peak (`scripts/diff_runs.py`, keyed on `sample_item_id|peak_id`):

- **server path** (`match_compounds`): ~35 min, 1774 M0 peaks — the baseline.
- **local path, first cut**: ~18 min but only 1591 M0 peaks, **0.874** agreement.
  Diagnosed a **−10% coverage regression**: 486 server-confident peaks (eff≥0.8) left
  `unexplained`, dominated by the `[M-H]-` channel (66 local vs **582** server M0s).
- **local path, after the fix**: ~8–18 min, 2149 M0 peaks, **0.932** agreement,
  `[M-H]-` fully recovered (594 M0s).

**The bug (fixed):** the server names deprotonation `-H+` ("remove a proton H⁺") with
`ionization_mechanism_polarity='-'`, but it yields an ANION. The dispatch reverse-maps
mechanism_ids straight to that name and hands it to `parse_ionization`, which reads the
trailing `+` as net charge **+1** — so every `[M-H]-` candidate was scored as a +1
cation, predicted the wrong m/z and matched nothing, silently dropping the whole
deprotonation channel (37% of Bromide assignments). The parity eval never hit this
because it used `adduct_to_mech` (correct `-H-`); only the pipeline dispatch uses
`_mechanism_names`. **Fix:** `_mechanism_names` normalises the trailing sign to the
mechanism's polarity (`-H+`→`-H-`); consistent adduct names are untouched. Guarded by
`tests/test_mechanism_names.py`. (Positive-mode Uronium — `+H+`/`+NH4+` — was never
affected, consistent with its 0.98 eval agreement.)

**Post-fix residual differences are benign:** local assigns *more* than the server
(2149 vs 1774), but the 599 local-only M0s are mostly real (median eff 0.845, 46%≥0.8;
server had them `unexplained`/`iso_child`), and the 224 local misses are the server's
*weak* ones (median eff 0.670) that local mostly re-reads as `iso_child`. The
disagreement set is now M0-vs-iso_child interpretation calls plus the systematic
**+0.052 score offset** (local runs optimistic, same as the eval's +0.07).

**Open call — threshold recalibration:** the +0.052 offset inflates the Assigned tier
(1739 vs 1361). Either de-bias the local `eff_score` by ~−0.05 to the server scale (one
spot, keeps all existing gates/margins) OR shift peaky's tier thresholds +0.05 for the
local path. Both are cosmetic for *which* peaks get found (the extra local peaks are
mostly real); they matter for tier-as-confidence parity with historical server runs.

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
5. **15N-labelled nitrate scoring** — DONE in `mascope_tools` >= 2026.06.25
   (custom-element `predict_isotopes` + per-reagent purity).
6. **Dependency** — add `mascope-tools` to peaky `pyproject.toml` (public PyPI); during
   dev use a `[tool.uv.sources]` editable override on the local checkout (like mascope-sdk).
7. **Optional memory win** — emit only matched isotopologues (+ M0) instead of the full
   predicted envelope to shrink the per-sample frame further.

## 15N-labelled nitrate scoring (resolved in mascope_tools >= 2026.06.25)

Previously broken (mass-only). `predict_isotopes` now handles labelled `^X` custom
elements self-contained: it computes the base (non-custom) envelope with IsoSpec, then
convolves it with the labelled element's distribution (multinomial for `^N2`, Cartesian
for several), yielding the **15N base + ~2% 14N satellite ~0.997 Da below**. Purity is a
parameter (`predict_isotopes(..., purity=...)`), carried per reagent
(`NO3_15N.purity = 0.98`); `combine_formula_and_ionization` forms the labelled ion via
the `^N`-aware `parse_composition`. A molmass parity test guards against drift.

The Mascope backend still has its own (molmass-based) `predict_isotopes`; converging it
onto this shared implementation and retiring the molmass fork is tracked as a separate
issue. (Historically the backend modelled the reagent N as natural abundance and peaky's
`_reanchor_labelled_reagent` re-anchored; the local path is now the more direct route.)

## Quick repro

```python
from mascope_tools.composition import CompositionSearchConfig
from mascope_tools.composition.finder import find_compositions
from peaky.local_scoring import score_candidates_local
# peaks: DataFrame with mz, height, peak_id
flat = score_candidates_local(peaks, neutral_formulas, ["[M+Br]-", "[M-H]-"])
```

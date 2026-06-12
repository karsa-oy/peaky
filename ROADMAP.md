# Roadmap — improving the low-quality / unassigned tail

Status 2026-06-11: the pipeline reaches ~95% "signal explained" on the Br-CIMS
test sample, but a large share of the tail is **Low-confidence, high-ppm
(3–4 ppm) heteroatom-rich assignments that are probably wrong** (e.g. deprotonated
bromo-organics). The user flagged this. The conclusion below reframes the goal
and lists the concrete methods to fix it. **Nothing here is implemented yet** —
this is the agreed direction for the next session.

## The core problem

At high m/z with heteroatoms (Br/N/S/O/F) enabled, the density of candidate
formulas within any fixed ±ppm window is large, so "best fit within tolerance"
stops being an identification — it's the closest of many. High ppm error is the
*symptom*; formula-space density is the disease. And the objective we optimised,
**`% signal explained`, actively rewards false positives** because it pushes the
pipeline to formula-fit the residual at any cost.

## The reframe (most important)

Stop maximising coverage. Optimise for a **defensible identification core** plus
an **honest map of what cannot be known**:
- `identified` — unique, sub-σ mass, isotope-consistent, ideally series-supported
- `ambiguous (N candidates)` — report the alternatives, pick none
- `isotope-partner only` — explained as a satellite, no independent formula
- `below assignability` — m/z too high for a unique formula at this accuracy
- `noise / unexplained`

A trustworthy ~70–80% beats a padded 95%. Expect the headline number to DROP;
that is the progress.

## Methods, in priority order

1. **Self-calibration + z-score mass gate (highest leverage, single-spectrum).**
   Fit μ, σ of the ppm error over the High/Good CHO/CHON backbone (real
   instrument accuracy is a tight cluster, e.g. +0.2 ± 0.4 ppm). Judge every
   candidate by `z = |ppm − μ| / σ`, not raw ppm. Accept within ~2–3σ of the
   *calibration*. A 4-ppm peak when μ=0.2, σ=0.5 is ~7σ out → reject. This alone
   deletes the bromo-organic tail while keeping the backbone. Slots into
   `arbitrate()` — no pass-structure change.

2. **Isotope envelope as a carbon clamp (biggest density-killer).** The ¹³C₁/M
   intensity ratio is a carbon counter (±1 at high m/z). Use the quantitative
   envelope (already in Mascope's `match_score_isotope` / abundances) as a HARD
   constraint that shrinks the candidate grid *before* scoring, not as a post-hoc
   boost. Clamping carbon collapses candidate counts ~10× and removes most
   heteroatom guessing because the mass budget closes.

3. **Candidate density = confidence.** Count formulas surviving the calibrated
   window per peak; 1 → confident, several → `ambiguous` (report the list). This
   is the honest confidence we currently fake with score thresholds.

4. **Spectrum-as-network / belief propagation.** Don't score peaks
   independently. Start from the backbone; let a formula propagate to a neighbour
   only via an exact mass difference (CH₂, O, …) with consistent mass error and
   agreeing isotope pattern. Peaks that can't connect stay unassigned. Refuses
   isolated high-ppm heteroatom formulas structurally. (Generalises Pass 2/4.)

5. **Assignability ceiling.** For the calibrated accuracy + element set, compute
   the m/z above which no formula is unique, and declare peaks above it
   `below-assignability` rather than fitting through them.

6. **Time-resolved / multi-sample correlation (strongest for the tail; needs
   data).** The SDK has `get_peak_timeseries`. Real compounds' isotopologues and
   homologs co-vary across time/samples; false formula-fits on noise don't.
   Co-variation is orthogonal evidence mass accuracy can't provide — it can
   confirm or kill the Low-confidence tail independent of ppm. **User is sourcing
   time-resolved data.** When available, build a correlation-based confirmer:
   require a candidate's isotopologues (and ideally series neighbours) to
   correlate above a threshold before promoting confidence.

## Two-tier reporting

Never present a candidate as an identification. The Excel/markdown should split
`Identifications` (tier-1) from `Candidates` (tier-2, with full alternative lists
and the reason none wins). Add a `tier` and `n_candidates` column to the ledger.

## Suggested next-session order

1. Self-calibration + z-score gate (1) — fast, immediately cleans the tail; show
   before/after on `<sample-id>`, expect ~95% → ~80% *real*.
2. Carbon clamp (2) — compounding density reduction.
3. Density/tier reporting (3 + two-tier) — make confidence honest.
4. If time-resolved data arrives: correlation confirmer (6) — the real fix for
   the junk that survives mass-based filtering.
5. Network propagation (4) and assignability ceiling (5) as the deeper rebuild.

Treat contaminant families (Pass 3) as opt-in / flagged, not coverage-seeking.

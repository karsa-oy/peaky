# Roadmap — state after v28 (2026-06-12) and what comes next

## Where the pipeline stands

v28 on the reference sample `<sample-id>` (Br-CIMS, ambient air):
**261 M0 / 56.5% peaks / 89.3% signal explained, 21/21 flagships, 0 junk,
ledger clean, ~180 s.** Run `python3 scripts/check_flagships.py <ledger.csv>`
after ANY change — it asserts the validated identifications (TFA both
channels, hydroxy-acid ladder, monoterpene rungs + Br partners, silanediol
series, HO2/HNO3/HNO2/HNO4) and the banned junk classes.

Pipeline shape: pass 0 (known species: contaminant series + atmospheric
acids/radicals, locked, twin-gated) → pass 1 (CHO/CHON backbone + calibration)
→ pass 2 (GKA series) → pass 3 (evidence-opened contaminant families) →
pass 4 (residual: iso-pairs incl. BrCl, deep series, carbon-clamped) →
pass 5 (known-neutral completion: cross-channel + series gaps) →
isotope-physics audit → calibrated mass-gate audit.

## What the residual IS (eliminated dimensions — do not re-chase)

- **Sulfate (SO3)**: dead, 0.06x decoys. **Nitrate (NO3/HNO3/NO2)**: dead,
  0.0–0.1x decoys (the sample's N lives in the small acids, now assigned).
- **Fluorine**: real but small (TFA + telomer); O>6 junk class banned.
- **Silicon**: solved (silanediol ladder, ~115k cps, locked contaminant).
- **Chlorine**: faint (Cl-H at 6x decoys), consistent with the known BrCl
  mixed-halogen family — a thread, not a dimension.
- **What remains**: ~40 bright unknowns in multi-halogen C/H lattices
  (CH2 93x, clean C2H4/H2/C2 ladders, Br/BrCl-twinned, heavy defect).
  PROVEN unsolvable from the sum spectrum alone: no CHNOSPSi+I formula closes
  under the carbon clamp in 7 adduct frames; candidates ambiguous.

## Next big improvements, in priority order

1. **Time-resolved correlation confirmer (THE unlock; user sourcing data).**
   SDK hook exists (`get_peak_timeseries`). Build: co-variation clustering of
   the residual; members of one lattice family must correlate. Resolves
   composites (393/395-style), confirms/kills the 36 Low/Suspect, separates
   inlet contamination (flat) from ambient chemistry (variable), and names
   whole lattices from any single identified member. The ledger already
   carries halogen counts, carbon brackets, and lattice memberships as
   correlation groups.

2. **Two-tier reporting: Identified / Candidates / Below-assignability.**
   Stop presenting one formula per peak. Tier rules already de facto exist
   (uniqueness in calibrated window, isotope consistency, alternatives list).
   The BrCl family (ClN2O10-type vs I-bearing candidates) and the surviving
   O15 monsters (C19H12O15, C23H18O15 — they sit ON the unexplained C/H
   lattice, so they are family members wearing CHO fantasies) belong in
   Candidates. Includes: candidate-density as the confidence currency.

3. **Below-assignability certificates (automate the manual proof).** For each
   bright has-constraints residual peak, run the clamped frame x element
   enumeration and stamp the result into the Unassigned sheet ("searched
   N frames x CHNOSPSi+I, C clamped a–b: zero fits → composite or exotic").

4. **Composite-doublet detection.** Twin ratio deviating >15% from isotope
   prediction → flag "composite: ~N cps hidden component"; marks peaks whose
   clamps are unreliable and feeds the time-series target list.

5. **mu(m/z) mass-dependent calibration.** ~30 z>2 stragglers remain; settle
   drift-vs-wrong-formula (evidence so far says wrong formulas: mixed signs
   at same m/z).

6. **Robustness niceties.** One batch-level re-retry before the (new,
   correct) fail-loud raise in score_candidates bites a whole pass on flaky
   server days; pass-4 `residual_ppm_*` could inherit the calibrated sigma;
   pin the X=76.5 half-integer "C/H lattice" view in the GKA widget; write
   lattice memberships into ledger commentary.

## Standing lessons (encode-don't-remember)

- Every evidence GATE needs a complementary RECOVERY path (TFA: chain gate
  stopped the flood and silently dropped the one real fluorochemical).
- Negative evidence (absent satellite) never overrides independent positive
  evidence (agreeing channel, twin satellite) — peak pickers lose peaks.
- A locked claim must pass self-consistency (own-twin gate) — composites.
- The server's scoring is authoritative; never trust hand ion-mass arithmetic
  (electron-mass sign error cost us TFA's [M-H]- twice).
- Coverage metrics reward fiction; the flagship list + count-based coverage
  are the honest scoreboard.
- GKA rotations: trust rows at integer A(R); verify any non-integer rotation
  (defect-quantum coincidences); half-integer X is a feature (C/H lattice).

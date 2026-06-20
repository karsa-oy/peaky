# mascope-peak-assign — developer / iteration guide

From-scratch, test-driven successor to `mascope-formula-assignment`. This dir is
the canonical home: edit here, run tests here, ship here. (An earlier scratch
copy may still exist at `~/mascope-assign` — ignore/delete it; this is the one.)

## Layout

```
mascope-peak-assign/
  SKILL.md            invocable manifest + full usage (read this first)
  ROADMAP.md          the agreed next-step quality work (calibration, carbon clamp, ...)
  README.md           this file
  mascope_assign/     the package (see SKILL.md module map) — single-sample assign
                      + the batch pipeline: sampling / assign_batch / cluster /
                      analyte_viz (full VK) / pdf_report
  tests/              one test_<module>.py per module (26 files) + fixtures/match_tree.json
  scripts/
    run_assignment.py one-shot: pipeline -> csv/xlsx/md/json/html
    gka_widget.py     standalone interactive rotating-GKA from a ledger CSV
```

For a whole batch (5 time-spaced + max-TIC samples merged → clustering → PDF
report), see the "Representative-sample batch pipeline" section of SKILL.md and the
reference drivers in `~/mascope-output/orange-assign/` (run_orange / run_clusters /
run_vankrevelen / run_report).

## Install & run

```bash
pip install -e .                   # pulls mascope-sdk + deps; registers `mascope-assign`
cp .env.example ~/.mascope/.env    # then fill in MASCOPE_URL + MASCOPE_ACCESS_TOKEN

mascope-assign list datasets                              # discover your data
mascope-assign list batches  --dataset "<workspace>"
mascope-assign list samples  --batch "<batch>" --dataset "<workspace>"
mascope-assign assign --sample-id <ID> --reagent Br \
    --height-cutoff 100 --output-dir ~/mascope-output/<name>
```
`--reagent {auto,Br,Ur,…}` forces the analyte channels (a positive/sparse sample
otherwise mis-detects as negative). Heavy work runs on the host Python; a Mascope
token is read from `~/.mascope/.env` (or `--env` / `$MASCOPE_ENV`). `~5 min` for a
~1000-peak sample at cutoff 100. (`python3 scripts/run_assignment.py …` still works
as a thin forwarder; `python3 -m mascope_assign …` is equivalent to the script.)

## Test loop

```bash
python3 tests/test_smoke.py          # 2s "install OK" check (no creds, no network)
pytest tests/                        # or: for t in tests/test_*.py; do python3 "$t"; done
```
833 offline assertions across 30 files, no network. Live smoke for io_mascope:
`MASCOPE_LIVE=1 python3 tests/test_io_mascope.py`. **Rule: every code change
ships with a test; keep the suite green.** Tests use plain asserts and run as
scripts (exit non-zero on failure); each also exposes a validating `test_all`
so `pytest tests/` collects and passes them. CI runs the suite with no creds.

## Design invariants (don't regress)

- One ledger DataFrame, one row per peak; passes only fill/annotate. `ledger.py`
  enforces structural invariants on commit; `ledger.validate()` must return `[]`.
- Mascope is the only scorer (`io_mascope`). Other modules never call the network.
- Chemistry gates are structural (integer-DBE-on-neutral, Senior, O-cap,
  halogens-as-H). See SKILL.md "Chemistry rules".
- Heteroatoms enter the neutral only with positive evidence; relaxed filtering is
  "earned by evidence" (chain membership / isotope confirmation), never default.

## Current status (2026-06-16 — see ROADMAP.md for full state, lessons + next steps)

Two validated modes now:
- **Negative Br-CIMS** (reference) — `<sample-id>` 404'd on <server>; the same
  physical sample lives at `<sample-id>` (<dataset>,
  08:21, −1.99 ppm). Full run: **22/22 flagships + 0 junk** (offset-aware
  check_flagships), 263 M0 / 177 Identified.
- **Positive urea-CIMS** (`uronium` context, NEW) — `<sample-id>`. Full
  pipeline at the −2.45 ppm offset: **719 M0 (604 Id), signal explained 64%→81%
  after the PDMS/siloxane-ladder pass** (`siloxane.py`). Outputs in
  `~/mascope-output/uronium-v5/`. The session added: positive context + urea
  reagent library, offset-tolerant calibration (calibrate/confidence/relabel/
  tiers/arbitration + `estimate_offset` pre-cal), Br-pass polarity guards, the
  dedicated siloxane pass, and a multi-file "experiment" merge. **560+ tests green.**

Earlier reference numbers (v44, the −0.6 ppm Br copy), cutoff 100:
- 269 M0, **65.7% peaks / 91.6% signal explained**, 21/21 flagships, 0 junk,
  ledger clean, **446 offline tests green**. **Tiered: 170 Identified / 99
  Candidate** (`tiers.py`; now also a mass-error-distribution gate, a
  CO₃-background-channel gate, and degeneracy-awareness). Roles: M0 / iso_child /
  reagent / **artifact** (ringing) / unexplained. New modules `degeneracy.py`
  (honest cross-family degeneracy) + `cleanup.py` (isotope-confirmed recovery,
  bromide-cluster labelling, ringing-artifact flagging). Outputs archived
  per-version in `~/mascope-output/assign-dev/v*/`.
- ⚠️ A tested reagent-formula / `[81BrO]-` fix is in the working tree but NOT yet
  in an output — **re-run to v45 next session** (ROADMAP "0. PICK UP HERE").
- **Composite-peak detection (`detect_composites`)**: the M+1 region (13C/29Si)
  is halogen-free, so it scales only with the assigned compound; if observed M0
  exceeds the M+1-implied intensity, an unresolved co-eluting compound shares
  the m/z. Flags (does not demote) + reads the co-component halogen off the
  even-shift residual. The silanediol C8H26O5Si4 (393) scores only 35% in
  Mascope because it is ~45% co-eluting BrCl -- formula (Si4, proven by +74.019
  rung spacing) and prediction (binomial) are BOTH correct; the peak is mixed.
  n=2 is clean. New `composite_note` column, surfaced in Identified + ownership.
- **Isotope-envelope completion (`complete_isotope_envelopes`)**: predicts each
  committed ion's full M+1/M+2/M+4 envelope (`isotopes.isotope_pattern`,
  per-element convolution incl. Si/Br/Cl combos) and claims it — ~44% of the
  bright "residual"/Candidate peaks were isotope satellites of brighter peaks,
  not independent compounds. The 393/395 case: the silanediol Si4+Br M+2 at
  395 was mis-assigned a phantom `C8H12ClF6NO2S` because its M+4/M+2 ratio
  (~0.26) mimics a Cl doublet; now 395/397 are its envelope. Runs before pass 4
  (so satellites never reach the iso-pair stage) and post-audit. iso_child rose
  276→304; peaks-explained 58.6→60.6%.
- **Pass 6 ladder gap-fill (`ladders.py`)**: walks homolog/oxidation diagonals
  out from committed anchors (+O/+CO/+CO2/+CH2O/±CH2/+C2H4/-H2O, same adduct)
  and fills the gaps, gated hard against the false positives the adversarial
  diagonal analysis surfaced (fluorinated/Si contaminant ladders excluded,
  81Br/13C isotope-satellite guard, O<=min(C,9), unexplained-only, Candidate
  tier). Di-bromide [M+HBr+Br]- gaps need the +HBr pairing (bromine-free
  neutral also at [M+Br]-) as corroboration. Adds 7 SOA homolog/oxidation
  rungs. The diagonals are MOSTLY contaminants + isotope satellites, not SOA —
  only a handful of genuine biogenic-SOA ladders survive verification.
- **Di-bromide SOA clusters (the ex-"unsolvable C/H lattice")**: the bright
  n_Br=2 residual is biogenic SOA (mono-/sesquiterpene oxidation products)
  detected as `[M+HBr+Br]-` reagent clusters, NOT exotic organohalogens. 10
  now assigned (e.g. 409.0015 = C15H22O3). Enabled by registering the
  di-bromide frames + the user's server-side `+Br2-` mechanism, and a
  pre-pass-4 carbon-clamp that frees O15-monster-occupied peaks for re-claim.
- **Ambient acids un-buried**: `[Br1+acid]-` is the `[M+Br]-` analyte channel,
  so formic acid (232k cps), acetic, pinic etc. are now Identified analytes,
  not "reagent". Reagent role bucket 24 -> 11; even-n bare clusters (Br2-.)
  labelled.
- Pipeline is now 6 passes (0: known species, 1: backbone+calibration, 2: GKA,
  3: evidence-opened families, 4: residual iso-pairs/series, 5: completion)
  plus two post-run audits (isotope physics, calibrated mass gate).
- **Regression protection**: `python3 scripts/check_flagships.py <ledger.csv>`
  after every change. Git history in this directory is the change log.
- **Remaining frontier**: multi-halogen C/H-lattice families, unsolvable from
  the sum spectrum — blocked on time-resolved data. See `ROADMAP.md`.

## Performance notes

- Cost = `match_compounds` ≈ 3.8 s / 200-formula batch; batches run concurrently
  (`io_mascope.MATCH_WORKERS`). `chemistry._grid_cached` memoises grid
  enumeration (a missing cache once caused a 60× regression).
- Per-pass timing is logged by `assign._safe` and stored in the manifest.
- `cheminfo` is off by default (flaky/slow; grid is the primary enumerator).

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
  mascope_assign/     the 13-module package (see SKILL.md module map)
  tests/              one test_<module>.py per module + fixtures/match_tree.json
  scripts/
    run_assignment.py one-shot: pipeline -> csv/xlsx/md/json/html
    gka_widget.py     standalone interactive rotating-GKA from a ledger CSV
```

## Run it

```bash
cd ~/.claude/skills/mascope-peak-assign
python3 scripts/run_assignment.py --sample-id <ID> --context ambient-air \
    --height-cutoff 100 --output-dir ~/mascope-output/<name>
```
Needs host Python (has `mascope-sdk`) + `~/mascope-mcp/.env`. Run via the shell
MCP. `~5 min` for a ~1000-peak sample at cutoff 100.

## Test loop

```bash
for t in chemistry contexts ledger isotopes series_gka io_mascope reagents \
         passes residual report series_detect; do
  echo "== $t =="; python3 tests/test_$t.py || break
done
```
229 offline assertions, no network. Live smoke for io_mascope:
`MASCOPE_LIVE=1 python3 tests/test_io_mascope.py`. **Rule: every code change
ships with a test; keep the suite green.** Tests use plain asserts (no pytest),
exit non-zero on failure.

## Design invariants (don't regress)

- One ledger DataFrame, one row per peak; passes only fill/annotate. `ledger.py`
  enforces structural invariants on commit; `ledger.validate()` must return `[]`.
- Mascope is the only scorer (`io_mascope`). Other modules never call the network.
- Chemistry gates are structural (integer-DBE-on-neutral, Senior, O-cap,
  halogens-as-H). See SKILL.md "Chemistry rules".
- Heteroatoms enter the neutral only with positive evidence; relaxed filtering is
  "earned by evidence" (chain membership / isotope confirmation), never default.

## Current status (2026-06-13, v37/v38)

Test sample `<sample-id>` (Br-CIMS, atmospheric), cutoff 100:
- 267 M0, **60.6% peaks / 91.0% signal explained**, 21/21 flagships, ledger
  clean. **Tiered: 179 Identified / 88 Candidate** (`tiers.py`; mechanical
  rules, candidate-density currency, lattice-monster + BrCl demotions, same-ion
  decomposition aliases excluded). Excel is an 11-sheet styled workbook
  (Identified / Candidates-per-formula / evidence-characterized Unassigned +
  legend). Outputs archived per-version in `~/mascope-output/assign-dev/v*/`.
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

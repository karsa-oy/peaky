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

## Current status (2026-06-11)

Test sample `<sample-id>` (Br-CIMS, atmospheric), cutoff 100:
- ~95% signal explained, ~290 s, ledger clean. Outputs archived in
  `~/mascope-output/assign-dev/50BkQ_c100v13_*`.
- **Known limitation**: the Low-confidence / high-ppm tail contains likely-wrong
  heteroatom assignments. The fix is a reframe + new methods — see `ROADMAP.md`.
  That is the next session's work; the user is sourcing time-resolved data to
  enable the correlation-based confirmer.

## Performance notes

- Cost = `match_compounds` ≈ 3.8 s / 200-formula batch; batches run concurrently
  (`io_mascope.MATCH_WORKERS`). `chemistry._grid_cached` memoises grid
  enumeration (a missing cache once caused a 60× regression).
- Per-pass timing is logged by `assign._safe` and stored in the manifest.
- `cheminfo` is off by default (flaky/slow; grid is the primary enumerator).

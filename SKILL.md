---
name: mascope-peak-assign
description: >-
  Multi-pass chemical-formula assignment for high-resolution mass-spec peaks
  stored in Mascope. Use when a user asks to assign formulas, annotate a
  spectrum, identify compounds, build a target list, or explain unassigned /
  contaminant / homolog peaks for a Mascope sample_id. SDK-native, runs locally
  via the shell MCP; defers all mass + isotope scoring to Mascope's
  match_compounds; produces an 11-sheet tiered Excel (Identified / Candidates /
  below-assignability) with commentary, close alternatives, per-isotopologue
  scores, a peak-ownership audit, and an interactive rotating-GKA widget. Triggers: "assign formulas", "peak
  assignment", "what's in sample X", "annotate spectrum", "unassigned peaks",
  "Kendrick / GKA", "homologous series", "CIMS", "HOM", "PFAS / contaminants".
---

# Mascope multi-pass peak assignment

A reproducible, test-driven, SDK-native pipeline. **All heavy work runs locally
via `mcp__shell__run_command`** against the host Python (which has `mascope-sdk`);
the Mascope MCP is never used to transport peak tables through context.

This is the from-scratch successor to `mascope-formula-assignment`. Canonical
home and iteration repo: `~/.claude/skills/mascope-peak-assign/`.

## Operating principle

Division of labor is the core design decision:

- **Mascope is the scoring oracle.** `match_compounds` returns a
  compound → ion → isotopologue tree; every node carries its own `match_score`
  and the attributed `sample_peak_id`. We never invent a match score; we hand it
  candidate neutral formulas and read its per-isotopologue verdict.
- **We own candidate generation, chemistry plausibility, series logic, and
  arbitration** — i.e. *which* formulas are worth asking about and *which*
  answer wins each peak.

State lives in one mutable **ledger DataFrame** (one row per physical peak).
Passes are functions ledger→ledger; they fill and annotate, never drop rows. The
ledger's commit API enforces structural invariants so no pass can corrupt it.

## Pre-flight

1. Use `mcp__shell__run_command` (host Python), not the cowork sandbox.
2. `.env` at `~/mascope-mcp/.env` has `MASCOPE_URL` + `MASCOPE_ACCESS_TOKEN`
   (auto-loaded). `openpyxl` must import (`pip3 install --user openpyxl`).
3. Pick `--context` from the sample's setting (see below).

## Running

```bash
cd ~/.claude/skills/mascope-peak-assign
python3 scripts/run_assignment.py \
    --sample-id <ID> --context ambient-air \
    --height-cutoff 100 --output-dir ~/mascope-output/<name>
```

Writes `<ID>_<UTC>_{ledger.csv, assignments.xlsx, summary.md, manifest.json,
gka.html}` plus per-pass ledger checkpoints. Full run on a ~1000-peak Br-CIMS
sample is ~5 min (cutoff 100).

### Contexts
`ambient-air` (= atmospheric), `chamber`, `indoor-air`, `object-headspace`,
`combustion`, `water`, `food`, `none`. Context sets plausibility bounds + which
Pass-3 contaminant families are eligible. Reagent adducts are NOT set by context
— they are **detected from the sample** (`ionization_mechanism` column).

### Key flags
`--ppm` (m/z trust, default 1.0) · `--search-ppm` (enumeration tol, 5.0) ·
`--height-cutoff` (cps, 100) · `--no-pass2/3/4` · `--no-cache`.

## The pipeline

| pass | what it does |
|---|---|
| Pre | detect reagent adducts from the sample; prescan isotope fingerprint; **label reagent-ion clusters** (Brₙ, Brₙ·neutral, BrO) so they are never assignment candidates |
| 1 | lock the high-confidence **CHO/CHON backbone**: grid-enumerate, score with `match_compounds`, arbitrate (complexity-penalised, isotopologue-gated), commit M0 owners + attach Mascope's isotopologue children, lock High peaks |
| 2 | **iterative GKA series** expansion from locked anchors (CH₂/O/H₂O/CO/CO₂/C₂H₂O + siloxane/CF₂), chaining confirmed members as new anchors |
| 3 | **automatic series detection** (the machine "rotating plot") opens contaminant families on decoy-controlled evidence; HBr-cluster ladder; organosulfate/nitrate/siloxane/amine + iso-gated bromo/chloro-organics |
| 4 | **residual explainer**: isotope-pair resolution of ~1.998-Da doublets + deep 2-step series; DBE-only plausibility; ppm-disciplined acceptance |
| 5 | **known-neutral completion**: cross-channel partners + series gaps of passes 1–4 (no new formula space) |
| 6 | **anchored ladder gap-fill** (`ladders.py`): walk homolog/oxidation diagonals (+O/+CH₂/+CO₂/−H₂O, constant-DBE for carbon growth) out from Identified anchors, satellite-guarded, Candidate tier |
| iso-env | **isotope-envelope completion** (`complete_isotope_envelopes`, before pass 4 + post-audit): claim every committed peak's full predicted M+2/M+4 envelope (Si/Br/Cl combos), attaching unexplained satellites and **displacing weak M0s that are really a parent's satellite** — kills the ~44% of "residual" peaks that are isotope lines (the silanediol M+2 mis-read as a Cl-F-S organic) |
| audits | 13C carbon-clamp (pre-pass-4 + post), Br-doublet repair, calibrated mass gate |

## Chemistry rules (enforced + regression-tested)

- **DBE on the neutral**: non-negative INTEGER + Senior's rule. Half-integer is
  an ion-only artifact (deprotonation), so organic nitrates pass as neutrals.
- **Halogens count as hydrogens** everywhere: DBE, and the Van Krevelen ratio
  uses `(H+F+Cl+Br+I)/(C+Si)` (Si is a C-equivalent backbone atom).
- **Structural oxygen cap**: `O ≤ 2·(C+N+S+P)+4` (valence, not a Van Krevelen
  prior — kills `C3H5ClO17`-type mass-fits; real HOMs pass).
- **Reagent-halogen policy**: the reagent halogen's ION isotope can't prove it
  sits in the NEUTRAL (covalent `X(Br)[M-H]⁻` ≡ `Y·HBr·Br⁻`). The complexity
  prior on the reagent element is never waived by isotope confirmation, and
  anchor·HX clusters are resolved as clusters of the bare analyte.
- **Isotopologue-gated heteroatoms**: a neutral S/Cl/Br needs its Mascope-
  confirmed ³⁴S/³⁷Cl/⁸¹Br or its complexity skepticism is not waived.

## Outputs

| file | contents |
|---|---|
| `_ledger.csv` | every peak: role, formula, adduct, scores (incl. arbitration `eff_score`/`eff_margin`/`tied`), ppm, confidence, **tier + tier_reason + candidate_density**, provenance, commentary, alternatives, isotopologues |
| `_assignments.xlsx` | Summary · Read me (legend) · **Identified** · **Candidates** (one row per candidate formula) · Unassigned (evidence-characterized) · By class · Unique formulas · Isotopologues · Peak ownership (all peaks) · Target list · Reagent ions — styled: frozen headers, autofilters, number formats, tier/confidence color chips |
| `_summary.md` | narrative + top assignments + coverage |
| `_manifest.json` | module versions, prescan, series evidence table, per-pass timing |
| `_gka.html` | interactive rotating-GKA widget (see below) |

## Interactive rotating-GKA widget

`python3 scripts/gka_widget.py LEDGER.csv [-o out.html] [--ppm 2]` → a
self-contained HTML (no server). Slider rotates the scaling factor X so
homologous series flatten into horizontal rows; peaks colored by status
(backbone / low / unassigned). Band detector uses the mass-accuracy-derived
tolerance `δGKA ≈ (X/mass(R))·δm`, `δm = ppm·(m/z)·1e-6`. Use it to spot
structure (CF₂ contaminant ladders, oxidation series) the auto-detector did not
already open. `run_assignment.py` emits one per run.

## Module map (`mascope_assign/`)

| module | role |
|---|---|
| `chemistry.py` | masses, formula algebra, grid (integer-DBE / Senior / O-cap), complexity penalty, grid cache |
| `contexts.py` | context profiles + plausibility filter + contaminant families |
| `ledger.py` | the peak DataFrame + invariants + commit API |
| `io_mascope.py` | the ONLY Mascope I/O: peaks, cheminfo, parallel `match_compounds` + per-isotopologue parser, adduct detection |
| `isotopes.py` | prescan fingerprint → grid constraints; **`isotope_pattern()`** envelope predictor (per-element convolution) |
| `series_gka.py` | GKA/Kendrick math, repeat units, propagation |
| `ladders.py` | pass-6 anchored homolog/oxidation ladder gap-fill |
| `series_detect.py` | automatic decoy-controlled series detection + chain extraction |
| `reagents.py` | reagent-cluster library + labeler |
| `passes.py` | arbitration + the 4-pass director |
| `residual.py` | Pass 4 residual explainer |
| `tiers.py` | Identified/Candidate tiering (margin, density, lattice/BrCl demotions) |
| `ladders.py` | Pass 6 anchored homolog/oxidation-ladder gap-fill (diagonal SOA series) |
| `report.py` | Excel / markdown / sheets |
| `assign.py` | orchestrator + `PassConfig` + module manifest |

## Testing & iteration

`for t in chemistry contexts ledger isotopes series_gka io_mascope reagents
passes residual ladders tiers report series_detect; do python3 tests/test_$t.py; done`
— 390 offline assertions, no network (io_mascope live smoke gated behind
`MASCOPE_LIVE=1`). Every module has a matching `tests/test_<module>.py`. Add a
test with each change; keep the suite green. See `README.md` for the dev loop and
`ROADMAP.md` for the open quality work.

## Gotchas

- `match_compounds` (plural); integer `mz_tolerance`; batched at 200, now scored
  concurrently (5 workers).
- `cheminfo` is flaky/slow and OFF by default (`cfg.use_cheminfo`) — the local
  grid is the primary, complete enumerator. It only adds compound names.
- Extra channels (e.g. `+CO3-`) must be passed as explicit `mechanism_ids`;
  the server's auto-select only covers the sample's own channels.
- `% signal explained` is a coverage metric, not a quality metric — see ROADMAP.

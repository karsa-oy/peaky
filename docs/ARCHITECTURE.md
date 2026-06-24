# Peaky — Architecture

This document explains **how Peaky works** for someone reading or extending the
code. Companion docs: [`README.md`](../README.md) (install + dev loop),
[`QUICKSTART.md`](../QUICKSTART.md) (5-minute run), [`SKILL.md`](../SKILL.md)
(Claude-Code operating instructions), [`ASSIGNMENT.md`](ASSIGNMENT.md) (what
assignment produces, for a scientist), [`OUTPUTS.md`](OUTPUTS.md) (every artifact,
where + what), and [`ROADMAP.md`](ROADMAP.md) (development history).

---

## 1. Where Peaky sits

```
Mascope   →  data platform   (the app + database; system of record, the scorer)
Peaky     →  analysis layer   (this package: assignment, clustering, reporting)
Claude    →  the interface    (drives Peaky in plain language; never does chemistry)
```

Peaky reads high-resolution CIMS mass-spec peaks from **Mascope** (via
`mascope-sdk`), assigns chemical formulas to them, and turns the result into
tiered spreadsheets, figures, and a PDF report. Two hard rules shape everything
below:

1. **Mascope's `match_compounds` is the only scorer.** Peaky never invents a
   mass/isotope score. It enumerates candidate formulas, hands them to Mascope,
   and arbitrates among the scores Mascope returns.
2. **The chemistry gates are structural, not statistical.** Integer DBE, Senior's
   rule, the oxygen cap, evidence-gated heteroatoms/halogens — these are valence
   facts, applied identically every run. **No LLM is in the assignment loop**, which
   is why a run is reproducible and auditable.

---

## 2. The core design: one ledger

All pipeline state lives in **a single mutable pandas DataFrame — the ledger**
(`ledger.py`). One **row per observed peak**. Every stage is a
`ledger → ledger` function that only *fills or annotates* columns; nothing else
holds authoritative state.

A peak's row carries its **role**, and once it is explained, the formula/adduct
and the evidence behind it:

| column group        | meaning                                                              |
| ------------------- | ------------------------------------------------------------------- |
| `role`              | `M0` (monoisotopic owner) · `iso_child` · `reagent` · `artifact` · `unexplained` |
| formula / adduct    | the committed neutral + ionization adduct (`ion_formula` for reagent ions) |
| scores              | Mascope `match_score` plus arbitration `eff_score` / `eff_margin` / `tied` |
| evidence            | ppm error, confidence label, isotopologue children, alternatives, commentary |
| verdict             | `tier` + `tier_reason`, `candidate_density`, `degeneracy_density` / `degeneracy_note` |

Because state is one table, the invariant "every peak is in exactly one role and
M0 ownership is unique" is checkable at any point (`ledger.py` enforces the commit
API). Passes are ordered but communicate only through the ledger — there is no
hidden cross-pass state to reason about.

---

## 3. End-to-end data flow

### Single sample (`assign.run`)

```
Mascope sample_id
   │  io_mascope.fetch_peaks            pull peaks (+ Mascope's own matches), cache as parquet
   ▼
reagent + prescan                       profiles.resolve / reagents.label / isotopes.prescan
   │                                     detect adducts, label reagent-ion clusters, build grid constraints
   ▼
candidate generation                    chemistry.py grid (integer-DBE / Senior / O-cap) per peak
   │
   ▼
SCORING  ── match_compounds ──►         io_mascope: parallel call, flatten_match_tree → per-isotopologue table
   │  (Mascope is the oracle)
   ▼
arbitration + multi-pass commit         passes.py director (passes 1–6 + sweeps; §4)
   │                                     complexity-penalised, isotopologue-gated; commit M0 owners
   ▼
cleanup → degeneracy → tiers            cleanup.py · degeneracy.py · tiers.py
   │                                     recover/relabel, measure mass-degeneracy, assign Identified/Candidate
   ▼
report                                  report.py / pdf — _ledger.csv, _assignments.xlsx, _summary.md, _gka.html
```

### Whole batch (`pipeline.run_batch`)

A batch is many samples over time. `match_compounds` scores against one real
server sample, so assigning one averaged spectrum would miss analytes that appear
only briefly. The batch path assigns a **subset of samples** and merges by m/z.
Two selection strategies (both feed the identical assign → merge → report chain):

- **`representative`** (default) — `sampling.select_representative_samples`: 5
  evenly time-spaced samples + the max-TIC sample.
- **`brightest`** (`--select brightest`) — `sampling.select_brightest_coverage_samples`:
  bin *all* batch peaks by m/z and assign each significant bin's *brightest*
  sample (greedy set-cover). Coverage tracks analyte signal, not a fixed time grid
  — on reagent-CIMS data the max-TIC pick is dominated by the reagent ion and is
  the brightest sample for a small fraction of analyte peaks. A coverage play, not
  a speed play.

```
batch + reagent
   │  pipeline.run_batch
   ▼
sampling.select_*_samples                pick the sample subset (representative | brightest)
   │
   ▼
assign each rep (assign.run)             per-file ledgers kept
   │
   ▼
assign_batch.run  → merge                offset-aware align + jitter table; positive amine gate at merge
   │  (merged_ledger.csv)
   ▼
generate_report (offline, no network):
   ├─ clustering.cluster_batch           correlation clusters of the batch time-series → A4 figure panels
   ├─ analyte_viz.van_krevelen_batch     Van Krevelen (every assigned peak, backbone-colored)
   └─ pdf_report.build                   iterable A4 PDF (cover…methods), + optional compressed companion
   │
   ▼
provenance.record_run                    run_manifest.json + append cross-run index.jsonl
```

The orchestration spine is `pipeline.py`. A `RunContext` (created by
`make_run_context`) carries the run's identity — `out_dir`, `run_id`, the single
`when` timestamp, profile, dataset — into every stage as a plain argument, instead
of threading environment variables or re-deriving `now()` per stage.

---

## 4. The assignment pass sequence

The director lives in `passes.py`. Passes run in order; each only adds
commitments the previous ones justify. (Condensed; the authoritative table is in
[`SKILL.md`](../SKILL.md#the-pipeline).)

| stage        | module            | what it commits                                                                 |
| ------------ | ----------------- | ------------------------------------------------------------------------------- |
| **Pre**      | `reagents`/`isotopes` | detect reagent adducts; prescan isotope fingerprint; label reagent-ion clusters (Brₙ, BrO/BrO₂/BrO₃, ⁷⁹/⁸¹Br) so they are never candidates |
| **0**        | `passes`          | **known species** (committed + locked, first): atmospheric acids, nitroaromatics, PFCAs, ³⁷Cl-confirmed chlorinated paraffins, silanediols, +mode organophosphates — families the generic grid would miss |
| **1**        | `passes`          | lock the high-confidence **CHO/CHON backbone**: enumerate → score → arbitrate → commit M0 owners + isotopologue children |
| **2**        | `passes` (`series_gka`) | **iterative GKA series** expansion from locked anchors (CH₂/O/H₂O/CO/CO₂/…) |
| **3**        | `passes` (`series_detect`) | **automatic series detection** ("rotating plot") opens contaminant families on decoy-controlled evidence |
| **4**        | `residual`        | **residual explainer**: ~1.998-Da isotope doublets, deep 2-step series, ppm-disciplined |
| **5**        | `passes`          | **known-neutral completion**: cross-channel partners + series gaps (no new formula space) |
| **6**        | `ladders`         | **anchored ladder gap-fill**: walk +O/+CH₂/+CO₂/−H₂O diagonals out from Identified anchors (Candidate tier) |
| iso-env      | `isotopes`        | claim each committed peak's full predicted M+2/M+4 envelope; displaces weak M0s that are really a parent's satellite |
| siloxane     | `siloxane`        | dedicated PDMS/siloxane ladder on spacing + ²⁹Si/³⁰Si envelope (CHON monsters out-score the true Si formula otherwise) |
| cleanup      | `cleanup`         | isotope-confirmed recovery, bromide-cluster labelling, ringing/sidelobe artifact flagging |
| degeneracy   | `degeneracy`      | honest cross-family mass-degeneracy density; an uncorroborated mass-degenerate commit is capped at Candidate |
| tiers        | `tiers`           | final **Identified / Candidate** verdict (margin, density, mass-error gate, degeneracy-aware) |

Also interleaved: **composite detection** (`cleanup`/`degeneracy`) flags an M0
whose intensity exceeds its M+1-implied owner — **halide-CIMS only, a no-op in
positive urea mode**. CLI toggles: `--no-pass2/3/4`; `--no-pass5` disables **both**
Pass 5 and the Pass-6 ladder gap-fill.

---

## 5. Chemistry gates (structural invariants)

Enforced in `chemistry.py` and regression-tested; they are the reason a "mass
fit" alone never wins:

- **Integer DBE on the neutral + Senior's rule.** Half-integer DBE is an ion-only
  artifact (deprotonation), so organic nitrates pass as neutrals.
- **Halogens count as hydrogens** in DBE; the Van Krevelen ratio uses
  `(H+F+Cl+Br+I)/(C+Si)` (Si is a carbon-equivalent backbone atom).
- **Structural oxygen cap** `O ≤ 2·(C+N+S+P)+4` — valence, not a Van Krevelen
  prior. Kills `C3H5ClO17`-type mass-fits while real HOMs pass.
- **Reagent-halogen policy.** A reagent halogen's *ion* isotope can't prove it
  sits in the *neutral* (covalent `X(Br)[M-H]⁻` ≡ `Y·HBr·Br⁻`). The complexity
  prior on the reagent element is never waived by isotope confirmation.
- **Isotopologue-gated heteroatoms.** A neutral S/Cl/Br needs its Mascope-confirmed
  ³⁴S/³⁷Cl/⁸¹Br, or its complexity skepticism stands.

---

## 6. Reagent profiles

`profiles.py` defines one `ReagentProfile` per reagent system (bromide `Br⁻`,
urea/uronium `Ur⁺`, nitrate `NO3` / `NO3_15N`), bundling polarity, adducts,
normaliser ion, and the chemistry context. `resolve('auto', peaks)` detects the
reagent from the spectrum; passing `--reagent` forces it (a sparse positive sample
otherwise mis-detects as negative). New reagents are added from a JSON/TOML file
via `register()` / `load_config()` (`--reagent-config`) **without forking the
package**. The reagent is the single switch that makes the same pipeline
negative- or positive-mode.

---

## 7. Reproducibility & provenance

**The scientific content of a run is a pure function of its inputs.** Every
figure's pixels, the ledger, and every cluster/Van-Krevelen table is byte-identical
whenever you re-run the same data — *regardless of when you run it*. The **only**
thing the run timestamp changes is the PDF cover's "generated" line, the **Report
ID**, the run-folder name, and the `run_manifest.json` provenance. So a re-run
reproduces the analysis exactly; only *which run it was* is stamped on top.

How that's enforced:

- **Fixed content epoch.** `pipeline.stamp_source_date_epoch()` pins
  `SOURCE_DATE_EPOCH` to a constant (`CONTENT_EPOCH`, 1980-01-01Z), **not** the run
  time. matplotlib reads it for PNG/PDF metadata and the xlsx writer reads it via
  `cluster._resolve_when`, so the embedded timestamps are constant and the bytes
  depend only on the data. No run-time value ever leaks into a figure or a data
  table (the assignment xlsx carries no "generated" cell, by design).
- **Run time as text only.** `RunContext.when` (UTC) names the run folder
  (`run_id = slug + YYYY-MM-DDTHHMMSSZ`) and is rendered as *visible text* on the
  PDF cover (the "generated" line + Report ID) — never through `SOURCE_DATE_EPOCH`.
  A report PDF therefore differs run-to-run only because that visible text differs.
- **`run_manifest.json`** (`provenance.py`) pins each run to its exact code
  (package + per-module version + content hash + git commit), input-data hash
  (`ts_sha1`), resolved config, and output hash (`merged_ledger_sha1`).
- **`index.jsonl`** at the `--out-dir` base (not inside the run folder) is the
  cross-run registry — one compact row per run, loadable with
  `pandas.read_json(lines=True)` to find or diff runs.

The one genuine source of cross-run variation is the live `match_compounds` call
(server-side), not Peaky. `test_determinism.py` locks the contract: two runs at
different times over the same inputs → identical figure/xlsx/csv bytes, with the
PDF differing only by its visible cover timestamp.

---

## 8. Outputs

A batch run writes one **versioned folder** per run (`<slug>_<timestamp>Z/`) so a
re-run never overwrites a previous one. The layout is the single source of truth
in `paths.py` (`RunPaths`), shared by the writers and the report reader so their
filename contract can't drift:

```
<run>/
  merged_ledger.csv     the result (provenance anchor)        — ROOT
  run_manifest.json     reproducibility manifest              — ROOT
  batch_summary.json    counts / per-file offsets             — ROOT
  per_file/             per-sample ledgers
  figures/              all .png (cluster panels, GKA, Van Krevelen)
  tables/               all .csv / .xlsx (cluster tables, jitter, channel QC)
  report/               the PDF report (+ compressed companion)
  data/                 a fetched time-series, kept with the run (only when no
                        on-disk source exists — a parquet passed by path is
                        referenced, never copied)
```

(`index.jsonl` — the cross-run registry — lives at the `--out-dir` base, *not*
inside the run folder.) The canonical artifacts are `merged_ledger.csv` (the
result) and `run_manifest.json` (the provenance anchor); both stay at the run root
because several modules + the registry read them. **[OUTPUTS.md](OUTPUTS.md) is the
full per-artifact reference** (batch + single-sample, one line each).

---

## 9. Module map

Grouped by responsibility (full table in [`SKILL.md`](../SKILL.md#module-map-peaky)):

**I/O & scoring** — `io_mascope.py` (the only Mascope I/O: peaks, cheminfo,
parallel `match_compounds` + per-isotopologue parser, offset estimation; legacy
workspace-server fallback).

**Chemistry & candidates** — `chemistry.py` (masses, formula algebra, grid,
complexity penalty), `isotopes.py` (prescan → grid constraints, envelope
predictor), `reagents.py` / `profiles.py` (reagent library + per-reagent config).

**Assignment** — `ledger.py` (state + invariants + commit API), `passes.py`
(arbitration + pass director + calibration), `series_gka.py` / `series_detect.py`
/ `ladders.py` (series math, detection, ladder gap-fill), `residual.py` (pass 4),
`siloxane.py` (PDMS ladder), `cleanup.py` (residual cleanup), `degeneracy.py`
(mass-degeneracy), `tiers.py` (Identified/Candidate verdict), `plausibility.py`
(QC), `assign.py` (orchestrator + `PassConfig`).

**Batch** — `sampling.py` (sample selection: representative subset OR
brightest-coverage), `assign_batch.py`
(assign reps + offset-aware merge), `timeseries.py` (time-resolved disposition),
`clustering.py` + `cluster.py` (correlation-cluster figures), `composition.py`
(signal-weighted composition accounting).

**Reporting & orchestration** — `analyte_viz.py` (Van Krevelen + time-series),
`gka_figure.py` (static GKA page), `gka_widget.py` (interactive HTML widget),
`report.py` (Excel/markdown), `pdf_report.py` (iterable A4 PDF),
`pipeline.py` (the spine + `RunContext`), `provenance.py` (manifest + registry),
`cli.py` (the `peaky` console entry point).

---

## 10. Testing

The offline suite (`tests/`, 850+ assertions, no network) is the contract. Pure
functions (`flatten_match_tree`, the chemistry gates, the merge, the cluster
math) are unit-tested against captured fixtures; determinism is locked by
`test_determinism.py`; the no-network install guarantee by `test_smoke.py`. CI
runs the suite on Python 3.11–3.13 with **no credentials**. Every change ships
with a test, and the suite must stay green — that is what keeps the pipeline
trustworthy as it grows.

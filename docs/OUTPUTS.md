# Peaky — Outputs reference

Where Peaky writes things, and what each artifact is. The batch run-folder layout
is the single source of truth in [`peaky/paths.py`](../peaky/paths.py) (`RunPaths`),
shared by the writers and the report reader so the filenames can't drift.

---

## Batch run — one versioned folder per run

`peaky batch` creates **one timestamped folder** under `--out-dir`
(default `~/peaky-output`), so a re-run never overwrites a previous one:

```
<batch-slug>_<YYYY-MM-DDTHHMMSSZ>/     ← folder name == Report ID (UTC stamp)
```

### Run root (flat — read by several modules + the cross-run registry)

| Artifact | What it is / what it's for |
|---|---|
| `merged_ledger.csv` | **The result.** Every merged peak (one row each): role, neutral formula + adduct, scores, ppm, confidence, tier, provenance. The provenance anchor. |
| `run_manifest.json` | **Reproducibility manifest.** Pins the run to its exact code (package + per-module version + content hash + git commit), input-data hash (`ts_sha1`), resolved config (incl. `select` / `coverage_target`), and output hash (`merged_ledger_sha1`). |
| `batch_summary.json` | Run counts + per-file calibration offsets (M0/tier counts, n_files, offsets) and the selection strategy used (`select`, `coverage_target`). |
| `per_file/<sid>_ledger.csv` | The full single-sample ledger for **each** assigned sample, kept for audit / re-merge. |
| `index.jsonl` | **At the `--out-dir` base, not inside the run folder.** Cross-run registry — one compact row per run, loadable with `pandas.read_json(lines=True)` to find or diff runs. |

### `figures/` — all `.png`

| Artifact | What it is |
|---|---|
| `van_krevelen_<tag>.png` | Van Krevelen of the assigned analytes (Si excluded — the clean atmospheric view). |
| `van_krevelen_full_<tag>.png` | Full Van Krevelen — every assigned peak by CHO/CHON/CHOS backbone (Si/F/halogen folded in). |
| `clusters_changing_<tag>_p*.png` | Correlation-cluster panels (A4 portrait, paginated) of the dynamic, co-varying analyte families. |
| `clusters_flat_<tag>_p1.png` | The uncorrelated/flat remainder + Si contamination, bunched into one overview. |
| `clusters_changers_<tag>_p*.png` | Big standalone changers — single channels that move ≥~5–10× on their own. |
| `clusters_unassigned_<tag>_p*.png` | The same clustering applied to the **unexplained** residual. |
| `gka_<tag>.png` | The static GKA findings page embedded in the report. |

### `tables/` — all `.csv` / `.xlsx`

| Artifact | What it is |
|---|---|
| `selected_samples.csv` | Which samples were assigned + why: a `role` (`time-grid` / `max-TIC`, or `coverage-winner` for `--select brightest`) and `bins_won` (significant m/z bins the sample is brightest for). |
| `jitter.csv` | Per-(cluster, file) mass-jitter table — raw vs calibration-adjusted ppm spread of each merged assignment. |
| `van_krevelen_full_<tag>.csv` | The full-VK data behind the figure (one row per assigned neutral). |
| `clusters_changing_<tag>.csv` / `.xlsx` | Cluster membership; the XLSX has one tab per cluster (formula / channel / m/z / match_score / tier). |
| `clusters_flat_<tag>.csv`, `clusters_changers_<tag>.csv`, `clusters_unassigned_<tag>.csv` | Membership for the flat / changers / unassigned figures. |
| `channel_agreement_<tag>.csv` | QC: how often a multi-channel neutral's ion channels agree in time. |

### `report/` — the PDF

| Artifact | What it is |
|---|---|
| `report_<run-id>.pdf` | **The standard iterable A4 report** (cover · findings · coverage · composition · scrutiny · GKA · families · changers · clusters · methods). The cover shows the Report ID + a date+time "generated" line. |
| `report_<run-id>_compressed.pdf` | Optional size-reduced companion for emailing (needs `pip install mascope-peaky[compress]`). The full report is left byte-for-byte untouched. |

### `data/` — bulky inputs kept with the run

| Artifact | What it is |
|---|---|
| `<tag>_ts.parquet` | The full-batch per-sample peak time series — written here **only** when fetched live (no on-disk source). A parquet passed by `--ts` is *referenced*, never copied. |

---

## Single-sample run — `peaky assign`

Writes into `--output-dir` with the prefix `<sample-id>_<YYYYMMDD-HHMM>`:

| Artifact | What it is |
|---|---|
| `<prefix>_ledger.csv` | Every peak: role (`M0` / `iso_child` / `reagent` / `artifact` / `unexplained`), formula, adduct, all scores (incl. arbitration `eff_score`/`eff_margin`/`tied`), ppm, confidence, `tier` + `tier_reason`, candidate/degeneracy density, provenance, commentary, alternatives, isotopologues. |
| `<prefix>_assignments.xlsx` | The styled multi-sheet workbook: Summary · Read-me legend · **Identified** · **Candidates** · Unassigned (evidence-characterized) · By class · Unique formulas · Isotopologues · Peak ownership · Target list · Reagent ions. Frozen headers, autofilters, tier/confidence color chips. |
| `<prefix>_summary.md` | Narrative + top assignments + coverage. |
| `<prefix>_manifest.json` | Module versions, prescan fingerprint, series-evidence table, per-pass timing. |
| `<prefix>_gka.html` | Interactive rotating-GKA widget (self-contained, no server). |
| `<prefix>_gka_unexplained.html` | The same widget over the **unexplained residual only** — the place to hunt for missed homologous structure. |
| `checkpoints/` | Per-pass ledger checkpoints (an audit trail of what each pass committed). |

---

## Reproducibility note

Every figure, table, and ledger above is a **pure function of the input data** —
byte-identical whenever you re-run the same data. The **only** thing the run
timestamp changes is the PDF cover's "generated" line + the Report ID + the
run-folder name + `run_manifest.json`. See
[ARCHITECTURE.md §7](ARCHITECTURE.md#7-reproducibility--provenance).

> `peaky report --run-dir <folder> ...` regenerates the `figures/` + `report/`
> artifacts of an existing run **offline** (no assignment, no network) from the
> ledgers already on disk + the TS parquet.

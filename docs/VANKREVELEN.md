# Peaky — Van Krevelen & analyte traces

This document explains **how a committed ledger becomes the Van Krevelen scatter
and the per-ion/per-compound time traces** — the analyte definition, the H/C–O/C
convention, the RAW-intensity choice, and the "changing" threshold shared by both
instruments. It is a module deep-dive companion to
[`ARCHITECTURE.md`](ARCHITECTURE.md) (the whole pipeline),
[`CLUSTERING.md`](CLUSTERING.md) (which clusters the `ion_traces` built here), and
[`TIMESERIES.md`](TIMESERIES.md) (whose `build_matrix` these traces sit on).

**Code:** `peaky/reporting/analyte_viz.py`. Pure data functions (no plotting deps)
+ lazy-matplotlib renderers + a `widget_payload()` for the interactive chat widget.
Negative Br-CIMS and positive urea-CIMS figures are computed **identically** —
same analyte definition, same convention, same threshold.

> Keep this in sync with the code. Every threshold below is a named constant or a
> literal in `analyte_viz.py`; if you change one there, change it here.

---

## 1. What this stage does

Two figures from one ledger (+ an optional batch time series): a **Van Krevelen**
(every analyte placed by H/C vs O/C of its *neutral*) and **time traces** (each
compound's RAW intensity over the run). The same module also runs a **channel
agreement** QC (do a neutral's ion channels track in time?) and emits the
`widget_payload` arrays the interactive plot consumes.

```
committed ledger (M0 rows)                 batch time series (optional)
   │ analyte_table: organic neutral M0          │ build_matrix → samples × m/z bins
   │   add H/C=nH/nC, O/C=nO/nC, klass          │
   │   one row per NEUTRAL (best tier)          │
   ▼                                            ▼
 Van Krevelen (x=O/C, y=H/C)        time_traces (sum adducts) · ion_traces (per ion)
   │  organic-only (Si excluded)      │  RAW intensity (mode='reagent'/'tic' opt-in)
   │  full (every peak, by backbone)  │  attach_dynamics: median_cps, cv, changing
   ▼                                  ▼
 van_krevelen_{tag}.png  +  _full.png + .csv      widget_payload / channel_agreement
```

---

## 2. Inputs

- A **committed ledger** (the merged batch ledger or a single-file ledger; only
  `role`, `neutral_formula`, `adduct`, `tier`, `mz` are read).
- An optional **batch time-series peak table** (`sample_item_id`, `mz`, `height`,
  `datetime_utc`) for the dynamics + traces.

---

## 3. The transformation, stage by stage

1. **Analyte table** (`analyte_table`). Keep `role == M0` with a non-empty
   `neutral_formula`; **exclude contamination** (`Si`, the
   siloxane/silanediol set) when `exclude_contaminant`. Parse counts → `nC nH nO
   nN nS`; require `nC > 0`; compute `hc = nH/nC`, `oc = nO/nC`; class
   `CHOS`(S) / `CHON`(N) / `CHO`. **One row per neutral** — a compound seen in two
   channels (`[M+Br]⁻` and `[M-H]⁻`, or `[M+H]⁺` and `[M+urea+H]⁺`) is a single
   analyte; the best tier wins. `channels` records the adduct(s) it was seen in.

2. **The Van Krevelen convention.** Always on the **neutral**: `x = O/C`,
   `y = H/C`. Halogens are not in these neutrals (the reagent halogen lives in the
   adduct), so the axes are clean. Polarity/reagent-agnostic — it reads the neutral
   formula, never the adduct.

3. **Compound traces** (`time_traces`). Per neutral, the **RAW summed** intensity
   of its bins across the given adducts (`_bins_for_formula`, ≤ 8 ppm). Default
   `bin_minutes=None` = **native per-sample resolution** (one point per sample at
   its real hour, no re-gridding); an int time-bins (legacy, median per bin).

4. **Per-ion traces** (`ion_traces`). One trace **per ion** (no summing across
   adducts), each from the single TS bin matching that ion's m/z (≤ 8 ppm). Keeping
   channels separate lets divergent channels of one neutral cluster apart — this is
   the input [`CLUSTERING.md`](CLUSTERING.md) consumes.

5. **Dynamics** (`attach_dynamics`). From the traces, add `median_cps`,
   `cv = std/mean`, and `changing = cv ≥ CHANGING_CV (0.30)` — the **same
   threshold both instruments**, so the figures are comparable.

6. **Channel agreement QC** (`channel_agreement`). For every neutral with ≥ 2
   testable channels (≥ `min_points` 8 finite points **and** median ≥ `floor` 150
   cps), correlate each channel pair on `log10`; report `worst_r` (least-agreeing
   pair), `top2_r` (the two brightest channels — the pair that dominates the sum),
   and a `verdict`: **agree** `r ≥ 0.7` / **marginal** `≥ 0.4` / **disagree**
   `< 0.4` on `worst_r`.

7. **Full Van Krevelen** (`render_van_krevelen_full`). EVERY assigned peak by
   **CHO/CHON/CHOS backbone** (`backbone_class` folds Si/F/halogen *into* their
   backbone class, not a separate class); changing analytes solid with a white
   edge, flat ones dimmed, size = log intensity. (The organic-only
   `render_van_krevelen` excludes Si and colours the changing CHO/CHON/CHOS.)

8. **Widget payload** (`widget_payload`). `vk` rows `[oc, hc, hasN, changing, logI,
   formula, channels, tier]`; `ts` = the top `top_ts` (28) changing series by
   `median_cps`.

---

## 4. Constants reference

All in `peaky/reporting/analyte_viz.py`.

| constant | value | role |
| --- | --- | --- |
| `CHANGING_CV` | 0.30 | `cv` at/above which a compound is "changing" (both instruments) |
| `DEFAULT_BIN_MIN` | 30 | legacy time-bin width (minutes); native cadence is the default |
| `_CONTAM_ELEMENTS` | `(Si,)` | siloxane/silanediol = instrument contamination (excluded from organic VK) |
| `_bins_for_formula` / `_bin_for_mz` tol | 8 ppm | m/z window matching a formula/ion to a TS bin |
| `channel_agreement` `floor` / `min_points` | 150 cps / 8 | min brightness + finite points for a testable channel |
| channel verdict | agree ≥ 0.7 · marginal ≥ 0.4 · disagree < 0.4 | on `worst_r` (log10 Pearson) |
| `render_van_krevelen` limits | x∈[0,1.3], y∈[0.3,3.0] | organic-only VK axes |
| `render_van_krevelen_full` `xmax`/`ymax`/`dpi` | 1.4 / 4.2 / 150 | full VK axes + resolution |
| `widget_payload` / `render_timeseries` top | 28 | brightest changing series rendered |

---

## 5. Metrics, defined

- **H/C, O/C** — `nH/nC`, `nO/nC` of the *neutral*; the Van Krevelen coordinates.
- **`median_cps`** — median of a compound's (RAW) detected trace; the size/brightness
  axis.
- **`cv`** — `std/mean` of the (RAW) trace; `changing` iff ≥ 0.30.
- **`worst_r` / `top2_r`** — least-agreeing / brightest-pair log10 Pearson r among a
  neutral's channels; the channel-agreement QC.
- **backbone class vs full class** — `backbone_class` folds Si/F/halogen into
  CHO/CHON/CHOS; `full_class` splits them out (siloxane / F-containing /
  halogenated) for the full-VK legend.

---

## 6. Outputs

| artifact | content |
| --- | --- |
| `figures/van_krevelen_<tag>.png` | organic-only VK (Si excluded), changing analytes coloured |
| `figures/van_krevelen_full_<tag>.png` | every assigned peak by CHO/CHON/CHOS backbone |
| `tables/van_krevelen_full_<tag>.csv` | one row per assigned neutral: formula, adduct, tier, oc, hc, fclass, median_cps, cv, changing |
| `time_traces` / `ion_traces` | `(x-hours, DataFrame)` per-compound / per-ion traces |
| `channel_agreement` | one row per multi-channel neutral: n_channels, worst_r, top2_r, channels, verdict |
| `widget_payload` | `{vk:[…], ts:{grid, series}}` for the interactive chat widget |

---

## 7. Properties, invariants & gotchas

- **ANALYTE = committed organic M0**, excluding Si contamination and the
  reagent/artifact roles; defined on the **neutral**, so the same code gives
  comparable negative- and positive-mode figures.
- **RAW intensity, not normalised, by default.** A positive urea-CIMS spectrum
  often excludes the reagent ions from its mass range and its TIC is
  analyte-dominated, so reagent/TIC normalisation introduces *closure artifacts*.
  `mode='reagent'/'tic'` exists for when a real normaliser is present.
- **Native cadence by default** (`bin_minutes=None`): one point per sample at its
  real hour, no re-gridding (which aliases into spurious empty bins). The legacy
  int path medians into a uniform grid.
- **One row per neutral.** Multi-channel detections collapse to a single analyte
  with the best tier — the VK never double-plots the same compound.
- **`ion_traces` is the clustering input**: it deliberately does **not** sum across
  adducts, so anti-phase channels of one neutral can split into different families.
- **`backbone_class` ≠ `full_class`.** The full VK colours by backbone (Si folded
  in); the legend's `full_class` is a separate split for showing contamination.
- **Channel agreement is QC only** — it never changes an assignment; it reports
  whether summing channels (`time_traces`) is safe for a given neutral.

---

## 8. Code map

| function | role |
| --- | --- |
| `analyte_table` | committed organic M0 → one row per neutral with H/C, O/C, class, channels |
| `full_class` / `backbone_class` | composition class for the full VK legend / backbone colouring |
| `time_traces` | per-compound RAW trace (summed across adducts) |
| `ion_traces` | per-ion trace (channels kept separate) — the clustering input |
| `attach_dynamics` | add `median_cps` / `cv` / `changing` from the traces |
| `channel_agreement` / `_logcorr` | QC: do a neutral's channels track in time? |
| `render_van_krevelen` / `render_van_krevelen_full` | organic-only / full VK PNGs |
| `render_timeseries` | the changing-analyte trace plot |
| `widget_payload` | arrays for the interactive chat widget |
| `van_krevelen_batch` | the batch entry point: both figures + the full-VK CSV |

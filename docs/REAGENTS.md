# Peaky — Reagents & profiles (the single mode switch)

This document explains **how one reagent choice configures the whole pipeline** —
the analyte adduct channels, the grid element box, the time-series normaliser, and
the labeled reagent-ion clusters — and **how reagent-cluster m/z are enumerated
and matched**. It is a module deep-dive companion to
[`ARCHITECTURE.md`](ARCHITECTURE.md) (the whole pipeline, §6 *Reagent profiles*),
[`CHEMISTRY.md`](CHEMISTRY.md) (the adduct shifts + grid the profile drives), and
[`ASSIGNMENT.md`](ASSIGNMENT.md) (which consumes the reagent labels).

**Code:** `peaky/chem/profiles.py` (`ReagentProfile` + the registry + `resolve`)
and `peaky/chem/reagents.py` (the cluster-ion library + `label_reagents`).

> Keep this in sync with the code. Every threshold below is a named constant or a
> literal in `reagents.py` / `profiles.py`; if you change one there, change it here.

---

## 1. What this stage does

In chemical-ionization MS the reagent anion/cation forms **bright cluster ions
that are not sample chemistry** (bare Rₙ, R·water, R·HBr, protonated-urea
clusters). Left alone they dominate the unexplained residual by signal. This layer
does two things:

1. **Profile selection** — pick (or auto-detect) one `ReagentProfile`, which is
   *everything the pipeline needs to treat the reagent correctly*: polarity, the
   analyte adduct channels, the grid ranges, the normaliser (`reagent`/`tic`), the
   reagent-ion regex, the assignment context, and any isotopic purity. **The
   reagent is the single switch that flips the same pipeline negative- or
   positive-mode.**
2. **Reagent labeling** — enumerate the cluster m/z (with halogen isotopologue
   combinations) and mark matching ledger peaks `role='reagent'`, recording the
   *known* cluster ion formula as the assignment.

```
peaks  ──► resolve('auto', peaks)              name/alias ──► resolve(name)
              │ detect_adduct in sample? → profile      │
              │ else polarity → profile                 ▼
              ▼                              ReagentProfile {polarity, adducts,
        ReagentProfile  ────────────────────► ranges, normaliser, reagent_ion_re,
              │                                 detect_adduct, context, purity}
              ▼  reagent_for_adducts → library key ("Br"/"I"/"Cl"/"urea")
        build_library:
          halide  → bare Rₙ⁻ (odd closed-shell, even radical) · Rₙ·neutralₖ · ROₙ⁻
                    (all halogen isotopologue combos)
          urea    → [Rₙ+H]⁺ protonated series (n=1..6)
              ▼  label_reagents (±ppm 15, nearest label, known ion_formula)
        ledger peaks marked role='reagent'
```

---

## 2. Inputs

- A **reagent selector**: a name/alias (`"Br"`, `"uronium"`, `"15no3"`), or
  `"auto"` + a loaded peak table.
- Optional **`--reagent-config`** JSON/TOML registering extra/override profiles.
- For labeling: the **ledger** (peaks with `role`/`mz`/`peak_id`).

---

## 3. The transformation, stage by stage

1. **Resolve a profile** (`profiles.resolve`). A name/alias hits `_BY_ALIAS`
   directly. `"auto"` detects from the sample: the **diagnostic `detect_adduct`**
   among the server's own adduct mechanisms (`io_mascope.detect_adducts`) wins
   first; failing that, **polarity** (`_detect_polarity`) picks the first
   profile of that sign. A `config` path is loaded (registered) before resolving.

2. **The profile configures everything else.** `ReagentProfile` (frozen) carries:
   `polarity`, `adducts` (analyte channels), `ranges` (the grid box string fed to
   [`CHEMISTRY.md`](CHEMISTRY.md)), `normaliser` (`reagent` or `tic`, for the
   TS/correlation layer), `reagent_ion_re` (regex on `ion_formula` picking reagent
   ions), `detect_adduct`, `context` (the assign-time mode + VK priors + caps), and
   `purity` (a labelled reagent's isotopic purity, threaded to
   `predict_isotopes`), and — for labelled reagents — `label_isotope` /
   `label_max` (the caret heavy isotope a *covalent product* can carry and its max
   count; drives the heavy-isotope rescue in `labeled.py`). Built-ins: **`BR`**
   (Br⁻, neg, normalise on reagent), **`UR`** (urea/uronium, pos, normalise on
   TIC), **`NO3`** (nitrate, neg, reagent), **`NO3_15N`** (¹⁵N nitrate, neg, TIC,
   `purity 0.98`, `label_isotope='^N'`, `label_max=2`).

   > **¹⁵N-nitrate ¹⁴NO₃-cluster hazard.** In a NOx-oxidation run the chamber holds
   > abundant *unlabelled* ¹⁴NO₃⁻, so a highly-oxygenated analyte X forms
   > `[X+¹⁴NO₃]⁻` — the **exact isobar** of the deprotonated covalent organonitrate
   > `[Y−H]⁻` (Y = X + HNO₃). `[M+NO3]-` (¹⁴N) is therefore **deliberately kept OFF
   > `NO3_15N.adducts`** (adding it would let the scorer arbitrate the isobar
   > arbitrarily); instead the post-tier `relabel_nitrate_clusters` pass re-reads a
   > covalent organonitrate as the cluster only when the parent X is independently
   > corroborated (ASSIGNMENT_DETAIL §3.9). The ¹⁵N cluster channel is the ordinary
   > `[M+^NO3]-`.

3. **Pick the cluster-library key** (`reagent_for_adducts`). From the analyte
   adducts: `CH4N2O` → `"urea"`; `Br` → `"Br"`; `I` → `"I"`; `Cl` → `"Cl"`. This
   is the **cluster-library key, not** the arbitration `reagent_element` (a
   molecular reagent puts no halogen in the neutral, so `assign.run` sets
   `reagent_element` only for the halogen keys).

4. **Build the cluster library** (`build_library`). For a **halide** (Br/Cl/I):
   - **bare Rₙ⁻**, `n = 1..max_n (4)`, every isotopologue combination
     (`combinations_with_replacement`); **odd n** are closed-shell anions (R⁻,
     R₃⁻), **even n** are radical anions (R₂⁻·, R₄⁻·); each anion adds `M_E`.
   - **Rₙ⁻·(neutral)ₖ**, `k = 1..max_neutral (1)`, over `_CLUSTER_NEUTRALS`
     (`H2O`, `HBr`, `HF` only).
   - **reagent-oxide anions** RO⁻/RO₂⁻/RO₃⁻, *both* halogen isotopologues.

   For a **molecular reagent** (urea): the protonated series **[Rₙ+H]⁺**,
   `n = 1..max_n (6)` (`_build_positive_library`) — a cation (lose an electron),
   repeat unit `CH4N2O`, no isotope branching.

5. **Label matching peaks** (`label_reagents`, `ppm 15`, `only_unexplained=True`).
   `bisect` the sorted library masses within `±(mz·ppm·1e-6)`, take the **nearest**
   label, and `L.mark_reagent` records the **known cluster ion formula** as the
   assignment (a reagent cluster has a known formula → it is *assigned*, just a
   different class — never left blank/red in the residual). Returns the count.

6. **The positive cluster adducts carry N/O into every analyte ion** (uronium/
   ammonium mode). The analyte channels of a positive N-reagent are not clean
   protonation: the reagent contributes atoms to the observed ion. `[M+NH4]⁺` adds
   `NH3` over `[M+H]⁺`, and the uronium/urea cluster `[M+(CH4N2O)H]⁺` adds a whole
   `CH4N2O`. Two downstream consequences fall out of that extra N/O and are handled
   in assignment, not here:
   - **CHO-via-N-cluster is isobaric with protonated CHON.** A CHO neutral seen as
     `[M+NH4]⁺` or `[M+(CH4N2O)H]⁺` yields the *same* ion formula as a CHON neutral
     seen as `[M+H]⁺` (e.g. `C12H14O4[M+NH4]⁺` and `C12H17NO4[M+H]⁺` are both
     `C12H18NO4⁺`). Mass and isotopes cannot separate them — the tier layer demotes
     such a commit to Candidate unless an extra-spectral discriminator corroborates
     it. See the reagent-N isobar gate in
     [`ASSIGNMENT_DETAIL.md`](ASSIGNMENT_DETAIL.md) (tiers).
   - **A pure hydrocarbon "seen via an N-cluster" is re-read as an N-heterocycle.**
     A hydrocarbon has no site to bind `NH4⁺`/uronium and would show `[M+H]⁺`, so an
     `[M+NH4]⁺` / `[M+(CH4N2O)H]⁺` reading of a bare CₓHᵧ is implausible; assignment
     re-reads it as `[M+H]⁺` of the N-heterocycle `M' = M + (cluster − H)` (see
     [`ASSIGNMENT_DETAIL.md`](ASSIGNMENT_DETAIL.md), reagent-N re-read).

---

## 4. Constants reference

`reagents.py` (library) + `profiles.py` (profiles).

| constant | value | role |
| --- | --- | --- |
| `_HALOGEN_ISO` | Br: ⁷⁹/⁸¹ · Cl: ³⁵/³⁷ · I: ¹²⁷ | reagent-halogen isotope masses/labels |
| `_CLUSTER_NEUTRALS` | `H2O`, `HBr`, `HF` | neutrals that cluster on a halide core (organic acids deliberately removed) |
| `_POSITIVE_REAGENTS` | `{urea: CH4N2O}` | molecular positive reagents (protonated series) |
| reagent-N cluster `[M+NH4]⁺` | `{N:1, H:3}` over `[M+H]⁺` | ammonium cluster adduct N/O added to the ion (isobar gate / re-read) |
| reagent-N cluster `[M+(CH4N2O)H]⁺` | `{C:1, H:4, N:2, O:1}` | uronium/urea cluster adduct atoms added to the ion (isobar gate / re-read) |
| `build_library` `max_n` | 4 (halide) / 6 (positive) | largest Rₙ cluster enumerated |
| `build_library` `max_neutral` | 1 | max neutral adducts per halide core |
| `label_reagents` `ppm` | 15.0 | reagent-cluster mass-match window |
| `label_reagents` `only_unexplained` | True | only relabels still-unexplained peaks |
| `ReagentProfile.normaliser` | `reagent` / `tic` | TS/correlation normalisation basis |
| `BR.ranges` | `C0-40 H0-80 N0-3 O0-18 S0-2 Cl0-2 Br0-2` | bromide grid box |
| `UR.ranges` | `C0-40 H0-90 N0-8 O0-15 S0-2` | uronium grid box |
| `NO3.ranges` / `NO3_15N.ranges` | `C0-40 H0-60 N0-3 O0-25 S0-2` | nitrate grid box |
| `NO3_15N.purity` | 0.98 | ~98 % ¹⁵N reagent (→ `predict_isotopes`) |

---

## 5. Metrics, defined

- **cluster ion m/z** — built from `_HALOGEN_ISO` masses ± `M_E` (anion +e,
  cation −e), summed over the isotopologue combination + any clustered neutral.
- **reagent ppm match** — `(mz − ion_mz)/ion_mz · 1e6`; within ±15 ppm of the
  *nearest* library entry the peak is labeled reagent (the ppm is recorded in the
  note).
- **normaliser** — `reagent` divides traces by a reagent-ion signal; `tic` divides
  by total ion current (used when the reagent ions sit below the acquisition
  window, e.g. ¹⁵N-nitrate / uronium).

---

## 6. Outputs

| artifact | content |
| --- | --- |
| `ReagentProfile` | the run's mode config: polarity, adducts, ranges, normaliser, reagent_ion_re, detect_adduct, context, purity |
| `build_library` | `[(label, ion_mz, ion_formula)]` — the enumerated reagent-cluster ions |
| `label_reagents` | count of ledger peaks set to `role='reagent'` (each with a known `ion_formula`) |
| `reagent_for_adducts` | the cluster-library key (`"Br"`/`"I"`/`"Cl"`/`"urea"`/`None`) |

---

## 7. Properties, invariants & gotchas

- **Organic acids are NOT reagent neutrals.** A `[Br+acid]⁻` ion *is* the
  `[acid+Br]⁻` primary `[M+Br]⁻` analyte channel — the labeler used to steal real
  ambient acids (formic acid's 232k-cps line among them) and bury them as reagent.
  `_CLUSTER_NEUTRALS` is now water + the reagent's own HBr + background HF only.
- **Both Rₙ parities are real** reagent ions: odd = closed-shell anion, even =
  radical anion (e.g. the Br₂⁻· di-bromide). All are pure reagent → labeled, not
  left in the residual.
- **A molecular reagent puts no halogen in the neutral.** `reagent_for_adducts`
  returns the *library key*; `reagent_element` (the arbitration complexity element)
  is set only for halogen reagents — never for urea.
- **Charge bookkeeping:** halide clusters *gain* an electron (`+M_E`); the
  protonated positive series *loses* one (`−M_E`).
- **Positive N-reagents make CHO readings isobaric with CHON.** Because `[M+NH4]⁺`
  and `[M+(CH4N2O)H]⁺` fold reagent N/O into the observed ion, a CHO neutral via one
  of these clusters is the *same ion formula* as a protonated CHON neutral — neither
  mass nor isotopes distinguish them. This is a positive-mode-only hazard (the
  halide anions add no N); the tier layer resolves it (isobar gate,
  [`ASSIGNMENT_DETAIL.md`](ASSIGNMENT_DETAIL.md)). It is why `UR.ranges` opens N to
  `N0-8` — the analyte space genuinely spans both readings.
- **A bare hydrocarbon has no N-cluster site.** `NH4⁺`/uronium bind at a polar/basic
  site a pure CₓHᵧ lacks; a hydrocarbon ionises as `[M+H]⁺`. An `[M+NH4]⁺` /
  `[M+(CH4N2O)H]⁺` reading of a hydrocarbon is therefore re-read as the N-heterocycle
  `[M+(cluster−H)+H]⁺` (assignment, not this layer) — *unless* that hydrocarbon also
  has a genuine `[M+H]⁺` row (a real terpene forming `[M+NH4]⁺`), which is left
  alone.
- **`detect_adduct` disambiguates isotopic twins.** `NO3` vs `NO3_15N` differ only
  by their diagnostic adduct (`[M+NO3]⁻` vs `[M+^NO3]⁻`), so auto-detect picks the
  right one; `purity` then flows to the labelled-reagent envelope predictor.
- **A reagent cluster is an assignment, not a blank.** Its formula is known, so it
  is committed with that `ion_formula` (a distinct class), never red in the report.
- **New reagent = a `ReagentProfile`, no fork.** `register` / `from_dict` /
  `load_config` add reagents from JSON/TOML (`--reagent-config`); a top-level list,
  a `{"reagents": [...]}` wrapper, or a `{name: {fields}}` mapping all parse.
- **Auto-detect needs peaks**; with no diagnostic adduct it falls back to polarity,
  and a sparse positive sample can mis-detect as negative — pass `--reagent` to
  force it.

---

## 8. Code map

| function | role |
| --- | --- |
| `profiles.ReagentProfile` | the frozen per-reagent config dataclass |
| `profiles.resolve` | name/alias or `auto` (detect_adduct → polarity) → a profile |
| `profiles.register` / `from_dict` / `load_config` | registry + JSON/TOML reagent loading |
| `profiles._detect_polarity` | infer `+`/`−` from the peak table |
| `reagents.reagent_for_adducts` | analyte adducts → cluster-library key |
| `reagents.build_library` | enumerate the reagent-cluster ions (halide + positive) |
| `reagents._build_positive_library` | the `[Rₙ+H]⁺` protonated-reagent series |
| `reagents.label_reagents` | mass-match + mark ledger peaks `role='reagent'` |

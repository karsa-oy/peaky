# Peaky Assignment Pipeline — Detailed Reference

## 1. Overview, Ledger Lifecycle, and Role/Tier Vocabulary

### 1.1 What Peaky does

Peaky is an AI-native analysis toolbox for the Mascope mass-spectrometry platform. Its assignment pipeline transforms a raw high-resolution MS peak table into a chemical-formula inventory with quantified confidence. The core abstraction is the **ledger**: a pandas DataFrame with **one row per physical peak**, mutated in place by a sequence of assignment passes. Each pass scores candidate neutral formulas against the **Mascope server's `match_compounds`** (an isotope-scored oracle), arbitrates competing formulas per peak, and commits winners with metadata (role, lock status, confidence, tier, provenance).

Once peaks are fetched, the pipeline is offline (no further network I/O is needed beyond `match_compounds` scoring calls), deterministic (content epoch pinned), and SDK-native.

### 1.2 The ledger

`ledger.new_ledger(peaks: pd.DataFrame) -> pd.DataFrame` creates a fresh ledger from the raw peaks table:
- Deduplicates by `peak_id` (keeps highest intensity), so every row is a unique physical peak.
- Initializes all assignment columns to defaults: `role=ROLE_UNEXPLAINED`, `locked=False`, `neutral_formula=NA`, etc.

A ledger row carries (among others): `peak_id`, `mz`, `height`, `area`, `role`, `locked`, `neutral_formula`, `adduct`, `ion_formula`, `ion_score`, `compound_score`, `ppm_error`, `confidence`, `tier`, `tier_reason`, `candidate_density`, `pass_no`, `method`, `commentary`, `alternatives` (JSON), `isotopologues` (JSON), `parent_peak_id`, `iso_label`, `iso_match_score`, `synthetic`, `host_peak_id`, `assigned_fraction`, `below_assignability`, plus time-series stamps (`ts_cv_norm`, `ts_r_mono`, `ts_r_formic`, `ts_disposition`).

### 1.3 Roles (mutually exclusive, enforced by invariants)

| Role | Constant | Meaning |
|---|---|---|
| Unexplained | `ROLE_UNEXPLAINED` | Initial state; no assignment yet. Only this role may be claimed by passes 2–5 (`claim_unexplained_only=True`). |
| M0 | `ROLE_M0` | Monoisotopic owner; owns `neutral_formula` + `adduct`; may own iso_children. |
| Iso child | `ROLE_ISO` | Isotopologue child; points to a parent ROLE_M0 via `parent_peak_id`; carries `iso_label` / `iso_match_score`. |
| Reagent | `ROLE_REAGENT` | Reagent-ion cluster (e.g. `[Br3]-`, `[urea_n+H]+`); labeled BEFORE passes so never a candidate. |
| Artifact | `ROLE_ARTIFACT` | Instrumental noise (ringing/shoulder); only unexplained peaks may be reclassified. |

### 1.4 Ledger invariants

- **I2**: every iso_child must have a parent that exists and owns M0 (enforced in `attach_isotopologue`, ledger.py:212).
- **I3/I4**: locked peaks are immutable; commit/clear/displace on a locked peak is refused. An iso_child of a locked parent refuses re-parenting.
- **I5**: every M0 must carry provenance.
- `ledger.validate(ledger) -> list[str]` runs the post-run audit (no duplicate peak_ids, valid roles, I2, I5); returns an empty list when healthy.

### 1.5 Ledger mutation API (ledger.py)

- `commit_assignment(...)` — atomic commit of an M0 to a peak. Enforces I3/I4. Records `pass_no`, `method`, `confidence`, `commentary`, `alternatives`, `isotopologues`. Sets `role=ROLE_M0`, clears `iso_label`/`parent_peak_id`. `overwrite=False` by default.
- `attach_isotopologue(...)` — marks `child_peak_id` as an isotopologue of a parent owning M0 (I2). Sets `role=ROLE_ISO`, `parent_peak_id`, `iso_label` (`13C` / `81Br` / `13C+81Br` / …), `iso_match_score`. Child must be unexplained or `overwrite=True`.
- `clear_assignment(...)` — demotes an M0 owner back to unexplained, orphaning its iso_children. Forbidden on locked peaks. Records reason in commentary.
- `displace_to_isotopologue(...)` — converts a peak that owns M0 into an isotopologue of a stronger parent; re-parents the child's former iso_children with combined labels (e.g. `13C+81Br`). Forbidden on locked peaks. Used in M0-vs-iso arbitration.
- `mark_reagent(...)` / `mark_artifact(...)` / `lock_peaks(...)` — role/lock setters.
- `role_of(...)` / `is_locked(...)` — queries (raise `LedgerError` if peak_id absent).
- `stats(...)` — coverage summary `{n_peaks, by_role, signal_by_role, count_frac_by_role, by_confidence, by_tier}`. Synthetic sub-peaks are excluded from count but included in signal (via `assigned_fraction`).

### 1.6 Tiers (per M0, assigned post-cleanup)

- **Identified** — corroborated, on-trend, plausible chemistry.
- **Candidate** — plausible but uncorroborated / tied / degenerate / off-trend / heteroatom-coincidence-risk.
- **Below assignability** — flagged monsters (e.g. O≥11 + mass-saturated; unconfirmed high fluorine).

### 1.7 The orchestrator `assign.run()`

```
run(sample_id, context='ambient-air', cfg=None, use_cache=True,
    do_pass2/3/4/5=True, ts_peaks=None, adducts=None, log=print,
    checkpoint_dir=None) -> dict
```
Chains: fetch peaks → `new_ledger` → isotope prescan → reagent labeling → pass 0 → pass 1 → calibrate → relabel → iso-envelope (pre-pass-4) → demote_carbon/massgate → pass 2 → pass 3 → pass 4 → pass 5 → reagent sweep → audits → iso-envelope (2nd) → composites → pass 6 (ladder) → iso-envelope (3rd) → cleanup → siloxane → degeneracy → tiers → fluorine demotion → carbon-cluster demotion → ionization-plausibility demotion → speculative-residual demotion → reference-list rescue-verify → optional time-series → `validate` + `stats`. Returns `{ledger, stats, summaries, prescan, problems, module_versions, module_hashes, context, sample_id}`.

---

## 2. Candidate Enumeration and the Scoring Model

### 2.1 Chemistry layer (chemistry.py)

The lowest layer provides exact monoisotopic masses, Hill-notation formula parsing/formatting, neutral-mass and ion-m/z calculation across 11 adducts, DBE accounting, and the candidate grid.

- `neutral_mass(formula)` — monoisotopic mass from `M[]` constants.
- `ion_mz(neutral, adduct)` — `neutral_mass + ADDUCT_SHIFTS[adduct]` (chemistry.py:159–164).
- `dbe(cnt)` — `1 + (C+Si) + (N+P)/2 − (H+F+Cl+Br+I)/2`; O,S divalent contribute 0.
- `seniors_cap` — `C + Si + N/2 + 1` (chemistry.py:184–186).
- `oxygen_ok` — `O ≤ 2·(C+N+S+P) + 4` (chemistry.py:189–205); each O needs 2 skeleton bonds + 4 headroom.
- `dbe_ok` — hard gate: DBE non-negative integer (±1e-9) and ≤ Senior cap (chemistry.py:208–220).

### 2.2 The grid

`enumerate_grid(ranges, mass_min=30, mass_max=900)` emits all valid `(mass, formula)` tuples in an element box:
- Loops C,Si,N,P,F,Cl,Br,I; for each computes Senior cap; loops `DBE∈[0,cap]` and **derives** `H = 2·(1+C+Si) + (N+P) − 2·DBE − halogens` (chemistry.py:266–268) — only integer-DBE neutrals are produced.
- Checks H in bounds; loops S; caps O at `min(grid.O_hi, oxygen_ok limit)`; loops O; filters by mass bounds (30–900 Da).
- Cached via `_grid_cached` (LRU, max 16 boxes), keyed on `(sorted ranges, round(mass_min,3), round(mass_max,3))`.

`candidates_for_peaks(peak_mzs, ranges, adducts, ppm=search_ppm)` pre-filters the grid: for each peak m/z × adduct, `m_neu = peak_mz − shift`, tolerance `m_neu·ppm·1e-6`, binary-search the sorted grid masses, accumulate matched formula strings.

### 2.3 Enumeration parameters

| Gate | Value | Effect / ref |
|---|---|---|
| `search_ppm` | **3.0 ppm** | Grid enumeration window (~8σ of 0.35 ppm instrument accuracy). Reduced 5→3 (v0.8.0), ~1.7× fewer candidates. `match_compounds` independently keeps a **5 ppm** window so real 29Si/81Br satellites are still found; the z-gate owns final rejection. (PassConfig, passes.py:47) |
| Default grid C-max (ambient) | 40, auto-scaled `min(40, max(12, est_max_C+4))` | contexts.py:52; build_ranges |
| Default grid O-max (ambient) | 30 (uronium: 32) | contexts.py:53,205 |
| Uronium grid | C-max 46, O-max 32, max_Si 12, min_C_for{Si}=2 | contexts.py:186–210 |
| Grid mass bounds | 30–900 Da | chemistry.py:238–239,286–287 |
| Pass 1/2 grid | CHO(N) only; S/P/Cl/Br/F/I = [0,0] | build_ranges; heteroatoms enter only via Pass-3 families |

### 2.4 Complexity penalty (heteroatom skepticism)

`complexity_penalty(formula, scale=0.01, cap=0.20)` (chemistry.py:300–308):
- `_COMPLEXITY_WEIGHT = {N:3, S:8, P:25, Cl:50, Br:50, Si:80, I:80, F:30}`.
- Penalty = `min(Σ weight[el]·count[el] · 0.01, 0.20)`. CHO forms face zero prior; Br/Cl/Si/I need a 0.05–0.20 eff_score margin over CHO to win.

### 2.5 The arbitration scoring model (`arbitrate`, passes.py:242–384)

A **pure function** selecting a single best M0 owner per observed peak from the flat `match_compounds` table.

**Input** `scored` rows: one per (compound, ion, isotopologue) triplet — `compound_formula`, `ion_formula`, `sample_peak_id`, `ion_score`, `compound_score`, `iso_score`, `iso_label`, `is_base`, `theo_mz`, `ppm_error`.

**Per-peak computation:**
1. Filter anchors: `is_base=True ∧ sample_peak_id notna ∧ ion_score notna` (lines 258–262).
2. `raw_score = min(ion_score, compound_score)`.
3. `eff_score = raw_score − penalty − adduct_penalty − cal_penalty + reflist_prior` where:
   - **penalty** = `_evidence_penalty(...)` (lines 292–314): the complexity prior, **waived** for a heteroatom only if its diagnostic isotope (37Cl/81Br/34S) is Mascope-confirmed AND it is not the reagent element. Monoisotopic elements (N, P, F, Si, I) keep the plain prior. An unconfirmed Cl/Br adds the gate penalty `het_iso_penalty_halogen=0.30` on top of the prior; unconfirmed S adds `het_iso_penalty_S=0.12`.
   - **adduct_penalty** = `minor_channel_penalty=0.12` if the winning adduct is in `minor_channels` (`[M+CO3]-`, `[M+O2]-`, `[M]-.`) (line 335).
   - **cal_penalty** = `0` if uncalibrated, else `max(0, (z − cal_z_accept)) · CAL_ARB_WEIGHT` with `CAL_ARB_WEIGHT=0.04` per σ (line 149, `_cal_offtrend` line 332).
   - **reflist_prior** (the SELECTION PRIOR) = `+reflist_prior` (default `0.04`) ADDED when `compound_formula ∈ cfg.reflist_formulas` — the union of the run's context-active reference peaklists (`reflists.active_lists`, set by `assign.run` from `reflists_active`). A small tie-break: it flips a near-tie (gap < the 0.05 margin) toward a published HOM / known contaminant over a mass-coincidence monster, but **cannot override a clearly-better fit**. Empty set → no-op (the default, so assignment is unchanged unless a list is unlocked).
4. **Winner** = highest eff_score; **tied** if `(winner_eff − runner_eff) < 0.05` (`TIE_MARGIN`, line 345).
5. **Alternatives**: up to 6 runners-up recorded with `eff_score`, `raw_score`, `ppm`.

**Isotopologue children** (lines 372–382): non-base rows (`is_base=False`, `iso_score > 0.4`) attributed to the winning (compound, ion) pair are emitted as `iso_children` with `iso_label` and `iso_score`.

**Output**: `{winners DataFrame[peak_id, neutral, ion_formula, adduct, ion_score, compound_score, raw_score, eff_score, eff_margin, ppm_error, n_iso, tied, alternatives], iso_children DataFrame[peak_id, parent_peak_id, iso_label, iso_score]}`.

### 2.6 Confidence labelling (`confidence_label`, passes.py:119–138)

Grades a winner on raw score + offset-aware mass proximity + isotope count + tie status:

| Label | Condition |
|---|---|
| **High** | `score ≥ tau_high (0.90)` AND `|ppm − cal_mu| ≤ 1.5·ppm_user (≤1.5 ppm)` AND `n_iso ≥ 1` (`require_iso_for_high=True`) AND **not tied** |
| **Good** | `score ≥ tau_good (0.80)` AND `|ppm − cal_mu| ≤ 2·ppm_user (≤2 ppm)` |
| **Low** | `score ≥ tau_low (0.70)` |
| **Suspect** | `score ≥ tau_suspect (0.50)` |
| **Reject** | else (no commit) |

Pre-calibration (`cal_mu=None`) the center is 0 ppm. The method suffix (e.g. `series`, `siloxane`) is appended to the label and preserved by `relabel_confidence`.

---

## 3. The Passes, In Order

`commit_winners(ledger, arb, *, pass_no, method, context, cfg, lock, min_raw_score, confidence_suffix="", claim_unexplained_only=False, only_peaks=None)` (passes.py:437–577) is the shared commit engine for all passes. Gates at commit time:
- Reject if `raw_score < min_raw_score` (line 461).
- Reject if `ppm_error` is NaN (no mass evidence) (lines 465–467).
- **Calibrated mass gate** (lines 472–477): if `cal_mu` set, `z = |ppm − cal_mu|/cal_sigma`; reject if `z > cal_z_pattern (4.0)`; if `cal_z_accept (2.0) < z ≤ 4.0` require pattern evidence (`n_iso ≥ 1` OR `series_like`). Uncalibrated → gate off.
- **Minor-channel corroboration** (lines 481–488): a minor-channel winner with `raw_score < tau_good` commits only if `series_like` OR the same neutral is independently assigned on a primary channel elsewhere.
- **Reagent-halogen decomposition** (`_prefer_adduct_reading`): relabel adduct so the reagent halogen sits in the adduct/cluster, not the neutral.
- Attach iso_children; displace weak non-locked child M0s onto stronger parents within the pass (lines 534–570).
- `lock=True` only in pass 1 (and pass 0 / siloxane via their own `lock_peaks`).

### 3.0 Pass 0 — Known-Species Registry (passes.py:686–830)

**Purpose**: assign high-confidence contaminant series and atmospheric radicals BEFORE the organic CHO/CHON pass, and **lock** them so pass 1 grid fits cannot displace them.

**`_known_species(polarity)` registry** (passes.py:595–672), polarity-specific:
- **Negative-mode Br-CIMS**: `atmospheric {HO2, HNO3, HNO2, HNO4}`, `nitroaromatic {C6H4N2O5, C6H5NO3, C6H5NO4, C7H6N2O5}`, `perfluoroacid {C2HF3O2 … C12HF23O2}`, `chlorinated_paraffin {C10H17Cl3 … C30H3Cl15}`, `contaminant:silanediol {C2H8O2Si … C16H50O9Si8}`.
- **Positive-mode urea-CIMS**: ONLY `organophosphate {C6H15O4P … C21H21O4P}` — the grid itself reaches N-bases and oxygenated VOCs, so positive mode largely skips pass 0 (polarity check at passes.py:603).

**Driver `run_pass0_known`**: scores all known-species formulas via Mascope (isotope model covers 29Si/30Si + reagent halogen). For each base ion, commits those passing **three independent gates**:

| Gate | Value | Ref |
|---|---|---|
| Mass offset-aware | `|ppm_error − prior_offset| ≤ 2.0 ppm` | passes.py:731 |
| Own-81Br twin ratio (`[M+Br]-` adducts only) | `0.5 ≤ h_twin/h_M0 ≤ 1.7` (twin at m/z+1.9979535) — composite detector | passes.py:743–753 |
| Organophosphate channel gate | `≥ 2` distinct ion mechanism_ids on-cal (`|ppm−prior_offset|≤2.0`) — P is monoisotopic, no twin | passes.py:755–761 |
| Chlorinated-paraffin isotope gate | `n_kids ≥ 2` matched 37Cl satellites — Cl is off-grid, isotope envelope alone proves it | passes.py:770–779 |

**Confidence**: Good if `ion_score ≥ 0.7` or `n_kids ≥ 2`, else Low. **All pass-0 commits are locked** (`lock_peaks`, passes.py:808). Iso children attached from `is_base=False` rows with `iso_score > 0.4`. Method = `known:{family}`; commentary records the tag (atmospheric / nitroaromatic / organophosphate / perfluoroacid / chlorinated-paraffin / contaminant). `relabel_confidence` skips locked M0s, so a deliberate Low stays Low — the lock is the intent.

**Recovery pass `_recover_isotope_locked_known`** (passes.py:833–944): for `_RECOVERABLE_KNOWN_FAMS = {"chlorinated_paraffin"}` only (Cl/Br/S families with isotope diagnostics). Re-anchors families the server scored too low to anchor (`sample_peak_id` NaN — e.g. 15N-labelled poly-Cl whose aggregate score collapsed under 14N phantoms + wide envelope). Commits **only when BOTH**: (1) a real unexplained ledger peak within `anchor_tol=2.0 ppm` (offset-aware, `theo·(1+prior_offset·1e-6)`, line 877) of the theoretical M0, AND (2) `≥ min_sats=2` confirmed 37Cl satellites within `sat_ppm=7.0 ppm` (looser than M+1/M+2's 5 ppm). Height floors: M0 `≥ height_floor=20 cps`, satellites `≥ 10 cps`. Cannot fabricate — no peak or no envelope means no commit. Monoisotopic F (perfluoroacid) and P (organophosphate) are NOT recoverable.

### Pre-Pass-1 — Offset Estimation

`io_mascope.estimate_offset(peaks, min_n=8)` computes the median ppm of base-ion server matches (skipping heavy-isotope rows), returning None if `< 8` matches. `cfg.prior_offset = prior if prior is not None else estimate_offset(raw)` (default 0.0). This seeds the offset-aware pass-0 / pass-1 pre-commit gates so a large systematic offset (e.g. −2.45 ppm uronium source) doesn't blind pass 0 to on-trend contaminants. `prescan(ledger)` detects isotope patterns (Br/Cl/S/13C) and caps the grid C-range.

### 3.1 Pass 1 — CHO/CHON Backbone (passes.py:1880–1901)

**Purpose**: assign and **lock** the high-confidence CHO/CHON backbone before calibration.

**Flow:**
1. `_target_peaks` = unexplained peaks with `height ≥ height_cutoff (100.0 cps)`.
2. `build_ranges(profile, pre, include_N=True)` → CHO(N) box; C capped at `min(grid_c_max, max(12, est_max_C+4))`.
3. `_enumerate` → candidates via `candidates_for_peaks` (grid at `search_ppm=3.0`) + optional cheminfo (`use_cheminfo`, default False).
4. `_context_filter` prunes against the profile's O/DBE/H/N/Cl/Br/I ceilings.
5. `IO.score_candidates(mechanism_ids=cfg.mechanism_ids)` → scored per-isotopologue rows.
6. `arbitrate(scored, cfg)` → winners + iso_children.
7. `commit_winners(pass_no=1, method="cheminfo+grid", lock=True, min_raw_score=tau_low=0.70)`. High M0s are **locked**; others unlocked. Pre-calibration confidence is judged against center 0.

**Gates** beyond the shared commit engine: `tau_high=0.90`, `tau_good=0.80`, `tau_low=0.70` (commit floor), `require_iso_for_high=True`, `tied < 0.05 margin` prevents High promotion.

### Pre/Pass-1 Self-Calibration

`calibrate(ledger, cfg)` (passes.py:152–183) runs AFTER pass 1:
1. Select backbone: M0 rows with `ppm_error` notna AND `score ≥ tau_good (0.80)` AND **CHO-CHON only** (`set(parse_formula) ⊆ {C,H,O,N}`). Backbone is **score-selected**, not confidence-selected, so a large offset doesn't exclude the real backbone.
2. `mu = median(ppm_error)`, `sigma = max(1.4826·MAD, cal_sigma_floor=0.25 ppm)`.
3. Require `n ≥ cal_min_n=20`, else return None (calibration stays off, mass gate disabled).
4. Store `cfg.cal_mu`, `cfg.cal_sigma`.

`relabel_confidence(ledger, cfg)` (passes.py:201–237) then re-grades **unlocked** pass-1 M0s against `cal_mu` (vs 0 pre-calibration), preserving the method suffix. At a large offset the whole backbone reads Low pre-calibration; this recovers true High/Good. Locked commits (pass-0 known, siloxane, pass-1 High) are immune. `z_of(ppm, cfg) = |ppm − cal_mu|/cal_sigma` (passes.py:185–192).

### Pre-Pass-4 Demotions and First Iso-Envelope

- `complete_isotope_envelopes` (1st run, see §3.5) — claim full patterns, displace weak M0s.
- `demote_carbon_inconsistent` (passes.py:1284–1342): clears unlocked M0s whose 13C satellite contradicts the carbon count — `C ≥ 8`, reliable 13C `≥ height_cutoff`, `|c_est − n_c| > max(2.5, 0.35·n_c)`, **no Si** (Si skipped at 1306–1307: 29Si M+1 overwhelms 13C). Frees bright peaks for re-assignment.
- `demote_massgate_monsters` (passes.py:1345–1371): clears unlocked M0s with `z > cal_z_pattern (4.0)`.

Both run BEFORE pass 4 so freed peaks are re-offered.

### 3.1b Siloxane ladder (siloxane.py, runs after cleanup)

`assign_siloxane_ladder` claims the PDMS `+C2H6OSi` (+74.019) oligomer ladder, mass-degenerate per rung with high-O CHON fits, using the series spacing + the 29Si/30Si envelope as decisive evidence, and **locks** its commits so the CHON-centric audits can't undo them. **Si-count intensity gate** (`_m1_ratio`, audit rule): the 29Si M+1 must not only be *matched* by the oracle but its **observed (M+1)/(M0) ratio must be ≥ `SI_M1_MIN_FRAC` (0.6) × the predicted `nSi·4.68% + nC·1.07%`** — otherwise the Si count is over-claimed and the commit is skipped. This stopped C₈H₂₆O₅Si₄ being locked over the real HOM C₁₀H₁₈O₁₁ at m/z 393.004 (M+1 ~13% where Si₄ needs ~27%).

### 3.2 Pass 2 — GKA Series Expansion (passes.py:1906–1947)

**Purpose**: iterative greatest-common-addition series expansion from locked M0 anchors. Walks homologous chains (CH2, PDMS `C2H6OSi`, CF2) outward step-by-step, re-anchoring on each round's confirmed members; each proposal is Mascope-scored.

| Gate | Value |
|---|---|
| `series_ppm` | 3.0 ppm (proposal must be within this of anchor ± unit) |
| `series_min_score` | 0.60 (min_raw_score for series proposals) |
| `series_max_iter` | 3 (max outward steps per iteration; stops on no new proposals or 0 commits) |

Commits with `method="gka-series"`, `confidence_suffix="series"`, **not locked**, `claim_unexplained_only=True` (only fills gaps; never displaces a prior commit).

### 3.3 Pass 3 — Contaminant Families & HX-Clusters (passes.py:1953–2091)

**Purpose**: low-quality recovery. Opens the context's contaminant families and resolves HX-cluster artifacts.

**Stage 1 — HX clusters** (`_resolve_hx_clusters`, passes.py:1688–1813): for each locked anchor Y (± 1 CH2 GKA homologs), proposes `Y+HX` scored under `[M+X]-` (e.g. `Y·HBr·Br-` on Br-CIMS — the identical ion to the covalent alias), commits with `neutral=Y`, `adduct=[M+HX+X]-`, method `cluster:Br`. Accepts the base line OR the +2 heavy-isotope line. The `cluster_claimed` set excludes these compositions from the covalent family below.

**Stage 2 — families** (`pass3_families` per context + auto-detected GKA evidence): family-specific element budgets override context caps (e.g. sulfate S 1-1/O 3-4; organosulfate S 1-1/O 3-6; nitrate N 1-2/O 3-8; siloxane Si 1-6/O 1-6/C 2-12/H 6-36; PDMS Si 4-12/O 3-14/C 8-26/H 18-78). **Chain-based enumeration** (passes.py:2020–2050): detected repeat-unit chains (CF2 links, Si-O-Si rungs) open a family and bypass context elemental caps (CF2 chains open `fluorinated` even on ambient `max_F=0`), guarded by series consistency + arbitration priors. Bromo/chloro: drop covalent-X where `X−HX` is an existing anchor or is `cluster_claimed`.

Commits with `method="contaminant:{family}"`, `confidence_suffix=family`, **not locked**, `claim_unexplained_only=True`, `min_raw_score=tau_suspect=0.50`. `bromo_organic`/`chloro_organic` auto-added on Br/Cl reagent.

### 3.4 Pass 4 — Residual (referenced in PassConfig:64–65)

Isotope pairs + series chains, DBE-only plausibility (no `match_compounds`, just formula validation). Gates: `residual_ppm_strict=1.0` on score alone; `residual_ppm_pattern=4.0` only with pattern evidence (isotope partner or ≥2 series anchors). Runs after the pre-pass-4 demotions free bright peaks.

### 3.5 Pass 5 — Completion & Isotope Envelopes

**`run_pass5_completion`** (passes.py:947–1027). **Purpose**: open the known-neutral space after passes 1–4 lock. Two mechanisms:
- **(a) Cross-channel partners**: for each adduct, `ion_mz(neutral, adduct)` of a High/Good assigned neutral, `_peak_near` within `search_ppm` → targets (e.g. `[M+Br]-` partner of a Good `[M-H]-`).
- **(b) Series-gap members**: `formula_add(neutral, "CH2", k)` to find ladder anchors, interpolate missing rungs, `ion_mz`, `_peak_near`. Malformed/invalid additions return None and are skipped.

`score_fn` scores all targets → `arbitrate` → `commit_winners(claim_unexplained_only=True, only_peaks=union)`. The **`completion` method tag grants the pattern-evidence z-band** (z up to `cal_z_pattern=4.0`) because the neutral is already independently assigned (lines 455–456).

**`complete_isotope_envelopes(ledger, cfg, min_rel=0.06, ppm=12.0)`** (passes.py:1045–1140). Runs **3 times** (before pass 4, after audits, after pass 6). Claims the FULL predicted isotope envelope of every committed M0:
1. `isotope_pattern(ion_formula, min_rel=0.06, max_shift=12.0)` (isotopes.py:90–163) predicts `(dmass, rel_intensity, label)` lines ≥ 6% via per-atom convolution, merging within ~3 mDa.
2. Process M0s in **ascending m/z** (line 1070–1071) so a satellite cannot claim a lighter parent.
3. For each predicted line: `line_ppm = 5.0 if dmass < 2.5 else cfg.ppm (12.0)` (line 1096) — tight for M+1/M+2 to separate 13C (+1.00335) from 29Si (+0.99957, 3.8 mDa apart), loose for multi-isotope M+4+ centroids.
4. **Attach** an unexplained peak as iso_child if `0.3 ≤ ratio ≤ 3.5` where `ratio = h_sat/(h_parent·rel)` (line 1110).
5. **Displace** a committed M0 onto the parent only if (line 1118–1128): not locked, confidence not `High`, standalone `score < tau_high`, and `0.45 ≤ ratio ≤ 2.2` (tighter window protects strong victims). Displaced victims' own iso_children are re-parented with combined labels.

| Gate | Value | Ref |
|---|---|---|
| `min_rel` | 0.06 (only claim lines ≥ 6% of M0) | passes.py:1046 |
| M+1/M+2 ppm | 5.0 | passes.py:1096 |
| M+4+ ppm | 12.0 (`cfg.ppm`) | passes.py:1096 |
| Attach ratio | 0.3–3.5 | passes.py:1110 |
| Displace ratio | 0.45–2.2 | passes.py:1128 |
| Locked immunity | locked M0 never displaced | passes.py:1118 |

### Post-Run Audits

**`audit_isotopes`** (passes.py:1374–1527): on unlocked M0s.
- **Br-doublet repair**: two M0s 1.9980–1.9988 apart at height ratio 0.6–1.45 — one is the 81Br isotopologue. If the lighter carries Br → attach the heavier as `81Br` child. If **neither** carries Br → **clear both, ONLY if `cfg.reagent_element == "Br"`** (lines 1424–1440). On non-Br reagents (e.g. 15NO3-) unrelated CHON pairs routinely sit 1.998 apart; clearing both there destroyed 54 real `[M+15NO3]-` M0s (fixed ed2001a).
- **13C sweeper**: attaches obvious unclaimed 13C satellites. Includes twin-satellite fallback (a 13C+81Br satellite counts as carbon evidence, lines 1487–1497) and cross-channel fallback (same neutral assigned High/Good on another channel spares a missing-13C clear, lines 1498–1510).
- **13C carbon-clamp**: same threshold as `demote_carbon_inconsistent`; only fires on a 13C satellite `≥ height_cutoff` (line 1475); Si-bearing skipped.
- **13C completeness**: formula predicting a bright 13C with no peak → cleared.

**`audit_mass_gate`** (passes.py:1530–1569): applies the calibrated mass gate to pre-calibration (pass-1) commits. Clears (never rewrites) M0s with `z > cal_z_pattern (4.0)`, or `cal_z_accept (2.0) < z ≤ 4.0` with no pattern evidence. No-op when uncalibrated. Returns `{cleared_z, cleared_z_noiso, cleared_nan}`.

### Composites and Pass 6

- `detect_composites` (passes.py:1150–1234): flags (does not demote) M0s whose intensity exceeds what their halogen-free M+1 (13C/29Si/15N) implies — `min_m1_rel≥0.06`, `excess_frac≥0.25`, `min_excess≥400 cps`, `ppm=8.0`. Halogen content guessed from the even-shift M+2/M+4 residual. Runs only when `has_halogen_adduct` (in positive mode an even shift is isotope structure, not co-component).
- `split_composites` (passes.py:1237–1281): de-blends — owner keeps `assigned_fraction` of measured height; a synthetic sub-peak `<id>.2` (same m/z, `synthetic=True`, `host_peak_id`) carries the co-component share + halogen guess. Signal conserved.
- **Pass 6 (ladder)**: gapfill homolog/oxidation diagonals; then the 3rd `complete_isotope_envelopes`.

---

## 4. Arbitration & Tiering Rules

### 4.1 Arbitration (recap of §2.5)

Per peak: `eff_score = raw_score − complexity/iso penalty − minor-channel penalty − calibration off-trend penalty + reference-list prior`. Winner = max eff_score; tied if margin `< 0.05`. The isotope-gating gotcha: an unconfirmed heteroatom pays BOTH the prior AND the gate (e.g. S without 34S pays 0.08 + 0.12 = 0.20); a confirmed one pays zero; a reagent element keeps the prior even when ion-confirmed (covalent vs cluster ambiguity), with `_prefer_adduct_reading` recovering the cluster reading post-arbitration. The **reference-list prior** (`+0.04` when the neutral is on an active reference peaklist, §8.4) is a tie-break only — enough to flip a near-tie toward a known literature/contaminant formula, never enough to beat a clearly-better fit.

### 4.2 Tier classification (tiers.py)

`compute_tiers(ledger)` (tiers.py:293–419) / `apply_tiers(ledger)` (tiers.py:422–445). The tier engine **re-calibrates independently** (`tiers._calibrate`, median-centered, offset-aware outlier rejection) on the corroborated CHO/CHON backbone (High/Good + isotope, excluding halogen/Si/S/F, N≤1) — avoiding circular logic; its fit wins the tier verdict if it disagrees with the pass fit.

**Tier gates:**

| Gate | Value | Ref |
|---|---|---|
| `CLOSE_MARGIN` | 0.10 (alternatives within this eff_score count toward density; `candidate_density = 1 + #close`) | tiers.py:64 |
| `O_MAX_IDENTIFIED` | 11 (O≥12 = lattice-monster → Candidate / below-assignability) | tiers.py:65–67 |
| `Z_TAIL_DEMOTE` | 2.6 σ (uncorroborated M0 with `|z|>2.6` → Candidate) | tiers.py:77 |
| `DEGEN_DEMOTE_DENSITY` | 2 (degenerate if `>2` distinct cross-family plausible ions, i.e. ≥3, OR MASS-SATURATED) | tiers.py:105 |

**Identified** (default) when: known/locked species; OR unique in the calibrated window (density=1) with isotope/cross-channel/series support or no close alternatives; OR `O ≤ 11`, mass on-trend (`|z| ≤ 2.6` or corroborated), not mass-degenerate or corroborated.

**Candidate** when any of: base confidence Low/Suspect; `O ≥ 12`; mixed Br/Cl backbone ambiguity; tied without cross-channel/series corroboration; close alternatives (density>1) uncorroborated; background air-ion channel without primary status or corroboration; `|z| > 2.6` uncorroborated; mass-degenerate uncorroborated.

**Below assignability** (flag): `O ≥ 11` AND mass-saturated.

The **degeneracy audit** (degeneracy.py) re-counts distinct plausible ions across ALL chemical families (not just the narrow pass box), catching honest cross-family ambiguity. Corroboration (isotopologue / cross-channel / series) is exactly the evidence that breaks degeneracy because it pins the specific ion. The confidence suffix lifecycle (`Low (series)`, `Low (recovered)`) is preserved through re-grading so the report knows the evidence type.

---

## 5. Cleanup Rules (cleanup.py — Pass 7)

`run_cleanup` runs, in order: (1) `recover_isotope_gated`, (2) `label_bromide_clusters`, (2b) `relabel_reagent_halocarbons`, (3) `flag_ringing_artifacts`, (4) `reclaim_satellites`, (5) `reclaim_envelope_tails`. `demote_unconfirmed_fluorine` is **excluded** here — it must run AFTER `tiers.apply_tiers` (which re-promotes), called by `assign.run` post-tiering. All functions mutate the ledger in place (role/tier/commentary only); no new peaks.

### 5.1 Ringing artifacts (cleanup.py:54–87)

`flag_ringing_artifacts(factor=100.0, dmz=0.012, min_parent=50000.0)`: marks an unexplained peak as `ROLE_ARTIFACT` when a saturating parent (`≥ 50000 cps`) sits within `±0.012 Da` (~4.4 ppm at m/z 400, sub-resolution) and is `≥ 100×` brighter (satellite <1% of parent). Below 50k cps a peak cannot ring; the 100× gate distinguishes a true sidelobe from a resolved neighbor.

### 5.2 Bromide clusters (cleanup.py:100–186)

`label_bromide_clusters(defect_max=-0.16, covalent_tau=0.70)`: labels strongly negative mass-defect peaks (`defect < -0.16`, i.e. ≥2 Br in the adduct region) carrying a Br isotope partner (1.998 Da, ratio 0.4–3.0) as `ROLE_REAGENT`. If an oracle is available, scores degenerate covalent di-/tri-bromo organics over box `C0-12 H0-22 N0-1 O0-8 S0-1 Cl0-1 Br1-3`; a fit `≥ 0.70` is recorded as a commentary alternative (reagent-adduct reading preferred per the fewest-halogens-in-neutral policy). Three-way honesty: tie found → record; oracle present, no tie → "not above threshold"; offline → defect-only note.

### 5.2b Reagent-precursor / brominated-background halocarbons (`relabel_reagent_halocarbons`)

A bromomethane reagent-precursor fragment (CH₂Br₂ → CHBr₂⁻, m/z 170.845) is **mass-degenerate** with an absurd bare-element + reagent-cluster reading (neutral `C` via `[M+HBr+Br]-`): the *same ion at the same mass*. Scoring ties exactly and the neutral-halogen complexity penalty then hands the win to the bare-element cluster, so the report names "neutral C". This step catches them on the **invariant ion composition** (parsed element counts, independent of the guessed neutral and robust to ion-formula string ordering): `CHBr2`/`CBr3` → `ROLE_REAGENT` (out of the analyte pool); `C2HBr2O2` → renamed to the real neutral **C₂H₂Br₂O₂ (dibromoacetic acid)** as `[M-H]-` with a background note. The ≥2 bromines are isotope-confirmable and the match is exact-composition, so this is safe. **Br-CIMS only** (no-op for other reagents). Reference: the F-monster / carbon-cluster "background" ions are *not* registered here — their composition (monoisotopic F, or bare carbon clusters) is unconfirmable, so they are left to the plausibility scan / `demote_unconfirmed_fluorine` rather than asserted as named species.

### 5.3 Isotope-gated recovery (cleanup.py:223–311)

`recover_isotope_gated(score_floor=0.65, z_max=2.5)`: revives low-complexity CHO (± ≤1 covalent Br/Cl) molecules dropped by an aggressive score gate, **only when the measured halogen isotope envelope confirms the halogen count**. Enumerates `RECOVERY_BOX='C0-20 H0-36 O0-12 Cl0-1 Br0-2'`, scores via oracle, filters by calibration (`|z| ≤ 2.5`, default σ=0.5 if uncalibrated), then `_pattern_ok`:
- 1Br: `0.78 ≤ M+2/M0 ≤ 1.20`; 2Br: `1.55 ≤ r2 ≤ 2.35 ∧ 0.55 ≤ r4 ≤ 1.35`; 1Cl: `0.20 ≤ r2 ≤ 0.48`; 1Br1Cl: `1.10 ≤ r2 ≤ 1.55`; `(0,0)` → False (no corroboration).
- Rejects any fit with N/S/P (`_het_types ≤ 2`, CHON/CHOS rejected — halogen isotope confirms the adduct halogen, not N/S). Commits `tier=Identified`, confidence "Good (recovered)".
- `_decompose` reverse-maps the ion to `(neutral, adduct)`, deterministically preferring fewest halogens in the neutral (Br-in-adduct `[M+Br]-` over covalent `[M-H]-`), iterating `RECOVERY_ADDUCTS=['[M-H]-','[M+Br]-']`.

### 5.4 Satellite reclaim (cleanup.py:342–395)

`reclaim_satellites(ppm=6.0)`: attaches clean monoisotopic satellites (13C, 81Br, 37Cl) of assigned M0s as iso_children. 13C gate carbon-aware: `0.3·(nC·0.0107) ≤ ratio ≤ 2.5·(nC·0.0107) ∧ ratio < 1.0`; 81Br `0.55 ≤ ratio ≤ 1.4·nBr`; 37Cl `0.18 ≤ ratio ≤ 0.5·nCl`. Touches only unexplained rows; never demotes M0s. Deltas: 13C=1.003355, 81Br=1.9979521, 37Cl=1.997050.

### 5.5 Envelope tails (cleanup.py:398–460)

`reclaim_envelope_tails(ppm=6.0)`: attaches deep multi-halogen envelopes (`k=2..10` of 37Cl/81Br) via binomial gate `0.35·C(nX,k)·p^k ≤ ratio ≤ 2.8·...` (`p_Cl=0.3199`, `p_Br=0.9728`). **KNOWN LIMITATION**: a no-op on real batches (the deep-tail leak is absorbed upstream by `reclaim_satellites` + isotope-locked known-species/CP recovery). Kept harmless; only synthetic tests exercise it.

### 5.6 Fluorine demotion (cleanup.py:466–502)

`demote_unconfirmed_fluorine(f_min=4)`: demotes M0s on unconfirmed high fluorine (`F ≥ 4`) from Identified → Candidate + sets `below_assignability` (19F is monoisotopic, no twin). Exempts known PFCAs (`CnH F(2n-1) O2`, n≥2) and any fit with a Cl/Br/S anchor whose diagnostic isotope is **CONFIRMED** — `34S`/`37Cl`/`81Br` present in the row's `isotopologues` (not merely the element in the formula); a *reagent* Br adduct's `81Br` does **not** count (it confirms the adduct, not the neutral). On the merged ledger (no isotopologues) it falls back to element-presence. **Must run after `apply_tiers`** so demotion sticks. (Audit rule-gap 1 — caught the F7+S monster C₈H₁₃F₇N₂O₄S that was exempted only because S was in the formula.)

### 5.6b Carbon-cluster demotion (`demote_implausible_carbon`)

`demote_implausible_carbon(hc_max=0.35)`: the F-free counterpart of the fluorine demotion. An M0 whose neutral is **F-free with H/C below 0.35** (e.g. C₂₇H₈ at 0.30, C₃₆H₆O at 0.17) is a high-mass coincidence, not a real organic-aerosol molecule (real SOA sits at H/C ≈ 1–2) → Identified → Candidate + `below_assignability`. Same arithmetic as the plausibility "carbon-rich" flag, applied to the tier so the demotion is consistent; the fluorine-rich low-H/C case is left to `demote_unconfirmed_fluorine`. **Runs after `apply_tiers`** (called by `assign.run` right after the fluorine demotion).

### 5.6c Ionization-plausibility demotion (`demote_implausible_ionization`)

A neutral can only be detected on a channel its chemistry supports: `[M-H]-` needs an **acidic proton**, the anion-cluster adducts (`[M+Br]-`/`[M+CO3]-`/`[M+NO3]-`/`[M+HSO4]-`/`[M+CHO2]-`/…) need an **H-bond donor / polar site**. A **pure hydrocarbon** (no O/N/S/P/halogen/Si) has neither, so it cannot ionize on these channels — a high exact-mass + isotope score on a C/H ion just confirms the carbon count, not a real analyte. Such M0s (e.g. C₇H₁₀/C₇H₁₂ `[M-H]-`, C₂H₂ `[M+CO3]-`) → Identified → Candidate + `below_assignability`. **Electron attachment (`[M]-.`/`[M+O2]-`) is exempt** — the one route a heteroatom-free, electron-poor species has. Negative-mode anion channels only. **Runs after `apply_tiers`**, between the carbon-cluster demotion and rescue-verify.

### 5.6d Speculative-residual demotion (`demote_speculative_residual`)

Targets **only** `method` starting `residual` (the pass-4 residual explainer / series gap-fill — the speculative tail; pass-0/1/2 grid analytes and known-species are untouched). Demotes Identified → Candidate + `below_assignability` when the commit reached the top tier on weak evidence (audit rule-gaps 2–4): **off-calibration** (`|z| > cal_z_accept`); **uncorroborated multi-N** (`n_iso==0` AND N≥3 — a Br doublet confirms the adduct Br, not the neutral's C/N backbone); a **series gap-fill with no anchors** ("0 supporting anchors" in the commentary); or a **sole minor-background channel** (`n_iso==0`, adduct in `minor_channels`, no primary-channel partner for the neutral). Caught C₆H₅N₃ and C₁₂H₆O. **Runs after `apply_tiers`**, just before rescue-verify; takes `cfg` for the calibration band + minor-channel set.

### 5.6e Reference-list rescue-verify (`reflists.rescue_unexplained_by_reflist`)

Runs **last** in `assign.run` (after the demotions, so it sets its own tier). Matches the still-`unexplained` residual **by mass** against the run's active reference peaklists, then SCORES those specific formulas with the server (`match_compounds`) — turning a literature lead into a verified ID or a refutation. Decision per matched peak (mass gate: server `ion_score ≥ tau_low` AND on-cal `z ≤ cal_z_accept`):
- **isotope-confirmed** → commit literature-anchored M0, tier Identified, `confidence="Good (literature)"`;
- **too dim to confirm** (the predicted ¹³C M+1 `0.011·nC·height` falls below `height_cutoff`, so no satellite *could* show) → commit a low-quality Candidate + `below_assignability`, `confidence="Candidate (literature, dim)"` — the lead is never lost back to `unexplained` (the small-peak rule);
- **bright enough but isotopes absent**, or off-cal / poor score → **left unexplained** (a real mass coincidence).

Soft and provenance-tagged (every commit records the source list); only `ROLE_UNEXPLAINED` peaks are touched, so it never overrides an existing assignment. Active only when `reflists_active` is supplied (the batch pipeline resolves it from metadata); a no-op otherwise.

### 5.7 Amine re-read (cleanup.py:522–618, uronium-only)

`prefer_amine_over_ammonium(ts_peaks=None, r_min=0.7)`: `[M+NH4]+` of CHO neutral X is mass/isotope-identical to `[M+H]+` of amine X+NH3. Re-reads as the protonated amine UNLESS (1) X is corroborated — present as `[M+H]+` or `[M+(CH4N2O)H]+`, or time-series log-Pearson `r ≥ 0.7` co-variation (≥6 points) — OR (2) the amine is valence-impossible (negative DBE → NH4 forced). Touches only `neutral_formula`/`adduct`.

---

## 6. Batch Pipeline (assign_batch.py + sampling.py)

`assign_batch.run(...)` assigns a representative subset SEPARATELY, then offset-aware merges their M0 peaks into a merged ledger.

### 6.1 Sample selection (sampling.py)

**THE RULE — `select_representative_samples(n_time=N_TIME=5, include_max_tic=True)`** (sampling.py:69–117): `N_TIME=5` evenly TIME-spaced samples (both endpoints always included; nearest distinct sample to each `linspace` target) + the max-TIC sample. Falls back to all samples if `n_samples ≤ n_time` or no `datetime_utc`. Adds `role` = `time-grid` / `max-TIC` / `time-grid+max-TIC`.

**COVERAGE — `select_brightest_coverage_samples(coverage_target=0.85, k_max=10, k_min=N_TIME+1=6, height_floor=1000.0)`** (sampling.py:146–213): bins all batch peaks by m/z; a bin is significant if its max height across samples `≥ 1000 cps`; greedily picks samples that are the brightest for the most significant bins until `0.85` of significant bins are covered, bounded `[6, 10]`; pads with richest remaining; adds time-grid endpoints. Role = `coverage-winner` / `time-grid` / `coverage+time-grid` + `bins_won`.

### 6.2 Offset-aware merge (assign_batch.align, assign_batch.py:52–107)

For each selected sample: `assign.run` → save `per_file/<sid>_ledger.csv` → `estimate_offset` (median ppm of the sample's own matches, None if <8) → extract M0 rows (`_M0_COLS = [mz, neutral_formula, adduct, tier, ion_score]`).

`align(per_file, tol_ppm=DEFAULT_TOL_PPM=6.0, offsets)` (pure):
1. Per file `_mz_adj = mz·(1 − offset/1e6)` (offset SUBTRACTED; a +3 ppm file shifts down).
2. Single-linkage gap clustering (`_cluster_mz`) on `_mz_adj` at 6.0 ppm.
3. Per cluster: pick best assignment by `TIER_RANK = {Identified:2, Candidate:1}` then `ion_score`; consensus `mz` = mean of raw; `n_files`; `mz_jitter_ppm_raw` and `mz_jitter_ppm_caldj` (offset-corrected residual); `formula_agree` (True if all files report the same neutral — different adducts of one neutral still agree).

**Polarity-aware cleanup**: for positive uronium, `cleanup.prefer_amine_over_ammonium(merged, ts_peaks, r_min=0.7)` at the merged level.

**Outputs**: `merged_ledger.csv` (root), `batch_summary.json` (offsets, tier counts, `n_in_all_files`, `formula_disagreement_count`), `tables/jitter.csv` (one row per cluster×file), `tables/selected_samples.csv`, `per_file/*.csv`.

---

## 7. Time-Series Clustering & the Unexplained Funnel (clustering.py, cluster.py, timeseries.py, analyte_viz.py)

`cluster_batch(...)` (clustering.py:39–315) reads `merged_ledger` (M0 only) + `per_file/*` (all roles) + batch TS; produces three figure sets + QC + `clusters_summary.json` (which documents every threshold and the funnel counts — required by the PDF report).

### 7.1 Binning (timeseries.py)

`build_matrix(peaks, tol_ppm=5.0)` pivots per-sample peaks into a samples × m/z-bin intensity matrix (gap-cluster bins by 5 ppm, weighted-average m/z, summed height). `reagent_total(tol_ppm=8.0)` sums reagent-ion bins; `normalize` divides by per-sample reagent total (0→NaN). `auto_bin_minutes(target_bins=50)` now defaults to NATIVE per-sample cadence (`ceil(median spacing)`, floored at 1 min) — the old 29-min binning aliased into spurious empty bins.

### 7.2 Assigned-analyte path

1. Build ion-mz map from M0 (formula|adduct → mz); extract reagent mz from per-file ledgers.
2. Median + CV per channel; **brightness gate `FLOOR_DEFAULT=200.0 cps`** (median) and `≥8` finite points.
3. `correlate` (log10, Pearson, `MIN_POINTS=8`) on RAW traces (preserves multi-channel sums).
4. `cluster(dist_t=DIST_T=0.40, link='complete', min_members=MIN_MEMBERS=3)` — cut at `r > 0.60`.
5. `merge_similar(merge_r=MERGE_R=0.85, complete linkage on centroids)` — fold near-duplicate families.
6. `split_flat_clusters(range_min=FLAT_CLUSTER_RANGE=1.4)` — demote families whose member-mean doesn't move (smoothed max/median < 1.4) to background.
7. `big_changers(fold_min=BIG_CHANGE_FOLD=3.0, baseline=10th-percentile)` — standalone large changes regardless of co-variation.
8. `panel_median` fills below-detection NaN to the detection floor (NOT dropped — nanmedian would have survivorship bias through zero-air dips, cluster.py:67).

### 7.3 Unexplained funnel (the gates)

A TS bin enters unassigned clustering only if: **median `≥ 50.0 cps`** AND `≥8` finite samples AND NOT within **`8.0 ppm`** of any explained peak (M0 + iso_child + reagent + artifact from per-file ledgers — using all roles, else satellites would falsely look unassigned). Then `split_varying(cv_min=FLAT_CV=0.30, range_min=PEAK_RANGE=1.7, smooth_w=SMOOTH_W=3)` partitions into **varying** (CV ≥ 0.30 OR smoothed max/median ≥ 1.7, catching transient bursts) and **flat**; only varying traces are correlated and clustered (`dist_t=0.40`, `min_members=3`) on REAGENT-NORMALISED traces. Flat bins are bunched, not clustered (their shape is noise).

### 7.4 Channel-agreement QC

`channel_agreement(floor=150.0, min_points=8)` (analyte_viz.py:258–305): for each neutral with ≥2 testable channels (median ≥ 150 cps, ≥8 points), correlates all channel pairs (log10). Verdict: `agree` (worst_r ≥ 0.7), `marginal` (≥ 0.4), `disagree` (< 0.4) — QC for the summing assumption.

### 7.5 Time-series annotation

`apply_timeseries(reagent_mzs, mono_anchor_mzs, formic_mz, tol_ppm=5.0, demote=True)` (timeseries.py:185–260) stamps M0s: `ts_cv_norm`, `ts_r_mono`, `ts_r_formic`, `ts_disposition`. Gates (on reagent-NORMALISED traces): `background` if `cv_norm < 0.25` (FLAT_CV); `ambient:biogenic-SOA` if `r_mono ≥ 0.70` (COVARY_R); `ambient:acid/oxygenate pool` if `r_formic ≥ 0.90`. If `demote=True` and an Identified di-bromide/CO3 commit is flat → **tier demoted to Candidate** (formula unchanged). `trace(run_dir, query, tol_ppm=5.0)` (timeseries.py:275–330) is the reproducible single-compound query (reads the run's own parquet + merged ledger).

---

## 8. Reagent / Context / Reference-List / Plausibility Model

### 8.1 Reagents (reagents.py)

`build_library(reagent='Br', max_n=4, max_neutral=1)` (reagents.py:83–143) enumerates cluster-ion m/z: bare halide `R_n-` (odd-n closed-shell, even-n radical), `R_n·(H2O/HBr/HF)_k` (`_CLUSTER_NEUTRALS = {H2O, HBr, HF}` — organic acids and HNO3/HNO2 were removed because `[Br+acid]-` IS the analyte `[M+Br]-` channel), and halide oxides. **Isotopologues** enumerated via `combinations_with_replacement` over `_HALOGEN_ISO` (Br: 79Br/81Br; Cl: 35Cl/37Cl; I: 127I) — both BrO twins now in the library. Positive `_build_positive_library` (reagents.py:65–80) makes `[urea_n+H]+` (`_POSITIVE_REAGENTS = {urea: CH4N2O}`, n=1..6, mass = `neutral_mass(R_n+H) − electron_mass`).

`label_reagents(reagent='Br', ppm=15.0, only_unexplained=True)` (reagents.py:146–174): binary-search the library within ±15 ppm; sets `role=ROLE_REAGENT`, records the known `ion_formula`. `reagent_for_adducts(adducts)` (reagents.py:177–195) infers the library key from detected adducts (returns None for `[M+NO3]-` since +NO3- is both reagent and analyte adduct — `resolve(peaks=df)` decides).

### 8.2 Profiles (profiles.py)

`ReagentProfile` (profiles.py:15–26): `name, label, polarity, adducts, normaliser ('reagent'|'tic'), reagent_ion_re, ranges, detect_adduct, context, aliases`. Built-ins BR, UR, NO3, NO3_15N (profiles.py:29–70). Br/NO3 use `normaliser='reagent'` ([Br3]- dominates); NO3_15N uses `'tic'` (15NO3 clusters below the acquisition window); UR uses `'tic'` (positive mode). `register`/`from_dict`/`load_config` (JSON/TOML) support user reagents. `resolve(reagent='auto', peaks, config)` (profiles.py:126–148) looks up by name/alias or auto-detects via `detect_adducts` then polarity.

### 8.3 Contexts (contexts.py)

`ContextProfile` (contexts.py:25–60): Van Krevelen windows, heteroatom caps, grid bounds, `min_C_for`, `reagent_adducts`, `pass3_families`. `filter_by_profile(formula, profile)` (contexts.py:243–301) gate sequence: `dbe_ok` → no-C inorganic allowlist → heteroatom caps → `min_C_for` reagent-alias guard → Van Krevelen windows (C≥3 only; `Heff = H+F+Cl+Br+I`, `Ceff = C+Si`; C1-C2 special caps).

| Context | Key bounds |
|---|---|
| ambient-air | H/C [0.7,2.75], O/C [0,1.5], N/C [0,0.4], DBE/C [0,0.75]; N≤3, S≤1, P=0, F=0, Cl≤2, Br≤2, Si≤1; min_C_for {Br:5, Cl:5, F:3}; pass3 (organosulfate, nitrate, siloxane, amine) |
| uronium | H/C 0.4–2.6, N≤5, Si≤12, max_Cl/Br=0; pass3 (amine, siloxane, pdms, glycol_peg, phthalate) |

Context caps always win over family expansion (e.g. ambient `max_Si=1` clamps any pdms family `Si(4,12)`). `classify_compound(formula)` (contexts.py:356–374) returns class band / oxidation level / heteroatom tags for reports.

### 8.4 Reference lists (reflists.py)

`load_catalog(directory)` loads packaged `*.json` peaklists into `{id: ReferenceList}` (self-describing: id/system/label/data_version/references/provenance/`always_active`/`applies_to_contexts`, `species[{formula, conditions, radical}]`). `resolve_context_tags(*texts)` infers experimental-context tags from run metadata (batch name / label) via `CONTEXT_KEYWORDS` — the "unlock with metadata" step. `active_lists(catalog, context_tags)` returns `always_active` lists (e.g. lab contaminants) + any whose `applies_to_contexts` intersect the tags (polarity NOT filtered — neutral formulas transfer across reagent ion forms). Current lists: `monoterpene_hom_kang2024` (830 HOM, context-gated) + `contaminants_keller2008` (59 organic neutrals, always-active).

A reference list is used in **three** places, all soft and provenance-tagged (a list is never a measurement and never overrides an isotope-scored Identified):
1. **Selection prior** (`arbitrate`, §2.5/§4.1) — a candidate neutral on an active list gets a `+0.04` tie-break in `eff_score`. `assign.run` sets `cfg.reflist_formulas` from `reflists_active`.
2. **Rescue-verify** (`rescue_unexplained_by_reflist`, §5.6c) — match `unexplained` peaks by mass, score the matched formula with the server, and commit (confirmed → M0; too-dim → tentative Candidate; else leave). The assign-time mass→score→commit path. `match_by_mass(tol_ppm)` is the matcher.
3. **Report annotation** (pdf_report) — `match_assigned` corroborates assigned Candidate neutrals by formula; `match_by_mass` lists unexplained leads on the "Reference-list corroboration & rescue" page + `tables/reflist_matches_*.csv`. Read-only, offline, deterministic.

`assign_batch.run` resolves the active lists once from batch metadata and passes `reflists_active` to every per-sample `assign.run` (enabling 1 + 2); the report (3) re-resolves them at build time.

### 8.5 Plausibility (plausibility.py)

`implausible(formula, tier, polarity)` (plausibility.py:35–63) flags **Candidate-tier only** (Identified never second-guessed): `N≥3 ∧ O/C≥1.0` (heteroatom coincidence); `N≥4 ∧ O≥8`; `F≥4` (no isotope twin); `F=0 ∧ H/C<0.35` (carbon-rich); `polarity='+' ∧ (Br>0 ∨ Cl>0)`. `scan(merged, polarity)` (plausibility.py:66–87) groups by neutral, skips any neutral Identified in any channel, returns flagged set for the report Scrutiny sheet. Flagged formulas are KEPT, not removed.

---

## 9. Report, Output Layout, and Determinism

### 9.1 Pipeline orchestration (pipeline.py)

- `run_batch(...)` (pipeline.py:232–286): one-call full pipeline — fetch/reuse TS → resolve profile → select samples → `assign_batch.run` per sample → merge → cluster → Van Krevelen → PDF. Returns `{ctx, assign, cluster, vk, report_pdf}`.
- `make_run_context` / `run_id(batch, when)` (pipeline.py:64–66): one `when` per run → folder name = cover Report-ID = `batch_slug + YYYY-MM-DDTHHMMSSZ`.
- `generate_report(ctx, ts, ...)` (pipeline.py:183–229): offline — pins `SOURCE_DATE_EPOCH` then `cluster_batch` → `van_krevelen_batch` → `pdf_report.build` → `provenance.record_run`.

### 9.2 Determinism

`stamp_source_date_epoch(when=None)` (pipeline.py:69–80) exports `CONTENT_EPOCH = 315532800` (1980-01-01Z) as `SOURCE_DATE_EPOCH`. matplotlib (PNG/PDF metadata) and the xlsx writer stamp this fixed epoch, so figures/tables are a pure function of input data — byte-identical re-runs. Run time appears only as visible cover text + folder name. If unset, matplotlib uses the system clock (non-reproducible).

### 9.3 I/O (io_mascope.py)

- `connect(env_path)` (io_mascope.py:79–96): builds `MascopeClient` from `MASCOPE_URL` + `MASCOPE_ACCESS_TOKEN`; .env search precedence repo-root → cwd → `$MASCOPE_ENV` → `~/.mascope/.env`. Reads token from disk each call (avoids stale-token 401). Legacy-server patch degrades missing `/api/datasets` gracefully.
- `fetch_batch_peaks` / `fetch_peaks(use_cache, CACHE_ROOT=~/.mascope-assign-cache)`.
- `score_candidates(...)` (io_mascope.py:565–626): batches formulas (`MATCH_BATCH=200`, >500 timeouts) scored concurrently (`MATCH_WORKERS=5`); raises on partial unless `allow_partial=True`. `DEFAULT_MATCH_PARAMS = {mz_tolerance:5 ppm, isotope_ratio_tolerance:0.2, peak_min_intensity:0.0, min_isotope_abundance:0.15, min_isotope_correlation:0.7, probable_match_threshold:0.8, possible_match_threshold:0.4}`.
- `flatten_match_tree(tree)` (io_mascope.py:483–539): pure; flattens compound→ion→isotope to one row per triplet; emits `ppm_error` only for genuinely matched peaks. **15N-labelled reagent re-anchor** (`_reanchor_labelled_reagent`, delta=0.997035 Da, label='15N', line 537): moves `is_base` from the phantom all-light M0 to the actual 15N monoisotopic line (the reagent is 100% 15N) — without it the `[M+15NO3]-` channel is dropped.

### 9.4 PDF report (pdf_report.py)

`build(out_dir, ...)` (pdf_report.py:1058–1084) iterates `SECTIONS = [cover, findings, coverage, composition, scrutiny, reference_lists, gka, families, changers, clusters, methods, assignments_table]`. `load_context` (pdf_report.py:57–276) reads all artifacts (merged ledger, per-file ledgers, figures, summaries, reflists, composition stats, plausibility flags), degrading silently on missing artifacts; one failed section renders an error page, the rest continue. Report-flag thresholds: reagent-signal note if `role_signal['reagent'] ≥ 0.05`; amine caveat if `[M+NH4]+`/`[M+(CH4N2O)H]+` present. `compress_pdf(max_px=850, quality=58, min_mb=2.0)` writes an optional companion (primary stays byte-deterministic).

### 9.5 Output layout (paths.RunPaths)

`RunPaths(out_dir)` (paths.py:1–68): `root` (run folder), `figures/` (PNGs), `tables/` (CSV/XLSX), `report/` (PDF), `data/` (bulk TS), `per_file/` (per-sample ledgers). `place(filename)` routes by extension/role; `ROOT_ANCHORS` (paths.py:27) keep `merged_ledger.csv`, `run_manifest.json`, `batch_summary.json` at root.

---

## 10. Module Reference & End-to-End Flow

### 10.1 Module reference table

| File | Role |
|---|---|
| `chemistry.py` | Element masses, formula parse/format, neutral/ion mass, DBE + Senior + oxygen gates, grid enumeration, complexity penalty |
| `ledger.py` | The per-peak DataFrame: create, commit/attach/clear/displace/mark/lock, role/lock queries, validate, stats; invariants I2–I5 |
| `passes.py` | All assignment passes (0–6), arbitration, calibration, confidence labelling, commit engine, iso-envelope completion, audits, composites, pre-pass-4 demotions |
| `assign.py` | `run()` orchestrator wiring the full per-sample multi-pass chain |
| `contexts.py` | ContextProfile (Van Krevelen windows, heteroatom caps, grid bounds, pass3 families), filter_by_profile, classify_compound, contaminant family budgets |
| `reagents.py` | Cluster-ion library (isotopologue-enumerated), `label_reagents`, reagent auto-inference |
| `profiles.py` | ReagentProfile (channels, normaliser, context); register / load_config / resolve (auto-detect) |
| `reflists.py` | Curated reference peaklists + the catalog loader/context-unlock; the **selection prior** set, the **rescue-verify** pass (`rescue_unexplained_by_reflist`), report corroboration (`match_assigned`) and mass-match (`match_by_mass`) |
| `plausibility.py` | Candidate-tier scrutiny flags (heteroatom coincidence, carbon-rich, wrong-mode halogen) |
| `tiers.py` | Tier classification (Identified / Candidate / below-assignability), independent re-calibration, degeneracy demotion |
| `degeneracy.py` | Cross-family ion-density audit and heteroatom-type counting |
| `cleanup.py` | Pass-7 residual reclassification: ringing artifacts, bromide clusters, reagent-halocarbon relabel, isotope-gated recovery, satellite/envelope reclaim, fluorine demotion, carbon-cluster demotion, amine re-read |
| `isotopes.py` | Per-atom isotope-distribution convolution → predicted envelope `(dmass, rel, label)` |
| `assign_batch.py` | Per-file assignment + offset-aware m/z merge into the merged ledger; jitter report |
| `sampling.py` | Sample selection: time-grid (5 + max-TIC) and brightest-coverage |
| `timeseries.py` | Matrix binning, reagent normalization, native cadence, TS annotation/demotion, reproducible single-compound trace |
| `cluster.py` | Correlation, complete-linkage clustering, merge, flat-split, big-changers, panel median |
| `clustering.py` | Cluster orchestrator: assigned + unassigned funnel, channel-agreement QC, `clusters_summary.json` |
| `analyte_viz.py` | Van Krevelen batch figures, channel-agreement |
| `pipeline.py` | Top-level orchestration, RunContext, determinism epoch, report generation |
| `io_mascope.py` | Mascope connect/fetch/score; flatten match tree; 15N re-anchor; offset estimate; adduct detect |
| `paths.py` | RunPaths output-folder layout |
| `pdf_report.py` | 12-section PDF assembly, context loader, optional compression |
| `provenance.py` | Run metadata + hashing into a cross-run registry |

### 10.2 End-to-end flow (text diagram)

```
RAW PEAKS
  │ io_mascope.fetch_peaks → new_ledger()  [role=UNEXPLAINED, locked=False]
  ▼
PRESCAN + estimate_offset → cfg.prior_offset       (offset-aware gates seeded)
  │
LABEL REAGENTS  → ROLE_REAGENT (never candidates)
  ▼
PASS 0  known-species (polarity-specific): mass-offset + 81Br-twin + per-family
        gates → commit + LOCK; recovery re-anchor for chlorinated paraffins
  ▼
PASS 1  CHO/CHON grid → arbitrate → commit_winners(lock=True)  [High M0 LOCKED]
  ▼
CALIBRATE (n≥20 CHO-CHON, score≥0.80) → cal_mu, cal_sigma
RELABEL unlocked pass-1 M0s vs cal_mu
  ▼
ISO-ENVELOPE #1  +  demote_carbon_inconsistent + demote_massgate_monsters
  ▼
PASS 2  GKA series (CH2/PDMS/CF2), iter≤3, score≥0.60, claim_unexplained_only
  ▼
PASS 3  HX-clusters (Y·HX·X-) then contaminant families (chain-opened),
        min_raw_score=0.50, claim_unexplained_only
  ▼
PASS 4  residual (isotope pairs + series, DBE-only plausibility)
  ▼
PASS 5  completion: cross-channel partners + series-gap rungs (pattern z-band)
  ▼
REAGENT SWEEP → AUDITS (audit_isotopes: Br-doublet/13C; audit_mass_gate)
  ▼
ISO-ENVELOPE #2 → COMPOSITES (detect + split synthetic <id>.2)
  ▼
PASS 6  ladder gapfill → ISO-ENVELOPE #3
  ▼
CLEANUP (pass 7): recover_isotope_gated → bromide clusters → ringing artifacts
        → reclaim_satellites → reclaim_envelope_tails(no-op) → SILOXANE(locked)
  ▼
DEGENERACY → TIERS (apply_tiers) → demote_unconfirmed_fluorine (post-tier)
  ▼
[opt] TIME-SERIES annotate/demote
  ▼
VALIDATE invariants + STATS  →  per-sample ledger
```

**Batch layer** (`assign_batch.run`): select samples (5 time-grid + max-TIC, or brightest-coverage) → run the per-sample chain above for each → estimate per-file offset → `align()` clusters M0s by offset-corrected m/z (6 ppm), best by tier→score, flags formula disagreement → merged ledger → (positive uronium) amine re-read → cluster_batch + Van Krevelen + 12-section PDF, all under a deterministic content epoch.

**Pass count**: pass 0 (known) + pass 1 (backbone) + pass 2 (GKA) + pass 3 (contaminants/clusters) + pass 4 (residual) + pass 5 (completion) + pass 6 (ladder) + pass 7 (cleanup) = **8 numbered stages**, with the isotope-envelope completion running 3 times and the post-run audits, composites, degeneracy, and tiering interleaved as shown.
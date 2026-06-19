# AGENT PEAKY — RESUME HERE (updated 2026-06-19)

**Project:** consolidate this Mascope pipeline (assign → bin/TS → cluster → local
isotope-validate → figures) into ONE scalable, shareable Claude Code skill ("peaky").
Memory: `agent-peaky` (+ `mascope-sdk-knowledge`, `mascope-assign-package`).

**Data path:** SDK-over-shell, credentials at `~/.mascope/.env` (MASCOPE_URL +
MASCOPE_ACCESS_TOKEN). The `mcp__mascope__*` MCP tools 401 (stale token) — don't use
them; `io_mascope.connect()` reads the live .env.

**Spine built (2026-06-19):** `mascope_assign/profiles.py` (ReagentProfile Br/Ur +
`resolve(auto, peaks)`), `mascope_assign/pipeline.py` (`run(batch|dataset|peaks,
reagent='auto', stages)` — load + profile + 'matrix' wired + representative-sample
selection), `mascope_assign/sampling.py` (THE RULE, below), `io_mascope` (canonical
.env search + `fetch_batch_peaks`). Tested on the orange batches.

**SAMPLE-SELECTION RULE (set 2026-06-19, `sampling.py`, LIVE in pipeline):** we do
NOT assign a single averaged file — a peak only present part of the run is then
invisible (the uronium 08:20 snapshot at hour 11.3 missed the hour-18-22 event
analytes). Instead assign a representative subset and MERGE by m/z: **5 samples
evenly spaced in TIME** (nearest distinct sample to each of 5 equally-spaced target
times; endpoints always in) **+ the 1 max-TIC sample** (richest spectrum). `run()`
always computes `out['assign_samples']` (table w/ role: time-grid / max-TIC /
both) + `out['assign_sample_ids']`. Selecting in TIME not row-index = a lone late
file in an irregular run still gets a pick. 19 tests in test_sampling.py.

**NEXT (in priority order):**
1. Wire **assign** stage: the sample SELECTION is done — what remains is to loop
   `assign.run` over `out['assign_sample_ids']` and MERGE the ledgers by m/z (the
   merge logic exists in scratch `ts_uro/merge_experiment.py` — fold it in).
2. Fold the **validate** stage into pipeline.py — a THIN layer on the existing
   `isotopes.py` (`isotope_pattern` envelope + `prescan`), NOT the scratch
   `~/mascope-output/assign-dev/isotope_validate.py` (which duplicates it). The only
   new part = scoring predicted M+2 vs the OBSERVED spectrum over the brightest samples.
3. Fold the **cluster** stage (scratch: cluster_analytes/fold_orphans/*_raw) → module,
   profile-driven; keep the <12-sample guard.
4. Add binning **max-width split guard** to `timeseries.build_matrix` (single-linkage
   chains on dense/drifting data; negligible on orange but harden for sharing).
5. SKILL.md entry + tests + de-hardcode → push to GitHub (remote not yet created).

**Test batches** ("Aleksei's workspace"): Orange peeling Br `NH7D3KHzoGcXCycw`,
Orange peeling Ur `WcEpq37OUtyzkwlP` (23 samples each). Parquets cached in
`~/mascope-output/orange/`.

---

# Roadmap — state after v48 (2026-06-16) and what comes next

## ✅ DONE (2026-06-16): URONIUM run through the FULL pipeline (positive mode)
The positive-mode **uronium MODE** is built and the full `assign.run` ran on
`<sample-id>` (2025-10-02 08:20, batch "<batch>", dataset
"<dataset>", <server>). Outputs: `~/mascope-output/uronium-v2/`
(ledger.csv + 11-sheet assignments.xlsx + summary.md + gka.html + FINDINGS.md).
**Result: 701 M0 (606 Identified / 95 Candidate), 165 iso_child, 13 artifact,
398 unexplained; 0 ledger problems; 540 offline tests green.** Channels 474
[M+H]+ / 227 [M+(CH4N2O)H]+ (bright analytes appear in BOTH -> cross-channel
corroboration). What was built:
- **`contexts.uronium`** (alias `urea-cims`/`urea`): `polarity="positive"`,
  N-heavy VK priors (h/c 0.4-2.6, o/c 0-1.5, n/c 0-0.6, dbe/c 0-1.1), grid
  **C46/O32** (new `grid_c_max`/`grid_o_max` ContextProfile fields, read by
  `build_ranges`), max_Si 4 / max_Br=Cl=F 0, pass3 amine/siloxane/glycol/phthalate.
- **urea reagent library** (`reagents._build_positive_library`): [urea_n+H]+ @
  61/121/181/241...; `reagent_for_adducts` returns "urea"; `label_reagents` works
  for it. (No-reagent sample -> ~0 clusters present, as expected.)
- **`assign.run` wiring**: polarity detection from adducts; positive opportunistic
  channels [M+Na]+/[M+NH4]+; `reagent_element` set ONLY for halogens (urea puts no
  halogen in the neutral -> di-bromide/iso-pair/_prefer_adduct inert); urea-cluster
  TS normaliser. `passes._mech_to_adduct`: added the urea diff -> the channel was
  being MISLABELED [M-H]- (it derives the adduct from the ion-vs-neutral element
  diff and fell through to the [M-H]- default). `detect_composites` GATED on a
  halogen adduct (misfires without one). pass-0 known-species is a no-op in
  positive mode.
- **-2 ppm "be aware" (user)**: the source sits at a uniform **-2.45 ppm (sigma
  0.27)**. Reported `ppm_error` stays RAW (you SEE -2.4); the pipeline re-centers
  every quality gate on -2.45. THREE offset bugs fixed (a 0-centered ppm gate
  broke calibration/tiers at a large offset): (1) `passes.calibrate` selects the
  backbone by SCORE not the |ppm|<=2 confidence label; (2) `confidence_label`
  judges ppm vs `cal_mu`; (3) new `passes.relabel_confidence` re-grades pass-1's
  pre-calibration labels post-calibrate; (4) `tiers._calibrate` outlier guard is
  median-relative not 0-centered. Without these the offset collapsed everything to
  Candidate (75/690) and let off-trend monsters survive (e.g. 462=C15H27NO15 at
  +0.32 ppm); fixed -> 606/95 and the monster is correctly cleared.
- **isotopologue attachment** now works (165 iso_child incl. 5 15N / 6 18O / 138
  13C) -> the shortcut's "79/127 unassigned = satellites" artefact is gone.
- **time-series** re-derived on the full ledger (urea-cluster normalised, 1.52 M
  batch peaks): **580 flat-background / 91 intermediate / 30 ambient** (~83% flat
  -- the no-reagent background character the ROADMAP predicted).
- **Br-reference regression: PASSED (live, 2026-06-16).** The v48 sample
  `<sample-id>` is 404 on <server>, but the user pointed out the SAME physical
  sample lives in dataset "<dataset>" / batch "<batch>" at **08:21:59 = `<sample-id>`** (just a different server copy).
  Full run on it: **22/22 flagships present + 0 junk** (offset-aware
  check_flagships). This sample sits at **-1.99 ppm** (vs the old -0.6 ppm copy) --
  a great stress test of the offset-tolerance work, which it passed. It exposed
  TWO more 0-centered gaps, now fixed: (a) **pass-0's |ppm|<=2 known-species gate
  was offset-blind** -> at -1.9 ppm it dropped the silanediol ladder and pass 1
  grabbed 244.97 with the off-trend C5H10O6 junk; fixed with a rough pre-cal
  (`io_mascope.estimate_offset` from the sample's own matches -> `cfg.prior_offset`
  seeds the pass-0 gate); recovered all 5 silanediols + killed C5H10O6. (b)
  **`relabel_confidence` was re-grading LOCKED pass-0 commits** -> the known
  composite silanediol n=4 (35% Mascope score) became "Reject"; now skips locked
  rows (deliberate grades preserved). `check_flagships` is now **offset-aware**
  (bounds judged vs the median-ppm backbone center, so a different-server copy
  passes). +3 tests (io_mascope estimate_offset, passes pass-0 offset gate).

### ✅ DONE (2026-06-16, user-requested): siloxane ladder ASSIGNED + analyte cluster
Outputs now `~/mascope-output/uronium-v5/`. Two follow-ups the user asked for:
**(A) Assign the characterized PDMS/siloxane ladder.** The 36% unexplained was
dominated by a PDMS ladder spaced by **+74.0186 = C2H6OSi** (the server even fit
684 as a Si10 monster). It's mass-DEGENERATE per peak: at the -2.45 ppm offset a
CHON O-monster out-scores the true Si formula (O is free in the complexity prior),
wins arbitration, then a CHON-centric audit clears it -> stuck unexplained. Three
fixes: (1) new uronium-scoped **`pdms` Pass-3 family** + max_Si 4->12 (shared
siloxane family untouched); (2) **calibration-aware arbitration** (`CAL_ARB_WEIGHT`
in passes.arbitrate: penalise off-trend |z|>accept so an off-trend coincidence no
longer wins-then-z-rejects); (3) **carbon clamp skips Si** (29Si dominates the M+1,
not 13C, so it was clearing real siloxanes). The decisive piece is the new
**`siloxane.py` (assign_siloxane_ladder)**: a LATE, LOCKED pass that finds the
+C2H6OSi ladder and commits each member on the **29Si/30Si isotope envelope**
(oracle-confirmed) as a PDMS oligomer, displacing UNLOCKED monsters, Candidate
tier -- bypassing the CHON heuristics that run earlier. Result: **22 ladder
members committed (+44 Si satellites), Si-bearing M0 ~10->29, unexplained
398->342, signal explained 64%->80.7%** (M0 69.6 + iso 11.1). Tier 604 Id / 115
Cand. 560+ offline tests green (new test_siloxane.py).
**(B) Analyte cluster (the peaks changing over the experiment).** The batch is a
**24-h run** (1188 samples, Oct 1 21:01 -> Oct 2 20:59). Urea-cluster-normalised,
cv_norm + hierarchical clustering: **138 changing peaks (cv>=0.3) vs 1084 flat**.
The analyte cluster (`cluster 12`, ~97 members) is **flat for ~17 h then spikes
~4x at hours 18-22**, cleanly separated from the constant flat background. Led by
monoterpene oxidation products **C10H14O (cv 1.36) / C10H16O2 (cv 1.28)** + a large
oxygenated-CHON/amine pool (C12-15 NOx). Tooling `assign-dev/ts_uro/
analyte_cluster.py` -> plot `uronium_analyte_cluster.png` + `uronium_changing_
peaks.csv`. (The siloxane background is correctly NOT in the spiking analyte cluster.)
Residual gaps (minor): 610/653/684 ladder rungs still unassigned (weak 29Si
confirmation); 684 keeps a C44 amine monster. Characterisation in uronium-v2/
FINDINGS.md; superseded outputs uronium-v2/v3/v4/.

## ⭐ LESSONS from the uronium/positive-mode session (2026-06-16)

1. **A large systematic mass offset exposes EVERY 0-centered ppm gate.** The
   pipeline was built for a ~0-offset instrument; the uronium source sits at
   -2.45 ppm and the Br low-temp copy at -1.99 ppm. SIX independent gates were
   silently hard-wired to 0 ppm and each broke differently: (a) `passes.calibrate`
   selected its backbone by the confidence LABEL (which embeds |ppm|<=2) ->
   biased mu toward 0; (b) `confidence_label` judged |ppm| vs 0 -> whole backbone
   read 'Low' -> tier collapsed to Candidate; (c) `tiers._calibrate` had a
   0-centered |ppm|<=2 outlier cut -> tier engine never calibrated; (d) the
   pass-0 known-species gate |ppm|<=2 dropped on-trend contaminants -> the
   silanediol-vs-C5H10O6 collision; (e) **arbitration** scores |ppm| vs theoretical
   (offset-BLIND) so an off-trend mass-coincidence out-scores the true formula,
   WINS, then is z-rejected -> peak left unexplained; (f) `check_flagships` fixed
   bound. RULE: every ppm judgement is relative to the calibrated/rough center,
   never 0. Reported ppm stays RAW (the user still SEES the offset).
2. **Self-calibration runs too late for the early passes.** pass-0/pass-1 commit
   BEFORE the pass-1 self-cal. Fix = a rough pre-calibration from the sample's OWN
   server matches (`io_mascope.estimate_offset` -> `cfg.prior_offset`) to seed the
   early gates; the real fit refines it. `relabel_confidence` then re-grades pass-1
   post-cal (but must SKIP locked commits -- they carry deliberate grades).
3. **Br-specific machinery misfires in positive/no-halogen mode** -- each needs
   gating on the actual reagent/polarity: composite test (M+1), carbon-clamp
   (reads 29Si as 13C), di-bromide/iso-pair, reagent_element complexity prior.
4. **A mass-degenerate contaminant ladder is NOT assignable by the general
   per-peak machinery** -- CHON O-monsters (O is free in the complexity prior)
   out-score the true Si formula, and CHON-centric audits then clear it. It needs
   a DEDICATED, evidence-decisive pass (series spacing + the actual isotope
   envelope) that runs LATE and LOCKS (`siloxane.py`). Generalisable pattern for
   any homologous contaminant family.
5. **Per-file assignment misses the experiment.** `assign.run` assigns ONE sample;
   a peak weak in that snapshot (e.g. an event analyte while the snapshot is
   pre-event) is never assigned. Assign representative files + MERGE by m/z (one
   event-peak file recovered +146 peaks / +21 changing analytes here).
6. **Reagent-normalise the time series ONLY if the reagent is in the mass range.**
   The uronium spectrum starts at 122 m/z, above the main `[urea_n+H]+` ions
   (61/121); normalising to the weak in-range cluster (or analyte-dominated TIC)
   gives closure artifacts. cv-based SELECTION survived (the divisor was flat) but
   the timing did not. Check first; use RAW when no good normaliser exists.
7. **Cross-CIMS comparison is ionisation-selective** -- Br⁻ and urea⁺ detect
   different compound sets (Br⁻ won't ionise low-O monoterpene products); a
   matching formula need not be the same molecule. Don't assume two reagents on
   the same day = the same air without checking shared compounds co-vary.
8. **A reference sample can move servers / change calibration.** The Br v48
   reference (<sample-id>, -0.6 ppm) 404'd on <server>; the same physical sample
   exists in "<dataset>"/08:21 (<sample-id>) at -1.99
   ppm. Regression tooling MUST be offset-aware to survive this.

## ⭐ NEXT TIME (improvements, priority order)
1. **Finish the siloxane ladder** -- 610/653/684 still unassigned (weak 29Si
   confirmation) and 684 keeps a C44 monster. Relax `siloxane.py`'s 29Si gate for
   INTERIOR ladder rungs (membership between two confirmed rungs is itself strong
   evidence) and let the series spacing carry them.
2. **"Experiment-assignment" mode** -- a flag that auto-picks representative files
   (background + TS event-extremes), runs `assign.run` on each, and emits ONE
   merged ledger. Formalise `assign-dev/ts_uro/merge_experiment.py`. (match_compounds
   is per-sample, so a synthetic union can't be scored -- merge real files.)
3. **Time-series-DRIVEN assignment** -- use the TS to find the changing peaks
   FIRST, then target-assign them across files, so experiment analytes are
   prioritised over flat background.
4. **Robust no-reagent normaliser** -- detect when the reagent ion is absent from
   the mass range and fall back to a median-of-flat-bins normaliser (not TIC).
5. **Positive-mode flagship set** -- establish + pin a validated uronium flagship
   list (like the 22 Br flagships) so positive runs have a regression guard.
6. **Generalise offset-awareness** -- audit for any remaining 0-centered ppm
   reference; consider a single `calibrated_ppm()` helper used everywhere.
7. **ladders.py positive mode** -- the pass-6 gap-fill is Br-anchor-specific; make
   the oxidation/homolog diagonals fire on positive adducts too.

## Where the pipeline stands

**LATEST GOOD: v48** (Br-CIMS) on `<sample-id>` (ambient air, <instrument>
2025.10.02). v48 = v46 + time-series findings converted to standing pipeline code
(new timeseries.py, [BrHF]- inorganic cluster, reclaim_satellites, below-
assignability flag+sheet) + the terminal steps applied:
**270 M0 (164 Identified / 106 Candidate — 7 flat di-bromide/CO3 demoted by TS),
iso_child 320, reagent 32, artifact 17, unexplained 330, 26 below-assignability,
22/22 flagships, 0 junk, 540 offline tests green** (was 497; +43 for the positive
uronium mode + offset-tolerant calibration). NB the Br sample is now 404 on <server>
(only the local cache remains) — re-derive flagships from the cached v48 ledger. Outputs `assign-dev/v48/`
(v48_ledger.csv + v48_assignments.xlsx, built by make_v48.py). Run the full
pipeline time-resolved with `run_assignment.py --ts-batch '<batch>'
--ts-dataset '<dataset>'`. (Earlier good: v46 = 471 tests,
21/21 flagships; v47 = consolidation/Orbitool audit.) Run `python3
scripts/check_flagships.py <ledger.csv>` after ANY change — it asserts the 22
validated Br identifications (Br-specific; N/A to uronium) and the banned junk.

### v45 -> v46: Orbitool cross-check fixes (chemistry-gate gaps)
Triggered by comparing our v45 to an Orbitool peak list of the SAME sample (the
Orbitool author flagged that DBE filtering can reject real species). Two real,
narrow gaps fixed (we were NOT making the core error — we threshold DBE on the
NEUTRAL, already assign glycerol, keep isoprene polyols):
- **h_to_c CEILING 2.6 -> 2.75** (ambient/chamber/indoor, contexts.py): admits C3
  saturated polyols (glycerol C3H8O3 H/C 2.67, propylene glycol) in candidate
  generation. They previously slipped in ONLY via the pass-4 iso-pair bypass;
  DBE>=0 already caps H/C at 2+2/Ceff so the ceiling only clipped C3 glycols.
  Effect: glycerol now a pass-1 **High** commit (was "Good (iso-pair)").
- **Nitroaromatics added to pass-0** (passes.py `_known_species` "nitroaromatic"
  family: C6H4N2O5 dinitrophenol, C6H5NO3 nitrophenol, C6H5NO4 nitrocatechol,
  C7H6N2O5 dinitrocresol): H-poor / high-DBE BrC tracers the ambient VK floor +
  DBE/C ceiling block from the grid. Only those present + |ppm|<=2 commit. Effect:
  **dinitrophenol C6H4N2O5 [M-H]- @183.0047 now Identified** (the only M0 add).
- Result: +1 M0 (dinitrophenol), 0 lost, 21/21 flagships, +7 tests (contexts
  31->36, passes 110->112). Comparison + finds reports in `assign-dev/v45/`
  (ORBITOOL_COMPARISON.md, ORBITOOL_GENUINE_FINDS.md, TOP10_DBE_ISOTOPE.md).
- **Orbitool genuine-finds verdict (oracle-confirmed):** Orbitool has NO bright,
  well-scoring organic peak we wrongly drop. The strong SOA (succinic, C6H10O3,
  C8H10O4, C9H14O4) we already assign; the genuinely-missed leads (C5H8O4,
  C10H16O8/O18O8, C9H14O5, C13H20O5, C5H6O3) score only 0.43-0.47 on the oracle
  (mass-only, no 81Br-twin corroboration) -> real compounds but unconfirmed here;
  time-series (#3) is the principled way to claim them, not a looser bar.

### v44 → v45: reagent-formula / `[81BrO]-` fix MATERIALISED (was uncommitted)
The post-v44 reagent fix is now in an output. `reagents.build_library` returns
`(label, mass, ion_formula)` and enumerates BOTH 79/81Br isotopologues for the
BrO/BrO2/BrO3 oxide anions; `ledger.mark_reagent(ion_formula=)` records every
reagent ion's formula. PRINCIPLE: known formula -> assigned, regardless of class
(analyte/contaminant/ion-source). Verified in v45: the `[81BrO]-` twin at 96.91
(2350 cps) is now `reagent`/`BrO-` (was unexplained); **17/30 reagent rows carry
an `ion_formula`** (Br3-, Br2-, BrO-, C2H2Br3O2-, CH2Br3O2-; was 0/29) — the 13
without one are the generic multi-Br clusters where no confident formula scores
above threshold. **Zero regression**: the M0 set is byte-identical to v44 (0
lost, 0 new, 0 tier changes); the only delta is reagent 29→30 / unexplained
330→329 (the one [81BrO]- peak). Outputs in `assign-dev/v44/` retained for the
diff.

### v45 step 0(b): top-10 unexplained <400 redone with DBE + isotope enforced
Done (analysis only, no locked IDs added — see `assign-dev/v45/TOP10_DBE_ISOTOPE.md`
and `assign-dev/dbe_iso_top10_v2.py`). Method validated PRINCIPLE 2 live: the
isotope-envelope check (predicted `isotope_pattern` vs observed M+2/M+4) decisively
re-ranks candidates. KEY RESULTS: (1) **386.98 = C17H12N2O2S [M+Br]- (0.73,
+0.20 ppm on-cal, DBE 13), isotope-confirmed single Br** — a benzothiazole/azo-dye
N2S aromatic *contaminant*; the higher-scoring Br-free C9H12N2O11S2 (0.79) is
ruled out by the 81Br twin (exactly the lesson). Tentative Candidate-grade, NOT a
pass-0 lock (it's a one-off contaminant at 0.73). (2) **365.10**: isotope FLIPS
the old guess — C20H17N2O3S [M-H]- (0.75, on-cal) is Br-free and contradicted by
the 367 twin; the iso-consistent C15H26O5 [M+Br]- (0.77) is -1.75 ppm off-cal, so
1-Br is confirmed but the formula is marginal. (3) **HALF the top-12 are isotope
SATELLITES** of brighter neighbors (388.98=81Br twin of 386.98; 335.07=81Br twin
of the recovered 333.07 C14H22O4 = the deferred M1; 300.75=37Cl twin of a reagent;
398.01=13C twin of an iso_child; the 138.94/140.94/142.93 Br lattice). CONCLUSION
unchanged: sub-400 residual is reagent-bromide lattice + satellites of assigned
peaks, not missed organics; DBE+isotope narrows a couple of "mass-saturated"
peaks to a best guess but adds no locked IDs.

### v43->v44: adversarial-review fixes (workflow wf_51c774a8, 12 confirmed-real)
- **H3** (tiers.py v0.3.0, in pipeline): tier engine consumes `degeneracy_density`;
  uncorroborated mass-degenerate commits capped at Candidate.
- **H1/H2** (cleanup.recover_isotope_gated): recovery restricted to CHO +
  isotope-confirmed covalent halogen (RECOVERY_BOX dropped N/S/P). C14H22O4
  [M+Br]- KEPT (->Candidate), C18H14N2O3S REJECTED (N2+S uncorroborated). A
  broad-competitor-score veto was tried first but over-rejected C14H22O4 (an
  impossible DBE<0 fluorinated fit out-scores it) -> the CHO-only +
  heteroatom-corroboration rule is the right discriminator.
- **H4** (cleanup.flag_ringing_artifacts): MIN_RING_PARENT=50000 cps + 100x ratio
  floor -> the physically-impossible 273.05 flag (parent 6385 cps) dropped; 17
  real artifacts kept (FT sidelobe amplitude scales UP with parent intensity).
- **M4** (cleanup.label_bromide_clusters): covalent-fit oracle check before
  commentary -> 5/18 clusters get "reagent-adduct reading preferred over
  degenerate di-bromo organic (ion C2H2Br3O2- 0.96 / CH2Br3O2- 0.83)".
- DEFERRED: M1 (wire recovery 81Br twins as iso_child -- would make recovery
  read "corroborated" and H3 would then SPARE it from demotion; the adduct-Br
  twin is non-discriminating, so leaving it unattached is correct for tiering);
  M3 (mostly resolved by H4); M5 (report molecular vs non-molecular separately).

### KEY METHOD LESSON (user, 2026-06-13): always enforce DBE + isotope count
The server `cheminfo` enumerator does NOT enforce DBE -- a raw cheminfo re-score
returns valence-impossible "fits" (C5H12F13O DBE=-6.5, F14H11N3O no-carbon) that
look high-scoring but are nonsense. Use the LOCAL grid (DBE/Senior/O-cap enforced)
or filter with `chemistry.dbe_ok`. AND constrain by the isotope-confirmed halogen
count: 386.98/388.98 are a 1-Br pair (M+2/M0=1.09, no M+4) -> only 2 DBE-valid
1-Br candidates, best guess **C17H12N2O2S [M+Br]- (0.73, +0.2 ppm, DBE 13** --
likely a benzothiazole/azo-dye N2S aromatic contaminant). The Br-free 0.79 fit
(C9H12N2O11S2) is RULED OUT by the 81Br twin. So several "mass-saturated" verdicts
were partly an artifact of skipping DBE -- some heavy peaks DO narrow to a real
best guess. TODO next session: redo the top-10 unexplained <400 with DBE+isotope
enforced (top-3 DBE-valid isotope-consistent formulas each).

Pipeline shape: pass 0 (known species, locked, twin-gated) → pass 1 (CHO/CHON
backbone + self-calibration) → pass 2 (GKA series) → pass 3 (evidence-opened
contaminant families) → **iso-envelope completion** → **pre-pass-4 carbon
clamp** → pass 4 (residual iso-pairs/deep series) → pass 5 (known-neutral
completion) → isotope-physics audit → calibrated mass-gate audit →
**iso-envelope completion (2nd sweep)** → **composite detection + de-blending**
→ pass 6 (anchored ladder gap-fill) → **iso-envelope completion (3rd sweep,
post-pass-6: claims the di-bromide [M+HBr+Br]- cores' M+2/M+4 satellites that
pass 6 commits too late for the earlier sweeps)** → **residual cleanup
(cleanup.py: isotope-confirmed low-complexity recovery + bromide-cluster
labelling + ringing/shoulder artifact flagging; new ROLE_ARTIFACT)** →
**degeneracy audit (degeneracy.apply_degeneracy)** → tiers.apply_tiers
**(degeneracy-aware: caps uncorroborated mass-degenerate commits at Candidate)**.

Quality work added 2026-06-13 (single-file, all tested, 21/21 flagships, 0 junk):
mass-error-distribution test + background-CO3-channel gate (tiers.py v0.2.0,
de-risk 179->170 Identified); honest degeneracy measurement (degeneracy.py,
note: unique / mass-degenerate+tie-set / mass-saturated) now CONSUMED by the
tier engine (tiers.py v0.3.0): an uncorroborated commit whose cross-family
degeneracy_density>2 (or MASS-SATURATED) is capped at Candidate -- degeneracy
runs before tiers, corroborated backbones (iso/cross-channel/series) are spared
(v43: 2 self-contradictory "unique"+MASS-SATURATED recovered ions demoted,
172->170 Identified, all 108 corroborated high-degeneracy rows kept);
3rd envelope sweep (+13 di-bromide isotopologues, 60.6->61.4% peaks explained).
INVESTIGATED + EXHAUSTED for sub-400 unexplained: the residual there is the
reagent-bromide lattice, not missed organics. The CH2-ladder "no-fit" peaks =
real homologous families: Family A (di-bromide SOA M+2/M+4) now attached;
Family B (heavy-halogen CH2 ladder, defect -0.10..-0.16) is mass-degenerate AND
its isotope envelopes are owned by neighbors (halogen-counting fails) -> NOT
crackable from one averaged spectrum, needs time-series. NEXT real unlock =
time-series co-variation.

## Built this session (v36–v41) — the isotope machinery

- **Isotope-envelope completion** (`isotopes.isotope_pattern` per-element
  convolution + `passes.complete_isotope_envelopes`, before pass 4 + post-audit).
  ~44% of the bright "residual"/Candidate peaks were ISOTOPE SATELLITES (M+2/13C)
  of brighter peaks, not new compounds. Attaches unexplained satellites +
  DISPLACES weak (non-High, weak-score) M0s onto their true parent. Fixed the
  393/395 bug (silanediol Si4+Br M+2 mis-read as a phantom Cl-F-S). iso_child
  276→304; peaks-explained 58.6→60.6%.
- **Composite-peak detection** (`passes.detect_composites`). The M+1 region
  (13C/29Si) is halogen-FREE so it scales only with the assigned compound; if
  observed M0 exceeds the M+1-implied intensity, an unresolved CO-ELUTING
  compound shares the m/z. Reads the co-component's halogen off the even-shift
  residual (M+2/M0, M+4/M+2 ~ Br/BrCl/Br2). The silanediol n>=3 rungs are
  ~30–45% co-eluting BrCl/Br — formula (Si4, proven by +74.0188 rung spacing)
  AND prediction (binomial, Mascope agrees) are BOTH correct; the PEAK is mixed.
- **Composite de-blending** (`passes.split_composites`, user-designed). The
  owner keeps `assigned_fraction` of its measured height; the co-eluting
  compound becomes a SYNTHETIC `<id>.2` sub-peak at the same m/z (synthetic=True,
  host_peak_id->host, carries co_height + co_halogen). Signal conserved:
  silanediol 393 → host eff 10994 + sub 9092 = 20086 measured. stats() uses
  effective height + excludes synthetic from n_peaks. Signal 91.0→90.2% (the
  honest drop = co-component fractions no longer mis-credited).
- 6 adversarial-review findings fixed in the isotope code (81Br mass constant
  was 0.16 mDa off + inconsistent; dead tier-guard; max_shift truncation; flat
  ppm window mislabeling 29Si/13C). See commit c882594.

## What the residual IS (eliminated — do not re-chase)

- **Sulfate / nitrate**: dead (0.06x / 0.0–0.1x decoys; the sample's N lives in
  the small acids + dinitrophenol, see below).
- **The C/H lattice is biogenic SOA** (v30-v32): bright n_Br=2 peaks are
  mono-/sesquiterpene oxidation products as `[M+HBr+Br]-` reagent clusters
  (covalent-Br2 and a Br2 reagent adduct are the SAME ion). 12 di-bromide CHO
  cores assigned (C9–C15 H_xO_4 ladders, e.g. 409=C15H22O3). All CHO-only.
- **The diagonals are MOSTLY NOT new compounds.** Repeatedly characterized
  (workflows + live scoring): the rotating-GKA diagonals are dominated by
  (a) 81Br/13C ISOTOPE SATELLITES of brighter peaks marching as parallel CH2
  ladders (an isotope envelope of a CH2 series IS a CH2 series), and
  (b) fluorinated-contaminant CH2/C2H4O ladders. A non-integer GKA rotation
  makes contaminant ladders + satellites look like SOA — adversarial
  verification (constant-DBE, same-adduct, satellite guard, dedicated isotope
  envelope, live score) is MANDATORY before encoding any diagonal as chemistry.
- **Di-bromide naming campaign on the remaining nBr=2 residual: ZERO new
  names** (workflow wf_94552dc0). All 14 candidates rejected on adversarial
  verification — including 424.99=C15H22O4, whose M+2/M+4 are OWNED by a
  confirmed single-Br fluorinated contaminant (C14H12F8O [M+Br]-), so its +HBr
  pairing only proves a SINGLE-Br adduct. The remaining nBr=2 residual is:
  tribromide/3-Br reagent clusters (194.95/208.96/280.76/282.76/199.85/238.76),
  multi-Br mass-domain-ambiguous (124.92/140.92/192.95/194.94), one
  multiply-charged (239.03), and uncorroborated N2/S mass-fits with FABRICATED
  isotope envelopes (their apparent M+2/M+4 are crowding from adjacent assigned
  di-bromides). The clean di-bromide SOA already got assigned in v41; what
  remains is not nameable from the sum spectrum.
- **Reagent = dibromomethane (CH2Br2)** → Br-/Br2-/Br3- clusters (Br3- DOMINANT
  ~129k cps) + trace HBr. Server has +Br-/+Br2-/+Br3- mechanisms. [M+Br3]- must
  NOT be a blanket scoring channel (timed out batches, lost 40 base M0s incl.
  TFA). OPEN: the [M+HBr+Br]- vs covalent-monobromo [M+Br]- LABEL choice for
  di-bromide cores is the user's call (same ion; 480-cps HBr.Br- argues for
  covalent-monobromo). Left as-is pending decision.

## Next steps, in priority order

0. **DONE (v45, 2026-06-15):** (a) re-ran to v45 — reagent-formula / `[81BrO]-`
   fix materialised, zero regression, 21/21 flagships; (b) redid the top-10
   unexplained <400 with DBE + isotope enforced — 386.98 -> C17H12N2O2S [M+Br]-
   confirmed, no locked IDs added (see the v45 block at the top + `assign-dev/
   v45/TOP10_DBE_ISOTOPE.md`). **Still open from 0(c):** **M1** — wire the
   recovered-ion 81Br twins (e.g. 335.07 = the twin of 333.07 = C14H22O4 [M+Br]-)
   as `iso_child`. STILL DELIBERATELY DEFERRED: attaching corroborates the
   recovery, so H3 would spare it from the Candidate cap, contradicting the
   intended Candidate outcome. Only do it if the user wants those twins off the
   residual AND accepts the tier consequence. **M3** (symmetric artifact split)
   is largely resolved by H4. Nothing here is blocking; the real remaining unlock
   is time-series (step 3).
1. **Name the composite co-components** (the open de-blending step). The
   `<id>.2` synthetic sub-peaks now sit in the ledger with their m/z + measured
   halogen constraint (the silanediol rungs host ~9k-cps BrCl sub-peaks). Build
   a constrained match on them (halogen pinned, mass fixed) to NAME the
   co-component → this formally allows TWO M0 owners per peak (the data model
   supports it: distinct peak_id, host link, fractional signal). DESIGN CALL the
   user flagged: do we commit the `.2` peak as a second M0?
2. **Add dinitrophenol C6H4N2O5 [M-H]- (183.0047, 808 cps, -0.24 ppm, DBE 6)**
   to the pass-0 known-atmospheric list. A genuine missed nitroaromatic tracer
   (verified credible during the v41 row analysis; the ONE solid new name from
   the diagonal hunting — it is CHON, immune to the di-bromide envelope trap).
   Candidate (tentative): 293.0579 = C15H15ClO4 [M-H]- (respects measured Cl,
   1215 cps, but -2.4 ppm — marginal).
3. **Time-resolved correlation confirmer — DONE (2026-06-16, THE unlock realised).**
   Data arrived on the **<server>** server: dataset "<dataset>",
   batch **"<batch>"** (= Oct-2, the SAME day as the v46 reference
   sample), **1199 samples ~1/min over 24 h**. Pipeline: `assign-dev/ts_*.py`
   (load -> bin 5 ppm -> sample×bin matrix -> **reagent-normalize** by Br3- to
   kill sensitivity drift -> cv_norm + Pearson corr + hierarchical clusters ->
   interpret/verify workflow). FULL FINDINGS: `assign-dev/ts/TIMESERIES_FINDINGS.md`
   + plot `family_timeseries.png`. KEY RESULTS (adversarially verified, wf agents):
   - **Inlet vs ambient cleanly separated**: median |diel_ratio-1| = 0.013 (flat
     background) vs 0.246 (biogenic). Fluorinated (61/64 flat) + silanediol (9)
     "Identified" are inlet/instrument background -> re-tier to contaminant.
   - **Monoterpene biogenic SOA CONFIRMED** (cluster 6, 25 CHO members C10H16/18Ox,
     r_mono 0.94, true intra-family co-variation not shared-diel). v46's biogenic
     assignments validated. BrO-/HO2 are HOx hitchhikers on the morning envelope,
     not SOA members.
   - **C15H24O3 [M+Br]- (331.09): promote Candidate->Identified-Good** (r_formic
     0.96, presence 0.998, unique formula). Cap at Good (mass -2..-3 ppm).
   - **DI-BROMIDE "SOA CORES" OVERTURNED**: every di-Br core (398.98/409.00/440.99/
     342.92/370.95/304.90/358.95) is BRIGHT + FLAT (cv_norm 0.06-0.22, r_mono~0,
     diel~1.0) = stable inlet/instrument background, NOT biogenic SOA. (Flatness
     is not low-S/N: the flattest are the brightest, 33-51k cps.) Keep the mass/
     isotope tier, REMOVE the biogenic-SOA-core label. This corrects the v30-v41
     "C/H lattice = di-bromide biogenic SOA" theme below.
   - The "anti-correlated nighttime" cluster is mostly a **Br3-normalization
     closure artifact** (anti-PHASE, diel~1.0), not night chemistry. Its 11
     unassigned co-varying peaks are single-Br (1.9979 Da) inorganic/halogen
     condition ions, NOT missed organics; 169.954/171.952 are 13C/13C+81Br
     satellites of the Identified C3H6BrO3-.
   - 13 concrete pipeline actions in TIMESERIES_FINDINGS.md §5 (re-tier, promote,
     downgrade sparse, name/decline). NEXT: ingest time-series into the pipeline
     (a `timeseries.py` module: cv_norm + family-correlation -> auto re-tier
     contaminants, promote co-varying Candidates, gate di-bromide labels). Needs
     oracle confirmation on the flagged items (server was offline during the wf).
4. **Below-assignability certificates.** For each bright has-constraints
   residual peak, run the clamped frame×element enumeration and stamp the result
   into the Unassigned sheet ("searched N frames, zero fits → composite/exotic").
5. **Robustness niceties.** Batch-level re-retry before fail-loud in
   score_candidates; the [M+Br]- vs [M+HBr+Br]- label decision (#OPEN above).

## Standing lessons (encode-don't-remember)

- **TIME-SERIES: always reagent-normalize (analyte/Br3-) BEFORE correlating** —
  raw heights carry an instrument-sensitivity/reagent common-mode that inflates
  every correlation (a 234-bin "cluster" of everything). Normalization isolates
  real chemistry. BUT it introduces a CLOSURE artifact: anything that doesn't
  co-vary with the (SOA-tracking) Br3- denominator shows up ANTI-correlated —
  report anti-correlation, not a night/RH mechanism, unless the ABSOLUTE diel
  confirms it. **cv_norm < ~0.22 + flat diel = inlet/instrument background**
  (bright+flat is the textbook stable-source signature); **high cv_norm + high
  intra-family r = real ambient**. This is the inlet-vs-ambient discriminator.
- **Di-bromide "SOA cores" are flat inlet background, not biogenic SOA** (time-
  series, 2026-06-16): they fail the ambient benchmark (cv_norm 0.06-0.22 vs
  0.36-0.64, r_mono~0 vs ~0.95). The v30-v41 "C/H lattice = di-bromide biogenic
  SOA" reading is OVERTURNED for the flat cores — keep their mass/isotope tier
  but drop the biogenic-SOA label. (The mono-Br monoterpene SOA C10H16/18Ox IS
  confirmed ambient.)
- A di-bromide (or any multi-isotope) name needs its OWN DEDICATED isotope
  envelope — M+2/M+4 at the predicted intensity that isn't already owned by an
  adjacent assigned species. A good ppm + an apparent M+2/M+4 is NOT enough; the
  envelope is routinely fabricated by crowding from neighbors (the v41 campaign:
  all 14 di-bromide candidates failed exactly here).
- The diagonals' "well-aligning rows" are mostly ISOTOPE SATELLITES (an isotope
  envelope of a CH2 series is itself a CH2 ladder) + the multi-halogen lattice.
  Good ppm ≠ real; high-DBE / O-monster / het-ignoring mass-fits are fictions.
- The M+1 region (13C/29Si) is halogen-free → it gives a compound's true
  intensity independent of any coincident halogen co-eluter (the composite test).
- Every evidence GATE needs a complementary RECOVERY path.
- Negative evidence (absent satellite) never overrides independent positive
  evidence (agreeing channel, twin satellite).
- Server scoring is authoritative; never trust hand ion-mass arithmetic.
- Coverage metrics reward fiction; flagships + count-based coverage are honest.
- GKA: trust integer-A(R) rows; verify non-integer rotations; half-integer X is
  a feature (C/H lattice / fixed-heteroatom view).

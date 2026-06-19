"""Top-level orchestrator + CLI.

Wires the spine (chemistry, contexts, ledger), the oracle (io_mascope), the
prescan (isotopes), and the three-pass director (passes) into one run and
records a reproducibility manifest with every locked module version.
"""
from __future__ import annotations

import argparse
import json

from . import (chemistry, cleanup, contexts, degeneracy, io_mascope, isotopes,
               ladders, ledger, passes, reagents, residual, series_gka,
               siloxane, tiers, timeseries)

__version__ = "0.3.0"

MODULE_VERSIONS = {
    "assign": __version__,
    "chemistry": chemistry.__version__,
    "contexts": contexts.__version__,
    "ledger": ledger.__version__,
    "io_mascope": io_mascope.__version__,
    "isotopes": isotopes.__version__,
    "series_gka": series_gka.__version__,
    "reagents": reagents.__version__,
    "passes": passes.__version__,
    "residual": residual.__version__,
    "ladders": ladders.__version__,
    "tiers": tiers.__version__,
    "degeneracy": degeneracy.__version__,
    "cleanup": cleanup.__version__,
    "siloxane": siloxane.__version__,
    "timeseries": timeseries.__version__,
}


def _degen_summary(led) -> dict:
    """Compact manifest summary of the degeneracy audit."""
    m0 = led[led["role"] == ledger.ROLE_M0]
    d = m0["degeneracy_density"].dropna()
    if not len(d):
        return {"measured": 0}
    d = d.astype(int)
    return {"measured": int(len(d)), "degenerate_ge2": int((d >= 2).sum()),
            "max_density": int(d.max())}


def _module_hashes() -> dict:
    """sha1 of each module file -- the manifest must pin the EXACT code of a
    run; static version strings are not bumped on every edit (v15-vs-v16
    lesson: two runs differed only via an un-versioned passes.py edit)."""
    import hashlib
    from pathlib import Path
    d = Path(__file__).parent
    return {p.name: hashlib.sha1(p.read_bytes()).hexdigest()[:12]
            for p in sorted(d.glob("*.py"))}


def run(sample_id: str, context: str = "ambient-air", *,
        cfg: passes.PassConfig | None = None, use_cache: bool = True,
        do_pass2: bool = True, do_pass3: bool = True, do_pass4: bool = True,
        do_pass5: bool = True, ts_peaks=None,
        log=print, checkpoint_dir=None) -> dict:
    cfg = cfg or passes.PassConfig()
    profile = contexts.get_context(context)
    client = io_mascope.connect()

    raw = io_mascope.fetch_peaks(client, sample_id, use_cache=use_cache)
    led = ledger.new_ledger(raw)
    adducts = io_mascope.detect_adducts(raw)
    # Polarity is read from the detected adducts (cation forms end with "+"); the
    # context's declared polarity is a cross-check, not the authority (adducts are
    # always detected from the sample -- the SKILL design rule).
    polarity = "positive" if any(str(a).rstrip().endswith("+") for a in adducts) \
        else "negative"
    # opportunistic extra channels, scored only if the server has the mechanism
    # registered (auto-selection covers only the sample's own channels, so the
    # ids must be passed explicitly). Negative (halide-CIMS):
    #   +CO3-  carbonate adducts of aldehydes from lingering air ions
    #   +Br2-  di-bromide reagent-cluster adducts of analytes -- the n_Br=2
    #          "C/H lattice" residual is biogenic SOA seen this way (2026-06-12)
    # NB: +Br3- is intentionally NOT a blanket scoring channel -- adding it lost
    # 40 base M0s (incl. TFA) for 0 gains (heavier per-formula scoring times out
    # batches; Br3- is the dominant reagent ion). Positive (urea-CIMS): the
    # alkali / ammonium adducts the source also produces.
    opportunistic = (["[M+Na]+", "[M+NH4]+"] if polarity == "positive"
                     else ["[M+CO3]-", "[M+Br2]-"])
    extra_channels = [a for a in opportunistic
                      if io_mascope.resolve_mechanism_ids(
                          client, [io_mascope.ADDUCT_TO_MECH[a]])]
    mech_names = [io_mascope.ADDUCT_TO_MECH[a]
                  for a in adducts + extra_channels
                  if a in io_mascope.ADDUCT_TO_MECH]
    mech_map = io_mascope.resolve_mechanism_ids(client, mech_names)
    cfg.mechanism_ids = list(mech_map.values()) or None
    adducts = adducts + [a for a in extra_channels if a not in adducts]
    has_halogen_adduct = any(h in str(a) for a in adducts
                             for h in ("Br", "Cl", "I"))
    # rough mass offset from the sample's own matches -> seeds the pre-calibration
    # pass-0 gate (the pass-1 self-calibration refines it). Without it a large
    # systematic offset is invisible to pass 0.
    prior = io_mascope.estimate_offset(raw)
    cfg.prior_offset = prior if prior is not None else 0.0
    log(f"[run] {len(led)} unique peaks; context={profile.label}; "
        f"polarity={polarity}; prior_offset={cfg.prior_offset:+.2f} ppm; "
        f"adducts={adducts}; mechanisms={sorted(mech_map)}")

    pre = isotopes.prescan(led)
    log(f"[run] prescan {pre.as_dict()}")

    # Label reagent-ion clusters BEFORE the passes so they are never assignment
    # candidates (e.g. [Br3]-, [Br+HBr]-, BrO- in a Br-CIMS sample; [urea_n+H]+
    # in a uronium sample).
    reagent = reagents.reagent_for_adducts(adducts)
    # The arbitration complexity prior is kept on a NEUTRAL halogen only -- a
    # molecular positive reagent (urea) puts no halogen in the neutral, so it has
    # no reagent_element (its [urea_n+H]+ clusters are still labelled via
    # `reagent`). This also keeps _prefer_adduct_reading / the di-bromide iso-pair
    # logic inert in positive mode.
    cfg.reagent_element = reagent if reagent in ("Br", "Cl", "I") else None
    if reagent:
        n_reag = reagents.label_reagents(led, reagent, ppm=12.0)
        log(f"[run] pre-labeled {n_reag} reagent-cluster peaks ({reagent})")

    def _checkpoint(tag):
        if checkpoint_dir:
            from pathlib import Path
            Path(checkpoint_dir).mkdir(parents=True, exist_ok=True)
            led.to_csv(Path(checkpoint_dir) / f"{sample_id}_ledger_{tag}.csv", index=False)

    def _safe(tag, fn):
        # a pass failure (e.g. server 500) must not lose prior passes' work
        import time as _time
        t0 = _time.time()
        try:
            s = fn()
        except Exception as e:
            log(f"[run] {tag} FAILED: {type(e).__name__}: {e}")
            s = {"committed": 0, "locked": 0, "iso_attached": 0, "error": str(e)}
        s["elapsed_s"] = round(_time.time() - t0, 1)
        log(f"[run] {tag} took {s['elapsed_s']}s")
        _checkpoint(tag)
        return s

    # Pass 0: explicit KNOWN species -- instrument contaminants (silanediol)
    # AND small atmospheric acids/radicals (HO2/HNO3/HNO2/HNO4) that the
    # integer-DBE / C>=1 organic grid cannot reach -- locked before pass 1
    summaries = {"pass0": _safe("pass0", lambda: passes.run_pass0_known(
        client, sample_id, led, profile, cfg, adducts, log=log))}
    summaries["pass1"] = _safe("pass1", lambda: passes.run_pass1(
        client, sample_id, led, profile, pre, cfg, adducts, log=log))
    # self-calibrate the mass gate on the pass-1 backbone; all later commits
    # are judged by calibrated z-score instead of a fixed ppm window. Then
    # re-grade pass-1's pre-calibration confidence labels against the fitted
    # center -- otherwise a uniform instrument offset (the uronium source sits
    # at ~-2.4 ppm) leaves the whole backbone mislabeled 'Low' and capped at
    # Candidate.
    passes.calibrate(led, cfg, log=log)
    passes.relabel_confidence(led, cfg, log=log)
    if do_pass2:
        summaries["pass2"] = _safe("pass2", lambda: passes.run_pass2(
            client, sample_id, led, profile, cfg, adducts, log=log))
    if do_pass3:
        summaries["pass3"] = _safe("pass3", lambda: passes.run_pass3(
            client, sample_id, led, profile, pre, cfg, adducts, log=log))
    if do_pass4:
        # claim the full isotope envelope (M+2/M+4 incl. Si/Br/Cl combos) of
        # every committed peak BEFORE pass 4, so multi-isotope satellites (the
        # silanediol Si4+Br M+2 at 395) are attached to their parents instead of
        # being mis-assigned by pass 4's iso-pair stage (a Si4+Br M+4/M+2 ratio
        # of ~0.26 otherwise reads as a phantom Cl doublet).
        summaries["iso_env_pre4"] = _safe(
            "iso_env_pre4",
            lambda: passes.complete_isotope_envelopes(led, cfg, log=log))
        # free the bright lattice peaks that pass 1 grabbed with low-carbon
        # CHON mass-fits (O15 monsters): their 13C satellite contradicts the
        # formula's carbon count, so clear them HERE -- pass 4's carbon-clamped
        # iso-pair enumeration then re-claims them as di-bromide SOA clusters.
        summaries["carbon_clamp_pre4"] = _safe(
            "carbon_clamp_pre4",
            lambda: {"demoted": passes.demote_carbon_inconsistent(led, cfg, log=log)})
        summaries["pass4"] = _safe("pass4", lambda: residual.explain_residual(
            client, sample_id, led, profile, pre, cfg, adducts,
            reagent=reagent, log=log))
    if do_pass5:
        # known-neutral completion: cross-channel partners + series gaps of
        # the assignments made by passes 1-4 (no new formula space)
        summaries["pass5"] = _safe("pass5", lambda: passes.run_pass5_completion(
            client, sample_id, led, profile, cfg, adducts, log=log))

    # post-run audit: apply the calibrated mass gate to pre-calibration commits
    # (pass 1 Lows) and anything that slipped through with no mass evidence
    summaries["audit_iso"] = passes.audit_isotopes(led, cfg, log=log)
    summaries["audit"] = passes.audit_mass_gate(led, cfg, log=log)
    # second envelope sweep: displace any satellite that a later pass still
    # mis-assigned, and attach satellites of the pass-4/5 commits
    summaries["iso_env_post"] = _safe(
        "iso_env_post",
        lambda: passes.complete_isotope_envelopes(led, cfg, log=log))
    # composite detection: flag M0s whose even-shift peaks (M0/M+2) are inflated
    # beyond what their halogen-free M+1 satellite implies -> an unresolved
    # co-eluting compound shares the m/z (the silanediol n>=3 rungs sit on a
    # coincident BrCl/Br compound; formula + prediction are both correct).
    # GATED ON A HALOGEN ADDUCT: the discriminator reads the co-component's
    # identity off the EVEN-shift (M+2/M+4) residual, which is only a halogen
    # signature when a halogen reagent is in play. In positive urea-CIMS (no
    # halogen adduct) the even-shift residual is ordinary 13C2 / 18O / 34S
    # isotope structure, so the test would flag every C-rich peak as composite.
    if has_halogen_adduct:
        summaries["composite"] = _safe(
            "composite", lambda: passes.detect_composites(led, cfg, log=log))
        # de-blend: owner keeps assigned_fraction; the co-eluting compound becomes
        # a synthetic '<id>.2' sub-peak (a target for constrained naming later).
        # split_composites appends rows, so rebind led.
        n_before = len(led)
        led = passes.split_composites(led, cfg, log=log)
        summaries["composite_split"] = {"split": len(led) - n_before}
    else:
        summaries["composite"] = {"flagged": 0,
                                  "skipped": "no halogen adduct"}
        summaries["composite_split"] = {"split": 0}
        log("[run] composite test skipped (no halogen adduct -- even-shift "
            "residual is isotope structure, not a co-component signature)")
    _checkpoint("audit")

    # final reagent sweep: catch any cluster peaks the passes left unexplained
    if reagent:
        n_reag2 = reagents.label_reagents(led, reagent, ppm=12.0)
        if n_reag2:
            log(f"[run] post-labeled {n_reag2} more reagent-cluster peaks")

    # Pass 6: anchored ladder gap-fill -- walk the homolog/oxidation diagonals
    # out from committed anchors and fill the gaps (the C15H22O3->O4->O5 type
    # series). Runs AFTER the audits so their mass/isotope gates don't clear its
    # pattern-evidenced (ladder-membership) completions. Candidate tier only.
    if do_pass5:
        summaries["pass6_ladder"] = _safe("pass6_ladder", lambda: ladders.run_ladder_gapfill(
            client, sample_id, led, profile, cfg, adducts, log=log))

    # Isotope-envelope completion, 3rd sweep -- AFTER pass 6. The di-bromide SOA
    # cores ([M+HBr+Br]-, 2 Br in the ion) commit in pass 6, so the pre-pass-4
    # and post-audit sweeps never see them; their M+2/M+4 satellites (a 1:1.95:
    # 0.95 envelope) were left in the residual. This claims them so the
    # di-bromide families don't leak their own isotope lines into "unexplained".
    summaries["iso_env_post6"] = _safe(
        "iso_env_post6",
        lambda: passes.complete_isotope_envelopes(led, cfg, log=log))

    # Residual cleanup (single-file, honest reclassification of the leftover
    # 'unexplained' residual): isotope-confirmed low-complexity recovery, bromide
    # reagent-cluster labelling, and ringing/shoulder artifact flagging.
    summaries["cleanup"] = _safe(
        "cleanup",
        lambda: cleanup.run_cleanup(client, sample_id, led, profile, cfg, log=log))

    # Dedicated PDMS/siloxane-ladder assignment: the +C2H6OSi (74.019) oligomer
    # ladder is mass-degenerate per peak (CHON O-monsters out-score the true Si
    # formula at the calibrated offset), so the general passes miss it or commit
    # a monster a CHON-centric audit then clears. This claims the ladder using the
    # series spacing + the 29Si/30Si isotope envelope as decisive evidence, and
    # LOCKS it so the audits (already run) can't undo it. Positive sources with a
    # silicone background (uronium); inert where the context forbids Si.
    summaries["siloxane"] = _safe(
        "siloxane", lambda: siloxane.assign_siloxane_ladder(
            client, sample_id, led, profile, cfg, adducts=adducts, log=log))

    # honest mass-degeneracy measurement: stamps degeneracy_density /
    # degeneracy_note with the cross-family competing tie set so a reader sees
    # how identifiable each mass really is. MUST run BEFORE tiers -- the tier
    # engine reads these columns to cap uncorroborated mass-degenerate commits at
    # Candidate (a per-pass "unique in the window" claim is only unique inside
    # its narrow box). Disk-cached grid, so this is ~free after the first build.
    summaries["degeneracy"] = _safe("degeneracy", lambda: _degen_summary(
        degeneracy.apply_degeneracy(led, context=profile.label, log=log)))

    # report tier on every committed assignment: Identified vs Candidate
    # (mechanical rules over evidence columns; ROADMAP 2). Degeneracy-aware: a
    # high degeneracy_density with no isotope / cross-channel / series
    # corroboration is capped at Candidate.
    tiers.apply_tiers(led)
    tc = led.loc[led["role"] == ledger.ROLE_M0, "tier"].value_counts().to_dict()
    log(f"[run] tiers {tc}")

    # OPTIONAL time-resolved disposition: when a batch time series is supplied,
    # reagent-normalise it and stamp each M0 with its temporal behaviour
    # (inlet-flat background vs ambient analyte), demoting flat di-bromide/CO3
    # background commits. Runs LAST so it refines the final tiers. (TS unlock.)
    if ts_peaks is not None and len(ts_peaks):
        # reagent-normalise to the source's reagent ion. For a halide source the
        # timeseries module finds the Br_n- rows in the ledger itself (default);
        # a positive molecular reagent (urea) has no such labelled isotopologue
        # rows, so pass the [urea_n+H]+ cluster m/z explicitly as the normaliser.
        ts_reagent_mzs = None
        if reagent in reagents._POSITIVE_REAGENTS:
            ts_reagent_mzs = [m for (_l, m, _f) in reagents.build_library(reagent)]
        summaries["timeseries"] = _safe(
            "timeseries", lambda: timeseries.apply_timeseries(
                led, ts_peaks, reagent_mzs=ts_reagent_mzs, log=log))

    problems = ledger.validate(led)
    if problems:
        log(f"[run] LEDGER VALIDATION PROBLEMS: {problems}")
    st = ledger.stats(led)
    log(f"[run] stats {json.dumps(st)}")
    return {"ledger": led, "stats": st, "summaries": summaries,
            "prescan": pre.as_dict(), "problems": problems,
            "module_versions": MODULE_VERSIONS,
            "module_hashes": _module_hashes(), "context": profile.label,
            "sample_id": sample_id}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Mascope multi-pass peak assignment")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--context", default="ambient-air")
    ap.add_argument("--ppm", type=float, default=1.0)
    ap.add_argument("--search-ppm", type=float, default=3.0)
    ap.add_argument("--height-cutoff", type=float, default=500.0)
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--no-pass2", action="store_true")
    ap.add_argument("--no-pass3", action="store_true")
    ap.add_argument("--output-dir", default=".")
    args = ap.parse_args(argv)

    cfg = passes.PassConfig(ppm=args.ppm, search_ppm=args.search_ppm,
                            height_cutoff=args.height_cutoff)
    out = run(args.sample_id, args.context, cfg=cfg, use_cache=not args.no_cache,
              do_pass2=not args.no_pass2, do_pass3=not args.no_pass3)
    # report.py will own file outputs; for now write the ledger + manifest
    from pathlib import Path
    od = Path(args.output_dir)
    od.mkdir(parents=True, exist_ok=True)
    out["ledger"].to_csv(od / f"{args.sample_id}_ledger.csv", index=False)
    manifest = {k: v for k, v in out.items() if k != "ledger"}
    (od / f"{args.sample_id}_manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    print(f"wrote {args.sample_id}_ledger.csv + manifest")


if __name__ == "__main__":
    main()

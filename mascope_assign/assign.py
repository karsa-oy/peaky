"""Top-level orchestrator + CLI.

Wires the spine (chemistry, contexts, ledger), the oracle (io_mascope), the
prescan (isotopes), and the three-pass director (passes) into one run and
records a reproducibility manifest with every locked module version.
"""
from __future__ import annotations

import argparse
import json

from . import (chemistry, contexts, io_mascope, isotopes, ledger, passes,
               reagents, residual, series_gka)

__version__ = "0.2.0"

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
}


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
        do_pass5: bool = True,
        log=print, checkpoint_dir=None) -> dict:
    cfg = cfg or passes.PassConfig()
    profile = contexts.get_context(context)
    client = io_mascope.connect()

    raw = io_mascope.fetch_peaks(client, sample_id, use_cache=use_cache)
    led = ledger.new_ledger(raw)
    adducts = io_mascope.detect_adducts(raw)
    # opportunistic background air-ion channels: if the server has +CO3-
    # registered, search and SCORE the carbonate channel too (aldehyde adducts
    # from lingering air ions). Mechanism ids must be passed explicitly --
    # the server's auto-selection only covers the sample's own channels.
    extra_channels = (["[M+CO3]-"] if
                      io_mascope.resolve_mechanism_ids(client, ["+CO3-"]) else [])
    mech_names = [io_mascope.ADDUCT_TO_MECH[a]
                  for a in adducts + extra_channels
                  if a in io_mascope.ADDUCT_TO_MECH]
    mech_map = io_mascope.resolve_mechanism_ids(client, mech_names)
    cfg.mechanism_ids = list(mech_map.values()) or None
    adducts = adducts + [a for a in extra_channels if a not in adducts]
    log(f"[run] {len(led)} unique peaks; context={profile.label}; "
        f"adducts={adducts}; mechanisms={sorted(mech_map)}")

    pre = isotopes.prescan(led)
    log(f"[run] prescan {pre.as_dict()}")

    # Label reagent-ion clusters BEFORE the passes so they are never assignment
    # candidates (e.g. [Br3]-, [Br+HBr]-, BrO- in a Br-CIMS sample).
    reagent = reagents.reagent_for_adducts(adducts)
    cfg.reagent_element = reagent   # arbitration keeps the prior on this element
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

    # Pass 0: known instrument-contaminant series (silanediol/PDMS ladder),
    # locked before the organic grid can mis-claim their peaks
    summaries = {"pass0": _safe("pass0", lambda: passes.run_pass0_contaminants(
        client, sample_id, led, profile, cfg, adducts, log=log))}
    summaries["pass1"] = _safe("pass1", lambda: passes.run_pass1(
        client, sample_id, led, profile, pre, cfg, adducts, log=log))
    # self-calibrate the mass gate on the pass-1 backbone; all later commits
    # are judged by calibrated z-score instead of a fixed ppm window
    passes.calibrate(led, cfg, log=log)
    if do_pass2:
        summaries["pass2"] = _safe("pass2", lambda: passes.run_pass2(
            client, sample_id, led, profile, cfg, adducts, log=log))
    if do_pass3:
        summaries["pass3"] = _safe("pass3", lambda: passes.run_pass3(
            client, sample_id, led, profile, pre, cfg, adducts, log=log))
    if do_pass4:
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
    _checkpoint("audit")

    # final reagent sweep: catch any cluster peaks the passes left unexplained
    if reagent:
        n_reag2 = reagents.label_reagents(led, reagent, ppm=12.0)
        if n_reag2:
            log(f"[run] post-labeled {n_reag2} more reagent-cluster peaks")

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
    ap = argparse.ArgumentParser(description="Mascope 3-pass peak assignment")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--context", default="ambient-air")
    ap.add_argument("--ppm", type=float, default=1.0)
    ap.add_argument("--search-ppm", type=float, default=5.0)
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

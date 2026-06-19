"""One-shot entry point: run the full pipeline and write every output.

    python3 scripts/run_assignment.py --sample-id <ID> --context ambient-air \
        --height-cutoff 100 --output-dir ~/mascope-output/<name>

Produces in --output-dir:
    <ID>_<UTC>_ledger.csv      full per-peak ledger (the source of truth)
    <ID>_<UTC>_assignments.xlsx 11-sheet tiered workbook (Identified /
                                Candidates / Unassigned + legend, styled)
    <ID>_<UTC>_summary.md       narrative summary
    <ID>_<UTC>_manifest.json    reproducibility manifest (module versions, prescan, series evidence, timing)
    <ID>_<UTC>_gka.html         interactive rotating-GKA widget over the residual

Heavy work runs locally against the host Python (mascope-sdk). Run via the
shell MCP, not the cowork sandbox.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))   # the skill dir holds the package

from mascope_assign import assign, passes, report  # noqa: E402
from scripts import gka_widget  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description="Mascope multi-pass peak assignment")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--context", default="ambient-air")
    ap.add_argument("--ppm", type=float, default=1.0)
    ap.add_argument("--search-ppm", type=float, default=3.0)
    ap.add_argument("--height-cutoff", type=float, default=100.0)
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--no-pass2", action="store_true")
    ap.add_argument("--no-pass3", action="store_true")
    ap.add_argument("--no-pass4", action="store_true")
    ap.add_argument("--no-pass5", action="store_true")
    ap.add_argument("--output-dir", default=".")
    # OPTIONAL time-resolved disposition: a batch (+dataset) whose per-sample peaks
    # are loaded and used to classify each assignment as inlet-flat background vs
    # ambient analyte (and demote flat di-bromide/CO3 background). The TS unlock.
    ap.add_argument("--ts-batch", default=None,
                    help="batch name to load as the time series (e.g. '<batch>')")
    ap.add_argument("--ts-dataset", default=None,
                    help="dataset for --ts-batch (defaults to the sample's dataset)")
    args = ap.parse_args(argv)

    cfg = passes.PassConfig(ppm=args.ppm, search_ppm=args.search_ppm,
                            height_cutoff=args.height_cutoff)
    od = Path(args.output_dir).expanduser()
    od.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M")
    base = od / f"{args.sample_id}_{stamp}"

    ts_peaks = None
    if args.ts_batch:
        from mascope_assign import io_mascope  # noqa: E402
        client = io_mascope.connect()
        ts_peaks = client.load_peaks(
            dataset=args.ts_dataset, batches=args.ts_batch,
            matches=False, areas=False, heights=True, average=False,
            confirm_above=None)
        print(f"[ts] loaded {len(ts_peaks)} peaks across "
              f"{ts_peaks['sample_item_id'].nunique()} samples for the time-series step")

    out = assign.run(args.sample_id, args.context, cfg=cfg,
                     use_cache=not args.no_cache, do_pass2=not args.no_pass2,
                     do_pass3=not args.no_pass3, do_pass4=not args.no_pass4, do_pass5=not args.no_pass5,
                     ts_peaks=ts_peaks, checkpoint_dir=str(od / "checkpoints"))
    led = out["ledger"]

    led.to_csv(f"{base}_ledger.csv", index=False)
    report.write_excel(led, f"{base}_assignments.xlsx", out["context"],
                       sample_id=args.sample_id)
    report.write_markdown(out, f"{base}_summary.md")
    manifest = {k: v for k, v in out.items() if k != "ledger"}
    Path(f"{base}_manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    pts = gka_widget.build_points(led)
    Path(f"{base}_gka.html").write_text(
        gka_widget.render_html(pts, args.sample_id, args.ppm))
    # second widget over the UNEXPLAINED residual only -- the honest place to
    # hunt for missed homologous structure (most of the unexplained signal in a
    # halide-CIMS sample sits in Br/Cl isotope doublets and polymer ladders)
    un = led[led["role"] == "unexplained"]
    un_pts = gka_widget.build_points(un)
    Path(f"{base}_gka_unexplained.html").write_text(gka_widget.render_html(
        un_pts, f"{args.sample_id} — UNEXPLAINED residual "
        f"({len(un)} peaks)", args.ppm))

    st = out["stats"]
    expl = 100 * (st["signal_by_role"]["M0"] + st["signal_by_role"]["iso_child"]
                  + st["signal_by_role"]["reagent"])
    cf = st.get("count_frac_by_role", {})
    print(f"\nwrote {base}_*.{{csv,xlsx,md,json,html}} (+ _gka_unexplained.html)")
    print(f"assigned {st['by_role']['M0']} | iso {st['by_role']['iso_child']} | "
          f"reagent {st['by_role']['reagent']} | unexplained {st['by_role']['unexplained']}")
    head = (f"peaks explained {100*(1-cf['unexplained']):.1f}% "
            f"({st['by_role']['unexplained']}/{st['n_peaks']} unexplained)  | "
            if cf else "")
    print(head + f"signal explained {expl:.1f}%"
          + f"  | ledger problems: {out['problems'] or 'none'}")


if __name__ == "__main__":
    main()

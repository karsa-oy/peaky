"""Console entry point for `peaky` (alias: `mascope-assign`).

Subcommands:
  list     discover data on the server (datasets / batches / samples)
  assign   single-sample multi-pass assignment -> ledger/xlsx/md/json/gka.html
  gka      build the interactive rotating-GKA HTML from a ledger CSV (offline)

Run `peaky <cmd> --help` for each. Heavy work runs on the host Python
(this package + mascope-sdk). A Mascope account/token is read from ~/.mascope/.env
(or --env / $MASCOPE_ENV / the process environment).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


# --------------------------------------------------------------------------- #
# shared helpers
# --------------------------------------------------------------------------- #
def _require_creds() -> None:
    """Fail fast with an actionable message if no Mascope creds are resolvable
    (before any expensive work). Process env vars satisfy this without a file."""
    from . import io_mascope as IO

    if os.environ.get("MASCOPE_URL") and os.environ.get("MASCOPE_ACCESS_TOKEN"):
        return
    path = IO._find_env()
    if not os.path.exists(path):
        sys.exit(
            "No Mascope credentials found.\n"
            f"  Looked for: {path}\n"
            "  Fix: copy .env.example to ~/.mascope/.env and fill in MASCOPE_URL +\n"
            "       MASCOPE_ACCESS_TOKEN, pass --env PATH, or export the two vars.")


def _friendly_server_error(e: Exception) -> str | None:
    """Map a raw SDK/HTTP exception to an actionable hint, or None if unrecognised."""
    s = f"{type(e).__name__}: {e}".lower()
    if "403" in s or "attention required" in s or "cloudflare" in s:
        return ("Mascope is rate-limiting you (Cloudflare WAF 403). Wait 15-30 min "
                "with NO traffic, then retry — polling extends the block.")
    if "401" in s or "unauthorized" in s or "forbidden token" in s:
        return ("Authorization failed (401). MASCOPE_ACCESS_TOKEN is likely expired "
                "— refresh it in ~/.mascope/.env.")
    if "no peaks" in s or "no samples" in s or "no batches" in s or "404" in s \
            or "not found" in s:
        return ("Not found on the server. IDs can go stale when a server copy is "
                "renamed — re-fetch fresh names/ids with `peaky list`.")
    return None


def _run_guarded(fn) -> int:
    try:
        fn()
        return 0
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001 — CLI boundary
        hint = _friendly_server_error(e)
        msg = f"\nERROR: {hint}\n  (raw: {type(e).__name__}: {e})" if hint \
            else f"\nERROR: {type(e).__name__}: {e}"
        print(msg, file=sys.stderr)
        return 1


# --------------------------------------------------------------------------- #
# subcommands
# --------------------------------------------------------------------------- #
def cmd_list(args) -> None:
    _require_creds()
    from . import io_mascope as IO

    client = IO.connect()
    if args.what == "datasets":
        ds = IO.list_datasets(client)
        col = "dataset_name" if "dataset_name" in ds.columns else ds.columns[0]
        print(f"{len(ds)} datasets:")
        for v in ds[col].tolist():
            print("  ", v)
    elif args.what == "batches":
        if not args.dataset:
            sys.exit("`list batches` needs --dataset NAME "
                     "(see `peaky list datasets`)")
        bs = IO.list_batches(client, args.dataset)
        cols = [c for c in ("sample_batch_name", "polarity", "status") if c in bs.columns]
        print(f"{len(bs)} batches in {args.dataset!r}:")
        print(bs[cols].to_string(index=False) if cols else bs.to_string(index=False))
    elif args.what == "samples":
        if not (args.batch and args.dataset):
            sys.exit("`list samples` needs --batch NAME --dataset NAME")
        sl = IO.fetch_batch_samples(client, args.batch, dataset=args.dataset)
        cols = [c for c in ("sample_item_id", "sample_item_name", "datetime_utc",
                            "tic", "polarity") if c in sl.columns]
        print(f"{len(sl)} samples in {args.batch!r}:")
        print(sl[cols].to_string(index=False) if cols else sl.to_string(index=False))


def _resolve_reagent(args):
    """Return (adducts, context, note). Forces the analyte channels so a positive
    or sparse-match sample never silently falls back to [M-H]- (wrong polarity).
    adducts=None means 'let assign.run auto-detect from the sample'."""
    from . import profiles

    config = getattr(args, "reagent_config", None)
    if args.adducts:
        return list(args.adducts), (args.context or "ambient-air"), \
            f"forced adducts={list(args.adducts)}"
    if args.reagent and args.reagent.lower() != "auto":
        prof = profiles.resolve(args.reagent, config=config)   # name/alias, no peaks needed
        return list(prof.adducts), (args.context or prof.context), \
            f"{prof.name} ({prof.label})"
    # auto: detect from the sample's own peaks (cached, so assign.run reuses it)
    from . import io_mascope as IO

    client = IO.connect()
    raw = IO.fetch_peaks(client, args.sample_id, use_cache=not args.no_cache)
    try:
        prof = profiles.resolve("auto", raw, config=config)
        return list(prof.adducts), (args.context or prof.context), \
            f"auto-detected {prof.name} ({prof.label})"
    except Exception as e:                           # noqa: BLE001
        return None, (args.context or "ambient-air"), \
            (f"auto-detect found no known profile ({e}); using per-sample adduct "
             "detection — pass --reagent explicitly for a positive/sparse sample")


def cmd_assign(args) -> None:
    _require_creds()
    from . import assign, gka_widget, io_mascope, passes, report

    cfg = passes.PassConfig(ppm=args.ppm, search_ppm=args.search_ppm,
                            height_cutoff=args.height_cutoff)
    od = Path(args.output_dir).expanduser()
    od.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M")
    base = od / f"{args.sample_id}_{stamp}"

    adducts, context, note = _resolve_reagent(args)
    print(f"[reagent] {note}; adducts={adducts}; context={context}")

    ts_peaks = None
    if args.ts_batch:
        client = io_mascope.connect()
        ts_peaks = client.load_peaks(dataset=args.ts_dataset, batches=args.ts_batch,
                                     matches=False, areas=False, heights=True,
                                     average=False, confirm_above=None)
        print(f"[ts] loaded {len(ts_peaks)} peaks across "
              f"{ts_peaks['sample_item_id'].nunique()} samples")

    out = assign.run(args.sample_id, context, cfg=cfg, use_cache=not args.no_cache,
                     do_pass2=not args.no_pass2, do_pass3=not args.no_pass3,
                     do_pass4=not args.no_pass4, do_pass5=not args.no_pass5,
                     adducts=adducts, ts_peaks=ts_peaks,
                     checkpoint_dir=str(od / "checkpoints"))
    led = out["ledger"]

    led.to_csv(f"{base}_ledger.csv", index=False)
    report.write_excel(led, f"{base}_assignments.xlsx", out["context"],
                       sample_id=args.sample_id)
    report.write_markdown(out, f"{base}_summary.md")
    manifest = {k: v for k, v in out.items() if k != "ledger"}
    Path(f"{base}_manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    Path(f"{base}_gka.html").write_text(
        gka_widget.render_html(gka_widget.build_points(led), args.sample_id, args.ppm))
    # second widget over the UNEXPLAINED residual only — the honest place to hunt
    # for missed homologous structure
    un = led[led["role"] == "unexplained"]
    Path(f"{base}_gka_unexplained.html").write_text(gka_widget.render_html(
        gka_widget.build_points(un),
        f"{args.sample_id} — UNEXPLAINED residual ({len(un)} peaks)", args.ppm))

    st = out["stats"]
    expl = 100 * (st["signal_by_role"]["M0"] + st["signal_by_role"]["iso_child"]
                  + st["signal_by_role"]["reagent"])
    cf = st.get("count_frac_by_role", {})
    print(f"\nwrote {base}_*.{{csv,xlsx,md,json,html}} (+ _gka_unexplained.html)")
    print(f"assigned {st['by_role']['M0']} | iso {st['by_role']['iso_child']} | "
          f"reagent {st['by_role']['reagent']} | unexplained {st['by_role']['unexplained']}")
    head = (f"peaks explained {100*(1-cf['unexplained']):.1f}% | " if cf else "")
    print(head + f"signal explained {expl:.1f}%  | "
          f"ledger problems: {out['problems'] or 'none'}")


def cmd_batch(args) -> None:
    _require_creds()
    from . import pipeline as PL

    res = PL.run_batch(batch=args.batch, dataset=args.dataset, reagent=args.reagent,
                       base_out=os.path.expanduser(args.out_dir), ts=args.ts,
                       subject=args.subject, do_report=not args.no_report,
                       config=args.reagent_config, select=args.select,
                       coverage_target=args.coverage_target, k_max=args.k_max,
                       height_floor=args.height_floor)
    ctx = res["ctx"]
    print(f"\n[batch] done -> {ctx.out_dir}")
    if res.get("report_pdf"):
        print(f"  report: {res['report_pdf']}")
    if res.get("report_pdf_small"):
        print(f"  report (small): {res['report_pdf_small']}")


def cmd_report(args) -> None:
    # offline: regenerate cluster figures + Van Krevelen + the PDF report from an
    # existing run folder's ledgers (no assignment, no network).
    from . import pipeline as PL
    from . import profiles as P

    prof = P.resolve(args.reagent)
    run_dir = os.path.expanduser(args.run_dir)
    ctx = PL.RunContext(
        out_dir=run_dir, batch_name=(args.batch or prof.label),
        tag=(args.tag or prof.name), label=prof.label, when=None,
        run_id=(args.run_id or os.path.basename(run_dir.rstrip("/"))),
        generated=(args.generated or ""), profile=prof)
    out = PL.generate_report(ctx, os.path.expanduser(args.ts), subject=args.subject)
    print("wrote", out.get("report_pdf"))
    if out.get("report_pdf_small"):
        print("wrote", out.get("report_pdf_small"), "(compressed)")


def cmd_gka(args) -> None:
    import pandas as pd

    from . import gka_widget

    led = pd.read_csv(args.ledger_csv)
    pts = gka_widget.build_points(led)
    out = args.out or (Path(args.ledger_csv).with_suffix("").as_posix() + "_gka.html")
    Path(out).write_text(gka_widget.render_html(pts, Path(args.ledger_csv).stem, args.ppm))
    print(f"wrote {out}  ({len(pts)} points)")


# --------------------------------------------------------------------------- #
# parser + entry point
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="peaky",
        description="Peaky — reproducible multi-pass formula assignment for Mascope peaks.")
    ap.add_argument("--env", default=None,
                    help="path to a Mascope .env (else ~/.mascope/.env or $MASCOPE_ENV)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    pl = sub.add_parser("list", help="discover datasets / batches / samples")
    pl.add_argument("what", choices=["datasets", "batches", "samples"])
    pl.add_argument("--dataset", default=None, help="dataset (workspace) name")
    pl.add_argument("--batch", default=None, help="sample-batch name (for `samples`)")
    pl.set_defaults(func=cmd_list)

    pa = sub.add_parser("assign", help="assign one sample")
    pa.add_argument("--sample-id", required=True)
    pa.add_argument("--reagent", default="auto",
                    help="reagent profile: auto | Br | Ur | ... — forces the analyte "
                         "channels + default context ('auto' detects from the sample)")
    pa.add_argument("--adducts", nargs="+", default=None,
                    help="explicit analyte adduct channels (overrides --reagent)")
    pa.add_argument("--context", default=None,
                    help="plausibility context (default = the reagent profile's context)")
    pa.add_argument("--reagent-config", default=None,
                    help="JSON/TOML file registering extra reagent profiles")
    pa.add_argument("--ppm", type=float, default=1.0)
    pa.add_argument("--search-ppm", type=float, default=3.0)
    pa.add_argument("--height-cutoff", type=float, default=100.0)
    pa.add_argument("--no-cache", action="store_true")
    pa.add_argument("--no-pass2", action="store_true")
    pa.add_argument("--no-pass3", action="store_true")
    pa.add_argument("--no-pass4", action="store_true")
    pa.add_argument("--no-pass5", action="store_true")
    pa.add_argument("--output-dir", default=".")
    pa.add_argument("--ts-batch", default=None,
                    help="batch name to load as the time series (optional TS step)")
    pa.add_argument("--ts-dataset", default=None, help="dataset for --ts-batch")
    pa.set_defaults(func=cmd_assign)

    pb = sub.add_parser("batch", help="assign + cluster + Van Krevelen + report for a whole batch")
    pb.add_argument("--batch", required=True, help="sample-batch name")
    pb.add_argument("--dataset", default=None, help="dataset (workspace) name")
    pb.add_argument("--reagent", default="auto", help="auto | Br | Ur | NO3 | ...")
    pb.add_argument("--reagent-config", default=None,
                    help="JSON/TOML file registering extra reagent profiles")
    pb.add_argument("--out-dir", default="~/mascope-output",
                    help="base output dir (a versioned run folder is created under it)")
    pb.add_argument("--ts", default=None,
                    help="cached full-batch TS parquet (else fetched live from the server)")
    pb.add_argument("--subject", default=None, help="optional subject phrase for the VK title")
    pb.add_argument("--no-report", action="store_true", help="skip the PDF report")
    pb.add_argument("--select", choices=["representative", "brightest"],
                    default="representative",
                    help="sample-selection strategy: 'representative' (5 time-spaced + "
                         "max-TIC) or 'brightest' (bin all peaks, assign each significant "
                         "m/z bin's brightest sample — better analyte coverage)")
    pb.add_argument("--coverage-target", type=float, default=0.85,
                    help="brightest: fraction of significant m/z bins to cover (default 0.85)")
    pb.add_argument("--k-max", type=int, default=10,
                    help="brightest: max number of winner samples to assign (default 10)")
    pb.add_argument("--height-floor", type=float, default=1000.0,
                    help="brightest: a bin is significant if its max height >= this (cps)")
    pb.set_defaults(func=cmd_batch)

    pr = sub.add_parser("report",
                        help="regenerate figures + PDF report from an existing run folder (offline)")
    pr.add_argument("--run-dir", required=True,
                    help="run folder holding merged_ledger.csv + per_file/")
    pr.add_argument("--reagent", required=True, help="Br | Ur | NO3 | ...")
    pr.add_argument("--ts", required=True, help="full-batch TS parquet")
    pr.add_argument("--batch", default=None, help="batch name for the report title")
    pr.add_argument("--tag", default=None, help="filename token (default: reagent name)")
    pr.add_argument("--run-id", default=None, help="Report ID (default: run-dir basename)")
    pr.add_argument("--generated", default=None, help="generated stamp for the cover")
    pr.add_argument("--subject", default=None)
    pr.set_defaults(func=cmd_report)

    pg = sub.add_parser("gka", help="interactive rotating-GKA HTML from a ledger CSV")
    pg.add_argument("ledger_csv")
    pg.add_argument("-o", "--out", default=None)
    pg.add_argument("--ppm", type=float, default=2.0, help="mass accuracy for band width")
    pg.set_defaults(func=cmd_gka)

    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.env:
        os.environ["MASCOPE_ENV"] = os.path.expanduser(args.env)
    return _run_guarded(lambda: args.func(args))


if __name__ == "__main__":
    sys.exit(main())

"""Van Krevelen + raw time-series plots for a run's analytes, computed
identically for any instrument via peaky.analyte_viz.

    python3 scripts/analyte_plots.py \
        --ledger <LEDGER.csv> --ts-parquet <BATCH_peaks.parquet> \
        --adducts '[M+Br]-,[M-H]-,[M+HBr+Br]-' \
        --out-prefix ~/peaky-output/<name>/<name> --label 'Br⁻ CIMS' --batch <SAMPLE_ID>

    --label is the reagent (e.g. 'Br⁻ CIMS', 'Ur⁺ CIMS'); --batch is the batch /
    sample id. Both are shown in the figure titles as "<label> · <batch> — ...".

Writes <prefix>_{vankrevelen.png, timeseries.png, viz.json}. The json is the
array payload the interactive chat widget consumes (vk + ts). RAW intensity by
default (`--mode raw`); `--mode reagent --reagent-mzs a,b,c` or `--mode tic`
only when a real normaliser exists.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from peaky import analyte_viz as V  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description="analyte Van Krevelen + time series")
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--ts-parquet", required=True)
    ap.add_argument("--adducts", required=True, help="comma-separated, e.g. '[M+Br]-,[M-H]-'")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--batch", default="", help="batch / sample id, shown in figure titles")
    ap.add_argument("--mode", default="raw", choices=["raw", "reagent", "tic"])
    ap.add_argument("--reagent-mzs", default=None, help="comma-separated m/z for --mode reagent")
    ap.add_argument("--top-ts", type=int, default=28)
    args = ap.parse_args(argv)

    adducts = [a.strip() for a in args.adducts.split(",") if a.strip()]
    reagent_mzs = ([float(x) for x in args.reagent_mzs.split(",")]
                   if args.reagent_mzs else None)
    ledger = pd.read_csv(args.ledger)
    ts = pd.read_parquet(args.ts_parquet)

    an = V.analyte_table(ledger)
    an = V.attach_dynamics(an, ts, adducts, mode=args.mode, reagent_mzs=reagent_mzs)
    grid, traces = V.time_traces(ts, an["neutral_formula"].tolist(), adducts,
                                 mode=args.mode, reagent_mzs=reagent_mzs)

    pre = Path(args.out_prefix).expanduser()
    pre.parent.mkdir(parents=True, exist_ok=True)
    lab = " · ".join(x for x in (args.label, args.batch) if x)
    V.render_van_krevelen(an, f"{pre}_vankrevelen.png",
                          title=f"{lab} analytes — Van Krevelen".strip())
    V.render_timeseries(grid, traces, an, f"{pre}_timeseries.png", top=args.top_ts,
                        title=f"{lab} changing analytes — raw log time series".strip())
    payload = V.widget_payload(an, grid, traces, top_ts=args.top_ts)
    Path(f"{pre}_viz.json").write_text(json.dumps(payload))

    nchg = int(an["changing"].sum())
    print(f"[{args.label or 'analytes'}] {len(an)} analytes ({nchg} changing, "
          f"{int((an.nN>0).sum())} CHON); wrote {pre}_{{vankrevelen.png,timeseries.png,viz.json}}")


if __name__ == "__main__":
    main()

"""Diel (time-of-day) composition figures for a batch run — class-resolved.

Turns a merged ledger + full-batch time series into two figures that show HOW
each compound class breathes over the day:

  * <prefix>_diel_composite.png   one line per class (median signal / class mean)
  * <prefix>_diel_individual.png  small-multiples: every ion's own diel profile
                                  (thin) + the class median (bold), one panel/class

The scientific point (measured on a nitrate-CIMS field dataset): a molecule's oxidation state
predicts its time of day — daytime photochemistry BUILDS oxygen (oxygenated CHO,
organonitrates, diacids peak ~13-15h) while nighttime EMISSION accumulates reduced
species (fatty-acid ladder peaks ~03-05h), and anthropogenic semivolatiles
(chloroparaffins) lag to a ~17h temperature-driven volatilization peak.

Everything is computed as a per-ion HOUR-OF-DAY profile (median within each local
hour, on log intensity) normalised to the ion's own mean, so shape is comparable
across ions of very different brightness. Classes come from the backbone formula
(peaky.batch.cluster.formula_class), so this works for any reagent/polarity.

Usage
-----
    # from a batch run folder (auto-finds merged_ledger.csv + data/*_ts.parquet)
    python3 scripts/diel_classes.py --run-dir <RUN_FOLDER> --label 'NO3 CIMS'

    # or explicitly
    python3 scripts/diel_classes.py \
        --ledger  <merged_ledger.csv> \
        --ts-parquet <BATCH_ts.parquet> \
        --out-prefix <dir>/<name> --label 'Br CIMS' --tz-offset -5

Key options: --tz-offset H (local = UTC+H; CDT = -5, default 0 = UTC) ·
--tier {assigned,all} (default assigned) · --min-detect (min finite points per
ion, default 80) · --classes a,b,c to restrict/order the panels.

Importable: `load_ion_diel()` / `class_of()` / `diel_profile()` are pure and can
be reused from a notebook.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from peaky import analyte_viz as V                     # noqa: E402
from peaky.batch.cluster import formula_class          # noqa: E402  (backbone classer)

# day/night shading + the sub-class panels the individual figure breaks out. Each
# panel is (title, predicate over the per-ion counts+class, colour). Predicates are
# evaluated per row so the set is easy to extend without touching the plot code.
DAY = (11, 17)
NIGHT_A, NIGHT_B = (0, 6), (21, 24)


def _counts(formula: str) -> dict:
    d: dict = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", str(formula)):
        if el:
            d[el] = d.get(el, 0) + (int(n) if n else 1)
    return d


def class_of(formula: str) -> str:
    """Backbone class of a neutral formula (reuses the report's classer)."""
    return formula_class(formula)


# richer, chemistry-aware panels for the small-multiples view (supersets of the
# coarse backbone classes; each takes the ions whose formula matches)
def _panel_members(m: pd.DataFrame):
    C, H, O = m["C"], m["H"], m["O"]
    halo = m["Cl"] + m["Br"]
    panels = [
        ("CHO oxygenated (O>=3)", (m.cls == "CHO") & (O >= 3), "#1D9E75"),
        ("CHON organonitrates", m.cls == "CHON", "#7F77DD"),
        ("fatty acids  CnH2nO2", (m.cls == "CHO") & (H == 2 * C) & (O == 2), "#0F6E56"),
        ("dicarboxylic  CnH(2n-2)O4", (m.cls == "CHO") & (O == 4) & (H == 2 * C - 2), "#378ADD"),
        ("chloroparaffins  CnClx", (halo >= 4) & (O == 0), "#D85A30"),
        ("fluorinated", m.cls == "F-containing", "#BA7517"),
        ("CHOS / CHONS organosulfur", m.cls.isin(["CHOS", "CHONS"]), "#BA1750"),
        ("Si / siloxane", m.cls == "Si / siloxane", "#888780"),
    ]
    return panels


def load_ion_diel(ledger, ts, *, tier="assigned", tz_offset=0.0, min_detect=80,
                  bin_minutes=None):
    """Return (ions_df, profiles) where ions_df is one row per ion channel with its
    class + element counts, and profiles[key] is a length-24 array = that ion's
    hour-of-day profile (local hour = UTC + tz_offset) normalised to its own mean.
    `ledger`/`ts` may be DataFrames or paths."""
    m = pd.read_csv(ledger) if isinstance(ledger, str) else ledger.copy()
    if tier == "assigned" and "tier" in m.columns:
        m = m[m["tier"] == "Assigned"].copy()
    m = m.dropna(subset=["neutral_formula", "mz", "adduct"]).copy()
    for el in ("C", "H", "O", "N", "S", "F", "Cl", "Br", "Si", "P"):
        m[el] = m["neutral_formula"].map(lambda f, e=el: _counts(f).get(e, 0))
    m["cls"] = m["neutral_formula"].map(class_of)
    m["key"] = m["neutral_formula"].astype(str) + "|" + m["adduct"].astype(str)
    m = m.drop_duplicates("key")

    ts = pd.read_parquet(os.path.expanduser(ts)) if isinstance(ts, str) else ts
    ts = ts.copy()
    ts["datetime_utc"] = pd.to_datetime(ts["datetime_utc"], utc=True)
    t0 = ts["datetime_utc"].min()
    grid, tr = V.ion_traces(ts, dict(zip(m["key"], m["mz"])), mode="raw",
                            bin_minutes=bin_minutes)
    # local hour-of-day for each trace point
    hod = ((np.asarray(grid, float) + t0.hour + t0.minute / 60.0 + tz_offset) % 24)
    hbin = np.floor(hod).astype(int)

    prof: dict = {}
    for k in m["key"]:
        if k not in tr:
            continue
        y = tr[k].to_numpy(float)
        y = np.where(y > 0, y, np.nan)
        if np.isfinite(y).sum() < min_detect:
            continue
        s = pd.Series(np.log10(y)).groupby(hbin).median().reindex(range(24))
        s = s.interpolate(limit_direction="both")
        lin = 10 ** s.to_numpy()
        mu = np.nanmean(lin)
        if mu > 0 and np.isfinite(mu):
            prof[k] = lin / mu
    m = m[m["key"].isin(prof)].reset_index(drop=True)
    return m, prof


def peak_hour(profile) -> float:
    """Circular-mean peak hour of a length-24 diel profile (amplitude-weighted)."""
    x = np.arange(24)
    w = np.clip(profile - np.nanmin(profile), 0, None)
    ang = (np.nansum(w * np.exp(1j * x * 2 * np.pi / 24)))
    return float((np.angle(ang) % (2 * np.pi)) * 24 / (2 * np.pi))


def diel_profile(m, prof, mask):
    """Class median profile over the ions selected by boolean `mask` (aligned to m)."""
    keys = m.loc[mask, "key"]
    P = [prof[k] for k in keys if k in prof]
    return (np.nanmedian(np.vstack(P), axis=0), len(P)) if P else (None, 0)


def _bands(ax):
    ax.axvspan(*DAY, color="gold", alpha=0.11)
    ax.axvspan(*NIGHT_A, color="navy", alpha=0.06)
    ax.axvspan(*NIGHT_B, color="navy", alpha=0.06)


def render(m, prof, out_prefix, *, label="", classes=None):
    """Write the composite + individual small-multiples figures. Returns their paths."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.arange(24)
    panels = _panel_members(m)
    if classes:
        want = [c.strip().lower() for c in classes]
        panels = [p for p in panels if any(w in p[0].lower() for w in want)]
    panels = [(t, mask, c) for (t, mask, c) in panels if mask.sum()]
    title_pfx = f"{label} · " if label else ""
    paths = []

    # ---- composite: one line per class ----------------------------------------
    figc, ax = plt.subplots(figsize=(9.8, 5.6), constrained_layout=True)
    _bands(ax)
    for t, mask, c in panels:
        med, n = diel_profile(m, prof, mask)
        if med is None:
            continue
        xs = np.r_[x, 24]
        ax.plot(xs, np.r_[med, med[0]], "-o", ms=3, lw=2, color=c, label=f"{t} (n={n})")
    ax.set_xlabel("hour of day (local)"); ax.set_ylabel("median signal / class mean")
    ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 3)); ax.grid(alpha=0.2)
    ax.legend(fontsize=8.5, loc="upper left")
    ax.set_title(f"{title_pfx}class-resolved diel composite", fontsize=12, weight="bold")
    pc = f"{out_prefix}_diel_composite.png"
    figc.savefig(pc, dpi=150, bbox_inches="tight"); plt.close(figc); paths.append(pc)

    # ---- individual small-multiples -------------------------------------------
    ncol = 3
    nrow = int(np.ceil(len(panels) / ncol))
    figi, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4.3 * nrow),
                              sharex=True, constrained_layout=True, squeeze=False)
    for ax, (t, mask, c) in zip(axes.ravel(), panels):
        _bands(ax)
        P = [prof[k] for k in m.loc[mask, "key"] if k in prof]
        for p in P:
            ax.plot(x, p, color=c, lw=0.6, alpha=0.28)
        if P:
            med = np.nanmedian(np.vstack(P), axis=0)
            ax.plot(x, med, color="#111", lw=2.4, zorder=5)
            ph = peak_hour(med)
            ax.axvline(ph, color="#111", ls=":", lw=1)
            ax.set_title(f"{t}   (n={len(P)}, peak ~{ph:.0f}h)", fontsize=10.5,
                         weight="bold", loc="left")
        ax.set_xlim(0, 23); ax.set_xticks(range(0, 24, 6)); ax.grid(alpha=0.18)
        ax.set_ylim(0.4, 2.2)
    for ax in axes.ravel()[len(panels):]:      # hide unused cells
        ax.set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("signal / own mean")
    for ax in axes[-1, :]:
        ax.set_xlabel("hour of day (local)")
    figi.suptitle(f"{title_pfx}individual ion diel profiles by class  "
                  "(thin = each ion, bold = class median; gold=day, blue=night)",
                  fontsize=12.5, weight="bold")
    pi = f"{out_prefix}_diel_individual.png"
    figi.savefig(pi, dpi=145, bbox_inches="tight"); plt.close(figi); paths.append(pi)
    return paths


def _resolve_run(run_dir):
    """(ledger, ts, out_prefix, label) from a batch run folder."""
    run_dir = os.path.expanduser(run_dir)
    ledger = os.path.join(run_dir, "merged_ledger.csv")
    tss = glob.glob(os.path.join(run_dir, "data", "*_ts.parquet"))
    if not os.path.exists(ledger):
        raise SystemExit(f"no merged_ledger.csv in {run_dir}")
    if not tss:
        raise SystemExit(f"no data/*_ts.parquet in {run_dir}")
    tag = os.path.basename(tss[0]).replace("_ts.parquet", "")
    figs = os.path.join(run_dir, "figures"); os.makedirs(figs, exist_ok=True)
    return ledger, tss[0], os.path.join(figs, f"story_{tag}"), tag


def main(argv=None):
    ap = argparse.ArgumentParser(description="class-resolved diel figures for a batch")
    ap.add_argument("--run-dir", help="batch run folder (auto-finds ledger + TS)")
    ap.add_argument("--ledger", help="merged_ledger.csv (if not --run-dir)")
    ap.add_argument("--ts-parquet", help="full-batch time series parquet")
    ap.add_argument("--out-prefix", help="output path prefix")
    ap.add_argument("--label", default="", help="reagent/title label, e.g. 'NO3 CIMS'")
    ap.add_argument("--tz-offset", type=float, default=0.0,
                    help="local = UTC + this many hours (CDT = -5; default 0 = UTC)")
    ap.add_argument("--tier", default="assigned", choices=["assigned", "all"])
    ap.add_argument("--min-detect", type=int, default=80,
                    help="min finite trace points for an ion to be plotted")
    ap.add_argument("--bin-minutes", type=int, default=None,
                    help="time-bin width (default: native per-sample)")
    ap.add_argument("--classes", default=None,
                    help="comma-separated substrings to restrict/order panels")
    args = ap.parse_args(argv)

    if args.run_dir:
        ledger, ts, out_prefix, tag = _resolve_run(args.run_dir)
        label = args.label or tag
        out_prefix = args.out_prefix or out_prefix
    else:
        if not (args.ledger and args.ts_parquet and args.out_prefix):
            raise SystemExit("need --run-dir, or all of --ledger/--ts-parquet/--out-prefix")
        ledger, ts, out_prefix, label = args.ledger, args.ts_parquet, args.out_prefix, args.label

    m, prof = load_ion_diel(ledger, ts, tier=args.tier, tz_offset=args.tz_offset,
                            min_detect=args.min_detect, bin_minutes=args.bin_minutes)
    print(f"[diel] {len(m)} ion channels with a diel profile "
          f"(tier={args.tier}, tz_offset={args.tz_offset:+g}h)")
    classes = args.classes.split(",") if args.classes else None
    paths = render(m, prof, out_prefix, label=label, classes=classes)
    for p in paths:
        print("wrote", p)


if __name__ == "__main__":
    main()

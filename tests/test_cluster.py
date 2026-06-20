"""Offline tests for cluster.py. Run: python3 tests/test_cluster.py"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import cluster as CL  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


# --- two anti-phase families: signed distance must keep them apart -----------
n = 40
t = np.linspace(0, 1, n)
rise = 10 ** (1 + 2 * t)          # rising
fall = 10 ** (3 - 2 * t)          # falling (anti-phase)
rng_noise = np.cos(np.arange(n))  # deterministic wiggle (no Math.random)
traces = {}
for i in range(4):
    traces[f"R{i}"] = rise * (1 + 0.02 * np.roll(rng_noise, i))
for i in range(4):
    traces[f"F{i}"] = fall * (1 + 0.02 * np.roll(rng_noise, i))
df = pd.DataFrame(traces)
cols = list(df.columns)

Lg, cm = CL.correlate(df, cols)
lab, big = CL.cluster(cm)
check("rising and falling families separate into 2 clusters", len(big) == 2, dict(lab))
rcl = {lab[c] for c in cols if c.startswith("R")}
fcl = {lab[c] for c in cols if c.startswith("F")}
check("all rising share one cluster", len(rcl) == 1, rcl)
check("all falling share one cluster", len(fcl) == 1, fcl)
check("rising and falling are DIFFERENT clusters (anti-phase not merged)",
      rcl.isdisjoint(fcl), (rcl, fcl))

# --- shape_of ----------------------------------------------------------------
check("shape_of rise", CL.shape_of(np.linspace(-2, 2, n)) == "rise")
check("shape_of fall", CL.shape_of(np.linspace(2, -2, n)) == "fall")
check("shape_of peak", CL.shape_of(np.concatenate([np.linspace(-1, 2, n // 2),
                                                   np.linspace(2, -1, n // 2)])) == "peak")

# --- threshold_scan returns sizes --------------------------------------------
ts = CL.threshold_scan(cm)
check("threshold_scan returns a dict keyed by cut", set(ts) >= {0.3, 0.4, 0.5})

# --- graceful on tiny / empty inputs -----------------------------------------
Lg0, cm0 = CL.correlate(df, [])
check("correlate([]) -> empty frames", len(cm0.columns) == 0)
lab0, big0 = CL.cluster(cm0)
check("cluster(empty) -> no clusters, no crash", big0 == [] and len(lab0) == 0)
lab1, big1 = CL.cluster(cm.iloc[:1, :1])     # single item
check("cluster(1 item) -> no >=3 cluster, no crash", big1 == [])
check("threshold_scan(empty) -> {}", CL.threshold_scan(cm0) == {})

# --- cluster_rows ------------------------------------------------------------
grid = t * 1.6
rows, Z = CL.cluster_rows(cols, lab, big, cm, df, grid)
check("cluster_rows builds one row per big cluster", len(rows) == len(big), len(rows))
check("cluster_rows carries (cid, members, rbar, shape, peak_hr)",
      all(len(r) == 5 for r in rows) and all(r[2] > 0.9 for r in rows),
      [(r[0], round(r[2], 2), r[3]) for r in rows])

# --- flatness gate: split varying vs flat, render the flat bunch --------------
# rise/fall traces vary strongly; add two genuinely flat traces (cv ~ 0).
flat_df = df.copy()
flat_df["FLAT0"] = 1000.0 + 0.0 * t                  # dead flat
flat_df["FLAT1"] = 500.0 * (1 + 0.001 * rng_noise)   # ~flat (tiny wiggle)
varying, flat = CL.split_varying(flat_df, list(flat_df.columns))
check("split_varying: rise/fall go to varying", set(cols) <= set(varying), varying)
check("split_varying: flat traces pulled out", set(flat) == {"FLAT0", "FLAT1"}, flat)
check("trace_cv: a varying trace clears the cv gate", CL.trace_cv(flat_df, "R0") >= CL.CHANGING)
check("trace_cv: a flat trace is below the gate", CL.trace_cv(flat_df, "FLAT0") < CL.CHANGING)
check("split_varying: a no-points trace is treated flat",
      "Z" in CL.split_varying(pd.DataFrame({"Z": [np.nan] * n}), ["Z"])[1])

import os, tempfile  # noqa: E402
with tempfile.TemporaryDirectory() as d:
    out = CL.render_flat_panel(flat, flat_df, grid, f"{d}/flat_p1.png", lambda c: str(c),
                               label="flat", title="flat test")
    check("render_flat_panel: writes one PNG", bool(out) and os.path.exists(out) and os.path.getsize(out) > 4000)
    check("render_flat_panel(empty) -> None", CL.render_flat_panel([], flat_df, grid, f"{d}/x.png", str) is None)

# --- per-cluster workbook (one tab per cluster) ------------------------------
wb_rows = [(1, ["A", "B", "C"], 0.9, "rise", 0.5),
           ("remaining (singletons)", ["D", "E"], float("nan"), "n/a", 0.0)]
wb_meta = {"A": {"neutral_formula": "C6H14O4", "channel": "+H⁺", "match_score": 0.98, "tier": "Identified"},
           "B": {"neutral_formula": "C6H14O4", "channel": "+Ur⁺", "match_score": 0.91, "tier": "Candidate"}}
with tempfile.TemporaryDirectory() as d:
    out = CL.write_cluster_workbook(wb_rows, f"{d}/wb.xlsx", meta=wb_meta,
                                    item_label=lambda k: f"ion-{k}",
                                    member_cols=["neutral_formula", "channel", "match_score", "tier"])
    check("write_cluster_workbook: file created", bool(out) and os.path.exists(out))
    import openpyxl  # noqa: E402
    names = openpyxl.load_workbook(out).sheetnames
    check("workbook: summary sheet + one sheet per cluster",
          names[0] == "summary" and "c1" in names and len(names) == 3, names)
    s1 = pd.read_excel(out, "c1")
    check("cluster sheet carries member + formula/channel/score/tier",
          "member" in s1.columns and "channel" in s1.columns and "tier" in s1.columns and len(s1) == 3,
          list(s1.columns))
    check("write_cluster_workbook([]) -> None", CL.write_cluster_workbook([], f"{d}/x.xlsx") is None)

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

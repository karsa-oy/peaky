"""Offline tests for analyte_viz.py. Run: python3 tests/test_analyte_viz.py"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import analyte_viz as V  # noqa: E402
from mascope_assign import chemistry as C    # noqa: E402

PASS = FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


# --- analyte_table: M0 organics, excludes Si contamination, VK coords ---
led = pd.DataFrame([
    dict(role="M0", neutral_formula="C9H19NO", adduct="[M+H]+", tier="Identified"),
    dict(role="M0", neutral_formula="C9H19NO", adduct="[M+Br]-", tier="Candidate"),  # SAME neutral, 2nd channel
    dict(role="M0", neutral_formula="C7H16O3", adduct="[M+H]+", tier="Identified"),  # CHO analyte
    dict(role="M0", neutral_formula="C8H26O5Si4", adduct="[M+H]+", tier="Candidate"),  # siloxane CONTAM
    dict(role="reagent", neutral_formula="", adduct="", tier=""),             # not an analyte
    dict(role="unexplained", neutral_formula=np.nan, adduct=np.nan, tier=""),
])
an = V.analyte_table(led)
check("analyte_table = one row per neutral (dedup cross-channel), drops Si + non-M0",
      len(an) == 2, an["neutral_formula"].tolist())
check("dedup keeps the Identified tier of a cross-channel neutral",
      an[an.neutral_formula == "C9H19NO"].iloc[0].tier == "Identified")
check("no Si-bearing analyte survives",
      not an["neutral_formula"].str.contains("Si").any())
r = an[an.neutral_formula == "C7H16O3"].iloc[0]
check("H/C and O/C on the neutral", abs(r.hc - 16/7) < 1e-6 and abs(r.oc - 3/7) < 1e-6,
      (r.hc, r.oc))
check("CHON vs CHO class",
      an[an.neutral_formula == "C9H19NO"].iloc[0].klass == "CHON"
      and r.klass == "CHO")

# --- time_traces + attach_dynamics: changing vs flat from a synthetic series ---
mzA = C.ion_mz("C9H19NO", "[M+H]+")     # make this one VARY -> high cv
mzB = C.ion_mz("C7H16O3", "[M+H]+")     # keep this one FLAT  -> low cv
rows = []
heights_A = [1e3, 1e3, 1e3, 5e4, 1e3, 1e3]   # one big spike
for i, hA in enumerate(heights_A):
    t = f"2025-10-02 0{i}:00:00"
    rows.append(dict(sample_item_id=f"s{i}", datetime_utc=t, mz=mzA, height=hA))
    rows.append(dict(sample_item_id=f"s{i}", datetime_utc=t, mz=mzB, height=2e4))
ts = pd.DataFrame(rows)
grid, traces = V.time_traces(ts, ["C9H19NO", "C7H16O3"], ["[M+H]+"], bin_minutes=60)
check("time_traces returns a trace per formula matched to its m/z",
      "C9H19NO" in traces and traces["C9H19NO"].notna().sum() >= 5,
      traces["C9H19NO"].tolist())
dyn = V.attach_dynamics(an, ts, ["[M+H]+"])
cvA = dyn[dyn.neutral_formula == "C9H19NO"].iloc[0].cv
cvB = dyn[dyn.neutral_formula == "C7H16O3"].iloc[0].cv
check("the spiking compound is 'changing' (high cv)",
      cvA >= V.CHANGING_CV and bool(dyn[dyn.neutral_formula == "C9H19NO"].iloc[0].changing), cvA)
check("the flat compound is NOT changing (low cv)",
      cvB < V.CHANGING_CV and not bool(dyn[dyn.neutral_formula == "C7H16O3"].iloc[0].changing), cvB)

# --- widget_payload structure ---
pl = V.widget_payload(dyn, grid, traces)
check("vk payload: one row per analyte, carries formula + channels for hover",
      len(pl["vk"]) == 2 and len(pl["vk"][0]) == 8
      and pl["vk"][0][5] and pl["vk"][0][6], pl["vk"])
check("ts payload lists only the changing analyte(s), with channel(s)",
      pl["ts"] is not None and [s["f"] for s in pl["ts"]["series"]] == ["C9H19NO"]
      and "ch" in pl["ts"]["series"][0], pl["ts"])

# --- a missing-formula trace degrades to NaN, not a crash ---
g2, tr2 = V.time_traces(ts, ["C40H80O2"], ["[M+H]+"])   # not in the synthetic ts
check("unmatched formula -> all-NaN trace (no crash)", tr2["C40H80O2"].isna().all())

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

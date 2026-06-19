"""Offline tests for timeseries.py. Run: python3 tests/test_timeseries.py"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import timeseries as TS  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


# --- synthetic 20-sample time series ---------------------------------------
N = 20
k = np.arange(N)
var = 1.0 + 0.8 * np.sin(2 * np.pi * k / N)      # variable (diel-like), cv ~0.57
flat = np.ones(N)                                  # flat / background
PEAKS = {
    236.7555: ("reagent", 1e6 * flat),             # Br3- reagent (normaliser)
    409.0011: ("dibromide", 5000 * flat),          # flat di-bromide cluster
    463.0000: ("contaminant", 8000 * flat),        # flat fluorinated inlet
    279.0236: ("mono1", 5000 * var),               # monoterpene anchor (variable)
    311.0134: ("mono2", 4000 * var),               # monoterpene anchor (variable)
    265.0080: ("ambient", 3000 * var),             # co-varies with monoterpenes
    124.9243: ("formic", 10000 * var),             # formic ref
}
rows = []
for mz, (_tag, tr) in PEAKS.items():
    for s in range(N):
        rows.append(dict(sample_item_id=f"s{s}", mz=mz, height=float(tr[s])))
peaks = pd.DataFrame(rows)

# --- matching ledger -------------------------------------------------------
led = pd.DataFrame([
    dict(peak_id="r1", mz=236.7555, height=1e6, role="reagent",
         neutral_formula=None, ion_formula="Br3-", adduct=None, tier=None, tier_reason=""),
    dict(peak_id="db", mz=409.0011, height=5000, role="M0",
         neutral_formula="C15H22O3", ion_formula="C15H23Br2O3-", adduct="[M+HBr+Br]-",
         tier="Identified", tier_reason="series"),
    dict(peak_id="ct", mz=463.0000, height=8000, role="M0",
         neutral_formula="C12H12F12", ion_formula="C12H12BrF12-", adduct="[M+Br]-",
         tier="Identified", tier_reason="unique"),
    dict(peak_id="m1", mz=279.0236, height=5000, role="M0",
         neutral_formula="C10H16O4", ion_formula="C10H16BrO4-", adduct="[M+Br]-",
         tier="Identified", tier_reason="iso"),
    dict(peak_id="m2", mz=311.0134, height=4000, role="M0",
         neutral_formula="C10H16O6", ion_formula="C10H16BrO6-", adduct="[M+Br]-",
         tier="Identified", tier_reason="iso"),
    dict(peak_id="am", mz=265.0080, height=3000, role="M0",
         neutral_formula="C9H14O4", ion_formula="C9H14BrO4-", adduct="[M+Br]-",
         tier="Candidate", tier_reason="ladder"),
    dict(peak_id="fo", mz=124.9243, height=10000, role="M0",
         neutral_formula="CH2O2", ion_formula="CH2BrO2-", adduct="[M+Br]-",
         tier="Identified", tier_reason="iso"),
])

summ = TS.apply_timeseries(led, peaks, log=lambda *a: None)

# --- assertions ------------------------------------------------------------
def disp(pid): return led.loc[led.peak_id == pid, "ts_disposition"].iloc[0]
def cv(pid): return led.loc[led.peak_id == pid, "ts_cv_norm"].iloc[0]
def tier(pid): return led.loc[led.peak_id == pid, "tier"].iloc[0]

check("matrix built + annotated all M0", summ["annotated"] == 6, summ)
check("di-bromide flat -> background disposition", disp("db").startswith("background:di-bromide"), disp("db"))
check("di-bromide flat -> DEMOTED Identified->Candidate", tier("db") == "Candidate", tier("db"))
check("demote count >=1", summ["demoted"] >= 1, summ)
check("fluorinated flat -> inlet contaminant", "inlet/instrument" in disp("ct"), disp("ct"))
check("contaminant NOT demoted (not di-bromide/CO3)", tier("ct") == "Identified", tier("ct"))
check("flat peaks have low cv_norm", cv("db") < 0.25 and cv("ct") < 0.25, (cv("db"), cv("ct")))
check("variable ambient peak high cv_norm", cv("am") > 0.4, cv("am"))
check("co-varying peak -> ambient disposition", disp("am").startswith("ambient"), disp("am"))
check("monoterpene anchor -> ambient", disp("m1").startswith("ambient"), disp("m1"))

# CO3-channel flat demotion
led2 = pd.DataFrame([
    dict(peak_id="r1", mz=236.7555, height=1e6, role="reagent", neutral_formula=None,
         ion_formula="Br3-", adduct=None, tier=None, tier_reason=""),
    dict(peak_id="co3", mz=361.0653, height=4000, role="M0", neutral_formula="C12H15NO8",
         ion_formula="C13H15NO11-", adduct="[M+CO3]-", tier="Identified", tier_reason="x"),
])
peaks2 = pd.DataFrame([dict(sample_item_id=f"s{s}", mz=mz, height=float(h[s]))
                       for mz, h in {236.7555: 1e6*flat, 361.0653: 4000*flat}.items() for s in range(N)])
TS.apply_timeseries(led2, peaks2, log=lambda *a: None)
check("flat CO3-channel -> demoted", led2.loc[led2.peak_id=="co3","tier"].iloc[0] == "Candidate",
      led2.loc[led2.peak_id=="co3","tier"].iloc[0])

# no reagent -> graceful (cv still computed on raw)
led3 = led.copy()
s3 = TS.apply_timeseries(led3, peaks, reagent_mzs=[], log=lambda *a: None)
check("runs without a reagent normaliser", s3["annotated"] == 6, s3)

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

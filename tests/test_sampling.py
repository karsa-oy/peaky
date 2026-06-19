"""Offline tests for sampling.py — the representative-sample rule.
Run: python3 tests/test_sampling.py"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import sampling as SS  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


def make_peaks(times, tics, *, sample_ids=None, peaks_per_sample=4):
    """Build a synthetic batch peak frame. `times` = list of datetimes, `tics` =
    per-sample total height (split across peaks_per_sample peaks)."""
    sample_ids = sample_ids or [f"s{i:02d}" for i in range(len(times))]
    t0 = pd.Timestamp("2025-10-01 21:00:00", tz="UTC")
    rows = []
    for sid, t, tic in zip(sample_ids, times, tics):
        for k in range(peaks_per_sample):
            rows.append(dict(sample_item_id=sid, sample_item_name=str(t),
                             datetime_utc=t, mz=100.0 + k, height=tic / peaks_per_sample))
    return pd.DataFrame(rows)


# --- regularly sampled 24-h run, evening TIC spike ---------------------------
t0 = pd.Timestamp("2025-10-01 21:00:00", tz="UTC")
times = [t0 + pd.Timedelta(minutes=72 * i) for i in range(20)]   # 20 samples / 24h
tics = np.ones(20) * 1e5
tics[15] = 9e5                                                    # evening spike
peaks = make_peaks(times, tics)

sel = SS.select_representative_samples(peaks)
ids = list(sel["sample_item_id"])
roles = dict(zip(sel["sample_item_id"], sel["role"]))

check("returns n_time grid + max-TIC = 6 rows", len(sel) == 6, f"got {len(sel)}")
check("time endpoints included (first+last sample)",
      "s00" in ids and "s19" in ids, ids)
check("max-TIC sample selected (s15, the evening spike)", "s15" in ids, ids)
check("max-TIC role flagged", roles.get("s15", "").endswith("max-TIC"), roles)
check("exactly 5 grid roles (one may also be max-TIC)",
      sum("time-grid" in r for r in roles.values()) == 5, roles)
check("rows are time-ordered",
      list(sel["datetime_utc"]) == sorted(sel["datetime_utc"]), ids)
# grid picks evenly spread: with 20 samples and 5 targets -> ~indices 0,5,9,14,19
check("grid spans the range (max gap between consecutive grid picks <= ~6)",
      max(np.diff(sorted(int(s[1:]) for s, r in roles.items() if r == "time-grid"))) <= 7,
      sorted(roles.items()))

# --- ids convenience matches the table --------------------------------------
check("select_representative_sample_ids == table ids",
      SS.select_representative_sample_ids(peaks) == ids)

# --- irregular sampling: dense early + one lone late file --------------------
dense = [t0 + pd.Timedelta(minutes=2 * i) for i in range(15)]    # 15 in first 30 min
late = [t0 + pd.Timedelta(hours=20)]                             # 1 file 20h later
itimes = dense + late
itics = list(np.ones(16) * 1e5)
ipeaks = make_peaks(itimes, itics, sample_ids=[f"d{i:02d}" for i in range(16)])
isel = SS.select_representative_samples(ipeaks)
iids = list(isel["sample_item_id"])
check("TIME-based (not index): lone late file is selected", "d15" in iids, iids)
check("irregular: still 5 distinct grid picks",
      sum("time-grid" in r for r in isel["role"]) == 5, list(isel["role"]))

# --- fewer than n_time samples -> all returned ------------------------------
few = make_peaks(times[:3], [1e5, 5e5, 2e5], sample_ids=["a", "b", "c"])
fsel = SS.select_representative_samples(few)
check("n < n_time -> all 3 returned", len(fsel) == 3, len(fsel))
check("n < n_time -> max-TIC (b) still flagged",
      fsel.set_index("sample_item_id").at["b", "role"].endswith("max-TIC"),
      fsel[["sample_item_id", "role"]].to_dict("records"))

# --- max-TIC coincides with a grid pick -> combined role --------------------
ctics = np.ones(20) * 1e5
ctics[0] = 9e5                                                    # first sample brightest
cpeaks = make_peaks(times, ctics)
csel = SS.select_representative_samples(cpeaks)
check("max-TIC == grid pick -> 'time-grid+max-TIC' role, 5 rows total",
      len(csel) == 5 and "time-grid+max-TIC" in set(csel["role"]),
      f"{len(csel)} rows, roles={list(csel['role'])}")

# --- include_max_tic=False, custom n_time -----------------------------------
nsel = SS.select_representative_samples(peaks, n_time=3, include_max_tic=False)
check("n_time=3, no max-TIC -> 3 rows, all time-grid",
      len(nsel) == 3 and set(nsel["role"]) == {"time-grid"},
      f"{len(nsel)} rows, roles={list(nsel['role'])}")

# --- sample_table basics ----------------------------------------------------
tab = SS.sample_table(peaks)
check("sample_table: one row per sample", len(tab) == 20, len(tab))
check("sample_table: tic = sum of heights", np.isclose(tab["tic"].max(), 9e5), tab["tic"].max())
check("sample_table: n_peaks counted", set(tab["n_peaks"]) == {4}, set(tab["n_peaks"]))

# --- no clock column -> falls back to all (no crash) ------------------------
noclock = peaks.drop(columns=["datetime_utc"])
ncsel = SS.select_representative_samples(noclock, n_time=5)
check("no datetime column -> returns all samples (no crash)", len(ncsel) == 20, len(ncsel))

# --- empty input ------------------------------------------------------------
empty = pd.DataFrame(columns=["sample_item_id", "datetime_utc", "height", "mz"])
esel = SS.select_representative_samples(empty)
check("empty peaks -> empty selection with role column",
      len(esel) == 0 and "role" in esel.columns, list(esel.columns))

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

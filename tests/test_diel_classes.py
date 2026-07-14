"""Offline smoke test for scripts/diel_classes.py — a tiny synthetic batch runs
through load_ion_diel + render and writes both figures. Run:
    python3 tests/test_diel_classes.py
"""
import os
import sys
import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import diel_classes as D  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


check("class_of routes backbones", [D.class_of(f) for f in
      ("C6H10O4", "C6H11NO6", "C10H16Cl6", "C3HF5", "C9H8OSi", "C8H10O5S")]
      == ["CHO", "CHON", "halogenated", "F-containing", "Si / siloxane", "CHOS"])
check("peak_hour of a noon-peaked profile is ~12",
      abs(D.peak_hour(1 + 0.5 * np.cos((np.arange(24) - 12) * 2 * np.pi / 24)) - 12) < 1.5)

# synthetic batch: a daytime CHO acid family + a nocturnal fatty-acid family, 4 days hourly
N = 96
t0 = datetime(2026, 6, 1, tzinfo=timezone.utc)
times = [t0 + timedelta(hours=i) for i in range(N)]
hod = np.array([(t.hour) for t in times])
day = 1 + 0.5 * np.cos((hod - 13) * 2 * np.pi / 24)      # peaks 13h UTC
night = 1 + 0.4 * np.cos((hod - 3) * 2 * np.pi / 24)     # peaks 03h UTC
chans = [("C5H8O4", "[M-H]-", 131.03, day, 3000), ("C6H10O4", "[M-H]-", 145.05, day, 2500),
         ("C7H12O4", "[M-H]-", 159.07, day, 2000), ("C8H16O2", "[M-H]-", 143.11, night, 4000),
         ("C10H20O2", "[M-H]-", 171.14, night, 3500), ("C12H24O2", "[M-H]-", 199.17, night, 3000)]
rows = []
for i, t in enumerate(times):
    for f, a, mz, shp, base in chans:
        rows.append((f"s{i:03d}", t, mz, base * shp[i]))
ts = pd.DataFrame(rows, columns=["sample_item_id", "datetime_utc", "mz", "height"])
led = pd.DataFrame([{"neutral_formula": f, "adduct": a, "mz": mz, "tier": "Assigned"}
                    for f, a, mz, _, _ in chans])

m, prof = D.load_ion_diel(led, ts, tz_offset=0.0, min_detect=24)
check("load_ion_diel returns a profile per ion", len(m) == 6 and len(prof) == 6, len(m))
check("each profile is length-24 and ~unit-mean", all(len(p) == 24 for p in prof.values())
      and all(abs(np.nanmean(p) - 1) < 0.05 for p in prof.values()))
# the daytime acids should peak in daytime, fatty acids at night
fa_mask = (m.cls == "CHO") & (m.H == 2 * m.C) & (m.O == 2)
ac_mask = (m.cls == "CHO") & (m.O == 4)
famed, fan = D.diel_profile(m, prof, fa_mask)
acmed, acn = D.diel_profile(m, prof, ac_mask)
check("fatty-acid family peaks at night (~3h)", fan == 3 and abs(D.peak_hour(famed) - 3) < 3,
      D.peak_hour(famed))
check("diacid family peaks in daytime (~13h)", acn == 3 and abs(D.peak_hour(acmed) - 13) < 3,
      D.peak_hour(acmed))

with tempfile.TemporaryDirectory() as d:
    paths = D.render(m, prof, os.path.join(d, "t"), label="TEST")
    check("render writes composite + individual PNGs",
          len(paths) == 2 and all(os.path.exists(p) for p in paths), paths)
    # run-dir resolution + full CLI path
    os.makedirs(os.path.join(d, "data")); os.makedirs(os.path.join(d, "figures"))
    led.to_csv(os.path.join(d, "merged_ledger.csv"), index=False)
    ts.to_parquet(os.path.join(d, "data", "SMOKE_ts.parquet"))
    D.main(["--run-dir", d, "--label", "smoke", "--min-detect", "24"])
    check("CLI --run-dir wrote figures",
          bool(list(Path(d, "figures").glob("*_diel_*.png"))))


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

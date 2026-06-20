"""Offline tests for series_detect.py. Run: python3 tests/test_series_detect.py"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import series_detect as SD  # noqa: E402
from mascope_assign import ledger as L  # noqa: E402

PASS = FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


rng = np.random.default_rng(7)


def ledger_with(mzs):
    peaks = pd.DataFrame({"peak_id": [f"p{i}" for i in range(len(mzs))],
                          "mz": mzs, "height": [1000.0] * len(mzs)})
    return L.new_ledger(peaks)


# planted CF2 ladder (8 members) + 60 random background peaks
cf2 = SD.UNIT_LIBRARY["CF2"][0]
ladder = [150.0123 + k * cf2 for k in range(8)]
ladder += [201.4567 + k * cf2 for k in range(6)]
background = list(rng.uniform(100, 700, 60))
led = ledger_with(ladder + background)
ev = SD.detect_series(led, ppm=5, min_links=8, min_enrichment=3.0)
cf2row = ev[ev.unit == "CF2"].iloc[0]
check("CF2 ladder detected with links >= 12", cf2row.n_links >= 12, cf2row.n_links)
check("CF2 chains3 detected", cf2row.n_chains3 >= 8, cf2row.n_chains3)
check("CF2 significant", bool(cf2row.significant), cf2row.to_dict())
check("CF2 enrichment >> decoys", cf2row.enrichment >= 3, cf2row.enrichment)
other = ev[(ev.unit != "CF2") & (ev.unit != "C2F4")]
check("no other unit significant on this data", not other.significant.any(),
      other[other.significant].unit.tolist())

# pure random background -> nothing significant
led2 = ledger_with(list(rng.uniform(100, 700, 80)))
ev2 = SD.detect_series(led2, ppm=5)
check("random background: nothing significant", not ev2.significant.any(),
      ev2[ev2.significant].unit.tolist())

# families_from_evidence maps CF2 -> fluorinated
fams = SD.families_from_evidence(ev)
check("CF2 evidence opens 'fluorinated'", "fluorinated" in fams, fams)
check("no spurious families", set(fams) <= {"fluorinated"}, fams)
nan_ev = pd.DataFrame([{"significant": True, "action": np.nan},
                       {"significant": True, "action": None},
                       {"significant": True, "action": "fluorinated"}])
check("NaN/None actions do not open fake families",
      SD.families_from_evidence(nan_ev) == ["fluorinated"],
      SD.families_from_evidence(nan_ev))

# low-confidence M0 peaks are included in the scanned population
led3 = ledger_with(ladder)
for pid in led3.peak_id[:4]:
    L.commit_assignment(led3, pid, neutral_formula="C5H8O2", adduct="[M-H]-",
                        ion_score=0.55, pass_no=3, method="x",
                        confidence="Suspect", commentary="weak")
ev3 = SD.detect_series(led3, ppm=5, min_links=8)
check("low-confidence peaks count toward links",
      ev3[ev3.unit == "CF2"].iloc[0].n_links >= 10,
      ev3[ev3.unit == "CF2"].iloc[0].n_links)

# --- unit_chains: ladder splits into chains with correct heads ---
chains = SD.unit_chains(led, cf2, ppm=5, min_len=3)
check("two CF2 chains found", len(chains) == 2, [len(c) for c in chains])
check("chain lengths 8 and 6", sorted(len(c) for c in chains) == [6, 8],
      sorted(len(c) for c in chains))
ch8 = max(chains, key=len)
check("chain head is lowest mz", abs(ch8[0][1] - 150.0123) < 1e-3, ch8[0])
check("chain steps are exactly CF2",
      all(abs(ch8[i + 1][1] - ch8[i][1] - cf2) < 1e-3 for i in range(len(ch8) - 1)))

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

"""Offline smoke test for clustering.cluster_batch — a tiny synthetic batch runs
end to end and writes the cluster / flat / unassigned CSVs + a figure without
crashing. (Faithful byte-equivalence vs the old scratch driver is checked
separately against real Orange data.) Run: python3 tests/test_clustering.py"""
import os
import sys
import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky import clustering as CLU, profiles as P  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


N = 14
t0 = datetime(2026, 1, 1, tzinfo=timezone.utc)
times = [t0 + timedelta(minutes=8 * i) for i in range(N)]      # ~1.7 h span
burst = np.array([1, 1, 1, 5, 9, 7, 4, 2, 1, 1, 1, 1, 1, 1], float)
flat = np.ones(N)
# (neutral_formula, adduct, mz, ion_score, tier, shape, base_cps)
channels = [
    ("C10H16O2", "[M+Br]-", 286.9, 0.95, "Identified", burst, 4000),
    ("C10H16O3", "[M+Br]-", 302.9, 0.92, "Identified", burst, 3000),
    ("C9H14O2",  "[M+Br]-", 272.9, 0.90, "Identified", burst, 2500),
    ("C6H14O4",  "[M+Br]-", 230.0, 0.88, "Candidate",  flat,  1500),
    ("C8H18O3",  "[M+Br]-", 244.0, 0.85, "Candidate",  flat,  1200),
    ("C8H26O5Si4", "[M+Br]-", 472.9, 0.80, "Candidate", flat,  800),   # Si contaminant
]
rows = []
for i, t in enumerate(times):
    sid = f"s{i:02d}"
    for (nf, ad, mz, sc, tier, shape, base) in channels:
        rows.append((sid, t, mz, base * shape[i]))
    rows.append((sid, t, 199.1234, 1000 * burst[i]))          # unassigned, co-varying
    rows.append((sid, t, 188.0000, 900.0))                    # unassigned, flat
ts = pd.DataFrame(rows, columns=["sample_item_id", "datetime_utc", "mz", "height"])
merged = pd.DataFrame([{"neutral_formula": nf, "adduct": ad, "mz": mz,
                        "ion_score": sc, "tier": tier}
                       for (nf, ad, mz, sc, tier, _, _) in channels])

with tempfile.TemporaryDirectory() as d:
    merged.to_csv(os.path.join(d, "merged_ledger.csv"), index=False)
    os.makedirs(os.path.join(d, "per_file"))
    pf = merged.copy(); pf["role"] = "M0"; pf["ion_formula"] = ""
    pf.to_csv(os.path.join(d, "per_file", "s00_ledger.csv"), index=False)

    res = CLU.cluster_batch(d, ts, P.resolve("Br"), tag="T", log=lambda *a: None)

    check("returns a dict with the expected keys",
          all(k in res for k in ("changing", "flat_clusters", "changers",
                                 "unassigned", "bin_minutes", "out_dir", "tag")))
    check("bin_minutes is None by default (native per-sample resolution)",
          res["bin_minutes"] is None, res["bin_minutes"])
    # an explicit bin width still works (legacy time-binning)
    res2 = CLU.cluster_batch(d, ts, P.resolve("Br"), tag="T2", bin_minutes=3, log=lambda *a: None)
    check("explicit bin_minutes is honored (legacy binning still available)",
          res2["bin_minutes"] == 3, res2["bin_minutes"])
    # tables -> tables/, figures -> figures/ (see paths.RunPaths)
    for fn in ("clusters_changing_T.csv", "clusters_flat_T.csv",
               "clusters_unassigned_T.csv", "channel_agreement_T.csv"):
        check(f"wrote tables/{fn}", os.path.exists(os.path.join(d, "tables", fn)))
    check("wrote a changing-cluster figure png under figures/",
          bool(list(Path(d, "figures").glob("clusters_changing_T_p*.png"))))
    check("wrote the per-cluster workbook under tables/",
          os.path.exists(os.path.join(d, "tables", "clusters_changing_T.xlsx")))

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

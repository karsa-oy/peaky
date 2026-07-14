"""Offline smoke test for clustering.cluster_batch — a tiny synthetic batch runs
end to end and writes the cluster / flat / unassigned CSVs + a figure without
crashing. (Faithful byte-equivalence vs the old scratch driver is checked
separately against real batch data.) Run: python3 tests/test_clustering.py"""
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
    ("C10H16O2", "[M+Br]-", 286.9, 0.95, "Assigned", burst, 4000),
    ("C10H16O3", "[M+Br]-", 302.9, 0.92, "Assigned", burst, 3000),
    ("C9H14O2",  "[M+Br]-", 272.9, 0.90, "Assigned", burst, 2500),
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
    # a ~1.7 h batch has no diel common mode to remove — corr_space must fall
    # back to raw (the panel median of a short run IS the event; removing it
    # would erase the main family)
    check("short batch falls back to raw correlation space",
          res["summary"].get("corr_space") == "raw", res["summary"].get("corr_space"))

# ---- LONG batch: residual space + unified unassigned union ------------------
N2 = 169                                                       # 7 days, hourly
times2 = [t0 + timedelta(hours=i) for i in range(N2)]
hrs = np.arange(N2, dtype=float)
day = (hrs // 24).astype(int)
daypat = np.array([0.0, 0.35, -0.25, 0.3, -0.35, 0.15, -0.2, 0.25])[day]   # family day-to-day
wave = 0.25 * np.sin(2 * np.pi * hrs / 24)                     # shared diel wave
fam = 10 ** (wave + 0.15 * daypat)                             # family shape
# a population of background channels with DIVERSE day-to-day patterns, so the
# panel-median common mode is not the family's own pattern (with a tiny panel the
# 3-member family would BE the median and residualization would erase it —
# documented small-panel trade-off; production panels have hundreds of channels)
bg_pats = {j: np.array([np.cos(1.3 * j + 2.7 * dd) for dd in range(8)])[day] for j in range(10)}
bg_mz = [230.0 + 3.7 * j for j in range(10)]
rows2 = []
for i, t in enumerate(times2):
    sid = f"L{i:03d}"
    rows2 += [(sid, t, 286.9, 4000 * fam[i]), (sid, t, 302.9, 3000 * fam[i]),
              (sid, t, 272.9, 2500 * fam[i])]                  # assigned family
    for j, bmz in enumerate(bg_mz):                            # assigned background
        rows2.append((sid, t, bmz, 1000 * 10 ** (wave[i] + 0.15 * bg_pats[j][i])))
    rows2.append((sid, t, 199.1234, 1500 * fam[i]))            # unknown, co-varies w/ family
    rows2.append((sid, t, 287.903355, 300 * fam[i]))           # 13C satellite of 286.9 (dimmer)
ts2 = pd.DataFrame(rows2, columns=["sample_item_id", "datetime_utc", "mz", "height"])
merged2 = pd.DataFrame(
    [{"neutral_formula": "C10H16O2", "adduct": "[M+Br]-", "mz": 286.9, "ion_score": .95, "tier": "Assigned"},
     {"neutral_formula": "C10H16O3", "adduct": "[M+Br]-", "mz": 302.9, "ion_score": .92, "tier": "Assigned"},
     {"neutral_formula": "C9H14O2", "adduct": "[M+Br]-", "mz": 272.9, "ion_score": .90, "tier": "Assigned"}]
    + [{"neutral_formula": f"C{8 + j}H{18 + 2 * j}O3", "adduct": "[M+Br]-", "mz": bmz,
        "ion_score": .85, "tier": "Candidate"} for j, bmz in enumerate(bg_mz)])
with tempfile.TemporaryDirectory() as d2:
    merged2.to_csv(os.path.join(d2, "merged_ledger.csv"), index=False)
    os.makedirs(os.path.join(d2, "per_file"))
    pf2 = merged2.copy(); pf2["role"] = "M0"; pf2["ion_formula"] = ""
    pf2.to_csv(os.path.join(d2, "per_file", "s00_ledger.csv"), index=False)
    res3 = CLU.cluster_batch(d2, ts2, P.resolve("Br"), tag="L", floor=100.0,
                             log=lambda *a: None)
    s3 = res3["summary"]
    check("long batch clusters in residual space", s3["corr_space"] == "residual")
    check("the 13C satellite bin is isotope-rejected from the union",
          s3["unassigned"]["n_isotope_rejected"] >= 1, s3["unassigned"])
    check("the co-varying unknown entered the union",
          s3["unassigned"]["n_entered_union"] >= 1, s3["unassigned"])
    cc = pd.read_csv(os.path.join(d2, "tables", "clusters_changing_L.csv"))
    unk = cc[cc.member_type == "unassigned"]
    check("unknown union member landed in a family with the assigned channels",
          len(unk) >= 1 and s3["unassigned"]["n_in_families"] >= 1, len(unk))
    check("unknown member's neutral_formula is a non-empty '?<mz>' (report-join safe)",
          bool(len(unk)) and unk["neutral_formula"].notna().all()
          and unk["neutral_formula"].astype(str).str.startswith("?").all(),
          list(unk.get("neutral_formula", [])))
    # the exact expression the PDF families() page runs must not raise
    for _cid, _g in cc.groupby("cluster"):
        ", ".join(_g.sort_values("median_cps", ascending=False)["neutral_formula"]
                  .dropna().astype(str).drop_duplicates().head(4))
    check("families()-style join over the changing CSV does not crash", True)
    check("family label gives unknowns chemical context",
          bool(len(unk)) and unk["cluster_label"].astype(str).str.contains("co-varies with").all(),
          list(unk.get("cluster_label", [])))
    # the rejected satellite must land in its own labeled bunch (not the flat-
    # unexplained page): a satellite tracks its bright parent's time series, so
    # it would otherwise sit there looking like an unexplained diel compound
    check("rejected satellite goes to the satellite bunch",
          s3["unassigned"]["n_isotope_satellites_bunched"] >= 1, s3["unassigned"])
    check("satellite panel page written",
          bool(list(Path(d2, "figures").glob("clusters_unassigned_L_p*.png"))))
    # gain vs family: fully-modulated members ~1; the column must be present and
    # finite for clustered members (flags baseline-dominated bright members)
    check("changing CSV carries finite gain_vs_family",
          "gain_vs_family" in cc.columns and cc["gain_vs_family"].notna().all(),
          list(cc.columns))

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

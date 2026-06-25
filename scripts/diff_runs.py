"""Diff two peaky run dirs' final assignments, peak-by-peak.

Validates the local-scoring pipeline against a committed backend run: same batch,
same samples, same peaks -> did the local pipeline assign the SAME neutral formula?

Keyed on (sample_item_id, peak_id) over the primary role=="M0" rows of each
merged_ledger.csv. Reports assignment agreement, coverage, tier shifts and the
score offset (local score_pattern is ~+0.07 vs the backend, per the parity eval).

    python scripts/diff_runs.py <backend_run_dir> <local_run_dir>
"""

from __future__ import annotations

import sys

import pandas as pd


def _m0(run_dir: str) -> pd.DataFrame:
    import glob

    frames = [
        pd.read_csv(p, low_memory=False)
        for p in glob.glob(f"{run_dir}/per_file/*_ledger.csv")
    ]
    df = pd.concat(frames, ignore_index=True)
    m0 = df[df["role"] == "M0"].copy()
    m0["key"] = m0["sample_item_id"].astype(str) + "|" + m0["peak_id"].astype(str)
    # one assignment per peak (defensive: keep the highest eff_score)
    m0 = m0.sort_values("eff_score", ascending=False).drop_duplicates("key")
    return m0.set_index("key")


def main(backend_dir: str, local_dir: str) -> None:
    b = _m0(backend_dir)
    l = _m0(local_dir)
    both = b.index.intersection(l.index)
    only_b = b.index.difference(l.index)
    only_l = l.index.difference(b.index)

    bf = b.loc[both, "neutral_formula"].fillna("")
    lf = l.loc[both, "neutral_formula"].fillna("")
    agree = bf == lf

    print(f"backend M0 peaks: {len(b)}   local M0 peaks: {len(l)}")
    print(f"assigned in both: {len(both)}")
    print(f"  agreement: {agree.mean():.3f}  ({agree.sum()}/{len(both)})")
    print(f"only backend assigned: {len(only_b)}   only local assigned: {len(only_l)}")

    # tier distribution shift
    print("\ntier counts (backend -> local):")
    bt = b["tier"].value_counts()
    lt = l["tier"].value_counts()
    for t in sorted(set(bt.index) | set(lt.index)):
        print(f"  {t:<10} {int(bt.get(t, 0)):>5} -> {int(lt.get(t, 0)):>5}")

    # score offset on agreed peaks
    bs = b.loc[both, "eff_score"].astype(float)
    ls = l.loc[both, "eff_score"].astype(float)
    d = (ls - bs)[agree.values]
    print(f"\neff_score offset (local - backend) on agreed peaks:")
    print(f"  mean {d.mean():+.4f}  median {d.median():+.4f}  std {d.std():.4f}")

    # show a sample of disagreements
    dis = both[~agree.values]
    if len(dis):
        print(f"\nsample disagreements (up to 15 of {len(dis)}):")
        cols = ["neutral_formula", "adduct", "eff_score", "tier"]
        for k in dis[:15]:
            bb = b.loc[k]
            ll = l.loc[k]
            print(
                f"  mz={bb['mz']:.4f}  backend {bb['neutral_formula']}/{bb['adduct']}"
                f" ({bb['eff_score']:.3f})  ->  local {ll['neutral_formula']}/{ll['adduct']}"
                f" ({ll['eff_score']:.3f})"
            )


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

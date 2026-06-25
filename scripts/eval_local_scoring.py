"""Parity eval: local mascope_tools scoring vs the backend `match_compounds`.

Runs over the per-file ledgers of a backend batch run (which carry both the sample
peaks and peaky's backend assignment) and reports:

  A. SCORE parity (isolated): score each backend-assigned M0 ion locally and compare
     the local ion score to the backend ion_score. Answers "is the score close enough".
  B. ASSIGNMENT agreement (end-to-end): local enumerate+score+argmax per peak vs the
     backend's assigned neutral formula. Answers "does it reach the same answer".

Usage:  python scripts/eval_local_scoring.py <run_dir>...  (Bromide and/or Uronium)
"""
from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from peaky.io.local_scoring import adduct_to_mech, score_candidates_local  # noqa: E402

# reagent -> (adducts, element ranges) from peaky/profiles.py
REAGENTS = {
    "Br": (["[M+Br]-", "[M-H]-", "[M+CO3]-", "[M+HBr+Br]-", "[M+HBr+CO3]-"],
           "C0-40 H0-80 N0-3 O0-18 S0-2 Cl0-2 Br0-2"),
    "Ur": (["[M+H]+", "[M+(CH4N2O)H]+", "[M+NH4]+"], "C0-40 H0-90 N0-8 O0-15 S0-2"),
}


def reagent_of(run_dir: str) -> str:
    return "Ur" if "uronium" in run_dir.lower() else "Br"


def eval_ledger(path: Path, adducts, ranges) -> dict:
    L = pd.read_csv(path)
    peaks = L[["mz", "height", "peak_id"]].dropna(subset=["mz"]).copy()
    m0 = L[(L["role"] == "M0") & L["neutral_formula"].notna()].copy()
    if m0.empty:
        return {}

    # ---- A. score parity: score the backend-assigned neutrals locally ----
    t = time.time()
    flat = score_candidates_local(peaks, sorted(m0["neutral_formula"].unique()), adducts)
    t_score = time.time() - t
    flat["mech_adduct"] = flat["mechanism_id"]
    base = flat[flat["is_base"] & flat["sample_peak_id"].notna()]
    # local score keyed by (neutral, mech, peak_id)
    local_score = {
        (r.compound_formula, r.mechanism_id, r.sample_peak_id): r.ion_score
        for r in base.itertuples()
    }
    pairs = []
    for r in m0.itertuples():
        mech = adduct_to_mech(r.adduct) if isinstance(r.adduct, str) else None
        ls = local_score.get((r.neutral_formula, mech, r.peak_id))
        if ls is not None and pd.notna(r.ion_score):
            pairs.append((float(r.ion_score), float(ls)))
    A = pd.DataFrame(pairs, columns=["backend", "local"]) if pairs else pd.DataFrame()

    # ---- B. assignment agreement: local argmax per backend M0 peak ----
    # best local M0 (any candidate neutral) per peak, restricted to backend M0 peaks
    bm0_ids = set(m0["peak_id"])
    bestloc = (base[base["sample_peak_id"].isin(bm0_ids)]
               .sort_values("ion_score", ascending=False)
               .drop_duplicates("sample_peak_id"))
    loc_formula = dict(zip(bestloc["sample_peak_id"], bestloc["compound_formula"]))
    agree = sum(1 for r in m0.itertuples()
                if loc_formula.get(r.peak_id) == r.neutral_formula)
    covered = sum(1 for r in m0.itertuples() if r.peak_id in loc_formula)

    return {
        "sample": path.stem.split("_")[0], "backend_M0": len(m0),
        "scored_pairs": len(A), "t_score": t_score,
        "score_mae": A["local"].sub(A["backend"]).abs().mean() if len(A) else np.nan,
        "score_corr": A["backend"].corr(A["local"]) if len(A) > 2 else np.nan,
        "within_0.1": (A["local"].sub(A["backend"]).abs() <= 0.1).mean() if len(A) else np.nan,
        "covered_frac": covered / len(m0),
        "agree_frac": agree / len(m0),
        "_A": A,
    }


def main(run_dirs: list[str]):
    rows, allA = [], []
    for rd in run_dirs:
        reagent = reagent_of(rd)
        adducts, ranges = REAGENTS[reagent]
        for led in sorted(Path(rd, "per_file").glob("*_ledger.csv")):
            r = eval_ledger(led, adducts, ranges)
            if r:
                r["reagent"] = reagent
                allA.append(r.pop("_A").assign(reagent=reagent))
                rows.append(r)
                print(f"[{reagent}] {r['sample']:18} M0={r['backend_M0']:4} "
                      f"pairs={r['scored_pairs']:4} MAE={r['score_mae']:.3f} "
                      f"corr={r['score_corr']:.3f} w0.1={r['within_0.1']:.2f} "
                      f"cover={r['covered_frac']:.2f} agree={r['agree_frac']:.2f} "
                      f"({r['t_score']:.1f}s)")
    df = pd.DataFrame(rows)
    A = pd.concat(allA, ignore_index=True)
    print("\n=== AGGREGATE ===")
    for reagent, g in df.groupby("reagent"):
        a = A[A["reagent"] == reagent]
        print(f"{reagent}: {len(g)} samples, {len(a)} scored ion-pairs | "
              f"score MAE={a['local'].sub(a['backend']).abs().mean():.3f} "
              f"corr={a['backend'].corr(a['local']):.3f} "
              f"within0.1={(a['local'].sub(a['backend']).abs()<=0.1).mean():.2f} | "
              f"cover={g['covered_frac'].mean():.2f} agree={g['agree_frac'].mean():.2f}")
    print(f"ALL: score MAE={A['local'].sub(A['backend']).abs().mean():.3f} "
          f"corr={A['backend'].corr(A['local']):.3f} "
          f"within0.1={(A['local'].sub(A['backend']).abs()<=0.1).mean():.2f}")
    return df, A


if __name__ == "__main__":
    main(sys.argv[1:])

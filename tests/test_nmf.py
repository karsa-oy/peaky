"""Offline tests for nmf.py (error-weighted NMF). Run: python3 tests/test_nmf.py"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky.batch import nmf as N  # noqa: E402

PASS = FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


def _synthetic(seed=1, n=200, m=60, k=3):
    """Two anti-phase diel factors + one flat factor, over ions whose brightness
    spans 5 orders of magnitude (the case that breaks unweighted NMF)."""
    rng = np.random.default_rng(seed)
    t = np.arange(n)
    Gt = np.stack([
        1 + 0.8 * np.sin(2 * np.pi * (t / 24 - 0.2)),   # midday
        1 + 0.8 * np.sin(2 * np.pi * (t / 24 - 0.7)),   # nocturnal (anti-phase)
        1 + 0.3 * rng.random(n),                        # flat-ish
    ], axis=1).clip(0, None)
    Ht = np.zeros((k, m))
    for f in range(k):
        Ht[f, rng.choice(m, 12, replace=False)] = rng.random(12)
    bright = 10 ** rng.uniform(1, 5, m)                  # 10 .. 100000 cps
    X = (Gt @ Ht) * bright
    X = (X * (1 + 0.05 * rng.standard_normal((n, m)))).clip(0, None)
    sigma = np.sqrt((N.REL_ERR * X) ** 2 + (N.FLOOR_FRAC * np.median(X[X > 0])) ** 2)
    return X, sigma, Gt, Ht


def _best_corr(G, Gt):
    return [max(abs(np.corrcoef(G[:, g], Gt[:, f])[0, 1]) for g in range(G.shape[1]))
            for f in range(Gt.shape[1])]


X, sigma, Gt, Ht = _synthetic()

# 1) weighted NMF recovers the two structured (anti-phase) factors
G, H, q, hist = N.wnmf(X, sigma, 3, seed=0, n_iter=800)
rec = _best_corr(G, Gt)
check("weighted NMF recovers midday factor (|corr|>0.95)", rec[0] > 0.95, rec)
check("weighted NMF recovers nocturnal factor (|corr|>0.95)", rec[1] > 0.95, rec)

# 2) weighting beats unweighting on brightness-skewed data (the whole point)
G2, H2, _, _ = N.wnmf(X, np.ones_like(X), 3, seed=0, n_iter=800)
rec_uw = _best_corr(G2, Gt)
check("error-weighting improves recovery over unweighted",
      min(rec[:2]) > max(rec_uw[:2]) + 0.05, f"weighted={rec[:2]} unweighted={rec_uw[:2]}")

# 3) non-negativity + normalisation convention (H rows sum to 1)
check("G, H non-negative", (G >= 0).all() and (H >= 0).all())
check("H rows sum to 1 (fractional profiles)", np.allclose(H.sum(axis=1), 1.0, atol=1e-6),
      H.sum(axis=1))

# 4) monotone Q history + rank scan Q decreases with k
check("Q decreases monotonically over iterations", all(np.diff(hist) <= 1e-6),
      hist)
rs = N.rank_scan(X, sigma, [2, 3, 4], n_iter=200)
check("rank_scan Q decreases with more factors", rs.Q.is_monotonic_decreasing, rs.Q.tolist())
check("rank_scan reports Q/Qexp", "Q_over_Qexp" in rs.columns and (rs.Q_over_Qexp > 0).all())

# 5) build_input from a tiny peaks frame: shape, non-detect->0, missing down-weighted
rows = []
for s in range(4):                       # 4 samples
    for mz, h in [(153.1274, 1000 * (s + 1)), (151.1118, 50)]:
        rows.append(dict(sample_item_id=f"s{s}", mz=mz, height=h,
                         datetime_utc=pd.Timestamp("2026-06-07", tz="UTC")))
# one sample is missing the dim ion -> a non-detection cell
peaks = pd.DataFrame([r for r in rows if not (r["sample_item_id"] == "s2" and r["mz"] == 151.1118)])
Xb, sb, bmz, tt, cc = N.build_input(peaks, tol_ppm=5.0)
check("build_input matrix shape = samples x bins", Xb.shape == (4, 2), Xb.shape)
check("non-detection cell filled with 0", (Xb == 0).sum() == 1, (Xb == 0).sum())
check("non-detection cell is down-weighted (max sigma)",
      sb.max() == sb[Xb == 0].max(), (sb.max(), sb[Xb == 0]))
check("build_input returns one datetime per sample", len(tt) == 4)

# 6) presence filter drops sparse bins
Xp, sp, bmzp, ttp, ccp = N.build_input(peaks, tol_ppm=5.0, min_presence=0.9)
check("min_presence drops the sparse (3/4) bin", Xp.shape[1] == 1, Xp.shape)

# 7) run_nmf packaging + factor helpers
res = N.run_nmf(peaks, 1, tol_ppm=5.0, n_iter=100)
check("run_nmf returns G,H,ts,top_ions", all(k in res for k in ("G", "H", "ts", "top_ions")))
check("factor timeseries has datetime + one col per factor", list(res["ts"].columns) == ["datetime_utc", "f0"],
      list(res["ts"].columns))


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

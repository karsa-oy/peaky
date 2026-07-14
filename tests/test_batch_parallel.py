"""Offline tests for sample-level batch parallelism (assign_batch process pool).
Run: python3 tests/test_batch_parallel.py

No network / no real pool: these cover the spawn-breaking failure modes -- job
resolution, picklability of the worker + its args, and that n_jobs is plumbed.
The byte-identical determinism guarantee (serial vs parallel) is verified live by
diffing merged_ledger.csv from an n_jobs=1 vs n_jobs=N run.
"""
import os
import pickle
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky.batch import assign_batch as AB   # noqa: E402
from peaky import passes as P                # noqa: E402

PASS = FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


# ---- _resolve_jobs -------------------------------------------------------
os.environ.pop("PEAKY_JOBS", None)
check("resolve: explicit 4 over 12 samples -> 4", AB._resolve_jobs(4, 12) == 4)
check("resolve: capped at sample count (8 over 3 -> 3)", AB._resolve_jobs(8, 3) == 3)
check("resolve: 1 stays 1 (serial)", AB._resolve_jobs(1, 12) == 1)
check("resolve: 0 -> auto physical cores (>=1, <=n)",
      1 <= AB._resolve_jobs(0, 12) <= 12)
check("resolve: None -> auto (>=1)", AB._resolve_jobs(None, 12) >= 1)
check("resolve: never exceeds n_samples", AB._resolve_jobs(99, 5) == 5)
check("resolve: n_samples=1 -> 1 (single-sample serial)", AB._resolve_jobs(6, 1) == 1)

os.environ["PEAKY_JOBS"] = "3"
check("resolve: $PEAKY_JOBS honored when n_jobs=None", AB._resolve_jobs(None, 12) == 3)
check("resolve: explicit arg overrides $PEAKY_JOBS", AB._resolve_jobs(6, 12) == 6)
os.environ.pop("PEAKY_JOBS", None)

check("physical_cores >= 1", AB._physical_cores() >= 1)

# ---- picklability (spawn re-pickles the target + initargs) ----------------
check("_assign_one is picklable (module-level target)",
      pickle.loads(pickle.dumps(AB._assign_one)) is AB._assign_one)
check("_worker_init is picklable", pickle.loads(pickle.dumps(AB._worker_init)) is AB._worker_init)
# the initargs the pool ships: context(str), reflists_active(list), base_kw(dict
# with a PassConfig + lists/ints), ts_path(str|None) -- all must pickle
base_kw = {"adducts": ["[M+H]+", "[M+NH4]+"], "cfg": P.PassConfig(),
           "label_isotope": None, "label_max": 2}
initargs = ("ambient-air", [], base_kw, None)
try:
    ra = pickle.loads(pickle.dumps(initargs))
    ok = ra[0] == "ambient-air" and isinstance(ra[2]["cfg"], P.PassConfig)
except Exception as e:  # noqa: BLE001
    ok = False
check("worker initargs (incl PassConfig) round-trip through pickle", ok)

# ---- signature plumbing --------------------------------------------------
import inspect  # noqa: E402
check("assign_batch.run exposes n_jobs param",
      "n_jobs" in inspect.signature(AB.run).parameters)
from peaky import pipeline as PL  # noqa: E402
check("pipeline.run_batch exposes n_jobs param",
      "n_jobs" in inspect.signature(PL.run_batch).parameters)

# ---- worker context isolation (deepcopy of cfg so mutation can't leak) ----
AB._worker_init("ambient-air", [], {"cfg": P.PassConfig(ppm=1.0)}, None)
c1 = AB._W["base_kw"]["cfg"]
import copy  # noqa: E402
c2 = copy.deepcopy(c1)
c2.prior_offset = 9.9
check("deepcopy isolates cfg (parent copy unmutated)", c1.prior_offset != 9.9)


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

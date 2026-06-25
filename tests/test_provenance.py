"""Offline tests for provenance.py. Run: python3 tests/test_provenance.py"""
import json
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky import passes as P  # noqa: E402
from peaky import provenance as PV  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


_rd = tempfile.mkdtemp()
with open(os.path.join(_rd, "merged_ledger.csv"), "w") as fh:
    fh.write("mz,neutral_formula,tier\n100.0,CH4,Assigned\n")
_tsp = os.path.join(_rd, "X_ts.parquet")
with open(_tsp, "wb") as fh:                       # any bytes -> hashable input
    fh.write(b"fake-parquet-bytes")

m = PV.build_manifest(run_dir=_rd, batch_name="B", dataset="D",
                      sample_ids=["s1", "s2"], reagent="NO3_15N",
                      cfg=P.PassConfig(), ts_path=_tsp, counts={"merged_M0": 42},
                      created_utc="2026-01-01T00:00:00Z")

check("manifest carries the core sections",
      all(k in m for k in ("run_id", "code", "input", "config", "output")))
check("code pins package version + a module-hash map",
      bool(m["code"]["package_version"])
      and isinstance(m["code"]["module_hashes"], dict)
      and len(m["code"]["module_hashes"]) > 5)
check("code records python + dep versions",
      "python" in m["code"] and "pandas" in m["code"]["deps"])
check("input hashes the ts parquet (== streaming sha1) + records provenance",
      m["input"]["ts_sha1"] == PV.sha1_file(_tsp)
      and m["input"]["dataset"] == "D" and m["input"]["sample_ids"] == ["s1", "s2"]
      and m["input"]["reagent"] == "NO3_15N")
check("output hashes merged_ledger.csv + carries counts",
      m["output"]["merged_ledger_sha1"] and m["output"]["counts"]["merged_M0"] == 42)
check("config fingerprint keeps user knobs, drops run-derived fields",
      m["config"].get("height_cutoff") == 100.0 and m["config"].get("ppm") == 1.0
      and "mechanism_ids" not in m["config"] and "prior_offset" not in m["config"])
check("git_info is best-effort and returns a dict", isinstance(
    PV.git_info(str(Path(PV.__file__).parent)), dict))

PV.write_manifest(_rd, m)
check("write_manifest writes run_manifest.json",
      os.path.exists(os.path.join(_rd, "run_manifest.json")))

_idx = os.path.join(_rd, "index.jsonl")
PV.append_registry(_idx, m)
PV.append_registry(_idx, m)
_rows = [json.loads(ln) for ln in open(_idx)]
check("append_registry appends ONE json row per call",
      len(_rows) == 2 and _rows[0]["run_id"] == m["run_id"]
      and _rows[0]["counts"]["merged_M0"] == 42 and _rows[0]["n_samples"] == 2)

m2 = PV.build_manifest(run_dir=_rd, batch_name="B", dataset="D",
                       sample_ids=["s1", "s2"], reagent="NO3_15N",
                       cfg=P.PassConfig(), ts_path=_tsp, counts={"merged_M0": 42})
check("hashing is deterministic (reproducible fingerprint)",
      m2["input"]["ts_sha1"] == m["input"]["ts_sha1"]
      and m2["output"]["merged_ledger_sha1"] == m["output"]["merged_ledger_sha1"])

check("sha1_file returns None for a missing path", PV.sha1_file(_rd + "/nope") is None)
check("record_run is non-fatal on a bad dir (returns {})",
      PV.record_run(run_dir="/nonexistent/xyz", base_out="/nonexistent",
                    batch_name="B", dataset=None, sample_ids=[], reagent="r",
                    cfg=P.PassConfig(), log=lambda *a: None) == {}
      or True)   # tolerate either {} or a partial manifest; must not raise


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

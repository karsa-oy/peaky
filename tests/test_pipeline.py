"""Offline tests for pipeline.py run-versioning helpers (deterministic via a
fixed datetime). Run: python3 tests/test_pipeline.py"""
import os
import sys
import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import pipeline as PL  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


WHEN = datetime(2026, 6, 20, 14, 35, 12)   # naive -> assumed UTC

check("slugify: spaces/punct -> single dashes, trimmed",
      PL.slugify("Orange peeling (Ur+ CIMS)") == "Orange-peeling-Ur-CIMS",
      PL.slugify("Orange peeling (Ur+ CIMS)"))
check("slugify: empty -> 'run'", PL.slugify("  ") == "run")

folder, human = PL.run_stamp(WHEN)
check("run_stamp: folder stamp is ISO-UTC to the second (…Z)", folder == "2026-06-20T143512Z", folder)
check("run_stamp: human stamp is date + HH:MM UTC", human == "2026-06-20 14:35 UTC", human)
# a tz-aware time in another zone is converted to UTC (the whole point of the fix)
aware = datetime(2026, 6, 20, 16, 35, 12, tzinfo=timezone(timedelta(hours=2)))   # = 14:35 UTC
check("run_stamp: tz-aware input converts to UTC", PL.run_stamp(aware)[0] == "2026-06-20T143512Z",
      PL.run_stamp(aware))

rid = PL.run_id("Orange peeling (Ur+ CIMS)", WHEN)
check("run_id: slug + ISO-UTC stamp", rid == "Orange-peeling-Ur-CIMS_2026-06-20T143512Z", rid)

with tempfile.TemporaryDirectory() as d:
    rd = PL.make_run_dir(d, "Orange peeling (Br- CIMS)", WHEN)
    check("make_run_dir: creates a fresh per-run folder", os.path.isdir(rd), rd)
    check("make_run_dir: folder basename == run_id (id locates the folder)",
          os.path.basename(rd) == PL.run_id("Orange peeling (Br- CIMS)", WHEN),
          os.path.basename(rd))
    check("make_run_dir: name carries the batch, the date AND the UTC time",
          "Orange-peeling-Br-CIMS" in os.path.basename(rd)
          and "2026-06-20" in os.path.basename(rd) and "143512Z" in os.path.basename(rd))

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

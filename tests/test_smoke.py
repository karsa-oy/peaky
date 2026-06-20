"""No-network 'install OK' smoke test — the fast first check a validator runs to
confirm the package imports, third-party deps resolve, and a tiny end-to-end
offline path works WITHOUT a Mascope account. Run: python3 tests/test_smoke.py
(or `pytest tests/test_smoke.py`). Should finish in a couple of seconds."""
import importlib
import json
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


# 1) public API surface + version
import mascope_assign as ma  # noqa: E402

check("package __version__ present", isinstance(ma.__version__, str) and bool(ma.__version__))
for name in ("run", "run_batch", "run_pipeline", "PassConfig", "get_context",
             "resolve_reagent", "ReagentProfile", "build_report"):
    check(f"public API exposes {name}", hasattr(ma, name))

# 2) every package submodule imports (NO network — io_mascope's SDK import is lazy)
SUBMODULES = [
    "chemistry", "contexts", "ledger", "io_mascope", "isotopes", "series_gka",
    "ladders", "series_detect", "reagents", "passes", "residual", "siloxane",
    "analyte_viz", "degeneracy", "cleanup", "sampling", "assign_batch", "cluster",
    "clustering", "composition", "plausibility", "pdf_report", "gka_figure",
    "profiles", "timeseries", "tiers", "report", "assign", "pipeline", "cli",
    "gka_widget", "__main__",
]
for m in SUBMODULES:
    try:
        importlib.import_module(f"mascope_assign.{m}")
        ok, detail = True, ""
    except Exception as e:                            # noqa: BLE001
        ok, detail = False, f"{type(e).__name__}: {e}"
    check(f"import mascope_assign.{m}", ok, detail)

# 3) third-party stack present (the deps not transitively guaranteed by mascope-sdk)
for lib in ("pandas", "numpy", "scipy", "matplotlib", "openpyxl", "dotenv"):
    try:
        importlib.import_module(lib); ok = True
    except Exception:                                 # noqa: BLE001
        ok = False
    check(f"dependency import: {lib}", ok)

# 4) a tiny PURE path: the match-tree fixture -> a flat per-isotopologue table
import pandas as pd  # noqa: E402

from mascope_assign import io_mascope as IO  # noqa: E402

tree = json.loads((Path(__file__).resolve().parent / "fixtures" / "match_tree.json").read_text())
flat = IO.flatten_match_tree(tree)
check("flatten_match_tree -> non-empty DataFrame",
      isinstance(flat, pd.DataFrame) and len(flat) > 0, len(flat) if hasattr(flat, "__len__") else "?")

# 5) reagent profile resolves offline
from mascope_assign import profiles as P  # noqa: E402

check("resolve_reagent('Br') works offline", P.resolve("Br").name == "Br")

# 6) openpyxl round-trips an xlsx (the report writer's dependency)
import openpyxl  # noqa: E402

with tempfile.TemporaryDirectory() as d:
    p = os.path.join(d, "t.xlsx")
    wb = openpyxl.Workbook(); wb.active["A1"] = "ok"; wb.save(p)
    check("xlsx write+reopen round-trips", openpyxl.load_workbook(p).active["A1"].value == "ok")


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

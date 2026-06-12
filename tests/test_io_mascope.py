"""Tests for io_mascope.py.

Offline: parser tests against tests/fixtures/match_tree.json (always run).
Live:    smoke test against Mascope, only when MASCOPE_LIVE=1.
Run: python3 tests/test_io_mascope.py   [MASCOPE_LIVE=1 for live]
"""
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import io_mascope as IO  # noqa: E402

PASS = FAIL = 0
HERE = Path(__file__).resolve().parent


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


# ---------- isotope label parsing ----------
check("base formula -> M0", IO.parse_isotope_label("C4H5O2-") == ("M0", True))
check("13C single", IO.parse_isotope_label("[13C]C3H5O2-") == ("13C", False))
check("13C2 multiple", IO.parse_isotope_label("[13C]2C2H5O2-") == ("13C2", False))
check("81Br", IO.parse_isotope_label("[81Br]CH2O2-") == ("81Br", False))
check("combined 13C+81Br",
      IO.parse_isotope_label("[13C]2[81Br]C2H6O2-") == ("13C2+81Br", False),
      IO.parse_isotope_label("[13C]2[81Br]C2H6O2-"))
check("18O", IO.parse_isotope_label("[18O]CH2BrO-") == ("18O", False))

# ---------- flatten_match_tree against the real fixture ----------
fix = HERE / "fixtures" / "match_tree.json"
if not fix.exists():
    check("fixture present", False, f"missing {fix}")
else:
    tree = json.load(open(fix))
    df = IO.flatten_match_tree(tree)
    check("flatten returns rows", len(df) > 0, len(df))
    cols = {"compound_formula", "compound_score", "ion_formula", "ion_score",
            "mechanism_id", "iso_label", "is_base", "theo_mz", "iso_score",
            "sample_peak_id", "ppm_error"}
    check("flatten has expected columns", cols <= set(df.columns),
          cols - set(df.columns))
    # exactly one base (M0) per ion
    per_ion = df.groupby(["compound_formula", "ion_formula"])["is_base"].sum()
    check("exactly one M0 base per ion", (per_ion == 1).all(), per_ion.to_dict())
    # C4H6O2 present with a base ion C4H5O2-
    sub = df[(df.compound_formula == "C4H6O2") & (df.is_base)]
    check("C4H6O2 base ion C4H5O2- present",
          "C4H5O2-" in set(sub.ion_formula), sorted(set(sub.ion_formula)))
    # base ion of a matched compound has a sample_peak_id and a finite iso_score
    base_matched = sub[sub.ion_formula == "C4H5O2-"]
    check("matched base ion has a sample_peak_id",
          base_matched["sample_peak_id"].notna().any(),
          base_matched[["sample_peak_id", "iso_score"]].to_dict("records"))
    check("base ion iso_score is high (>0.9)",
          (base_matched["iso_score"].fillna(0) > 0.9).any(),
          base_matched["iso_score"].tolist())
    # ppm error is defined only for genuinely matched isotopes, and is small
    real = df[df["ppm_error"].notna()]
    check("matched isotopes have a real sample_peak_id",
          real["sample_peak_id"].notna().all() and len(real) > 0, len(real))
    check("ppm errors |ppm|<10 for matched isotopes",
          (real["ppm_error"].abs() < 10).all(), real["ppm_error"].describe().to_dict())
    # brominated ion form present too ([M+Br]- => C4H6BrO2-)
    check("Br adduct ion form present",
          any("Br" in f for f in df[df.compound_formula == "C4H6O2"].ion_formula.dropna()))

# ---------- empty input safety ----------
import pandas as pd  # noqa: E402
check("flatten([]) -> empty df", len(IO.flatten_match_tree([])) == 0)

# ---------- live smoke (opt-in) ----------
if os.environ.get("MASCOPE_LIVE") == "1":
    print("\n-- live smoke --")
    SID = os.environ.get("MASCOPE_SID", "<sample-id>")
    cl = IO.connect()
    peaks = IO.fetch_peaks(cl, SID, use_cache=False)
    check("live: peaks fetched", peaks is not None and len(peaks) > 0, None if peaks is None else peaks.shape)
    mech = IO.resolve_mechanism_ids(cl, ["-H+", "+Br-"])
    check("live: mechanisms resolved", set(mech) == {"-H+", "+Br-"}, mech)
    scored = IO.score_candidates(cl, SID, ["C4H6O2", "C6H8O3"], mechanism_ids=None)
    check("live: score_candidates returns per-iso rows", len(scored) > 0, len(scored))
    check("live: at least one high ion_score",
          (scored["ion_score"].fillna(0) > 0.8).any())
else:
    print("\n(live smoke skipped; set MASCOPE_LIVE=1 to run)")

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

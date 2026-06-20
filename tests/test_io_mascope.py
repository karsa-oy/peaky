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


# ---------- score_candidates failure policy ----------
class _Matching:
    def match_compounds(self, sample_id, formulas, match_params, ionization_mechanism_ids):
        if "BAD" in formulas:
            raise RuntimeError("boom")
        return []


class _Client:
    matching = _Matching()


try:
    IO.score_candidates(_Client(), "SID", ["A", "BAD"], batch=1, workers=1)
    raised = False
except RuntimeError as e:
    raised = "match_compounds failed" in str(e) and "BAD" in str(e)
check("score_candidates raises on partial batch failure by default", raised)

partial = IO.score_candidates(_Client(), "SID", ["A", "BAD"], batch=1,
                              workers=1, allow_partial=True)
check("score_candidates allow_partial records failures",
      len(partial.attrs.get("match_batch_failures", [])) == 1,
      partial.attrs)

# ---------- estimate_offset: rough offset from the sample's own matches ----------
from mascope_assign import chemistry as _C  # noqa: E402
# build a synthetic match table at a uniform -1.9 ppm offset (Br-CIMS)
_off_rows = []
for f, mech in [("C10H16O4", "+Br-"), ("C10H16O3", "+Br-"), ("C5H10O3", "+Br-"),
                ("C2HF3O2", "+Br-"), ("C9H14O4", "-H+"), ("C6H12O3", "+Br-"),
                ("C4H6O4", "-H+"), ("HNO3", "+Br-"), ("C8H12O4", "+Br-"),
                ("C10H18O4", "+Br-")]:
    add = IO.MECH_TO_ADDUCT[mech]
    mz = _C.ion_mz(f, add) * (1 - 1.9e-6)            # observed = -1.9 ppm
    _off_rows.append({"target_compound_formula": f, "ionization_mechanism": mech,
                      "target_isotope_formula": f, "mz": mz})
# a heavy-isotope row + an unmatched row must be ignored
_off_rows.append({"target_compound_formula": "C10H16O4", "ionization_mechanism": "+Br-",
                  "target_isotope_formula": "[13C]C9H16BrO4-", "mz": 999.0})
_off_rows.append({"target_compound_formula": None, "ionization_mechanism": "+Br-",
                  "target_isotope_formula": None, "mz": 300.0})
off = IO.estimate_offset(pd.DataFrame(_off_rows))
check("estimate_offset recovers the -1.9 ppm instrument offset",
      off is not None and abs(off + 1.9) < 0.2, off)
check("estimate_offset None when too few matches",
      IO.estimate_offset(pd.DataFrame(_off_rows[:3])) is None)
check("estimate_offset None on a table with no match columns",
      IO.estimate_offset(pd.DataFrame({"mz": [1.0, 2.0]})) is None)

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

def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

"""Offline tests for profiles.py: built-in reagents (Br / Ur / NO3 / NO3_15N), alias
resolution, register(), and config-driven loading (JSON + TOML). The global
registry is snapshotted and restored so this file does not pollute other tests.
Run: python3 tests/test_profiles.py"""
import json
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky import profiles as P  # noqa: E402

PASS = FAIL = 0
def check(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  ok  {name}")
    else: FAIL += 1; print(f"FAIL  {name}  {detail}")


_SAVED = (dict(P.PROFILES), dict(P._BY_ALIAS))   # snapshot the registry

# ---- built-in reagents + alias resolution -----------------------------------
check("resolve('Br')", P.resolve("Br").name == "Br")
check("resolve('uronium' alias) -> Ur", P.resolve("uronium").name == "Ur")
check("resolve('NO3') built-in", P.resolve("NO3").name == "NO3")
check("NO3 is negative mode", P.resolve("nitrate").polarity == "-")
check("NO3 analyte channel is [M+NO3]-", "[M+NO3]-" in P.resolve("NO3").adducts)

# ---- 15N-labelled nitrate (distinct from the 14N NO3 profile) ----------------
check("resolve('NO3_15N') built-in", P.resolve("NO3_15N").name == "NO3_15N")
for _a in ("15no3", "^no3-", "nitrate-15n", "15n-nitrate"):
    check(f"alias {_a!r} -> NO3_15N", P.resolve(_a).name == "NO3_15N")
check("15N profile is negative mode", P.resolve("NO3_15N").polarity == "-")
check("15N analyte channel is [M+^NO3]-", "[M+^NO3]-" in P.resolve("NO3_15N").adducts)
check("15N keeps deprotonation channel", "[M-H]-" in P.resolve("NO3_15N").adducts)
check("15N detect_adduct distinguishes it from 14N",
      P.resolve("NO3_15N").detect_adduct == "[M+^NO3]-"
      and P.resolve("NO3").detect_adduct == "[M+NO3]-")
check("plain 'no3' still resolves to 14N NO3, not the 15N profile",
      P.resolve("no3").name == "NO3")
check("15N profile normalises on TIC (reagent ions out of window)",
      P.resolve("NO3_15N").normaliser == "tic")
try:
    P.resolve("xenon"); check("unknown reagent raises", False, "no raise")
except KeyError:
    check("unknown reagent raises KeyError", True)

# ---- register a new profile in code -----------------------------------------
acet = P.ReagentProfile(
    name="Ac", label="Acetate⁻", polarity="-", adducts=["[M+CH3COO]-", "[M-H]-"],
    normaliser="reagent", reagent_ion_re=r"C2H3O2-?$", ranges="C0-30 H0-50 O0-15",
    detect_adduct="[M+CH3COO]-", aliases=("acetate", "ac-"))
P.register(acet)
check("register() -> resolve by name", P.resolve("Ac").name == "Ac")
check("register() -> resolve by alias", P.resolve("acetate").name == "Ac")

# ---- config-driven loading (JSON + TOML) ------------------------------------
with tempfile.TemporaryDirectory() as d:
    cfgj = os.path.join(d, "r.json")
    json.dump([{"name": "Iod", "label": "I⁻ CIMS", "polarity": "-",
                "adducts": ["[M+I]-", "[M-H]-"], "normaliser": "reagent",
                "reagent_ion_re": "I-?$", "ranges": "C0-30 H0-50 O0-12",
                "detect_adduct": "[M+I]-", "aliases": ["iodide"]}], open(cfgj, "w"))
    P.load_config(cfgj)
    check("load_config(JSON list) registers", P.resolve("iodide").name == "Iod")
    check("loaded aliases are a tuple", isinstance(P.resolve("Iod").aliases, tuple))

    cfgw = os.path.join(d, "r2.json")
    json.dump({"reagents": [{"name": "Qx", "label": "Qx", "polarity": "-",
               "adducts": ["[M+Qx]-"], "normaliser": "tic", "reagent_ion_re": None,
               "ranges": "C0-10 H0-20", "detect_adduct": "[M+Qx]-"}]}, open(cfgw, "w"))
    check("resolve(config=) loads then resolves", P.resolve("Qx", config=cfgw).name == "Qx")

    cfgt = os.path.join(d, "r.toml")
    open(cfgt, "w").write(
        '[[reagents]]\nname="Tz"\nlabel="Tz"\npolarity="-"\nadducts=["[M+Tz]-"]\n'
        'normaliser="tic"\nreagent_ion_re=""\nranges="C0-5 H0-10"\ndetect_adduct="[M+Tz]-"\n')
    P.load_config(cfgt)
    check("load_config(TOML) registers", P.resolve("Tz").name == "Tz")

# ---- restore the registry (no cross-test pollution) -------------------------
P.PROFILES.clear(); P.PROFILES.update(_SAVED[0])
P._BY_ALIAS.clear(); P._BY_ALIAS.update(_SAVED[1])
check("registry restored (Ac gone after cleanup)", "Ac" not in P.PROFILES and "ac-" not in P._BY_ALIAS)
check("built-ins intact after restore", {"Br", "Ur", "NO3", "NO3_15N"} <= set(P.PROFILES))


def test_all():
    assert FAIL == 0, f"{FAIL} checks failed"


if __name__ == "__main__":
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)

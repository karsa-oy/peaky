"""Offline unit tests for contexts.py. Run: python3 tests/test_contexts.py"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mascope_assign import contexts as X  # noqa: E402

PASS = FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


# --- requested contexts all resolve, incl. aliases ---
for name in ("ambient-air", "chamber", "indoor-air", "object-headspace",
             "atmospheric", "smog-chamber", "indoor", "headspace"):
    try:
        X.get_context(name)
        check(f"context resolves: {name}", True)
    except Exception as e:
        check(f"context resolves: {name}", False, str(e))

def _raises(f):
    try:
        f(); return False
    except Exception:
        return True

check("unknown context raises", _raises(lambda: X.get_context("zzz")))


# --- structural gate composed in: half-integer ion formula rejected everywhere ---
keep, why = X.filter_by_context("C8H14NO12", "ambient-air")
check("half-integer neutral rejected by context filter", not keep, why)
keep, why = X.filter_by_context("C8H15NO12", "chamber")
# C8H15NO12: N/C = 0.125 ok, O/C = 1.5 ok for chamber(<=2.2), H/C=1.875, DBE integer
check("neutral nitrate C8H15NO12 plausible in chamber", keep, why)

# --- real atmospheric molecule passes ambient ---
keep, why = X.filter_by_context("C10H16O4", "ambient-air")
check("C10H16O4 plausible ambient", keep, why)

# --- PAH-like (low H/C) rejected in ambient but allowed in combustion ---
keep_a, _ = X.filter_by_context("C16H10", "ambient-air")   # H/C=0.625 < 0.7
keep_c, _ = X.filter_by_context("C16H10", "combustion")    # H/C ok, DBE/C ok
check("PAH C16H10 rejected ambient", not keep_a)
check("PAH C16H10 allowed combustion", keep_c)

# --- siloxane: rejected ambient (Si cap 1 + needs C scaffold), allowed indoor ---
# D4 siloxane octamethylcyclotetrasiloxane C8H24O4Si4
keep_amb, why_amb = X.filter_by_context("C8H24O4Si4", "ambient-air")
keep_ind, why_ind = X.filter_by_context("C8H24O4Si4", "indoor-air")
check("D4 siloxane rejected ambient (Si cap)", not keep_amb, why_amb)
check("D4 siloxane allowed indoor", keep_ind, why_ind)

# --- heteroatom-in-neutral guard: small Br formula is a reagent alias ---
keep, why = X.filter_by_context("C2H3O2Br", "ambient-air")  # C=2 < min_C_for Br=5
check("small-C brominated rejected as reagent alias", not keep, why)

# --- inorganic allowlist: H2SO4 ok in ambient, random C-free junk not ---
keep_so4, _ = X.filter_by_context("H2SO4", "ambient-air")
keep_junk, _ = X.filter_by_context("N4O2", "ambient-air")
check("H2SO4 allowed ambient (inorganic acid)", keep_so4)
check("N4O2 rejected ambient (not allowlisted)", not keep_junk)

# --- pass3 families wired per context ---
check("indoor opens siloxane+phthalate+glycol",
      {"siloxane", "phthalate", "glycol_peg"} <= set(X.get_context("indoor-air").pass3_families))
check("water opens halogen_dbp", "halogen_dbp" in X.get_context("water").pass3_families)
check("ambient opens organosulfate+nitrate",
      {"organosulfate", "nitrate"} <= set(X.get_context("ambient-air").pass3_families))

# --- contaminant family table sane ---
check("sulfate family adds S", X.CONTAMINANT_FAMILIES["organosulfate"]["add"].get("S") == (1, 1))
check("siloxane family adds Si", "Si" in X.CONTAMINANT_FAMILIES["siloxane"]["add"])
check("fluorinated family adds F", "F" in X.CONTAMINANT_FAMILIES["fluorinated"]["add"])

# --- classify_compound ---
cls, ox, tags = X.classify_compound("C10H16O4")
check("classify C10 monomer", cls == "C10 monomer", cls)
check("classify CHO only", tags == "CHO only", tags)
_, _, tags2 = X.classify_compound("C5H9NO4S")
check("classify organic-N + organosulfur", "organic-N" in tags2 and "organosulfur" in tags2, tags2)

# --- halogen-corrected H/C: halogens count as hydrogens in the ratio ---
# trichloroacetic acid C2HCl3O2: H/C=0.5 would fail, (H+X)/C=2.0 passes (water ctx)
keep_tca, why_tca = X.filter_by_context("C2HCl3O2", "water")
check("TCA C2HCl3O2 passes water ctx via (H+X)/C", keep_tca, why_tca)
# bromoform CHBr3 in water context (small-molecule branch, Heff=4 fine)
keep_bf, why_bf = X.filter_by_context("CHBr3", "water")
check("bromoform CHBr3 passes water ctx", keep_bf, why_bf)
# heavily brominated aromatic in water ctx: C6HBr5O (pentabromophenol-ish)
# (H+X)/C = 6/6 = 1.0 — inside bounds; plain H/C=0.17 would have failed
keep_pbp, why_pbp = X.filter_by_context("C6HBr5O", "water")
check("C6HBr5O rejected only by Br cap, not H/C", (not keep_pbp) and "Br=5" in str(why_pbp), why_pbp)

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

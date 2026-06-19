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

# --- saturated C3 polyols: glycerol/glycols are REAL (H/C 2.67); the h_to_c
#     ceiling was raised 2.6->2.75 so candidate generation no longer clips them
#     (they used to slip in only via the pass-4 iso-pair bypass). ---
keep_gly, why_gly = X.filter_by_context("C3H8O3", "ambient-air")   # glycerol, H/C 2.67
check("glycerol C3H8O3 plausible ambient (H/C 2.67)", keep_gly, why_gly)
keep_pg, _ = X.filter_by_context("C3H8O2", "ambient-air")          # propylene glycol
check("propylene glycol C3H8O2 plausible ambient", keep_pg)
keep_gly_in, _ = X.filter_by_context("C3H8O3", "indoor-air")
check("glycerol plausible indoor (glycols are real indoor signal)", keep_gly_in)
# DBE still bounds H/C: over-saturated junk stays rejected (these have DBE<0)
keep_j1, _ = X.filter_by_context("C4H22O4", "ambient-air")   # H/C 5.5, DBE -6
keep_j2, _ = X.filter_by_context("C3H10O3", "ambient-air")   # H/C 3.33, DBE -1
check("over-saturated C4H22O4 still rejected (DBE caps)", not keep_j1)
check("over-H C3H10O3 still rejected (DBE caps)", not keep_j2)

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

# --- POSITIVE uronium (urea-CIMS) context -----------------------------------
uro = X.get_context("uronium")
check("uronium resolves", uro.label == "uronium")
check("uronium aliases", X.get_context("urea-cims") is uro and X.get_context("urea") is uro)
check("uronium polarity positive", uro.polarity == "positive", uro.polarity)
check("negative contexts default polarity negative",
      X.get_context("ambient-air").polarity == "negative")
check("uronium widens the grid box (C46/O32)",
      uro.grid_c_max == 46 and uro.grid_o_max == 32, (uro.grid_c_max, uro.grid_o_max))
check("default grid box unchanged (C40/O30)",
      X.get_context("ambient-air").grid_c_max == 40
      and X.get_context("ambient-air").grid_o_max == 30)
check("uronium offers positive adducts incl. urea",
      "[M+(CH4N2O)H]+" in uro.reagent_adducts and "[M+H]+" in uro.reagent_adducts)
check("uronium bans neutral halogens (positive, no halide reagent)",
      uro.max_Br == 0 and uro.max_Cl == 0 and uro.max_F == 0)
check("uronium opens positive pass-3 families (amine/siloxane)",
      {"amine", "siloxane"} <= set(uro.pass3_families))
# an N-base analyte the urea source is selective for passes; a brominated
# neutral is rejected by the halogen cap (no Br reagent in positive mode)
keep_nbase, why_nbase = X.filter_by_context("C11H17NO2", "uronium")   # N/C 0.09, ok
check("N-base C11H17NO2 passes uronium", keep_nbase, why_nbase)
keep_amine, _ = X.filter_by_context("C6H15N", "uronium")              # tripropylamine-ish
check("amine C6H15N passes uronium", keep_amine)
keep_bromo, why_bromo = X.filter_by_context("C10H15BrO", "uronium")
check("brominated neutral rejected in uronium (Br cap 0)",
      (not keep_bromo) and "Br=" in str(why_bromo), why_bromo)
# a high-O HOM passes (urea-CIMS detects oxygenated VOC)
keep_hom, why_hom = X.filter_by_context("C10H16O6", "uronium")        # O/C 0.6
check("oxygenated VOC C10H16O6 passes uronium", keep_hom, why_hom)

# --- PDMS / long-siloxane family (assign the characterized siloxane ladder) ---
check("uronium opens the pdms family", "pdms" in uro.pass3_families)
check("uronium max_Si raised to 12 (reach the long PDMS ladder)", uro.max_Si == 12)
check("pdms family reaches Si>6 (the rungs siloxane can't)",
      X.CONTAMINANT_FAMILIES["pdms"]["add"]["Si"][1] >= 10)
check("pdms family offers positive Si adducts incl. NH4",
      "[M+NH4]+" in X.CONTAMINANT_FAMILIES["pdms"]["adducts"]
      and "[M+H]+" in X.CONTAMINANT_FAMILIES["pdms"]["adducts"])
# a high-Si PDMS oligomer passes the uronium VK filter (Heff/(C+Si)~2, O/(C+Si)~0.3)
keep_pdms, why_pdms = X.filter_by_context("C12H39NO6Si6", "uronium")  # 462.145 [M+H]+
check("PDMS oligomer C12H39NO6Si6 passes uronium", keep_pdms, why_pdms)
# the short siloxane family is UNCHANGED (no impact on the shared reference contexts)
check("shared siloxane family Si cap unchanged at 6",
      X.CONTAMINANT_FAMILIES["siloxane"]["add"]["Si"][1] == 6)

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

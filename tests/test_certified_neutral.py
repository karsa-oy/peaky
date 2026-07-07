"""Offline tests for certified_neutral.py (pure module -- no client, no network).
Run: python3 tests/test_certified_neutral.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peaky.assignment import certified_neutral as CN  # noqa: E402
from peaky.chem import chemistry as C  # noqa: E402
from peaky.chem import contexts as XC  # noqa: E402

PASS = FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


URO = XC.get_context("uronium")
ADDUCTS = ["[M+H]+", "[M+(CH4N2O)H]+", "[M+NH4]+"]

# ---------- channel_offsets ----------
offs = CN.channel_offsets(ADDUCTS, "urea", max_cluster_order=2)
keys = {(a, n) for a, n, _ in offs}
check("offsets include bare [M+H]+", ("[M+H]+", 0) in keys)
# the proton+2-urea rung (121.0720) survives as an OFFSET -- labeled as the
# registered urea adduct at order 1 after alias de-dup, not as ([M+H]+, 2)
check("offsets include the 2-urea proton rung (121.072, any label)",
      any(abs(off - 121.0720) < 1e-3 for _, _, off in offs))
# the alias: [M+H]+ @ order1 == [M+(CH4N2O)H]+ @ order0 -> exactly one survives
alias = [(a, n) for a, n, off in offs if abs(off - 61.0396392121) < 1e-3]
check("urea alias de-duplicated to ONE channel", len(alias) == 1, alias)
check("alias keeps the registered adduct form (order 0)",
      alias[0] == ("[M+(CH4N2O)H]+", 0), alias)
# no reagent -> no synthetic ladder
offs_nr = CN.channel_offsets(ADDUCTS, None)
check("no reagent -> order-0 only", all(n == 0 for _, n, _ in offs_nr))

# ---------- find_certificates: the NBBS ladder (real session ground truth) ----
nbbs = pd.DataFrame({
    "peak_id": ["p214", "p274", "p334", "noise1", "noise2"],
    "mz": [214.0896, 274.1220, 334.1544, 251.9070, 401.3337],
    "height": [71626.0, 496505.0, 3057.0, 500.0, 400.0],
})
certs = CN.find_certificates(nbbs, offs)
big = [c for c in certs if c.n_channels >= 3]
check("NBBS ladder certifies ONE 3-channel core", len(big) == 1, [c.core_mass for c in certs])
if big:
    c0 = big[0]
    check("NBBS core mass = 213.0823 (+-1 mDa)", abs(c0.core_mass - 213.0823) < 1e-3, c0.core_mass)
    check("NBBS spread tight (<0.5 mDa)", c0.spread_mda < 0.5, c0.spread_mda)
    check("NBBS members are the 3 ladder peaks",
          {h.peak_id for h in c0.hits} == {"p214", "p274", "p334"})
    check("cluster orders are 0/1/2 (after alias de-dup)",
          sorted(h.cluster_order for h in c0.hits) in ([0, 0, 1], [0, 1, 2]),
          [(h.adduct, h.cluster_order) for h in c0.hits])

# ---------- malathion cross-channel (H+ / urea / NH4) ----------
mal_mz = [C.ion_mz("C10H19O6PS2", a) for a in ADDUCTS]
mal = pd.DataFrame({"peak_id": [f"m{i}" for i in range(3)],
                    "mz": mal_mz, "height": [350128.0, 20000.0, 9000.0]})
mcerts = [c for c in CN.find_certificates(mal, offs) if c.n_channels >= 2]
best = max(mcerts, key=lambda c: c.n_channels) if mcerts else None
check("malathion certifies across 3 channels", best is not None and best.n_channels == 3,
      [(c.core_mass, c.n_channels) for c in mcerts])
if best:
    check("malathion core = 330.0360 (+-1 mDa)", abs(best.core_mass - 330.0360) < 1e-3, best.core_mass)

# a lone peak certifies nothing
lone = pd.DataFrame({"peak_id": ["x"], "mz": [331.0433], "height": [1.0]})
check("single peak -> no certificate", not CN.find_certificates(lone, offs))

# same peak cannot supply two channels to one core
dup = pd.DataFrame({"peak_id": ["a", "a2"], "mz": [331.0433, 331.0433], "height": [10.0, 9.0]})
dcerts = CN.find_certificates(dup, offs)
check("two peaks at the SAME m/z -> no multi-channel certificate",
      all(c.n_channels < 2 for c in dcerts), [(c.core_mass, c.n_channels) for c in dcerts])

# wide pair (outside tol) does not certify
wide = pd.DataFrame({"peak_id": ["w1", "w2"],
                     "mz": [331.0433, 331.0433 + 60.0324 + 0.008],  # 8 mDa off the rung
                     "height": [10.0, 10.0]})
check("8-mDa-off ladder step -> no certificate",
      all(c.n_channels < 2 for c in CN.find_certificates(wide, offs)))

# ---------- enumerate_certified: off-grid P licensed by a certificate ----------
forms = CN.enumerate_certified(330.0360, URO, force=True)
check("certified 330.0360 enumerates malathion C10H19O6PS2 (P+S2 open)",
      "C10H19O6PS2" in forms, forms[:8])
forms_noforce = CN.enumerate_certified(330.0360, URO, force=False,
                                       extra_elements={"P": (0, 0), "S": (0, 0), "Cl": (0, 0)})
check("P/S-closed box does NOT reach malathion", "C10H19O6PS2" not in forms_noforce)
forms_nbbs = CN.enumerate_certified(213.0823, URO, force=True)
check("certified 213.0823 enumerates NBBS C10H15NO2S", "C10H15NO2S" in forms_nbbs, forms_nbbs[:8])

# ---------- corroboration_count ----------
check("3 channels + iso = 3 corroborations",
      CN.corroboration_count(3, iso_confirmed=True) == 3)
check("2 channels alone = 1", CN.corroboration_count(2) == 1)
check("TS layers add", CN.corroboration_count(2, ts_covaries=True,
                                              order_lift_consistent=True) == 3)

# ---------- determinism ----------
c1 = CN.find_certificates(nbbs, offs)
c2 = CN.find_certificates(nbbs.sample(frac=1.0, random_state=7).reset_index(drop=True), offs)
check("certificates independent of input row order",
      [(round(a.core_mass, 6), a.n_channels) for a in c1]
      == [(round(b.core_mass, 6), b.n_channels) for b in c2])

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)

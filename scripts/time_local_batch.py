"""Measure local-scoring wall time per representative sample (enumerate + score)."""
import sys, time, glob, warnings
import pandas as pd
warnings.filterwarnings("ignore")
from mascope_tools.composition import CompositionSearchConfig, HeuristicFilterConfig
from mascope_tools.composition.finder import find_compositions
from peaky.io.local_scoring import score_candidates_local

CFG = {
 "Br": (CompositionSearchConfig(ionizations="+Br-,-H-",
        element_count_ranges="C0-40 H0-80 N0-3 O0-18 S0-2 Cl0-2 Br0-2",
        mass_range_ppm=5.0, use_unsaturation=True, only_integer_unsaturation=True),
        ["[M+Br]-","[M-H]-","[M+CO3]-","[M+HBr+Br]-","[M+HBr+CO3]-"]),
 "Ur": (CompositionSearchConfig(ionizations="+H+,+NH4+,+(CH4N2O)H+",
        element_count_ranges="C0-40 H0-90 N0-8 O0-15 S0-2",
        mass_range_ppm=5.0, use_unsaturation=True, only_integer_unsaturation=True),
        ["[M+H]+","[M+(CH4N2O)H]+","[M+NH4]+"]),
}
def reagent(d): return "Ur" if "uronium" in d.lower() else "Br"

for d in sys.argv[1:]:
    tag = reagent(d); cfg, adducts = CFG[tag]
    leds = sorted(glob.glob(d+"/per_file/*_ledger.csv"))
    tot_e=tot_s=0.0; n=0
    print(f"\n### {tag} batch ({len(leds)} representative samples available) ###")
    for f in leds:
        L = pd.read_csv(f); peaks = L[["mz","height","peak_id"]].dropna(subset=["mz"])
        t=time.time(); neutrals=set()
        for mz in peaks[peaks.height>=100].mz: neutrals|={r["formula"] for r in find_compositions(float(mz),cfg)}
        te=time.time()-t
        t=time.time(); flat=score_candidates_local(peaks,sorted(neutrals),adducts); ts=time.time()-t
        tot_e+=te; tot_s+=ts; n+=1
        print(f"  {f.split(chr(92))[-1].split('_')[0]:18} peaks={len(peaks):4} cand={len(neutrals):5} enum={te:5.1f}s score={ts:4.1f}s")
    per=(tot_e+tot_s)/n
    print(f"  per-sample avg: {per:.1f}s  | 6-sample subset: {per*6/60:.1f} min")

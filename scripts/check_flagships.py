"""Flagship-assignment regression check for a finished run's ledger.

Usage:  python3 scripts/check_flagships.py LEDGER.csv

Asserts that the assignments we have validated by hand (isotope physics,
cross-channel corroboration, GKA series membership) are present and sane in
the given ledger, and that the known junk classes are absent. Run this after
every pipeline change: TFA silently vanishing between v12 and v19 is the
failure mode this protects against.
"""
from __future__ import annotations

import re
import sys

import pandas as pd

# (neutral, adduct, max |ppm|, why it must be there)
FLAGSHIPS = [
    ("C2HF3O2",  "[M+Br]-", 1.5, "TFA: 0.978-ratio Br doublet at 192.91, Mascope 0.96"),
    ("C3H6O3",   "[M+Br]-", 1.5, "hydroxy-acid ladder C3 rung, 10.3k cps, twin-satellite confirmed"),
    ("C2H4O3",   "[M+Br]-", 1.5, "hydroxy-acid ladder C2 rung"),
    ("C5H10O3",  "[M+Br]-", 1.5, "hydroxy-acid ladder C5 rung"),
    ("C6H12O3",  "[M+Br]-", 1.5, "hydroxy-acid ladder C6 rung"),
    ("C10H16O5", "[M+Br]-", 1.5, "monoterpene ladder rung; 81Br twin at 297.016"),
    ("C10H16O4", None,      1.5, "monoterpene ladder rung, both channels"),
    ("C10H16O6", None,      1.5, "monoterpene ladder rung"),
    # cross-channel partners that the missing-13C audit falsely cleared in
    # v18-v21 (absent satellite = picker loss; agreeing second channel wins)
    ("C10H16O6", "[M+Br]-", 1.5, "Br partner of Good [M-H]-, 2.4k cps"),
    ("C9H14O3",  "[M+Br]-", 1.5, "Br partner of Good [M-H]-, 1.8k cps"),
    ("C10H18O4", "[M+Br]-", 1.5, "Br partner of Good [M-H]-, 2.0k cps"),
    ("C10H16O3", "[M+Br]-", 1.5, "Br partner of Good [M-H]-, 1.5k cps"),
    ("C2HF3O2",  "[M-H]-",  1.5, "TFA second channel (pass-5 completion)"),
    # silanediol/PDMS ladder (GKA-discovered 2026-06-12; Mascope-verified with
    # 29Si/30Si/81Br satellites; overturned C5H10O6@244.97 and C7H16O7Si@318.99)
    ("C2H8O2Si1",  "[M+Br]-", 1.5, "dimethylsilanediol n=1, 12k cps"),
    ("C4H14O3Si2", "[M+Br]-", 1.5, "silanediol oligomer n=2 (was bogus C5H10O6 'High', z=3.6)"),
    ("C6H20O4Si3", "[M+Br]-", 1.5, "silanediol oligomer n=3 (was bogus C7H16O7Si, z=3.0)"),
    ("C8H26O5Si4", "[M+Br]-", 1.5, "silanediol oligomer n=4: the 20k cps ex-#1 unknown"),
    ("C10H32O6Si5","[M+Br]-", 1.5, "silanediol oligomer n=5"),
]

# junk classes that must NEVER reappear (formula regexes on M0 neutrals)
JUNK = [
    (r"F\d*", lambda n: n.get("F", 0) >= 1 and n.get("O", 0) > 6,
     "fluorochemical with O>6 (v16/v19 flood class)"),
    (r"C5H10O6", lambda n: n == {"C": 5, "H": 10, "O": 6},
     "the CHO fantasy that mis-claimed silanediol n=2 at 244.9668"),
]


def parse(f: str) -> dict:
    d: dict = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", str(f)):
        if el:
            d[el] = d.get(el, 0) + (int(n) if n else 1)
    return d


def main() -> int:
    led = pd.read_csv(sys.argv[1])
    m0 = led[led["role"] == "M0"]
    failed = 0

    for nf, ad, max_ppm, why in FLAGSHIPS:
        rows = m0[m0["neutral_formula"] == nf]
        if ad is not None:
            rows = rows[rows["adduct"] == ad]
        if not len(rows):
            print(f"FAIL  {nf} {ad or '(any adduct)'} MISSING -- {why}")
            failed += 1
            continue
        r = rows.iloc[0]
        ppm = abs(r["ppm_error"]) if pd.notna(r["ppm_error"]) else 99.0
        conf = str(r["confidence"])
        if ppm > max_ppm:
            print(f"FAIL  {nf} present but |ppm|={ppm:.2f} > {max_ppm} -- {why}")
            failed += 1
        elif conf.startswith(("Suspect", "Reject")):
            print(f"FAIL  {nf} present but confidence '{conf}' -- {why}")
            failed += 1
        else:
            print(f"  ok  {nf:10s} {r['adduct']:9s} {conf} ({ppm:.2f} ppm)")

    n_junk = 0
    for _pat, pred, why in JUNK:
        for f in m0["neutral_formula"].dropna():
            if pred(parse(f)):
                print(f"FAIL  junk class present: {f} -- {why}")
                n_junk += 1
    failed += n_junk

    n = len(m0)
    print(f"\n{len(FLAGSHIPS) - failed + n_junk}/{len(FLAGSHIPS)} flagships ok, "
          f"{n_junk} junk hits, {n} M0 total")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

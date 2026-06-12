"""Reagent-ion library + labeler.

In chemical-ionization MS the reagent anion (Br-, I-, NO3-, ...) forms bright
cluster ions that are NOT sample chemistry: bare R_n clusters and R.(neutral)
clusters (water, HNO3, small acids). These otherwise land in 'unexplained' and
dominate the residual by signal. This module enumerates the cluster m/z (with
halogen isotopologue combinations) and labels matching ledger peaks as reagent.

It is keyed on the reagent detected from the sample's adducts, so a Br-CIMS run
gets the Br_n / Br.(acid) library and an I-CIMS run gets the I_n library.
"""
from __future__ import annotations

import itertools

import pandas as pd

from . import chemistry as C
from . import ledger as L

__version__ = "0.1.0"

# isotope masses for the reagent halogens (light, heavy, heavy abundance)
_HALOGEN_ISO = {
    "Br": [(78.9183371, "79Br"), (80.9162906, "81Br")],
    "Cl": [(34.96885268, "35Cl"), (36.96590259, "37Cl")],
    "I":  [(126.9044719, "127I")],
}
_M_E = C.M_E

# small neutrals that cluster onto the reagent anion (incl. organic acids that
# commonly form Br_n.acid clusters in Br-CIMS)
_CLUSTER_NEUTRALS = {
    "H2O": "H2O", "HNO3": "HNO3", "HBr": "HBr", "HNO2": "HNO2",
    "HCOOH": "CH2O2", "CH3COOH": "C2H4O2", "CO2": "CO2",
    "SO2": "O2S", "H2SO4": "H2O4S",
    "propionic": "C3H6O2", "pyruvic": "C3H4O3", "glyoxal": "C2H2O2",
    "oxalic": "C2H2O4", "MSA": "CH4O3S", "pinic": "C9H14O4",
}


def build_library(reagent: str = "Br", *, max_n: int = 3, max_neutral: int = 1
                  ) -> list[tuple[str, float]]:
    """Return [(label, ion_mz)] for negative-mode reagent clusters:
      * bare R_n^-  (odd n carries the charge: R-, R3-, R5-)
      * R^- . (neutral)_k  for the small neutral list
    All halogen isotopologue combinations are enumerated.
    """
    if reagent not in _HALOGEN_ISO:
        return []
    isos = _HALOGEN_ISO[reagent]
    out: list[tuple[str, float]] = []

    # bare R_n clusters (charge -1): odd n carries the charge (R-, R3-, R5-)
    core_masses: list[tuple[str, float, int]] = []   # (label, mass, n)
    for n in range(1, max_n * 2, 2):
        for combo in itertools.combinations_with_replacement(range(len(isos)), n):
            mass = sum(isos[i][0] for i in combo) + _M_E   # anion: +1 electron
            tag = "+".join(isos[i][1] for i in combo)
            label = f"[{reagent}{n}]- ({tag})"
            out.append((label, mass))
            core_masses.append((label, mass, n))

    # R_n^- . (neutral)_k clusters -- neutrals adduct onto each bare core
    for label, core_mz, n in core_masses:
        for name, formula in _CLUSTER_NEUTRALS.items():
            nm = C.neutral_mass(formula)
            for k in range(1, max_neutral + 1):
                out.append((f"[{reagent}{n}+{k}x{name}]-", core_mz + k * nm))

    # reagent-halogen oxide anions: RO-, RO2-, RO3-  (e.g. BrO-, BrO2-, BrO3-)
    light = isos[0][0]
    for no in (1, 2, 3):
        out.append((f"[{reagent}O{no if no > 1 else ''}]-",
                    light + no * C.M["O"] + _M_E))
    return out


def label_reagents(ledger: pd.DataFrame, reagent: str = "Br", *, ppm: float = 15.0,
                   only_unexplained: bool = True) -> int:
    """Mark ledger peaks matching a reagent-cluster m/z as role='reagent'.
    Returns the number of peaks labeled."""
    lib = build_library(reagent)
    if not lib:
        return 0
    lib_sorted = sorted(lib, key=lambda x: x[1])
    masses = [m for _, m in lib_sorted]
    import bisect
    n = 0
    for i, row in ledger.iterrows():
        if only_unexplained and row["role"] != L.ROLE_UNEXPLAINED:
            continue
        mz = row["mz"]
        tol = mz * ppm * 1e-6
        lo = bisect.bisect_left(masses, mz - tol)
        hi = bisect.bisect_right(masses, mz + tol)
        if hi > lo:
            # nearest label
            best = min(lib_sorted[lo:hi], key=lambda x: abs(x[1] - mz))
            try:
                L.mark_reagent(ledger, row["peak_id"],
                               f"reagent ion: {best[0]} ({(mz-best[1])/best[1]*1e6:+.1f} ppm)")
                n += 1
            except L.LedgerError:
                continue
    return n


def reagent_for_adducts(adducts: list[str]) -> str | None:
    """Pick the reagent halogen implied by the sample's adducts."""
    for a in adducts:
        if "Br" in a:
            return "Br"
        if a.endswith("I]-") or "+I" in a:
            return "I"
        if "Cl" in a:
            return "Cl"
    return None

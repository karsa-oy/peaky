"""Reagent profiles — ONE config per reagent system, replacing the adducts /
element-ranges / normaliser / label constants that were copy-pasted inline across
the time-series, clustering and validation scripts.

A profile is everything the pipeline needs to treat a batch's reagent correctly.
New reagent = add a ReagentProfile, not edit code. `resolve()` picks one by name
or auto-detects from a loaded peak table (polarity + the server's own adduct
mechanisms via io_mascope.detect_adducts).
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class ReagentProfile:
    name: str                       # short key, e.g. "Br" / "Ur"
    label: str                      # display label for figures, e.g. "Br⁻ CIMS"
    polarity: str                   # "-" or "+"
    adducts: list[str]              # analyte channels (mascope_assign adduct labels)
    normaliser: str                 # "reagent" | "tic"  (for the TS/correlation layer)
    reagent_ion_re: str | None      # regex on ion_formula picking the reagent ions
    ranges: str                     # grid element ranges for local enumeration
    detect_adduct: str | None       # presence of this adduct => this reagent (auto-detect)
    context: str = "ambient-air"    # default assign.run context (mode + VK priors + caps)
    aliases: tuple = field(default_factory=tuple)


BR = ReagentProfile(
    name="Br", label="Br⁻ CIMS", polarity="-",
    adducts=["[M+Br]-", "[M-H]-", "[M+HBr+Br]-"],
    normaliser="reagent", reagent_ion_re=r"Br\d-$",
    ranges="C0-40 H0-80 N0-3 O0-18 S0-2 Cl0-2 Br0-2",
    detect_adduct="[M+Br]-", context="ambient-air",
    aliases=("br", "bromide", "br-cims", "br-"))

UR = ReagentProfile(
    name="Ur", label="Ur⁺ CIMS", polarity="+",
    adducts=["[M+H]+", "[M+(CH4N2O)H]+"],
    normaliser="tic", reagent_ion_re=None,
    ranges="C0-40 H0-90 N0-8 O0-15 S0-2",
    detect_adduct="[M+(CH4N2O)H]+", context="uronium",
    aliases=("ur", "uronium", "urea", "urea-cims", "ur+"))

PROFILES: dict[str, ReagentProfile] = {BR.name: BR, UR.name: UR}
_BY_ALIAS = {a: p for p in PROFILES.values() for a in (p.name.lower(), *p.aliases)}


def resolve(reagent: str = "auto", peaks=None) -> ReagentProfile:
    """Return a ReagentProfile. `reagent` may be a name/alias, or 'auto' to detect
    from a loaded peak table (its server adduct mechanisms, then polarity)."""
    if reagent and reagent.lower() in _BY_ALIAS:
        return _BY_ALIAS[reagent.lower()]
    if reagent != "auto":
        raise KeyError(f"unknown reagent {reagent!r}; known: {sorted(_BY_ALIAS)}")
    if peaks is None:
        raise ValueError("reagent='auto' needs a peaks table to detect from")
    from . import io_mascope as IO
    seen = set(IO.detect_adducts(peaks))
    for p in PROFILES.values():
        if p.detect_adduct in seen:
            return p
    # fall back on polarity if no diagnostic adduct matched
    pol = _detect_polarity(peaks)
    for p in PROFILES.values():
        if p.polarity == pol:
            return p
    raise ValueError(f"could not auto-detect reagent (adducts={seen}, polarity={pol})")


def _detect_polarity(peaks) -> str | None:
    for col in ("polarity", "sample_batch_name", "ionization_mechanism"):
        if col in getattr(peaks, "columns", []):
            s = " ".join(map(str, peaks[col].dropna().unique()[:20]))
            if "+" in s and "-" not in s:
                return "+"
            if "-" in s and "+" not in s:
                return "-"
    return None

"""Chemical-plausibility QC for assigned formulas.

High-resolution accurate mass alone lets the fitter hit almost any target once it
is allowed many heteroatoms (every extra N/O/S/halogen adds free parameters), so a
small fraction of assignments are mass-coincidence "monsters" rather than real
molecules in an organic-aerosol matrix: extreme heteroatom counts, implausibly
carbon-rich (very low H/C) skeletons, or a covalent halogen in a positive-mode
spectrum whose reagent provides no halogen.

`scan` flags these for SCRUTINY — it does not delete or re-assign anything. It is
deliberately conservative and only inspects CANDIDATE-tier neutrals: an Assigned
assignment cleared the server's isotope-pattern score, which is independent
corroboration we don't second-guess from element ratios. A neutral seen as
Assigned in ANY channel is therefore never flagged.

Pure formula arithmetic; deterministic. Thresholds are intentionally loose so the
flagged set is small and defensible (the clear coincidences), not a dragnet.
"""
from __future__ import annotations

from . import chemistry as C

# thresholds (loose on purpose — flag the clear coincidences only)
N_HIGH_OC = 3       # N>=3 combined with...
OC_HIGH = 1.0       # ...O/C >= this  -> high-heteroatom coincidence
N_VERY_HIGH = 4     # N>=4 combined with...
O_HIGH = 8          # ...O>= this
HC_FLOOR = 0.35     # H/C below this -> implausibly carbon-rich (F-FREE formulas only)
F_HIGH = 4          # F>=this -> heavily fluorinated; F is monoisotopic, so the fit
                    # has NO isotope twin to confirm it (a fluorine mass coincidence).
                    # Flagged with the RIGHT reason: low H/C here is F displacing H
                    # (fluorine-rich), NOT a carbon-rich skeleton.


def implausible(neutral_formula: str, *, tier: str | None = None,
                polarity: str | None = None) -> str | None:
    """Return a short reason string if `neutral_formula` looks like a mass-coincidence
    fit rather than a real molecule, else None. Only Candidate-tier is scrutinised
    (pass tier=None to scrutinise regardless). `polarity` ('+'/'-') enables the
    wrong-mode-halogen check."""
    if tier is not None and str(tier) != "Candidate":
        return None
    c = C.parse_formula(str(neutral_formula))
    nc = c.get("C", 0)
    if nc == 0:
        return None                      # carbon-free handled elsewhere (reagent/inorganic)
    h, n, o = c.get("H", 0), c.get("N", 0), c.get("O", 0)
    f = c.get("F", 0)
    br, cl = c.get("Br", 0), c.get("Cl", 0)
    hc, oc = h / nc, o / nc
    # Terse labels (the full meaning is spelled out in the scrutiny-page legend);
    # keeping them short stops the table overflowing the page width.
    if n >= N_HIGH_OC and oc >= OC_HIGH:
        return f"N{n}, O/C {oc:.1f} (heteroatom coincidence)"
    if n >= N_VERY_HIGH and o >= O_HIGH:
        return f"N{n}O{o} (heteroatom coincidence)"
    if f >= F_HIGH:           # heavily fluorinated: 19F is 100% monoisotopic
        # NB any 13C/81Br satellites the row carries confirm the CARBON count / the
        # adduct halogen, NOT the fluorine -- 19F has no heavier stable isotope, so
        # the F COUNT is never isotope-confirmable (do NOT say "no isotope twin").
        return f"F{f}: 19F monoisotopic, fluorine count not isotope-confirmable"
    if f == 0 and hc < HC_FLOOR:     # genuine carbon-rich skeleton (F not displacing H)
        return f"H/C {hc:.2f} (carbon-rich)"
    if polarity == "+" and (br > 0 or cl > 0):
        return "halogen in neutral, +mode"
    return None


def scan(merged, *, polarity: str | None = None) -> list[dict]:
    """Flag Candidate-only neutrals that look implausible. Returns one dict per
    distinct neutral: {neutral_formula, reason, ion_score, tier}. A neutral that is
    Assigned in any ion channel is excluded (it is corroborated)."""
    if merged is None or "neutral_formula" not in getattr(merged, "columns", []):
        return []
    g = merged.dropna(subset=["neutral_formula"]).copy()
    if not len(g):
        return []
    g["neutral_formula"] = g["neutral_formula"].astype(str)
    has_tier = "tier" in g.columns
    out = []
    for f, sub in g.groupby("neutral_formula"):
        best = ("Assigned" if has_tier and (sub["tier"] == "Assigned").any()
                else "Candidate")
        reason = implausible(f, tier=best, polarity=polarity)
        if reason:
            sc = sub["ion_score"].max() if "ion_score" in sub.columns else None
            out.append({"neutral_formula": f, "reason": reason, "tier": best,
                        "ion_score": (float(sc) if sc is not None and sc == sc else None)})
    out.sort(key=lambda d: (d["reason"], d["neutral_formula"]))
    return out

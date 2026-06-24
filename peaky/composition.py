"""Composition accounting for the assignment report: signal-weighting and the
ammonium/amine degeneracy.

The composition page used to report distinct-neutral COUNTS by backbone only. For
a positive urea-CIMS batch that is misleading on two fronts, both addressed here:

  1. SIGNAL vs COUNT. A compound is counted once regardless of abundance, so a
     swarm of dim species can dominate the count while a few bright ones carry the
     chemistry. `signal_by_backbone` weights each distinct neutral by its summed
     M0 intensity, so the report can show e.g. "CHON is 47% by count but 6% by
     signal".

  2. THE AMMONIUM/AMINE SHADOW. `cleanup.prefer_amine_over_ammonium` re-reads an
     [M+NH4]+ adduct of a CHO neutral X as [M+H]+ of the amine X+NH3 (the two ions
     are mass- and isotope-identical). That conversion inflates the CHON count, and
     frequently the original CHO partner X is ALSO assigned (from its own [M+H]+),
     so the merged ledger carries BOTH X (CHO) and X+NH3 (CHON) as distinct
     neutrals for what may be a single molecule. `amine_shadow_stats` /
     `collapsed_composition` quantify and optionally collapse that degeneracy, so
     the page can report the count "two ways" (as-assigned vs ammonium-as-CHO).

All pure (formula arithmetic only); no I/O, no plotting. `neutral_signal` is a
{neutral_formula -> summed cps} map the caller builds from the per-file M0 rows.
"""
from __future__ import annotations

from collections import Counter

from . import chemistry as C


def backbone(formula: str) -> str:
    """CHOS if S present, else CHON if N present, else CHO. Si/F/Cl/Br are
    additions to the backbone, not a separate class (matches analyte_viz)."""
    c = C.parse_formula(str(formula))
    if c.get("S", 0):
        return "CHOS"
    if c.get("N", 0):
        return "CHON"
    return "CHO"


def minus_nh3(formula: str) -> str | None:
    """The CHO 'shadow' of an amine: X+NH3 -> X (remove one N and three H). The
    [M+H]+ of the returned neutral is the SAME ion as the [M+NH4]+ of `formula`'s
    de-aminated parent, i.e. the mass-degenerate partner. None if no NH3 to remove."""
    c = dict(C.parse_formula(str(formula)))
    if c.get("N", 0) < 1 or c.get("H", 0) < 3:
        return None
    c["N"] -= 1
    c["H"] -= 3
    if c["N"] == 0:
        del c["N"]
    return C.format_formula(c)


def _neutrals(merged) -> list[str]:
    return [str(f) for f in merged["neutral_formula"].dropna().unique() if str(f) != "nan"]


def count_by_backbone(merged) -> dict:
    """Distinct-neutral count per backbone class (the as-assigned composition)."""
    return dict(Counter(backbone(f) for f in _neutrals(merged)))


def signal_by_backbone(merged, neutral_signal: dict) -> tuple[dict, dict]:
    """Backbone composition weighted by summed M0 signal. Returns
    (fractions, absolute) where fractions sum to ~1 over classes with signal.
    A neutral missing from `neutral_signal` contributes 0 (it had no M0 height)."""
    absolute: dict = {}
    for f in _neutrals(merged):
        kl = backbone(f)
        absolute[kl] = absolute.get(kl, 0.0) + float(neutral_signal.get(f, 0.0) or 0.0)
    tot = sum(absolute.values()) or 1.0
    fractions = {k: v / tot for k, v in absolute.items()}
    return fractions, absolute


def amine_shadow_stats(merged) -> dict:
    """Quantify the ammonium/amine degeneracy in the distinct-neutral count.

    A 'shadowed' amine is an N-bearing neutral whose exact X-NH3 CHO twin is ALSO
    present in the ledger — the [M+H]+(amine) / [M+NH4]+(CHO) pair the re-read
    cannot distinguish, counted twice. Returns counts + a few examples.
    `collapsed_neutrals` is the distinct count after removing each shadowed amine
    (its CHO twin remains)."""
    neu = set(_neutrals(merged))
    amines = [f for f in neu if C.parse_formula(f).get("N", 0) >= 1]
    shadowed = []
    for f in amines:
        twin = minus_nh3(f)
        if twin is not None and twin in neu:
            shadowed.append((f, twin))
    return {
        "n_neutrals": len(neu),
        "n_amine": len(amines),
        "n_shadowed": len(shadowed),
        "collapsed_neutrals": len(neu) - len(shadowed),
        "examples": [f"{a}={t}+NH3" for a, t in sorted(shadowed)[:6]],
    }


def collapsed_composition(merged) -> tuple[dict, dict, int]:
    """Two-way backbone counts: (as_assigned, ammonium_as_cho, n_collapsed).

    `ammonium_as_cho` re-reads every shadowed amine (one with a present X-NH3 twin)
    back into its CHO twin's class — i.e. the composition if the parsimony NH4->amine
    re-read had NOT been applied to the cases where the bare CHO is independently
    seen. The twin already exists, so collapsing just removes the duplicate amine."""
    neu = set(_neutrals(merged))
    as_assigned = dict(Counter(backbone(f) for f in neu))
    dropped: Counter = Counter()
    for f in neu:
        if C.parse_formula(f).get("N", 0) >= 1:
            twin = minus_nh3(f)
            if twin is not None and twin in neu:
                dropped[backbone(f)] += 1
    collapsed = {k: as_assigned.get(k, 0) - dropped.get(k, 0) for k in as_assigned}
    return as_assigned, collapsed, int(sum(dropped.values()))


def top_species_by_signal(merged, neutral_signal: dict, *, n: int = 8) -> list[dict]:
    """Top-n distinct neutrals by summed M0 signal, with class + signal fraction.
    Useful for a findings page: the chemistry lives in a handful of bright peaks."""
    tot = sum(float(v or 0.0) for v in neutral_signal.values()) or 1.0
    rows = []
    seen = set()
    for f in _neutrals(merged):
        if f in seen:
            continue
        seen.add(f)
        rows.append({"neutral_formula": f, "signal": float(neutral_signal.get(f, 0.0) or 0.0),
                     "frac": float(neutral_signal.get(f, 0.0) or 0.0) / tot,
                     "klass": backbone(f)})
    rows.sort(key=lambda r: r["signal"], reverse=True)
    return rows[:n]


def oligomer_flag(merged, *, c_min: int = 18, c_max: int = 40, o_min: int = 7) -> list[str]:
    """Distinct neutrals that look like accretion / oligomer products (high carbon
    AND high oxygen) — the HOM dimers that are often the most event-specific signal.
    `c_max` excludes the absurdly large fits (C>40 in a monoterpene system is almost
    always a high-heteroatom mass coincidence, not a real oligomer). Returned sorted
    by carbon then oxygen; the caller may re-sort by signal."""
    out = []
    for f in _neutrals(merged):
        c = C.parse_formula(f)
        nc = c.get("C", 0)
        if c_min <= nc <= c_max and c.get("O", 0) >= o_min and not c.get("Si", 0):
            out.append((nc, c.get("O", 0), f))
    out.sort(key=lambda t: (-t[0], -t[1]))
    return [f for _, _, f in out]

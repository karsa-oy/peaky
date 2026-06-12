"""Automatic GKA series detection -- the machine version of 'rotating the plot'.

A homologous series is visible in a rotating GKA plot as a horizontal row at
some (base, X). Mathematically that is identical to peaks being linked by an
exact repeat-unit mass, which we can test directly: for each unit in a library,
count residual peaks whose partner at +1 unit also exists (within ppm).

Chance alignments are controlled with DECOY units: the same scan run at the
unit mass shifted by irrational offsets. A unit is significant when its link
count is well above the decoy mean AND absolutely large enough to matter.

The detector's output gates the pipeline: CF2 significant -> open the
'fluorinated' contaminant family; C2H6OSi -> 'siloxane'; HBr -> the cluster
ladder is already handled. The evidence table goes into the run manifest so
every data-driven decision is auditable.
"""
from __future__ import annotations

import bisect

import pandas as pd

from . import ledger as L

__version__ = "0.1.0"

# unit -> (exact mass, pipeline action when significant)
UNIT_LIBRARY: dict[str, tuple[float, str | None]] = {
    "CH2":     (14.015650, None),          # generic homology (pass-4 series)
    "O":       (15.994915, None),
    "H2O":     (18.010565, None),
    "CO2":     (43.989830, None),
    "C2H4O":   (44.026215, None),          # PEG / ethoxylate ladder
    "C3H6O":   (58.041865, None),          # PPG / propoxylate ladder
    "C2H4O2":  (60.021130, None),
    "CF2":     (49.996806, "fluorinated"),
    "C2F4":    (99.993612, "fluorinated"),
    "C2H6OSi": (74.018792, "siloxane"),
    "HBr":     (79.926160, None),          # cluster ladder handled separately
    "SO3":     (79.956815, "organosulfate"),
}

_DECOY_OFFSETS = (0.0317, -0.0473, 0.0689)   # Da; irrational-ish, off any unit


def _link_count(mz_sorted: list[float], unit: float, ppm: float) -> int:
    n = 0
    for m in mz_sorted:
        t = m + unit
        tol = t * ppm * 1e-6
        j = bisect.bisect_left(mz_sorted, t - tol)
        if j < len(mz_sorted) and mz_sorted[j] <= t + tol:
            n += 1
    return n


def _chain3_count(mz_sorted: list[float], unit: float, ppm: float) -> int:
    """Peaks that have partners at BOTH +1 and +2 units (chains of >=3)."""
    n = 0
    for m in mz_sorted:
        ok = True
        for k in (1, 2):
            t = m + k * unit
            tol = t * ppm * 1e-6
            j = bisect.bisect_left(mz_sorted, t - tol)
            if not (j < len(mz_sorted) and mz_sorted[j] <= t + tol):
                ok = False
                break
        if ok:
            n += 1
    return n


def detect_series(ledger: pd.DataFrame, *, ppm: float = 5.0,
                  min_height: float = 0.0,
                  include_low_confidence: bool = True,
                  min_links: int = 12, min_enrichment: float = 3.0,
                  units: dict | None = None) -> pd.DataFrame:
    """Scan the unexplained (+ optionally Low/Suspect) peaks for repeat-unit
    structure. Returns a DataFrame[unit, mass, n_links, n_chains3, decoy_mean,
    enrichment, significant, action].
    """
    units = units or UNIT_LIBRARY
    pop = ledger[ledger["role"] == L.ROLE_UNEXPLAINED]
    if include_low_confidence:
        low = ledger[(ledger["role"] == L.ROLE_M0)
                     & ledger["confidence"].astype(str).str.match(r"(Low|Suspect)")]
        pop = pd.concat([pop, low])
    pop = pop[pop["height"].fillna(0) >= min_height]
    mzs = sorted(pop["mz"].dropna().tolist())
    rows = []
    for name, (mass, action) in units.items():
        n = _link_count(mzs, mass, ppm)
        c3 = _chain3_count(mzs, mass, ppm)
        decoys = [_link_count(mzs, mass + d, ppm) for d in _DECOY_OFFSETS]
        dmean = sum(decoys) / len(decoys)
        enrich = n / dmean if dmean > 0 else (float(n) if n else 0.0)
        sig = bool(n >= min_links and enrich >= min_enrichment)
        rows.append({"unit": name, "mass": mass, "n_links": n, "n_chains3": c3,
                     "decoy_mean": round(dmean, 1), "enrichment": round(enrich, 2),
                     "significant": sig, "action": action})
    return pd.DataFrame(rows).sort_values("n_links", ascending=False).reset_index(drop=True)


def unit_members(ledger: pd.DataFrame, unit_mass: float, *, ppm: float = 5.0,
                 min_height: float = 0.0,
                 include_low_confidence: bool = True) -> set:
    """Peak_ids participating in at least one +/-1-unit link for this unit.
    Used to restrict an evidence-opened family's TARGETS to the chain members
    that justified opening it (prevents e.g. fluorine, a dense mass-filler,
    from claiming the whole residual)."""
    pop = ledger[ledger["role"] == L.ROLE_UNEXPLAINED]
    if include_low_confidence:
        low = ledger[(ledger["role"] == L.ROLE_M0)
                     & ledger["confidence"].astype(str).str.match(r"(Low|Suspect)")]
        pop = pd.concat([pop, low])
    pop = pop[pop["height"].fillna(0) >= min_height]
    pop = pop.sort_values("mz")
    mzs = pop["mz"].tolist()
    pids = pop["peak_id"].tolist()
    members: set = set()
    for i, m in enumerate(mzs):
        for sign in (+1, -1):
            t = m + sign * unit_mass
            tol = abs(t) * ppm * 1e-6
            j = bisect.bisect_left(mzs, t - tol)
            if j < len(mzs) and mzs[j] <= t + tol:
                members.add(pids[i])
                break
    return members


def unit_chains(ledger: pd.DataFrame, unit_mass: float, *, ppm: float = 5.0,
                min_height: float = 0.0, min_len: int = 2,
                include_low_confidence: bool = True) -> list[list[tuple]]:
    """Maximal chains of peaks spaced by exactly `unit_mass`. Each chain is a
    list of (peak_id, mz) sorted ascending; only chains with >= min_len members
    are returned. Chain heads have no -1-unit partner."""
    pop = ledger[ledger["role"] == L.ROLE_UNEXPLAINED]
    if include_low_confidence:
        low = ledger[(ledger["role"] == L.ROLE_M0)
                     & ledger["confidence"].astype(str).str.match(r"(Low|Suspect)")]
        pop = pd.concat([pop, low])
    pop = pop[pop["height"].fillna(0) >= min_height].sort_values("mz")
    mzs = pop["mz"].tolist()
    pids = pop["peak_id"].tolist()

    def _partner(m, direction):
        t = m + direction * unit_mass
        tol = abs(t) * ppm * 1e-6
        j = bisect.bisect_left(mzs, t - tol)
        if j < len(mzs) and mzs[j] <= t + tol:
            return j
        return None

    chains: list[list[tuple]] = []
    used: set = set()
    for i, m in enumerate(mzs):
        if pids[i] in used:
            continue
        if _partner(m, -1) is not None:
            continue   # not a head
        chain = [(pids[i], m)]
        cur = m
        while True:
            j = _partner(cur, +1)
            if j is None:
                break
            chain.append((pids[j], mzs[j]))
            cur = mzs[j]
        if len(chain) >= min_len:
            for pid, _ in chain:
                used.add(pid)
            chains.append(chain)
    return chains


def families_from_evidence(evidence: pd.DataFrame) -> list[str]:
    """Contaminant families to open based on significant units."""
    out: list[str] = []
    for _, r in evidence.iterrows():
        if r["significant"] and r["action"] and r["action"] not in out:
            out.append(r["action"])
    return out

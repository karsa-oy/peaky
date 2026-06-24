"""Reference peaklists — a curated catalog of known-molecule formula lists from
published chemical systems (oxidation HOM, contaminant families, ...), used as a
SOFT, context-gated prior.

Two uses, both decoupled from the expensive server scoring (one offline pass):
  * corroborate  -- a Candidate-tier neutral whose formula is on a list for the
                    sample's chemistry gets a literature-corroboration tag.
  * rescue/annotate -- an UNEXPLAINED peak (no formula) is matched BY MASS, under
                    the run's reagent adducts, against the list formulas, turning
                    the remainder into leads the reader can chase.

This is NOT the Pass-0 known-species lock (`passes._known_species`, for
contaminants the grid cannot reach). It is a post-hoc PRIOR: it never overrides an
isotope-scored Identified, never fabricates a peak, and every match carries the
source list id for provenance. Masses are recomputed for whatever reagent adduct
the run uses (`chemistry.ion_mz`), so one list serves Br-/NO3-/I-/urea+ runs alike.

A list is one self-describing JSON in data/peaklists/ (see README.md). Add a
credible source = drop in one JSON; nothing else couples.
"""
from __future__ import annotations

import bisect
import glob
import json
import os
from dataclasses import dataclass

from . import chemistry as C

__version__ = "0.1.0"

_DIR = os.path.join(os.path.dirname(__file__), "data", "peaklists")

# batch-name / label keywords -> experimental-context tag (the metadata "unlock").
# Conservative: only fire on clearly source-diagnostic words.
CONTEXT_KEYWORDS = {
    "monoterpene_ox": ("monoterpene", "pinene", "limonene", "terpene", "orange", "carene"),
    "limonene_ox": ("limonene", "orange", "d-limonene"),
    "ap_ox": ("pinene", "a-pinene", "apinene"),
}


@dataclass(frozen=True)
class ReferenceList:
    id: str
    system: str
    label: str
    data_version: str
    polarity: str
    native_detection: str
    applies_to_contexts: tuple
    references: tuple
    formulas: frozenset            # closed-shell neutral formulas (matchable)
    radicals: frozenset            # odd-H radical formulas (excluded by default)
    conditions_of: dict            # formula -> tuple(conditions)
    source_file: str
    always_active: bool = False    # universal lists (e.g. contaminants) ignore context gating
    meta_of: dict = None           # formula -> {name, origin, ...} (display extras)

    def pool(self, include_radicals: bool = False) -> frozenset:
        return self.formulas | self.radicals if include_radicals else self.formulas

    def cite(self) -> str:
        if not self.references:
            return self.label
        r = self.references[0]
        bits = [r.get("authors"), r.get("title"), r.get("publisher"),
                str(r.get("year", "")), r.get("section")]
        return ", ".join(b for b in bits if b)


# ---------------------------------------------------------------------------
def load_catalog(directory: str | None = None) -> dict:
    """Load every *.json reference list under `directory` (default: packaged
    data/peaklists). Returns {id: ReferenceList}."""
    directory = directory or _DIR
    out: dict = {}
    for p in sorted(glob.glob(os.path.join(directory, "*.json"))):
        with open(p, encoding="utf-8") as fh:
            d = json.load(fh)
        sp = d.get("species", [])
        closed, rad, cond, meta = set(), set(), {}, {}
        for s in sp:
            f = s.get("formula")
            if not f:
                continue
            (rad if s.get("radical") else closed).add(f)
            cond[f] = tuple(s.get("conditions", ()))
            extra = {k: s[k] for k in ("name", "origin", "modes") if k in s}
            if extra:
                meta[f] = extra
        out[d["id"]] = ReferenceList(
            id=d["id"], system=d.get("system", ""), label=d.get("label", d["id"]),
            data_version=str(d.get("data_version", "")),
            polarity=d.get("polarity", ""), native_detection=d.get("native_detection", ""),
            applies_to_contexts=tuple(d.get("applies_to_contexts", ())),
            references=tuple(d.get("references", ())),
            formulas=frozenset(closed), radicals=frozenset(rad), conditions_of=cond,
            source_file=os.path.basename(p),
            always_active=bool(d.get("always_active", False)), meta_of=meta)
    return out


def resolve_context_tags(*texts: str) -> set:
    """Infer experimental-context tags from run metadata (batch name / label).
    This is the 'unlock with metadata' step — only the matched contexts' lists
    become active, keeping precision high on unrelated samples."""
    blob = " ".join(t for t in texts if t).lower()
    return {tag for tag, kws in CONTEXT_KEYWORDS.items()
            if any(k.lower() in blob for k in kws)}


def active_lists(catalog: dict, *, context_tags=()) -> list:
    """Select active lists: every `always_active` list (ubiquitous, e.g. lab
    contaminants) PLUS any whose `applies_to_contexts` intersect the run's context
    tags (the metadata unlock). Polarity is NOT a filter: a neutral-formula library
    transfers across reagent ion forms."""
    tags = set(context_tags)
    return [L for L in catalog.values()
            if L.always_active or (tags and set(L.applies_to_contexts) & tags)]


# ---------------------------------------------------------------------------
def match_assigned(neutral_formulas, lists, *, include_radicals: bool = False) -> dict:
    """Formula-membership corroboration: which assigned neutrals appear on a list.
    Returns {formula: [{list, conditions}]}."""
    res: dict = {}
    for f in {str(x) for x in neutral_formulas if x and str(x) != "nan"}:
        hits = [{"list": L.id, "conditions": list(L.conditions_of.get(f, ()))}
                for L in lists if f in L.pool(include_radicals)]
        if hits:
            res[f] = hits
    return res


def _target_table(lists, adducts, include_radicals: bool):
    """Sorted [(mz, formula, adduct, list_id)] of every list formula under every
    reagent adduct — built once per match call."""
    rows = []
    for L in lists:
        for f in L.pool(include_radicals):
            for a in adducts:
                try:
                    rows.append((C.ion_mz(f, a), f, a, L.id))
                except Exception:
                    continue
    rows.sort(key=lambda r: r[0])
    return rows


def match_by_mass(mz_values, lists, adducts, *, tol_ppm: float = 5.0,
                  include_radicals: bool = False) -> list:
    """Rescue/annotate UNEXPLAINED peaks (which have no formula) BY MASS: for each
    observed m/z, find the closest list-formula ion (any reagent adduct) within
    `tol_ppm`. Returns one dict per matched observed peak (best match only):
    {obs_mz, formula, adduct, list, ppm, target_mz}."""
    targets = _target_table(lists, adducts, include_radicals)
    if not targets:
        return []
    tmz = [t[0] for t in targets]
    out = []
    for obs in mz_values:
        try:
            obs = float(obs)
        except (TypeError, ValueError):
            continue
        if obs <= 0:
            continue
        w = obs * tol_ppm * 1e-6
        lo = bisect.bisect_left(tmz, obs - w)
        hi = bisect.bisect_right(tmz, obs + w)
        best = None
        for j in range(lo, hi):
            m, f, a, lid = targets[j]
            ppm = (obs - m) / obs * 1e6
            if abs(ppm) <= tol_ppm and (best is None or abs(ppm) < abs(best["ppm"])):
                best = {"obs_mz": round(obs, 5), "formula": f, "adduct": a,
                        "list": lid, "ppm": round(ppm, 2), "target_mz": round(m, 5)}
        if best:
            out.append(best)
    return out

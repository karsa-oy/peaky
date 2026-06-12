"""The ONLY module that talks to Mascope.

It wraps the mascope-sdk MascopeClient and exposes exactly the operations the
pipeline needs:

  * connect()                  -- build a client from ~/mascope-mcp/.env
  * fetch_peaks()              -- pull + cache the raw peak table
  * resolve_mechanism_ids()    -- ionization name -> id
  * query_candidates()         -- cheminfo formula enumerator for one m/z
  * score_candidates()         -- match_compounds -> flat per-isotopologue table

The scoring oracle is Mascope: match_compounds returns a compound -> ion ->
isotopologue tree, every node carrying its own match_score and (for isotopes)
the attributed sample_peak_id. flatten_match_tree() turns that tree into a flat
table and is a PURE function, unit-tested offline against a captured fixture.
"""
from __future__ import annotations

import os
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

__version__ = "0.1.0"

DEFAULT_ENV = os.path.expanduser("~/mascope-mcp/.env")
CACHE_ROOT = Path(os.path.expanduser("~/.mascope-assign-cache"))
MATCH_BATCH = 200   # match_compounds times out above ~500

# Default server-side match parameters. mz_tolerance is INTEGER ppm.
DEFAULT_MATCH_PARAMS = {
    "mz_tolerance": 5,
    "isotope_ratio_tolerance": 0.2,
    "peak_min_intensity": 0.0,
    "min_isotope_abundance": 0.15,
    "min_isotope_correlation": 0.7,
    "probable_match_threshold": 0.8,
    "possible_match_threshold": 0.4,
}

# bracketed heavy-isotope tokens, e.g. [13C], [13C]2, [81Br], [18O], [37Cl]
_ISO_TOKEN = re.compile(r"\[(\d+[A-Z][a-z]?)\](\d*)")


# ---------------------------------------------------------------------------
# Connection
# ---------------------------------------------------------------------------
def connect(env_path: str = DEFAULT_ENV):
    """Build a MascopeClient from the .env (MASCOPE_URL, MASCOPE_ACCESS_TOKEN)."""
    from dotenv import load_dotenv
    load_dotenv(env_path)
    url = os.environ.get("MASCOPE_URL")
    tok = os.environ.get("MASCOPE_ACCESS_TOKEN")
    if not url or not tok:
        raise RuntimeError(f"MASCOPE_URL / MASCOPE_ACCESS_TOKEN not found in {env_path}")
    from mascope_sdk import MascopeClient
    return MascopeClient(url=url, access_token=tok)


# ---------------------------------------------------------------------------
# Peaks
# ---------------------------------------------------------------------------
def fetch_peaks(client, sample_id: str, *, use_cache: bool = True,
                cache_root: Path = CACHE_ROOT) -> pd.DataFrame:
    """Pull the raw peak table (with Mascope's own matches flattened in) and
    cache it. Returns the full multi-row-per-peak frame; dedup is the ledger's
    job."""
    cdir = Path(cache_root) / sample_id
    cfile = cdir / "peaks.parquet"
    if use_cache and cfile.exists():
        return pd.read_parquet(cfile)
    peaks = client.samples.get_peaks(sample_id=sample_id, matches=True)
    if peaks is None or len(peaks) == 0:
        raise RuntimeError(f"no peaks returned for sample {sample_id!r}")
    cdir.mkdir(parents=True, exist_ok=True)
    try:
        peaks.to_parquet(cfile)
    except Exception:
        peaks.to_csv(cdir / "peaks.csv", index=False)
    return peaks


def resolve_mechanism_ids(client, names: list[str]) -> dict[str, str]:
    """Map ionization-mechanism names (e.g. '-H+', '+Br-') to their ids."""
    table = client.ionization.list()
    by_name = {r.ionization_mechanism: r.ionization_mechanism_id
               for r in table.itertuples()}
    out: dict[str, str] = {}
    for n in names:
        if n in by_name:
            out[n] = by_name[n]
    return out


# Map our adduct labels to the server's ionization-mechanism names.
ADDUCT_TO_MECH = {
    "[M-H]-": "-H+",
    "[M+Br]-": "+Br-",
    "[M+Cl]-": "+Cl-",
    "[M+I]-": "+I-",
    "[M+NO3]-": "+NO3-",
    "[M+HSO4]-": "+HSO4-",
    "[M+H]+": "+H+",
    "[M+Na]+": "+Na+",
    "[M+NH4]+": "+NH4+",
    "[M+CO3]-": "+CO3-",
}
MECH_TO_ADDUCT = {v: k for k, v in ADDUCT_TO_MECH.items()}


def detect_adducts(peaks: pd.DataFrame) -> list[str]:
    """Infer the reagent/adduct system from the sample's own peak matches
    (the `ionization_mechanism` column). This is what makes a Br-CIMS sample
    get [M+Br]- offered as an interpretation instead of forcing Br into the
    neutral. Falls back to [M-H]- if nothing is recognised."""
    if peaks is None or "ionization_mechanism" not in peaks.columns:
        return ["[M-H]-"]
    out: list[str] = []
    for name in peaks["ionization_mechanism"].dropna().unique():
        a = MECH_TO_ADDUCT.get(str(name))
        if a and a not in out:
            out.append(a)
    return out or ["[M-H]-"]


# ---------------------------------------------------------------------------
# Candidate enumeration (cheminfo)
# ---------------------------------------------------------------------------
def query_candidates(client, mz: float, mechanism_ids: list[str], *,
                     formula_ranges: str, ppm: float = 5.0,
                     limit: int = 25) -> list[str]:
    """Return candidate NEUTRAL formulas for one m/z (deduped).

    cheminfo is a flaky endpoint (timeouts / 500s). A failure here must NOT kill
    the run: the local grid covers the same CHO/CHON formula space, so degrade
    to [] on any error and let the grid carry that m/z."""
    try:
        res = client.cheminfo.query_by_mz(
            mz=mz, ionization_mechanism_ids=mechanism_ids,
            formula_ranges=formula_ranges, mz_tolerance=float(ppm), limit=limit) or []
    except Exception:
        return []
    out = []
    seen = set()
    for r in res:
        f = r.get("target_compound_formula")
        if f and f not in seen:
            seen.add(f)
            out.append(f)
    return out


def query_candidates_bulk(client, mzs: list[float], mechanism_ids: list[str], *,
                          formula_ranges: str, ppm: float = 5.0, limit: int = 25,
                          workers: int = 12) -> dict[float, list[str]]:
    """Parallel cheminfo enumeration over many m/z."""
    def _one(mz):
        return mz, query_candidates(client, mz, mechanism_ids,
                                    formula_ranges=formula_ranges, ppm=ppm, limit=limit)
    out: dict[float, list[str]] = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for mz, cands in ex.map(_one, mzs):
            out[mz] = cands
    return out


# ---------------------------------------------------------------------------
# Scoring oracle (match_compounds) + PURE parser
# ---------------------------------------------------------------------------
def parse_isotope_label(isotope_formula: str) -> tuple[str, bool]:
    """('[13C]C3H5O2-') -> ('13C', False). Base (no heavy isotope) -> ('M0', True)."""
    toks = _ISO_TOKEN.findall(isotope_formula or "")
    if not toks:
        return "M0", True
    parts = []
    for sym, n in toks:
        parts.append(sym + (n if n and int(n) > 1 else ""))
    return "+".join(parts), False


def flatten_match_tree(tree: list[dict]) -> pd.DataFrame:
    """PURE: flatten match_compounds output into one row per
    (compound, ion, isotopologue). No network. Columns:

      compound_formula, compound_score, compound_category,
      ion_formula, ion_score, ion_category, mechanism_id,
      isotope_formula, iso_label, is_base, theo_mz, rel_abundance,
      iso_score, iso_category, sample_peak_id, sample_peak_mz,
      sample_peak_intensity, ppm_error, abundance_error
    """
    rows: list[dict] = []
    for comp in tree or []:
        cf = comp.get("target_compound_formula")
        cs = comp.get("match_score")
        cc = comp.get("match_category")
        for ion in comp.get("children", []) or []:
            ifl = ion.get("target_ion_formula")
            isc = ion.get("match_score")
            icat = ion.get("match_category")
            mech = ion.get("ionization_mechanism_id")
            for iso in ion.get("children", []) or []:
                iso_f = iso.get("target_isotope_formula")
                label, is_base = parse_isotope_label(iso_f)
                theo = iso.get("mz")
                spk_mz = iso.get("sample_peak_mz")
                spk_int = iso.get("sample_peak_intensity") or 0.0
                spid = iso.get("sample_peak_id") or None
                # ppm is only meaningful for a GENUINELY matched isotope (a real
                # attributed peak). Unmatched/forced nodes carry sample_peak_mz ==
                # theoretical mz and zero intensity -> leave ppm undefined.
                matched = bool(spid) and float(spk_int) > 0 and spk_mz and float(spk_mz) > 0
                if matched and theo:
                    ppm = (float(spk_mz) - float(theo)) / float(theo) * 1e6
                else:
                    ppm = None
                rows.append({
                    "compound_formula": cf, "compound_score": cs, "compound_category": cc,
                    "ion_formula": ifl, "ion_score": isc, "ion_category": icat,
                    "mechanism_id": mech,
                    "isotope_formula": iso_f, "iso_label": label, "is_base": is_base,
                    "theo_mz": theo, "rel_abundance": iso.get("relative_abundance"),
                    "iso_score": iso.get("match_score"), "iso_category": iso.get("match_category"),
                    "sample_peak_id": spid, "sample_peak_mz": spk_mz,
                    "sample_peak_intensity": iso.get("sample_peak_intensity"),
                    "ppm_error": ppm, "abundance_error": iso.get("match_abundance_error"),
                })
    return pd.DataFrame(rows)


MATCH_WORKERS = 5   # concurrent match_compounds batches (I/O-bound; server-safe)


def score_candidates(client, sample_id: str, formulas: list[str], *,
                     match_params: dict | None = None,
                     mechanism_ids: list[str] | None = None,
                     batch: int = MATCH_BATCH,
                     workers: int = MATCH_WORKERS) -> pd.DataFrame:
    """Score candidate NEUTRAL formulas against the sample. Returns the flat
    per-isotopologue table (see flatten_match_tree). Batches are scored
    CONCURRENTLY -- match_compounds is network-bound, so the wall-clock for a
    many-batch pass (e.g. a wide heteroatom family) scales with the worker
    count, not the batch count. Per-batch failures degrade to empty, never
    losing the other batches."""
    formulas = sorted({f for f in formulas if f})
    if not formulas:
        return pd.DataFrame()
    mp = dict(DEFAULT_MATCH_PARAMS)
    if match_params:
        mp.update(match_params)
    mp["mz_tolerance"] = int(round(mp.get("mz_tolerance", 5)))  # server needs int ppm
    chunks = [formulas[i:i + batch] for i in range(0, len(formulas), batch)]

    def _score(chunk):
        try:
            tree = client.matching.match_compounds(
                sample_id=sample_id, formulas=chunk, match_params=mp,
                ionization_mechanism_ids=mechanism_ids)
        except Exception:
            tree = None
        return flatten_match_tree(tree or [])

    if len(chunks) == 1:
        frames = [_score(chunks[0])]
    else:
        with ThreadPoolExecutor(max_workers=min(workers, len(chunks))) as ex:
            frames = list(ex.map(_score, chunks))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

"""Certified-neutral discovery: multi-channel / cluster-ladder mass certification.

THE MECHANISM. A single peak's neutral interpretation is mass-degenerate: at
m/z 300+ with a +-3 ppm window, dozens of formulas fit, and monoisotopic
elements (P) make it worse -- which is why the per-peak grid keeps P closed.
But when the SAME neutral core mass is reachable from >=2 DISTINCT ion channels
in one spectrum -- a different adduct ([M+H]+ vs [M+NH4]+) or a different
cluster order on the reagent ladder ([M+H]+, [M+urea+H]+, [M+2urea+H]+) -- the
channels impose N independent mass constraints on ONE unknown. Their
convergence *certifies* the neutral, and only for a certified mass do we open
the expanded element box (P, higher S, Cl) that the per-peak grid forbids.

This INVERTS pass-5 (completion): pass-5 walks from KNOWN neutrals to their
missing cross-channel partners; this module walks from UNKNOWN peak groups to
a licensed neutral mass. Validated by hand on real data before being built:

  * NBBS ladder -- 214.0896 [M+H]+, 274.1220 [M+urea+H]+, 334.1544
    [M+2urea+H]+ all back-calculate to core 213.0823 within 0.1 mDa
    (C10H15NO2S, N-butylbenzenesulfonamide instrument background).
  * malathion -- 331.0433 [M+H]+ / 391.0756 [M+urea+H]+ / 348.0699 [M+NH4]+
    -> core 330.0360 (C10H19O6PS2), invisible to the P-free grid.

Everything here is PURE: deterministic functions of their inputs, no client,
no ledger mutation, no I/O. The pipeline entry (run_pass_certified) and the
optional time-series corroboration layer live elsewhere.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from ..chem import chemistry as C
from ..chem import reagents as RG

__all__ = [
    "ChannelHit", "Certificate",
    "channel_offsets", "find_certificates", "enumerate_certified",
    "ts_covariation", "corroboration_count",
]

# A certificate's channels must converge this tightly in the NEUTRAL-MASS
# domain. The NBBS ground truth converges within 0.1 mDa; 3 mDa accommodates
# real calibration scatter while keeping random co-convergence rare (two
# unrelated peaks landing the same back-calc mass within 3 mDa under a
# specific offset pair is ~1e-2 per pair even in a dense spectrum; demanding
# DISTINCT channel types cuts it further).
CORE_TOL_MDA = 3.0
# Offsets closer than this are the same physical channel (the urea alias:
# ["[M+H]+" + 1 urea] == "[M+(CH4N2O)H]+" exactly) -> de-duplicated.
ALIAS_TOL = 1e-4
MIN_CORE_MASS = 30.0


@dataclass(frozen=True)
class ChannelHit:
    """One (peak, channel) interpretation contributing to a certified core."""
    peak_id: str
    mz: float
    adduct: str          # peaky adduct label (a chemistry.ADDUCT_SHIFTS key)
    cluster_order: int   # 0 = bare adduct; n>=1 = [M + n*reagent + adduct]
    back_calc_mass: float
    height: float = 0.0


@dataclass(frozen=True)
class Certificate:
    """A neutral core mass certified by >=2 distinct converging channels."""
    core_mass: float               # intensity-weighted mean of back_calc_mass
    hits: tuple[ChannelHit, ...]   # sorted by m/z
    n_channels: int                # distinct (adduct, cluster_order) pairs
    spread_mda: float              # (max-min back_calc_mass) * 1000

    @property
    def total_order(self) -> int:
        """Total cluster units assumed across the member interpretations --
        the parsimony key: a ladder of k peaks certifies SEVERAL cores shifted
        by whole repeat units (each peak re-read with one more/fewer urea);
        the reading that assumes the FEWEST attached reagent molecules is the
        physical one (its lowest rung is the bare adduct)."""
        return sum(h.cluster_order for h in self.hits)


def _cluster_unit(reagent: str | None) -> float | None:
    """Neutral repeat-unit mass for reagents that form molecular cluster
    ladders ([M + n*R + adduct]). Only molecular reagents qualify (urea);
    halide reagents' extra-X channels are already distinct ADDUCT_SHIFTS keys
    ([M+Br2]-, [M+HBr+Br]- ...) -- building a synthetic ladder for them would
    double-count those channels."""
    if reagent and reagent in RG._POSITIVE_REAGENTS:
        return C.neutral_mass(RG._POSITIVE_REAGENTS[reagent])
    return None


def channel_offsets(
    adducts: list[str],
    reagent: str | None = None,
    *,
    max_cluster_order: int = 2,
) -> list[tuple[str, int, float]]:
    """Enumerate (adduct, cluster_order, total_offset) triples observable in a
    single spectrum for this profile: every adduct at order 0, plus -- for a
    molecular reagent (urea) -- the same adducts shifted by n repeat units.

    Aliased offsets (identical within ALIAS_TOL, e.g. [M+H]+ at order 1 vs the
    registered urea adduct at order 0) are de-duplicated keeping the REGISTERED
    adduct form (lowest cluster_order), so one physical ion can never supply
    two "distinct" channels to a certificate.
    """
    unit = _cluster_unit(reagent)
    triples: list[tuple[str, int, float]] = []
    for ad in adducts:
        shift = C.ADDUCT_SHIFTS.get(ad)
        if shift is None:
            continue
        orders = range(max_cluster_order + 1) if unit else (0,)
        for n in orders:
            triples.append((ad, n, shift + (n * unit if unit else 0.0)))
    # de-dup aliases: sort so the lowest cluster_order wins its offset bin
    triples.sort(key=lambda t: (round(t[2] / ALIAS_TOL), t[1], t[0]))
    out: list[tuple[str, int, float]] = []
    for ad, n, off in triples:
        if out and abs(off - out[-1][2]) < ALIAS_TOL:
            continue
        out.append((ad, n, off))
    return out


def find_certificates(
    peaks: pd.DataFrame,
    offsets: list[tuple[str, int, float]],
    *,
    tol_mda: float = CORE_TOL_MDA,
    min_channels: int = 2,
    mz_col: str = "mz",
    id_col: str = "peak_id",
    intensity_col: str = "height",
) -> list[Certificate]:
    """Group peaks whose back-calculated neutral masses converge.

    Every (peak, offset) pair yields a candidate core mass = mz - offset
    (dropped below MIN_CORE_MASS). Candidates are single-linkage clustered
    within tol_mda; a cluster qualifies as a Certificate when it holds
    >= min_channels DISTINCT (adduct, cluster_order) channels from DISTINCT
    peaks. Deterministic: stable sorts throughout.
    """
    hits: list[ChannelHit] = []
    heights = peaks[intensity_col] if intensity_col in peaks.columns else None
    for row_i, (pid, mz) in enumerate(zip(peaks[id_col], peaks[mz_col])):
        h = float(heights.iloc[row_i]) if heights is not None else 0.0
        for ad, n, off in offsets:
            core = float(mz) - off
            if core < MIN_CORE_MASS:
                continue
            hits.append(ChannelHit(str(pid), float(mz), ad, n, core, h))
    hits.sort(key=lambda x: (x.back_calc_mass, x.peak_id, x.adduct, x.cluster_order))

    certs: list[Certificate] = []
    tol = tol_mda / 1000.0
    i = 0
    while i < len(hits):
        # greedy single-linkage: extend while consecutive gaps stay inside tol
        j = i + 1
        while j < len(hits) and hits[j].back_calc_mass - hits[j - 1].back_calc_mass <= tol:
            j += 1
        group = hits[i:j]
        i = j
        # one peak may appear once per cluster; one channel type once per cert
        by_channel: dict[tuple[str, int], ChannelHit] = {}
        seen_peaks: set[str] = set()
        for hcand in sorted(group, key=lambda x: -x.height):
            key = (hcand.adduct, hcand.cluster_order)
            if key in by_channel or hcand.peak_id in seen_peaks:
                continue
            by_channel[key] = hcand
            seen_peaks.add(hcand.peak_id)
        if len(by_channel) < min_channels:
            continue
        members = tuple(sorted(by_channel.values(), key=lambda x: x.mz))
        spread = (max(m.back_calc_mass for m in members)
                  - min(m.back_calc_mass for m in members)) * 1000.0
        if spread > tol_mda:
            continue
        wsum = sum(m.height for m in members)
        core = (sum(m.back_calc_mass * m.height for m in members) / wsum
                if wsum > 0 else
                sum(m.back_calc_mass for m in members) / len(members))
        certs.append(Certificate(core, members, len(members), spread))
    return _select_certificates(certs)


def _select_certificates(certs: list[Certificate]) -> list[Certificate]:
    """Resolve the ladder-shift ambiguity. A k-rung cluster ladder converges at
    SEVERAL cores offset by whole repeat units (the same peaks re-read with
    +-n reagent molecules each). Greedily keep, per overlapping peak set, the
    reading with the MOST channels and then the FEWEST assumed cluster units
    (total_order) -- for the NBBS ladder that picks core 213.0823 (lowest rung
    = bare [M+H]+) over the shifted 153.05/273.11 readings. Deterministic."""
    ranked = sorted(certs, key=lambda c: (-c.n_channels, c.total_order,
                                          -c.core_mass, c.hits[0].peak_id))
    taken: set[str] = set()
    out: list[Certificate] = []
    for c in ranked:
        pids = {h.peak_id for h in c.hits}
        if pids & taken:
            continue
        taken |= pids
        out.append(c)
    out.sort(key=lambda c: c.core_mass)
    return out


# elements a certificate is allowed to open beyond the per-peak grid, with the
# hard caps used when the caller FORCES them open (the whole point: the
# per-peak grid keeps P shut because a lone mass can't earn it; a
# multi-channel certificate can).
_CERT_EXTRA = {"P": (0, 2), "S": (0, 3), "Cl": (0, 2)}


def enumerate_certified(
    core_mass: float,
    profile,
    *,
    tol_mda: float = CORE_TOL_MDA,
    extra_elements: dict[str, tuple[int, int]] | None = None,
    force: bool = True,
) -> list[str]:
    """Enumerate neutral formulas within +-tol_mda of a CERTIFIED core mass
    under the expanded element box.

    Base CHON box from the context profile (grid_c_max/grid_o_max/max_N); the
    extra elements (default P/S/Cl per _CERT_EXTRA) are opened at their
    certificate caps when force=True (bypassing the profile's max_<el>, which
    exists to protect the PER-PEAK grid), else clamped to the profile caps.
    """
    extra = dict(_CERT_EXTRA if extra_elements is None else extra_elements)
    ranges: dict[str, tuple[int, int]] = {
        "C": (1, int(getattr(profile, "grid_c_max", 40))),
        "H": (0, 2 * int(getattr(profile, "grid_c_max", 40)) + 2 + 8),
        "N": (0, min(6, int(getattr(profile, "max_N", 6)))),
        "O": (0, int(getattr(profile, "grid_o_max", 30))),
    }
    for el, (lo, hi) in extra.items():
        cap = int(getattr(profile, f"max_{el}", hi))
        ranges[el] = (lo, hi if force else min(hi, cap))
    tol = tol_mda / 1000.0
    grid = C.enumerate_grid(ranges, core_mass - tol, core_mass + tol)
    return sorted({f for m, f in grid if abs(m - core_mass) <= tol})


def ts_covariation(
    ts_peaks: pd.DataFrame,
    mzs: list[float],
    *,
    tol: float = 0.006,
    mz_col: str = "mz",
    time_col: str = "datetime_utc",
    intensity_col: str = "height",
    min_points: int = 8,
) -> float | None:
    """OPTIONAL time-series corroboration: the minimum pairwise log-intensity
    correlation among a certificate's member channels across a batch TS.

    Channels carrying the SAME neutral must co-vary; a strongly negative or
    near-zero minimum says the pairing is coincidental. Returns None when the
    TS cannot support the test (missing channels / too few shared points), so
    the caller treats it as no-evidence rather than refutation. Pure.
    """
    import numpy as np

    traces = []
    for t in mzs:
        s = ts_peaks[(ts_peaks[mz_col] - t).abs() < tol]
        if not len(s):
            return None
        tr = s.groupby(time_col)[intensity_col].max()
        if len(tr) < min_points:
            return None
        traces.append(np.log10(tr.clip(lower=1)))
    r_min = 1.0
    for i in range(len(traces)):
        for j in range(i + 1, len(traces)):
            joined = pd.concat([traces[i], traces[j]], axis=1, join="inner").dropna()
            if len(joined) < min_points:
                return None
            r = float(joined.iloc[:, 0].corr(joined.iloc[:, 1]))
            r_min = min(r_min, r)
    return r_min


def corroboration_count(
    n_channels: int,
    *,
    iso_confirmed: bool = False,
    ts_covaries: bool = False,
    order_lift_consistent: bool = False,
) -> int:
    """Count INDEPENDENT corroborations backing a certified commit.

    Sources: extra ion channels beyond the first (n_channels - 1), a confirmed
    diagnostic heavy-isotope envelope (34S/37Cl/81Br -- never 13C, which every
    C-bearing formula has), time-series co-variation of the member channels,
    and a cluster-order zeroing-lift gradient consistent with the adduct
    stoichiometry (reagent batches only). The pipeline maps counts to tiers:
    >=2 with isotope evidence -> Assigned-grade, else Candidate.
    """
    n = max(0, n_channels - 1)
    n += 1 if iso_confirmed else 0
    n += 1 if ts_covaries else 0
    n += 1 if order_lift_consistent else 0
    return n

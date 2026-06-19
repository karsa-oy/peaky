"""Pipeline spine — `run(batch, reagent)` instead of a copy-pasted script per
dataset. Resolves a reagent profile, loads the batch (SDK or a cached parquet),
and dispatches the requested stages, all parameterised by the profile.

This is the thin orchestrator the rest of agent-peaky hangs off. Stage functions
live in their own modules and are called with profile params — nothing reagent-
or batch-specific is hardcoded here.

Stages: 'matrix' (TS intensity matrix) is wired; 'assign' / 'cluster' / 'validate'
are folded in as their scratch logic is consolidated into the package.
"""
from __future__ import annotations

import os
import pandas as pd

from . import io_mascope as IO
from . import profiles as P
from . import sampling as SS
from . import timeseries as TS

__version__ = "0.2.0"  # + representative-sample selection (5 time-grid + max-TIC)

STAGES = ("matrix", "assign", "cluster", "validate")


def load(*, batch: str | None = None, dataset: str | None = None,
         peaks: "str | pd.DataFrame | None" = None, save_path: str | None = None
         ) -> pd.DataFrame:
    """Get the batch peak time-series — from a parquet/DataFrame if given (offline,
    cached), else fetched from Mascope via the SDK."""
    if peaks is not None:
        return pd.read_parquet(os.path.expanduser(peaks)) if isinstance(peaks, str) else peaks
    if not (batch and dataset):
        raise ValueError("need peaks=, or both batch= and dataset=")
    return IO.fetch_batch_peaks(IO.connect(), dataset, batch, save_path=save_path)


def run(*, batch: str | None = None, dataset: str | None = None,
        peaks: "str | pd.DataFrame | None" = None, reagent: str = "auto",
        stages: tuple = ("matrix",), out_dir: str | None = None,
        n_time: int = SS.N_TIME, include_max_tic: bool = True) -> dict:
    """Run the pipeline on one batch.

    Returns a dict with at least {profile, peaks, n_samples, assign_samples}
    plus per-stage outputs. Pass reagent by name ('Br'/'Ur') or 'auto' to detect
    from the data.

    THE RULE (always computed, regardless of stages): `assign_samples` is the
    representative subset to assign + merge — `n_time` samples evenly spaced in
    TIME plus the max-TIC sample (see sampling.py). `assign_sample_ids` is the
    bare id list. Assignment runs on these, not on a single averaged file, so the
    merged peak list covers analytes that only appear at part of the run.
    """
    pk = load(batch=batch, dataset=dataset, peaks=peaks)
    prof = P.resolve(reagent, pk)
    n_samples = pk["sample_item_id"].nunique() if "sample_item_id" in pk.columns else None
    assign_samples = SS.select_representative_samples(
        pk, n_time=n_time, include_max_tic=include_max_tic)
    out: dict = {"profile": prof, "peaks": pk, "n_samples": n_samples,
                 "assign_samples": assign_samples,
                 "assign_sample_ids": assign_samples["sample_item_id"].tolist()
                 if "sample_item_id" in assign_samples.columns else [],
                 "stages": tuple(stages)}

    if "matrix" in stages:
        mat, bin_mz = TS.build_matrix(pk)
        out["matrix"] = mat
        out["bin_mz"] = bin_mz

    # 'assign' / 'cluster' / 'validate' are added as the scratch modules
    # (cluster_*, isotope_validate) are consolidated into the package against
    # isotopes.py. Flag clearly until then rather than silently no-op. NOTE the
    # 'assign' stage's sample SELECTION is already live (assign_samples above);
    # what remains is to run assign.run per selected id and merge the ledgers.
    todo = [s for s in stages if s in ("assign", "cluster", "validate")]
    if todo:
        out["pending_stages"] = todo
        if n_samples and n_samples < 12 and "cluster" in todo:
            out["cluster_warning"] = (f"{n_samples} samples < 12: temporal clustering "
                                      "unreliable; isotope/coverage layer preferred")

    if out_dir:
        os.makedirs(os.path.expanduser(out_dir), exist_ok=True)
    return out

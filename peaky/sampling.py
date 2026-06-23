"""Representative-sample selection for batch assignment — THE RULE.

A single averaged / single-file spectrum misses peaks that are only present at
certain times: `assign.run` scores ONE sample, so an analyte that spikes during
an event window is never in the candidate list when the chosen file predates the
spike. This bit us on the uronium 24-h run — the 08:20 snapshot sat at hour 11.3,
BEFORE the hour 18-22 event, so the event-only analytes (C10H19NO2, C9H14O2, ...)
were invisible until event files were assigned and merged.

THE RULE (agent-peaky, set 2026-06-19): to get a peak list representative of the
WHOLE batch, assign a small fixed subset and merge it, where the subset is

  * N_TIME (=5) samples evenly spaced across the batch's TIME range — catches
    peaks that appear / disappear over the run; the time endpoints are always
    included, so the start and end of the experiment are covered; PLUS
  * the single MAX-TIC sample — the richest spectrum, the most peaks above the
    detection floor in one shot.

Selecting in TIME (not by row index) matters: an irregularly-sampled run (dense
early, a lone late file, or a gap) must still place a pick on the sparse late
region. The downstream merge (by m/z) of the per-sample assignments is the
representative peak matrix.

Pure pandas/numpy; no network. The selected `sample_item_id`s feed `assign.run`
one at a time (match_compounds is per-sample), then the ledgers are merged.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__version__ = "0.1.0"

N_TIME = 5  # evenly-time-spaced samples per batch (the rule)


def sample_table(peaks: pd.DataFrame, *, sample_col: str = "sample_item_id",
                 time_col: str = "datetime_utc", name_col: str = "sample_item_name",
                 height_col: str = "height", tic_col: str = "tic") -> pd.DataFrame:
    """One row per sample: id, time, name, `tic`, `n_peaks`. Accepts EITHER

      * a per-peak batch table (has `height`) -> tic = sum of peak heights, OR
      * a per-sample table such as `samples.list` (already has a `tic` column,
        one row per sample) -> tic taken directly, no peak load needed.

    Sorted by time when a clock is present. Only `sample_col` is required."""
    if peaks is None or sample_col not in getattr(peaks, "columns", []):
        raise ValueError(f"peaks needs a {sample_col!r} column")
    cols = peaks.columns
    if tic_col in cols and height_col not in cols:
        # already per-sample (e.g. samples.list)
        keep = [c for c in (sample_col, time_col, name_col, tic_col) if c in cols]
        tab = peaks[keep].drop_duplicates(subset=sample_col).reset_index(drop=True)
        tab = tab.rename(columns={tic_col: "tic", time_col: "datetime_utc",
                                  name_col: "sample_item_name"})
    else:
        parts: dict[str, tuple] = {"n_peaks": (sample_col, "size")}
        if height_col in cols:
            parts["tic"] = (height_col, "sum")
        if time_col in cols:
            parts["datetime_utc"] = (time_col, "first")
        if name_col in cols:
            parts["sample_item_name"] = (name_col, "first")
        tab = peaks.groupby(sample_col).agg(**parts).reset_index()
    sort_key = "datetime_utc" if "datetime_utc" in tab.columns else sample_col
    return tab.sort_values(sort_key).reset_index(drop=True)


def select_representative_samples(peaks: pd.DataFrame, *, n_time: int = N_TIME,
                                  include_max_tic: bool = True,
                                  sample_col: str = "sample_item_id",
                                  **table_kw) -> pd.DataFrame:
    """Pick the representative sample subset for assignment (THE RULE).

    Returns the selected rows of `sample_table()`, time-ordered, with an added
    `role` column: 'time-grid', 'max-TIC', or 'time-grid+max-TIC' (when the
    max-TIC sample is also a time-grid pick). Assign each id, then merge by m/z.

    - `n_time` samples are chosen evenly across the batch TIME range: the nearest
      DISTINCT sample to each of `n_time` equally-spaced target times (so the two
      endpoints are always included and dense regions are not over-sampled).
    - The single max-TIC sample is added (union; deduped against the grid).
    - Fewer than `n_time` samples -> all are returned (max-TIC still flagged).
    """
    tab = sample_table(peaks, sample_col=sample_col, **table_kw)
    n = len(tab)
    if n == 0:
        return tab.assign(role=pd.Series(dtype=str))
    has_time = "datetime_utc" in tab.columns
    has_tic = "tic" in tab.columns

    # --- n_time evenly TIME-spaced picks (nearest distinct sample per target) ---
    if not has_time or n <= n_time:
        grid = list(range(n))                       # too few to thin, or no clock
    else:
        t = tab["datetime_utc"].astype("int64").to_numpy()   # ns since epoch
        targets = np.linspace(t[0], t[-1], n_time)
        grid, used = [], set()
        for tg in targets:
            for j in np.argsort(np.abs(t - tg)):     # nearest sample not yet taken
                j = int(j)
                if j not in used:
                    used.add(j)
                    grid.append(j)
                    break
    roles = {i: "time-grid" for i in grid}

    # --- max-TIC (the richest single spectrum) ---
    if include_max_tic and has_tic:
        mi = int(np.asarray(tab["tic"]).argmax())
        roles[mi] = "time-grid+max-TIC" if mi in roles else "max-TIC"

    keep = sorted(roles)
    sel = tab.loc[keep].copy()
    sel["role"] = [roles[i] for i in keep]
    sort_key = "datetime_utc" if has_time else sample_col
    return sel.sort_values(sort_key).reset_index(drop=True)


def select_representative_sample_ids(peaks: pd.DataFrame, *,
                                     sample_col: str = "sample_item_id",
                                     **kw) -> list:
    """Convenience: just the selected `sample_item_id`s (time order)."""
    sel = select_representative_samples(peaks, sample_col=sample_col, **kw)
    return sel[sample_col].tolist()

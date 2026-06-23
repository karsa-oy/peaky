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
import re
from dataclasses import dataclass
from datetime import datetime, timezone

import pandas as pd

from . import io_mascope as IO
from . import profiles as P
from . import sampling as SS
from . import timeseries as TS

__version__ = "0.2.0"  # + representative-sample selection (5 time-grid + max-TIC)

STAGES = ("matrix", "assign", "cluster", "validate")


# ---------------------------------------------------------------------------
# run versioning: every set of outputs goes in its OWN timestamped folder so a
# re-run never overwrites a previous one. The folder name = the run id = the
# batch slug + date + time, and that same id is stamped on the report cover.
# Timestamps are in UTC (explicitly marked Z / "UTC") — the run environment's
# local clock can differ from the user's, so UTC is unambiguous for sharing.
# Pass `when` explicitly (one now() per run) so folder/id/cover all agree.
# ---------------------------------------------------------------------------
def slugify(name: str) -> str:
    """Filesystem-safe batch slug: 'Orange peeling (Ur+ CIMS)' -> 'Orange-peeling-Ur-CIMS'."""
    return re.sub(r"[^0-9A-Za-z]+", "-", str(name)).strip("-") or "run"


def run_stamp(when: datetime | None = None) -> tuple[str, str]:
    """(folder stamp 'YYYY-MM-DDTHHMMSSZ', human 'YYYY-MM-DD HH:MM UTC') in UTC.
    `when`: None -> now (UTC); a tz-aware datetime is converted to UTC; a naive
    datetime is assumed to already be UTC."""
    if when is None:
        when = datetime.now(timezone.utc)
    elif when.tzinfo is not None:
        when = when.astimezone(timezone.utc)
    return when.strftime("%Y-%m-%dT%H%M%SZ"), when.strftime("%Y-%m-%d %H:%M UTC")


def run_id(batch_name: str, when: datetime | None = None) -> str:
    """Versioned run id = batch slug + timestamp; also the output folder's name."""
    return f"{slugify(batch_name)}_{run_stamp(when)[0]}"


def make_run_dir(base: str, batch_name: str, when: datetime | None = None) -> str:
    """Create and return a fresh per-run output folder `base/<run_id>`."""
    d = os.path.join(os.path.expanduser(base), run_id(batch_name, when))
    os.makedirs(d, exist_ok=True)
    return d


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


# ---------------------------------------------------------------------------
# RunContext + the batch pipeline (assign -> cluster -> Van Krevelen -> report).
# Replaces the out-of-repo run_peaky.py orchestrator and its PEAKY_* env-var +
# subprocess threading: ONE object carries the run's identity (out_dir / run_id /
# generated stamp / batch / profile) into every stage as a plain argument.
# ---------------------------------------------------------------------------
@dataclass
class RunContext:
    out_dir: str                 # the versioned run folder (holds all artifacts)
    batch_name: str
    tag: str                     # filename token (e.g. the profile name)
    label: str                   # display label for figure/report titles
    when: datetime
    run_id: str                  # == folder basename; stamped on the report cover
    generated: str               # human stamp 'YYYY-MM-DD HH:MM UTC'
    profile: object = None       # ReagentProfile
    dataset: str | None = None
    ts_path: str | None = None   # parquet the report reads for the event TIC


def make_run_context(base_out: str, batch_name: str, profile, *, when=None,
                     tag=None, label=None, dataset=None) -> RunContext:
    """Create a fresh timestamped run folder and the context that identifies it.
    Pass ONE `when` per run so folder / run_id / cover stamp all agree."""
    when = when or datetime.now(timezone.utc)
    return RunContext(
        out_dir=make_run_dir(base_out, batch_name, when), batch_name=batch_name,
        tag=tag or profile.name, label=label or profile.label, when=when,
        run_id=run_id(batch_name, when), generated=run_stamp(when)[1],
        profile=profile, dataset=dataset)


def generate_report(ctx: RunContext, ts, *, subject: str | None = None,
                    do_cluster=True, do_vk=True, do_report=True, log=print) -> dict:
    """Offline generation half: cluster figures + Van Krevelen + the PDF report,
    from the merged / per-file ledgers already in `ctx.out_dir` and the batch time
    series. No network. Same artifacts the run_clusters / run_vankrevelen /
    run_report scratch chain produced, but in-process via the RunContext."""
    from . import analyte_viz as V
    from . import clustering as CLU
    from . import pdf_report as R

    if isinstance(ts, str):
        ctx.ts_path = os.path.expanduser(ts)
        ts = pd.read_parquet(ctx.ts_path)
    elif ctx.ts_path is None:
        ctx.ts_path = os.path.join(ctx.out_dir, f"{ctx.tag}_ts.parquet")
        ts.to_parquet(ctx.ts_path)

    out: dict = {"ctx": ctx}
    if do_cluster:
        out["cluster"] = CLU.cluster_batch(ctx.out_dir, ts, ctx.profile,
                                           tag=ctx.tag, label=ctx.label, log=log)
    if do_vk:
        out["vk"] = V.van_krevelen_batch(ctx.out_dir, ts, ctx.profile, tag=ctx.tag,
                                         label=ctx.label, batch_name=ctx.batch_name,
                                         subject=subject, log=log)
    if do_report:
        out["report_pdf"] = R.build(ctx.out_dir, tag=ctx.tag, label=ctx.label,
                                     ts_path=ctx.ts_path, batch_name=ctx.batch_name,
                                     run_id=ctx.run_id, generated=ctx.generated)
        log(f"[report] wrote {out['report_pdf']}")
        # also emit a size-reduced companion for emailing (optional deps; no-op if
        # absent or already small). The full report above is left untouched so it
        # stays byte-for-byte deterministic.
        small = R.compress_pdf(out["report_pdf"], log=log)
        if small:
            out["report_pdf_small"] = small
            log(f"[report] compressed -> {small} ({os.path.getsize(small) / 1e6:.1f} MB)")
    return out


def run_batch(*, batch: str, dataset: str | None = None, reagent: str = "auto",
              base_out: str, ts=None, when=None, subject: str | None = None,
              amine_r_min: float = 0.7, do_report=True, config: str | None = None,
              log=print, **assign_kw) -> dict:
    """Full batch pipeline in ONE call: representative-sample ASSIGN (live
    match_compounds) -> merge -> cluster figures -> Van Krevelen -> PDF report,
    into one versioned run folder. `ts` is the full-batch per-sample peak time
    series (DataFrame or parquet path); if None it is fetched live and reused for
    the amine gate + clustering. Returns {ctx, assign, cluster, vk, report_pdf}."""
    from . import assign_batch as AB

    if isinstance(ts, str):
        ts = pd.read_parquet(os.path.expanduser(ts))
    if ts is None:
        log(f"[batch] fetching full-batch time series for {batch!r} ...")
        ts = load(batch=batch, dataset=dataset)
    prof = P.resolve(reagent, ts, config=config)
    ctx = make_run_context(base_out, batch, prof, when=when, dataset=dataset)
    log(f"[batch] {ctx.run_id} -> {ctx.out_dir}")

    res = AB.run(batch=batch, dataset=dataset, reagent=prof.name,
                 out_dir=ctx.out_dir, ts_peaks=ts, amine_r_min=amine_r_min,
                 log=log, **assign_kw)
    gen = generate_report(ctx, ts, subject=subject, do_report=do_report, log=log)

    # provenance: pin this run to its exact code + input-data hash + config +
    # output hash, and append it to the cross-run registry. Best-effort (never
    # fatal). Runs LAST so ts_path / merged_ledger.csv exist to be hashed.
    from . import passes as PA
    from . import provenance as PV
    summ = res.get("summary", {}) if isinstance(res, dict) else {}
    PV.record_run(
        run_dir=ctx.out_dir, base_out=os.path.expanduser(base_out),
        batch_name=batch, dataset=dataset,
        sample_ids=(res.get("sample_ids") if isinstance(res, dict) else None),
        reagent=prof.name, cfg=assign_kw.get("cfg") or PA.PassConfig(),
        ts_path=ctx.ts_path,
        counts={"merged_M0": summ.get("merged_M0"),
                "merged_tiers": summ.get("merged_tiers"),
                "n_samples": summ.get("n_files")},
        created_utc=ctx.when.isoformat(), log=log)
    return {"ctx": ctx, "assign": res, **gen}

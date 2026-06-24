"""peaky — multi-pass chemical-formula assignment for high-resolution
mass-spec peaks stored in Mascope.

SDK-native, test-driven. Mascope's ``match_compounds`` is the scoring oracle; this
package owns candidate generation, chemistry plausibility, series logic and
arbitration. State lives in one mutable ledger DataFrame (one row per peak);
passes are ledger -> ledger functions that only fill/annotate.

Public API (import from the package root; internals may move):

    import peaky as ma

    ma.run(sample_id, context)        # single-sample assignment -> dict
    ma.run_batch(batch=..., ...)      # FULL batch pipeline: assign -> cluster -> VK -> PDF report
    ma.run_assign_batch(batch=..., ...) # the assign + merge half only (no figures/report)
    ma.PassConfig(...)                # assignment knobs
    ma.get_context("ambient-air")     # context (plausibility + contaminant families)
    ma.resolve_reagent("auto", peaks) # ReagentProfile (Br / Ur / NO3 ...)
    ma.build_report(out_dir, ...)     # standard PDF assignment report

Attributes are resolved lazily (PEP 562) so ``import peaky`` stays cheap
and does not pull matplotlib or the Mascope SDK until a heavy entry point is used.
"""
from __future__ import annotations

__version__ = "0.5.0"

# public name -> (submodule, attribute)
_LAZY = {
    "run": ("assign", "run"),                       # single-sample assignment
    "run_batch": ("pipeline", "run_batch"),         # FULL pipeline (assign->cluster->VK->report)
    "run_assign_batch": ("assign_batch", "run"),    # assign + merge half only
    "run_pipeline": ("pipeline", "run_batch"),      # back-compat alias of run_batch
    "PassConfig": ("passes", "PassConfig"),
    "get_context": ("contexts", "get_context"),
    "resolve_reagent": ("profiles", "resolve"),
    "ReagentProfile": ("profiles", "ReagentProfile"),
    "build_report": ("pdf_report", "build"),
}

__all__ = ["__version__", *sorted(_LAZY)]


def __getattr__(name: str):
    import importlib

    if name in _LAZY:
        mod, attr = _LAZY[name]
        return getattr(importlib.import_module(f".{mod}", __name__), attr)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(__all__)

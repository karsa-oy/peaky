"""Run-folder layout — the single source of truth for where a batch run's
artifacts live, shared by the WRITERS (assign_batch / clustering / analyte_viz /
pdf_report) and the READER (pdf_report.load_context). Both sides derive paths from
the same ``RunPaths(out_dir)``, so the writer/reader filename contract can never
drift: change a subdir here and both ends move together.

Layout under a run folder::

    <run>/                      run root
      merged_ledger.csv         the result (provenance anchor)   — ROOT
      run_manifest.json         reproducibility manifest         — ROOT
      batch_summary.json        run counts/offsets               — ROOT
      per_file/                 per-sample ledgers               — ROOT (own subdir)
      figures/                  all .png (cluster panels, GKA, Van Krevelen)
      tables/                   all .csv / .xlsx (cluster tables, jitter, channel QC)
      report/                   the PDF report (+ compressed companion)
      data/                     bulky inputs kept with the run (e.g. a fetched TS)

ROOT-level names stay flat on purpose: they are read by several modules + the
cross-run registry, so moving them is high blast radius for no gain.
"""
from __future__ import annotations

import os

# The installed package root (the `peaky/` directory). Resolve bundled resources
# (e.g. `data/peaklists/`) and the workspace root from HERE, not from each calling
# module's `__file__`, so modules stay correct regardless of which sub-package they
# live in. `paths.py` itself stays at the package root, so this anchor is stable.
PKG_ROOT = os.path.dirname(os.path.abspath(__file__))


def pkg_data(*parts: str) -> str:
    """Absolute path to a bundled data resource under `peaky/data/`."""
    return os.path.join(PKG_ROOT, "data", *parts)


# Names that intentionally stay at the run root (not routed into a subdir).
ROOT_ANCHORS = frozenset({"merged_ledger.csv", "run_manifest.json", "batch_summary.json"})


class RunPaths:
    """Subdir paths for one run folder. Attributes are absolute; ``ensure()``
    creates the subdirs and returns self."""

    __slots__ = ("root", "figures", "tables", "report", "data", "per_file")

    def __init__(self, out_dir: str):
        self.root = os.path.expanduser(out_dir)
        self.figures = os.path.join(self.root, "figures")
        self.tables = os.path.join(self.root, "tables")
        self.report = os.path.join(self.root, "report")
        self.data = os.path.join(self.root, "data")
        self.per_file = os.path.join(self.root, "per_file")

    def ensure(self) -> "RunPaths":
        for d in (self.figures, self.tables, self.report, self.data):
            os.makedirs(d, exist_ok=True)
        return self

    def place(self, filename: str) -> str:
        """Route a BARE filename to its subdir by role/extension (root anchors stay
        at root). Use this for single named files; use the subdir attributes
        directly for glob patterns and render-prefixes."""
        if filename in ROOT_ANCHORS:
            return os.path.join(self.root, filename)
        ext = filename.rsplit(".", 1)[-1].lower() if "." in filename else ""
        if ext == "png":
            return os.path.join(self.figures, filename)
        if ext in ("csv", "xlsx"):
            return os.path.join(self.tables, filename)
        if ext == "pdf":
            return os.path.join(self.report, filename)
        return os.path.join(self.root, filename)


def run_paths(out_dir: str) -> RunPaths:
    """`RunPaths` for a run folder (does NOT create subdirs; call `.ensure()`)."""
    return RunPaths(out_dir)

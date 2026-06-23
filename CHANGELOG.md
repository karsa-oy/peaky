# Changelog

All notable changes to Peaky are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] — 0.4.0 (public-release refactor)

A refactor pass preparing Peaky for the public `karsa-oy/peaky` repo: cleaner
install, enforced reproducibility, organized outputs, and a design doc.

### Added
- `docs/ARCHITECTURE.md` — the canonical design doc (ledger model, pass sequence,
  end-to-end data flow with diagram, reproducibility model, module map).
- `CHANGELOG.md` (this file).
- Legacy workspace-based Mascope server support (`io_mascope`): connects to older
  deployments where `/api/datasets` 404s, resolving workspaces/batches via the raw
  endpoints. Additive and gated — modern servers are unaffected.

### Changed
- Moved `ROADMAP.md` → `docs/ROADMAP.md` (kept as development history); README now
  points at `docs/ARCHITECTURE.md` as the entry point for how Peaky works.

<!-- Filled in as the remaining phases land:
### Changed (packaging)   — import package renamed mascope_assign → peaky (back-compat shim); single canonical lock.
### Fixed (reproducibility) — run driver now exports SOURCE_DATE_EPOCH from the run time so figures/PDF are byte-stable.
### Changed (outputs)     — run dir organized into figures/ tables/ report/; input time-series no longer copied per run.
-->

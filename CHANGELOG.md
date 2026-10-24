# Changelog

All notable changes to Tessera. Software versions follow semver; the science model has its own
version (see `docs/science.md`).

## [1.0.1] - 2026-10-24

### Fixed
- README images use absolute URLs so they render on PyPI (the 1.0.0 page showed broken images).
- Removed a leftover, unused partial copy of the AT2017gfo photometry from the package.

### Added
- Releases publish to PyPI through trusted publishing (`release.yml`, environment `pypi`),
  uploading the exact files attached to the GitHub release.

## [1.0.0] - 2026-10-24

First public release. Science model 1.1.0.

### Science core (`tessera.core`)
- Sky-map ingestion for flat (RING/NESTED, gzip) and multi-order LVK maps with strict validation
  and resource limits; `p_sky` and credible levels agree with `ligo.skymap` crossmatch to 1e-12.
- Galaxy association: `p_sky`, `p_dist` (closed-form marginalisation over galaxy distance error),
  astrophysical weights (K-band, B-band, stellar mass, flat, GLADE+ BNS rate), tiers and data-quality
  flags. Ranking score `S`; opt-in `host-posterior` strategy with an explicit uncataloged term.
- Catalog probability budget from a smoothed completeness cube (`P_cat`, `P_uncat`, footprint,
  depth), with versioned Parquet catalog snapshots indexed by HEALPix.
- Field generation (grid, tessellation, galaxy-targeted, hybrid), astroplan visibility, greedy-v1
  scheduler with decision log, optional filter blocks, and an independent 10 s feasibility validator.
- MILP upper bound to measure scheduler optimality (greedy reaches ≥ 0.91 on GW170817).
- Reproducibility manifests, `tessera rerun`, and `tessera compare` for independent validation.
- AT2017gfo light curve from Villar et al. 2017; reference densities from Kopparapu et al. 2008,
  Kochanek et al. 2001 and Bell et al. 2003.

### Service
- FastAPI REST API (`/v1`) with API keys, in-process OIDC, and trusted-gateway headers.
- Database job queue with isolated worker processes, wall-time and memory limits, retries.
- LVK alert listener (GCN Kafka, SCiMMA, GraceDB, directory) with revision-aware replanning,
  idempotency and retraction handling.
- Versioned schedule edits, approvals bound to an exact version, gated submission, audit log.
- Scheduler adapters: generic, ACP (unverified), LCO, Treasure Map, Rubin ToO footprint.
- Web UI, Prometheus metrics, structured logs, Docker Compose and a Helm chart.

### Known limitations
- ACP plan syntax is not yet confirmed against vendor documentation.
- LCO instrument codes and submission endpoint, Treasure Map cancellation, and Rubin's
  operational ToO triggers must be confirmed with each facility before submission is enabled.
- The stellar-mass reference density is an upper limit; its IMF consistency with GLADE+ masses is
  unverified. Kilonova diversity σ_M = 0.5 mag is an assumption.
- Packaged telescope profiles are illustrative, not facility-validated.
- The two-group independent rerun (docs/validation-protocol.md) has not yet been performed.

[1.0.1]: https://github.com/bibinthomas123/Tessera/releases/tag/v1.0.1
[1.0.0]: https://github.com/bibinthomas123/Tessera/releases/tag/v1.0.0

# Changelog

All notable changes to Tessera. Software versions follow semver; the science model has its own
version (see `docs/science.md`).

## [Unreleased]

### Fixed
- Hybrid field mode no longer leaves the telescope idle once the galaxy fields are done. The rest
  of the credible region is added as low-priority `fill` tiles, scheduled only when no galaxy or
  low-coverage field has gain. GW170817 / example 1 m: 25 → 69 visits, sky probability 16% → 72%,
  galaxy score and the first 25 visits unchanged. Found by comparing with gwemopt and tilepy.

- Very large maps no longer lose all their fields. The low-value cut (`fields.low_value_eps`, 1e-4)
  now also requires a field to be below 1% of the best field (`fields.low_value_rel`). Before, every
  tile of a map spread over tens of thousands of deg² fell under the absolute cut: S230529ay with
  the example 1 m telescope got 0 tiling visits, and now gets 70; the default mode now gets 0.0018
  expected detections. When the `max_fields` cap binds, tiles the site can never observe are
  dropped before the cap (`FIELDS_UNREACHABLE_FROM_SITE`).

### Changed
- Exposure times under the yield objective: new `scheduler.exposure_choice` (`best` by default,
  `yield_rate`, `depth`). `best` schedules with depth-solved exposures and with exposures chosen
  for detections per hour, and keeps the schedule with more expected detections, so it is never
  worse than before. Example 1 m: S190425z 0.013 → 0.017 and S200105ae 0.006 → 0.012 expected
  detections. The wide-field example keeps depth. Yield schedules take about twice as long to
  compute. `depth` restores the 1.1.0 behaviour.

## [1.1.0] - 2026-10-24

Science model 1.2.0.

### Changed
- Stellar-mass reference density is 2.3×10⁸ M☉ Mpc⁻³ (Baldry et al. 2012, GAMA, Chabrier IMF),
  on the same mass scale as GLADE+ stellar masses (W1 mass-to-light ratios of Kettlety et al. 2018,
  calibrated on GAMA). It was the Bell et al. 2003 upper limit, 3.9×10⁸. No longer provisional.
- Kilonova diversity σ_M for `at2017gfo` is 1.0 mag, sourced from Ascenzi et al. 2019 (was an
  assumed 0.5 mag). Expected yields of distant events go down; schedules change little.

### Added
- `tessera benchmark` reports the sensitivity of schedules and yields to σ_M (0.25 / 0.5 / 1.0).
- The LCO adapter checks instrument types, telescope classes, filter codes and guiding/acquisition
  modes against a dated snapshot of LCO's public instrument list, with `refresh_instruments()`.

### Known limitations
- Still open: ACP syntax, the LCO submission endpoint, Rubin's operational ToO triggers, facility
  validation of profiles, and the two-group rerun.

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

[1.1.0]: https://github.com/bibinthomas123/Tessera/releases/tag/v1.1.0
[1.0.1]: https://github.com/bibinthomas123/Tessera/releases/tag/v1.0.1
[1.0.0]: https://github.com/bibinthomas123/Tessera/releases/tag/v1.0.0

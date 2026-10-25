# Changelog

All notable changes to Tessera. Software versions follow semver; the science model has its own
version (see `docs/science.md`).

## [Unreleased]

### Added
- `ztf-approx` telescope profile built from published values: the 4 × 4 CCD mosaic with chip gaps,
  g/r/i depths, 8 s readout and 10.2 s median overhead (Bellm et al. 2019), drive speeds
  (Dekany et al. 2020), and observing limits from ZTF's scheduler. Planning with it is slower (S190425z: about 50 s) because each pointing is rasterised
  as 16 chips.

- The GLADE+ reader also reads the full fixed-width CDS file `gladep.dat` from VizieR (VII/291,
  23.2 million rows; about 3 minutes to convert), detected automatically.

- Mount `ha_limits_by_dec`: hour-angle windows that apply only south of a given declination,
  enforced in visibility, the fast refinement and the validator.

### Fixed
- `ztf-approx` used rules stricter than ZTF's own scheduler: declinations south of −22° are now
  allowed down to −45° with HA ≥ −17.6° (ztf_sim QueueManager), and the night uses 12° twilight
  (ztf_sim Scheduler.py) instead of astronomical twilight. Both were found while checking the
  comparison harness, and would otherwise have penalised other tools' valid ZTF visits.
- The scheduler ended the night early when nothing observable *at that moment* had gain, even if
  fields rose later. It now jumps to the next window opening, or the next due revisit, and stops
  only when none remain. Example: S240422ed on the wide-field example observed until 04:12 UTC
  and left 7.9 dark hours unused; it now observes until 12:07, 54 → 365 visits, sky probability
  0.757 → 0.791, 89% → 93% of the MILP upper bound. In the 8-event comparison, 5 cases improve and
  none get worse. Found while measuring the gap to the upper bound.
- `optimality.simple_upper_bound`: a provably valid upper bound for any number of fields (union of
  observable coverage, and the N best single fields), reported by `greedy_vs_bound` next to the
  MILP bound.
- A start time the server could not parse (for example `2017-08-17 22:00 UTC`) caused an HTTP 500.
  Times now accept a space instead of `T`, missing seconds, and a trailing `Z` or `UTC`; anything
  else is rejected with HTTP 422 and code `INVALID_TIME`.

### Changed
- Web UI: a running plan's progress (Uploading…, Queued…, Running · stage) is shown inside the Plan
  button, which is disabled until the plan finishes; messages and errors appear below it.

## [1.3.0] - 2026-10-24

Scheduling improvements from the comparison with gwemopt and tilepy. Science model 1.2.0 (unchanged).

### Changed
- Expected detections use the best look per pixel instead of treating repeated or overlapping looks
  as independent chances (they are not: a counterpart too faint for one exposure stays too faint).
  Reported yields of revisit patterns go down; the GW170817 golden schedule is unchanged.
- Half-tile-shifted grids (plus quarter shifts for fields ≥ 2°) are added as lowest-priority
  `offset` fields (`fields.offset_fill`), and the scheduler compares complete schedules with and
  without them (and, for the yield objective, with each exposure policy), keeping the best by the
  plan's objective (`summary.variants_compared`).
- In the comparison with gwemopt and tilepy (8 events, 2 telescopes, 30 s exposures for all
  tools), Tessera now gives the most expected detections in 14 of 14 observable cases and the most
  sky probability in 13 (tilepy leads S250206dm wide-field by 0.001). With the default settings,
  wide-field expected detections rise, for example S240422ed 0.48 → 0.59 and S190425z 0.45 → 0.51.

- Faster visibility windows: edge bisection uses Sun and Moon positions interpolated from the
  coarse grid, and every edge is then re-checked with the exact constraints. This is 2.5–3.4×
  faster, with identical windows to within one refinement step.
- Offset variants run only when the plain schedule had usable idle time, and offset fields'
  visibility is computed only then. Planning time is back near 1.2.0: S190425z with two
  telescopes 18 s (1.2.0: 18 s), S200105ae 14 s (14 s), GW170817 on the 1 m 4 s (2 s).

### Known regressions
- With the wide-field example's g-band revisit pattern, expected detections are slightly lower
  than in 1.2.0 when scored by the same independent evaluator: S190425z 0.454 → 0.452 and
  S200105ae 0.453 → 0.444. The greedy choices under the corrected accounting land slightly worse
  there.

## [1.2.0] - 2026-10-24

Scheduling fixes found by comparing Tessera with gwemopt and tilepy on 8 events (3 historical,
5 public O4 alerts) and two telescope classes. Science model 1.2.0 (unchanged).

### Fixed
- Hybrid field mode no longer leaves the telescope idle once the galaxy fields are done. The rest
  of the credible region is added as low-priority `fill` tiles, scheduled only when no galaxy or
  low-coverage field has gain. GW170817 / example 1 m: 25 → 69 visits, sky probability 16% → 72%,
  galaxy score and the first 25 visits unchanged.
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

[1.3.0]: https://github.com/bibinthomas123/Tessera/releases/tag/v1.3.0
[1.2.0]: https://github.com/bibinthomas123/Tessera/releases/tag/v1.2.0
[1.1.0]: https://github.com/bibinthomas123/Tessera/releases/tag/v1.1.0
[1.0.1]: https://github.com/bibinthomas123/Tessera/releases/tag/v1.0.1
[1.0.0]: https://github.com/bibinthomas123/Tessera/releases/tag/v1.0.0

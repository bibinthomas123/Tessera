# Independent operational validation protocol (spec 11.3)

Goal: two or more groups, on different institutions' hardware, run the **same** inputs and show
that tessera's outputs agree exactly, or within the documented tolerance, with every difference
explained.

## Inputs (fixed by the coordinating group)

| Input | How it is pinned |
| --- | --- |
| Sky map | file + sha256 (e.g. the GWTC-1 `GW170817_skymap.fits.gz` in `tests/data`) |
| Catalog snapshot | `snapshot_id` (sha256 of `snapshot.json`); share the snapshot directory or the ingestion command + upstream file checksum |
| Telescope profiles | `id@version` + profile sha256 |
| Config | the `request_config` block of a reference manifest |
| Planning start | `t_start` (ISO UTC) |
| Software | container image digest **and** a source install of the same tag |

## Procedure

1. Each group runs, in both the published container and a source install:

   ```bash
   tessera plan --skymap MAP --catalog SNAPSHOT --profile example-1m --config config.yaml \
                --t-start 2017-08-17T22:00:00 --out out/<group>-<container|source>
   ```

2. Groups exchange their `out/` directories (they contain `manifest.json` with every input hash).
3. Each pair is compared:

   ```bash
   tessera compare out/groupA-container out/groupB-container --out compare-AB.json
   ```

   Verdicts: `identical` (every output hash equal), `within_tolerance` (same inputs; top-100
   galaxies in the same order except ties; scores within a relative 1e-9; same field sequence;
   visit times within 1 s), or `different`.
4. Every `different` result is triaged using `likely_causes`: input differences, software
   versions, IERS table date, ephemeris, platform. A difference is closed only when a rerun with
   the cause removed gives `identical` or `within_tolerance`.
5. The comparison reports and conclusions go into the benchmark report for the release.

## Acceptance (v1.0 gate)

* Same container, any host: `identical`.
* Source installs on different platforms: `within_tolerance`.
* Any `different` verdict has a documented cause and fix.

This protocol needs at least one external group; tessera provides the tooling but cannot
perform the cross-institution run itself.

# Roadmap

| Phase | Content | Gate |
| --- | --- | --- |
| 1 Offline engine | sky map + snapshot + profile → galaxies, fields, schedules, manifest | crossmatch agreement, injection tests |
| 2 API and UI | authenticated /v1, Postgres queue, Docker Compose, observability | latency and fuzz targets |
| 3 Live alerts, v1.0 | GCN/SCiMMA/GraceDB, revisions, retractions, adapters | O4 alert replay clean, two-group rerun |
| 4 Coordination, v1.1 | multi-facility allocation, observation feedback, candidate association | — |

Status: v1.0.0 released (phases 1–3 implemented); phase 4 has sequential coordination, feedback replanning and
transient association. Done since: tabulated AT2017gfo light curve, sourced reference densities,
greedy-vs-MILP upper bound, filter-block sequencing, in-process OIDC, Helm chart, `tessera compare`
and the validation protocol, Rubin/LCO/Treasure Map format checks.

Open items:

* ACP directive syntax: the vendor documentation could not be reached; confirm with an ACP site.
* LCO instrument types, modes and filters (configdb) and the submission endpoint, with LCO.
* Rubin ToO operational trigger criteria and transport, with Rubin's scheduler team.
* Stellar-mass reference density: IMF consistency with GLADE+ masses.
* The two-group independent rerun itself (needs an external group; see docs/validation-protocol.md).

# Roadmap

| Phase | Content | Gate |
| --- | --- | --- |
| 1 Offline engine | sky map + snapshot + profile → galaxies, fields, schedules, manifest | crossmatch agreement, injection tests |
| 2 API and UI | authenticated /v1, Postgres queue, Docker Compose, observability | latency and fuzz targets |
| 3 Live alerts, v1.0 | GCN/SCiMMA/GraceDB, revisions, retractions, adapters | O4 alert replay clean, two-group rerun |
| 4 Coordination, v1.1 | multi-facility allocation, observation feedback, candidate association | — |

Status: v1.3.0 released (phases 1–3 implemented); phase 4 has sequential coordination, feedback replanning and
transient association. Done since: tabulated AT2017gfo light curve, sourced reference densities,
greedy-vs-MILP upper bound, filter-block sequencing, in-process OIDC, Helm chart, `tessera compare`
and the validation protocol, Rubin/LCO/Treasure Map format checks, GAMA-matched stellar-mass density,
sourced kilonova diversity with a sensitivity benchmark, LCO codes checked against LCO's instrument list.

Open items:

* ACP directive syntax: the vendor documentation could not be reached; confirm with an ACP site.
* LCO submission endpoint and token header, via LCO's validate endpoint with a real account
  (instrument types, filters and modes are now checked against LCO's public instrument list).
* Rubin ToO operational trigger criteria and transport, with Rubin's scheduler team.
* The two-group independent rerun itself (needs an external group; see docs/validation-protocol.md).

# tessera — gravitational-wave follow-up engine

An open-source, telescope-agnostic engine that turns an LVK HEALPix sky map (flat or multi-order)
into:

1. **Ranked host galaxies.** Every score component is exposed (`p_sky`, `p_dist`, `w`, credible
   level, distance z-score, data-quality flags). `S` is a *ranking score*, not a host probability.
   The opt-in `host-posterior` strategy gives normalised probabilities under a stated generative
   model.
2. **A catalog probability budget.** It reports how much probability the catalog can explain and how
   much is uncataloged, outside the footprint, or beyond reliable depth.
3. **Feasible schedules.** Astroplan hard constraints, per-axis slew and filter-change overheads, a
   deterministic greedy optimiser with a decision log, and an independent validator. Infeasible
   schedules cannot be exported.
4. **A reproducibility manifest.** `plan_key`, input checksums, catalog snapshot ID, profile
   hashes, full config, software and IERS versions, and output hashes. `tessera rerun` diffs a rerun
   against it.

The core (`tessera.core`) is pure Python with no network, database or web imports, so pipelines can
import it directly. The service adds a REST API, a Postgres job queue, isolated workers, an LVK alert
listener with revision-aware replanning, a schedule edit/approve workflow, scheduler adapters (ACP,
LCO, Treasure Map, Rubin ToO footprint) and a small web UI.

## Quick start (library / CLI)

```bash
python -m venv .venv && .venv/bin/pip install -e ".[dev]"
tessera ingest-catalog glade+ tests/data/glade_plus_gw170817.tsv --out catalogs/glade-excerpt \
     --upstream-version 2022-10-31
tessera plan --skymap tests/data/GW170817_skymap.fits.gz --catalog catalogs/glade-excerpt \
     --profile example-1m --preset bns --nights 1 --t-start 2017-08-17T22:00:00 --out out/gw170817
tessera validate-schedule "out/gw170817/schedules/example-1m@1.json"
tessera export "out/gw170817/schedules/example-1m@1.json" --adapter acp
```

```python
from tessera import plan, PlanConfig
r = plan("bayestar.multiorder.fits", "catalogs/gladeplus-2022.10", ["example-1m"],
         PlanConfig(preset="bns"), t_start="2017-08-17T22:00:00")
r.galaxies[:10]; r.budget["p_cat"]; r.schedules["example-1m@1"]["visits"]; r.write("out/")
```

## Service

`docker compose up` starts api, worker, listener and Postgres. The UI is at `/ui/` and the OpenAPI
docs at `/v1/docs`. See [docs/deployment.md](docs/deployment.md) and [docs/api.md](docs/api.md).
Schedules stay `proposed` until a person approves an exact version.

## Validation

`tessera benchmark` produces [docs/benchmark-report.md](docs/benchmark-report.md). Current results:

- **GW170817:** NGC 4993 ranks 1st with `vol3d-lumK` (14th with sky position only). The example 1 m
  profile observes its field on the 6th visit, 10.7 h after merger (SSS17a was found at ~10.9 h).
- **Synthetic injections:** top-10 host recovery is 0.43 with 3D scoring vs 0.07 for the 2D
  baseline (200 injections).
- **Budget accuracy:** for a randomly halved catalog, the predicted completeness (relative to the
  full catalog) is within 0.01 of the truth. Absolute P_cat also absorbs local over/under-density
  relative to the cosmic mean; see the report.
- **Scheduler optimality:** greedy reaches ≥ 0.91 of a provable MILP upper bound on GW170817.
- **Engineering:** every schedule passes the validator, outputs are deterministic, malformed and
  oversized inputs are rejected with typed errors, and edge-case maps (multi-modal, poles, RA = 0,
  all-sky, no distance) are handled.
- **Reference check:** p_sky and credible levels agree with `ligo.skymap` crossmatch to 1e-12.

`pytest` runs 96 tests, including service integration and a GW170817 golden regression.

## Installing

The distribution is published as **`tessera-gw`** (the PyPI name `tessera` belongs to an unrelated
Graphite dashboard); the import package and the command are `tessera`. Do not install both in one
environment. Extras: `tessera-gw[service]`, `[postgres]`, `[alerts]`, `[rubin]`, `[validation]`.

Kubernetes: `helm install tessera deploy/helm/tessera` (see [docs/deployment.md](docs/deployment.md)).
Independent two-group validation: [docs/validation-protocol.md](docs/validation-protocol.md) and
`tessera compare`.

## Data and caveats

- **Test data** is public: `tests/data` holds the GWTC-1 GW170817 map (LIGO DCC P1800381), the
  GraceDB BAYESTAR maps for S190425z and S200105ae, and a GLADE+ excerpt from VizieR VII/291 (Dalya
  et al. 2022; cite when used). The AT2017gfo light curve is derived from Villar et al. 2017, Table 3.
- **Telescope profiles** packaged here are *illustrative or approximate*, not facility-validated.
- **Science defaults (model 1.1.0):** the AT2017gfo light curve and the B- and K-band densities are
  sourced; the stellar-mass density is an upper limit whose IMF match to GLADE+ masses is unverified,
  and σ_M = 0.5 mag kilonova diversity is an assumption.
- **Adapter formats:** Rubin is checked against `rubin-scheduler`; LCO against the open-source
  Observation Portal; Treasure Map against its official client; ACP is unverified. Each facility must
  still confirm before submission is enabled. See [docs/adapters.md](docs/adapters.md).

Science reference: [docs/science.md](docs/science.md). Contributing and SMEP process:
[CONTRIBUTING.md](CONTRIBUTING.md). Roadmap: [ROADMAP.md](ROADMAP.md). Licence: BSD-3-Clause.

# Tessera

[![CI](https://github.com/bibinthomas123/Tessera/actions/workflows/ci.yml/badge.svg)](https://github.com/bibinthomas123/Tessera/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)
![License](https://img.shields.io/badge/license-BSD--3--Clause-green)

**Tessera is an open-source, telescope-agnostic engine for gravitational-wave follow-up.** It turns an
LVK HEALPix sky map into ranked host galaxies and feasible, validated observing schedules for your
telescopes. Every result comes with a reproducibility manifest.

![Tessera planning GW170817: credible region, scheduled fields, top galaxies, probability budget and strategy comparison](https://raw.githubusercontent.com/bibinthomas123/Tessera/main/docs/images/ui-plan-overview.png)

*GW170817 planned with the public GWTC-1 sky map, a GLADE+ excerpt, and an example 1 m telescope
profile. NGC 4993 ranks first, and its field is observed on the sixth visit, 10.7 h after merger.*

## What Tessera does

1. **Ranks host galaxies transparently.** Every score component is shown: sky density `p_sky`, the
   distance term `p_dist`, the astrophysical weight `w`, the credible level, the distance z-score and
   data-quality flags. `S` is a *ranking score*, not a host probability. The opt-in `host-posterior`
   strategy gives normalised probabilities under a stated generative model.
2. **Reports what the catalog cannot explain.** The probability budget shows how much is cataloged
   (`P_cat`), how much is uncataloged, how much lies outside the footprint, and how much lies beyond
   reliable depth.
3. **Builds schedules that can actually be executed.** Hard constraints come from astroplan:
   twilight, airmass, horizon mask, Moon, hour-angle limits and operational windows. Slews and filter
   changes are modelled per axis. A deterministic greedy optimiser records why each visit was chosen,
   and an independent validator rechecks every visit at 10 s resolution. Infeasible schedules cannot
   be exported.
4. **Makes runs reproducible.** Each run records a `plan_key`, input checksums, the catalog snapshot
   ID, profile hashes, the full config, and software and IERS versions. `tessera rerun` and
   `tessera compare` diff any two runs.

## Screenshots

<table>
<tr>
<td width="50%"><b>Start a plan</b>: upload a sky map (or give a GraceDB URL), pick a catalog, preset and telescopes.<br><img src="https://raw.githubusercontent.com/bibinthomas123/Tessera/main/docs/images/ui-new-plan.png" alt="New plan form"></td>
<td width="50%"><b>REST API</b>: versioned <code>/v1</code> endpoints with OpenAPI docs.<br><img src="https://raw.githubusercontent.com/bibinthomas123/Tessera/main/docs/images/api-docs.png" alt="OpenAPI documentation"></td>
</tr>
</table>

**Ranked host galaxies**, with every score component and flag (downloadable as CSV, ECSV, JSON,
Parquet or VOTable):

![Ranked galaxy table](https://raw.githubusercontent.com/bibinthomas123/Tessera/main/docs/images/ui-galaxies.png)

**A validated schedule**, with the decision behind each visit and a per-visit feasibility check.
The schedule stays `proposed` until a person approves one exact version. It exports to generic
JSON/CSV, ACP, LCO, Treasure Map, or a Rubin ToO package:

![Schedule with decision log](https://raw.githubusercontent.com/bibinthomas123/Tessera/main/docs/images/ui-schedule.png)

## Quick start

```bash
git clone https://github.com/bibinthomas123/Tessera.git && cd Tessera
python -m venv .venv && .venv/bin/pip install -e ".[dev]"

# 1. build a versioned catalog snapshot (here: the GLADE+ excerpt shipped for tests)
.venv/bin/tessera ingest-catalog glade+ tests/data/glade_plus_gw170817.tsv \
    --out catalogs/glade-excerpt --upstream-version 2022-10-31

# 2. plan GW170817 for an example 1 m telescope
.venv/bin/tessera plan --skymap tests/data/GW170817_skymap.fits.gz --catalog catalogs/glade-excerpt \
    --profile example-1m --preset bns --nights 1 --t-start 2017-08-17T22:00:00 --out out/gw170817

# 3. re-check feasibility independently, and export for an ACP-controlled telescope
.venv/bin/tessera validate-schedule "out/gw170817/schedules/example-1m@1.json"
.venv/bin/tessera export "out/gw170817/schedules/example-1m@1.json" --adapter acp
```

As a library, inside an existing pipeline:

```python
from tessera import plan, PlanConfig

r = plan("bayestar.multiorder.fits", "catalogs/gladeplus-2022.10", ["example-1m"],
         PlanConfig(preset="bns"), t_start="2017-08-17T22:00:00")
r.galaxies[:10]            # ranked galaxies (astropy Table)
r.budget["p_cat"]          # share of probability the catalog can explain
r.schedules["example-1m@1"]["visits"]
r.write("out/")            # all outputs + manifest.json
```

## Running the service

```bash
cp .env.example .env      # then set POSTGRES_PASSWORD (e.g. openssl rand -hex 24)
docker compose up -d      # api, worker, alert listener, PostgreSQL
docker compose exec api tessera keys create --name ops --scopes read,plan,approve
```

The UI is at `http://localhost:8000/ui/` and the API docs at `/v1/docs`. For a multi-user deployment
behind an institutional gateway, there is a Helm chart in `deploy/helm/tessera`, and the service
supports OIDC and gateway authentication; see [docs/deployment.md](docs/deployment.md). Live LVK alerts
arrive from GCN Kafka, SCiMMA or GraceDB. New sky maps are treated as revisions of the same event, and
retractions withdraw schedules automatically.

## Validation

`tessera benchmark` regenerates [docs/benchmark-report.md](docs/benchmark-report.md):

| Check | Result |
| --- | --- |
| NGC 4993 rank for GW170817 (3D scoring with luminosity weight / sky position only) | 1st / 14th |
| Top-10 host recovery, 200 synthetic injections (3D / 2D baseline) | 0.43 / 0.07 |
| Catalog-completeness prediction vs truth (randomly halved catalog) | within 0.01 |
| Greedy scheduler vs provable optimum (MILP upper bound) | ≥ 0.91 |
| `p_sky` and credible levels vs `ligo.skymap` crossmatch | agree to 1e-12 |
| Scheduled visits passing the independent validator | 100% |
| Deterministic outputs; malformed and oversized inputs rejected with typed errors | yes |

There are 105 tests, covering service integration, alert revision handling and a GW170817 golden
regression. CI runs them on Python 3.10, 3.11 and 3.12 and builds the Docker image.

## Installing

The PyPI distribution name is **`tessera-gw`**, because the name `tessera` on PyPI belongs to an
unrelated project. The Python package and the command are both `tessera`. Optional extras:
`[service]`, `[postgres]`, `[alerts]`, `[rubin]`, `[validation]`.

## Data and caveats

- **Test data is public.** It includes the GWTC-1 GW170817 sky map (LIGO DCC P1800381), the GraceDB
  BAYESTAR maps for S190425z and S200105ae, and a GLADE+ excerpt from VizieR VII/291 (Dalya et al.
  2022; please cite it). The AT2017gfo light curve is derived from Villar et al. 2017, Table 3.
- **Telescope profiles** shipped here are illustrative or approximate, not facility-validated.
- **Science model 1.2.0.** The light curve, all reference densities (B-band, K-band, and stellar
  mass on the same GAMA/Chabrier scale as GLADE+ masses) and the kilonova diversity σ_M = 1.0 mag
  are sourced from the literature; see [docs/science.md](docs/science.md).
- **Adapters.**
  - The Rubin adapter is checked against `rubin-scheduler`.
  - LCO is checked against the open-source Observation Portal and a dated snapshot of LCO's
    public instrument list (instrument types, telescope classes, filters, modes).
  - Treasure Map is checked against its official client.
  - ACP is unverified.
  - Each facility must confirm the format before submission is enabled. See
    [docs/adapters.md](docs/adapters.md).

Further reading: [science model](docs/science.md) · [API](docs/api.md) ·
[validation protocol](docs/validation-protocol.md) · [contributing](CONTRIBUTING.md) ·
[roadmap](ROADMAP.md). Licensed under BSD-3-Clause.

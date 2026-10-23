# Deployment

## Single group (Docker Compose)

```bash
mkdir catalogs
tessera ingest-catalog glade+ GLADE+.txt --out catalogs/gladeplus-2022.10 --upstream-version 2022-10-31 \
     --source-url http://elysium.elte.hu/~dalyag/GLADE+.txt --ebv-map SFD_dust_4096.fits
export TESSERA_DEFAULT_CATALOG=gladeplus-2022.10
docker compose up -d
docker compose exec api tessera keys create --name ops --scopes read,plan,approve
```

UI: `http://localhost:8000/ui/` · OpenAPI: `http://localhost:8000/v1/docs` · metrics: `/metrics`
(scope `metrics`; worker metrics on port 9101).

The listener in the compose file watches `/data/alerts` (drop LVK JSON alerts there). For live alerts
replace its command with `tessera listen --source gcn` and set `GCN_CLIENT_ID` / `GCN_CLIENT_SECRET`
(needs the `alerts` extra), or `--source scimma` with hop credentials.

## Behind an institutional gateway

Terminate TLS and authenticate at the gateway (e.g. oauth2-proxy). Forward identity as
`X-Tessera-User` and scopes as `X-Tessera-Scopes` (comma-separated), and list the gateway's address in
`TESSERA_TRUSTED_PROXIES` (CIDRs); these headers are ignored from any other source. Set
`TESSERA_ROOT_PATH` if served under a sub-path. Scale workers horizontally; they share Postgres.

## Kubernetes (Helm)

```bash
kubectl create secret generic tessera-secrets \
  --from-literal=TESSERA_DATABASE_URL="$TESSERA_DATABASE_URL"
helm install tessera deploy/helm/tessera \
  --set ingress.enabled=true --set ingress.host=tessera.example.org \
  --set config.defaultCatalog=gladeplus-2022.10 \
  --set config.oidc.issuer=https://idp.example.org --set config.oidc.audience=tessera
```

The chart runs api, worker and listener Deployments as non-root with a read-only root filesystem,
mounts catalogs read-only, and (by default) applies NetworkPolicies so workers accept no inbound
traffic except metrics. The data volume is shared by api and workers, so it needs `ReadWriteMany`
storage. PostgreSQL is external. The chart passes `helm lint --strict` and its rendered manifests
validate against the Kubernetes 1.29 schemas (kubeconform, strict).

## Authentication

* **API keys** (`tessera keys create`): `Authorization: Bearer tessera_<prefix>_<secret>`.
* **OIDC, in process**: set `TESSERA_OIDC_ISSUER` and `TESSERA_OIDC_AUDIENCE`. Keys come from
  `TESSERA_OIDC_JWKS_URL`, the issuer's discovery document, or `TESSERA_OIDC_JWKS_FILE` (air-gapped).
  Signature, issuer, audience and expiry are checked. Scopes are read from `TESSERA_OIDC_SCOPES_CLAIM`
  (default `scope`); `plan`, `tessera:plan` and `tessera.plan` are all accepted.
* **Gateway headers**: as below, only from `TESSERA_TRUSTED_PROXIES`.

## Configuration

| Variable | Default | Meaning |
| --- | --- | --- |
| `TESSERA_DATABASE_URL` | `sqlite:///./tessera-data/tessera.db` | SQLAlchemy URL (Postgres in production) |
| `TESSERA_DATA_DIR` | `./tessera-data` | artifacts, uploaded profiles |
| `TESSERA_CATALOG_DIRS` | `./catalogs` | where snapshots are found |
| `TESSERA_DEFAULT_CATALOG` | — | snapshot used when a request names none |
| `TESSERA_TRIGGER_POLICY` | — | YAML rules mapping alerts to plans (below) |
| `TESSERA_MODE` | `production` | `test` also plans MS/TS mock events |
| `TESSERA_MAX_UPLOAD_BYTES` | 100 MiB | upload limit (decompressed limit 1 GiB) |
| `TESSERA_MAX_QUEUE_DEPTH` / `TESSERA_MAX_JOBS_PER_KEY` | 50 / 4 | admission control (HTTP 429) |
| `TESSERA_JOB_TIMEOUT_S` / `TESSERA_JOB_MEMORY_BYTES` | 300 / 4 GiB | per-job limits |
| `TESSERA_SKYMAP_URL_ALLOWLIST` | GraceDB hosts | hosts allowed for `skymap_url` |
| `TESSERA_IERS_A_PATH` | — | local IERS-A table (auto-download is off) |
| `TESSERA_AUTH_DISABLED` | 0 | single-user local use only |

Trigger policy example:

```yaml
require_significant: true
alert_types: [PRELIMINARY, INITIAL, UPDATE]
rules:
  - preset: bns
    profiles: [example-1m]
    when: {min_has_ns: 0.1}
  - preset: bbh
    profiles: [example-wide-north]
    when: {max_has_ns: 0.1, max_p_terrestrial: 0.5}
```

## Safety defaults

Schedules are `proposed` until a person with the `approve` scope approves one exact version.
Submission additionally requires `submission_enabled: true` in the profile (admin-published) and
adapter credentials (`TESSERA_LCO_TOKEN`, `TESSERA_TREASUREMAP_TOKEN`). Rubin ToO packages are never
submitted automatically.

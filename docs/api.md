# REST API v1

Authenticate with `Authorization: Bearer tessera_<prefix>_<secret>` (or `X-API-Key`). Errors are
RFC 9457 `application/problem+json` with a stable `code`. Full schema: `/v1/openapi.json`.

| Method and path | Scope | Purpose |
| --- | --- | --- |
| `POST /v1/skymaps` (multipart `file`, or form `url`) | plan | store a sky map; returns `skymap_id` |
| `POST /v1/plans` (+ `Idempotency-Key`) | plan | `{skymap_id or skymap_url, catalog, profiles, config, t_start, event, prior_observations}` → job |
| `GET/DELETE /v1/jobs/{id}` | read/plan | status (stage, attempts, error) / cancel |
| `GET /v1/plans/{key}` | read | summary, budget, strategy comparison, diagnostics, schedules |
| `GET /v1/plans/{key}/galaxies?format=csv\|ecsv\|json\|parquet\|votable` | read | ranked galaxies |
| `GET /v1/plans/{key}/fields/{profile}?format=csv\|json` | read | candidate fields and reasons |
| `GET /v1/plans/{key}/manifest`, `/files/{path}` | read | manifest, any output file |
| `POST /v1/plans/{key}/rerun` | plan | rerun and diff output hashes |
| `GET /v1/events`, `/v1/events/{sid}` | read | revisions with deltas between them |
| `POST /v1/alerts` | admin | inject an LVK JSON alert (replay, testing) |
| `GET /v1/schedules/{id}[/versions/{v}]` | read | schedule document and history |
| `POST /v1/schedules/{id}/versions` | plan | `{base_version, visits}` → validated new version |
| `POST /v1/schedules/{id}/versions/{v}/approve` | approve | bind approval to that version's sha256 |
| `GET /v1/schedules/{id}/versions/{v}/export?adapter=` | read | generic, acp, lco, treasuremap, rubin_too |
| `POST /v1/schedules/{id}/versions/{v}/submit?adapter=` | approve | approved + submission-enabled only |
| `POST /v1/observations` | plan | completed/failed/weathered report; `replan: true` credits it |
| `POST /v1/transients/associate` | read | galaxies near a transient + covering fields |
| `GET/POST /v1/profiles`, `GET /v1/catalogs`, `GET /v1/science` | read/admin | reference data |
| `GET /healthz`, `/readyz`, `/metrics`, `/v1/version` | public/metrics | operations |

Example:

```bash
KEY=tessera_...; H="Authorization: Bearer $KEY"
SID=$(curl -s -H "$H" -F file=@bayestar.multiorder.fits localhost:8000/v1/skymaps | jq -r .skymap_id)
curl -s -H "$H" -H 'content-type: application/json' localhost:8000/v1/plans \
  -d "{\"skymap_id\":\"$SID\",\"profiles\":[\"example-1m\"],\"config\":{\"preset\":\"bns\"}}"
```

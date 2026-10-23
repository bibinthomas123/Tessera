"""REST API, versioned under ``/v1`` (spec 10.1).  OpenAPI at ``/v1/openapi.json``."""

from __future__ import annotations

import json
import os
import uuid
from pathlib import Path
from typing import Any

import structlog
from fastapi import Depends, FastAPI, File, Form, Header, Query, Request, Response, UploadFile
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from prometheus_client import CONTENT_TYPE_LATEST, generate_latest
from pydantic import BaseModel, Field
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from tessera import SCIENCE_MODEL_VERSION, __version__
from tessera.adapters.alerts.envelope import parse_lvk_json
from tessera.adapters.schedulers import available as available_adapters
from tessera.core import science
from tessera.core.errors import TesseraError
from tessera.core.profile import list_profiles, load_profile, parse_profile
from tessera.service import obs
from tessera.service.auth import AuthError, Principal, verify_key
from tessera.service.db import (
    AuditLog,
    Event,
    Job,
    Observation,
    Plan,
    Revision,
    Schedule,
    ScheduleVersion,
    audit,
    new_id,
    session_factory,
)
from tessera.service.events import TriggerPolicy, handle_alert
from tessera.service.jobs import NotFoundError, enqueue, list_catalogs, request_cancel, store_skymap
from tessera.service.settings import Settings, get_settings
from tessera.service.store import ArtifactStore

STATIC = Path(__file__).parent / "static"


# ---------------------------------------------------------------------------
# request models
# ---------------------------------------------------------------------------

class PlanRequest(BaseModel):
    skymap_id: str | None = None
    skymap_url: str | None = None
    catalog: str | None = None
    profiles: list[str] = Field(default_factory=list)
    config: dict[str, Any] = Field(default_factory=dict)
    t_start: str | None = None
    event: dict[str, Any] | None = None
    prior_observations: list[dict[str, Any]] | None = None


class EditRequest(BaseModel):
    base_version: int
    visits: list[dict[str, Any]]


class ObservationReport(BaseModel):
    superevent_id: str | None = None
    schedule_id: str | None = None
    profile_id: str | None = None
    field_id: str | None = None
    ra: float
    dec: float
    pa_deg: float | None = None
    filter: str | None = None
    start: str
    end: str | None = None
    status: str = Field(pattern="^(completed|failed|weathered|aborted)$")
    m_lim: float | None = None
    replan: bool = False


class AssociateRequest(BaseModel):
    plan_key: str
    ra: float
    dec: float
    max_offset_kpc: float = 50.0


# ---------------------------------------------------------------------------
# app factory
# ---------------------------------------------------------------------------

def create_app(settings: Settings | None = None) -> FastAPI:
    settings = settings or get_settings()
    obs.configure_logging(os.environ.get("TESSERA_LOG_LEVEL", "INFO"))
    prof_dir = settings.data_dir / "profiles"
    prof_dir.mkdir(parents=True, exist_ok=True)
    dirs = [d for d in os.environ.get("TESSERA_PROFILE_DIRS", "").split(os.pathsep) if d]
    if str(prof_dir) not in dirs:
        os.environ["TESSERA_PROFILE_DIRS"] = os.pathsep.join([str(prof_dir), *dirs])
    S = session_factory(settings.database_url)
    store = ArtifactStore(settings.artifacts_dir)
    policy = TriggerPolicy.load(settings.trigger_policy, settings.mode)
    app = FastAPI(title="tessera", version=__version__, root_path=settings.root_path,
                  openapi_url="/v1/openapi.json", docs_url="/v1/docs", redoc_url=None,
                  description="Telescope-agnostic gravitational-wave follow-up engine.")
    app.state.settings, app.state.sessions, app.state.store, app.state.policy = settings, S, store, policy
    from tessera.service.oidc import OIDCConfig, OIDCVerifier

    oidc_cfg = OIDCConfig.from_env()
    oidc = OIDCVerifier(oidc_cfg) if oidc_cfg else None
    if settings.auth_disabled:
        obs.log.warning("auth_disabled", message="TESSERA_AUTH_DISABLED=1: every request has admin rights")

    def db():
        with S() as s:
            yield s

    def principal(request: Request, s: Session = Depends(db),
                  authorization: str | None = Header(None), x_api_key: str | None = Header(None)) -> Principal:
        if settings.auth_disabled:
            return Principal("local", frozenset({"admin"}), "disabled")
        client = request.client.host if request.client else None
        if settings.is_trusted_proxy(client) and request.headers.get("x-tessera-user"):
            scopes = frozenset(x.strip() for x in request.headers.get("x-tessera-scopes", "read").split(",") if x.strip())
            return Principal(request.headers["x-tessera-user"], scopes, "gateway")
        token = x_api_key or (authorization[7:] if authorization and authorization.lower().startswith("bearer ") else None)
        if not token:
            raise AuthError("missing credentials (Authorization: Bearer <API key or OIDC token>, or X-API-Key)")
        if token.startswith("tessera_"):
            return verify_key(s, token)
        if oidc is not None:
            return oidc.verify(token)
        raise AuthError("bearer token is not a tessera API key and OIDC is not configured")

    def need(scope: str):
        def dep(p: Principal = Depends(principal)) -> Principal:
            p.require(scope)
            return p
        return dep

    # -- middleware & errors -----------------------------------------------
    @app.middleware("http")
    async def request_context(request: Request, call_next):
        rid = request.headers.get("x-request-id") or uuid.uuid4().hex
        structlog.contextvars.clear_contextvars()
        structlog.contextvars.bind_contextvars(request_id=rid, path=request.url.path, method=request.method)
        cl = request.headers.get("content-length")
        if cl and cl.isdigit() and int(cl) > settings.max_upload_bytes + 1024 * 1024:
            return JSONResponse(TesseraError("request body too large", code="SKYMAP_TOO_LARGE").to_problem() |
                                {"status": 413}, status_code=413, media_type="application/problem+json")
        response = await call_next(request)
        response.headers["x-request-id"] = rid
        obs.log.info("request", status=response.status_code)
        return response

    @app.exception_handler(TesseraError)
    async def tessera_error(_: Request, exc: TesseraError):
        return JSONResponse(exc.to_problem(), status_code=exc.http_status, media_type="application/problem+json")

    @app.exception_handler(RequestValidationError)
    async def validation_error(_: Request, exc: RequestValidationError):
        body = TesseraError("request validation failed", code="REQUEST_INVALID",
                         details={"errors": json.loads(json.dumps(exc.errors(), default=str))}).to_problem()
        body["status"] = 422
        return JSONResponse(body, status_code=422, media_type="application/problem+json")

    # -- operations ---------------------------------------------------------
    @app.get("/healthz", include_in_schema=False)
    def healthz():
        return {"status": "ok"}

    @app.get("/readyz", include_in_schema=False)
    def readyz(s: Session = Depends(db)):
        s.execute(select(func.count()).select_from(Job))
        return {"status": "ready", "catalogs": len(list_catalogs(settings))}

    @app.get("/v1/version")
    def version():
        return {"tessera": __version__, "science_model": SCIENCE_MODEL_VERSION, "api": "v1",
                "adapters": available_adapters()}

    @app.get("/metrics", include_in_schema=False)
    def metrics(s: Session = Depends(db), _: Principal = Depends(need("metrics"))):
        depth = s.scalar(select(func.count()).select_from(Job).where(Job.state.in_(("queued", "running")))) or 0
        obs.QUEUE_DEPTH.set(depth)
        return Response(generate_latest(), media_type=CONTENT_TYPE_LATEST)

    # -- reference data ------------------------------------------------------
    @app.get("/v1/profiles")
    def profiles(_: Principal = Depends(need("read"))):
        out = []
        for name in list_profiles():
            try:
                p = load_profile(name)
                out.append({"id": p.id, "version": p.version, "key": p.key, "description": p.description,
                            "validation": p.validation.model_dump(), "adapters": p.adapters,
                            "footprint": p.instrument.footprint.kind, "filters": sorted(p.instrument.filters)})
            except TesseraError as exc:
                out.append({"id": name, "error": exc.code})
        return {"profiles": out}

    @app.get("/v1/profiles/{pid}")
    def profile(pid: str, _: Principal = Depends(need("read"))):
        p = load_profile(pid)
        return {**p.canonical(), "sha256": p.sha256}

    @app.post("/v1/profiles", status_code=201)
    def add_profile(doc: dict[str, Any], s: Session = Depends(db), who: Principal = Depends(need("admin"))):
        p = parse_profile(doc)
        target = prof_dir / f"{p.id}.yaml"
        try:
            cur = load_profile(p.id)
            if p.version <= cur.version:
                raise TesseraError(f"profile {p.id} version must exceed {cur.version} (profiles are immutable)",
                                code="PROFILE_VERSION_CONFLICT")
        except TesseraError as exc:
            if exc.code != "PROFILE_NOT_FOUND":
                raise
        import yaml

        target.write_text(yaml.safe_dump(p.canonical(), sort_keys=False))
        (prof_dir / "history").mkdir(exist_ok=True)
        (prof_dir / "history" / f"{p.id}@{p.version}.yaml").write_text(yaml.safe_dump(p.canonical(), sort_keys=False))
        audit(s, who.name, "profile_publish", p.key, sha256=p.sha256)
        s.commit()
        return {"key": p.key, "sha256": p.sha256}

    @app.get("/v1/catalogs")
    def catalogs(_: Principal = Depends(need("read"))):
        return {"catalogs": list_catalogs(settings)}

    @app.get("/v1/science")
    def science_defaults(_: Principal = Depends(need("read"))):
        d = science.defaults()
        return {k: d[k] for k in ("science_model", "presets", "strategies", "light_curves", "reference_density",
                                  "cosmology", "peculiar_velocity_kms")}

    # -- sky maps and plans ----------------------------------------------------
    @app.post("/v1/skymaps", status_code=201)
    async def upload_skymap(file: UploadFile | None = File(None), url: str | None = Form(None),
                            s: Session = Depends(db), who: Principal = Depends(need("plan"))):
        if file is None and not url:
            raise TesseraError("provide a file or a url", code="REQUEST_INVALID")
        if file is not None:
            buf = bytearray()
            while chunk := await file.read(1 << 20):
                buf += chunk
                if len(buf) > settings.max_upload_bytes:
                    raise TesseraError("upload exceeds size limit", code="SKYMAP_TOO_LARGE")
            raw, name = bytes(buf), file.filename or "upload.fits"
        else:
            from tessera.adapters.alerts.sources import fetch_url

            raw, name = fetch_url(url, max_bytes=settings.max_upload_bytes), url.rsplit("/", 1)[-1]
        row = store_skymap(s, store, raw, name, who.name)
        return {"skymap_id": row.id, "sha256_data": row.sha256_data, "filename": row.filename,
                "note": "validated by the worker when planned"}

    @app.post("/v1/plans", status_code=202)
    def create_plan(body: PlanRequest, response: Response, s: Session = Depends(db),
                    who: Principal = Depends(need("plan")), idempotency_key: str | None = Header(None)):
        req = body.model_dump()
        if not req.get("skymap_id"):
            if not req.get("skymap_url"):
                raise TesseraError("skymap_id or skymap_url is required", code="REQUEST_INVALID")
            from tessera.adapters.alerts.sources import fetch_url

            raw = fetch_url(req["skymap_url"], max_bytes=settings.max_upload_bytes)
            req["skymap_id"] = store_skymap(s, store, raw, req["skymap_url"].rsplit("/", 1)[-1], who.name).id
        req.pop("skymap_url", None)
        job, plan, created = enqueue(s, settings, req, owner=who.name, idempotency_key=idempotency_key)
        if plan is not None:
            response.status_code = 200
            return {"state": "succeeded", "plan_key": plan.plan_key, "job_id": plan.job_id, "existing": True}
        response.status_code = 202 if created else 200
        if settings.inline_worker and created:
            from tessera.service.worker import run_one

            run_one(settings, "inline")
            s.expire_all()
            job = s.get(Job, job.id)
        return _job_view(job) | {"existing": not created}

    def _job_view(j: Job) -> dict[str, Any]:
        return {"job_id": j.id, "plan_key": j.plan_key, "state": j.state, "stage": j.stage, "attempts": j.attempts,
                "error": j.error, "result": j.result, "created": j.created.isoformat() if j.created else None,
                "duration_s": j.duration_s, "revision_id": j.revision_id}

    @app.get("/v1/jobs/{job_id}")
    def get_job(job_id: str, s: Session = Depends(db), _: Principal = Depends(need("read"))):
        j = s.get(Job, job_id)
        if j is None:
            raise NotFoundError(f"job {job_id} not found")
        return _job_view(j)

    @app.delete("/v1/jobs/{job_id}")
    def cancel_job(job_id: str, s: Session = Depends(db), who: Principal = Depends(need("plan"))):
        return _job_view(request_cancel(s, job_id, who.name))

    @app.get("/v1/jobs")
    def list_jobs(state: str | None = None, limit: int = Query(50, le=500), s: Session = Depends(db),
                  _: Principal = Depends(need("read"))):
        q = select(Job).order_by(Job.created.desc()).limit(limit)
        if state:
            q = q.where(Job.state == state)
        return {"jobs": [_job_view(j) for j in s.scalars(q).all()]}

    def _plan(s: Session, key: str) -> Plan:
        key = key if key.startswith("sha256:") else f"sha256:{key}"
        p = s.get(Plan, key)
        if p is None:
            raise NotFoundError(f"plan {key} not found")
        return p

    def _file(p: Plan, rel: str) -> bytes:
        path = Path(p.artifact_dir) / rel
        if not path.resolve().is_relative_to(Path(p.artifact_dir).resolve()) or not path.is_file():
            raise NotFoundError(f"{rel} not found in plan")
        return path.read_bytes()

    @app.get("/v1/plans")
    def list_plans(limit: int = Query(50, le=500), s: Session = Depends(db), _: Principal = Depends(need("read"))):
        rows = s.scalars(select(Plan).order_by(Plan.created.desc()).limit(limit)).all()
        return {"plans": [{"plan_key": p.plan_key, "created": p.created.isoformat(), "event": p.summary.get("event"),
                           "budget": p.summary.get("budget")} for p in rows]}

    @app.get("/v1/plans/{key}")
    def get_plan(key: str, s: Session = Depends(db), _: Principal = Depends(need("read"))):
        p = _plan(s, key)
        schedules = s.scalars(select(Schedule).where(Schedule.plan_key == p.plan_key)).all()
        return {"plan_key": p.plan_key, "summary": p.summary,
                "budget": json.loads(_file(p, "budget.json")),
                "strategy_comparison": json.loads(_file(p, "strategy_comparison.json")),
                "diagnostics": json.loads(_file(p, "diagnostics.json"))["diagnostics"],
                "schedules": [{"schedule_id": x.id, "profile": x.profile_key, "status": x.status,
                               "current_version": x.current_version} for x in schedules],
                "files": sorted(p.outputs)}

    @app.get("/v1/plans/{key}/galaxies")
    def plan_galaxies(key: str, format: str = Query("csv", pattern="^(csv|ecsv|json|parquet|votable)$"),
                      s: Session = Depends(db), _: Principal = Depends(need("read"))):
        p = _plan(s, key)
        if format in ("csv", "ecsv", "json"):
            data = _file(p, f"galaxies.{format}")
        else:
            from astropy.table import Table

            from tessera.core.io import table_to_bytes

            data = table_to_bytes(Table.read(Path(p.artifact_dir) / "galaxies.ecsv", format="ascii.ecsv"), format)
        from tessera.core.io import MEDIA_TYPES

        return Response(data, media_type=MEDIA_TYPES[format],
                        headers={"content-disposition": f'attachment; filename="galaxies.{format}"'})

    @app.get("/v1/plans/{key}/fields/{profile}")
    def plan_fields(key: str, profile: str, format: str = Query("csv", pattern="^(csv|json)$"),
                    s: Session = Depends(db), _: Principal = Depends(need("read"))):
        p = _plan(s, key)
        from tessera.core.io import MEDIA_TYPES

        return Response(_file(p, f"fields/{profile}.{format}"), media_type=MEDIA_TYPES[format],
                        headers={"content-disposition": f'attachment; filename="fields-{profile}.{format}"'})

    @app.get("/v1/plans/{key}/files/{name:path}")
    def plan_file(key: str, name: str, s: Session = Depends(db), _: Principal = Depends(need("read"))):
        p = _plan(s, key)
        mt = "application/json" if name.endswith(".json") else "text/plain"
        return Response(_file(p, name), media_type=mt)

    @app.get("/v1/plans/{key}/manifest")
    def plan_manifest(key: str, s: Session = Depends(db), _: Principal = Depends(need("read"))):
        return json.loads(_file(_plan(s, key), "manifest.json"))

    @app.post("/v1/plans/{key}/rerun", status_code=202)
    def rerun(key: str, s: Session = Depends(db), who: Principal = Depends(need("plan"))):
        p = _plan(s, key)
        orig = s.get(Job, p.job_id)
        req = {**orig.request, "rerun_of": p.plan_key}
        j = Job(id=new_id("job"), plan_key=p.plan_key, owner=who.name, state="queued", request=req,
                idempotency_key=f"rerun:{p.plan_key}:{uuid.uuid4().hex}")
        s.add(j)
        audit(s, who.name, "rerun", p.plan_key)
        s.commit()
        if settings.inline_worker:
            from tessera.service.worker import run_one

            run_one(settings, "inline")
            s.expire_all()
            j = s.get(Job, j.id)
        return _job_view(j)

    # -- events, alerts ---------------------------------------------------------
    @app.get("/v1/events")
    def list_events(s: Session = Depends(db), _: Principal = Depends(need("read"))):
        rows = s.scalars(select(Event).order_by(Event.updated.desc()).limit(200)).all()
        return {"events": [{"superevent_id": e.superevent_id, "state": e.state,
                            "active_revision_id": e.active_revision_id} for e in rows]}

    @app.get("/v1/events/{sid}")
    def get_event(sid: str, s: Session = Depends(db), _: Principal = Depends(need("read"))):
        e = s.get(Event, sid)
        if e is None:
            raise NotFoundError(f"event {sid} not found")
        revs = s.scalars(select(Revision).where(Revision.superevent_id == sid).order_by(Revision.number)).all()
        out = []
        prev_dir = None
        from tessera.core.delta import revision_delta

        for r in revs:
            plans = s.scalars(select(Plan).where(Plan.revision_id == r.id)).all()
            item = {"revision_id": r.id, "number": r.number, "state": r.state, "alert_types": r.alert_types,
                    "time_created": r.time_created, "pinned": r.pinned, "skymap_id": r.skymap_id,
                    "plans": [p.plan_key for p in plans]}
            if plans and prev_dir is not None:
                try:
                    item["delta_from_previous"] = revision_delta(prev_dir, plans[0].artifact_dir)
                except Exception as exc:  # pragma: no cover - missing artifacts
                    item["delta_from_previous"] = {"error": str(exc)}
            if plans:
                prev_dir = plans[0].artifact_dir
            out.append(item)
        return {"superevent_id": sid, "state": e.state, "active_revision_id": e.active_revision_id,
                "latest_alert": (e.meta or {}).get("latest_alert"), "revisions": out}

    @app.post("/v1/events/{sid}/revisions/{rid}/pin")
    def pin_revision(sid: str, rid: str, pinned: bool = True, s: Session = Depends(db),
                     who: Principal = Depends(need("admin"))):
        r = s.get(Revision, rid)
        if r is None or r.superevent_id != sid:
            raise NotFoundError("revision not found")
        r.pinned = pinned
        audit(s, who.name, "revision_pin", rid, pinned=pinned)
        s.commit()
        return {"revision_id": rid, "pinned": pinned}

    @app.post("/v1/alerts")
    async def ingest_alert(request: Request, s: Session = Depends(db), who: Principal = Depends(need("admin"))):
        raw = await request.body()
        env = parse_lvk_json(raw, source=f"api:{who.name}")
        out = handle_alert(s, settings, store, env, policy)
        obs.ALERTS.labels(source="api", type=env.alert_type, outcome=out.outcome).inc()
        if settings.inline_worker:
            from tessera.service.worker import run_one

            while run_one(settings, "inline"):
                pass
        return out.to_dict()

    # -- schedules ----------------------------------------------------------------
    from tessera.service import schedules as sched

    def _sched_view(s: Session, sid: str, version: int | None = None) -> dict[str, Any]:
        row = sched.get_schedule(s, sid)
        ver = sched.get_version(s, sid, version)
        appr = sched.latest_approval(s, sid, ver.version)
        return {"schedule_id": sid, "status": row.status, "current_version": row.current_version,
                "version": ver.version, "sha256": ver.sha256, "diff_from_parent": ver.diff,
                "approval": None if appr is None else {"approver": appr.approver, "created": appr.created.isoformat(),
                                                       "auto": appr.auto, "sha256": appr.sha256},
                "versions": [v.version for v in s.scalars(select(ScheduleVersion).where(
                    ScheduleVersion.schedule_id == sid).order_by(ScheduleVersion.version)).all()],
                "doc": ver.doc}

    @app.get("/v1/schedules/{sid}")
    def get_schedule(sid: str, s: Session = Depends(db), _: Principal = Depends(need("read"))):
        return _sched_view(s, sid)

    @app.get("/v1/schedules/{sid}/versions/{v}")
    def get_schedule_version(sid: str, v: int, s: Session = Depends(db), _: Principal = Depends(need("read"))):
        return _sched_view(s, sid, v)

    @app.post("/v1/schedules/{sid}/versions", status_code=201)
    def edit_schedule(sid: str, body: EditRequest, s: Session = Depends(db), who: Principal = Depends(need("plan"))):
        row = sched.edit(s, sid, body.base_version, body.visits, who.name)
        return _sched_view(s, sid, row.version)

    @app.post("/v1/schedules/{sid}/versions/{v}/approve")
    def approve_schedule(sid: str, v: int, s: Session = Depends(db), who: Principal = Depends(need("approve"))):
        sched.approve(s, sid, v, who.name)
        return _sched_view(s, sid, v)

    @app.get("/v1/schedules/{sid}/versions/{v}/export")
    def export_schedule(sid: str, v: int, adapter: str = "generic", fmt: str | None = None,
                        allow_infeasible: bool = False, s: Session = Depends(db), _: Principal = Depends(need("read"))):
        opts = {"fmt": fmt} if fmt and adapter == "generic" else {}
        exp, issues = sched.export(s, sid, v, adapter, store, allow_infeasible=allow_infeasible, **opts)
        data, mt, name = exp.as_bytes()
        headers = {"content-disposition": f'attachment; filename="{name}"', "x-tessera-adapter": f"{exp.adapter}/{exp.adapter_version}"}
        if issues:
            headers["x-tessera-issues"] = "; ".join(f"{i.level}:{i.message}" for i in issues)[:1000]
        return Response(data, media_type=mt, headers=headers)

    @app.post("/v1/schedules/{sid}/versions/{v}/submit")
    def submit_schedule(sid: str, v: int, adapter: str, s: Session = Depends(db), who: Principal = Depends(need("approve"))):
        sub = sched.submit(s, sid, v, adapter, store, who.name)
        return {"submission_id": sub.id, "state": sub.state, "receipt": sub.receipt}

    # -- feedback, candidates -------------------------------------------------------
    @app.post("/v1/observations", status_code=201)
    def report_observation(body: ObservationReport, s: Session = Depends(db), who: Principal = Depends(need("plan"))):
        o = Observation(id=new_id("obs"), superevent_id=body.superevent_id, schedule_id=body.schedule_id,
                        profile_id=body.profile_id, status=body.status, payload=body.model_dump(), reported_by=who.name)
        s.add(o)
        audit(s, who.name, "observation", o.id, status=body.status)
        s.commit()
        out: dict[str, Any] = {"observation_id": o.id}
        if body.replan and body.superevent_id:
            ev = s.get(Event, body.superevent_id)
            if ev is None or ev.active_revision_id is None:
                raise NotFoundError("no active revision to replan")
            rev = s.get(Revision, ev.active_revision_id)
            base = s.scalars(select(Plan).where(Plan.revision_id == rev.id).order_by(Plan.created.desc())).first()
            if base is None:
                raise NotFoundError("no plan exists for the active revision")
            prior = [x.payload for x in s.scalars(select(Observation).where(
                Observation.superevent_id == body.superevent_id)).all()]
            from astropy.time import Time

            from tessera.core.timeutil import iso

            req = {**s.get(Job, base.job_id).request, "prior_observations": prior, "t_start": iso(Time.now())}
            for k in ("plan_key", "catalog_path"):
                req.pop(k, None)
            job, plan, _ = enqueue(s, settings, req, owner=who.name, revision_id=rev.id)
            out["replan"] = {"job_id": job.id if job else None, "plan_key": (job or plan).plan_key}
        return out

    @app.post("/v1/transients/associate")
    def associate(body: AssociateRequest, s: Session = Depends(db), _: Principal = Depends(need("read"))):
        from tessera.core.candidates import associate_transient

        p = _plan(s, body.plan_key)
        gal = json.loads(_file(p, "galaxies.json"))["rows"]
        fields = {}
        for rel in p.outputs:
            if rel.startswith("fields/") and rel.endswith(".json"):
                prof_key = rel[len("fields/"):-5]
                rows = json.loads(_file(p, rel))["rows"]
                try:
                    from tessera.core.footprint import InstrumentFootprint

                    rad = InstrumentFootprint(load_profile(prof_key.split("@")[0]).instrument.footprint).bounding_radius_deg
                except TesseraError:
                    rad = 1.0
                fields[prof_key] = [{**r, "_radius_deg": rad} for r in rows]
        return associate_transient(body.ra, body.dec, gal, fields, max_offset_kpc=body.max_offset_kpc)

    @app.get("/v1/audit")
    def audit_log(limit: int = Query(100, le=1000), s: Session = Depends(db), _: Principal = Depends(need("admin"))):
        rows = s.scalars(select(AuditLog).order_by(AuditLog.id.desc()).limit(limit)).all()
        return {"entries": [{"ts": r.ts.isoformat(), "actor": r.actor, "action": r.action, "target": r.target,
                             "details": r.details} for r in rows]}

    # -- UI -------------------------------------------------------------------------
    app.mount("/ui", StaticFiles(directory=STATIC, html=True), name="ui")

    @app.get("/", include_in_schema=False)
    def root():
        return RedirectResponse(url=f"{settings.root_path}/ui/")

    return app


def app_from_env() -> FastAPI:  # uvicorn --factory tessera.service.api:app_from_env
    return create_app()

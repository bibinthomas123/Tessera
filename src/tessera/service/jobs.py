"""Job queue on the database (spec 10.2) and job execution.

* The plan key is computed *before* enqueueing (same function as the core),
  so a repeated request returns the existing plan or job (idempotency).
* Claims use compare-and-set (``UPDATE ... WHERE state='queued'``), which is
  safe on both PostgreSQL and SQLite; expired leases are reclaimed.
* Execution happens in a separate process with wall-time and memory limits
  (see :mod:`tessera.service.worker`).
"""

from __future__ import annotations

import zlib
from datetime import timedelta
from pathlib import Path
from typing import Any

from sqlalchemy import func, select, update
from sqlalchemy.orm import Session

from tessera.core import timeutil
from tessera.core.catalog import open_snapshot
from tessera.core.config import PlanConfig
from tessera.core.errors import CatalogError, ConfigError, TesseraError
from tessera.core.hashing import sha256_bytes
from tessera.core.manifest import plan_identity
from tessera.core.profile import load_profile
from tessera.service.db import Job, Plan, Schedule, ScheduleVersion, Skymap, audit, new_id, utcnow
from tessera.service.settings import Settings
from tessera.service.store import ArtifactStore


class QueueFullError(TesseraError):
    code = "QUEUE_FULL"
    http_status = 429


class NotFoundError(TesseraError):
    code = "NOT_FOUND"
    http_status = 404


# ---------------------------------------------------------------------------
# Sky-map intake (the API never parses FITS; it only hashes, spec 10.3)
# ---------------------------------------------------------------------------

def data_sha256(raw: bytes, max_decompressed: int) -> str:
    """sha256 of the decompressed bytes (gzip streamed with a cap) or of the raw bytes."""
    if raw[:2] != b"\x1f\x8b":
        return sha256_bytes(raw)
    import hashlib

    h = hashlib.sha256()
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    total = 0
    data = raw
    try:
        while data:
            chunk = d.decompress(data, 1 << 22)
            total += len(chunk)
            if total > max_decompressed:
                raise TesseraError("decompressed sky map exceeds limit", code="SKYMAP_TOO_LARGE")
            h.update(chunk)
            data = d.unconsumed_tail
            if d.eof:
                break
    except zlib.error as exc:
        raise TesseraError(f"corrupt gzip stream: {exc}", code="SKYMAP_INVALID_FORMAT") from exc
    return h.hexdigest()


def store_skymap(s: Session, store: ArtifactStore, raw: bytes, filename: str, uploaded_by: str,
                 max_decompressed: int = 1024**3) -> Skymap:
    if not raw.startswith(b"SIMPLE  =") and raw[:2] != b"\x1f\x8b":
        raise TesseraError("upload is not a FITS file (or gzip-compressed FITS)", code="SKYMAP_INVALID_FORMAT")
    digest = store.put(raw)
    row = s.get(Skymap, digest)
    if row is None:
        row = Skymap(id=digest, filename=filename[:255], sha256_data=data_sha256(raw, max_decompressed),
                     uploaded_by=uploaded_by, summary={})
        s.add(row)
        s.commit()
    return row


# ---------------------------------------------------------------------------
# Input resolution (shared by API and worker so plan keys always agree)
# ---------------------------------------------------------------------------

def find_catalog(settings: Settings, name: str | None) -> Path:
    name = name or settings.default_catalog
    if not name:
        raise CatalogError("no catalog given and TESSERA_DEFAULT_CATALOG is not set", code="CATALOG_SNAPSHOT_MISSING")
    p = Path(name)
    if (p / "snapshot.json").is_file():
        return p
    for d in settings.catalog_dirs:
        if (d / name / "snapshot.json").is_file():
            return d / name
    raise CatalogError(f"catalog snapshot {name!r} not installed", code="CATALOG_SNAPSHOT_MISSING",
                       details={"searched": [str(d) for d in settings.catalog_dirs]})


def list_catalogs(settings: Settings) -> list[dict[str, Any]]:
    out = []
    for d in settings.catalog_dirs:
        if not d.is_dir():
            continue
        for sub in sorted(d.iterdir()):
            if (sub / "snapshot.json").is_file():
                try:
                    d = open_snapshot(sub).describe()
                    # "name" is what plan requests use to find the snapshot: its directory name.
                    out.append({**d, "snapshot_name": d["name"], "name": sub.name})
                except TesseraError:
                    continue
    return out


def bucket_time(iso: str, minutes: int = 5) -> str:
    t = timeutil.to_time(iso)
    step = minutes / 1440.0
    import math

    from astropy.time import Time

    return timeutil.iso(Time(math.floor(t.jd / step) * step, format="jd", scale="utc"))


def resolve_request(s: Session, settings: Settings, req: dict[str, Any]) -> dict[str, Any]:
    """Validate a plan request and compute its plan key without parsing FITS."""
    sk = s.get(Skymap, req["skymap_id"])
    if sk is None:
        raise NotFoundError(f"sky map {req['skymap_id']} not found")
    cat_path = find_catalog(settings, req.get("catalog"))
    snap = open_snapshot(cat_path)
    profiles = [load_profile(p) for p in req.get("profiles", [])]
    try:
        cfg = PlanConfig.model_validate(req.get("config") or {})
    except Exception as exc:
        raise ConfigError(f"invalid config: {exc}") from exc
    t_start = req.get("t_start")
    if not t_start:
        from astropy.time import Time

        t_start = bucket_time(timeutil.iso(Time.now()), 5)
    t_start = timeutil.iso(timeutil.to_time(t_start))
    resolved = cfg.resolved()
    key, _ = plan_identity(skymap_sha256=sk.sha256_data, snapshot_id=snap.id,
                           profiles=[{"key": p.key, "sha256": p.sha256} for p in profiles],
                           config=resolved.model_dump(mode="json"), t_start=t_start, event=req.get("event"),
                           prior_observations=req.get("prior_observations"))
    return {**req, "catalog_path": str(cat_path), "t_start": t_start, "plan_key": key,
            "profiles": [p.key.split("@")[0] for p in profiles]}


def enqueue(s: Session, settings: Settings, req: dict[str, Any], *, owner: str, idempotency_key: str | None = None,
            revision_id: str | None = None) -> tuple[Job | None, Plan | None, bool]:
    """Returns (job, existing_plan, created)."""
    if idempotency_key:
        j = s.scalar(select(Job).where(Job.idempotency_key == idempotency_key, Job.owner == owner))
        if j is not None:
            return j, s.get(Plan, j.plan_key), False
    resolved = resolve_request(s, settings, req)
    key = resolved["plan_key"]
    existing = s.get(Plan, key)
    if existing is not None:
        return None, existing, False
    j = s.scalar(select(Job).where(Job.plan_key == key, Job.state.in_(("queued", "running"))))
    if j is not None:
        return j, None, False
    depth = s.scalar(select(func.count()).select_from(Job).where(Job.state.in_(("queued", "running")))) or 0
    if depth >= settings.max_queue_depth:
        raise QueueFullError("job queue is full", details={"depth": depth})
    mine = s.scalar(select(func.count()).select_from(Job).where(Job.owner == owner,
                                                                 Job.state.in_(("queued", "running")))) or 0
    if owner != "system" and mine >= settings.max_jobs_per_key:
        raise QueueFullError("too many concurrent jobs for this key", details={"limit": settings.max_jobs_per_key})
    job = Job(id=new_id("job"), plan_key=key, idempotency_key=idempotency_key, owner=owner, state="queued",
              request=resolved, revision_id=revision_id)
    s.add(job)
    audit(s, owner, "enqueue", job.id, plan_key=key)
    s.commit()
    return job, None, True


def claim(s: Session, worker_id: str, lease_s: float) -> Job | None:
    now = utcnow()
    cands = s.scalars(select(Job.id).where(
        ((Job.state == "queued") & ((Job.not_before.is_(None)) | (Job.not_before <= now)))
        | ((Job.state == "running") & (Job.lease_expires < now))
    ).order_by(Job.created).limit(5)).all()
    for jid in cands:
        res = s.execute(update(Job).where(Job.id == jid, (Job.state == "queued") | (
            (Job.state == "running") & (Job.lease_expires < now))).values(
            state="running", lease_owner=worker_id, lease_expires=now + timedelta(seconds=lease_s),
            started=now, attempts=Job.attempts + 1))
        s.commit()
        if res.rowcount == 1:
            return s.get(Job, jid)
    return None


def heartbeat(s: Session, job_id: str, worker_id: str, lease_s: float) -> bool:
    res = s.execute(update(Job).where(Job.id == job_id, Job.lease_owner == worker_id, Job.state == "running")
                    .values(lease_expires=utcnow() + timedelta(seconds=lease_s)))
    s.commit()
    return res.rowcount == 1


def request_cancel(s: Session, job_id: str, actor: str, reason: str = "cancelled") -> Job:
    j = s.get(Job, job_id)
    if j is None:
        raise NotFoundError(f"job {job_id} not found")
    if j.state == "queued":
        j.state = "superseded" if reason == "superseded" else "cancelled"
        j.finished = utcnow()
    elif j.state == "running":
        j.cancel_requested = True
        j.error = {"code": "JOB_CANCELLED" if reason != "superseded" else "REVISION_SUPERSEDED"}
    audit(s, actor, "cancel", job_id, reason=reason)
    s.commit()
    return j


def finish(s: Session, job_id: str, state: str, *, error: dict[str, Any] | None = None, retry_in_s: float | None = None
           ) -> None:
    j = s.get(Job, job_id)
    if j is None:
        return
    now = utcnow()
    if retry_in_s is not None:
        j.state = "queued"
        j.not_before = now + timedelta(seconds=retry_in_s)
        j.lease_owner = None
        j.error = error
    else:
        j.state = state
        j.finished = now
        j.error = error
        if j.started:
            started = j.started if j.started.tzinfo else j.started.replace(tzinfo=now.tzinfo)
            j.duration_s = (now - started).total_seconds()
    s.commit()


# ---------------------------------------------------------------------------
# Execution (runs inside the worker's child process)
# ---------------------------------------------------------------------------

def execute(s: Session, settings: Settings, job_id: str) -> str:
    """Run the plan for ``job_id`` and persist outputs.  Returns the plan key."""
    import time

    from tessera.core.plan import plan
    from tessera.core.skymap import load_skymap

    job = s.get(Job, job_id)
    req = job.request
    store = ArtifactStore(settings.artifacts_dir)
    last = {"t": 0.0, "cancel": False}

    def cancel_check() -> bool:
        now = time.monotonic()
        if now - last["t"] > 1.0:
            s.expire(job)
            last["cancel"] = bool(s.get(Job, job_id).cancel_requested)
            last["t"] = now
        return last["cancel"]

    def progress(stage: str) -> None:
        s.execute(update(Job).where(Job.id == job_id).values(stage=stage))
        s.commit()

    sk = s.get(Skymap, req["skymap_id"])
    raw = store.get(sk.id)
    if req.get("rerun_of"):
        return _rerun(s, settings, job, req, raw, sk.filename, store)
    m = load_skymap(raw, filename=sk.filename)
    sk.summary = m.summary()
    s.commit()
    result = plan(m, req["catalog_path"], req["profiles"], req.get("config") or {}, req["t_start"],
                  event=req.get("event"), prior_observations=req.get("prior_observations"),
                  cancel_check=cancel_check, progress=progress)
    if result.plan_key != job.plan_key:
        raise TesseraError("plan key mismatch between API and worker (input drift)", code="PLAN_KEY_MISMATCH",
                        details={"api": job.plan_key, "worker": result.plan_key})
    out_dir = store.plan_dir(result.plan_key)
    result.write(out_dir)
    s.add(Plan(plan_key=result.plan_key, job_id=job_id, revision_id=job.revision_id, artifact_dir=str(out_dir),
               summary=result.summary, manifest=result.manifest, outputs=result.manifest["outputs"]))
    from tessera.adapters.schedulers.base import schedule_sha256

    for key, doc in result.schedules.items():
        sid = doc["schedule_id"]
        s.add(Schedule(id=sid, plan_key=result.plan_key, profile_key=key, revision_id=job.revision_id,
                       status=doc["status"], current_version=0))
        s.add(ScheduleVersion(schedule_id=sid, version=0, parent_version=None, doc=doc, sha256=schedule_sha256(doc),
                              created_by="tessera"))
    s.commit()
    from tessera.service.schedules import maybe_auto_approve

    for key, doc in result.schedules.items():
        maybe_auto_approve(s, doc["schedule_id"], key)
    return result.plan_key



def _rerun(s: Session, settings: Settings, job: Job, req: dict[str, Any], raw: bytes, filename: str,
           store: ArtifactStore) -> str:
    """Re-execute a stored plan and diff output hashes against its manifest (spec 7.5)."""
    from tessera.core.hashing import sha256_bytes as _h
    from tessera.core.plan import plan
    from tessera.core.skymap import load_skymap

    orig = s.get(Plan, req["rerun_of"])
    if orig is None:
        raise NotFoundError(f"plan {req['rerun_of']} not found")
    r = plan(load_skymap(raw, filename=filename), req["catalog_path"], req["profiles"], req.get("config") or {},
             req["t_start"], event=req.get("event"), prior_observations=req.get("prior_observations"))
    new = {k: _h(v) for k, v in r.outputs().items() if k != "manifest.json"}
    old = {k: v for k, v in orig.outputs.items()}
    diff = sorted(k for k in set(new) | set(old) if new.get(k) != old.get(k))
    job.result = {"rerun_of": orig.plan_key, "plan_key": r.plan_key, "identical": not diff, "differences": diff,
                  "versions_then": orig.manifest.get("versions"), "versions_now": r.manifest.get("versions")}
    s.commit()
    return orig.plan_key

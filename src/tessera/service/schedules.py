"""Schedule lifecycle with human-in-the-loop safety (spec 10.5).

proposed (v0) -> edited (v1..n) -> approved -> submitted -> executing -> completed
any state -> withdrawn | superseded

Every edit is a new immutable version, validated by the independent
validator.  v0 is never altered.  An approval binds approver, time and the
sha256 of one exact version.
"""

from __future__ import annotations

import json
from typing import Any

from sqlalchemy import select
from sqlalchemy.orm import Session

from tessera.adapters.schedulers import get_adapter
from tessera.adapters.schedulers.base import Approval as ApprovalDC
from tessera.adapters.schedulers.base import schedule_sha256
from tessera.core.errors import AdapterError, TesseraError
from tessera.core.profile import load_profile
from tessera.core.validate import validate_visits
from tessera.service.db import (
    Approval,
    Job,
    Plan,
    Schedule,
    ScheduleVersion,
    Submission,
    audit,
    new_id,
    utcnow,
)
from tessera.service.jobs import NotFoundError

ACTIVE = ("proposed", "edited", "approved")


class ConflictError(TesseraError):
    code = "CONFLICT"
    http_status = 409


def get_schedule(s: Session, sid: str) -> Schedule:
    row = s.get(Schedule, sid)
    if row is None:
        raise NotFoundError(f"schedule {sid} not found")
    return row


def get_version(s: Session, sid: str, version: int | None = None) -> ScheduleVersion:
    sch = get_schedule(s, sid)
    v = sch.current_version if version is None else version
    row = s.scalar(select(ScheduleVersion).where(ScheduleVersion.schedule_id == sid, ScheduleVersion.version == v))
    if row is None:
        raise NotFoundError(f"schedule {sid} version {v} not found")
    return row


def _diff(old: list[dict[str, Any]], new: list[dict[str, Any]]) -> dict[str, Any]:
    def key(v):
        return (v["field_id"], v.get("pattern_step", 0))

    o = {key(v): v for v in old}
    n = {key(v): v for v in new}
    changed = []
    for k in sorted(set(o) & set(n)):
        fields = sorted(f for f in ("start", "end", "filter", "n_exp", "t_exp_s", "ra", "dec")
                        if o[k].get(f) != n[k].get(f))
        if fields:
            changed.append({"field_id": k[0], "pattern_step": k[1], "fields": fields})
    return {"added": [k[0] for k in sorted(set(n) - set(o))], "removed": [k[0] for k in sorted(set(o) - set(n))],
            "changed": changed}


def edit(s: Session, sid: str, base_version: int, visits: list[dict[str, Any]], actor: str) -> ScheduleVersion:
    sch = get_schedule(s, sid)
    if sch.status not in ACTIVE:
        raise ConflictError(f"schedule is {sch.status}; edits are not allowed")
    if base_version != sch.current_version:
        raise ConflictError("edit is based on a stale version", details={"current": sch.current_version})
    base = get_version(s, sid, base_version)
    profile = load_profile(sch.profile_key.split("@")[0])
    if profile.key != sch.profile_key:
        raise ConflictError("profile version changed since the plan; re-plan instead",
                            details={"plan": sch.profile_key, "installed": profile.key})
    clean = []
    for i, v in enumerate(sorted(visits, key=lambda x: x["start"])):
        req = {k: v[k] for k in ("field_id", "ra", "dec", "filter", "n_exp", "t_exp_s", "start", "end")}
        clean.append({**v, **req, "seq": i + 1})
    hz = base.doc["horizon"]
    from tessera.core.timeutil import to_time

    report = validate_visits(clean, profile, horizon=(to_time(hz["start"]).jd, to_time(hz["end"]).jd))
    for v, r in zip(clean, report["visits"]):
        v["feasible"], v["violations"], v["airmass"] = r["feasible"], r["violations"], r["airmass"]
        v["edited"] = True
    doc = json.loads(json.dumps(base.doc))
    doc.update({"version": base_version + 1, "parent_version": base_version, "visits": clean,
                "feasible": report["feasible"],
                "validation": {k: report[k] for k in ("feasible", "n_visits", "n_infeasible", "violations_by_constraint",
                                                      "step_s")},
                "status": "edited" if report["feasible"] else "infeasible", "edited_by": actor})
    row = ScheduleVersion(schedule_id=sid, version=base_version + 1, parent_version=base_version, doc=doc,
                          sha256=schedule_sha256(doc), diff=_diff(base.doc["visits"], clean), created_by=actor)
    s.add(row)
    sch.current_version = base_version + 1
    sch.status = doc["status"] if report["feasible"] else "edited"
    audit(s, actor, "schedule_edit", sid, version=row.version, feasible=report["feasible"])
    s.commit()
    return row


def approve(s: Session, sid: str, version: int, actor: str, *, auto: bool = False) -> Approval:
    sch = get_schedule(s, sid)
    if sch.status not in ACTIVE:
        raise ConflictError(f"schedule is {sch.status}; cannot approve")
    if version != sch.current_version:
        raise ConflictError("only the current version can be approved", details={"current": sch.current_version})
    ver = get_version(s, sid, version)
    if not ver.doc.get("feasible"):
        raise ConflictError("infeasible schedules cannot be approved", code="SCHEDULE_INFEASIBLE")
    a = Approval(id=new_id("apr"), schedule_id=sid, version=version, sha256=ver.sha256, approver=actor, auto=auto)
    s.add(a)
    sch.status = "approved"
    audit(s, actor, "schedule_auto_approve" if auto else "schedule_approve", sid, version=version, sha256=ver.sha256)
    s.commit()
    return a


def maybe_auto_approve(s: Session, sid: str, profile_key: str) -> Approval | None:
    """Explicit per-profile admin policy only (spec 10.5); off by default."""
    try:
        profile = load_profile(profile_key.split("@")[0])
    except TesseraError:
        return None
    pol = profile.auto_approve
    if not pol:
        return None
    ver = get_version(s, sid)
    doc = ver.doc
    if not doc.get("feasible"):
        return None
    if len(doc["visits"]) > int(pol.get("max_visits", 0)):
        return None
    plan = s.get(Plan, s.get(Schedule, sid).plan_key)
    preset = (plan.manifest.get("config") or {}).get("preset") if plan else None
    if pol.get("presets") and preset not in pol["presets"]:
        return None
    return approve(s, sid, ver.version, f"auto-approve:{profile.key}", auto=True)


def latest_approval(s: Session, sid: str, version: int) -> Approval | None:
    return s.scalar(select(Approval).where(Approval.schedule_id == sid, Approval.version == version)
                    .order_by(Approval.created.desc()))


def _export_context(s: Session, sch: Schedule, store) -> dict[str, Any]:
    """Extra inputs some adapters need (Rubin: the sky map and galaxy table)."""
    plan = s.get(Plan, sch.plan_key)
    job = s.get(Job, plan.job_id) if plan else None
    ctx: dict[str, Any] = {"event": (plan.manifest.get("event") if plan else None) or {}}
    if job is not None:
        from astropy.table import Table

        from tessera.core.skymap import load_skymap

        ctx["skymap"] = load_skymap(store.get(job.request["skymap_id"]))
        gpath = f"{plan.artifact_dir}/galaxies.ecsv"
        try:
            ctx["galaxies"] = Table.read(gpath, format="ascii.ecsv")
        except Exception:
            ctx["galaxies"] = None
        ctx["graceid"] = ctx["event"].get("superevent_id")
    return ctx


def export(s: Session, sid: str, version: int | None, adapter_name: str, store, *, allow_infeasible: bool = False,
           **opts: Any):
    sch = get_schedule(s, sid)
    ver = get_version(s, sid, version)
    profile = load_profile(sch.profile_key.split("@")[0])
    adapter = get_adapter(adapter_name)
    extra = _export_context(s, sch, store) if adapter_name in ("rubin_too", "treasuremap") else {}
    if adapter_name == "treasuremap":
        extra = {"graceid": extra.get("graceid")}
    exp = adapter.export({**ver.doc, "status": sch.status}, profile, allow_infeasible=allow_infeasible,
                         **extra, **opts)
    issues = adapter.validate(exp)
    if any(i.level == "error" for i in issues):
        raise AdapterError("export failed format validation", details={"issues": [i.__dict__ for i in issues]})
    return exp, issues


def submit(s: Session, sid: str, version: int, adapter_name: str, store, actor: str) -> Submission:
    sch = get_schedule(s, sid)
    if sch.status != "approved":
        raise ConflictError("only approved schedules can be submitted", code="NOT_APPROVED")
    ver = get_version(s, sid, version)
    appr = latest_approval(s, sid, version)
    if appr is None:
        raise ConflictError("this version has no approval", code="NOT_APPROVED")
    profile = load_profile(sch.profile_key.split("@")[0])
    exp, _ = export(s, sid, version, adapter_name, store)
    adapter = get_adapter(adapter_name)
    approval = ApprovalDC(schedule_id=sid, version=version, schedule_sha256=appr.sha256, approver=appr.approver,
                          approved_utc=appr.created.isoformat())
    receipt = adapter.submit(exp, approval, ver.doc, profile)
    sub = Submission(id=new_id("sub"), schedule_id=sid, version=version, adapter=adapter_name, state="submitted",
                     receipt=receipt.__dict__)
    s.add(sub)
    sch.status = "submitted"
    audit(s, actor, "schedule_submit", sid, version=version, adapter=adapter_name, remote_id=receipt.remote_id)
    s.commit()
    return sub


def retire_for_revision(s: Session, revision_id: str, new_status: str, actor: str) -> dict[str, int]:
    """Supersede or withdraw a revision's schedules.  Already-submitted schedules
    are never withdrawn silently: the operator is notified (audit + diagnostics)."""
    counts = {"retired": 0, "operator_notified": 0}
    for sch in s.scalars(select(Schedule).where(Schedule.revision_id == revision_id)).all():
        if sch.status in ACTIVE or sch.status == "infeasible":
            sch.status = new_status
            sch.updated = utcnow()
            counts["retired"] += 1
            audit(s, actor, f"schedule_{new_status}", sch.id)
        elif sch.status in ("submitted", "executing"):
            counts["operator_notified"] += 1
            audit(s, actor, "operator_notify", sch.id, reason=f"revision {new_status}; schedule already {sch.status}")
    s.commit()
    return counts

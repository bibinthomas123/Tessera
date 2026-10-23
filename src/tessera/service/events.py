"""Alert handling and revision-aware replanning (spec 8.2-8.5).

One LVK superevent = one Event with ordered Revisions.  Only the newest active
revision has live schedules.  Duplicates, out-of-order arrivals and
retractions are idempotent.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from tessera.adapters.alerts.envelope import AlertEnvelope
from tessera.core.errors import TesseraError
from tessera.service.db import AlertRaw, Event, Job, Revision, audit, new_id
from tessera.service.jobs import bucket_time, enqueue, request_cancel, store_skymap
from tessera.service.schedules import retire_for_revision
from tessera.service.settings import Settings
from tessera.service.store import ArtifactStore


@dataclass
class PlanRule:
    preset: str
    profiles: list[str]
    catalog: str | None = None
    when: dict[str, Any] = field(default_factory=dict)
    config: dict[str, Any] = field(default_factory=dict)

    def matches(self, env: AlertEnvelope) -> bool:
        w = self.when
        props, cls = env.properties, env.classification
        if "min_has_ns" in w and props.get("HasNS", 0.0) < w["min_has_ns"]:
            return False
        if "max_has_ns" in w and props.get("HasNS", 0.0) > w["max_has_ns"]:
            return False
        if "min_has_remnant" in w and props.get("HasRemnant", 0.0) < w["min_has_remnant"]:
            return False
        if "min_p_bbh" in w and cls.get("BBH", 0.0) < w["min_p_bbh"]:
            return False
        if "max_p_terrestrial" in w and cls.get("Terrestrial", 0.0) > w["max_p_terrestrial"]:
            return False
        if "groups" in w and env.group not in w["groups"]:
            return False
        return True


@dataclass
class TriggerPolicy:
    require_significant: bool = True
    include_mock: bool = False
    alert_types: tuple[str, ...] = ("PRELIMINARY", "INITIAL", "UPDATE")
    rules: list[PlanRule] = field(default_factory=list)

    @classmethod
    def load(cls, path: Path | None, mode: str) -> TriggerPolicy:
        if path is None:
            return cls(include_mock=(mode == "test"))
        doc = yaml.safe_load(Path(path).read_text()) or {}
        rules = [PlanRule(**r) for r in doc.get("rules", [])]
        return cls(require_significant=doc.get("require_significant", True),
                   include_mock=doc.get("include_mock", mode == "test"),
                   alert_types=tuple(doc.get("alert_types", cls.alert_types)), rules=rules)

    def rules_for(self, env: AlertEnvelope) -> tuple[list[PlanRule], str | None]:
        if env.is_mock and not self.include_mock:
            return [], "mock event ignored outside test mode"
        if env.alert_type not in self.alert_types:
            return [], f"alert type {env.alert_type} not planned by policy"
        if self.require_significant and env.significant is False:
            return [], "event not significant"
        rules = [r for r in self.rules if r.matches(env)]
        return rules, None if rules else "no trigger rule matched"


@dataclass
class AlertOutcome:
    outcome: str
    superevent_id: str
    revision_id: str | None = None
    jobs: list[str] = field(default_factory=list)
    detail: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return self.__dict__


def _t(iso: str) -> float:
    from tessera.core.timeutil import to_time

    return float(to_time(iso).jd)


def _cancel_revision_jobs(s: Session, revision_id: str, reason: str) -> int:
    n = 0
    for j in s.scalars(select(Job).where(Job.revision_id == revision_id, Job.state.in_(("queued", "running")))).all():
        request_cancel(s, j.id, "system", reason)
        n += 1
    return n


def handle_alert(s: Session, settings: Settings, store: ArtifactStore, env: AlertEnvelope,
                 policy: TriggerPolicy) -> AlertOutcome:
    sid = env.superevent_id
    if s.get(AlertRaw, env.alert_key) is not None:
        return AlertOutcome("duplicate", sid)
    raw_digest = store.put(env.raw)
    alert = AlertRaw(alert_key=env.alert_key, superevent_id=sid, alert_type=env.alert_type,
                     time_created=env.time_created, source=env.source, payload_sha256=raw_digest)
    s.add(alert)
    ev = s.get(Event, sid)
    if ev is None:
        ev = Event(superevent_id=sid, state="active", meta={})
        s.add(ev)
        s.flush()

    if env.alert_type == "RETRACTION":
        if ev.state == "retracted":
            alert.outcome = "retraction_duplicate"
            s.commit()
            return AlertOutcome("retraction_duplicate", sid)
        ev.state = "retracted"
        n_jobs = 0
        retired: dict[str, int] = {"retired": 0, "operator_notified": 0}
        for rev in s.scalars(select(Revision).where(Revision.superevent_id == sid)).all():
            n_jobs += _cancel_revision_jobs(s, rev.id, "retracted")
            c = retire_for_revision(s, rev.id, "withdrawn", "system")
            retired = {k: retired[k] + c[k] for k in retired}
        alert.outcome = "retracted"
        audit(s, "system", "event_retracted", sid, jobs_cancelled=n_jobs, **retired)
        s.commit()
        return AlertOutcome("retracted", sid, detail=f"cancelled {n_jobs} jobs; {retired}")

    if not env.skymap_bytes:
        alert.outcome = "no_skymap"
        s.commit()
        return AlertOutcome("no_skymap", sid, detail="alert carried no sky map")

    sk = store_skymap(s, store, env.skymap_bytes, env.skymap_filename or f"{sid}.fits", f"alert:{env.source}")
    rev = s.scalar(select(Revision).where(Revision.superevent_id == sid, Revision.skymap_sha256 == sk.id))
    if rev is not None:
        types = list(rev.alert_types)
        if env.alert_type not in types:
            rev.alert_types = types + [env.alert_type]
        rev.meta = {**rev.meta, "latest_alert": env.metadata()}
        alert.outcome = "metadata_updated"
        s.commit()
        return AlertOutcome("metadata_updated", sid, rev.id, detail="same sky map under a new alert; no replan")

    number = (s.scalar(select(func.count()).select_from(Revision).where(Revision.superevent_id == sid)) or 0) + 1
    rev = Revision(id=new_id("rev"), superevent_id=sid, number=number, skymap_sha256=sk.id, skymap_id=sk.id,
                   alert_types=[env.alert_type], time_created=env.time_created, state="pending",
                   meta={"alert": env.metadata()})
    s.add(rev)
    if ev.state == "retracted":
        rev.state = "after_retraction"
        alert.outcome = "after_retraction"
        audit(s, "system", "operator_alert", sid, reason="revision received after retraction", revision=rev.id)
        s.commit()
        return AlertOutcome("after_retraction", sid, rev.id, detail="stored, not planned; operator alerted")

    active = s.get(Revision, ev.active_revision_id) if ev.active_revision_id else None
    if active is not None and _t(env.time_created) <= _t(active.time_created):
        rev.state = "superseded"
        alert.outcome = "late"
        s.commit()
        return AlertOutcome("late", sid, rev.id, detail="older than the active revision; stored as superseded")

    if active is not None:
        active.state = "superseded"
        if not active.pinned:
            _cancel_revision_jobs(s, active.id, "superseded")
        retire_for_revision(s, active.id, "superseded", "system")
    rev.state = "active"
    ev.active_revision_id = rev.id
    ev.meta = {**ev.meta, "latest_alert": env.metadata()}
    s.commit()

    rules, why_not = policy.rules_for(env)
    if not rules:
        alert.outcome = "not_planned"
        s.commit()
        return AlertOutcome("not_planned", sid, rev.id, detail=why_not)
    jobs: list[str] = []
    t_start = bucket_time(env.time_created, 5)
    event = {"superevent_id": sid, "revision": rev.number, "alert_type": env.alert_type,
             "time_created": env.time_created, "classification": env.classification, "properties": env.properties}
    for r in rules:
        req = {"skymap_id": sk.id, "catalog": r.catalog, "profiles": r.profiles,
               "config": {"preset": r.preset, **r.config}, "t_start": t_start, "event": event}
        try:
            job, plan, _ = enqueue(s, settings, req, owner="system", revision_id=rev.id,
                                   idempotency_key=f"{env.alert_key}:{r.preset}:{','.join(r.profiles)}")
        except TesseraError as exc:
            audit(s, "system", "enqueue_failed", sid, error=exc.code, detail=exc.message)
            s.commit()
            continue
        jobs.append(job.id if job else plan.plan_key)
    alert.outcome = "planned"
    s.commit()
    return AlertOutcome("planned", sid, rev.id, jobs)

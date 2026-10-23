"""Persistence (SQLAlchemy 2.0).  PostgreSQL in production; SQLite for a
single-node install and for tests.  Plans and schedule versions are immutable."""

from __future__ import annotations

import uuid
from datetime import datetime, timezone
from typing import Any

from sqlalchemy import (
    JSON,
    Boolean,
    DateTime,
    Float,
    ForeignKey,
    Integer,
    String,
    Text,
    UniqueConstraint,
    create_engine,
    event,
    select,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, Session, mapped_column, sessionmaker


def utcnow() -> datetime:
    return datetime.now(timezone.utc)


def new_id(prefix: str) -> str:
    return f"{prefix}_{uuid.uuid4().hex[:20]}"


class Base(DeclarativeBase):
    type_annotation_map = {dict[str, Any]: JSON, list[Any]: JSON}


class ApiKey(Base):
    __tablename__ = "api_keys"
    id: Mapped[str] = mapped_column(String(40), primary_key=True)
    name: Mapped[str] = mapped_column(String(200))
    prefix: Mapped[str] = mapped_column(String(16), index=True, unique=True)
    hash: Mapped[str] = mapped_column(Text)
    scopes: Mapped[str] = mapped_column(String(200))
    created: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    revoked: Mapped[bool] = mapped_column(Boolean, default=False)


class Event(Base):
    __tablename__ = "events"
    superevent_id: Mapped[str] = mapped_column(String(64), primary_key=True)
    state: Mapped[str] = mapped_column(String(16), default="active")         # active | retracted
    active_revision_id: Mapped[str | None] = mapped_column(String(40), nullable=True)
    created: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow, onupdate=utcnow)
    meta: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)


class AlertRaw(Base):
    __tablename__ = "alerts_raw"
    alert_key: Mapped[str] = mapped_column(String(64), primary_key=True)
    superevent_id: Mapped[str] = mapped_column(String(64), index=True)
    alert_type: Mapped[str] = mapped_column(String(16))
    time_created: Mapped[str] = mapped_column(String(40))
    received: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    source: Mapped[str] = mapped_column(String(64))
    payload_sha256: Mapped[str] = mapped_column(String(64))
    outcome: Mapped[str] = mapped_column(String(40), default="")


class Revision(Base):
    __tablename__ = "revisions"
    __table_args__ = (UniqueConstraint("superevent_id", "skymap_sha256"),)
    id: Mapped[str] = mapped_column(String(40), primary_key=True)
    superevent_id: Mapped[str] = mapped_column(String(64), ForeignKey("events.superevent_id"), index=True)
    number: Mapped[int] = mapped_column(Integer)
    skymap_sha256: Mapped[str] = mapped_column(String(64))
    skymap_id: Mapped[str] = mapped_column(String(64))
    alert_types: Mapped[list[Any]] = mapped_column(JSON, default=list)
    time_created: Mapped[str] = mapped_column(String(40))
    state: Mapped[str] = mapped_column(String(24))      # active | superseded | after_retraction
    pinned: Mapped[bool] = mapped_column(Boolean, default=False)
    meta: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    created: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class Skymap(Base):
    __tablename__ = "skymaps"
    id: Mapped[str] = mapped_column(String(64), primary_key=True)        # sha256 of raw bytes
    filename: Mapped[str] = mapped_column(String(255))
    sha256_data: Mapped[str] = mapped_column(String(64))
    summary: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    uploaded_by: Mapped[str] = mapped_column(String(200), default="")
    created: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class Job(Base):
    __tablename__ = "jobs"
    id: Mapped[str] = mapped_column(String(40), primary_key=True)
    plan_key: Mapped[str] = mapped_column(String(80), index=True)
    idempotency_key: Mapped[str | None] = mapped_column(String(200), nullable=True, index=True)
    owner: Mapped[str] = mapped_column(String(200), default="")
    state: Mapped[str] = mapped_column(String(16), index=True, default="queued")
    stage: Mapped[str | None] = mapped_column(String(32), nullable=True)
    request: Mapped[dict[str, Any]] = mapped_column(JSON)
    revision_id: Mapped[str | None] = mapped_column(String(40), nullable=True, index=True)
    attempts: Mapped[int] = mapped_column(Integer, default=0)
    not_before: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    lease_owner: Mapped[str | None] = mapped_column(String(80), nullable=True)
    lease_expires: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    cancel_requested: Mapped[bool] = mapped_column(Boolean, default=False)
    error: Mapped[dict[str, Any] | None] = mapped_column(JSON, nullable=True)
    result: Mapped[dict[str, Any] | None] = mapped_column(JSON, nullable=True)
    created: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    started: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    finished: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    duration_s: Mapped[float | None] = mapped_column(Float, nullable=True)


class Plan(Base):
    __tablename__ = "plans"
    plan_key: Mapped[str] = mapped_column(String(80), primary_key=True)
    job_id: Mapped[str] = mapped_column(String(40))
    revision_id: Mapped[str | None] = mapped_column(String(40), nullable=True, index=True)
    artifact_dir: Mapped[str] = mapped_column(Text)
    summary: Mapped[dict[str, Any]] = mapped_column(JSON)
    manifest: Mapped[dict[str, Any]] = mapped_column(JSON)
    outputs: Mapped[dict[str, Any]] = mapped_column(JSON)          # relpath -> sha256
    created: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class Schedule(Base):
    __tablename__ = "schedules"
    id: Mapped[str] = mapped_column(String(40), primary_key=True)
    plan_key: Mapped[str] = mapped_column(String(80), index=True)
    profile_key: Mapped[str] = mapped_column(String(80))
    revision_id: Mapped[str | None] = mapped_column(String(40), nullable=True, index=True)
    status: Mapped[str] = mapped_column(String(16))   # proposed|edited|approved|submitted|executing|completed|withdrawn|superseded|infeasible
    current_version: Mapped[int] = mapped_column(Integer, default=0)
    created: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow, onupdate=utcnow)


class ScheduleVersion(Base):
    __tablename__ = "schedule_versions"
    __table_args__ = (UniqueConstraint("schedule_id", "version"),)
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    schedule_id: Mapped[str] = mapped_column(String(40), ForeignKey("schedules.id"), index=True)
    version: Mapped[int] = mapped_column(Integer)
    parent_version: Mapped[int | None] = mapped_column(Integer, nullable=True)
    doc: Mapped[dict[str, Any]] = mapped_column(JSON)
    sha256: Mapped[str] = mapped_column(String(64))
    diff: Mapped[dict[str, Any] | None] = mapped_column(JSON, nullable=True)
    created_by: Mapped[str] = mapped_column(String(200), default="tessera")
    created: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class Approval(Base):
    __tablename__ = "approvals"
    id: Mapped[str] = mapped_column(String(40), primary_key=True)
    schedule_id: Mapped[str] = mapped_column(String(40), index=True)
    version: Mapped[int] = mapped_column(Integer)
    sha256: Mapped[str] = mapped_column(String(64))
    approver: Mapped[str] = mapped_column(String(200))
    auto: Mapped[bool] = mapped_column(Boolean, default=False)
    created: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class Submission(Base):
    __tablename__ = "submissions"
    id: Mapped[str] = mapped_column(String(40), primary_key=True)
    schedule_id: Mapped[str] = mapped_column(String(40), index=True)
    version: Mapped[int] = mapped_column(Integer)
    adapter: Mapped[str] = mapped_column(String(40))
    state: Mapped[str] = mapped_column(String(16))
    receipt: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    created: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class Observation(Base):
    __tablename__ = "observations"
    id: Mapped[str] = mapped_column(String(40), primary_key=True)
    superevent_id: Mapped[str | None] = mapped_column(String(64), index=True, nullable=True)
    schedule_id: Mapped[str | None] = mapped_column(String(40), nullable=True)
    profile_id: Mapped[str | None] = mapped_column(String(80), nullable=True)
    status: Mapped[str] = mapped_column(String(16))
    payload: Mapped[dict[str, Any]] = mapped_column(JSON)
    reported_by: Mapped[str] = mapped_column(String(200), default="")
    created: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class AuditLog(Base):
    __tablename__ = "audit_log"
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    ts: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow, index=True)
    actor: Mapped[str] = mapped_column(String(200))
    action: Mapped[str] = mapped_column(String(64))
    target: Mapped[str] = mapped_column(String(200))
    details: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)


_engines: dict[str, Any] = {}


def get_engine(url: str):
    if url not in _engines:
        kwargs: dict[str, Any] = {"future": True, "pool_pre_ping": True}
        if url.startswith("sqlite"):
            kwargs["connect_args"] = {"check_same_thread": False, "timeout": 30}
        eng = create_engine(url, **kwargs)
        if url.startswith("sqlite"):
            @event.listens_for(eng, "connect")
            def _pragma(dbapi_con, _):  # pragma: no cover - trivial
                cur = dbapi_con.cursor()
                cur.execute("PRAGMA journal_mode=WAL")
                cur.execute("PRAGMA busy_timeout=30000")
                cur.close()
        Base.metadata.create_all(eng)
        _engines[url] = eng
    return _engines[url]


def session_factory(url: str) -> sessionmaker[Session]:
    return sessionmaker(bind=get_engine(url), expire_on_commit=False, future=True)


def audit(s: Session, actor: str, action: str, target: str, **details: Any) -> None:
    s.add(AuditLog(actor=actor, action=action, target=target, details=details))


__all__ = ["Base", "ApiKey", "Event", "AlertRaw", "Revision", "Skymap", "Job", "Plan", "Schedule", "ScheduleVersion",
           "Approval", "Submission", "Observation", "AuditLog", "get_engine", "session_factory", "audit", "select",
           "new_id", "utcnow"]

"""Scheduler adapter interface (spec 9.1).

Adapters translate the facility-agnostic ``tessera.schedule/1`` document.  They
refuse infeasible schedules, and ``submit`` refuses anything without an
approval bound to the exact schedule version (spec 10.5).
"""

from __future__ import annotations

import io
import zipfile
from dataclasses import dataclass, field
from typing import Any, Protocol

from tessera.core.errors import AdapterError
from tessera.core.hashing import sha256_json
from tessera.core.profile import TelescopeProfile


@dataclass
class Export:
    adapter: str
    adapter_version: str
    files: dict[str, bytes]
    media_type: str

    def as_bytes(self) -> tuple[bytes, str, str]:
        """(payload, media_type, filename) — a zip when there are several files."""
        if len(self.files) == 1:
            name, data = next(iter(self.files.items()))
            return data, self.media_type, name
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for name in sorted(self.files):
                zi = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                z.writestr(zi, self.files[name])
        return buf.getvalue(), "application/zip", f"{self.adapter}.zip"


@dataclass
class Approval:
    schedule_id: str
    version: int
    schedule_sha256: str
    approver: str
    approved_utc: str


@dataclass
class SubmitReceipt:
    adapter: str
    remote_id: str | None
    submitted_utc: str
    details: dict[str, Any] = field(default_factory=dict)


@dataclass
class Issue:
    level: str
    message: str


class SchedulerAdapter(Protocol):
    name: str
    version: str
    requires: set[str]

    def export(self, schedule: dict[str, Any], profile: TelescopeProfile, **opts: Any) -> Export: ...

    def validate(self, export: Export) -> list[Issue]: ...


class BaseAdapter:
    name = "base"
    version = "1"
    requires: set[str] = set()
    media_type = "application/octet-stream"
    can_submit = False

    def check_exportable(self, schedule: dict[str, Any], *, allow_infeasible: bool = False) -> None:
        if schedule.get("schema") != "tessera.schedule/1":
            raise AdapterError("not a tessera.schedule/1 document", code="ADAPTER_BAD_SCHEDULE")
        if not schedule.get("feasible", False) and not allow_infeasible:
            raise AdapterError("schedule is infeasible; export refused (spec 5.4)", code="SCHEDULE_INFEASIBLE")
        if schedule.get("status") in ("withdrawn", "superseded"):
            raise AdapterError(f"schedule is {schedule['status']}; export refused", code="SCHEDULE_NOT_ACTIVE")

    def options(self, profile: TelescopeProfile, **opts: Any) -> dict[str, Any]:
        merged = dict(profile.adapter_options.get(self.name, {}))
        merged.update({k: v for k, v in opts.items() if v is not None})
        missing = [r for r in self.requires if r not in merged]
        if missing:
            raise AdapterError(f"adapter {self.name!r} needs profile.adapter_options.{self.name}: {missing}",
                               code="ADAPTER_OPTIONS_MISSING")
        return merged

    def header_comment(self, schedule: dict[str, Any]) -> str:
        return (f"tessera adapter={self.name} v{self.version} plan_key={schedule.get('plan_key')} "
                f"profile={schedule.get('profile')} schedule_version={schedule.get('version')}")

    def validate(self, export: Export) -> list[Issue]:
        return []

    def submit(self, export: Export, approval: Approval | None, schedule: dict[str, Any],
               profile: TelescopeProfile) -> SubmitReceipt:
        check_submission(approval, schedule, profile)
        raise AdapterError(f"adapter {self.name!r} does not support submission; deliver the exported file",
                           code="ADAPTER_NO_SUBMIT")


def check_submission(approval: Approval | None, schedule: dict[str, Any], profile: TelescopeProfile) -> None:
    if not profile.submission_enabled:
        raise AdapterError(f"submission is disabled for profile {profile.key}", code="SUBMISSION_DISABLED")
    if approval is None:
        raise AdapterError("submission requires a human approval", code="NOT_APPROVED")
    if approval.schedule_sha256 != schedule_sha256(schedule):
        raise AdapterError("approval does not match this schedule version", code="APPROVAL_MISMATCH")
    if not schedule.get("feasible"):
        raise AdapterError("schedule is infeasible", code="SCHEDULE_INFEASIBLE")


def schedule_sha256(schedule: dict[str, Any]) -> str:
    """Hash of a schedule version, excluding mutable lifecycle fields."""
    doc = {k: v for k, v in schedule.items() if k not in ("status", "approval", "submission")}
    return sha256_json(doc)

"""Typed errors and diagnostics.

Every failure carries a stable machine-readable ``code`` (documented in
``docs/errors.md``).  Warnings never abort a run; they travel with the result
as :class:`Diagnostic` records so nothing is dropped silently.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal


class TesseraError(Exception):
    """Base class.  ``code`` is stable across releases."""

    code = "TESSERA_ERROR"
    http_status = 500

    def __init__(self, message: str, *, code: str | None = None, details: dict[str, Any] | None = None):
        super().__init__(message)
        if code is not None:
            self.code = code
        self.message = message
        self.details = details or {}

    def to_problem(self) -> dict[str, Any]:
        """RFC 9457 problem+json body."""
        return {
            "type": f"https://tessera.readthedocs.io/errors#{self.code.lower()}",
            "title": self.code,
            "status": self.http_status,
            "detail": self.message,
            "code": self.code,
            **({"details": self.details} if self.details else {}),
        }


class SkymapError(TesseraError):
    code = "SKYMAP_INVALID_FORMAT"
    http_status = 422


class CatalogError(TesseraError):
    code = "CATALOG_ERROR"
    http_status = 422


class ProfileError(TesseraError):
    code = "PROFILE_INVALID"
    http_status = 422


class ConfigError(TesseraError):
    code = "CONFIG_INVALID"
    http_status = 422


class ScheduleError(TesseraError):
    code = "SCHEDULE_INFEASIBLE"
    http_status = 409


class AdapterError(TesseraError):
    code = "ADAPTER_ERROR"
    http_status = 422


class CancelledError(TesseraError):
    code = "JOB_CANCELLED"
    http_status = 409


@dataclass(frozen=True)
class Diagnostic:
    """A warning or informational note attached to a result."""

    code: str
    message: str
    level: Literal["info", "warning", "error"] = "warning"
    details: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {"code": self.code, "level": self.level, "message": self.message, "details": self.details}


class Diagnostics(list):
    """List of :class:`Diagnostic` with a convenience adder."""

    def add(self, code: str, message: str, level: Literal["info", "warning", "error"] = "warning", **details: Any) -> None:
        self.append(Diagnostic(code=code, message=message, level=level, details=details))

    def codes(self) -> list[str]:
        return [d.code for d in self]

    def to_list(self) -> list[dict[str, Any]]:
        return [d.to_dict() for d in self]

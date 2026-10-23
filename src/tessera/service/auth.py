"""Authentication and scopes (spec 10.3).

* API keys ``tessera_<prefix>_<secret>``: stored as argon2 hashes, looked up by
  prefix, carried in ``Authorization: Bearer`` or ``X-API-Key``.
* Institutional gateway mode: identity and scopes forwarded in
  ``X-Tessera-User`` / ``X-Tessera-Scopes``, honoured only from ``TESSERA_TRUSTED_PROXIES``.
* ``TESSERA_AUTH_DISABLED=1`` for a single-user local install (logged loudly).
"""

from __future__ import annotations

import secrets
from dataclasses import dataclass

from argon2 import PasswordHasher
from argon2.exceptions import VerificationError
from sqlalchemy import select
from sqlalchemy.orm import Session

from tessera.core.errors import TesseraError
from tessera.service.db import ApiKey, new_id

SCOPES = ("read", "plan", "approve", "admin", "metrics")
_ph = PasswordHasher()


class AuthError(TesseraError):
    code = "UNAUTHORIZED"
    http_status = 401


class ForbiddenError(TesseraError):
    code = "FORBIDDEN"
    http_status = 403


@dataclass(frozen=True)
class Principal:
    name: str
    scopes: frozenset[str]
    via: str

    def require(self, scope: str) -> None:
        if scope not in self.scopes and "admin" not in self.scopes:
            raise ForbiddenError(f"scope {scope!r} required", details={"have": sorted(self.scopes)})


def create_key(s: Session, name: str, scopes: list[str]) -> str:
    bad = [x for x in scopes if x not in SCOPES]
    if bad:
        raise ValueError(f"unknown scopes {bad}; valid: {SCOPES}")
    prefix = secrets.token_hex(6)
    secret = secrets.token_urlsafe(32)
    s.add(ApiKey(id=new_id("key"), name=name, prefix=prefix, hash=_ph.hash(secret), scopes=",".join(sorted(scopes))))
    s.commit()
    return f"tessera_{prefix}_{secret}"


def verify_key(s: Session, token: str) -> Principal:
    try:
        _, prefix, secret = token.split("_", 2)
    except ValueError as exc:
        raise AuthError("malformed API key") from exc
    row = s.scalar(select(ApiKey).where(ApiKey.prefix == prefix))
    if row is None or row.revoked:
        raise AuthError("unknown or revoked API key")
    try:
        _ph.verify(row.hash, secret)
    except VerificationError as exc:
        raise AuthError("invalid API key") from exc
    return Principal(name=row.name, scopes=frozenset(row.scopes.split(",")), via="api_key")


def revoke_key(s: Session, prefix: str) -> bool:
    row = s.scalar(select(ApiKey).where(ApiKey.prefix == prefix))
    if row is None:
        return False
    row.revoked = True
    s.commit()
    return True

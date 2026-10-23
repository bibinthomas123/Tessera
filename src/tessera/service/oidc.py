"""In-process OIDC bearer-token validation (spec 10.3).

Enabled by ``TESSERA_OIDC_ISSUER`` and ``TESSERA_OIDC_AUDIENCE``.  Signing keys come
from ``TESSERA_OIDC_JWKS_URL``, from the issuer's discovery document
(``<issuer>/.well-known/openid-configuration``), or from a static
``TESSERA_OIDC_JWKS_FILE`` for air-gapped sites.  Tessera scopes are read from
``TESSERA_OIDC_SCOPES_CLAIM`` (default ``scope``, space-separated, or a list);
``read``, ``tessera:read`` and ``tessera.read`` are all accepted.  Signature,
issuer, audience, expiry and not-before are always checked.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from functools import cached_property

import jwt

from tessera.service.auth import SCOPES, AuthError, Principal

ALGORITHMS = ("RS256", "RS384", "RS512", "ES256", "ES384", "PS256")


@dataclass
class OIDCConfig:
    issuer: str
    audience: str
    jwks_url: str | None = None
    jwks_file: str | None = None
    scopes_claim: str = "scope"
    user_claim: str = "preferred_username"
    leeway_s: int = 30

    @classmethod
    def from_env(cls) -> OIDCConfig | None:
        iss, aud = os.environ.get("TESSERA_OIDC_ISSUER"), os.environ.get("TESSERA_OIDC_AUDIENCE")
        if not iss or not aud:
            return None
        return cls(issuer=iss, audience=aud, jwks_url=os.environ.get("TESSERA_OIDC_JWKS_URL"),
                   jwks_file=os.environ.get("TESSERA_OIDC_JWKS_FILE"),
                   scopes_claim=os.environ.get("TESSERA_OIDC_SCOPES_CLAIM", "scope"),
                   user_claim=os.environ.get("TESSERA_OIDC_USER_CLAIM", "preferred_username"))


class OIDCVerifier:
    def __init__(self, cfg: OIDCConfig):
        self.cfg = cfg

    @cached_property
    def _static_keys(self) -> jwt.PyJWKSet | None:
        if not self.cfg.jwks_file:
            return None
        with open(self.cfg.jwks_file) as f:
            return jwt.PyJWKSet.from_dict(json.load(f))

    @cached_property
    def _client(self) -> jwt.PyJWKClient:
        url = self.cfg.jwks_url
        if not url:
            import httpx

            disc = httpx.get(self.cfg.issuer.rstrip("/") + "/.well-known/openid-configuration", timeout=10)
            disc.raise_for_status()
            url = disc.json()["jwks_uri"]
        return jwt.PyJWKClient(url, cache_keys=True, lifespan=3600)

    def _key(self, token: str):
        if self._static_keys is not None:
            kid = jwt.get_unverified_header(token).get("kid")
            for k in self._static_keys.keys:
                if kid is None or k.key_id == kid:
                    return k.key
            raise AuthError("no JWKS key matches the token's kid")
        return self._client.get_signing_key_from_jwt(token).key

    def verify(self, token: str) -> Principal:
        try:
            claims = jwt.decode(token, self._key(token), algorithms=list(ALGORITHMS), audience=self.cfg.audience,
                                issuer=self.cfg.issuer, leeway=self.cfg.leeway_s,
                                options={"require": ["exp", "iss", "aud"]})
        except AuthError:
            raise
        except jwt.PyJWTError as exc:
            raise AuthError(f"invalid bearer token: {exc}") from exc
        except Exception as exc:  # JWKS fetch failures etc.
            raise AuthError(f"cannot validate bearer token: {exc}") from exc
        raw = claims.get(self.cfg.scopes_claim, [])
        items = raw.split() if isinstance(raw, str) else list(raw or [])
        scopes = set()
        for it in items:
            name = str(it).split(":", 1)[-1].split(".", 1)[-1] if str(it).startswith(("tessera:", "tessera.")) else it
            if name in SCOPES:
                scopes.add(name)
        user = str(claims.get(self.cfg.user_claim) or claims.get("sub") or "unknown")
        return Principal(name=f"oidc:{user}", scopes=frozenset(scopes), via="oidc")

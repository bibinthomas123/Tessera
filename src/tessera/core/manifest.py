"""Reproducibility manifest and plan identity (spec 7.5)."""

from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone
from importlib import metadata
from typing import Any

from tessera.core.hashing import sha256_json

TOLERANCE = {
    "scores_rel": 1e-9,
    "ranks": "identical except ties within scores_rel",
    "schedules": "identical (same inputs, same t_start, same container)",
    "cross_platform": "scores within scores_rel; schedule times within 1 s",
}

_PACKAGES = ("tessera", "numpy", "scipy", "astropy", "astroplan", "astropy-healpix", "pyarrow", "pydantic")


def software_versions() -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    for p in _PACKAGES:
        try:
            out[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            out[p] = None
    from tessera import SCIENCE_MODEL_VERSION

    out["science_model"] = SCIENCE_MODEL_VERSION
    out["python"] = sys.version.split()[0]
    import os

    out["container_image"] = os.environ.get("TESSERA_IMAGE_DIGEST")
    return out


def plan_identity(*, skymap_sha256: str, snapshot_id: str, profiles: list[dict[str, str]], config: dict[str, Any],
                  t_start: str, event: dict[str, Any] | None, prior_observations: list[dict[str, Any]] | None
                  ) -> tuple[str, dict[str, Any]]:
    from tessera import SCIENCE_MODEL_VERSION, __version__

    ident = {
        "skymap_sha256": skymap_sha256,
        "catalog_snapshot_id": snapshot_id,
        "profiles": sorted(profiles, key=lambda p: p["key"]),
        "config_sha256": sha256_json(config),
        "t_start": t_start,
        "event": {k: (event or {}).get(k) for k in ("superevent_id", "revision")},
        "prior_observations_sha256": sha256_json(prior_observations or []),
        "tessera": __version__,
        "science_model": SCIENCE_MODEL_VERSION,
    }
    return "sha256:" + sha256_json(ident), ident


def build_manifest(*, plan_key: str, identity: dict[str, Any], event: dict[str, Any] | None, skymap: dict[str, Any],
                   catalog: dict[str, Any], profiles: list[dict[str, Any]], config: dict[str, Any], t_start: str,
                   t_end: str, environment: dict[str, Any], outputs: dict[str, str], timings: dict[str, float]
                   ) -> dict[str, Any]:
    return {
        "schema": "tessera.manifest/1",
        "plan_key": plan_key,
        "identity": identity,
        "event": event or {},
        "skymap": skymap,
        "catalog": catalog,
        "profiles": profiles,
        "config": config,
        "config_sha256": identity["config_sha256"],
        "t_start": t_start,
        "t_end": t_end,
        "versions": software_versions(),
        "environment": {**environment, "platform": platform.platform()},
        "outputs": outputs,
        "tolerance": TOLERANCE,
        "timings_s": {k: round(v, 4) for k, v in timings.items()},
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
    }

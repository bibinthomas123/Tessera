"""Service settings from environment variables (12-factor).  See docs/deployment.md."""

from __future__ import annotations

import ipaddress
import os
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path


def _bool(name: str, default: bool = False) -> bool:
    return os.environ.get(name, "1" if default else "0").lower() in ("1", "true", "yes", "on")


@dataclass
class Settings:
    database_url: str = field(default_factory=lambda: os.environ.get("TESSERA_DATABASE_URL", "sqlite:///./tessera-data/tessera.db"))
    data_dir: Path = field(default_factory=lambda: Path(os.environ.get("TESSERA_DATA_DIR", "./tessera-data")))
    catalog_dirs: list[Path] = field(default_factory=lambda: [Path(p) for p in os.environ.get(
        "TESSERA_CATALOG_DIRS", "./catalogs").split(os.pathsep) if p])
    default_catalog: str | None = field(default_factory=lambda: os.environ.get("TESSERA_DEFAULT_CATALOG"))
    auth_disabled: bool = field(default_factory=lambda: _bool("TESSERA_AUTH_DISABLED"))
    trusted_proxies: list[str] = field(default_factory=lambda: [p for p in os.environ.get(
        "TESSERA_TRUSTED_PROXIES", "").split(",") if p])
    root_path: str = field(default_factory=lambda: os.environ.get("TESSERA_ROOT_PATH", ""))
    max_upload_bytes: int = field(default_factory=lambda: int(os.environ.get("TESSERA_MAX_UPLOAD_BYTES", 100 * 1024**2)))
    max_queue_depth: int = field(default_factory=lambda: int(os.environ.get("TESSERA_MAX_QUEUE_DEPTH", 50)))
    max_jobs_per_key: int = field(default_factory=lambda: int(os.environ.get("TESSERA_MAX_JOBS_PER_KEY", 4)))
    job_timeout_s: float = field(default_factory=lambda: float(os.environ.get("TESSERA_JOB_TIMEOUT_S", 300)))
    job_memory_bytes: int = field(default_factory=lambda: int(os.environ.get("TESSERA_JOB_MEMORY_BYTES", 4 * 1024**3)))
    job_max_attempts: int = field(default_factory=lambda: int(os.environ.get("TESSERA_JOB_MAX_ATTEMPTS", 3)))
    lease_s: float = field(default_factory=lambda: float(os.environ.get("TESSERA_LEASE_S", 60)))
    trigger_policy: Path | None = field(default_factory=lambda: (Path(os.environ["TESSERA_TRIGGER_POLICY"])
                                                                if os.environ.get("TESSERA_TRIGGER_POLICY") else None))
    mode: str = field(default_factory=lambda: os.environ.get("TESSERA_MODE", "production"))  # production | test
    worker_metrics_port: int = field(default_factory=lambda: int(os.environ.get("TESSERA_WORKER_METRICS_PORT", 9101)))
    inline_worker: bool = field(default_factory=lambda: _bool("TESSERA_INLINE_WORKER"))

    def is_trusted_proxy(self, host: str | None) -> bool:
        if not host or not self.trusted_proxies:
            return False
        try:
            ip = ipaddress.ip_address(host)
        except ValueError:
            return False
        return any(ip in ipaddress.ip_network(c, strict=False) for c in self.trusted_proxies)

    @property
    def artifacts_dir(self) -> Path:
        return self.data_dir / "artifacts"


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    s = Settings()
    s.data_dir.mkdir(parents=True, exist_ok=True)
    s.artifacts_dir.mkdir(parents=True, exist_ok=True)
    return s

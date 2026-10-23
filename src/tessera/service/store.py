"""Content-addressed artifact store on a filesystem (or a mounted S3-compatible volume)."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

from tessera.core.hashing import sha256_bytes


class ArtifactStore:
    def __init__(self, root: Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def _path(self, digest: str) -> Path:
        return self.root / "blobs" / digest[:2] / digest

    def put(self, data: bytes) -> str:
        digest = sha256_bytes(data)
        p = self._path(digest)
        if not p.exists():
            p.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp = tempfile.mkstemp(dir=p.parent)
            with os.fdopen(fd, "wb") as f:
                f.write(data)
            os.replace(tmp, p)
        return digest

    def get(self, digest: str) -> bytes:
        p = self._path(digest)
        if not p.is_file():
            raise FileNotFoundError(digest)
        return p.read_bytes()

    def path(self, digest: str) -> Path:
        return self._path(digest)

    def exists(self, digest: str) -> bool:
        return self._path(digest).is_file()

    def plan_dir(self, plan_key: str) -> Path:
        d = self.root / "plans" / plan_key.replace("sha256:", "")
        d.mkdir(parents=True, exist_ok=True)
        return d

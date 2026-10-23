"""Canonical JSON and content hashing used for identities and manifests."""

from __future__ import annotations

import hashlib
import json
import math
from datetime import date, datetime
from pathlib import Path
from typing import Any

import numpy as np


def _default(obj: Any) -> Any:
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (datetime, date)):
        return obj.isoformat()
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "model_dump"):
        return obj.model_dump(mode="json")
    if isinstance(obj, (set, frozenset)):
        return sorted(obj)
    raise TypeError(f"not JSON serialisable: {type(obj)!r}")


def _clean(obj: Any) -> Any:
    """Replace non-finite floats (not valid JSON) with strings."""
    if isinstance(obj, float) and not math.isfinite(obj):
        return "inf" if obj > 0 else ("-inf" if obj < 0 else "nan")
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    return obj


def canonical_json(obj: Any) -> str:
    """Deterministic JSON: sorted keys, no whitespace, finite floats only."""
    return json.dumps(_clean(json.loads(json.dumps(obj, default=_default))), sort_keys=True, separators=(",", ":"))


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_json(obj: Any) -> str:
    return sha256_bytes(canonical_json(obj).encode())


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            block = f.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()

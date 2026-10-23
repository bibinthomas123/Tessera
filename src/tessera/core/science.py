"""Access to the versioned science-model defaults (``data/science/defaults.yaml``)."""

from __future__ import annotations

import copy
from functools import lru_cache
from importlib import resources
from typing import Any

import yaml

from tessera.core.errors import ConfigError


@lru_cache(maxsize=1)
def _load() -> dict[str, Any]:
    text = resources.files("tessera.data.science").joinpath("defaults.yaml").read_text()
    return yaml.safe_load(text)


def defaults() -> dict[str, Any]:
    """A deep copy of the science defaults (safe to mutate)."""
    return copy.deepcopy(_load())


def preset(name: str) -> dict[str, Any]:
    p = _load()["presets"].get(name)
    if p is None:
        raise ConfigError(f"unknown preset {name!r}", details={"known": sorted(_load()["presets"])})
    return copy.deepcopy(p)


def strategy(name: str) -> dict[str, Any]:
    s = _load()["strategies"].get(name)
    if s is None:
        raise ConfigError(f"unknown strategy {name!r}", details={"known": sorted(_load()["strategies"])})
    return copy.deepcopy(s)


def light_curve(name: str) -> dict[str, Any]:
    lc = _load()["light_curves"].get(name)
    if lc is None:
        raise ConfigError(f"unknown light curve {name!r}", details={"known": sorted(_load()["light_curves"])})
    return copy.deepcopy(lc)


def reference_density(weight: str) -> dict[str, Any] | None:
    return copy.deepcopy(_load()["reference_density"].get(weight))

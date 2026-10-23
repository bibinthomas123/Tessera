"""Scheduler adapter registry.  Third-party adapters register under the
``tessera.adapters`` entry-point group."""

from __future__ import annotations

from importlib.metadata import entry_points

from tessera.core.errors import AdapterError

_BUILTIN = {
    "generic": "tessera.adapters.schedulers.generic:GenericAdapter",
    "acp": "tessera.adapters.schedulers.acp:ACPAdapter",
    "lco": "tessera.adapters.schedulers.lco:LCOAdapter",
    "treasuremap": "tessera.adapters.schedulers.treasuremap:TreasureMapAdapter",
    "rubin_too": "tessera.adapters.schedulers.rubin_too:RubinToOAdapter",
}


def _load(spec: str):
    import importlib

    mod, cls = spec.split(":")
    return getattr(importlib.import_module(mod), cls)


def available() -> list[str]:
    names = set(_BUILTIN)
    try:
        names.update(ep.name for ep in entry_points(group="tessera.adapters"))
    except Exception:
        pass
    return sorted(names)


def get_adapter(name: str):
    if name in _BUILTIN:
        return _load(_BUILTIN[name])()
    for ep in entry_points(group="tessera.adapters"):
        if ep.name == name:
            return ep.load()()
    raise AdapterError(f"unknown adapter {name!r}", code="ADAPTER_UNKNOWN", details={"available": available()})

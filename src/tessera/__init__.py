"""tessera: a telescope-agnostic gravitational-wave follow-up engine.

The scientific core (:mod:`tessera.core`) is pure and deterministic: it has no
network, database or web dependencies and can be imported directly by
follow-up pipelines.  The service (:mod:`tessera.service`) and adapters
(:mod:`tessera.adapters`) wrap it.

Typical use::

    from tessera import plan, load_skymap, open_snapshot, load_profile, PlanConfig

    result = plan(
        skymap=load_skymap("bayestar.multiorder.fits"),
        catalog=open_snapshot("catalogs/gladeplus-2022.10"),
        profiles=[load_profile("example-1m")],
        config=PlanConfig(preset="bns"),
        t_start="2017-08-17T22:00:00",
    )
    result.galaxies.write("galaxies.csv")
"""

__version__ = "1.2.0"

#: Version of the science model (formulas, priors, defaults).  Major bumps
#: change rankings; minor bumps add optional strategies; patch bumps are
#: numerical fixes within the documented tolerance.
SCIENCE_MODEL_VERSION = "1.2.0"

_LAZY = {
    "plan": ("tessera.core.plan", "plan"),
    "PlanResult": ("tessera.core.plan", "PlanResult"),
    "PlanConfig": ("tessera.core.config", "PlanConfig"),
    "load_skymap": ("tessera.core.skymap", "load_skymap"),
    "open_snapshot": ("tessera.core.catalog", "open_snapshot"),
    "load_profile": ("tessera.core.profile", "load_profile"),
}


def __getattr__(name: str):
    # Lazy imports keep ``import tessera`` cheap for pipelines that only need one module.
    if name in _LAZY:
        import importlib

        mod, attr = _LAZY[name]
        return getattr(importlib.import_module(mod), attr)
    raise AttributeError(name)


__all__ = ["__version__", "SCIENCE_MODEL_VERSION", *_LAZY]

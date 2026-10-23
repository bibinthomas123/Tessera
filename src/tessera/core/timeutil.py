"""Time, Earth-orientation and ephemeris configuration.

The critical path must never block on the network: IERS auto-download is
disabled unless ``TESSERA_IERS_AUTO_DOWNLOAD=1``.  A local IERS-A file can be
supplied with ``TESSERA_IERS_A_PATH`` (the container refreshes it daily).  The
table used is recorded in every manifest.
"""

from __future__ import annotations

import os
import warnings
from functools import lru_cache
from typing import Any

import numpy as np
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time
from astropy.utils import iers

_CONFIGURED = False


def configure() -> dict[str, Any]:
    global _CONFIGURED
    info: dict[str, Any] = {}
    if not _CONFIGURED:
        if os.environ.get("TESSERA_IERS_AUTO_DOWNLOAD", "0") != "1":
            iers.conf.auto_download = False
            iers.conf.auto_max_age = None
            iers.conf.iers_degraded_accuracy = "warn"
        path = os.environ.get("TESSERA_IERS_A_PATH")
        if path and os.path.isfile(path):
            iers.earth_orientation_table.set(iers.IERS_A.open(path))
        eph = os.environ.get("TESSERA_EPHEMERIS", "builtin")
        solar_system_ephemeris.set(eph)
        from astropy.coordinates.baseframe import NonRotationTransformationWarning

        warnings.filterwarnings("ignore", category=NonRotationTransformationWarning)
        warnings.filterwarnings("ignore", message=".*IERS.*")
        warnings.filterwarnings("ignore", message=".*polar motion.*")
        _CONFIGURED = True
    return environment_info() | info


@lru_cache(maxsize=1)
def environment_info() -> dict[str, Any]:
    table = iers.earth_orientation_table.get()
    try:
        mjd = np.asarray(table["MJD"])
        last = Time(float(mjd.max()), format="mjd").iso[:10]
    except Exception:
        last = None
    return {
        "iers_table": type(table).__name__,
        "iers_last_date": last,
        "iers_auto_download": bool(iers.conf.auto_download),
        "iers_a_path": os.environ.get("TESSERA_IERS_A_PATH"),
        "ephemeris": str(solar_system_ephemeris.get()),
    }


def to_time(t: str | float | Time) -> Time:
    if isinstance(t, Time):
        return t
    if isinstance(t, (int, float)):
        return Time(float(t), format="mjd", scale="utc")
    return Time(str(t).replace("Z", ""), format="isot" if "T" in str(t) else None, scale="utc")


def iso(t: Time | float) -> str:
    """ISO-8601 UTC with a trailing Z, second resolution (stable for hashing)."""
    tt = t if isinstance(t, Time) else Time(t, format="jd", scale="utc")
    return tt.utc.isot[:19] + "Z"


configure()

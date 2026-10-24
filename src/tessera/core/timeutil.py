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
    """Parse a UTC time: an astropy Time, an MJD number, or an ISO-8601 string.

    Strings may use a space instead of ``T``, omit the seconds, and end in ``Z`` or ``UTC``
    (for example ``2017-08-17 22:00 UTC``).  Anything else raises ConfigError (HTTP 422).
    """
    if isinstance(t, Time):
        return t
    if isinstance(t, (int, float)):
        return Time(float(t), format="mjd", scale="utc")
    import re

    from tessera.core.errors import ConfigError

    s = str(t).strip()
    s = re.sub(r"\s*(Z|UTC|UT)$", "", s, flags=re.IGNORECASE).strip()
    m = re.fullmatch(r"(\d{4}-\d{2}-\d{2})(?:[T ](\d{1,2}:\d{2})(:\d{2}(?:\.\d+)?)?)?", s)
    if not m:
        raise ConfigError(f"invalid time {t!r}: expected UTC ISO-8601 such as 2017-08-17T22:00:00",
                          code="INVALID_TIME")
    date, hm, sec = m.groups()
    iso_s = f"{date}T{hm or '00:00'}{sec or ':00'}"
    try:
        return Time(iso_s, format="isot", scale="utc")
    except ValueError as exc:
        raise ConfigError(f"invalid time {t!r}: {exc}", code="INVALID_TIME") from exc


def iso(t: Time | float) -> str:
    """ISO-8601 UTC with a trailing Z, second resolution (stable for hashing)."""
    tt = t if isinstance(t, Time) else Time(t, format="jd", scale="utc")
    return tt.utc.isot[:19] + "Z"


configure()

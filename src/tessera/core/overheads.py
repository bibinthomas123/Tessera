"""Slew, settle, acquisition, readout and filter-change overheads (spec 5.3)."""

from __future__ import annotations

import numpy as np

from tessera.core.geometry import angsep
from tessera.core.profile import Axis, TelescopeProfile


def axis_time(delta_deg, axis: Axis) -> np.ndarray:
    """Trapezoidal velocity profile: 2 sqrt(D/a) if D < v^2/a, else D/v + v/a."""
    d = np.abs(np.asarray(delta_deg, float))
    v, a = axis.vmax_dps, axis.acc_dps2
    return np.where(d < v * v / a, 2.0 * np.sqrt(d / a), d / v + v / a)


def hadec_to_altaz(ha_deg, dec_deg, lat_deg) -> tuple[np.ndarray, np.ndarray]:
    """Geometric alt/az (no refraction) from hour angle and declination; az east of north."""
    h, d, p = (np.radians(np.asarray(x, float)) for x in (ha_deg, dec_deg, lat_deg))
    sin_alt = np.sin(d) * np.sin(p) + np.cos(d) * np.cos(p) * np.cos(h)
    alt = np.arcsin(np.clip(sin_alt, -1, 1))
    az = np.arctan2(-np.cos(d) * np.sin(h), np.sin(d) * np.cos(p) - np.cos(d) * np.sin(p) * np.cos(h))
    return np.degrees(alt), np.mod(np.degrees(az), 360.0)


def wrap180(x) -> np.ndarray:
    return (np.asarray(x, float) + 180.0) % 360.0 - 180.0


class OverheadModel:
    def __init__(self, profile: TelescopeProfile):
        self.p = profile
        self.lat = profile.site.lat_deg

    def slew_time(self, ra1, dec1, ra2, dec2, lst_deg) -> np.ndarray:
        """Seconds to move from (ra1, dec1) to (ra2, dec2) at local sidereal time ``lst_deg``.

        Also includes dome rotation (max with mount) and settle.  Acquisition,
        readout and filter changes are added by :meth:`gap`.
        """
        s = self.p.slew
        ra1, dec1, ra2, dec2 = (np.asarray(x, float) for x in (ra1, dec1, ra2, dec2))
        ha1, ha2 = wrap180(lst_deg - ra1), wrap180(lst_deg - ra2)
        alt1, az1 = hadec_to_altaz(ha1, dec1, self.lat)
        alt2, az2 = hadec_to_altaz(ha2, dec2, self.lat)
        if s.axis1 is not None:
            if self.p.mount.type == "equatorial":
                d1, d2 = np.abs(wrap180(ha2 - ha1)), np.abs(dec2 - dec1)
            else:
                d1, d2 = np.abs(wrap180(az2 - az1)), np.abs(alt2 - alt1)
            t1, t2 = axis_time(d1, s.axis1), axis_time(d2, s.axis2)
            mount = np.maximum(t1, t2) if s.axes == "simultaneous" else t1 + t2
        else:
            mount = s.constant_s + angsep(ra1, dec1, ra2, dec2) / s.rate_dps
        if s.dome_dps:
            mount = np.maximum(mount, np.abs(wrap180(az2 - az1)) / s.dome_dps)
        same = (np.abs(ra1 - ra2) < 1e-9) & (np.abs(dec1 - dec2) < 1e-9)
        return np.where(same, 0.0, mount + s.settle_s)

    def gap(self, ra1, dec1, ra2, dec2, lst_deg, filter_change) -> np.ndarray:
        """Time between the end of one visit and the start of the next (seconds)."""
        slew = self.slew_time(ra1, dec1, ra2, dec2, lst_deg)
        fc = np.where(np.asarray(filter_change, bool), self.p.instrument.filter_change_s, 0.0)
        same = slew == 0
        acq = np.where(same, 0.0, self.p.slew.acquisition_s)
        # Spec 5.3: t_gap = max(t_slew, t_dome) + t_settle + t_acq + 1[filter change] t_filter.
        # Filter changes are charged in full (not overlapped with the slew): conservative.
        return slew + acq + fc

    def visit_duration(self, n_exp: int, t_single_s: float) -> float:
        return n_exp * (t_single_s + self.p.instrument.readout_s)

    def park_radec(self, lst_deg: float) -> tuple[float, float]:
        park = self.p.mount.park
        ha = park.get("ha_h", 0.0) * 15.0
        return float(np.mod(lst_deg - ha, 360.0)), float(park.get("dec_deg", self.lat))

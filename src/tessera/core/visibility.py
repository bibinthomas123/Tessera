"""Observability windows from hard constraints (spec 5.1, 5.2), using astroplan.

1. Evaluate the Sun and operational windows on a coarse grid to find
   candidate night times (cheap: one ephemeris per grid time).
2. Evaluate every target constraint on fields x night times (astroplan,
   gridded).
3. Refine every window edge by vectorised bisection to ``refine_s``, keeping
   only times verified to pass, so windows never include unchecked time.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import astropy.units as u
import numpy as np
from astroplan import (
    AirmassConstraint,
    AltitudeConstraint,
    AtNightConstraint,
    MoonIlluminationConstraint,
    MoonSeparationConstraint,
    Observer,
    SunSeparationConstraint,
)
from astroplan.constraints import Constraint, _get_altaz
from astroplan.moon import moon_illumination
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from astropy.time import Time

from tessera.core import timeutil  # noqa: F401  (configures IERS before any transform)
from tessera.core.profile import TelescopeProfile

# ---------------------------------------------------------------------------
# Custom astroplan constraints
# ---------------------------------------------------------------------------


class HorizonMaskConstraint(Constraint):
    """Piecewise-linear minimum altitude as a function of azimuth."""

    def __init__(self, mask: list[tuple[float, float]]):
        pts = np.asarray(sorted(mask), float)
        if pts[0, 0] > 0:
            pts = np.vstack([[0.0, pts[-1, 1]], pts])
        if pts[-1, 0] < 360:
            pts = np.vstack([pts, [360.0, pts[0, 1]]])
        self.az, self.alt = pts[:, 0], pts[:, 1]

    def min_alt(self, az_deg):
        return np.interp(np.mod(az_deg, 360.0), self.az, self.alt)

    def compute_constraint(self, times, observer, targets):
        aa = _get_altaz(times, observer, targets)["altaz"]
        return aa.alt.deg >= self.min_alt(aa.az.deg)


class HourAngleConstraint(Constraint):
    def __init__(self, lo_h: float, hi_h: float):
        self.lo, self.hi = lo_h, hi_h

    def compute_constraint(self, times, observer, targets):
        lst = observer.local_sidereal_time(times).hour
        ra = targets.ra.hour
        ha = (lst - ra + 12.0) % 24.0 - 12.0
        return (ha >= self.lo) & (ha <= self.hi)


class DeclinationConstraint(Constraint):
    def __init__(self, lo: float, hi: float):
        self.lo, self.hi = lo, hi

    def compute_constraint(self, times, observer, targets):
        d = targets.dec.deg
        res = (d >= self.lo) & (d <= self.hi)
        return np.broadcast_to(res, np.broadcast(np.empty(times.shape), np.empty(targets.shape)).shape)


class MoonSeparationByIllumination(Constraint):
    """Minimum Moon separation that depends on lunar illumination (table)."""

    def __init__(self, table: list[tuple[float, float]]):
        tab = np.asarray(sorted(table), float)
        self.x, self.y = tab[:, 0], tab[:, 1]

    def compute_constraint(self, times, observer, targets):
        moon = get_body("moon", times, location=observer.location)
        sep = moon.separation(targets).deg
        illum = moon_illumination(times)
        return sep >= np.interp(illum, self.x, self.y)


class MoonIlluminationAlways(Constraint):
    """Illumination limit applied whether or not the Moon is up."""

    def __init__(self, max_illum: float):
        self.max = max_illum

    def compute_constraint(self, times, observer, targets):
        return moon_illumination(times) <= self.max


class TimeWindowsConstraint(Constraint):
    def __init__(self, windows_jd: list[tuple[float, float]]):
        self.w = np.asarray(windows_jd, float).reshape(-1, 2)

    def compute_constraint(self, times, observer, targets):
        jd = times.jd
        ok = np.zeros(jd.shape, bool)
        for a, b in self.w:
            ok |= (jd >= a) & (jd <= b)
        return ok


def make_observer(profile: TelescopeProfile) -> Observer:
    s = profile.site
    loc = EarthLocation.from_geodetic(s.lon_deg * u.deg, s.lat_deg * u.deg, s.elevation_m * u.m)
    return Observer(location=loc, name=profile.id, timezone="UTC")


def build_constraints(profile: TelescopeProfile) -> list[tuple[str, Constraint]]:
    c = profile.constraints
    out: list[tuple[str, Constraint]] = [
        ("night", AtNightConstraint(max_solar_altitude=c.sun_alt_max_deg * u.deg)),
        ("altitude", AltitudeConstraint(min=c.min_alt() * u.deg, max=c.alt_max_deg * u.deg)),
        ("airmass", AirmassConstraint(max=c.airmass_max)),
    ]
    if c.horizon_mask:
        out.append(("horizon_mask", HorizonMaskConstraint(c.horizon_mask)))
    m = profile.mount
    if m.type == "equatorial" and (m.ha_limits_h[0] > -12 or m.ha_limits_h[1] < 12):
        out.append(("hour_angle", HourAngleConstraint(*m.ha_limits_h)))
    if m.dec_limits_deg[0] > -90 or m.dec_limits_deg[1] < 90:
        out.append(("declination", DeclinationConstraint(*m.dec_limits_deg)))
    if c.sun_sep_min_deg > 0:
        out.append(("sun_separation", SunSeparationConstraint(min=c.sun_sep_min_deg * u.deg)))
    ms = c.moon_sep()
    if ms.by_illumination:
        out.append(("moon_separation", MoonSeparationByIllumination(ms.by_illumination)))
    elif ms.constant:
        out.append(("moon_separation", MoonSeparationConstraint(min=ms.constant * u.deg)))
    if c.moon_illum_max < 1.0:
        out.append(("moon_illumination", MoonIlluminationConstraint(max=c.moon_illum_max)
                    if c.moon_illum_only_when_up else MoonIlluminationAlways(c.moon_illum_max)))
    if profile.windows:
        wins = [(Time(w.start.replace("Z", "")).jd, Time(w.end.replace("Z", "")).jd) for w in profile.windows]
        out.append(("operational_window", TimeWindowsConstraint(wins)))
    return out


# ---------------------------------------------------------------------------
# Window computation
# ---------------------------------------------------------------------------


@dataclass
class Visibility:
    profile_key: str
    jd: np.ndarray                    # full time grid (JD, UTC)
    step_s: float
    windows: list[list[tuple[float, float]]]
    alt: np.ndarray                   # (n_fields, n_times) deg, NaN outside night
    az: np.ndarray
    lst_deg: np.ndarray               # (n_times,)
    sun_alt: np.ndarray
    moon_alt: np.ndarray
    moon_illum: np.ndarray
    never_reason: list[str | None] = field(default_factory=list)
    constraint_names: list[str] = field(default_factory=list)
    timings: dict[str, float] = field(default_factory=dict)

    def subset(self, idx: np.ndarray) -> Visibility:
        """Visibility for a subset of fields (shared time grid and ephemerides)."""
        idx = np.asarray(idx, int)
        return Visibility(profile_key=self.profile_key, jd=self.jd, step_s=self.step_s,
                          windows=[self.windows[i] for i in idx], alt=self.alt[idx], az=self.az[idx],
                          lst_deg=self.lst_deg, sun_alt=self.sun_alt, moon_alt=self.moon_alt,
                          moon_illum=self.moon_illum, never_reason=[self.never_reason[i] for i in idx],
                          constraint_names=self.constraint_names, timings=self.timings)

    def airmass(self, i_field: int, jd: float) -> float:
        k = int(np.clip(np.round((jd - self.jd[0]) * 86400.0 / self.step_s), 0, len(self.jd) - 1))
        a = self.alt[i_field, k]
        if not np.isfinite(a) or a <= 0:
            return 10.0
        return float(1.0 / np.sin(np.radians(a)))

    def time_index(self, jd: float) -> int:
        return int(np.clip(np.round((jd - self.jd[0]) * 86400.0 / self.step_s), 0, len(self.jd) - 1))

    def lst_at(self, jd: float) -> float:
        # LST advances 360.9856 deg per solar day; interpolate from the nearest grid node
        k = self.time_index(jd)
        return float((self.lst_deg[k] + 360.98564736629 * (jd - self.jd[k])) % 360.0)

    def window_containing(self, i_field: int, jd0: float, jd1: float) -> tuple[float, float] | None:
        for a, b in self.windows[i_field]:
            if a <= jd0 + 1e-9 and jd1 <= b + 1e-9:
                return (a, b)
        return None

    def next_open(self, i_field: int, jd: float) -> float | None:
        for a, b in self.windows[i_field]:
            if b > jd:
                return max(a, jd)
        return None


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Inclusive (start, end) index runs of True."""
    if not mask.any():
        return []
    d = np.diff(np.concatenate([[0], mask.view(np.int8), [0]]))
    starts = np.flatnonzero(d == 1)
    ends = np.flatnonzero(d == -1) - 1
    return list(zip(starts.tolist(), ends.tolist()))


def compute_visibility(
    profile: TelescopeProfile,
    ra: np.ndarray,
    dec: np.ndarray,
    t_start: Time,
    t_end: Time,
    *,
    step_s: float = 300.0,
    refine_s: float = 10.0,
) -> Visibility:
    import time as _time

    tm: dict[str, float] = {}
    t0 = _time.perf_counter()
    observer = make_observer(profile)
    constraints = build_constraints(profile)
    n = len(ra)
    nstep = int(np.ceil((t_end - t_start).to_value(u.s) / step_s)) + 1
    jd = t_start.jd + np.arange(nstep) * step_s / 86400.0
    times = Time(jd, format="jd", scale="utc")
    frame = AltAz(obstime=times, location=observer.location)
    sun_alt = get_sun(times).transform_to(frame).alt.deg
    moon = get_body("moon", times, location=observer.location)
    moon_alt = moon.transform_to(frame).alt.deg
    moon_illum = np.asarray(moon_illumination(times), float)
    lst = observer.local_sidereal_time(times).deg
    # time-only constraints first: the Sun and operational windows
    time_ok = sun_alt <= profile.constraints.sun_alt_max_deg + 0.5  # +0.5 deg margin; exact test below
    for _, c in constraints:
        if isinstance(c, TimeWindowsConstraint):
            time_ok &= c.compute_constraint(times, observer, None)
    tm["ephemerides"] = _time.perf_counter() - t0

    alt = np.full((n, nstep), np.nan, dtype=np.float32)
    az = np.full((n, nstep), np.nan, dtype=np.float32)
    ok = np.zeros((n, nstep), bool)
    idx_t = np.flatnonzero(time_ok)
    targets = SkyCoord(ra=np.asarray(ra) * u.deg, dec=np.asarray(dec) * u.deg)
    t1 = _time.perf_counter()
    if n and len(idx_t):
        sub = times[idx_t]
        grid_ok = np.ones((n, len(idx_t)), bool)
        per_constraint = {}
        for name, c in constraints:
            r = np.asarray(c(observer, targets, sub, grid_times_targets=True), bool)
            per_constraint[name] = r
            grid_ok &= r
        tg = targets[:, np.newaxis]
        aa = _get_altaz(sub, observer, tg)["altaz"]
        alt[:, idx_t] = aa.alt.deg
        az[:, idx_t] = aa.az.deg
        ok[:, idx_t] = grid_ok
    else:
        per_constraint = {}
    tm["grid"] = _time.perf_counter() - t1

    # refine edges by vectorised bisection
    t2 = _time.perf_counter()
    edges: list[tuple[int, int, int]] = []        # (field, lo_idx(false side), hi_idx(true side)) rising/falling
    runs_per_field = [_runs(ok[i]) for i in range(n)]
    rise_f, rise_lo, rise_hi, fall_f, fall_lo, fall_hi = [], [], [], [], [], []
    for i, runs in enumerate(runs_per_field):
        for s, e in runs:
            if s > 0:
                rise_f.append(i)
                rise_lo.append(jd[s - 1])
                rise_hi.append(jd[s])
            if e < nstep - 1:
                fall_f.append(i)
                fall_lo.append(jd[e])
                fall_hi.append(jd[e + 1])
    del edges
    n_iter = max(0, int(np.ceil(np.log2(step_s / refine_s))))

    def evaluate(fidx: np.ndarray, tjd: np.ndarray) -> np.ndarray:
        if len(fidx) == 0:
            return np.zeros(0, bool)
        tt = Time(tjd, format="jd", scale="utc")
        tg = targets[fidx]
        res = np.ones(len(fidx), bool)
        for _, c in constraints:
            res &= np.asarray(c(observer, tg, tt, grid_times_targets=False), bool)
        return res

    rf, rlo, rhi = np.asarray(rise_f, int), np.asarray(rise_lo), np.asarray(rise_hi)   # lo false, hi true
    ff, flo, fhi = np.asarray(fall_f, int), np.asarray(fall_lo), np.asarray(fall_hi)   # lo true, hi false
    for _ in range(n_iter):
        if len(rf):
            mid = 0.5 * (rlo + rhi)
            r = evaluate(rf, mid)
            rhi = np.where(r, mid, rhi)
            rlo = np.where(r, rlo, mid)
        if len(ff):
            mid = 0.5 * (flo + fhi)
            r = evaluate(ff, mid)
            flo = np.where(r, mid, flo)
            fhi = np.where(r, fhi, mid)
    tm["refine"] = _time.perf_counter() - t2

    rise_map = {(int(f), float(h0)): float(h) for f, h0, h in zip(rf, np.asarray(rise_hi), rhi)}
    fall_map = {(int(f), float(l0)): float(lv) for f, l0, lv in zip(ff, np.asarray(fall_lo), flo)}
    windows: list[list[tuple[float, float]]] = []
    for i, runs in enumerate(runs_per_field):
        w = []
        for s, e in runs:
            a = rise_map.get((i, float(jd[s])), float(jd[s]))
            b = fall_map.get((i, float(jd[e])), float(jd[e]))
            if b > a:
                w.append((a, b))
        windows.append(w)

    never: list[str | None] = []
    for i in range(n):
        if windows[i]:
            never.append(None)
            continue
        if not len(idx_t):
            never.append("no_night_time")
            continue
        failing = [name for name, r in per_constraint.items() if not r[i].any()]
        never.append(failing[0] if failing else "no_common_window")
    return Visibility(profile_key=profile.key, jd=jd, step_s=step_s, windows=windows, alt=alt, az=az, lst_deg=lst,
                      sun_alt=sun_alt, moon_alt=moon_alt, moon_illum=moon_illum, never_reason=never,
                      constraint_names=[nm for nm, _ in constraints], timings=tm)


def describe_constraints(profile: TelescopeProfile) -> list[dict[str, Any]]:
    return [{"name": n, "type": type(c).__name__} for n, c in build_constraints(profile)]

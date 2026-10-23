"""Independent feasibility validator (spec 5.4).

Deliberately a separate code path from :mod:`tessera.core.visibility`: it uses
astropy directly (not astroplan, not the scheduler's cached windows) and
recomputes altitude, airmass, Sun and Moon at ``step_s`` resolution across
every visit.  Any violation marks the visit (and the schedule) infeasible;
there are no soft penalties.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
from astropy.time import Time

from tessera.core import timeutil  # noqa: F401
from tessera.core.overheads import OverheadModel
from tessera.core.profile import TelescopeProfile

TOL_S = 1.0          # timing tolerance (s)
TOL_DEG = 1e-3       # angular tolerance (deg) for floating noise at window edges


def _illumination(sun, moon) -> np.ndarray:
    """Fraction of the lunar disk illuminated, from the Sun-Moon-observer phase angle."""
    elong = sun.separation(moon).rad
    r_sun = sun.distance.to_value(u.km)
    r_moon = moon.distance.to_value(u.km)
    phase = np.arctan2(r_sun * np.sin(elong), r_moon - r_sun * np.cos(elong))
    return (1.0 + np.cos(phase)) / 2.0


def _to_jd(v: dict[str, Any], key: str) -> float:
    if f"{key}_jd" in v:
        return float(v[f"{key}_jd"])
    return Time(str(v[key]).replace("Z", ""), scale="utc").jd


def validate_visits(visits: list[dict[str, Any]], profile: TelescopeProfile, *, step_s: float = 10.0,
                    horizon: tuple[float, float] | None = None) -> dict[str, Any]:
    c = profile.constraints
    loc = EarthLocation.from_geodetic(profile.site.lon_deg * u.deg, profile.site.lat_deg * u.deg,
                                      profile.site.elevation_m * u.m)
    ovh = OverheadModel(profile)
    min_alt = c.min_alt()
    mask = None
    if c.horizon_mask:
        pts = np.asarray(sorted(c.horizon_mask), float)
        if pts[0, 0] > 0:
            pts = np.vstack([[0.0, pts[-1, 1]], pts])
        if pts[-1, 0] < 360:
            pts = np.vstack([pts, [360.0, pts[0, 1]]])
        mask = pts
    moon_sep = c.moon_sep()
    wins = None
    if profile.windows:
        wins = [(Time(w.start.replace("Z", "")).jd, Time(w.end.replace("Z", "")).jd) for w in profile.windows]

    # sample every visit at step_s (always including both ends)
    samples_jd, owner = [], []
    for k, v in enumerate(visits):
        a, b = _to_jd(v, "start"), _to_jd(v, "end")
        n = max(2, int(np.ceil((b - a) * 86400.0 / step_s)) + 1)
        samples_jd.append(np.linspace(a, b, n))
        owner.append(np.full(n, k))
    reports: list[dict[str, Any]] = []
    if not visits:
        return {"feasible": True, "n_visits": 0, "n_infeasible": 0, "violations_by_constraint": {}, "step_s": step_s,
                "visits": []}
    jd = np.concatenate(samples_jd)
    own = np.concatenate(owner)
    times = Time(jd, format="jd", scale="utc")
    ra = np.array([visits[k]["ra"] for k in own], float)
    dec = np.array([visits[k]["dec"] for k in own], float)
    frame = AltAz(obstime=times, location=loc)
    tgt = SkyCoord(ra=ra * u.deg, dec=dec * u.deg).transform_to(frame)
    sun = get_body("sun", times, location=loc)
    moon = get_body("moon", times, location=loc)
    sun_alt = sun.transform_to(frame).alt.deg
    moon_alt = moon.transform_to(frame).alt.deg
    alt = tgt.alt.deg
    az = tgt.az.deg
    with np.errstate(divide="ignore"):
        airmass = np.where(alt > 0, 1.0 / np.sin(np.radians(alt)), np.inf)
    sky = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    sun_sep = sun.separation(sky).deg
    m_sep = moon.separation(sky).deg
    illum = _illumination(sun, moon)
    lst = times.sidereal_time("apparent", longitude=profile.site.lon_deg * u.deg).deg
    ha_h = ((lst - ra) / 15.0 + 12.0) % 24.0 - 12.0

    checks: dict[str, np.ndarray] = {
        "night": sun_alt <= c.sun_alt_max_deg + TOL_DEG,
        "altitude": (alt >= min_alt - TOL_DEG) & (alt <= c.alt_max_deg + TOL_DEG),
        "airmass": airmass <= c.airmass_max + 1e-4,
    }
    if mask is not None:
        checks["horizon_mask"] = alt >= np.interp(np.mod(az, 360), mask[:, 0], mask[:, 1]) - TOL_DEG
    m = profile.mount
    if m.type == "equatorial":
        checks["hour_angle"] = (ha_h >= m.ha_limits_h[0] - 1e-4) & (ha_h <= m.ha_limits_h[1] + 1e-4)
    checks["declination"] = (dec >= m.dec_limits_deg[0]) & (dec <= m.dec_limits_deg[1])
    if c.sun_sep_min_deg > 0:
        checks["sun_separation"] = sun_sep >= c.sun_sep_min_deg - TOL_DEG
    checks["moon_separation"] = m_sep >= moon_sep.min_sep(illum) - 0.05
    if c.moon_illum_max < 1.0:
        lim = illum <= c.moon_illum_max + 1e-3
        checks["moon_illumination"] = (lim | (moon_alt < 0)) if c.moon_illum_only_when_up else lim
    if wins:
        okw = np.zeros(len(jd), bool)
        for a, b in wins:
            okw |= (jd >= a - TOL_S / 86400) & (jd <= b + TOL_S / 86400)
        checks["operational_window"] = okw

    by_constraint: dict[str, int] = {}
    for k, v in enumerate(visits):
        sel = own == k
        viol: list[dict[str, Any]] = []
        for name, okarr in checks.items():
            bad = ~okarr[sel]
            if bad.any():
                j = np.flatnonzero(sel)[np.flatnonzero(bad)[0]]
                worst = {"night": sun_alt, "altitude": alt, "airmass": airmass, "sun_separation": sun_sep,
                         "moon_separation": m_sep, "hour_angle": ha_h}.get(name)
                viol.append({"constraint": name, "at": timeutil.iso(float(jd[j])),
                             "value": None if worst is None else round(float(worst[j]), 4)})
        if v.get("filter") not in profile.instrument.filters:
            viol.append({"constraint": "filter_unavailable", "value": v.get("filter")})
        a, b = _to_jd(v, "start"), _to_jd(v, "end")
        dur_need = ovh.visit_duration(int(v.get("n_exp", 1)), float(v.get("t_exp_s", 0.0)))
        if (b - a) * 86400.0 + TOL_S < dur_need:
            viol.append({"constraint": "duration_shorter_than_exposures", "value": round((b - a) * 86400.0, 1),
                         "limit": round(dur_need, 1)})
        if horizon and (a < horizon[0] - TOL_S / 86400 or b > horizon[1] + TOL_S / 86400):
            viol.append({"constraint": "outside_horizon"})
        if k > 0:
            pv = visits[k - 1]
            pb = _to_jd(pv, "end")
            if a < pb - TOL_S / 86400:
                viol.append({"constraint": "overlap", "with_seq": pv.get("seq")})
            else:
                lst_b = float(Time(pb, format="jd").sidereal_time("apparent",
                                                                  longitude=profile.site.lon_deg * u.deg).deg)
                need = float(ovh.gap(pv["ra"], pv["dec"], v["ra"], v["dec"], lst_b, pv.get("filter") != v.get("filter")))
                if (a - pb) * 86400.0 + TOL_S < need:
                    viol.append({"constraint": "gap_shorter_than_overheads", "value": round((a - pb) * 86400, 1),
                                 "limit": round(need, 1)})
        for vv in viol:
            by_constraint[vv["constraint"]] = by_constraint.get(vv["constraint"], 0) + 1
        am = airmass[sel]
        reports.append({
            "seq": v.get("seq", k + 1),
            "feasible": not viol,
            "violations": viol,
            "airmass": {"start": round(float(am[0]), 3), "end": round(float(am[-1]), 3),
                        "max": round(float(np.max(am)), 3)},
            "moon_sep_min_deg": round(float(np.min(m_sep[sel])), 2),
            "moon_illumination": round(float(illum[sel][0]), 3),
            "sun_alt_max_deg": round(float(np.max(sun_alt[sel])), 2),
            "alt_min_deg": round(float(np.min(alt[sel])), 2),
        })
    # revisit gaps
    pattern = profile.visit_pattern
    seen: dict[str, list[int]] = {}
    for k, v in enumerate(visits):
        seen.setdefault(v["field_id"], []).append(k)
    for ks in seen.values():
        for j in range(1, len(ks)):
            step = visits[ks[j]].get("pattern_step", j)
            if step >= len(pattern):
                continue
            sd = pattern[step]
            gap_min = (_to_jd(visits[ks[j]], "start") - _to_jd(visits[ks[j - 1]], "end")) * 1440.0
            bad = gap_min + 1e-3 < sd.min_gap_min or (sd.max_gap_min is not None and gap_min > sd.max_gap_min + 1e-3)
            if bad:
                reports[ks[j]]["violations"].append({"constraint": "revisit_gap", "value": round(gap_min, 2),
                                                     "limit": [sd.min_gap_min, sd.max_gap_min]})
                reports[ks[j]]["feasible"] = False
                by_constraint["revisit_gap"] = by_constraint.get("revisit_gap", 0) + 1
    n_bad = sum(1 for r in reports if not r["feasible"])
    return {"feasible": n_bad == 0, "n_visits": len(visits), "n_infeasible": n_bad,
            "violations_by_constraint": by_constraint, "step_s": step_s, "visits": reports}

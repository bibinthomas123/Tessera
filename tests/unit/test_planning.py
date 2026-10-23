"""Footprints, overheads, visibility, scheduling and validation (spec sections 5 and 6)."""

from __future__ import annotations

import numpy as np
import pytest

from tessera.core.footprint import InstrumentFootprint
from tessera.core.geometry import angsep, inverse_gnomonic
from tessera.core.overheads import OverheadModel, axis_time, hadec_to_altaz
from tessera.core.profile import Axis, Footprint, load_profile
from tessera.core.validate import validate_visits

from ..conftest import T_GW170817

# -- footprints ---------------------------------------------------------------

def test_rect_containment_and_rotation():
    fp = InstrumentFootprint(Footprint(kind="rect", width_deg=2.0, height_deg=0.5))
    # 0.8 deg east of centre: inside at PA 0 (half-width 1), outside at PA 90 (half-width 0.25)
    ra, dec = inverse_gnomonic(100.0, 0.0, np.array([0.8]), np.array([0.0]))
    assert fp.contains(100.0, 0.0, 0.0, ra, dec)[0]
    assert not fp.contains(100.0, 0.0, 90.0, ra, dec)[0]
    assert fp.area_deg2 == pytest.approx(1.0)


def test_chip_gap_not_covered():
    chips = [[(-1.0, -0.5), (-0.05, -0.5), (-0.05, 0.5), (-1.0, 0.5)],
             [(0.05, -0.5), (1.0, -0.5), (1.0, 0.5), (0.05, 0.5)]]
    fp = InstrumentFootprint(Footprint(kind="chips", chips=chips))
    ra, dec = inverse_gnomonic(50.0, 20.0, np.array([0.0, 0.5, -0.5]), np.array([0.0, 0.0, 0.0]))
    assert list(fp.contains(50.0, 20.0, 0.0, ra, dec)) == [False, True, True]


def test_rotation_averaged_coverage_is_fractional():
    fp = InstrumentFootprint(Footprint(kind="rect", width_deg=2.0, height_deg=0.5, rotation_averaged=True))
    ra, dec = inverse_gnomonic(10.0, 0.0, np.array([0.0, 0.8]), np.array([0.0, 0.0]))
    w = fp.coverage(10.0, 0.0, 0.0, ra, dec)
    assert w[0] == 1.0 and 0.0 < w[1] < 1.0


def test_circle_raster_matches_area():
    fp = InstrumentFootprint(Footprint(kind="circle", radius_deg=1.75))
    r = fp.rasterize(30.0, -30.0, 0.0, 9)
    pix_area = 41252.96 / (12 * 4**9)
    assert len(r.ipix) * pix_area == pytest.approx(fp.area_deg2, rel=0.03)


def test_sky_polygon_roundtrip():
    fp = InstrumentFootprint(Footprint(kind="rect", width_deg=1.0, height_deg=1.0))
    (poly,) = fp.sky_polygons(200.0, -45.0, 30.0)
    ra, dec = np.array(poly).T
    assert np.all(angsep(200.0, -45.0, ra, dec) <= 0.71 + 1e-6)


# -- overheads -----------------------------------------------------------------

def test_axis_time_regimes():
    ax = Axis(vmax_dps=2.0, acc_dps2=1.0)
    assert float(axis_time(1.0, ax)) == pytest.approx(2.0)          # triangular: 2 sqrt(1/1)
    assert float(axis_time(10.0, ax)) == pytest.approx(10 / 2 + 2)  # D/v + v/a


def test_hadec_to_altaz_against_astropy():
    import astropy.units as u
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord
    from astropy.time import Time

    loc = EarthLocation.from_geodetic(-70.69 * u.deg, -29.01 * u.deg, 0 * u.m)
    t = Time("2017-08-17T23:30:00", scale="utc")
    lst = t.sidereal_time("apparent", longitude=loc.lon).deg
    c = SkyCoord(197.45 * u.deg, -23.38 * u.deg).transform_to(AltAz(obstime=t, location=loc))
    # compare with the target's apparent HA/Dec (precession makes ICRS differ by ~0.1 deg)
    from astropy.coordinates import TETE

    app = SkyCoord(197.45 * u.deg, -23.38 * u.deg).transform_to(TETE(obstime=t))
    alt, az = hadec_to_altaz(lst - app.ra.deg, app.dec.deg, -29.01)
    assert float(alt) == pytest.approx(c.alt.deg, abs=0.05)
    assert float(az) == pytest.approx(c.az.deg, abs=0.1)


def test_gap_charges_filter_change_and_acquisition():
    p = load_profile("example-1m")
    o = OverheadModel(p)
    same = float(o.gap(10, -20, 10, -20, 0.0, False))
    assert same == 0.0
    a = float(o.gap(10, -20, 11, -20, 0.0, False))
    b = float(o.gap(10, -20, 11, -20, 0.0, True))
    assert b - a == pytest.approx(p.instrument.filter_change_s)
    assert a >= p.slew.settle_s + p.slew.acquisition_s


# -- visibility, cross-checked with Skyfield ------------------------------------

def test_windows_satisfy_constraints_per_skyfield():
    """Window edges and midpoints are above the altitude limit according to Skyfield (independent library)."""
    pytest.importorskip("skyfield")
    from astropy.time import Time
    from skyfield.api import load, wgs84

    from tessera.core.geometry import radec_to_xyz
    from tessera.core.visibility import compute_visibility

    p = load_profile("example-1m")
    ra, dec = np.array([197.45, 210.0, 300.0]), np.array([-23.38, -40.0, -60.0])
    t0 = Time(T_GW170817, scale="utc")
    vis = compute_visibility(p, ra, dec, t0, t0 + 0.6, step_s=300, refine_s=10)
    ts = load.timescale(builtin=True)
    site = wgs84.latlon(p.site.lat_deg, p.site.lon_deg, elevation_m=p.site.elevation_m)
    checked = 0
    for i, wins in enumerate(vis.windows):
        v = radec_to_xyz(ra[i], dec[i])
        for a, b in wins:
            for jd in (a, 0.5 * (a + b), b):
                t = ts.from_astropy(Time(jd, format="jd", scale="utc"))
                alt = np.degrees(np.arcsin((site.rotation_at(t) @ v)[2]))
                # 0.02 deg covers annual aberration, which this geometric check omits
                assert alt >= p.constraints.min_alt() - 0.02, (i, jd, alt)
                checked += 1
    assert checked >= 3


# -- scheduling ------------------------------------------------------------------

def test_gw170817_schedule_is_feasible_and_ordered(gw170817_plan):
    s = gw170817_plan.schedules["example-1m@1"]
    assert s["status"] == "proposed" and s["feasible"] and s["validation"]["n_infeasible"] == 0
    v = s["visits"]
    assert len(v) > 10
    starts = [x["start"] for x in v]
    assert starts == sorted(starts)
    for a, b in zip(v, v[1:]):
        assert a["end"] <= b["start"]
    for x in v:
        assert x["window"][0] <= x["start"] and x["end"] <= x["window"][1]
        assert x["decision"]["rule"] and "cum" in x
    cums = [x["cum"]["p_sky"] for x in v]
    assert all(b >= a - 1e-12 for a, b in zip(cums, cums[1:]))


def test_three_coverage_metrics_reported_separately(gw170817_plan):
    sm = gw170817_plan.schedules["example-1m@1"]["summary"]
    for k in ("cum_p_sky", "cum_score_fraction_in_field", "cum_expected_detections", "expected_detections_per_hour",
              "p_uncat"):
        assert k in sm


def test_empty_schedule_explains_itself(glade_snapshot):
    """A target that never rises above the airmass limit from a northern site gets a reason."""
    from tessera import PlanConfig, plan
    from tessera.testing import gaussian_map

    m = gaussian_map([(0.0, -80.0, 1.0, 1.0)], order=8)
    r = plan(m, glade_snapshot, ["example-wide-north"], PlanConfig(preset="bbh"), "2017-08-17T22:00:00")
    s = r.schedules["example-wide-north@1"]
    assert s["visits"] == [] and s["empty_reason"] == "NO_FIELD_VISIBLE_IN_HORIZON"
    assert "EMPTY_SCHEDULE" in [d["code"] for d in r.diagnostics]
    assert any(k.startswith("REJECT") for k in s["unscheduled_counts"])


def test_revisit_pattern_respected(glade_snapshot):
    from tessera import PlanConfig, plan
    from tessera.testing import gaussian_map

    m = gaussian_map([(197.45, -23.38, 3.0, 1.0)], order=8, dist=(150, 40))
    r = plan(m, glade_snapshot, ["rubin-lsstcam-approx"], PlanConfig(preset="bns", horizon_nights=1,
                                                                      compare_strategies=False), T_GW170817)
    s = r.schedules["rubin-lsstcam-approx@1"]
    assert s["feasible"]
    by_field: dict[str, list] = {}
    for v in s["visits"]:
        by_field.setdefault(v["field_id"], []).append(v)
    pairs = [vs for vs in by_field.values() if len(vs) == 2]
    assert pairs
    from tessera.core.timeutil import to_time

    for a, b in pairs:
        assert (a["filter"], b["filter"]) == ("g", "i")
        gap_min = (to_time(b["start"]).jd - to_time(a["end"]).jd) * 1440
        assert 30 - 1e-6 <= gap_min <= 240 + 1e-6


def test_sequential_coordination_reduces_duplicate_coverage(glade_snapshot):
    from tessera import PlanConfig, plan
    from tessera.testing import gaussian_map

    m = gaussian_map([(197.45, -23.38, 2.0, 1.0)], order=8)
    cfg = dict(preset="bns", horizon_nights=1, compare_strategies=False, field_mode="tessellate")
    ind = plan(m, glade_snapshot, ["rubin-lsstcam-approx", "example-1m"], PlanConfig(**cfg), T_GW170817)
    seq = plan(m, glade_snapshot, ["rubin-lsstcam-approx", "example-1m"],
               PlanConfig(**cfg, coordination="sequential"), T_GW170817)
    a = ind.schedules["example-1m@1"]["summary"]["cum_p_sky"]
    b = seq.schedules["example-1m@1"]["summary"]["cum_p_sky"]
    assert b >= a  # the second facility's union coverage includes the first facility's fields
    assert seq.plan_key != ind.plan_key


def test_prior_observations_are_credited(gw170817_map, glade_snapshot, gw170817_plan):
    from tessera import PlanConfig, plan

    first = gw170817_plan.schedules["example-1m@1"]["visits"][:5]
    obs = [{"ra": v["ra"], "dec": v["dec"], "filter": v["filter"], "start": v["start"], "status": "completed",
            "profile_id": "example-1m"} for v in first]
    r = plan(gw170817_map, glade_snapshot, ["example-1m"], PlanConfig(preset="bns", horizon_nights=1), T_GW170817,
             prior_observations=obs)
    ids = {v["field_id"] for v in r.schedules["example-1m@1"]["visits"]}
    assert not ids & {v["field_id"] for v in first}
    assert "PRIOR_OBSERVATIONS_CREDITED" in [d["code"] for d in r.diagnostics]


# -- validator -----------------------------------------------------------------

def test_validator_flags_daylight_overlap_and_gaps(gw170817_plan):
    p = load_profile("example-1m")
    v = [dict(x) for x in gw170817_plan.schedules["example-1m@1"]["visits"][:3]]
    assert validate_visits(v, p)["feasible"]
    day = dict(v[0], start="2017-08-17T18:00:00Z", end="2017-08-17T18:01:00Z")
    rep = validate_visits([day], p)
    assert not rep["feasible"] and "night" in rep["violations_by_constraint"]
    ov = [v[0], dict(v[1], start=v[0]["start"], end=v[0]["end"])]
    rep = validate_visits(ov, p)
    assert "overlap" in rep["violations_by_constraint"]
    squeeze = [v[0], dict(v[1], start=v[0]["end"])]
    assert "gap_shorter_than_overheads" in validate_visits(squeeze, p)["violations_by_constraint"]
    bad_filter = [dict(v[0], filter="u")]
    assert "filter_unavailable" in validate_visits(bad_filter, p)["violations_by_constraint"]

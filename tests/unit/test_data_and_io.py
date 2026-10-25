"""Catalog snapshots, budget, adapters, alerts, manifests and the CLI."""

from __future__ import annotations

import json

import numpy as np
import pytest

from tessera.core.catalog import build_snapshot, find_duplicates, open_snapshot, pixels_to_index_ranges
from tessera.core.errors import AdapterError, CatalogError, ProfileError
from tessera.core.profile import load_profile, parse_profile

from ..conftest import DATA, GLADE_EXCERPT, GW170817_MAP, NGC4993, T_GW170817

# -- catalog -------------------------------------------------------------------

def test_glade_flags_mapped(glade_cols):
    ids = list(glade_cols["gal_id"])
    i = ids.index(NGC4993)
    assert glade_cols["dist_kind"][i] == "specz"           # f_dL = 3
    assert glade_cols["dist_err_includes_pv"][i]           # f_zcmb = 1
    assert "PGC:45657" in glade_cols["alt_ids"][i]
    assert set(glade_cols["dist_kind"]) <= {"none", "photoz", "z_indep", "specz"}
    assert 1e10 < 10 ** glade_cols["log_mstar"][i] < 1e11  # M* = 2.2 x 1e10


def test_snapshot_immutable_and_checksummed(glade_snapshot):
    snap = open_snapshot(glade_snapshot.path, verify=True)
    assert snap.id == glade_snapshot.id
    with pytest.raises(CatalogError):
        build_snapshot({"gal_id": np.array(["x"], dtype=object), "ra": np.zeros(1), "dec": np.zeros(1)},
                       glade_snapshot.path, name="x", upstream_version="1")


def test_region_extraction_equals_bruteforce(gw170817_map, glade_snapshot, glade_cols):
    m = gw170817_map
    reg = m.region_mask_inclusive(0.99)
    got = set(glade_snapshot.extract_for_map(m.order[reg], m.ipix[reg])["gal_id"])
    idx = m.lookup(glade_cols["ra"], glade_cols["dec"])
    want = set(np.asarray(glade_cols["gal_id"])[(idx >= 0) & reg[np.clip(idx, 0, None)]])
    assert got == want and len(got) > 50


def test_index_ranges_merge():
    # order-9 pixel 0 -> [0, 4); order-10 4, 5 -> [4, 6); order-11 24 -> parent 6: contiguous, one range
    r = pixels_to_index_ranges(np.array([9, 10, 10, 11]), np.array([0, 4, 5, 24]))
    assert r.tolist() == [[0, 7]]
    r = pixels_to_index_ranges(np.array([10, 10]), np.array([9, 3]))
    assert r.tolist() == [[3, 4], [9, 10]]


def test_duplicates_resolved_deterministically():
    cols = {"gal_id": np.array(["A", "B", "C"], dtype=object), "ra": np.array([10.0, 10.0 + 1 / 3600 / 2, 50.0]),
            "dec": np.zeros(3), "dist_mpc": np.array([40.0, 41.0, 40.0]), "dist_err_mpc": np.array([3.0, 1.0, 1.0]),
            "dist_kind": np.array(["specz", "z_indep", "specz"], dtype=object), "mag_K": np.full(3, 10.0)}
    dup = find_duplicates(cols)
    assert list(dup) == ["B", "", ""]   # B is redshift-independent, so it represents the pair


def test_generic_adapter_column_map(tmp_path):
    from tessera.adapters.catalogs import generic

    csv = tmp_path / "cat.csv"
    csv.write_text("TARGETID,RA,DEC,Z\n1,10.0,-5.0,0.01\n2,11.0,-6.0,\n")
    cols = generic.read(csv, {"columns": {"gal_id": "TARGETID", "ra": "RA", "dec": "DEC", "z": "Z"},
                              "id_prefix": "DESI:", "compute_distance_from_z": True})
    assert list(cols["gal_id"]) == ["DESI:1", "DESI:2"]
    assert 40 < cols["dist_mpc"][0] < 50 and np.isnan(cols["dist_mpc"][1])
    assert list(cols["dist_kind"]) == ["specz", "none"]


# -- budget --------------------------------------------------------------------

def test_budget_terms_are_consistent(gw170817_plan):
    b = gw170817_plan.budget
    assert b["p_cat"] + b["p_uncat"] == pytest.approx(b["p_in_footprint"], rel=1e-9)
    assert b["p_in_footprint"] + b["p_out_footprint"] + b["p_unassessed"] == pytest.approx(1.0, abs=1e-9)
    assert 0 <= b["p_cat"] <= 1 and b["by_distance"]


def test_gw170817_budget_matches_published_glade_completeness(gw170817_plan):
    """GLADE+ is complete to d_L = 47 Mpc in B-band luminosity (Dalya et al. 2022); GW170817 is ~40 Mpc."""
    assert gw170817_plan.budget["p_cat"] > 0.9


def test_budget_tracks_known_incompleteness(tmp_path, rng):
    """A halved synthetic catalog has roughly half the completeness of the full one."""
    from tessera.core.association import score_galaxies
    from tessera.core.budget import compute_budget
    from tessera.core.config import PlanConfig
    from tessera.testing import gaussian_map, synthetic_universe

    uni = synthetic_universe(rng, ra0=150, dec0=0, radius_deg=15, d_max=120, n_scale=1.0)
    keep = rng.uniform(size=len(uni.L)) < 0.5
    full = build_snapshot(uni.cols, tmp_path / "full", name="f", upstream_version="1")
    half = build_snapshot({k: v[keep] for k, v in uni.cols.items()}, tmp_path / "half", name="h", upstream_version="1")
    m = gaussian_map([(150, 0, 3.0, 1.0)], order=7, dist=(70, 12))
    cfg = PlanConfig().resolved()
    p = {}
    for name, snap in (("full", full), ("half", half)):
        reg = m.region_mask_inclusive(cfg.cr_max)
        a = score_galaxies(m, snap.extract_for_map(m.order[reg], m.ipix[reg]), cfg)
        p[name] = compute_budget(m, snap, cfg, assoc_summary=a.summary)["p_cat"]
    assert p["full"] > 0.6
    assert p["half"] / p["full"] == pytest.approx(0.5, abs=0.12)


def test_budget_without_distance_is_explicit(glade_snapshot):
    from tessera.core.budget import compute_budget
    from tessera.core.config import PlanConfig
    from tessera.core.errors import Diagnostics
    from tessera.testing import gaussian_map

    d = Diagnostics()
    b = compute_budget(gaussian_map([(197, -23, 2, 1)], order=7, dist=None), glade_snapshot,
                       PlanConfig().resolved(), diagnostics=d)
    assert b["p_cat"] is None and b["reason"] and "BUDGET_INCOMPLETE" in d.codes()


# -- profiles ------------------------------------------------------------------

def test_packaged_profiles_load():
    for pid in ("example-1m", "example-wide-north", "rubin-lsstcam-approx"):
        p = load_profile(pid)
        assert p.key.startswith(pid) and len(p.sha256) == 64


@pytest.mark.parametrize("mutate,msg", [
    (lambda d: d["site"].update(lat_deg=120), "lat"),
    (lambda d: d["instrument"]["footprint"].update(chips=[[[0, 0], [1, 1], [1, 0], [0, 1]]]), "convex"),
    (lambda d: d.update(visit_pattern=[{"filter": "z"}]), "filter"),
    (lambda d: d["instrument"]["filters"]["r"].update(m_ref=40), "m_ref"),
])
def test_profile_sanity_checks(mutate, msg):
    d = load_profile("example-1m").canonical()
    mutate(d)
    with pytest.raises(ProfileError) as e:
        parse_profile(d)
    assert msg in json.dumps(e.value.details)


# -- adapters ------------------------------------------------------------------

def test_adapters_export_and_validate(gw170817_plan, gw170817_map):
    from tessera.adapters.schedulers import available, get_adapter

    assert {"generic", "acp", "lco", "treasuremap", "rubin_too"} <= set(available())
    doc = gw170817_plan.schedules["example-1m@1"]
    p = load_profile("example-1m")
    for name, kw in (("generic", {}), ("acp", {}), ("lco", {}), ("treasuremap", {"graceid": "S170817a"}),
                     ("rubin_too", {"skymap": gw170817_map, "galaxies": gw170817_plan.galaxies,
                                    "weighting": "galaxy_boosted"})):
        ad = get_adapter(name)
        exp = ad.export(doc, p, **kw)
        assert not [i for i in ad.validate(exp) if i.level == "error"], name
        data, _, fname = exp.as_bytes()
        assert data and fname
    acp = get_adapter("acp").export(doc, p).files["tessera_plan.txt"].decode()
    assert acp.count("#waituntil") == len(doc["visits"]) and "#quitat" in acp
    csv = get_adapter("generic").export(doc, p, fmt="csv").files["schedule.csv"].decode()
    assert csv.splitlines()[1].startswith("seq,field_id")


def test_adapters_refuse_infeasible_and_unapproved(gw170817_plan):
    from tessera.adapters.schedulers import get_adapter

    p = load_profile("example-1m")
    doc = dict(gw170817_plan.schedules["example-1m@1"], feasible=False)
    with pytest.raises(AdapterError) as e:
        get_adapter("lco").export(doc, p)
    assert e.value.code == "SCHEDULE_INFEASIBLE"
    good = gw170817_plan.schedules["example-1m@1"]
    exp = get_adapter("lco").export(good, p)
    with pytest.raises(AdapterError) as e:
        get_adapter("lco").submit(exp, None, good, p)
    assert e.value.code == "SUBMISSION_DISABLED"
    with pytest.raises(AdapterError):
        get_adapter("rubin_too").submit()


# -- alerts --------------------------------------------------------------------

def test_lvk_alert_parsing():
    from tessera.adapters.alerts import AlertError, make_alert, parse_lvk_json

    raw = GW170817_MAP.read_bytes()[:100]
    env = parse_lvk_json(json.dumps(make_alert("S230518h", "PRELIMINARY", "2023-05-18T12:59:08Z", raw,
                                               classification={"NSBH": 0.86, "Terrestrial": 0.14})))
    assert env.skymap_bytes == raw and env.classification["NSBH"] == 0.86 and not env.is_mock
    r = parse_lvk_json(make_alert("S230518h", "RETRACTION", "2023-05-18T14:00:00Z"))
    assert r.skymap_bytes is None and r.alert_key != env.alert_key
    assert parse_lvk_json(make_alert("MS230518h", "PRELIMINARY", "2023-05-18T12:59:08Z", raw)).is_mock
    with pytest.raises(AlertError):
        parse_lvk_json(b"{not json")
    with pytest.raises(AlertError):
        parse_lvk_json({"superevent_id": "S1", "alert_type": "BOGUS", "time_created": "x"})


def test_url_fetch_allowlist():
    from tessera.adapters.alerts.sources import FetchError, fetch_url

    for url in ("http://gracedb.ligo.org/x", "https://example.com/x.fits", "file:///etc/passwd"):
        with pytest.raises(FetchError):
            fetch_url(url)


# -- reproducibility -----------------------------------------------------------

def test_plan_is_deterministic(gw170817_map, glade_snapshot, gw170817_plan):
    from tessera import PlanConfig, plan
    from tessera.core.hashing import sha256_bytes

    again = plan(gw170817_map, glade_snapshot, ["example-1m"], PlanConfig(preset="bns", horizon_nights=1),
                 T_GW170817, event={"superevent_id": "GW170817", "revision": 1})
    assert again.plan_key == gw170817_plan.plan_key
    a = {k: sha256_bytes(v) for k, v in gw170817_plan.outputs().items() if k != "manifest.json"}
    assert a == again.manifest["outputs"]


def test_plan_key_tracks_inputs(gw170817_map, glade_snapshot, gw170817_plan):
    from tessera import PlanConfig, plan

    other = plan(gw170817_map, glade_snapshot, [], PlanConfig(preset="bns", horizon_nights=1), T_GW170817,
                 event={"superevent_id": "GW170817", "revision": 1})
    assert other.plan_key != gw170817_plan.plan_key
    m = gw170817_plan.manifest
    for k in ("skymap", "catalog", "profiles", "config", "versions", "environment", "tolerance", "request_config"):
        assert k in m
    assert m["catalog"]["snapshot_id"] == glade_snapshot.id


def test_cli_plan_validate_export_rerun(tmp_path, glade_snapshot):
    from tessera.cli import main

    out = tmp_path / "out"
    assert main(["plan", "--skymap", str(GW170817_MAP), "--catalog", str(glade_snapshot.path), "--profile",
                 "example-1m", "--preset", "bns", "--nights", "1", "--t-start", T_GW170817, "--out", str(out)]) == 0
    sched = out / "schedules" / "example-1m@1.json"
    assert main(["validate-schedule", str(sched)]) == 0
    assert main(["export", str(sched), "--adapter", "acp", "--out", str(tmp_path / "plan.txt")]) == 0
    assert main(["rerun", str(out / "manifest.json"), "--skymap", str(GW170817_MAP), "--catalog",
                 str(glade_snapshot.path)]) == 0


def test_cli_ingest(tmp_path):
    from tessera.cli import main

    assert main(["ingest-catalog", "glade+", str(GLADE_EXCERPT), "--out", str(tmp_path / "snap"),
                 "--upstream-version", "2022-10-31"]) == 0
    assert (tmp_path / "snap" / "snapshot.json").exists()


def test_to_time_accepts_common_utc_forms_and_rejects_others():
    from tessera.core.errors import ConfigError
    from tessera.core.timeutil import iso, to_time

    for s in ("2017-08-17T22:00:00", "2017-08-17T22:00:00Z", "2017-08-17 22:00 UTC", "2017-08-17 22:00"):
        assert iso(to_time(s)) == "2017-08-17T22:00:00Z"
    for bad in ("next tuesday", "22:00", "2017-13-40T99:00"):
        with pytest.raises(ConfigError) as e:
            to_time(bad)
        assert e.value.code == "INVALID_TIME" and e.value.http_status == 422


def test_ztf_profile_matches_published_camera():
    from tessera.core.footprint import InstrumentFootprint
    from tessera.core.profile import load_profile

    p = load_profile("ztf-approx")
    fp = InstrumentFootprint(p.instrument.footprint)
    assert len(p.instrument.footprint.chips) == 16                 # 4 x 4 CCD mosaic
    assert abs(fp.area_deg2 - 47.7) < 0.2                          # Bellm et al. 2019: 47 deg2 with gaps
    assert p.mount.ha_limits_h == (-5.95, 5.95) and p.constraints.airmass_max == 2.5
    # ztf_sim positional rules: south of -22 deg only west of HA -17.6 deg
    assert p.mount.ha_ok(-2.0, -10.0) and not p.mount.ha_ok(-2.0, -30.0) and p.mount.ha_ok(0.5, -30.0)


def test_glade_plus_cds_fixed_width_reader():
    from tessera.adapters.catalogs import glade_plus

    c = glade_plus.read(DATA / "gladep_cds_sample.dat")     # first 3 records of VizieR VII/291 gladep.dat
    assert len(c["ra"]) == 3
    assert c["gal_id"][0] == "GLADE+:1" and "GWGC:NGC4736" in c["alt_ids"][0]
    assert abs(c["ra"][0] - 192.721451) < 1e-6 and abs(c["dec"][0] - 41.120152) < 1e-6
    assert abs(c["mag_B"][0] - 8.8) < 1e-9 and c["dist_kind"][0] == "specz"

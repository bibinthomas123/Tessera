"""End-to-end service tests: API + worker (inline) + SQLite + alert revisions."""

from __future__ import annotations

import json

import pytest

pytest.importorskip("fastapi")

from fastapi.testclient import TestClient  # noqa: E402

from tessera.testing import gaussian_map_bytes  # noqa: E402

from ..conftest import GLADE_EXCERPT  # noqa: E402

T0 = "2017-08-17T22:00:00"


@pytest.fixture()
def svc(tmp_path, monkeypatch):
    from tessera.adapters.catalogs import glade_plus
    from tessera.core.catalog import build_snapshot

    cats = tmp_path / "catalogs"
    build_snapshot(glade_plus.read(GLADE_EXCERPT), cats / "glade", name="glade", upstream_version="2022-10-31")
    env = {
        "TESSERA_DATABASE_URL": f"sqlite:///{tmp_path}/tessera.db",
        "TESSERA_DATA_DIR": str(tmp_path / "data"),
        "TESSERA_CATALOG_DIRS": str(cats),
        "TESSERA_DEFAULT_CATALOG": "glade",
        "TESSERA_INLINE_WORKER": "1",
        "TESSERA_MODE": "test",
        "TESSERA_MAX_UPLOAD_BYTES": str(5 * 1024**2),
    }
    policy = tmp_path / "policy.yaml"
    policy.write_text("include_mock: true\nrules:\n  - {preset: bns, profiles: [example-1m], "
                      "config: {horizon_nights: 1, compare_strategies: false}}\n")
    env["TESSERA_TRIGGER_POLICY"] = str(policy)
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    from tessera.service import settings as st

    st.get_settings.cache_clear()
    from tessera.service.api import create_app
    from tessera.service.auth import create_key
    from tessera.service.db import session_factory

    s = st.get_settings()
    app = create_app(s)
    with session_factory(s.database_url)() as sess:
        keys = {
            "admin": create_key(sess, "admin", ["admin"]),
            "reader": create_key(sess, "reader", ["read"]),
            "planner": create_key(sess, "planner", ["read", "plan"]),
            "approver": create_key(sess, "approver", ["read", "approve"]),
        }
    client = TestClient(app)
    yield client, keys
    st.get_settings.cache_clear()


def H(key):
    return {"Authorization": f"Bearer {key}"}


def _upload(client, key, data, name="map.fits"):
    r = client.post("/v1/skymaps", files={"file": (name, data, "application/fits")}, headers=H(key))
    assert r.status_code == 201, r.text
    return r.json()["skymap_id"]


def _plan(client, key, sid, **extra):
    body = {"skymap_id": sid, "profiles": ["example-1m"], "t_start": T0,
            "config": {"preset": "bns", "horizon_nights": 1, "compare_strategies": False}, **extra}
    return client.post("/v1/plans", json=body, headers=H(key))


def test_auth_and_scopes(svc):
    client, keys = svc
    assert client.get("/healthz").status_code == 200
    assert client.get("/v1/profiles").status_code == 401
    assert client.get("/v1/profiles", headers=H("tessera_bad_key")).status_code == 401
    assert client.get("/v1/profiles", headers=H(keys["reader"])).status_code == 200
    r = client.post("/v1/skymaps", files={"file": ("m.fits", b"SIMPLE  =", "x")}, headers=H(keys["reader"]))
    assert r.status_code == 403 and r.json()["code"] == "FORBIDDEN"


def test_upload_limits_and_url_allowlist(svc):
    client, keys = svc
    big = b"SIMPLE  =" + b"\0" * (6 * 1024**2)
    r = client.post("/v1/skymaps", files={"file": ("big.fits", big, "x")}, headers=H(keys["planner"]))
    assert r.status_code in (413, 422)
    r = client.post("/v1/skymaps", data={"url": "https://evil.example.com/x.fits"}, headers=H(keys["planner"]))
    assert r.status_code == 422 and r.json()["code"] == "SKYMAP_URL_REFUSED"
    r = client.post("/v1/skymaps", files={"file": ("x.txt", b"hello", "x")}, headers=H(keys["planner"]))
    assert r.json()["code"] == "SKYMAP_INVALID_FORMAT"


def test_catalog_listing_names_are_requestable(svc):
    """The UI sends back the listed name; it must resolve (regression: snapshot name shadowed dir name)."""
    client, keys = svc
    cats = client.get("/v1/catalogs", headers=H(keys["reader"])).json()["catalogs"]
    assert [c["name"] for c in cats] == ["glade"]
    sid = _upload(client, keys["planner"], gaussian_map_bytes([(197.45, -23.38, 2.0, 1.0)], order=8))
    assert _plan(client, keys["planner"], sid, catalog=cats[0]["name"]).status_code in (200, 202)


def test_plan_lifecycle(svc):
    client, keys = svc
    sid = _upload(client, keys["planner"], gaussian_map_bytes([(197.45, -23.38, 2.0, 1.0)], order=8))
    r = _plan(client, keys["planner"], sid)
    assert r.status_code == 202, r.text
    job = r.json()
    assert job["state"] == "succeeded", job
    key = job["plan_key"]
    # idempotent: same request returns the existing plan
    r2 = _plan(client, keys["planner"], sid)
    assert r2.status_code == 200 and r2.json()["plan_key"] == key and r2.json()["existing"]

    p = client.get(f"/v1/plans/{key}", headers=H(keys["reader"])).json()
    assert p["budget"]["p_cat"] is not None
    assert p["schedules"] and p["schedules"][0]["status"] == "proposed"
    csv = client.get(f"/v1/plans/{key}/galaxies?format=csv", headers=H(keys["reader"]))
    assert csv.status_code == 200 and csv.text.startswith("rank,tier,gal_id")
    for fmt in ("json", "ecsv", "parquet", "votable"):
        assert client.get(f"/v1/plans/{key}/galaxies?format={fmt}", headers=H(keys["reader"])).status_code == 200
    man = client.get(f"/v1/plans/{key}/manifest", headers=H(keys["reader"])).json()
    assert man["plan_key"] == key and man["catalog"]["snapshot_id"]
    assert client.get(f"/v1/plans/{key}/files/../../etc/passwd", headers=H(keys["reader"])).status_code == 404

    sch_id = p["schedules"][0]["schedule_id"]
    sch = client.get(f"/v1/schedules/{sch_id}", headers=H(keys["reader"])).json()
    assert sch["doc"]["feasible"] and sch["version"] == 0
    for adapter in ("generic", "acp", "treasuremap", "lco", "rubin_too"):
        r = client.get(f"/v1/schedules/{sch_id}/versions/0/export?adapter={adapter}", headers=H(keys["reader"]))
        assert r.status_code == 200, (adapter, r.text)

    # edit: move the first visit into daylight -> new version flagged infeasible; v0 untouched
    visits = sch["doc"]["visits"]
    bad = [dict(v) for v in visits]
    bad[0]["start"], bad[0]["end"] = "2017-08-17T18:00:00Z", "2017-08-17T18:01:00Z"
    r = client.post(f"/v1/schedules/{sch_id}/versions", json={"base_version": 0, "visits": bad},
                    headers=H(keys["planner"]))
    assert r.status_code == 201
    v1 = r.json()
    assert v1["version"] == 1 and not v1["doc"]["feasible"]
    assert v1["diff_from_parent"]["changed"]
    r = client.post(f"/v1/schedules/{sch_id}/versions/1/approve", headers=H(keys["approver"]))
    assert r.status_code == 409
    v0 = client.get(f"/v1/schedules/{sch_id}/versions/0", headers=H(keys["reader"])).json()
    assert v0["doc"]["visits"] == visits

    # a clean edit (drop the last visit) can be approved; stale base is refused
    r = client.post(f"/v1/schedules/{sch_id}/versions", json={"base_version": 1, "visits": visits[:-1]},
                    headers=H(keys["planner"]))
    assert r.status_code == 201 and r.json()["doc"]["feasible"]
    r = client.post(f"/v1/schedules/{sch_id}/versions", json={"base_version": 0, "visits": visits},
                    headers=H(keys["planner"]))
    assert r.status_code == 409
    assert client.post(f"/v1/schedules/{sch_id}/versions/2/approve", headers=H(keys["planner"])).status_code == 403
    r = client.post(f"/v1/schedules/{sch_id}/versions/2/approve", headers=H(keys["approver"]))
    assert r.status_code == 200 and r.json()["status"] == "approved"
    # submission stays disabled for the packaged profile
    r = client.post(f"/v1/schedules/{sch_id}/versions/2/submit?adapter=lco", headers=H(keys["approver"]))
    assert r.status_code in (409, 422)

    # rerun reproduces identical outputs
    r = client.post(f"/v1/plans/{key}/rerun", headers=H(keys["planner"])).json()
    assert r["state"] == "succeeded" and r["result"]["identical"], r

    # transient association near the top galaxy
    gal = client.get(f"/v1/plans/{key}/galaxies?format=json", headers=H(keys["reader"])).json()["rows"][0]
    a = client.post("/v1/transients/associate", json={"plan_key": key, "ra": gal["ra"] + 1e-3, "dec": gal["dec"]},
                    headers=H(keys["reader"])).json()
    assert a["galaxies"][0]["gal_id"] == gal["gal_id"]


def _alert(client, key, sid, atype, tc, data=None, **ev):
    from tessera.adapters.alerts.envelope import make_alert

    doc = make_alert(sid, atype, tc, data, **ev)
    r = client.post("/v1/alerts", content=json.dumps(doc), headers=H(key) | {"content-type": "application/json"})
    assert r.status_code == 200, r.text
    return r.json()


def test_alert_revisions_and_retraction(svc):
    client, keys = svc
    a = keys["admin"]
    m1 = gaussian_map_bytes([(197.45, -23.38, 2.0, 1.0)], order=8)
    m2 = gaussian_map_bytes([(197.0, -23.0, 1.5, 1.0)], order=8)
    sid = "MS170817a"
    o1 = _alert(client, a, sid, "PRELIMINARY", "2017-08-17T12:45:00Z", m1)
    assert o1["outcome"] == "planned" and o1["jobs"]
    assert _alert(client, a, sid, "PRELIMINARY", "2017-08-17T12:45:00Z", m1)["outcome"] == "duplicate"
    assert _alert(client, a, sid, "INITIAL", "2017-08-17T13:30:00Z", m1)["outcome"] == "metadata_updated"
    o3 = _alert(client, a, sid, "UPDATE", "2017-08-17T15:00:00Z", m2)
    assert o3["outcome"] == "planned" and o3["revision_id"] != o1["revision_id"]
    late = _alert(client, a, sid, "UPDATE", "2017-08-17T14:00:00Z",
                  gaussian_map_bytes([(197.2, -23.1, 1.8, 1.0)], order=8))
    assert late["outcome"] == "late"
    ev = client.get(f"/v1/events/{sid}", headers=H(keys["reader"])).json()
    states = {r["revision_id"]: r["state"] for r in ev["revisions"]}
    assert states[o1["revision_id"]] == "superseded" and states[o3["revision_id"]] == "active"
    assert "delta_from_previous" in ev["revisions"][1]
    # revision 1 schedules superseded, revision 2 proposed
    plans = client.get("/v1/plans", headers=H(keys["reader"])).json()["plans"]
    assert len(plans) == 2
    r = _alert(client, a, sid, "RETRACTION", "2017-08-17T16:00:00Z")
    assert r["outcome"] == "retracted"
    assert _alert(client, a, sid, "RETRACTION", "2017-08-17T16:05:00Z")["outcome"] == "retraction_duplicate"
    statuses = []
    for p in plans:
        detail = client.get(f"/v1/plans/{p['plan_key']}", headers=H(keys["reader"])).json()
        statuses += [s["status"] for s in detail["schedules"]]
    # the live revision's schedules are withdrawn; the earlier revision's stay superseded
    assert sorted(statuses) == ["superseded", "withdrawn"], statuses
    after = _alert(client, a, sid, "UPDATE", "2017-08-17T17:00:00Z",
                   gaussian_map_bytes([(197.3, -23.2, 1.0, 1.0)], order=8))
    assert after["outcome"] == "after_retraction"


def test_mock_events_ignored_in_production(svc, monkeypatch):
    from tessera.adapters.alerts.envelope import make_alert, parse_lvk_json
    from tessera.service.events import TriggerPolicy

    pol = TriggerPolicy.load(None, "production")
    env = parse_lvk_json(make_alert("MS1", "PRELIMINARY", "2020-01-01T00:00:00Z", b"x"))
    assert pol.rules_for(env)[1] == "mock event ignored outside test mode"


def test_observation_feedback_replans(svc):
    client, keys = svc
    a = keys["admin"]
    sid = "MS200105ae"
    o = _alert(client, a, sid, "PRELIMINARY", "2017-08-17T12:45:00Z",
               gaussian_map_bytes([(197.45, -23.38, 2.0, 1.0)], order=8))
    assert o["outcome"] == "planned"
    key = client.get("/v1/plans", headers=H(a)).json()["plans"][0]["plan_key"]
    sched = client.get(f"/v1/plans/{key}", headers=H(a)).json()["schedules"][0]
    v = client.get(f"/v1/schedules/{sched['schedule_id']}", headers=H(a)).json()["doc"]["visits"][0]
    r = client.post("/v1/observations", json={"superevent_id": sid, "profile_id": "example-1m", "ra": v["ra"],
                                              "dec": v["dec"], "filter": v["filter"], "start": v["start"],
                                              "end": v["end"], "status": "completed", "replan": True},
                    headers=H(a))
    assert r.status_code == 201, r.text
    rp = r.json()["replan"]
    assert rp["plan_key"] != key


def test_bad_start_time_is_a_422_not_a_500(svc):
    client, keys = svc
    sid = _upload(client, keys["planner"], gaussian_map_bytes([(197.45, -23.38, 2.0, 1.0)], order=8))
    r = _plan(client, keys["planner"], sid, t_start="next tuesday")
    assert r.status_code == 422, r.text
    assert r.json()["code"] == "INVALID_TIME"
    # human-friendly forms are accepted
    r = _plan(client, keys["planner"], sid, t_start="2017-08-17 22:00 UTC")
    assert r.status_code in (200, 201, 202), r.text

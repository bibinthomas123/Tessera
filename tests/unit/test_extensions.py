"""Filter blocks, optimality bound, OIDC, compare, facility-format conformance, light curves."""

from __future__ import annotations

import json
import shutil
import time

import numpy as np
import pytest

from tessera.core.profile import load_profile

from ..conftest import GW170817_MAP, T_GW170817

# -- light curve ------------------------------------------------------------------

def test_at2017gfo_table_is_default_and_sane():
    from tessera.core import science
    from tessera.core.detectability import LightCurve

    assert science.preset("bns")["light_curve"] == "at2017gfo"
    lc = LightCurve.named("at2017gfo")
    r = [float(lc.abs_mag(t, "r")) for t in (0.5, 1.0, 2.0, 4.0, 7.0)]
    assert -16.5 < r[0] < -15.5                     # AT2017gfo peaked near M_r ~ -16
    assert all(b > a for a, b in zip(r, r[1:]))     # fading
    # red bands fade more slowly than blue (the red kilonova component)
    assert float(lc.abs_mag(5.0, "g")) > float(lc.abs_mag(5.0, "z"))
    assert float(lc.abs_mag(12.0, "r")) > float(lc.abs_mag(7.5, "r"))   # extrapolated decline
    assert np.isinf(float(lc.abs_mag(-1.0, "r")))


def test_reference_densities_are_sourced():
    from tessera.core import science

    k = science.reference_density("lum_K")
    assert k["value"] == pytest.approx(4.01e8) and "Kochanek" in k["source"] and not k["provisional"]
    # reproduce the derivation from the published Schechter fit
    from scipy.special import gamma

    j = 0.0116 * 10 ** (-0.4 * (-23.39 - 3.28)) * gamma(-1.09 + 2) * 0.7
    assert j == pytest.approx(k["value"], rel=0.01)
    m = science.reference_density("mstar")
    assert m["value"] == pytest.approx(0.0020 * 2.775e11 * 0.7, rel=0.01) and m["provisional"]


# -- filter blocks ------------------------------------------------------------------

def test_filter_blocks_reduce_filter_changes(glade_snapshot):
    from tessera import PlanConfig, plan
    from tessera.testing import gaussian_map

    m = gaussian_map([(197.45, -23.38, 4.0, 1.0)], order=8, dist=(150, 40))
    out = {}
    for fb in (None, 60):
        r = plan(m, glade_snapshot, ["rubin-lsstcam-approx"],
                 PlanConfig(preset="bns", horizon_nights=1, compare_strategies=False,
                            scheduler={"filter_block_min": fb}), T_GW170817)
        s = r.schedules["rubin-lsstcam-approx@1"]
        assert s["feasible"]
        out[fb] = s["summary"]["n_filter_changes"]
        if fb:
            assert any(v["decision"]["rule"].startswith("filter_block") for v in s["visits"])
    assert out[60] < out[None]


# -- optimality ---------------------------------------------------------------------

def test_greedy_within_bound(gw170817_map, glade_snapshot):
    pytest.importorskip("ortools")
    from astropy.time import Time

    from tessera.core.association import score_galaxies
    from tessera.core.config import PlanConfig
    from tessera.core.fields import build_fields
    from tessera.core.optimality import greedy_vs_bound
    from tessera.core.visibility import compute_visibility

    m = gw170817_map
    cfg = PlanConfig(preset="bns", horizon_nights=1).resolved()
    reg = m.region_mask_inclusive(cfg.cr_max)
    a = score_galaxies(m, glade_snapshot.extract_for_map(m.order[reg], m.ipix[reg]), cfg)
    S = np.where(np.asarray(a.table["tier"]) == 1, np.asarray(a.table["S"]), 0.0)
    p = load_profile("example-1m")
    fs = build_fields(m, p, cfg, a.table, host_norm=float(S.sum()), host_metric="score_fraction_in_field",
                      mode="tessellate")
    t0 = Time(T_GW170817)
    vis = compute_visibility(p, fs.ra, fs.dec, t0, t0 + 1)
    r = greedy_vs_bound(fs, vis, cfg, t_start=t0, t_end=t0 + 1, t_merger=Time(m.meta.mjd_obs, format="mjd"),
                        galaxies_S=S, objective="sky_prob")
    assert r["status"] == "optimal"
    assert r["greedy"] <= r["upper_bound"] + 1e-9          # a valid bound is never beaten
    assert r["ratio_lower_bound"] >= 0.9                   # spec target: >= 90% of optimum


# -- OIDC ---------------------------------------------------------------------------

@pytest.fixture()
def oidc_env(tmp_path, monkeypatch):
    jwt = pytest.importorskip("jwt")
    from cryptography.hazmat.primitives.asymmetric import rsa

    key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
    jwk = json.loads(jwt.algorithms.RSAAlgorithm.to_jwk(key.public_key()))
    jwk.update({"kid": "k1", "use": "sig", "alg": "RS256"})
    (tmp_path / "jwks.json").write_text(json.dumps({"keys": [jwk]}))
    monkeypatch.setenv("TESSERA_OIDC_ISSUER", "https://idp.example")
    monkeypatch.setenv("TESSERA_OIDC_AUDIENCE", "tessera")
    monkeypatch.setenv("TESSERA_OIDC_JWKS_FILE", str(tmp_path / "jwks.json"))

    def token(scope="tessera:read tessera:plan", aud="tessera", exp=600, sub="alice"):
        now = int(time.time())
        return jwt.encode({"iss": "https://idp.example", "aud": aud, "sub": sub, "preferred_username": sub,
                           "scope": scope, "iat": now, "exp": now + exp}, key, algorithm="RS256",
                          headers={"kid": "k1"})
    return token


def test_oidc_verifier(oidc_env):
    from tessera.service.auth import AuthError
    from tessera.service.oidc import OIDCConfig, OIDCVerifier

    v = OIDCVerifier(OIDCConfig.from_env())
    p = v.verify(oidc_env())
    assert p.name == "oidc:alice" and p.scopes == {"read", "plan"} and p.via == "oidc"
    for bad in (oidc_env(aud="other"), oidc_env(exp=-3600), oidc_env()[:-4] + "AAAA"):
        with pytest.raises(AuthError):
            v.verify(bad)


def test_oidc_through_api(oidc_env, tmp_path, monkeypatch):
    pytest.importorskip("fastapi")
    from fastapi.testclient import TestClient

    from tessera.service import settings as st

    monkeypatch.setenv("TESSERA_DATABASE_URL", f"sqlite:///{tmp_path}/t.db")
    monkeypatch.setenv("TESSERA_DATA_DIR", str(tmp_path / "d"))
    st.get_settings.cache_clear()
    from tessera.service.api import create_app

    c = TestClient(create_app(st.get_settings()))
    assert c.get("/v1/profiles", headers={"Authorization": f"Bearer {oidc_env()}"}).status_code == 200
    r = c.get("/v1/audit", headers={"Authorization": f"Bearer {oidc_env()}"})
    assert r.status_code == 403
    assert c.get("/v1/profiles", headers={"Authorization": f"Bearer {oidc_env(aud='x')}"}).status_code == 401
    st.get_settings.cache_clear()


# -- compare ------------------------------------------------------------------------

def test_compare_identical_and_perturbed(gw170817_plan, tmp_path):
    from tessera.core.compare import compare_plans

    a, b = tmp_path / "a", tmp_path / "b"
    gw170817_plan.write(a)
    gw170817_plan.write(b)
    assert compare_plans(a, b)["verdict"] == "identical"
    # perturb group B: a different software version and a tiny score change
    man = json.loads((b / "manifest.json").read_text())
    man["versions"]["astropy"] = "0.0-other"
    man["outputs"]["galaxies.json"] = "0" * 64
    (b / "manifest.json").write_text(json.dumps(man))
    gal = json.loads((b / "galaxies.json").read_text())
    gal["rows"][0]["S"] *= 1 + 1e-12
    (b / "galaxies.json").write_text(json.dumps(gal))
    rep = compare_plans(a, b)
    assert rep["verdict"] == "within_tolerance"
    assert any("astropy" in c for c in rep["likely_causes"])
    gal["rows"][0]["S"] *= 1.01
    (b / "galaxies.json").write_text(json.dumps(gal))
    assert compare_plans(a, b)["verdict"] == "different"
    shutil.rmtree(b)


# -- facility formats -----------------------------------------------------------------

def test_rubin_case_mapping():
    from tessera.adapters.schedulers.rubin_too import classify

    assert classify(25, {"BNS": 0.99}, {"HasNS": 1})[0] == "GW_case_A"
    assert classify(240, {"NSBH": 0.9}, {"HasNS": 1})[0] == "GW_case_D"
    assert classify(800, {"BNS": 0.99}, {"HasNS": 1})[0] is None
    assert classify(15, {"BBH": 0.99}, {"HasNS": 0})[0] == "BBH_case_B"
    assert classify(50, {"BBH": 0.99}, {"HasNS": 0})[0] is None


def test_rubin_export_conforms_to_rubin_scheduler(gw170817_map, gw170817_plan):
    from tessera.adapters.schedulers import get_adapter

    ad = get_adapter("rubin_too")
    exp = ad.export(gw170817_plan.schedules["example-1m@1"], load_profile("rubin-lsstcam-approx"),
                    skymap=gw170817_map, event={"superevent_id": "S170817a", "revision": 1,
                                                "classification": {"BNS": 0.99}, "properties": {"HasNS": 1.0}})
    req = json.loads(exp.files["too_request.json"])
    assert req["targetoo"]["too_type"] == "GW_case_A" and req["eligible"]
    assert req["targetoo"]["posterior_distance"] == pytest.approx(38034.08, rel=1e-3)   # kpc
    issues = ad.validate(exp)
    assert not [i for i in issues if i.level == "error"], issues
    if pytest.importorskip("rubin_scheduler", reason="conformance needs rubin_scheduler"):
        assert any("followed by rubin_scheduler" in i.message for i in issues)


def test_lco_schema_rules(gw170817_plan):
    from tessera.adapters.schedulers import get_adapter

    p = load_profile("example-1m")
    ad = get_adapter("lco")
    doc = gw170817_plan.schedules["example-1m@1"]
    one = dict(doc, visits=doc["visits"][:1])
    exp = ad.export(one, p)
    rg = json.loads(exp.files["requestgroup.json"])
    assert rg["operator"] == "SINGLE"
    cfg = rg["requests"][0]["configurations"][0]
    assert cfg["acquisition_config"] == {"mode": "OFF"} and "guiding_config" in cfg
    assert rg["requests"][0]["optimization_type"] == "TIME"
    assert not [i for i in ad.validate(exp) if i.level == "error"]
    rg["ipp_value"] = 5.0
    rg["requests"][0]["configurations"][0]["constraints"]["max_airmass"] = 40
    exp.files["requestgroup.json"] = json.dumps(rg).encode()
    msgs = [i.message for i in ad.validate(exp)]
    assert any("ipp_value" in m for m in msgs) and any("max_airmass" in m for m in msgs)


def test_treasuremap_token_sent_as_header(gw170817_plan, monkeypatch):
    from tessera.adapters.schedulers import get_adapter
    from tessera.adapters.schedulers.base import Approval, schedule_sha256

    p = load_profile("example-1m").model_copy(update={"submission_enabled": True})
    doc = gw170817_plan.schedules["example-1m@1"]
    ad = get_adapter("treasuremap")
    exp = ad.export(doc, p, graceid="S170817a")
    body = json.loads(exp.files["treasuremap_pointings.json"])
    assert "api_token" not in body and "position" not in body["pointings"][0]
    assert isinstance(body["pointings"][0]["instrumentid"], int)
    seen = {}

    class R:
        status_code, text = 200, "{}"

    def fake_post(url, json=None, headers=None, timeout=None):
        seen.update(url=url, headers=headers, json=json)
        return R()

    import httpx

    monkeypatch.setattr(httpx, "post", fake_post)
    monkeypatch.setenv("TESSERA_TREASUREMAP_TOKEN", "fake-token-for-test")
    appr = Approval(doc["schedule_id"], 0, schedule_sha256(doc), "tester", "2017-08-17T22:00:00Z")
    ad.submit(exp, appr, doc, p)
    assert seen["headers"] == {"api_token": "fake-token-for-test"} and "api_token" not in seen["json"]
    assert seen["url"] == "https://treasuremap.space/api/v1/pointings"


def test_acp_flags_unconfirmed_syntax(gw170817_plan):
    from tessera.adapters.schedulers import get_adapter

    p = load_profile("example-1m")
    ad = get_adapter("acp")
    doc = gw170817_plan.schedules["example-1m@1"]
    assert any(i.level == "warning" and "not yet confirmed" in i.message for i in ad.validate(ad.export(doc, p)))
    ok = ad.export(doc, p, syntax_confirmed=True)
    assert not [i for i in ad.validate(ok) if i.level != "info"]
    assert GW170817_MAP.exists()

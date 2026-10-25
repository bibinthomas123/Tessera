"""Science-model tests: closed forms, scoring rules, flags (spec section 3)."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import integrate

from tessera.core.association import score_galaxies
from tessera.core.catalog import Galaxies, normalise_columns
from tessera.core.config import PlanConfig
from tessera.core.distance import conditional_pdf, gauss, shell_probability
from tessera.testing import distnorm, gaussian_map

from ..conftest import NGC4993


def test_distnorm_normalises_ansatz():
    for mu, sg in ((40, 8), (300, 150), (10, 30), (-20, 40)):
        n = distnorm(mu, sg)
        val, _ = integrate.quad(lambda r: conditional_pdf(r, mu, sg, n), 0, np.inf)
        assert val == pytest.approx(1.0, rel=1e-7)


@pytest.mark.parametrize("mu,sg,lo,hi", [(40, 8, 30, 50), (40, 8, 0, 10), (300, 150, 100, 1000), (5, 20, 0, 3)])
def test_shell_probability_closed_form(mu, sg, lo, hi):
    n = distnorm(mu, sg)
    num, _ = integrate.quad(lambda r: conditional_pdf(r, mu, sg, n), lo, hi)
    assert float(shell_probability(mu, sg, n, lo, hi)) == pytest.approx(num, rel=1e-9, abs=1e-14)


def test_pdist_marginalisation_closed_form():
    """N * Gauss(d; mu, sqrt(sigma^2 + s^2)) == integral of N Gauss(r; mu, sigma) Gauss(r; d, s) dr."""
    mu, sg, s, d = 41.0, 7.5, 4.4, 37.0
    n = distnorm(mu, sg)
    num, _ = integrate.quad(lambda r: n * gauss(r, mu, sg) * gauss(r, d, s), -np.inf, np.inf)
    assert float(n * gauss(d, mu, np.hypot(sg, s))) == pytest.approx(num, rel=1e-9)


def _cols(n, **over):
    base = {"gal_id": np.array([f"T:{i:03d}" for i in range(n)], dtype=object), "ra": np.full(n, 197.45),
            "dec": np.full(n, -23.38), "dist_mpc": np.full(n, 40.0), "dist_err_mpc": np.full(n, 2.0),
            "dist_kind": np.full(n, "z_indep", dtype=object), "mag_K": np.full(n, 10.0)}
    base.update(over)
    c = normalise_columns(base)
    from tessera.core.catalog import _hpx

    c["hpx10"] = _hpx(c["ra"], c["dec"], 10)
    return Galaxies(c)


def test_tie_break_is_by_gal_id(synthetic_map):
    g = _cols(5, gal_id=np.array(["T:e", "T:b", "T:d", "T:a", "T:c"], dtype=object))
    res = score_galaxies(synthetic_map, g, PlanConfig().resolved())
    assert list(res.table["gal_id"]) == ["T:a", "T:b", "T:c", "T:d", "T:e"]


def test_distance_consistency_orders_galaxies(synthetic_map):
    g = _cols(3, dist_mpc=np.array([40.0, 80.0, 15.0]))
    t = score_galaxies(synthetic_map, g, PlanConfig(strategy="vol3d-flat").resolved()).table
    assert list(t["gal_id"])[0] == "T:000"
    assert t["p_dist_mpc3"][0] > t["p_dist_mpc3"][1]


def test_luminosity_weight_orders_galaxies(synthetic_map):
    g = _cols(3, mag_K=np.array([12.0, 8.0, 10.0]))
    t = score_galaxies(synthetic_map, g, PlanConfig(strategy="vol3d-lumK").resolved()).table
    assert list(t["gal_id"]) == ["T:001", "T:002", "T:000"]


def test_no_distance_galaxies_are_tier2_never_interleaved(synthetic_map):
    g = _cols(4, dist_mpc=np.array([40.0, np.nan, 41.0, np.nan]),
              dist_kind=np.array(["z_indep", "none", "z_indep", "none"], dtype=object))
    t = score_galaxies(synthetic_map, g, PlanConfig().resolved()).table
    assert list(t["tier"]) == [1, 1, 2, 2]
    assert all("NO_DIST" in f for f in t["flags"][2:])
    t2 = score_galaxies(synthetic_map, g, PlanConfig(missing_distance="exclude").resolved()).table
    assert set(t2["tier"]) == {1}


def test_qso_and_duplicates_excluded(synthetic_map):
    g = _cols(3, obj_type=np.array(["galaxy", "qso", "galaxy"], dtype=object),
              duplicate_of=np.array(["", "", "T:000"], dtype=object))
    res = score_galaxies(synthetic_map, g, PlanConfig().resolved())
    assert list(res.table["gal_id"]) == ["T:000"]
    assert res.summary["n_qso_excluded"] == 1 and res.summary["n_duplicates_skipped"] == 1


def test_missing_weight_is_imputed_pessimistically(synthetic_map):
    g = _cols(4, mag_K=np.array([9.0, 10.0, 11.0, np.nan]))
    t = score_galaxies(synthetic_map, g, PlanConfig(strategy="vol3d-lumK").resolved()).table
    row = t[list(t["gal_id"]).index("T:003")]
    assert "W_IMPUTED" in row["flags"]
    assert row["w"] <= min(t["w"][[i for i, x in enumerate(t["gal_id"]) if x != "T:003"]])


def test_photoz_flag_and_wider_uncertainty(synthetic_map):
    g = _cols(2, dist_kind=np.array(["specz", "photoz"], dtype=object), dist_err_mpc=np.full(2, np.nan),
              z=np.full(2, 0.0093), dist_mpc=np.full(2, 40.0))
    t = score_galaxies(synthetic_map, g, PlanConfig().resolved()).table
    ph = t[list(t["gal_id"]).index("T:001")]
    sp = t[list(t["gal_id"]).index("T:000")]
    assert "PHOTOZ" in ph["flags"] and ph["dist_err_mpc"] > sp["dist_err_mpc"]


def test_sky2d_strategy_ignores_distance(synthetic_map):
    g = _cols(2, dist_mpc=np.array([40.0, 400.0]))
    t = score_galaxies(synthetic_map, g, PlanConfig(strategy="sky2d-flat").resolved()).table
    assert np.allclose(t["p_dist_mpc3"], 1.0)


def test_event_without_distance_falls_back_to_2d():
    m = gaussian_map([(197.45, -23.38, 2.0, 1.0)], order=7, dist=None)
    g = _cols(2, dist_mpc=np.array([40.0, 400.0]))
    res = score_galaxies(m, g, PlanConfig().resolved())
    assert res.summary["distance_used"] is False
    assert np.allclose(res.table["p_dist_mpc3"], 1.0)


def test_score_fraction_labelled_not_probability(gw170817_plan):
    t = gw170817_plan.galaxies
    assert "NOT a host probability" in t["score_fraction"].description
    assert "NOT a host probability" in t["S"].description
    assert "host_prob" not in t.colnames


def test_ngc4993_ranks_first(gw170817_plan):
    t = gw170817_plan.galaxies
    assert t["gal_id"][0] == NGC4993
    assert 0.3 < t["credible_level"][0] < 0.9


def test_host_posterior_normalisation(gw170817_map, glade_snapshot):
    from tessera import PlanConfig, plan

    r = plan(gw170817_map, glade_snapshot, [], PlanConfig(strategy="host-posterior"), "2017-08-17T22:00:00")
    t = r.galaxies
    t1 = np.asarray(t["tier"]) == 1
    total = float(np.nansum(np.asarray(t["host_prob"])[t1])) + r.budget["p_uncat_host_model"]
    assert total == pytest.approx(1.0, rel=1e-9)
    assert "STRATEGY_EXPERIMENTAL" in [d["code"] for d in r.diagnostics]


def test_crossmatch_agreement_with_ligo_skymap(s200105ae_map):
    """p_sky, credible level and dP/dV agree with ligo.skymap.postprocess.crossmatch (spec 11.4)."""
    pytest.importorskip("ligo.skymap")
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from ligo.skymap.io import read_sky_map
    from ligo.skymap.postprocess import crossmatch

    from ..conftest import S200105AE_MAP

    m = s200105ae_map
    rng = np.random.default_rng(7)
    # sample positions from the map itself so they fall in populated pixels
    idx = rng.choice(len(m), 300, p=m.prob / m.prob.sum())
    ra, dec = m.pixel_centers(idx)
    d = rng.uniform(150, 350, 300)
    tab = read_sky_map(str(S200105AE_MAP), moc=True)
    res = crossmatch(tab, SkyCoord(ra * u.deg, dec * u.deg, distance=d * u.Mpc))
    ours = m.lookup(ra, dec)
    np.testing.assert_allclose(m.probdensity[ours], res.probdensity, rtol=1e-12)
    np.testing.assert_allclose(m.credible_levels()[ours], res.searched_prob, rtol=1e-9, atol=1e-12)
    # crossmatch evaluates dP/dV on a 1000-bin radial grid (at the lower bin edge), whereas
    # tessera uses the exact closed form.  Reproduce its grid
    # node to check the formula itself, and bound the discretisation difference separately.
    from ligo.skymap import distance

    distmean, _ = distance.parameters_to_marginal_moments(m.prob, m.distmu, m.distsigma)
    d_r = 6 * distmean / 1000
    r = d_r * np.arange(1, 1000)
    r_i = r[np.digitize(d, r) - 1]
    exact_at_node = m.probdensity[ours] * m.distnorm[ours] * gauss(r_i, m.distmu[ours], m.distsigma[ours])
    np.testing.assert_allclose(exact_at_node, res.probdensity_vol, rtol=1e-9)
    exact = m.probdensity[ours] * m.distnorm[ours] * gauss(d, m.distmu[ours], m.distsigma[ours])
    assert np.median(np.abs(exact / res.probdensity_vol - 1)) < 0.05


def test_tabulated_luminosity_distance_matches_astropy():
    import numpy as np
    from astropy.cosmology import Planck18

    from tessera.core.association import _ddl_dz, _redshift_to_distance

    z = np.geomspace(1e-6, 2.99, 2000)
    assert np.max(np.abs(_redshift_to_distance(z) / Planck18.luminosity_distance(z).value - 1)) < 1e-9
    h = 1e-6
    num = (Planck18.luminosity_distance(z + h).value - Planck18.luminosity_distance(z - h / 2).value) / (1.5 * h)
    assert np.max(np.abs(_ddl_dz(z[z > 1e-3]) / num[z > 1e-3] - 1)) < 1e-5
    assert np.isnan(_redshift_to_distance(np.array([np.nan, -1.0]))).all()

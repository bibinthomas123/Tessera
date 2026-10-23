from __future__ import annotations

import gzip
import io

import numpy as np
import pytest
from astropy.io import fits

from tessera.core.errors import SkymapError
from tessera.core.skymap import SkymapLimits, from_arrays, load_skymap, to_fits_bytes
from tessera.testing import gaussian_map_arrays

from ..conftest import GW170817_MAP, S190425Z_MAP, S200105AE_MAP


def _flat_bytes(order: int, ordering: str, *, distance: bool = True) -> tuple[bytes, np.ndarray]:
    """Flat map written in RING or NESTED order from the same nested arrays."""
    from astropy_healpix.core import nested_to_ring

    uniq, rho, mu, sg, nm = gaussian_map_arrays([(40.0, 10.0, 3.0, 1.0)], order=order, dist=(100.0, 20.0))
    # expand to a full flat nested map at ``order``
    from tessera.core.skymap import _uniq_to_order_ipix

    o, ip = _uniq_to_order_ipix(uniq)
    npix = 12 * 4**order
    prob = np.zeros(npix)
    flat = {k: np.zeros(npix) for k in ("mu", "sg", "nm")}
    for k in range(len(uniq)):
        span = 4 ** (order - o[k])
        sl = slice(ip[k] * span, (ip[k] + 1) * span)
        prob[sl] = rho[k] * 4 * np.pi / npix
        flat["mu"][sl], flat["sg"][sl], flat["nm"][sl] = mu[k], sg[k], nm[k]
    arrays = {"PROB": prob, "DISTMU": flat["mu"], "DISTSIGMA": flat["sg"], "DISTNORM": flat["nm"]}
    if ordering == "RING":
        ring = np.asarray(nested_to_ring(np.arange(npix), 2**order))
        reordered = {}
        for name, v in arrays.items():
            out = np.empty_like(v)
            out[ring] = v
            reordered[name] = out
        arrays = reordered
    cols = [fits.Column(name=k, format="D", array=v) for k, v in arrays.items() if distance or k == "PROB"]
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.header.update({"PIXTYPE": "HEALPIX", "ORDERING": ordering, "NSIDE": 2**order, "INDXSCHM": "IMPLICIT",
                       "CREATOR": "BAYESTAR", "MJD-OBS": 58000.0})
    buf = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(buf)
    return buf.getvalue(), prob


def test_public_maps_load(gw170817_map, s200105ae_map):
    assert gw170817_map.meta.format == "flat" and gw170817_map.meta.nside == 1024
    assert s200105ae_map.meta.format == "multiorder" and s200105ae_map.meta.ordering == "NUNIQ"
    for m in (gw170817_map, s200105ae_map):
        assert abs(m.prob.sum() - 1) < 1e-3
        assert m.has_distance_columns
    # GWTC-1 GW170817: 90% area ~ 16 deg2 at NSIDE 1024
    assert 10 < gw170817_map.area_deg2(0.9) < 40
    assert gw170817_map.meta.sha256_raw != gw170817_map.meta.sha256_data  # gzip


def test_flat_ring_and_nested_agree():
    nest, _ = _flat_bytes(6, "NESTED")
    ring, _ = _flat_bytes(6, "RING")
    a, b = load_skymap(nest), load_skymap(ring)
    np.testing.assert_array_equal(a.uniq, b.uniq)
    np.testing.assert_allclose(a.probdensity, b.probdensity)
    np.testing.assert_allclose(a.distmu, b.distmu)
    ra, dec = np.array([40.0, 41.0, 200.0]), np.array([10.0, 12.0, -50.0])
    np.testing.assert_array_equal(a.lookup(ra, dec), b.lookup(ra, dec))


def test_multiorder_lookup_matches_flat_expansion():
    uniq, rho, mu, sg, nm = gaussian_map_arrays([(40.0, 10.0, 3.0, 1.0)], order=8)
    m = from_arrays(uniq, rho, mu, sg, nm)
    rng = np.random.default_rng(0)
    ra, dec = rng.uniform(0, 360, 2000), np.degrees(np.arcsin(rng.uniform(-1, 1, 2000)))
    idx = m.lookup(ra, dec)
    assert np.all(idx >= 0)
    # exact containment: re-index each point at its pixel's order
    import astropy.units as u
    from astropy_healpix import lonlat_to_healpix

    for o in np.unique(m.order[idx]):
        sel = m.order[idx] == o
        pix = lonlat_to_healpix(ra[sel] * u.deg, dec[sel] * u.deg, 2**int(o), order="nested")
        np.testing.assert_array_equal(np.asarray(pix), m.ipix[idx[sel]])


def test_credible_levels_bruteforce(s200105ae_map):
    m = s200105ae_map
    cl = m.credible_levels()
    order = np.argsort(-m.probdensity, kind="stable")
    np.testing.assert_allclose(np.sort(cl), np.sort(np.cumsum(m.prob[order])), rtol=1e-12)
    assert np.all((cl > 0) & (cl <= 1 + 1e-9))


def test_gzip_roundtrip(s200105ae_map):
    raw = S200105AE_MAP.read_bytes()
    m = load_skymap(gzip.compress(raw))
    np.testing.assert_array_equal(m.uniq, s200105ae_map.uniq)
    assert m.meta.sha256_data == s200105ae_map.meta.sha256_raw


def test_s190425z_flat_bayestar():
    m = load_skymap(S190425Z_MAP)
    assert m.meta.nside == 256 and m.meta.distance_prior == "r2"
    assert m.area_deg2(0.9) > 5000


@pytest.mark.parametrize("blob,code", [
    (b"", "SKYMAP_INVALID_FORMAT"),
    (b"not a fits file at all" * 10, "SKYMAP_INVALID_FORMAT"),
    (b"\x1f\x8b" + b"\x00" * 64, "SKYMAP_INVALID_FORMAT"),
])
def test_malformed_rejected(blob, code):
    with pytest.raises(SkymapError) as e:
        load_skymap(blob)
    assert e.value.code == code


def test_truncated_file_rejected():
    raw = S200105AE_MAP.read_bytes()
    with pytest.raises(SkymapError):
        load_skymap(raw[: len(raw) // 3])


def test_gzip_bomb_rejected():
    bomb = gzip.compress(b"\x00" * (40 * 1024**2))
    with pytest.raises(SkymapError) as e:
        load_skymap(bomb, limits=SkymapLimits(max_decompressed_bytes=8 * 1024**2))
    assert e.value.code == "SKYMAP_TOO_LARGE"


def test_upload_size_limit():
    with pytest.raises(SkymapError) as e:
        load_skymap(GW170817_MAP, limits=SkymapLimits(max_bytes=1000))
    assert e.value.code == "SKYMAP_TOO_LARGE"


def test_normalisation_rules():
    uniq, rho, mu, sg, nm = gaussian_map_arrays([(10.0, 10.0, 5.0, 1.0)], order=6)
    m = from_arrays(uniq, rho * 1.005, mu, sg, nm)
    assert "SKYMAP_RENORMALISED" in m.diagnostics.codes()
    assert abs(m.prob.sum() - 1) < 1e-12
    with pytest.raises(SkymapError) as e:
        from_arrays(uniq, rho * 1.5, mu, sg, nm)
    assert e.value.code == "SKYMAP_NOT_NORMALISED"


def test_invalid_values_rejected():
    uniq, rho, mu, sg, nm = gaussian_map_arrays([(10.0, 10.0, 5.0, 1.0)], order=6)
    bad = rho.copy()
    bad[3] = np.nan
    with pytest.raises(SkymapError):
        from_arrays(uniq, bad, mu, sg, nm)
    sg2 = sg.copy()
    sg2[5] = 0.0
    with pytest.raises(SkymapError):
        from_arrays(uniq, rho, mu, sg2, nm)


def test_partial_distance_columns_rejected():
    uniq, rho, mu, sg, nm = gaussian_map_arrays([(10.0, 10.0, 5.0, 1.0)], order=6)
    hdu = fits.BinTableHDU.from_columns([fits.Column(name="UNIQ", format="K", array=uniq),
                                         fits.Column(name="PROBDENSITY", format="D", array=rho),
                                         fits.Column(name="DISTMU", format="D", array=mu)])
    hdu.header.update({"PIXTYPE": "HEALPIX", "ORDERING": "NUNIQ"})
    buf = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(buf)
    with pytest.raises(SkymapError):
        load_skymap(buf.getvalue())


def test_overlapping_pixels_rejected():
    uniq, rho, mu, sg, nm = gaussian_map_arrays([(10.0, 10.0, 5.0, 1.0)], order=5, dist=None)
    # duplicate one pixel as its four children too
    o = 5
    first = uniq[0] - 4 * 4**o
    children = 4 * 4 ** (o + 1) + first * 4 + np.arange(4)
    u2 = np.concatenate([uniq, children])
    r2 = np.concatenate([rho, np.zeros(4)])
    with pytest.raises(SkymapError):
        load_skymap(to_fits_bytes(u2, r2))


def test_no_distance_map_classified(rng):
    from tessera.testing import gaussian_map

    m = gaussian_map([(10, 10, 2, 1)], order=7, dist=None)
    assert m.distance_class() == "none"
    assert "EVENT_NO_DIST" in m.diagnostics.codes()


def test_non_bayestar_prior_flagged(gw170817_map):
    # GWTC-1 PE map: CREATOR is ligo-skymap-from-samples
    assert "DIST_PRIOR_ASSUMED" in gw170817_map.diagnostics.codes()

"""Synthetic sky maps, catalogs and injections (validation level 1, spec 11.1).

Sky maps are analytic multi-order maps with the LVK distance ansatz; they are
written and re-read through the normal FITS path so the same validation
applies.  Galaxy universes are drawn uniformly in volume with a log-normal
K-band luminosity function whose mean luminosity density equals the
science-model reference density (so completeness estimates can be checked
against a known truth).
"""

from __future__ import annotations

from dataclasses import dataclass

import astropy.units as u
import numpy as np
from astropy_healpix import healpix_to_lonlat
from scipy.special import ndtr

from tessera.core import science
from tessera.core.geometry import angsep
from tessera.core.skymap import MultiOrderMap, from_arrays, to_fits_bytes


def distnorm(mu, sigma):
    """DISTNORM that normalises r^2 Gauss(r; mu, sigma) over r > 0."""
    mu, sigma = np.asarray(mu, float), np.asarray(sigma, float)
    a = mu / sigma
    phi = np.exp(-0.5 * a * a) / np.sqrt(2 * np.pi)
    return 1.0 / ((mu**2 + sigma**2) * ndtr(a) + mu * sigma * phi)


def _adaptive_pixels(components, order: int, coarse: int = 5):
    """Multi-order pixels: ``coarse`` everywhere, refined to ``order`` within 6 sigma of each component."""
    co = min(order, coarse)
    nside = 2**co
    ipix = np.arange(12 * nside * nside, dtype=np.int64)
    lon, lat = healpix_to_lonlat(ipix, nside, order="nested")
    ra, dec = lon.to_value(u.deg), lat.to_value(u.deg)
    pix_deg = np.degrees(np.sqrt(4 * np.pi / len(ipix)))
    refine = np.zeros(len(ipix), bool)
    if order > co:
        for (ra0, dec0, sig, _w) in components:
            refine |= angsep(ra0, dec0, ra, dec) <= 6 * sig + 1.5 * pix_deg
    orders = [np.full(int((~refine).sum()), co)]
    ips = [ipix[~refine]]
    if refine.any():
        k = 4 ** (order - co)
        child = (ipix[refine][:, None] * k + np.arange(k)[None, :]).ravel()
        orders.append(np.full(len(child), order))
        ips.append(child)
    o = np.concatenate(orders)
    ip = np.concatenate(ips)
    ra_all = np.empty(len(ip))
    dec_all = np.empty(len(ip))
    for oo in np.unique(o):
        m = o == oo
        lo, la = healpix_to_lonlat(ip[m], 2**int(oo), order="nested")
        ra_all[m], dec_all[m] = lo.to_value(u.deg), la.to_value(u.deg)
    return o, ip, ra_all, dec_all


def gaussian_map_arrays(components: list[tuple[float, float, float, float]], *, order: int = 8,
                        dist: tuple[float, float] | None = (40.0, 8.0), dist_gradient: float = 0.0):
    """Mixture of 2D Gaussians [(ra, dec, sigma_deg, weight)] with distance ansatz (mu, sigma).

    Returns multi-order arrays (UNIQ, PROBDENSITY, DISTMU, DISTSIGMA, DISTNORM)."""
    o, ip, ra, dec = _adaptive_pixels(components, order)
    rho = np.zeros(len(ip))
    for (ra0, dec0, sig, w) in components:
        th = angsep(ra0, dec0, ra, dec)
        rho += w * np.exp(-0.5 * (th / sig) ** 2) / (2 * np.pi * np.radians(sig) ** 2)
    area = 4 * np.pi / (12 * 4.0**o)
    rho /= np.sum(rho * area)
    uniq = 4 * 4**o + ip
    if dist is None:
        return uniq, rho, None, None, None
    mu = np.full(len(ip), dist[0]) + dist_gradient * (dec - components[0][1])
    sg = np.full(len(ip), dist[1])
    return uniq, rho, mu, sg, distnorm(mu, sg)


def gaussian_map(components, *, order: int = 8, dist=(40.0, 8.0), dist_gradient: float = 0.0,
                 mjd_obs: float = 57982.528523, creator: str = "BAYESTAR") -> MultiOrderMap:
    uniq, rho, mu, sg, nm = gaussian_map_arrays(components, order=order, dist=dist, dist_gradient=dist_gradient)
    return from_arrays(uniq, rho, mu, sg, nm, header={"OBJECT": "SYNTH", "MJD-OBS": mjd_obs, "CREATOR": creator,
                                                     "DATE-OBS": "2017-08-17T12:41:04"})


def gaussian_map_bytes(components, **kw) -> bytes:
    uniq, rho, mu, sg, nm = gaussian_map_arrays(components, order=kw.get("order", 8), dist=kw.get("dist", (40.0, 8.0)))
    return to_fits_bytes(uniq, rho, mu, sg, nm, header={"MJD-OBS": kw.get("mjd_obs", 57982.528523),
                                                        "CREATOR": "BAYESTAR", "OBJECT": "SYNTH"})


def edge_case_maps() -> dict[str, MultiOrderMap]:
    return {
        "bimodal": gaussian_map([(30, 10, 3, 1.0), (210, -10, 3, 0.7)], order=7),
        "four_modes": gaussian_map([(30, 10, 2, 1), (120, -40, 2, 1), (210, 30, 2, 1), (300, -5, 2, 1)], order=7),
        "ra_zero_wrap": gaussian_map([(0.2, -20, 2.0, 1.0)], order=8),
        "south_pole": gaussian_map([(0, -89.5, 3.0, 1.0)], order=7),
        "tiny_1deg2": gaussian_map([(150, 5, 0.3, 1.0)], order=10),
        "allsky_single_detector": gaussian_map([(100, 0, 60, 1.0)], order=5, dist=(150, 60)),
        "no_distance": gaussian_map([(197.45, -23.38, 2.0, 1.0)], order=8, dist=None),
        "beyond_catalog": gaussian_map([(197.45, -23.38, 2.0, 1.0)], order=8, dist=(2500, 400)),
    }


@dataclass
class Universe:
    cols: dict[str, np.ndarray]
    L: np.ndarray
    d_true: np.ndarray


def synthetic_universe(rng: np.random.Generator, *, ra0: float, dec0: float, radius_deg: float, d_max: float,
                       mean_L: float = 1e10, sigma_lnL: float = 1.0, n_scale: float = 1.0) -> Universe:
    """Galaxies uniform in volume inside a cone; luminosity density = rho_bar(lum_K) * n_scale."""
    rho_bar = float(science.reference_density("lum_K")["value"])
    omega = 2 * np.pi * (1 - np.cos(np.radians(radius_deg)))
    vol = omega / 3.0 * d_max**3
    # log-normal with E[L] = mean_L
    n = rng.poisson(n_scale * rho_bar * vol / mean_L)
    cos_t = rng.uniform(np.cos(np.radians(radius_deg)), 1.0, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    d = d_max * rng.uniform(0, 1, n) ** (1 / 3)
    # rotate cone around (ra0, dec0)
    st = np.sqrt(1 - cos_t**2)
    x, y, z = st * np.cos(phi), st * np.sin(phi), cos_t
    a, b = np.radians(ra0), np.radians(90 - dec0)
    # rotate about y by b then about z by a
    x1 = x * np.cos(b) + z * np.sin(b)
    z1 = -x * np.sin(b) + z * np.cos(b)
    x2 = x1 * np.cos(a) - y * np.sin(a)
    y2 = x1 * np.sin(a) + y * np.cos(a)
    ra = np.mod(np.degrees(np.arctan2(y2, x2)), 360)
    dec = np.degrees(np.arcsin(np.clip(z1, -1, 1)))
    mu_ln = np.log(mean_L) - 0.5 * sigma_lnL**2
    L = np.exp(rng.normal(mu_ln, sigma_lnL, n))
    msun = science.defaults()["solar_abs_mag"]["K"]
    mag_K = msun - 2.5 * np.log10(L) + 5 * np.log10(np.maximum(d, 1e-3)) + 25
    derr = 0.05 * d + 1.0
    cols = {
        "gal_id": np.array([f"SYN:{i:07d}" for i in range(n)], dtype=object),
        "ra": ra, "dec": dec,
        "dist_mpc": np.maximum(0.5, d + rng.normal(0, 1, n) * derr),
        "dist_err_mpc": derr,
        "dist_kind": np.full(n, "z_indep", dtype=object),
        "mag_K": mag_K,
    }
    return Universe(cols, L, d)


def injection_map(rng: np.random.Generator, ra: float, dec: float, d: float, *, sigma_deg: float = 2.0,
                  frac_dist_err: float = 0.2, order: int = 8, mjd_obs: float = 57982.528523) -> tuple[MultiOrderMap, tuple]:
    """GW posterior for a source at (ra, dec, d): centred on a noisy measurement, as a real localisation is."""
    off = rng.normal(0, sigma_deg / np.sqrt(2), 2)
    cra = ra + off[0] / max(np.cos(np.radians(dec)), 0.05)
    cdec = float(np.clip(dec + off[1], -89.9, 89.9))
    ds = frac_dist_err * d
    dmu = max(1.0, d + rng.normal(0, ds))
    return gaussian_map([(cra, cdec, sigma_deg, 1.0)], order=order, dist=(dmu, ds), mjd_obs=mjd_obs), (cra, cdec, dmu)

"""Closed-form expressions for the LVK distance ansatz (Singer et al. 2016).

Along each line of sight the conditional distance posterior is

    p(r | n) = N r^2 / (sqrt(2 pi) sigma) exp(-(r - mu)^2 / (2 sigma^2))

with (mu, sigma, N) = (DISTMU, DISTSIGMA, DISTNORM).  Pixels without usable
distance information have mu = +inf and N = 0 (ligo.skymap convention).
"""

from __future__ import annotations

import numpy as np
from scipy.special import ndtr

_SQRT2PI = np.sqrt(2.0 * np.pi)


def valid_distance(mu: np.ndarray, sigma: np.ndarray, norm: np.ndarray) -> np.ndarray:
    return np.isfinite(mu) & np.isfinite(sigma) & (sigma > 0) & np.isfinite(norm) & (norm > 0)


def gauss(x, mu, sigma) -> np.ndarray:
    x, mu, sigma = np.broadcast_arrays(np.asarray(x, float), np.asarray(mu, float), np.asarray(sigma, float))
    out = np.zeros(x.shape)
    ok = np.isfinite(mu) & np.isfinite(sigma) & (sigma > 0) & np.isfinite(x)
    z = (x[ok] - mu[ok]) / sigma[ok]
    out[ok] = np.exp(-0.5 * z * z) / (_SQRT2PI * sigma[ok])
    return out


def conditional_pdf(r, mu, sigma, norm) -> np.ndarray:
    """p(r | n) in Mpc^-1."""
    r = np.asarray(r, float)
    return np.where(valid_distance(np.asarray(mu), np.asarray(sigma), np.asarray(norm)),
                    np.asarray(norm) * r * r * gauss(r, mu, sigma), 0.0)


def volume_density_factor(r, mu, sigma, norm) -> np.ndarray:
    """N * Gauss(r; mu, sigma) [Mpc^-3]: dP/dV divided by the sky density."""
    return np.where(valid_distance(np.asarray(mu), np.asarray(sigma), np.asarray(norm)),
                    np.asarray(norm) * gauss(r, mu, sigma), 0.0)


def shell_probability(mu, sigma, norm, r_lo, r_hi) -> np.ndarray:
    """Integral of p(r | n) over [r_lo, r_hi] (closed form).

    With t = (r - mu)/sigma:
      int r^2 phi dt = (mu^2 + sigma^2)[Phi(b) - Phi(a)] + 2 mu sigma [phi(a) - phi(b)]
                       + sigma^2 [a phi(a) - b phi(b)]
    Broadcasts ``mu/sigma/norm`` (pixels) against ``r_lo/r_hi`` (shells).
    """
    mu = np.asarray(mu, float)
    sigma = np.asarray(sigma, float)
    norm = np.asarray(norm, float)
    ok = valid_distance(mu, sigma, norm)
    mu_s = np.where(ok, mu, 0.0)
    sg_s = np.where(ok, sigma, 1.0)
    a = (np.asarray(r_lo, float) - mu_s) / sg_s
    b = (np.asarray(r_hi, float) - mu_s) / sg_s
    phi_a = np.exp(-0.5 * a * a) / _SQRT2PI
    phi_b = np.exp(-0.5 * b * b) / _SQRT2PI
    # Guard a*phi(a) at infinities.
    aphi_a = np.where(np.isfinite(a), a * phi_a, 0.0)
    bphi_b = np.where(np.isfinite(b), b * phi_b, 0.0)
    val = (mu_s**2 + sg_s**2) * (ndtr(b) - ndtr(a)) + 2 * mu_s * sg_s * (phi_a - phi_b) + sg_s**2 * (aphi_a - bphi_b)
    return np.where(ok, norm * np.clip(val, 0.0, None), 0.0)

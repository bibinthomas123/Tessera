"""Detectability model: limiting magnitude, light curves, P_det (spec 6.1).

    m_lim = m_ref + 1.25 log10(t_exp / t_ref) - k (X - 1) - dm_full * illum * [Moon up]
    P_det = Phi((m_lim - [M(t) + DM(r) + A]) / sigma_M)

The distance-modulus distribution in each sky pixel is summarised by its
mean and standard deviation (computed by quadrature over p(r | n)); the
detection efficiency is then Phi((Delta - E[DM]) / sqrt(Var[DM] + sigma_M^2))
with Delta = m_lim - M - A.  This Gaussian-in-DM approximation is documented
in docs/science.md and is accurate to a few per cent for sigma/mu < 0.5.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.special import ndtr

from tessera.core import science
from tessera.core.distance import conditional_pdf, valid_distance
from tessera.core.profile import Exposure, FilterSpec

# A_band / E(B-V) for common optical filters (Schlafly & Finkbeiner 2011, R_V=3.1; LSST/SDSS-like).
R_FILTER = {"u": 4.145, "g": 3.237, "r": 2.273, "i": 1.684, "z": 1.323, "y": 1.088,
            "B": 3.626, "V": 2.742, "R": 2.169, "I": 1.505, "J": 0.709, "H": 0.449, "K": 0.302, "L": 2.5,
            "clear": 2.5, "w": 2.5}


def distance_modulus(r_mpc) -> np.ndarray:
    return 5.0 * np.log10(np.asarray(r_mpc, float)) + 25.0


@dataclass(frozen=True)
class LightCurve:
    """Absolute-magnitude light curve: parametric (peak + linear decline) or tabulated."""

    name: str
    t_peak_days: float
    sigma_mag: float
    filters: dict[str, dict[str, float]]
    provisional: bool
    table: dict[str, tuple[np.ndarray, np.ndarray]] | None = None
    filter_map: dict[str, str] | None = None

    @classmethod
    def named(cls, name: str) -> LightCurve:
        return _load_light_curve(name)

    def abs_mag(self, t_days, filt: str) -> np.ndarray:
        t = np.asarray(t_days, float)
        if self.table is not None:
            band = filt if filt in self.table else (self.filter_map or {}).get(filt, "r")
            tt, mm = self.table.get(band, self.table["r"])
            lt = np.log10(np.maximum(t, 1e-3))
            m = np.interp(lt, np.log10(tt), mm)              # holds first value before tt[0]
            if len(tt) >= 2:                                  # extrapolate the last decline
                slope = max(0.0, (mm[-1] - mm[-2]) / (tt[-1] - tt[-2]))
                m = np.where(t > tt[-1], mm[-1] + slope * (t - tt[-1]), m)
            return np.where(t > 0, m, np.inf)
        p = self.filters.get(filt, self.filters["default"])
        m = p["m_peak"] + p["decline"] * np.maximum(0.0, t - self.t_peak_days)
        return np.where(t > 0, m, np.inf)


@lru_cache(maxsize=16)
def _load_light_curve(name: str) -> LightCurve:
    lc = science.light_curve(name)
    table = None
    if lc.get("kind") == "table":
        import csv
        from importlib import resources

        rows: dict[str, list[tuple[float, float]]] = {}
        text = resources.files("tessera.data.science").joinpath(lc["file"]).read_text()
        for r in csv.DictReader(line for line in text.splitlines() if not line.startswith("#")):
            rows.setdefault(r["filter"], []).append((float(r["t_days"]), float(r["abs_mag"])))
        table = {}
        for b, pts in rows.items():
            pts.sort()
            table[b] = (np.array([p[0] for p in pts]), np.array([p[1] for p in pts]))
    return LightCurve(name=name, t_peak_days=float(lc.get("t_peak_days", 0.0)), sigma_mag=float(lc["sigma_mag"]),
                      filters=lc.get("filters", {}), provisional=bool(lc.get("provisional", False)), table=table,
                      filter_map=lc.get("filter_map"))


def m_lim(spec: FilterSpec, t_exp_s, airmass, moon_illum=0.0, moon_up=False) -> np.ndarray:
    t = np.asarray(t_exp_s, float)
    X = np.asarray(airmass, float)
    moon = spec.moon.dm_full * np.asarray(moon_illum, float) * np.asarray(moon_up, float)
    return spec.m_ref + 1.25 * np.log10(t / spec.t_ref_s) - spec.k_ext * (X - 1.0) - moon


def solve_exposure(spec: FilterSpec, exp: Exposure, m_needed: float, airmass: float, moon_illum: float,
                   moon_up: bool) -> tuple[int, float]:
    """(n_exp, t_single_s) reaching ``m_needed`` (clamped to [t_min, t_max])."""
    if exp.fixed_s is not None:
        total = exp.fixed_s * (exp.n_exp_fixed or 1)
    elif not np.isfinite(m_needed):
        total = exp.t_min_s
    else:
        loss = spec.k_ext * (airmass - 1.0) + spec.moon.dm_full * moon_illum * float(moon_up)
        total = spec.t_ref_s * 10 ** ((m_needed - spec.m_ref + loss) / 1.25)
        total = float(np.clip(total, exp.t_min_s, exp.t_max_s))
    if exp.n_exp_fixed:
        n = exp.n_exp_fixed
    else:
        n = int(np.ceil(total / exp.t_single_max_s - 1e-9))
    n = max(n, 1)
    return n, round(total / n, 1)


def pixel_dm_stats(mu, sigma, norm, n_nodes: int = 64) -> tuple[np.ndarray, np.ndarray]:
    """Mean and std of the distance modulus under p(r | n) for each pixel (NaN where invalid)."""
    mu = np.asarray(mu, float)
    sigma = np.asarray(sigma, float)
    norm = np.asarray(norm, float)
    ok = valid_distance(mu, sigma, norm)
    mean = np.full(mu.shape, np.nan)
    sd = np.full(mu.shape, np.nan)
    if not ok.any():
        return mean, sd
    m, s, nn = mu[ok], sigma[ok], norm[ok]
    # integrate on a grid spanning the bulk of r^2 Gauss(r; mu, sigma), r > 0
    hi = np.maximum(m, 0) + 6 * s
    lo = np.maximum(1e-3, np.maximum(m, 0) - 6 * s)
    t = np.linspace(0, 1, n_nodes)
    r = lo[:, None] + (hi - lo)[:, None] * t[None, :]
    p = conditional_pdf(r, m[:, None], s[:, None], nn[:, None])
    wsum = p.sum(axis=1)
    good = wsum > 0
    dm = distance_modulus(r)
    e1 = np.where(good, (p * dm).sum(axis=1) / np.where(good, wsum, 1), np.nan)
    e2 = np.where(good, (p * dm**2).sum(axis=1) / np.where(good, wsum, 1), np.nan)
    mean[ok] = e1
    sd[ok] = np.sqrt(np.maximum(e2 - e1**2, 0))
    return mean, sd


def efficiency(delta, dm_mean, dm_sd, sigma_M: float) -> np.ndarray:
    """P_det averaged over the pixel's distance distribution.  delta = m_lim - M - A."""
    dm_mean = np.asarray(dm_mean, float)
    tot = np.sqrt(np.nan_to_num(np.asarray(dm_sd, float)) ** 2 + sigma_M**2)
    out = ndtr((np.asarray(delta, float) - dm_mean) / tot)
    return np.where(np.isfinite(dm_mean), out, np.nan)

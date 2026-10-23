"""Astrophysical host weights w_i (spec 3.4).

One implementation serves both galaxy scoring and the completeness cube, so
the two can never disagree about what "luminosity" means.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from tessera.core import science
from tessera.core.errors import ConfigError

# Extinction coefficients A_band / E(B-V) (Schlafly & Finkbeiner 2011, R_V = 3.1).
R_BAND = {"B": 3.626, "K": 0.310, "W1": 0.18}

WEIGHT_KINDS = ("flat", "lum_B", "lum_K", "mstar", "bns_rate", "sfr")


@dataclass
class WeightResult:
    w: np.ndarray            # final weight (imputed where needed, ** alpha applied)
    raw: np.ndarray          # tracer value before imputation (NaN where missing)
    imputed: np.ndarray      # bool
    source: np.ndarray       # per-galaxy tracer used ("K", "W1", "B", "mstar", ...)
    impute_value: float | None


def _lum(mag: np.ndarray, dist_mpc: np.ndarray, ebv: np.ndarray, band: str) -> np.ndarray:
    msun = science.defaults()["solar_abs_mag"][band]
    with np.errstate(invalid="ignore", divide="ignore"):
        absmag = mag - 5.0 * np.log10(dist_mpc) - 25.0 - R_BAND[band] * ebv
        return 10.0 ** (-0.4 * (absmag - msun))


def tracer(cols: dict[str, np.ndarray], kind: str, dist_mpc: np.ndarray, ebv: np.ndarray | None = None
           ) -> tuple[np.ndarray, np.ndarray]:
    """Raw tracer value (solar units) and the band/column used, NaN where unavailable."""
    n = len(dist_mpc)
    ebv = np.zeros(n) if ebv is None else np.nan_to_num(np.asarray(ebv, float), nan=0.0)
    src = np.full(n, "", dtype=object)
    if kind == "flat":
        return np.ones(n), np.full(n, "flat", dtype=object)
    if kind == "lum_B":
        v = _lum(np.asarray(cols["mag_B"], float), dist_mpc, ebv, "B")
        src[np.isfinite(v)] = "B"
        return v, src
    if kind == "lum_K":
        v = _lum(np.asarray(cols["mag_K"], float), dist_mpc, ebv, "K")
        src[np.isfinite(v)] = "K"
        if "mag_W1" in cols:
            v1 = _lum(np.asarray(cols["mag_W1"], float), dist_mpc, ebv, "W1")
            use = ~np.isfinite(v) & np.isfinite(v1)
            v = np.where(use, v1, v)
            src[use] = "W1"
        return v, src
    if kind == "mstar":
        v = 10.0 ** np.asarray(cols["log_mstar"], float)
        src[np.isfinite(v)] = "mstar"
        return v, src
    if kind == "bns_rate":
        v = 10.0 ** np.asarray(cols["log_bns_rate"], float)
        src[np.isfinite(v)] = "bns_rate"
        return v, src
    if kind == "sfr":
        if "sfr" not in cols:
            raise ConfigError("weight 'sfr' requested but the catalog has no 'sfr' column", code="WEIGHT_UNAVAILABLE")
        v = np.asarray(cols["sfr"], float)
        src[np.isfinite(v)] = "sfr"
        return v, src
    raise ConfigError(f"unknown weight kind {kind!r}", details={"known": list(WEIGHT_KINDS)})


def compute_weights(cols: dict[str, np.ndarray], kind: str, dist_mpc: np.ndarray, *, ebv: np.ndarray | None = None,
                    alpha: float = 1.0, impute_percentile: float = 5.0) -> WeightResult:
    raw, src = tracer(cols, kind, dist_mpc, ebv)
    raw = np.where(np.isfinite(raw) & (raw > 0), raw, np.nan)
    missing = ~np.isfinite(raw)
    impute_value: float | None = None
    w = raw.copy()
    if missing.any():
        good = raw[~missing]
        impute_value = float(np.percentile(good, impute_percentile, method="lower")) if good.size else 1.0
        w[missing] = impute_value
    if alpha != 1.0:
        w = w**alpha
    return WeightResult(w=w, raw=raw, imputed=missing, source=src, impute_value=impute_value)

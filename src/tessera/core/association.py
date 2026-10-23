"""Galaxy association (science spec sections 3.2-3.8).

For every galaxy i in the candidate region::

    p_sky,i  = rho_j(i)                     [deg^-2], j(i) = containing pixel
    p_dist,i = N_j Gauss(d_i; mu_j, sqrt(sigma_j^2 + s_i^2))   [Mpc^-3]
    w_i      = tracer_i ** alpha
    S_i      = p_sky,i * p_dist,i * w_i     (a ranking score, not a probability)

``p_dist`` is the GW posterior density per unit volume at the galaxy,
marginalised over a Gaussian galaxy-distance uncertainty s_i.  It is valid for
maps whose distance prior is proportional to r^2 (BAYESTAR).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from astropy.cosmology import Planck18
from astropy.table import Table

from tessera.core.catalog import Galaxies
from tessera.core.config import ResolvedConfig
from tessera.core.distance import gauss, valid_distance
from tessera.core.errors import Diagnostics
from tessera.core.geometry import SR_PER_DEG2
from tessera.core.skymap import MultiOrderMap
from tessera.core.weights import compute_weights

C_KMS = 299792.458


@dataclass
class AssociationResult:
    table: Table                      # ranked galaxies (tier 1 then tier 2)
    summary: dict[str, Any] = field(default_factory=dict)

    @property
    def tier1(self) -> np.ndarray:
        return np.asarray(self.table["tier"]) == 1


def _redshift_to_distance(z: np.ndarray) -> np.ndarray:
    d = np.full(z.shape, np.nan)
    ok = np.isfinite(z) & (z > 0)
    if ok.any():
        d[ok] = Planck18.luminosity_distance(z[ok]).value
    return d


def _ddl_dz(z: np.ndarray) -> np.ndarray:
    """Numerical derivative of the luminosity distance with respect to redshift."""
    h = 1e-5
    zp = z + h
    zm = np.maximum(z - h, 1e-9)
    return (_redshift_to_distance(zp) - _redshift_to_distance(zm)) / (zp - zm)


def resolve_distances(g: Galaxies, cfg: ResolvedConfig) -> tuple[np.ndarray, np.ndarray, list[set[str]]]:
    """Galaxy distance d_i and total uncertainty s_i (spec 3.3) plus per-galaxy flags."""
    n = len(g)
    flags: list[set[str]] = [set() for _ in range(n)]
    kind = np.asarray(g["dist_kind"]).astype(str)
    d = np.asarray(g["dist_mpc"], float).copy()
    z = np.asarray(g["z"], float)
    from_z = ~np.isfinite(d) & np.isfinite(z) & (z > 0)
    if from_z.any():
        d[from_z] = _redshift_to_distance(z[from_z])
        kind = np.where(from_z & (kind == "none"), "specz", kind)
        for i in np.flatnonzero(from_z):
            flags[i].add("DIST_FROM_Z")

    s_meas = np.asarray(g["dist_err_mpc"], float).copy()
    s_meas[~(s_meas > 0)] = np.nan
    need = ~np.isfinite(s_meas) & np.isfinite(d)
    zerr = np.asarray(g["z_err"], float)
    photoz = kind == "photoz"
    zerr = np.where(photoz & ~(zerr > 0), cfg.photoz_sigma_frac * (1 + np.nan_to_num(z)), zerr)
    use_z = need & (zerr > 0) & np.isfinite(z)
    if use_z.any():
        s_meas[use_z] = _ddl_dz(z[use_z]) * zerr[use_z]
    still = ~np.isfinite(s_meas) & np.isfinite(d)
    s_meas[still] = cfg.dist_frac_floor * d[still]

    h0 = Planck18.H0.value
    s_pv = cfg.sigma_v_kms / h0
    includes = np.asarray(g["dist_err_includes_pv"], bool)
    z_derived = np.isin(kind, ("specz", "photoz"))
    add_pv = z_derived & ~includes
    s = np.sqrt(np.nan_to_num(s_meas) ** 2 + np.where(add_pv, s_pv**2, 0.0))
    s[~np.isfinite(d)] = np.nan
    lowz = z_derived & np.isfinite(z) & (z < 0.01)
    for i in np.flatnonzero(photoz):
        flags[i].add("PHOTOZ")
    for i in np.flatnonzero(lowz):
        flags[i].add("LOW_Z_PV")
    return d, s, flags


def score_galaxies(
    skymap: MultiOrderMap,
    gals: Galaxies,
    cfg: ResolvedConfig,
    *,
    ebv: np.ndarray | None = None,
    diagnostics: Diagnostics | None = None,
) -> AssociationResult:
    diags = diagnostics if diagnostics is not None else Diagnostics()
    n_in = len(gals)
    summary: dict[str, Any] = {"n_candidates": n_in, "strategy": cfg.strategy, "weight": cfg.weight,
                               "weight_alpha": cfg.weight_alpha}
    if n_in == 0:
        diags.add("NO_GALAXIES_IN_REGION", "the catalog has no galaxies inside the candidate region",
                  cr_max=cfg.cr_max)
        return AssociationResult(_empty_table(), {**summary, "n_scored": 0})

    keep = np.ones(n_in, bool)
    dup = np.asarray(gals["duplicate_of"]).astype(str) != ""
    qso = np.asarray(gals["obj_type"]).astype(str) == "qso"
    summary["n_duplicates_skipped"] = int(dup.sum())
    keep &= ~dup
    if cfg.exclude_qso:
        summary["n_qso_excluded"] = int((qso & keep).sum())
        keep &= ~qso
    g = gals.take(np.flatnonzero(keep))
    ebv_g = None if ebv is None else np.asarray(ebv)[keep]
    n = len(g)

    idx = skymap.lookup(g["ra"], g["dec"])
    inside = idx >= 0
    rho = np.where(inside, skymap.probdensity[np.clip(idx, 0, None)], 0.0)
    cl = np.where(inside, skymap.credible_levels()[np.clip(idx, 0, None)], 1.0)
    p_sky = rho * SR_PER_DEG2

    d, s, flags = resolve_distances(g, cfg)
    has_d = np.isfinite(d) & (d > 0)

    event_dist = skymap.distance_class()
    use_distance = cfg.uses_distance and event_dist != "none"
    summary["distance_used"] = use_distance
    if cfg.uses_distance and event_dist == "none" and "EVENT_NO_DIST" not in diags.codes():
        diags.add("EVENT_NO_DIST", "map has no usable distance; p_dist set to 1 (2D scoring)")

    p_dist = np.ones(n)
    zscore = np.full(n, np.nan)
    if use_distance:
        mu = np.where(inside, skymap.distmu[np.clip(idx, 0, None)], np.inf)
        sg = np.where(inside, skymap.distsigma[np.clip(idx, 0, None)], 1.0)
        nm = np.where(inside, skymap.distnorm[np.clip(idx, 0, None)], 0.0)
        pix_ok = valid_distance(mu, sg, nm)
        tot = np.sqrt(sg**2 + np.nan_to_num(s) ** 2)
        p_dist = np.where(pix_ok & has_d, nm * gauss(d, mu, tot), 0.0)
        with np.errstate(invalid="ignore"):
            zscore = np.where(pix_ok & has_d, (d - mu) / tot, np.nan)
        for i in np.flatnonzero(~pix_ok & has_d):
            flags[i].add("PIX_NO_DIST")
        p_dist[~has_d] = np.nan

    # tiers: distance-less galaxies never interleave with tier 1 (spec 3.7)
    tier = np.where(has_d | (not use_distance), 1, 2)
    for i in np.flatnonzero(~has_d):
        flags[i].add("NO_DIST")
    if cfg.missing_distance == "exclude":
        tier = np.where(tier == 2, 0, tier)

    wres = compute_weights(g.cols, cfg.weight, np.where(has_d, d, np.nan), ebv=ebv_g, alpha=cfg.weight_alpha,
                           impute_percentile=cfg.w_impute_percentile)
    w = wres.w
    if cfg.weight != "flat":
        t1 = tier == 1
        tier1_raw = wres.raw[t1]
        if t1.any() and np.isfinite(tier1_raw).any():
            iv = float(np.percentile(tier1_raw[np.isfinite(tier1_raw)], cfg.w_impute_percentile,
                                    method="lower")) ** cfg.weight_alpha
        else:
            iv = 1.0
        w = np.where(wres.imputed, iv, w)
        for i in np.flatnonzero(wres.imputed & t1):
            flags[i].add("W_IMPUTED")
        for i in np.flatnonzero(wres.source == "W1"):
            flags[i].add("W1_FOR_K")
        summary["w_impute_value"] = iv
    w = np.where(tier == 2, 1.0, w)

    S = np.where(tier == 1, p_sky * np.nan_to_num(p_dist) * w, np.where(tier == 2, p_sky, 0.0))
    for i in np.flatnonzero(cl > 0.9):
        flags[i].add("OUTSIDE_CR90")

    # pre-cut (efficiency only; its score mass is reported)
    precut = (cl > cfg.cr_max) | (np.isfinite(zscore) & (np.abs(zscore) > cfg.zscore_max)) | ~inside
    s1_total = float(S[tier == 1].sum())
    summary["p_precut"] = float(S[(tier == 1) & precut].sum() / s1_total) if s1_total > 0 else 0.0
    summary["n_precut"] = int(precut.sum())
    final = ~precut & (tier > 0)
    if cfg.missing_distance == "exclude":
        summary["n_no_dist_excluded"] = int((~has_d).sum())

    sel = np.flatnonzero(final)
    gid = np.asarray(g["gal_id"]).astype(str)
    order = sel[np.lexsort((gid[sel], -S[sel], tier[sel]))]
    t1sum = float(S[order][tier[order] == 1].sum())
    t2sum = float(S[order][tier[order] == 2].sum())
    frac = np.where(tier[order] == 1, S[order] / t1sum if t1sum > 0 else np.nan,
                    S[order] / t2sum if t2sum > 0 else np.nan)
    with np.errstate(divide="ignore"):
        log10_S = np.where(S[order] > 0, np.log10(np.where(S[order] > 0, S[order], 1.0)), -np.inf)

    tbl = Table()
    tbl["rank"] = np.arange(1, len(order) + 1, dtype=np.int64)
    tbl["tier"] = tier[order].astype(np.int16)
    tbl["gal_id"] = gid[order]
    tbl["alt_ids"] = np.asarray(g["alt_ids"]).astype(str)[order]
    tbl["ra"] = np.asarray(g["ra"], float)[order]
    tbl["dec"] = np.asarray(g["dec"], float)[order]
    tbl["dist_mpc"] = d[order]
    tbl["dist_err_mpc"] = s[order]
    tbl["dist_kind"] = np.asarray(g["dist_kind"]).astype(str)[order]
    tbl["z"] = np.asarray(g["z"], float)[order]
    tbl["p_sky_deg2"] = p_sky[order]
    tbl["credible_level"] = cl[order]
    tbl["p_dist_mpc3"] = p_dist[order]
    tbl["dist_zscore"] = zscore[order]
    tbl["w"] = w[order]
    tbl["w_kind"] = np.full(len(order), cfg.weight if cfg.weight else "flat")
    tbl["S"] = S[order]
    tbl["log10_S"] = log10_S
    tbl["score_fraction"] = frac
    tbl["flags"] = np.array(["|".join(sorted(flags[i])) for i in order], dtype=str) if len(order) else np.array([], str)
    _describe(tbl)
    summary.update({
        "n_scored": int(len(order)),
        "n_tier1": int((tier[order] == 1).sum()),
        "n_tier2": int((tier[order] == 2).sum()),
        "sum_S_tier1": t1sum,
        "sum_S_tier1_sr": t1sum / SR_PER_DEG2,
    })
    if len(order) and (tier[order] == 1).sum() == 0:
        diags.add("NO_TIER1_GALAXIES", "no galaxy with a distance survived the cuts")
    return AssociationResult(tbl, summary)


def add_host_probability(res: AssociationResult, U_deg: float) -> None:
    """Attach ``host_prob`` for the host-posterior strategy (spec 3.6)."""
    t = res.table
    S = np.asarray(t["S"], float)
    t1 = np.asarray(t["tier"]) == 1
    Z = float(S[t1].sum()) + U_deg
    t["host_prob"] = np.where(t1, S / Z if Z > 0 else 0.0, np.nan)
    t["host_prob"].description = ("Normalised host probability under the host-posterior generative model "
                                  "(experimental); sum over galaxies + p_uncat = 1")
    res.summary["p_uncat_host_model"] = U_deg / Z if Z > 0 else 1.0


_DESCRIPTIONS = {
    "rank": "Rank (tier 1 first, then tier 2); ties broken by gal_id",
    "tier": "1 = has a distance; 2 = no distance, ranked by sky density only (never interleaved)",
    "p_sky_deg2": "GW sky probability density at the galaxy (containing pixel), per deg^2",
    "credible_level": "Smallest credible region containing the galaxy (searched probability)",
    "p_dist_mpc3": "N_j Gauss(d; mu_j, sqrt(sigma_j^2 + s^2)): GW density per unit volume / sky density, Mpc^-3",
    "dist_zscore": "(d - mu) / sqrt(sigma^2 + s^2) in the galaxy's pixel",
    "w": "Astrophysical weight (tracer ** alpha), solar units",
    "S": "Ranking score p_sky * p_dist * w. NOT a host probability",
    "score_fraction": "S / sum(S) within the tier: conditional on the host being in this list; NOT a host probability",
    "dist_err_mpc": "Total distance uncertainty s_i including peculiar velocity where applicable",
}


def _describe(t: Table) -> None:
    for k, v in _DESCRIPTIONS.items():
        if k in t.colnames:
            t[k].description = v
    units = {"ra": "deg", "dec": "deg", "dist_mpc": "Mpc", "dist_err_mpc": "Mpc", "p_sky_deg2": "deg-2",
             "p_dist_mpc3": "Mpc-3"}
    for k, v in units.items():
        if k in t.colnames:
            t[k].unit = v


def _empty_table() -> Table:
    t = Table(names=["rank", "tier", "gal_id", "alt_ids", "ra", "dec", "dist_mpc", "dist_err_mpc", "dist_kind",
                     "z", "p_sky_deg2", "credible_level", "p_dist_mpc3", "dist_zscore", "w", "w_kind", "S",
                     "log10_S", "score_fraction", "flags"],
              dtype=["i8", "i2", "U32", "U64", "f8", "f8", "f8", "f8", "U8", "f8", "f8", "f8", "f8", "f8", "f8",
                     "U8", "f8", "f8", "f8", "U64"])
    _describe(t)
    return t

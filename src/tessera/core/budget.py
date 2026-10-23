"""Catalog incompleteness as a first-class output (spec 4.3).

    P_cat   = integral dP/dV * min(1, C_w) dV       (host is cataloged, under the weight model)
    P_uncat = P_in_footprint - P_cat
    P_out   = 1 - P_in_footprint

computed on the pixels of the ``budget_level`` credible region; the
remaining sky probability is reported as ``p_unassessed`` and is never
silently attributed to the catalog.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np
from astropy_healpix import lonlat_to_healpix

from tessera.core import science
from tessera.core.catalog import CUBE_ORDER, Snapshot
from tessera.core.config import ResolvedConfig
from tessera.core.distance import shell_probability
from tessera.core.errors import Diagnostics
from tessera.core.geometry import SR_PER_DEG2
from tessera.core.skymap import MultiOrderMap


def _pix6(m: MultiOrderMap, idx: np.ndarray) -> np.ndarray:
    ra, dec = m.pixel_centers(idx)
    return np.asarray(lonlat_to_healpix(ra * u.deg, dec * u.deg, 2**CUBE_ORDER, order="nested"), dtype=np.int64)


def low_coverage_mask(m: MultiOrderMap, snap: Snapshot, weight: str, level: float, threshold: float) -> np.ndarray:
    """Map pixels (in the ``level`` region) where C_w at the pixel's DISTMU is below ``threshold``.

    Pixels without distance information count as low coverage (they cannot be
    covered by galaxy targeting).  Used by the hybrid field strategy.
    """
    out = np.zeros(len(m), bool)
    cube = snap.completeness(weight) or snap.completeness("lum_K")
    region = np.flatnonzero(m.region_mask_inclusive(level))
    if cube is None or m.distmu is None:
        out[region] = True
        return out
    edges, C = cube
    mu = m.distmu[region]
    valid = m.distance_valid[region]
    shell = np.clip(np.searchsorted(edges, np.where(valid, mu, -1), side="right") - 1, 0, len(edges) - 2)
    c = C[_pix6(m, region), shell]
    out[region] = ~valid | (c < threshold)
    return out


def compute_budget(
    m: MultiOrderMap,
    snap: Snapshot,
    cfg: ResolvedConfig,
    *,
    assoc_summary: dict[str, Any] | None = None,
    diagnostics: Diagnostics | None = None,
) -> dict[str, Any]:
    diags = diagnostics if diagnostics is not None else Diagnostics()
    region = np.flatnonzero(m.region_mask_inclusive(cfg.budget_level))
    prob = m.prob[region]
    p_assessed = float(prob.sum())
    pix6 = _pix6(m, region)
    fp = snap.footprint_mask()[pix6]
    p_in_fp = float(prob[fp].sum())
    out: dict[str, Any] = {
        "level": cfg.budget_level,
        "p_unassessed": max(0.0, 1.0 - p_assessed),
        "p_in_footprint": p_in_fp,
        "p_out_footprint": max(0.0, p_assessed - p_in_fp),
        "p_precut": (assoc_summary or {}).get("p_precut", 0.0),
        "n_tier2_galaxies": (assoc_summary or {}).get("n_tier2", 0),
    }

    weight = cfg.weight
    cube = snap.completeness(weight) if weight != "flat" else None
    if cube is None:
        fallback = "lum_K" if snap.completeness("lum_K") is not None else None
        if fallback:
            diags.add("BUDGET_WEIGHT_FALLBACK",
                      f"no completeness cube for weight {weight!r}; budget uses {fallback!r}",
                      level="info", requested=weight, used=fallback)
            weight, cube = fallback, snap.completeness(fallback)
    out["weight"] = weight
    ref = science.reference_density(weight) if cube is not None else None
    if ref is not None:
        out["reference_density"] = ref
        if ref.get("provisional"):
            diags.add("REFERENCE_DENSITY_PROVISIONAL",
                      f"reference density for {weight!r} is provisional pending science review", source=ref["source"])

    if cube is None or not m.has_distance_columns or m.distance_class() == "none":
        out.update({"p_cat": None, "p_uncat": None, "p_beyond_dmax": None, "by_distance": [],
                    "low_coverage_regions": None,
                    "reason": "no completeness cube" if cube is None else "map has no distance information"})
        diags.add("BUDGET_INCOMPLETE", f"catalog budget not computed: {out['reason']}")
        return out

    edges, C = cube
    mu, sg, nm = m.distmu[region], m.distsigma[region], m.distnorm[region]
    valid = m.distance_valid[region]
    P_shell = shell_probability(mu[:, None], sg[:, None], nm[:, None], edges[None, :-1], edges[None, 1:])
    Cpix = C[pix6].astype(float)                       # (n, nshell)
    w_pix = prob[:, None] * P_shell                   # probability per pixel-shell
    p_cat = float(np.sum(w_pix[fp] * np.minimum(1.0, Cpix[fp])))
    p_no_dist_pix = float(prob[~valid].sum())
    p_beyond_cube = float(np.sum(prob[valid] * np.clip(1.0 - P_shell[valid].sum(axis=1), 0, 1)))
    low = Cpix < 0.5
    p_beyond_dmax = float(np.sum(w_pix[fp & valid][low[fp & valid]]))

    shells_tot = w_pix.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        c_mean = np.where(shells_tot > 0, (w_pix * np.minimum(1, Cpix)).sum(axis=0) / shells_tot, np.nan)
    by_distance = [
        {"r_lo_mpc": float(edges[k]), "r_hi_mpc": float(edges[k + 1]), "p": float(shells_tot[k]),
         "completeness": None if not np.isfinite(c_mean[k]) else float(c_mean[k])}
        for k in range(len(edges) - 1) if shells_tot[k] > 1e-6
    ]

    # low-coverage regions inside the 90% area (spec 4.3)
    r90 = m.region_mask_inclusive(0.9)[region]
    shell_mu = np.clip(np.searchsorted(edges, np.where(valid, mu, -1), side="right") - 1, 0, len(edges) - 2)
    c_at_mu = Cpix[np.arange(len(region)), shell_mu]
    lowcov = r90 & valid & (c_at_mu < 0.5)
    ebv = snap.ebv(*m.pixel_centers(region)) if snap.meta.get("has_ebv") else None
    if ebv is not None:
        lowcov |= r90 & (ebv > 0.3)
    moc: dict[str, list[int]] = {}
    if lowcov.any():
        for p6 in np.unique(pix6[lowcov]):
            moc.setdefault(str(CUBE_ORDER), []).append(int(p6))
    out.update({
        "p_cat": p_cat,
        "p_uncat": max(0.0, p_in_fp - p_cat),
        "p_no_distance_pixels": p_no_dist_pix,
        "p_beyond_cube": p_beyond_cube,
        "p_beyond_dmax": p_beyond_dmax,
        "by_distance": by_distance,
        "low_coverage_regions": {
            "moc": moc,
            "p": float(prob[lowcov].sum()),
            "area_deg2": float(m.area[region][lowcov].sum() / SR_PER_DEG2),
            "criteria": "C_w < 0.5 at the pixel's DISTMU" + (" or E(B-V) > 0.3" if ebv is not None else ""),
        },
    })
    if p_cat < 0.5:
        diags.add("CATALOG_INCOMPLETE",
                  f"only {p_cat:.2f} of the probability is expected to be cataloged under the {weight} model; "
                  "consider tiling or hybrid strategies", p_cat=p_cat)

    # host-posterior normalisation term (spec 3.6), in the same units as S (deg^-2 * Mpc^-3 * w)
    rho_bar = float(ref["value"]) if ref else None
    if rho_bar is not None:
        U_sr = rho_bar * max(0.0, 1.0 - p_cat)
        out["host_model"] = {"rho_bar_w": rho_bar, "U_deg": U_sr * SR_PER_DEG2}
        s_sr = (assoc_summary or {}).get("sum_S_tier1_sr")
        if s_sr and cfg.weight == weight and cfg.weight_alpha == 1.0:
            ratio = s_sr / (rho_bar * p_cat) if p_cat > 0 else float("inf")
            out["host_model"]["consistency_ratio"] = ratio
            # Sum(S) estimates integral(dP/dV rho_cat dV); C_w is capped at 1, so overdense
            # regions legitimately give ratio > 1.  A ratio well below 1 means the cube
            # claims more completeness than the galaxies present can support.
            if ratio < 1 / 3:
                diags.add("BUDGET_CONSISTENCY",
                          "sum of galaxy scores is below rho_bar * P_cat by more than a factor 3; "
                          "the completeness cube or reference density may be inaccurate", ratio=ratio)
            elif ratio > 3:
                diags.add("BUDGET_OVERDENSE",
                          "localisation volume is overdense relative to rho_bar (C_w capped at 1)",
                          level="info", ratio=ratio)
    return out

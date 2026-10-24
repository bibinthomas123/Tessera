"""Stage A: candidate fields and their time-independent metrics (spec 6.1).

Each field is rasterised onto a common HEALPix order chosen from the
footprint (pixel <= 1/4 of the smallest chip dimension, capped).  Every
raster pixel carries sky probability rho(centre) * A_R and distance-modulus
statistics for the yield objective.  Field-pixel and field-galaxy incidences
are stored as CSR arrays so Stage B can evaluate marginal gains for all
fields at once.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy_healpix import healpix_to_lonlat
from scipy.spatial import cKDTree

from tessera.core.config import ResolvedConfig
from tessera.core.detectability import pixel_dm_stats
from tessera.core.errors import Diagnostics, ProfileError
from tessera.core.footprint import InstrumentFootprint
from tessera.core.geometry import chord_for_angle, radec_to_xyz
from tessera.core.profile import TelescopeProfile
from tessera.core.skymap import MultiOrderMap


@dataclass
class FieldSet:
    profile: TelescopeProfile
    footprint: InstrumentFootprint
    order: int
    mode: str
    field_id: np.ndarray            # str
    ra: np.ndarray
    dec: np.ndarray
    pa: np.ndarray
    origin: np.ndarray              # "grid" | "tile" | "galaxy"
    # CSR incidence field -> raster pixel (index into upix)
    pix_ptr: np.ndarray
    pix_idx: np.ndarray
    pix_w: np.ndarray
    pix_edge: np.ndarray
    # unique raster pixels
    upix: np.ndarray
    uprob: np.ndarray
    u_dm_mean: np.ndarray
    u_dm_sd: np.ndarray
    u_ebv: np.ndarray
    u_ra: np.ndarray
    u_dec: np.ndarray
    # CSR incidence field -> galaxy (index into association table, tier-1 only)
    gal_ptr: np.ndarray
    gal_idx: np.ndarray
    gal_w: np.ndarray
    # metrics
    p_sky: np.ndarray
    p_host: np.ndarray
    n_galaxies: np.ndarray
    raster_err: np.ndarray
    dm90: np.ndarray
    ebv_mean: np.ndarray
    status: np.ndarray              # "candidate" | "rejected"
    reason: np.ndarray
    host_norm: float
    host_metric: str                # "p_host_cat" | "score_fraction_in_field"
    timings: dict[str, float] = field(default_factory=dict)
    n_unreachable: int = 0          # in-region tiles dropped because the site can never observe them

    def __len__(self) -> int:
        return len(self.field_id)

    def active(self) -> np.ndarray:
        return np.flatnonzero(self.status == "candidate")


def _read_grid(profile: TelescopeProfile) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    base = profile._base_dir or Path(".")
    path = Path(profile.field_grid)
    if not path.is_absolute():
        path = base / path
    if not path.is_file():
        raise ProfileError(f"field grid {path} not found", code="PROFILE_FIELD_GRID_MISSING")
    ids, ra, dec, pa = [], [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            ids.append(str(row["field_id"]))
            ra.append(float(row["ra"]))
            dec.append(float(row["dec"]))
            pa.append(float(row.get("pa_deg") or profile.instrument.footprint.pa_deg[0]))
    return ids, np.array(ra), np.array(dec), np.array(pa)


def _tessellation(fp: InstrumentFootprint, overlap: float, dec_lo: float = -90.0, dec_hi: float = 90.0,
                  phase: float = 0.0, phase_ra: float | None = None) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Deterministic all-sky tiling in declination bands.

    Tile IDs depend only on the footprint and overlap, so the same sky position
    always gets the same ID; only bands overlapping [dec_lo, dec_hi] are emitted.
    ``phase`` shifts the grid by that fraction of a tile in declination and ``phase_ra`` (default:
    the same) in RA; shifted grids get IDs ``o<dec phase><ra phase>...``.
    """
    pra = phase if phase_ra is None else phase_ra
    prefix = "t" if not (phase or pra) else f"o{int(round(phase * 100)):02d}{int(round(pra * 100)):02d}_"
    w, h = fp.tile_steps(overlap)
    nb = max(1, int(np.ceil(180.0 / h)))
    hb = 180.0 / nb
    ras, decs, ids = [], [], []
    for b in range(nb):
        dc = -90.0 + (b + 0.5 + phase) * hb
        if dc >= 90.0:
            continue
        if dc + hb / 2 + fp.bounding_radius_deg < dec_lo or dc - hb / 2 - fp.bounding_radius_deg > dec_hi:
            continue
        edge = min(89.999, abs(dc) + hb / 2)
        nr = max(1, int(np.ceil(360.0 * np.cos(np.radians(edge)) / w)))
        rr = np.mod((np.arange(nr) + 0.5 + pra) * 360.0 / nr, 360.0)
        ras.append(rr)
        decs.append(np.full(nr, dc))
        ids += [f"{prefix}{b:04d}_{k:05d}" for k in range(nr)]
    return np.concatenate(ras), np.concatenate(decs), ids


def _region_tree(m: MultiOrderMap, level: float, mask: np.ndarray | None = None):
    sel = m.region_mask_inclusive(level)
    if mask is not None:
        sel &= mask
    idx = np.flatnonzero(sel)
    ra, dec = m.pixel_centers(idx)
    # pixel radius, so coarse pixels still intersect small fields
    rad = np.degrees(np.sqrt(m.area[idx]))
    return idx, cKDTree(radec_to_xyz(ra, dec)) if len(idx) else None, m.prob[idx], float(rad.max() if len(idx) else 0)


def _prefilter(m: MultiOrderMap, level: float, ra: np.ndarray, dec: np.ndarray, fp: InstrumentFootprint,
               max_fields: int, region_mask: np.ndarray | None = None) -> np.ndarray:
    """Indices of candidate centres whose footprint can intersect the region, best first (capped)."""
    idx, tree, prob, prad = _region_tree(m, level, region_mask)
    if tree is None or len(ra) == 0:
        return np.zeros(0, int)
    xyz = radec_to_xyz(ra, dec)
    reach = chord_for_angle(fp.bounding_radius_deg + prad)
    hits = tree.query_ball_point(xyz, reach)
    est = np.array([prob[h].sum() if h else 0.0 for h in hits])
    keep = np.flatnonzero(est > 0)
    order = keep[np.lexsort((keep, -est[keep]))]
    return order[:max_fields]


def _galaxy_fields(fp: InstrumentFootprint, pa0: float, gra: np.ndarray, gdec: np.ndarray, S: np.ndarray,
                   gid: np.ndarray, cum: float, max_n: int) -> tuple[np.ndarray, np.ndarray, list[str]]:
    total = S.sum()
    if total <= 0 or len(S) == 0:
        return np.zeros(0), np.zeros(0), []
    order = np.lexsort((gid, -S))
    tree = cKDTree(radec_to_xyz(gra, gdec))
    reach = chord_for_angle(fp.bounding_radius_deg)
    assigned = np.zeros(len(S), bool)
    centres: list[tuple[float, float]] = []
    ids: list[str] = []
    acc = 0.0
    for i in order:
        if len(centres) >= max_n or acc >= cum * total:
            break
        if assigned[i]:
            continue
        near = np.asarray(tree.query_ball_point(radec_to_xyz(gra[i], gdec[i]), reach), int)
        near = near[~assigned[near]]
        inside = near[fp.contains(gra[i], gdec[i], pa0, gra[near], gdec[near])]
        c_ra, c_dec = gra[i], gdec[i]
        if len(inside) > 1:
            xyz = (radec_to_xyz(gra[inside], gdec[inside]) * S[inside, None]).sum(axis=0)
            from tessera.core.geometry import xyz_to_radec

            n_ra, n_dec = xyz_to_radec(xyz)
            if np.all(fp.contains(float(n_ra), float(n_dec), pa0, gra[inside], gdec[inside])):
                c_ra, c_dec = float(n_ra), float(n_dec)
        assigned[inside] = True
        acc += S[inside].sum()
        centres.append((c_ra, c_dec))
        ids.append(f"g:{gid[i]}")
    if not centres:
        return np.zeros(0), np.zeros(0), []
    arr = np.asarray(centres)
    return arr[:, 0], arr[:, 1], ids


def choose_mode(cfg: ResolvedConfig, profile: TelescopeProfile, fp: InstrumentFootprint, has_tier1: bool,
                distance_ok: bool) -> str:
    mode = cfg.field_mode
    if mode != "auto":
        if mode == "grid" and not profile.field_grid:
            raise ProfileError("field_mode 'grid' requires profile.field_grid", code="PROFILE_FIELD_GRID_MISSING")
        return mode
    tiled = "grid" if profile.field_grid else "tessellate"
    if not (cfg.galaxy_targeting and has_tier1 and distance_ok):
        return tiled
    return "hybrid" if fp.area_deg2 <= 1.0 else tiled


def build_fields(
    m: MultiOrderMap,
    profile: TelescopeProfile,
    cfg: ResolvedConfig,
    galaxies,                      # astropy Table from association (ranked)
    *,
    host_norm: float,
    host_metric: str,
    low_coverage_pixmask: np.ndarray | None = None,
    ebv_fn=None,
    mode: str | None = None,
    diagnostics: Diagnostics | None = None,
) -> FieldSet:
    import time as _time

    diags = diagnostics if diagnostics is not None else Diagnostics()
    t0 = _time.perf_counter()
    fp = InstrumentFootprint(profile.instrument.footprint)
    fc = cfg.fields
    pa0 = float(profile.instrument.footprint.pa_deg[0])
    tier = np.asarray(galaxies["tier"]) if len(galaxies) else np.zeros(0, int)
    t1 = np.flatnonzero(tier == 1)
    gra = np.asarray(galaxies["ra"], float)[t1] if len(t1) else np.zeros(0)
    gdec = np.asarray(galaxies["dec"], float)[t1] if len(t1) else np.zeros(0)
    gS = np.asarray(galaxies["S"], float)[t1] if len(t1) else np.zeros(0)
    gid = np.asarray(galaxies["gal_id"]).astype(str)[t1] if len(t1) else np.zeros(0, str)
    mode = mode or choose_mode(cfg, profile, fp, len(t1) > 0, m.distance_class() != "none")

    ids: list[str] = []
    ras, decs, pas, origin = [], [], [], []
    n_unreachable = 0

    def add(i_ids, i_ra, i_dec, i_pa, org):
        ids.extend(i_ids)
        ras.append(np.asarray(i_ra, float))
        decs.append(np.asarray(i_dec, float))
        pas.append(np.asarray(i_pa, float))
        origin.extend([org] * len(i_ids))

    if mode in ("grid", "tessellate", "hybrid"):
        if profile.field_grid:
            gids, gra_, gdec_, gpa = _read_grid(profile)
            org = "grid"
        else:
            reg = m.region_mask_inclusive(fc.region_level)
            _, rdec = m.pixel_centers(np.flatnonzero(reg))
            pad = float(np.degrees(np.sqrt(m.area[reg].max()))) if reg.any() else 0.0
            lo, hi = (float(rdec.min()) - pad, float(rdec.max()) + pad) if reg.any() else (0.0, -1.0)
            gra_, gdec_, gids = _tessellation(fp, fc.overlap, lo, hi)
            gpa = np.full(len(gids), pa0)
            org = "tile"
        # When the region has more tiles than max_fields, drop centres this site can never observe
        # (culmination below the minimum altitude, or outside the mount's declination limits)
        # before the cap, so the cap is not spent on fields the telescope cannot reach.  Applied
        # only when the cap binds, so smaller regions keep exactly the same field set.
        culm = 90.0 - np.abs(gdec_ - profile.site.lat_deg)
        dlo, dhi = profile.mount.dec_limits_deg
        reach = (culm >= profile.constraints.min_alt() - fp.bounding_radius_deg) & \
            (gdec_ >= dlo - fp.bounding_radius_deg) & (gdec_ <= dhi + fp.bounding_radius_deg)
        if not reach.all():
            n_region = len(_prefilter(m, fc.region_level, gra_, gdec_, fp, 10**9))
            if n_region > fc.max_fields:
                n_unreachable = len(_prefilter(m, fc.region_level, gra_[~reach], gdec_[~reach], fp, 10**9))
            elif not len(_prefilter(m, fc.region_level, gra_[reach], gdec_[reach], fp, 10**9)):
                n_unreachable = n_region   # nothing in the region is reachable: say so explicitly
            if n_unreachable:
                diags.add("FIELDS_UNREACHABLE_FROM_SITE",
                          f"{n_unreachable} tiles in the region never rise above the minimum altitude or lie "
                          f"outside the mount's declination limits at {profile.key}", level="info",
                          profile=profile.key, n=n_unreachable)
                gids = [g for g, k in zip(gids, reach) if k] if not isinstance(gids, np.ndarray) else gids[reach]
                gra_, gdec_, gpa = gra_[reach], gdec_[reach], gpa[reach]
        rmask = low_coverage_pixmask if mode == "hybrid" else None
        if mode == "hybrid" and rmask is None:
            rmask = np.zeros(len(m), bool)
        sel = _prefilter(m, fc.region_level, gra_, gdec_, fp, fc.max_fields, rmask)
        base = (gids, gra_, gdec_, gpa)
        groups = [(sel, org, base)]
        if mode == "hybrid":
            # The rest of the region as low-priority "fill" tiles: the scheduler takes them only
            # when no galaxy or low-coverage field has positive gain, so the telescope does not
            # sit idle once the galaxy fields are done.
            full = _prefilter(m, fc.region_level, gra_, gdec_, fp, fc.max_fields)
            groups.append((full[~np.isin(full, sel)][:max(0, fc.max_fields - len(sel))], "fill", base))
        if fc.offset_fill and org == "tile":
            # The grid shifted by half a tile, also at fill priority.  Fixed tiles whose centres
            # never become observable leave the visible edge of the region uncovered; shifted
            # pointings can reach it, using time that would otherwise be idle.
            # Origin "offset" ranks below "fill": it is used only when no other field has gain.
            # Quarter-tile shifts too for large fields (>= 2 deg), whose coarse lattice otherwise
            # leaves much of the visible edge unreachable; small fields are already fine-grained.
            phases = [(0.5, 0.5), (0.5, 0.0), (0.0, 0.5)]
            if fp.min_dimension_deg >= 2.0:
                phases += [(0.25, 0.25), (0.75, 0.75), (0.25, 0.75), (0.75, 0.25)]
            parts = [_tessellation(fp, fc.overlap, lo, hi, phase=pd, phase_ra=pr) for pd, pr in phases]
            ora = np.concatenate([p[0] for p in parts])
            odec = np.concatenate([p[1] for p in parts])
            oids = [i for p in parts for i in p[2]]
            oculm = 90.0 - np.abs(odec - profile.site.lat_deg)
            ok = (oculm >= profile.constraints.min_alt() - fp.bounding_radius_deg) & \
                (odec >= dlo - fp.bounding_radius_deg) & (odec <= dhi + fp.bounding_radius_deg)
            oids = [g for g, k in zip(oids, ok) if k]
            ora, odec = ora[ok], odec[ok]
            used = sum(len(g[0]) for g in groups)
            osel = _prefilter(m, fc.region_level, ora, odec, fp, max(0, fc.max_fields - used))
            groups.append((osel, "offset", (oids, ora, odec, np.full(len(oids), pa0))))
        for rows, kind, (g_ids, g_ra, g_dec, g_pa) in groups:
            for pa in (profile.instrument.footprint.pa_deg if org == "tile" else [None]):
                add([f"{g_ids[i]}" + ("" if pa is None or len(profile.instrument.footprint.pa_deg) == 1
                                       else f"@pa{pa:g}") for i in rows],
                    g_ra[rows], g_dec[rows], (g_pa[rows] if pa is None else np.full(len(rows), pa)), kind)
    if mode in ("galaxy", "hybrid"):
        cra, cdec, cids = _galaxy_fields(fp, pa0, gra, gdec, gS, gid, fc.galaxy_fields_cum_score,
                                         fc.galaxy_fields_max)
        add(cids, cra, cdec, np.full(len(cids), pa0), "galaxy")
    ra = np.concatenate(ras) if ras else np.zeros(0)
    dec = np.concatenate(decs) if decs else np.zeros(0)
    pa = np.concatenate(pas) if pas else np.zeros(0)
    t_gen = _time.perf_counter() - t0

    # -- rasterise -----------------------------------------------------------
    t1_ = _time.perf_counter()
    order = fp.raster_order(fc.max_raster_order)
    rasters = [fp.rasterize(ra[i], dec[i], pa[i], order) for i in range(len(ra))]
    allpix = np.concatenate([r.ipix for r in rasters]) if rasters else np.zeros(0, np.int64)
    upix, inv = np.unique(allpix, return_inverse=True)
    lon, lat = healpix_to_lonlat(upix, 2**order, order="nested")
    u_ra, u_dec = lon.to_value(u.deg), lat.to_value(u.deg)
    midx = m.lookup(u_ra, u_dec)
    area_R = 4 * np.pi / (12 * 4**order)
    rho = np.where(midx >= 0, m.probdensity[np.clip(midx, 0, None)], 0.0)
    uprob = rho * area_R
    # distance-modulus statistics of the containing map pixels (for yield)
    dm_mean = np.full(len(upix), np.nan)
    dm_sd = np.full(len(upix), np.nan)
    if m.has_distance_columns and len(upix):
        um, uinv = np.unique(np.clip(midx, 0, None), return_inverse=True)
        mm, ss = pixel_dm_stats(m.distmu[um], m.distsigma[um], m.distnorm[um])
        dm_mean, dm_sd = mm[uinv], ss[uinv]
        dm_mean[midx < 0] = np.nan
    if np.all(~np.isfinite(dm_mean)) and m.meta is not None and m.meta.header.get("DISTMEAN"):
        dmean, dstd = float(m.meta.header["DISTMEAN"]), float(m.meta.header.get("DISTSTD") or 0)
        mm, ss = pixel_dm_stats(np.array([dmean]), np.array([max(dstd, 1e-3)]), np.array([1.0 / max(dmean, 1) ** 2]))
        dm_mean[:] = mm[0]
        dm_sd[:] = ss[0]
        diags.add("YIELD_GLOBAL_DISTANCE", "pixel distances unavailable; yield uses the map's DISTMEAN/DISTSTD",
                  level="info")
    u_ebv = np.zeros(len(upix))
    if ebv_fn is not None and len(upix):
        e = ebv_fn(u_ra, u_dec)
        if e is not None:
            u_ebv = np.asarray(e, float)
    ptr = np.zeros(len(rasters) + 1, np.int64)
    ptr[1:] = np.cumsum([len(r.ipix) for r in rasters])
    pix_idx = inv.astype(np.int64)
    pix_w = np.concatenate([r.weight for r in rasters]) if rasters else np.zeros(0)
    pix_edge = np.concatenate([r.boundary for r in rasters]) if rasters else np.zeros(0, bool)
    t_ras = _time.perf_counter() - t1_

    # -- galaxies per field ----------------------------------------------------
    t2_ = _time.perf_counter()
    gptr = np.zeros(len(ra) + 1, np.int64)
    gidx_l: list[np.ndarray] = []
    gw_l: list[np.ndarray] = []
    if len(t1) and len(ra):
        tree = cKDTree(radec_to_xyz(gra, gdec))
        reach = chord_for_angle(fp.bounding_radius_deg)
        near = tree.query_ball_point(radec_to_xyz(ra, dec), reach)
        for i, cand in enumerate(near):
            cand = np.asarray(sorted(cand), int)
            w = fp.coverage(ra[i], dec[i], pa[i], gra[cand], gdec[cand]) if len(cand) else np.zeros(0)
            k = w > 0
            gidx_l.append(t1[cand[k]])          # index into the full galaxy table
            gw_l.append(w[k])
            gptr[i + 1] = gptr[i] + int(k.sum())
    gal_idx = np.concatenate(gidx_l) if gidx_l else np.zeros(0, np.int64)
    gal_w = np.concatenate(gw_l) if gw_l else np.zeros(0)
    t_gal = _time.perf_counter() - t2_

    # -- metrics ---------------------------------------------------------------
    nf = len(ra)
    Sfull = np.asarray(galaxies["S"], float) if len(galaxies) else np.zeros(0)
    p_sky = np.array([np.sum(uprob[pix_idx[ptr[i]:ptr[i + 1]]] * pix_w[ptr[i]:ptr[i + 1]]) for i in range(nf)])
    err = np.array([0.5 * np.sum(uprob[pix_idx[ptr[i]:ptr[i + 1]]][pix_edge[ptr[i]:ptr[i + 1]]])
                    for i in range(nf)])
    p_host = np.array([np.sum(Sfull[gal_idx[gptr[i]:gptr[i + 1]]] * gal_w[gptr[i]:gptr[i + 1]]) / host_norm
                       if host_norm > 0 else 0.0 for i in range(nf)])
    n_gal = np.diff(gptr)
    dm90 = np.full(nf, np.nan)
    ebv_mean = np.zeros(nf)
    for i in range(nf):
        sl = pix_idx[ptr[i]:ptr[i + 1]]
        w = uprob[sl] * pix_w[ptr[i]:ptr[i + 1]]
        ok = np.isfinite(dm_mean[sl]) & (w > 0)
        if ok.any():
            ww = w[ok] / w[ok].sum()
            mean = np.sum(ww * dm_mean[sl][ok])
            var = np.sum(ww * (dm_sd[sl][ok] ** 2 + (dm_mean[sl][ok] - mean) ** 2))
            dm90[i] = mean + 1.2816 * np.sqrt(var)
            ebv_mean[i] = float(np.sum(ww * u_ebv[sl][ok]))
    status = np.full(nf, "candidate", dtype=object)
    reason = np.full(nf, "", dtype=object)
    # Relative to the best field as well: on very large maps every field holds < low_value_eps,
    # and an absolute cut alone would reject them all.
    eps = min(fc.low_value_eps, fc.low_value_rel * float(p_sky.max())) if nf else fc.low_value_eps
    low = (p_sky < eps) & (n_gal == 0)
    status[low] = "rejected"
    reason[low] = "REJECT_LOW_VALUE"
    if nf == 0:
        diags.add("NO_CANDIDATE_FIELDS", f"no candidate fields generated in mode {mode!r}", level="warning")
    fs = FieldSet(profile=profile, footprint=fp, order=order, mode=mode, field_id=np.asarray(ids, dtype=object),
                  ra=ra, dec=dec, pa=pa, origin=np.asarray(origin, dtype=object), pix_ptr=ptr, pix_idx=pix_idx,
                  pix_w=pix_w, pix_edge=pix_edge, upix=upix, uprob=uprob, u_dm_mean=dm_mean, u_dm_sd=dm_sd,
                  u_ebv=u_ebv, u_ra=u_ra, u_dec=u_dec, gal_ptr=gptr, gal_idx=gal_idx, gal_w=gal_w, p_sky=p_sky,
                  p_host=p_host, n_galaxies=n_gal, raster_err=err, dm90=dm90, ebv_mean=ebv_mean, status=status,
                  reason=reason, host_norm=host_norm, host_metric=host_metric,
                  timings={"generate": t_gen, "rasterize": t_ras, "galaxies": t_gal}, n_unreachable=n_unreachable)
    return fs


def fields_summary(fs: FieldSet) -> dict[str, Any]:
    act = fs.active()
    return {"mode": fs.mode, "raster_order": fs.order, "n_fields": len(fs), "n_candidate": int(len(act)),
            "footprint_area_deg2": fs.footprint.area_deg2,
            "p_sky_union_all_candidates": float(_union_psky(fs, act))}


def _union_psky(fs: FieldSet, rows: np.ndarray) -> float:
    if len(rows) == 0:
        return 0.0
    cov = np.zeros(len(fs.upix))
    for i in rows:
        sl = slice(fs.pix_ptr[i], fs.pix_ptr[i + 1])
        np.maximum.at(cov, fs.pix_idx[sl], fs.pix_w[sl])
    return float(np.sum(cov * fs.uprob))

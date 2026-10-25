"""The engine entry point: ``plan()``.

Pure and deterministic: given the same sky map, catalog snapshot, profiles,
config and ``t_start``, it returns byte-identical outputs (spec section 1).
"""

from __future__ import annotations

import time as _time
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.table import Table
from astropy.time import Time
from astropy_healpix import lonlat_to_healpix

from tessera.core import timeutil
from tessera.core.association import add_host_probability, score_galaxies
from tessera.core.budget import compute_budget, low_coverage_mask
from tessera.core.catalog import Snapshot, open_snapshot
from tessera.core.config import PlanConfig, ResolvedConfig
from tessera.core.detectability import R_FILTER, efficiency
from tessera.core.errors import CancelledError, ConfigError, Diagnostics
from tessera.core.fields import FieldSet, build_fields, choose_mode, fields_summary
from tessera.core.footprint import InstrumentFootprint
from tessera.core.hashing import sha256_bytes
from tessera.core.io import doc_to_bytes, table_to_bytes
from tessera.core.manifest import build_manifest, plan_identity
from tessera.core.profile import TelescopeProfile, load_profile
from tessera.core.scheduler import Coverage, GreedyScheduler, ScheduleResult
from tessera.core.skymap import MultiOrderMap, load_skymap
from tessera.core.validate import validate_visits
from tessera.core.visibility import Visibility, compute_visibility

STAGES = ("parse", "catalog", "score", "budget", "fields", "visibility", "schedule", "validate")


@dataclass
class PlanResult:
    plan_key: str
    galaxies: Table
    fields: dict[str, Table]
    schedules: dict[str, dict[str, Any]]
    budget: dict[str, Any]
    strategy_comparison: dict[str, Any]
    diagnostics: list[dict[str, Any]]
    manifest: dict[str, Any]
    summary: dict[str, Any]
    timings: dict[str, float] = field(default_factory=dict)
    _outputs: dict[str, bytes] = field(default_factory=dict, repr=False)

    def outputs(self) -> dict[str, bytes]:
        """All output files (relative path -> bytes), including the manifest."""
        return {**self._outputs, "manifest.json": doc_to_bytes(self.manifest)}

    def write(self, out_dir: str | Path) -> list[Path]:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        written = []
        for rel, data in self.outputs().items():
            p = out / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(data)
            written.append(p)
        return written


def _vis_for(fs: FieldSet, uid: dict[str, int], vis_all: Visibility) -> Visibility:
    """Slice the shared visibility computation down to one field set."""
    idx = np.array([uid.get(str(f), -1) for f in fs.field_id], dtype=int)
    sub = vis_all.subset(np.clip(idx, 0, None))
    for k in np.flatnonzero(idx < 0):      # rejected at Stage A, or an offset field computed lazily
        sub.windows[k] = []
        sub.never_reason[k] = "deferred" if fs.origin[k] == "offset" else "rejected_stage_a"
    return sub


def _check(cancel_check: Callable[[], bool] | None, stage: str) -> None:
    if cancel_check is not None and cancel_check():
        raise CancelledError(f"plan cancelled before stage {stage!r}")


def _coerce_config(config: PlanConfig | dict[str, Any] | None) -> PlanConfig:
    if config is None:
        return PlanConfig()
    if isinstance(config, PlanConfig):
        return config
    try:
        return PlanConfig.model_validate(config)
    except Exception as exc:
        raise ConfigError(f"invalid plan config: {exc}") from exc


def _transfer_coverage(prev: list[tuple[FieldSet, Coverage]], fs: FieldSet, n_gal: int) -> Coverage:
    sky = np.zeros(len(fs.upix))
    miss = np.ones(len(fs.upix))
    gal = np.zeros(n_gal)
    for pfs, pcov in prev:
        if len(pfs.upix) and len(fs.upix):
            pix = np.asarray(lonlat_to_healpix(fs.u_ra * u.deg, fs.u_dec * u.deg, 2**pfs.order, order="nested"))
            pos = np.clip(np.searchsorted(pfs.upix, pix), 0, len(pfs.upix) - 1)
            hit = pfs.upix[pos] == pix
            sky = np.where(hit, np.maximum(sky, pcov.sky[pos]), sky)
            miss = np.where(hit, np.minimum(miss, pcov.miss[pos]), miss)
        gal = np.maximum(gal, pcov.gal)
    return Coverage(sky, miss, gal)


def _credit_observations(fs: FieldSet, cov: Coverage, obs: list[dict[str, Any]], sched: GreedyScheduler,
                         gal_ra: np.ndarray, gal_dec: np.ndarray) -> int:
    """Credit completed observations against the new map (incremental replanning)."""
    n = 0
    for o in obs:
        if str(o.get("status", "completed")) != "completed":
            continue
        if o.get("profile_id") and o["profile_id"] != fs.profile.id:
            continue
        pa = float(o.get("pa_deg", fs.profile.instrument.footprint.pa_deg[0]))
        r = fs.footprint.rasterize(float(o["ra"]), float(o["dec"]), pa, fs.order)
        pos = np.clip(np.searchsorted(fs.upix, r.ipix), 0, max(len(fs.upix) - 1, 0))
        hit = (fs.upix[pos] == r.ipix) if len(fs.upix) else np.zeros(0, bool)
        p, w = pos[hit], r.weight[hit]
        np.maximum.at(cov.sky, p, w)
        mlim = o.get("m_lim")
        if mlim is not None and sched.t_merger is not None:
            start = timeutil.to_time(o["start"]).jd
            filt = o.get("filter", fs.profile.visit_pattern[0].filter)
            M = sched.lc.abs_mag(np.full(len(p), start - sched.t_merger), filt)
            eff = np.nan_to_num(efficiency(float(mlim) - M - R_FILTER.get(filt, 2.5) * fs.u_ebv[p],
                                           fs.u_dm_mean[p], fs.u_dm_sd[p], sched.lc.sigma_mag))
        else:
            eff = np.ones(len(p))
        np.minimum.at(cov.miss, p, 1.0 - w * eff)
        if len(gal_ra):
            cov.gal = np.maximum(cov.gal, fs.footprint.coverage(float(o["ra"]), float(o["dec"]), pa, gal_ra, gal_dec))
        n += 1
    return n


def _fields_table(fs: FieldSet, vis: Visibility, sched: ScheduleResult, yields: np.ndarray | None,
                  nights: int) -> Table:
    scheduled: dict[str, int] = {}
    for v in sched.visits:
        scheduled[v["field_id"]] = scheduled.get(v["field_id"], 0) + 1
    n = len(fs)
    status, reason = [], []
    for i in range(n):
        fid = str(fs.field_id[i])
        if fid in scheduled:
            status.append("scheduled")
            reason.append(sched.unscheduled.get(fid + "#revisit", ""))
        elif fs.status[i] == "rejected":
            status.append("rejected")
            reason.append(str(fs.reason[i]))
        else:
            r = sched.unscheduled.get(fid, "LOWER_GAIN_THAN_SCHEDULED")
            status.append("rejected" if r.startswith("REJECT") else "not_selected")
            reason.append(r)
    t = Table()
    t["field_id"] = fs.field_id.astype(str)
    t["profile"] = np.full(n, fs.profile.key)
    t["origin"] = fs.origin.astype(str)
    t["ra"] = fs.ra
    t["dec"] = fs.dec
    t["pa_deg"] = fs.pa
    t["p_sky"] = fs.p_sky
    t[fs.host_metric] = fs.p_host
    t["n_galaxies"] = fs.n_galaxies
    t["raster_err"] = fs.raster_err
    t["dm90"] = fs.dm90
    if yields is not None:
        for k in range(nights):
            t[f"y_night{k + 1}"] = yields[:, k]
    t["windows"] = ["; ".join(f"{timeutil.iso(a)}/{timeutil.iso(b)}" for a, b in w) for w in vis.windows]
    t["n_visits"] = [scheduled.get(str(f), 0) for f in fs.field_id]
    t["status"] = status
    t["reason"] = reason
    t["ra"].unit = t["dec"].unit = t["pa_deg"].unit = "deg"
    t["p_sky"].description = "Sky probability inside the footprint (rasterised)"
    t[fs.host_metric].description = (
        "P(host in field) under the host-posterior model" if fs.host_metric == "p_host_cat"
        else "Score fraction of cataloged galaxies in the field: conditional on the host being cataloged")
    t["raster_err"].description = "Bound on P_sky rasterisation error (half the probability of boundary pixels)"
    t["dm90"].description = "90th-percentile distance modulus in the field (exposure-time solver)"
    order = np.lexsort((t["field_id"], -np.asarray(t["p_sky"])))
    return t[order]


def _night_yields(sched: GreedyScheduler, fs: FieldSet, vis: Visibility, t0: float, nights: int) -> np.ndarray | None:
    if sched.objective != "yield" and not np.isfinite(fs.u_dm_mean).any():
        return None
    if sched.t_merger is None or not np.isfinite(fs.u_dm_mean).any():
        return None
    out = np.full((len(fs), nights), 0.0)
    filt = fs.profile.visit_pattern[0].filter
    fresh = GreedyScheduler.__new__(GreedyScheduler)
    fresh.__dict__.update(sched.__dict__)
    fresh.cov = Coverage(np.zeros(len(fs.upix)), np.ones(len(fs.upix)), np.zeros(len(sched.S)))
    fresh.objective = "yield"
    for k in range(nights):
        lo, hi = t0 + k, t0 + k + 1
        rows, starts = [], []
        for i, ws in enumerate(vis.windows):
            for a, b in ws:
                if b > lo and a < hi:
                    rows.append(i)
                    starts.append(max(a, lo))
                    break
        if not rows:
            continue
        rows_a, st = np.asarray(rows), np.asarray(starts)
        n_exp, ts = fresh._exposures(rows_a, st, filt)
        end = st + n_exp * (ts + fs.profile.instrument.readout_s) / 86400.0
        ml = fresh._m_lim(rows_a, filt, n_exp, ts, st, end)
        out[rows_a, k] = fresh._gains(rows_a, filt, ml, st)
    return out


def plan(
    skymap: MultiOrderMap | str | Path | bytes,
    catalog: Snapshot | str | Path,
    profiles: Sequence[TelescopeProfile | str | Path | dict[str, Any]] = (),
    config: PlanConfig | dict[str, Any] | None = None,
    t_start: str | Time | None = None,
    *,
    event: dict[str, Any] | None = None,
    prior_observations: list[dict[str, Any]] | None = None,
    cancel_check: Callable[[], bool] | None = None,
    progress: Callable[[str], None] | None = None,
) -> PlanResult:
    """Rank host galaxies and produce validated schedules for each profile."""
    timings: dict[str, float] = {}
    diags = Diagnostics()

    def stage(name: str):
        _check(cancel_check, name)
        if progress:
            progress(name)
        return _time.perf_counter()

    t = stage("parse")
    m = skymap if isinstance(skymap, MultiOrderMap) else load_skymap(skymap)
    diags.extend(m.diagnostics)
    timings["parse"] = _time.perf_counter() - t
    snap = catalog if isinstance(catalog, Snapshot) else open_snapshot(catalog)
    profs = [p if isinstance(p, TelescopeProfile) else load_profile(p) for p in profiles]
    if len({p.id for p in profs}) != len(profs):
        raise ConfigError("each telescope profile may appear only once per plan")
    pcfg = _coerce_config(config)
    cfg: ResolvedConfig = pcfg.resolved()
    if cfg.preset_note:
        diags.add("PRESET_NOTE", cfg.preset_note, level="info")
    if cfg.strategy_def.get("experimental"):
        diags.add("STRATEGY_EXPERIMENTAL", f"strategy {cfg.strategy!r} is experimental", level="warning")
    if t_start is None:
        now = Time.now()
        t_start_t = Time(np.floor(now.jd * 1440.0) / 1440.0, format="jd", scale="utc")
        diags.add("PLAN_START_NOW", "t_start not given; using the current time (recorded in the manifest)",
                  level="info")
    else:
        t_start_t = timeutil.to_time(t_start)
    t_start_t = Time(round(t_start_t.jd * 86400.0) / 86400.0, format="jd", scale="utc")
    t_end_t = t_start_t + cfg.horizon_nights * u.day
    if cfg.t_merger:
        t_merger = timeutil.to_time(cfg.t_merger)
    elif m.meta is not None and m.meta.mjd_obs:
        t_merger = Time(m.meta.mjd_obs, format="mjd", scale="utc")
    else:
        t_merger = None
        diags.add("NO_MERGER_TIME", "no MJD-OBS in the map and no config.t_merger; yield objective unavailable")

    # -- catalog extraction ---------------------------------------------------
    t = stage("catalog")
    region = m.region_mask_inclusive(cfg.cr_max)
    gals = snap.extract_for_map(m.order[region], m.ipix[region])
    ebv = snap.ebv(gals["ra"], gals["dec"]) if len(gals) else None
    if not snap.meta.get("has_ebv"):
        diags.add("NO_DUST_MAP", "catalog snapshot has no E(B-V) map; extinction set to 0", level="info")
    timings["catalog"] = _time.perf_counter() - t

    t = stage("score")
    assoc = score_galaxies(m, gals, cfg, ebv=ebv, diagnostics=diags)
    timings["score"] = _time.perf_counter() - t

    t = stage("budget")
    budget = compute_budget(m, snap, cfg, assoc_summary=assoc.summary, diagnostics=diags)
    host_metric = "score_fraction_in_field"
    host_norm = float(assoc.summary.get("sum_S_tier1", 0.0))
    if cfg.host_posterior:
        hm = budget.get("host_model")
        if hm is None:
            raise ConfigError("host-posterior strategy needs a completeness cube and a distance map",
                              code="HOST_MODEL_UNAVAILABLE")
        add_host_probability(assoc, hm["U_deg"])
        host_norm += hm["U_deg"]
        host_metric = "p_host_cat"
        budget["p_uncat_host_model"] = assoc.summary.get("p_uncat_host_model")
    timings["budget"] = _time.perf_counter() - t

    lowcov = None
    nights = max(1, int(np.ceil(cfg.horizon_nights)))
    S_all = np.asarray(assoc.table["S"], float) if len(assoc.table) else np.zeros(0)
    S_all = np.where(np.asarray(assoc.table["tier"]) == 1, S_all, 0.0) if len(assoc.table) else S_all

    fields_out: dict[str, Table] = {}
    schedules: dict[str, dict[str, Any]] = {}
    comparison: dict[str, Any] = {}
    prev_cov: list[tuple[FieldSet, Coverage]] = []
    gal_fields: dict[int, list[str]] = {}
    for key in ("fields", "visibility", "schedule", "validate"):
        timings[key] = 0.0

    for prof in profs:
        t = stage("fields")
        fp = InstrumentFootprint(prof.instrument.footprint)
        has_t1 = bool(assoc.summary.get("n_tier1", 0))
        primary = choose_mode(cfg, prof, fp, has_t1, m.distance_class() != "none")
        tiled = "grid" if prof.field_grid else "tessellate"
        modes = [primary]
        if cfg.compare_strategies and has_t1:
            for alt in (tiled, "galaxy", "hybrid"):
                if alt not in modes:
                    modes.append(alt)
        if ("hybrid" in modes) and lowcov is None:
            lowcov = low_coverage_mask(m, snap, budget.get("weight") or "lum_K", cfg.fields.region_level,
                                       cfg.fields.hybrid_completeness_threshold)
        fsets: dict[str, FieldSet] = {}
        for md in modes:
            fsets[md] = build_fields(m, prof, cfg, assoc.table, host_norm=host_norm, host_metric=host_metric,
                                     low_coverage_pixmask=lowcov, ebv_fn=snap.ebv if snap.meta.get("has_ebv") else None,
                                     mode=md, diagnostics=diags if md == primary else Diagnostics())
        timings["fields"] += _time.perf_counter() - t

        # one visibility computation over the union of all candidate fields; shifted-grid
        # "offset" fields are computed later, only if a schedule variant uses them
        t = stage("visibility")
        uid: dict[str, int] = {}
        ura, udec = [], []
        off: dict[str, tuple[float, float]] = {}
        for fs in fsets.values():
            for i in fs.active():
                f = str(fs.field_id[i])
                if fs.origin[i] == "offset":
                    off.setdefault(f, (float(fs.ra[i]), float(fs.dec[i])))
                elif f not in uid:
                    uid[f] = len(ura)
                    ura.append(fs.ra[i])
                    udec.append(fs.dec[i])
        vis_all = compute_visibility(prof, np.asarray(ura), np.asarray(udec), t_start_t, t_end_t,
                                     step_s=cfg.scheduler.grid_step_s, refine_s=cfg.scheduler.refine_s)
        timings["visibility"] += _time.perf_counter() - t
        off_cache: dict[str, Any] = {}

        def offset_loader(fs: FieldSet, off=off, off_cache=off_cache, prof=prof):
            def load(target: Visibility) -> None:
                t_ = _time.perf_counter()
                if "vis" not in off_cache:
                    ids = list(off)
                    off_cache["uid"] = {f: k for k, f in enumerate(ids)}
                    off_cache["vis"] = compute_visibility(
                        prof, np.array([off[f][0] for f in ids]), np.array([off[f][1] for f in ids]),
                        t_start_t, t_end_t, step_s=cfg.scheduler.grid_step_s, refine_s=cfg.scheduler.refine_s)
                ov, ou = off_cache["vis"], off_cache["uid"]
                for k in np.flatnonzero(fs.origin == "offset"):
                    j = ou.get(str(fs.field_id[k]))
                    if j is None:
                        continue
                    target.windows[k] = ov.windows[j]
                    target.never_reason[k] = ov.never_reason[j]
                    target.alt[k] = ov.alt[j]
                    target.az[k] = ov.az[j]
                timings["visibility"] += _time.perf_counter() - t_
            return load if off else None

        t = stage("schedule")
        results: dict[str, tuple[FieldSet, Visibility, GreedyScheduler, ScheduleResult]] = {}
        for md, fs in fsets.items():
            vis = _vis_for(fs, uid, vis_all)
            cov = None
            if md == primary:
                if cfg.coordination == "sequential" and prev_cov:
                    cov = _transfer_coverage(prev_cov, fs, len(S_all))
            sch = GreedyScheduler(fs, vis, cfg, t_start=t_start_t, t_end=t_end_t, t_merger=t_merger,
                                  galaxies_S=S_all, coverage=cov,
                                  diagnostics=diags if md == primary else Diagnostics(),
                                  load_offset_visibility=offset_loader(fs))
            if md == primary and prior_observations:
                g_ra = np.asarray(assoc.table["ra"], float) if len(assoc.table) else np.zeros(0)
                g_dec = np.asarray(assoc.table["dec"], float) if len(assoc.table) else np.zeros(0)
                n_cred = _credit_observations(fs, sch.cov, prior_observations, sch, g_ra, g_dec)
                if n_cred:
                    diags.add("PRIOR_OBSERVATIONS_CREDITED", f"{n_cred} completed observations credited",
                              level="info", profile=prof.key)
            results[md] = (fs, vis, sch, sch.run())
        timings["schedule"] += _time.perf_counter() - t

        t = stage("validate")
        fs, vis, sch, res = results[primary]
        report = validate_visits(res.visits, prof, horizon=(t_start_t.jd, t_end_t.jd))
        by_seq = {r["seq"]: r for r in report["visits"]}
        for v in res.visits:
            r = by_seq[v["seq"]]
            v["feasible"] = r["feasible"]
            v["violations"] = r["violations"]
            v["airmass"] = r["airmass"]
            v["moon_sep_deg"] = r["moon_sep_min_deg"]
            v["moon_illumination"] = r["moon_illumination"]
            v["sun_alt_max_deg"] = r["sun_alt_max_deg"]
        if not report["feasible"]:
            diags.add("SCHEDULE_INFEASIBLE", f"validator rejected {report['n_infeasible']} visits for {prof.key}",
                      level="error", violations=report["violations_by_constraint"])
        timings["validate"] += _time.perf_counter() - t

        prev_cov.append((fs, res.coverage))
        yields = _night_yields(sch, fs, vis, t_start_t.jd, nights)
        fields_out[prof.key] = _fields_table(fs, vis, res, yields, nights)
        for i in range(len(fs)):
            for gi in fs.gal_idx[fs.gal_ptr[i]:fs.gal_ptr[i + 1]]:
                gal_fields.setdefault(int(gi), []).append(f"{prof.id}:{fs.field_id[i]}")
        for v in res.visits:
            v.pop("start_jd", None)
            v.pop("end_jd", None)
        schedules[prof.key] = {
            "schema": "tessera.schedule/1",
            "version": 0,
            "parent_version": None,
            "status": "proposed" if report["feasible"] else "infeasible",
            "profile": prof.key,
            "profile_sha256": prof.sha256,
            "profile_validation": prof.validation.model_dump(),
            "horizon": {"start": timeutil.iso(t_start_t), "end": timeutil.iso(t_end_t)},
            "objective": res.summary["objective"],
            "scheduler": cfg.scheduler.algorithm,
            "field_mode": primary,
            "light_curve": cfg.light_curve,
            "feasible": report["feasible"],
            "validation": {k: report[k] for k in ("feasible", "n_visits", "n_infeasible",
                                                  "violations_by_constraint", "step_s")},
            "visits": res.visits,
            "summary": {**res.summary, "p_uncat": budget.get("p_uncat"), "fields": fields_summary(fs)},
            "unscheduled_counts": {**_hist(res.unscheduled),
                                   **({"REJECT_UNREACHABLE_FROM_SITE": fs.n_unreachable} if fs.n_unreachable else {})},
            "empty_reason": res.empty_reason,
        }
        if res.empty_reason:
            diags.add("EMPTY_SCHEDULE", f"no visits scheduled for {prof.key}: {res.empty_reason}",
                      profile=prof.key, rejection_histogram=res.summary.get("rejection_histogram", {}))
        comp = {}
        for md, (cfs, _, _, cres) in results.items():
            comp[md] = {"n_visits": cres.summary["n_visits"], "cum_p_sky": cres.summary["cum_p_sky"],
                        f"cum_{cfs.host_metric}": cres.summary.get(f"cum_{cfs.host_metric}"),
                        "cum_expected_detections": cres.summary["cum_expected_detections"],
                        "expected_detections_per_hour": cres.summary["expected_detections_per_hour"],
                        "telescope_time_h": round((cres.summary["t_exposure_s"] + cres.summary["t_overhead_s"]) / 3600, 3),
                        "objective": cres.summary["objective"]}
        comparison[prof.key] = {"primary": primary, "strategies": comp,
                                "recommendation": _recommend(comp, budget, fp, results, host_metric, t_start_t)}

    profs_by_key = {p.key: p for p in profs}
    # galaxy -> field ids
    tbl = assoc.table
    if len(tbl):
        tbl["field_ids"] = ["|".join(gal_fields.get(i, [])) for i in range(len(tbl))]

    # -- identity, outputs, manifest --------------------------------------------
    cfg_dump = cfg.model_dump(mode="json")
    prof_ids = [{"key": p.key, "sha256": p.sha256} for p in profs]
    plan_key, ident = plan_identity(
        skymap_sha256=m.meta.sha256_data if m.meta else sha256_bytes(m.uniq.tobytes() + m.probdensity.tobytes()),
        snapshot_id=snap.id, profiles=prof_ids, config=cfg_dump, t_start=timeutil.iso(t_start_t), event=event,
        prior_observations=prior_observations)
    meta_common = {"plan_key": plan_key, "strategy": cfg.strategy, "science_model": cfg.science_model,
                   "catalog": snap.describe(), "skymap_sha256": ident["skymap_sha256"]}
    outputs: dict[str, bytes] = {}
    tbl_out = tbl
    if len(tbl) > pcfg.max_output_galaxies:
        # Distant events can put 10^6 catalogue galaxies in the region; writing them all as CSV,
        # ECSV and JSON dominated planning time.  Files keep the top-ranked rows.
        tbl_out = tbl[: pcfg.max_output_galaxies]
        diags.add("GALAXY_OUTPUT_TRUNCATED",
                  f"galaxies.* list the top {pcfg.max_output_galaxies:,} of {len(tbl):,} scored galaxies "
                  f"({float(np.sum(tbl_out['score_fraction'])):.4f} of the total score)", level="info",
                  n_scored=len(tbl), n_written=pcfg.max_output_galaxies)
    for fmt in ("csv", "ecsv", "json"):
        outputs[f"galaxies.{fmt}"] = table_to_bytes(tbl_out, fmt, meta={**meta_common, "summary": assoc.summary,
                                                                         "budget": budget})
    for key, ft in fields_out.items():
        for fmt in ("csv", "json"):
            outputs[f"fields/{key}.{fmt}"] = table_to_bytes(ft, fmt, meta={**meta_common, "profile": key})
    for key, sd in schedules.items():
        sd["plan_key"] = plan_key
        sd["schedule_id"] = "sch_" + sha256_bytes(f"{plan_key}|{key}".encode())[:24]
        outputs[f"schedules/{key}.json"] = doc_to_bytes(sd)
    outputs["region.json"] = doc_to_bytes({
        "levels": {str(lv): m.region_moc(lv, max_order=10) for lv in (0.5, 0.9)},
        "area_deg2": {str(lv): m.area_deg2(lv) for lv in (0.5, 0.9)},
        "preview": m.preview(0.9, order=min(7, m.max_order)),
        "footprints": {k: _footprint_outlines(profs_by_key[k], ft) for k, ft in fields_out.items()},
    })
    outputs["budget.json"] = doc_to_bytes(budget)
    outputs["strategy_comparison.json"] = doc_to_bytes(comparison)
    outputs["diagnostics.json"] = doc_to_bytes({"diagnostics": diags.to_list()})
    summary = {
        "plan_key": plan_key,
        "event": event or {},
        "skymap": m.summary(),
        "galaxies": {k: assoc.summary.get(k) for k in ("n_candidates", "n_scored", "n_tier1", "n_tier2",
                                                        "strategy", "weight")},
        "top_galaxies": [{"rank": int(r["rank"]), "gal_id": str(r["gal_id"]), "score_fraction": float(r["score_fraction"])}
                         for r in tbl[:10]] if len(tbl) else [],
        "budget": {k: budget.get(k) for k in ("p_cat", "p_uncat", "p_out_footprint", "p_unassessed", "weight")},
        "schedules": {k: {"status": v["status"], "n_visits": len(v["visits"]), "summary": v["summary"]}
                      for k, v in schedules.items()},
        "diagnostic_codes": diags.codes(),
    }
    outputs["summary.json"] = doc_to_bytes(summary)
    out_hashes = {k: sha256_bytes(v) for k, v in sorted(outputs.items())}
    manifest = build_manifest(
        plan_key=plan_key, identity=ident, event=event, skymap=m.meta.to_dict() if m.meta else {},
        catalog={**snap.describe(), "meta": {k: snap.meta.get(k) for k in ("source_url", "retrieved", "files")}},
        profiles=[{"key": p.key, "sha256": p.sha256, "validation": p.validation.model_dump()} for p in profs],
        config=cfg_dump, t_start=timeutil.iso(t_start_t), t_end=timeutil.iso(t_end_t),
        environment=timeutil.environment_info(), outputs=out_hashes, timings=timings)
    manifest["request_config"] = pcfg.model_dump(mode="json")
    manifest["prior_observations"] = prior_observations or []
    return PlanResult(plan_key=plan_key, galaxies=tbl, fields=fields_out, schedules=schedules, budget=budget,
                      strategy_comparison=comparison, diagnostics=diags.to_list(), manifest=manifest,
                      summary=summary, timings=timings, _outputs=outputs)


def _footprint_outlines(prof: TelescopeProfile, ft: Table, limit: int = 300) -> list[dict[str, Any]]:
    """Sky outlines of scheduled fields (for the UI and coordination tools)."""
    fp = InstrumentFootprint(prof.instrument.footprint)
    out = []
    for r in ft:
        if r["status"] != "scheduled":
            continue
        out.append({"field_id": str(r["field_id"]), "polygons": fp.sky_polygons(float(r["ra"]), float(r["dec"]),
                                                                               float(r["pa_deg"]))})
        if len(out) >= limit:
            break
    return out


def _hist(unscheduled: dict[str, str]) -> dict[str, int]:
    h: dict[str, int] = {}
    for r in unscheduled.values():
        k = r.split(":")[0]
        h[k] = h.get(k, 0) + 1
    return dict(sorted(h.items()))


def _recommend(comp: dict[str, Any], budget: dict[str, Any], fp: InstrumentFootprint, results, host_metric: str,
               t_start: Time) -> dict[str, Any]:
    """Spec 4.4: galaxy targeting when p_cat >= 0.5, FOV <= 1 deg2 and the fields holding
    50% of sum(S) fit in the first night; otherwise tiling or hybrid, with the failed conditions."""
    p_cat = budget.get("p_cat")
    reasons = []
    ok_cat = p_cat is not None and p_cat >= 0.5
    if not ok_cat:
        reasons.append(f"p_cat = {p_cat if p_cat is None else round(p_cat, 3)} < 0.5")
    ok_fov = fp.area_deg2 <= 1.0
    if not ok_fov:
        reasons.append(f"field of view {fp.area_deg2:.2f} deg2 > 1 deg2")
    ok_night = False
    if "galaxy" in results:
        _, _, _, gres = results["galaxy"]
        first_night_end = t_start.jd + 1.0
        for v in gres.visits:
            if Time(v["end"].replace("Z", ""), scale="utc").jd > first_night_end:
                break
            if v["cum"].get(host_metric, 0) >= 0.5 * (1.0 if host_metric == "score_fraction_in_field"
                                                       else max(1e-12, 1 - (budget.get("p_uncat_host_model") or 0))):
                ok_night = True
                break
    if not ok_night:
        reasons.append("galaxy fields holding 50% of the score do not fit in the first night")
    if ok_cat and ok_fov and ok_night:
        return {"strategy": "galaxy", "reasons": ["all galaxy-targeting conditions met"]}
    best = max((k for k in comp if k in ("tessellate", "grid", "hybrid")),
               key=lambda k: (comp[k]["cum_expected_detections"], comp[k]["cum_p_sky"]), default=None)
    return {"strategy": best or "tessellate", "reasons": reasons}

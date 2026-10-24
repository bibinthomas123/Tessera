"""Benchmark and validation report (spec section 11), produced by ``tessera benchmark``.

Level 1  synthetic injections: top-k host recovery for 3D strategies vs the 2D
         baseline, prior mismatch, budget accuracy under known incompleteness
Level 2  historical public events: GW170817 (NGC 4993 rank, earliest visit),
         S190425z (large single-detector map), S200105ae (multi-order NSBH)
Sensitivity: kilonova diversity sigma_M (0.25 / 0.5 / 1.0 mag) on yield-driven schedules
Engineering: determinism, latency per stage, feasibility of every schedule,
         malformed-input rejection, edge-case maps
"""

from __future__ import annotations

import gzip
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from tessera.core.association import score_galaxies
from tessera.core.budget import compute_budget
from tessera.core.catalog import build_snapshot, open_snapshot
from tessera.core.config import PlanConfig
from tessera.core.errors import SkymapError, TesseraError
from tessera.core.hashing import sha256_bytes
from tessera.core.plan import plan
from tessera.core.skymap import load_skymap

NGC4993 = "GLADE+:1151336"


def _snapshot(data: Path, out: Path):
    snap_dir = out / "catalogs" / "gladeplus-gw170817-excerpt"
    if (snap_dir / "snapshot.json").exists():
        return open_snapshot(snap_dir)
    from tessera.adapters.catalogs import glade_plus

    cols = glade_plus.read(data / "glade_plus_gw170817.tsv")
    return build_snapshot(cols, snap_dir, name="gladeplus-gw170817-excerpt", upstream_version="2022-10-31",
                          source_url="https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=VII/291",
                          citation=glade_plus.CITATION)


def _rank_of(table, gal_id: str) -> int | None:
    ids = list(np.asarray(table["gal_id"]).astype(str))
    return ids.index(gal_id) + 1 if gal_id in ids else None


def historical(data: Path, out: Path) -> dict[str, Any]:
    snap = _snapshot(data, out)
    res: dict[str, Any] = {}
    m = load_skymap(data / "GW170817_skymap.fits.gz")
    strategies = {}
    for st in ("vol3d-lumK", "vol3d-mstar", "vol3d-lumB", "vol3d-flat", "sky2d-flat"):
        cfg = PlanConfig(preset="bns", strategy=st).resolved()
        reg = m.region_mask_inclusive(cfg.cr_max)
        g = snap.extract_for_map(m.order[reg], m.ipix[reg])
        a = score_galaxies(m, g, cfg)
        r = _rank_of(a.table, NGC4993)
        strategies[st] = {"rank": r, "score_fraction": float(a.table["score_fraction"][r - 1]) if r else None,
                          "n_scored": len(a.table)}
    t0 = time.perf_counter()
    p = plan(m, snap, ["example-1m"], PlanConfig(preset="bns", horizon_nights=1), "2017-08-17T22:00:00",
             event={"superevent_id": "GW170817", "revision": 1})
    wall = time.perf_counter() - t0
    row = p.galaxies[np.asarray(p.galaxies["gal_id"]).astype(str) == NGC4993][0]
    fids = set(f.split(":", 1)[1] for f in str(row["field_ids"]).split("|") if f)
    sched = p.schedules["example-1m@1"]
    first = next((v for v in sched["visits"] if v["field_id"] in fids), None)
    res["GW170817"] = {
        "map": "GWTC-1 GW170817_skymap.fits.gz (flat NSIDE 1024)",
        "area_90_deg2": m.area_deg2(0.9),
        "ngc4993": {"credible_level": float(row["credible_level"]), "by_strategy": strategies},
        "example_1m_first_visit_covering_ngc4993": None if first is None else {
            "start": first["start"], "seq": first["seq"], "field_id": first["field_id"],
            "hours_after_merger": round((_jd(first["start"]) - (m.meta.mjd_obs + 2400000.5)) * 24, 2)},
        "reference": "SSS17a discovered ~10.9 h after merger by Swope (Coulter et al. 2017)",
        "schedule": {"status": sched["status"], "n_visits": len(sched["visits"]),
                     "summary": {k: sched["summary"][k] for k in ("cum_p_sky", "cum_score_fraction_in_field",
                                                                    "cum_expected_detections", "open_shutter_frac")}},
        "strategy_comparison": p.strategy_comparison["example-1m@1"],
        "budget": {k: p.budget.get(k) for k in ("p_cat", "p_uncat", "p_out_footprint", "p_unassessed")},
        "timings_s": {k: round(v, 3) for k, v in p.timings.items()},
        "wall_s": round(wall, 2),
        "diagnostics": [d["code"] for d in p.diagnostics],
    }
    for name, fn, preset, profs, t in (
        ("S190425z", "S190425z_bayestar.fits.gz", "bns", ["example-wide-north", "rubin-lsstcam-approx"],
         "2019-04-25T08:30:00"),
        ("S200105ae", "S200105ae_bayestar.multiorder.fits", "nsbh", ["example-wide-north", "rubin-lsstcam-approx"],
         "2020-01-05T16:40:00"),
    ):
        path = data / fn
        if not path.exists():
            continue
        t0 = time.perf_counter()
        mm = load_skymap(path)
        p = plan(mm, snap, profs, PlanConfig(preset=preset, horizon_nights=1, coordination="sequential"), t,
                 event={"superevent_id": name, "revision": 1})
        res[name] = {
            "area_90_deg2": mm.area_deg2(0.9), "distmean": mm.meta.header.get("DISTMEAN"),
            "wall_s": round(time.perf_counter() - t0, 2),
            "timings_s": {k: round(v, 3) for k, v in p.timings.items()},
            "schedules": {k: {"status": v["status"], "n_visits": len(v["visits"]),
                              "cum_p_sky": v["summary"]["cum_p_sky"], "field_mode": v["field_mode"],
                              "empty_reason": v["empty_reason"]} for k, v in p.schedules.items()},
            "budget": {k: p.budget.get(k) for k in ("p_cat", "p_uncat")},
            "diagnostics": [d["code"] for d in p.diagnostics],
        }
    return res


SIGMA_GRID = (0.25, 0.5, 1.0)


def sigma_sensitivity(data: Path, out: Path) -> dict[str, Any]:
    """Re-plan with the at2017gfo light curve at several sigma_M; compare to the adopted value."""
    import dataclasses

    from tessera.core import detectability

    snap = _snapshot(data, out)
    adopted = detectability.LightCurve.named("at2017gfo").sigma_mag
    cases = [("GW170817", "GW170817_skymap.fits.gz", "example-1m", "2017-08-17T22:00:00"),
             ("S190425z", "S190425z_bayestar.fits.gz", "example-wide-north", "2019-04-25T08:30:00")]
    original = detectability._load_light_curve
    res: dict[str, Any] = {"adopted_sigma_mag": adopted, "grid": list(SIGMA_GRID), "events": {}}
    try:
        for name, fn, prof, t in cases:
            if not (data / fn).exists():
                continue
            m = load_skymap(data / fn)
            runs = {}
            for sig in SIGMA_GRID:
                def patched(n, _sig=sig):
                    lc = original(n)
                    return dataclasses.replace(lc, sigma_mag=_sig) if n == "at2017gfo" else lc
                detectability._load_light_curve = patched
                p = plan(m, snap, [prof], PlanConfig(preset="bns", horizon_nights=1), t,
                         event={"superevent_id": name, "revision": 1})
                sch = p.schedules[f"{prof}@1"]
                runs[sig] = {"n_visits": len(sch["visits"]),
                             "cum_expected_detections": sch["summary"]["cum_expected_detections"],
                             "cum_p_sky": sch["summary"]["cum_p_sky"],
                             "fields": [v["field_id"] for v in sch["visits"]]}
            ref = set(runs[adopted]["fields"]) if adopted in runs else set()
            for r in runs.values():
                f = set(r.pop("fields"))
                r["field_overlap_with_adopted"] = (len(f & ref) / len(f | ref)) if (f | ref) else 1.0
            res["events"][name] = {"profile": prof, "distmean_mpc": m.meta.header.get("DISTMEAN"),
                                   "by_sigma": {str(k): v for k, v in runs.items()}}
    finally:
        detectability._load_light_curve = original
    return res


def _jd(iso: str) -> float:
    from tessera.core.timeutil import to_time

    return float(to_time(iso).jd)


def injections(out: Path, n: int, seed: int = 20170817) -> dict[str, Any]:
    from tessera.testing import injection_map, synthetic_universe

    rng = np.random.default_rng(seed)
    uni = synthetic_universe(rng, ra0=150.0, dec0=0.0, radius_deg=14.0, d_max=160.0)
    sdir = out / "catalogs" / f"synthetic-{seed}"
    snap = open_snapshot(sdir) if (sdir / "snapshot.json").exists() else build_snapshot(
        uni.cols, sdir, name="synthetic", upstream_version=str(seed), footprint="allsky")
    # incomplete version: keep each galaxy with probability 0.5
    keep = rng.uniform(size=len(uni.L)) < 0.5
    half_dir = out / "catalogs" / f"synthetic-{seed}-half"
    half_cols = {k: v[keep] for k, v in uni.cols.items()}
    half = open_snapshot(half_dir) if (half_dir / "snapshot.json").exists() else build_snapshot(
        half_cols, half_dir, name="synthetic-half", upstream_version=str(seed), footprint="allsky")
    half_ids = set(half_cols["gal_id"])
    inner = np.flatnonzero((uni.d_true < 100) & (uni.d_true > 10)
                           & (_sep(uni.cols["ra"], uni.cols["dec"], 150.0, 0.0) < 8.0))
    p_l = uni.L[inner] / uni.L[inner].sum()
    ks = (1, 10, 100)
    strategies = ("vol3d-lumK", "vol3d-flat", "sky2d-flat")
    hits = {s: {k: 0 for k in ks} for s in strategies}
    mismatch = {k: 0 for k in ks}
    pred_cat, in_cat, pred_full = [], [], []
    for i in range(n):
        host = int(rng.choice(inner, p=p_l)) if i % 2 == 0 else int(rng.choice(inner))
        m, _ = injection_map(rng, uni.cols["ra"][host], uni.cols["dec"][host], uni.d_true[host], sigma_deg=1.5,
                             order=7)
        hid = str(uni.cols["gal_id"][host])
        for st in strategies:
            cfg = PlanConfig(preset="bns", strategy=st).resolved()
            reg = m.region_mask_inclusive(cfg.cr_max)
            g = snap.extract_for_map(m.order[reg], m.ipix[reg])
            a = score_galaxies(m, g, cfg)
            r = _rank_of(a.table, hid)
            for k in ks:
                if r is not None and r <= k:
                    if i % 2 == 0:
                        hits[st][k] += 1
                    elif st == "vol3d-lumK":
                        mismatch[k] += 1
        if i % 2 == 0:
            cfg = PlanConfig(preset="bns").resolved()
            reg = m.region_mask_inclusive(cfg.cr_max)
            g = half.extract_for_map(m.order[reg], m.ipix[reg])
            a = score_galaxies(m, g, cfg)
            b = compute_budget(m, half, cfg, assoc_summary=a.summary)
            pred_cat.append(b["p_cat"])
            in_cat.append(hid in half_ids)
            g_full = snap.extract_for_map(m.order[reg], m.ipix[reg])
            a_full = score_galaxies(m, g_full, cfg)
            pred_full.append(compute_budget(m, snap, cfg, assoc_summary=a_full.summary)["p_cat"])
    truth = float(uni.L[inner][keep[inner]].sum() / uni.L[inner].sum())
    n_matched = (n + 1) // 2
    n_mis = n // 2
    rec = {s: {f"top{k}": hits[s][k] / n_matched for k in ks} for s in strategies}
    return {
        "n_injections": n, "seed": seed, "n_galaxies": int(len(uni.L)),
        "host_prior_matched": "hosts drawn proportional to L_K (same as ranking prior)",
        "top_k_recovery": rec,
        "prior_mismatch_flat_hosts_ranked_lumK": {f"top{k}": mismatch[k] / max(1, n_mis) for k in ks},
        "budget_accuracy": {
            "mean_predicted_p_cat_half": float(np.mean(pred_cat)) if pred_cat else None,
            "mean_predicted_p_cat_full": float(np.mean(pred_full)) if pred_full else None,
            "predicted_relative_completeness": (float(np.mean(pred_cat) / np.mean(pred_full))
                                                if pred_full and np.mean(pred_full) > 0 else None),
            "truth_weighted_kept_fraction": truth,
            "empirical_host_in_catalog": float(np.mean(in_cat)) if in_cat else None,
            "n_empirical": len(in_cat),
            "abs_difference": (float(abs(np.mean(pred_cat) / np.mean(pred_full) - truth))
                               if pred_full and np.mean(pred_full) > 0 else None),
            "abs_difference_absolute_p_cat": float(abs(np.mean(pred_cat) - truth)) if pred_cat else None,
            "note": ("catalog randomly halved. Completeness is measured against a fixed cosmic mean density, "
                     "so a locally under- or over-dense volume also moves P_cat (the full catalog scores "
                     "'mean_predicted_p_cat_full'). The target (|difference| <= 0.05) is applied to the "
                     "halved/full ratio, which isolates incompleteness; the absolute difference is also shown."),
        },
        "3d_beats_2d_at_all_k": all(rec["vol3d-lumK"][f"top{k}"] >= rec["sky2d-flat"][f"top{k}"] for k in ks),
    }


def _sep(ra, dec, ra0, dec0):
    from tessera.core.geometry import angsep

    return angsep(ra0, dec0, ra, dec)


def engineering(data: Path, out: Path) -> dict[str, Any]:
    snap = _snapshot(data, out)
    m = load_skymap(data / "GW170817_skymap.fits.gz")
    cfg = PlanConfig(preset="bns", horizon_nights=1)
    a = plan(m, snap, ["example-1m"], cfg, "2017-08-17T22:00:00")
    b = plan(m, snap, ["example-1m"], cfg, "2017-08-17T22:00:00")
    ha = {k: sha256_bytes(v) for k, v in a.outputs().items() if k != "manifest.json"}
    hb = {k: sha256_bytes(v) for k, v in b.outputs().items() if k != "manifest.json"}
    # malformed corpus
    good = (data / "S200105ae_bayestar.multiorder.fits").read_bytes()
    corpus = {
        "empty": b"",
        "not_fits": b"hello world" * 100,
        "truncated": good[: len(good) // 3],
        "gzip_garbage": b"\x1f\x8b" + b"\x00" * 100,
        "gzip_bomb": gzip.compress(b"\x00" * (64 * 1024**2)),
    }
    from tessera.core.skymap import SkymapLimits

    rejected = {}
    for name, blob in corpus.items():
        try:
            load_skymap(blob, limits=SkymapLimits(max_decompressed_bytes=32 * 1024**2))
            rejected[name] = "ACCEPTED (unexpected)"
        except SkymapError as exc:
            rejected[name] = exc.code
        except Exception as exc:  # unhandled = failure of this check
            rejected[name] = f"UNHANDLED {type(exc).__name__}"
    from tessera.testing import edge_case_maps

    edges = {}
    for name, em in edge_case_maps().items():
        try:
            r = plan(em, snap, ["example-1m"], PlanConfig(preset="bns", horizon_nights=1, compare_strategies=False),
                     "2017-08-17T22:00:00")
            s = r.schedules["example-1m@1"]
            edges[name] = {"ok": True, "status": s["status"], "n_visits": len(s["visits"]),
                           "empty_reason": s["empty_reason"], "diagnostics": [d["code"] for d in r.diagnostics]}
        except TesseraError as exc:
            edges[name] = {"ok": False, "error": exc.code}
        except Exception as exc:
            edges[name] = {"ok": False, "error": f"UNHANDLED {type(exc).__name__}: {exc}"}
    return {
        "optimality": _optimality(m, snap),
        "deterministic": ha == hb, "n_outputs_compared": len(ha),
        "malformed_inputs": rejected,
        "malformed_all_rejected_with_typed_error": all(v.startswith("SKYMAP_") for v in rejected.values()),
        "edge_case_maps": edges,
        "edge_cases_all_handled": all(v.get("ok") for v in edges.values()),
    }


def _optimality(m, snap) -> dict[str, Any]:
    """greedy-v1 against a provable MILP upper bound (tessera.core.optimality)."""
    from astropy.time import Time

    from tessera.core.fields import build_fields
    from tessera.core.optimality import greedy_vs_bound
    from tessera.core.profile import load_profile
    from tessera.core.visibility import compute_visibility

    cfg = PlanConfig(preset="bns", horizon_nights=1).resolved()
    reg = m.region_mask_inclusive(cfg.cr_max)
    a = score_galaxies(m, snap.extract_for_map(m.order[reg], m.ipix[reg]), cfg)
    S = np.where(np.asarray(a.table["tier"]) == 1, np.asarray(a.table["S"]), 0.0)
    p = load_profile("example-1m")
    t0 = Time("2017-08-17T22:00:00", scale="utc")
    out = {}
    for mode in ("tessellate", "galaxy"):
        fs = build_fields(m, p, cfg, a.table, host_norm=float(S.sum()), host_metric="score_fraction_in_field",
                          mode=mode)
        vis = compute_visibility(p, fs.ra, fs.dec, t0, t0 + 1)
        for obj in ("sky_prob", "host_prob"):
            r = greedy_vs_bound(fs, vis, cfg, t_start=t0, t_end=t0 + 1,
                                t_merger=Time(m.meta.mjd_obs, format="mjd"), galaxies_S=S, objective=obj)
            out[f"{mode}/{obj}"] = {k: r.get(k) for k in ("greedy", "upper_bound", "ratio_lower_bound", "status",
                                                          "n_fields")}
    return out


def run_benchmark(data: Path, out: Path, *, n_injections: int = 200, quick: bool = False) -> dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    hist = historical(data, out)
    inj = injections(out, 20 if quick else n_injections)
    eng = engineering(data, out)
    sens = sigma_sensitivity(data, out)
    feasible = all(v["status"] != "infeasible" for ev in hist.values()
                   for v in (ev.get("schedules") or {"_": ev.get("schedule", {"status": "proposed"})}).values())
    from tessera import SCIENCE_MODEL_VERSION, __version__
    from tessera.core.manifest import software_versions

    rep = {
        "tessera": __version__, "science_model": SCIENCE_MODEL_VERSION, "versions": software_versions(),
        "headline": {
            "ngc4993_rank_vol3d_lumK": hist["GW170817"]["ngc4993"]["by_strategy"]["vol3d-lumK"]["rank"],
            "ngc4993_rank_sky2d": hist["GW170817"]["ngc4993"]["by_strategy"]["sky2d-flat"]["rank"],
            "gw170817_plan_wall_s": hist["GW170817"]["wall_s"],
            "top10_recovery_vol3d_lumK": inj["top_k_recovery"]["vol3d-lumK"]["top10"],
            "top10_recovery_sky2d": inj["top_k_recovery"]["sky2d-flat"]["top10"],
            "3d_beats_2d_at_all_k": inj["3d_beats_2d_at_all_k"],
            "budget_abs_difference": inj["budget_accuracy"]["abs_difference"],
            "all_schedules_feasible": feasible,
            "deterministic": eng["deterministic"],
            "malformed_all_rejected": eng["malformed_all_rejected_with_typed_error"],
            "edge_cases_all_handled": eng["edge_cases_all_handled"],
            "greedy_vs_bound_min_ratio": min((v["ratio_lower_bound"] for v in eng["optimality"].values()
                                              if v.get("ratio_lower_bound") is not None), default=None),
        },
        "historical": hist, "injections": inj, "engineering": eng, "sigma_sensitivity": sens,
        "runtime_s": round(time.perf_counter() - t0, 1),
    }
    (out / "report.json").write_text(json.dumps(rep, indent=2, sort_keys=True, default=str))
    (out / "report.md").write_text(_markdown(rep))
    return rep


def _markdown(rep: dict[str, Any]) -> str:
    h = rep["headline"]
    inj = rep["injections"]
    lines = [f"# tessera benchmark report (tessera {rep['tessera']}, science model {rep['science_model']})", "",
             "## Headline", "", "| Check | Result |", "| --- | --- |"]
    lines += [f"| {k} | {v} |" for k, v in h.items()]
    lines += ["", "## Synthetic injections (top-k host recovery)", "",
              f"{inj['n_injections']} injections, seed {inj['seed']}, {inj['n_galaxies']} synthetic galaxies.", "",
              "| Strategy | top1 | top10 | top100 |", "| --- | --- | --- | --- |"]
    for s, r in inj["top_k_recovery"].items():
        lines.append(f"| {s} | {r['top1']:.2f} | {r['top10']:.2f} | {r['top100']:.2f} |")
    pm = inj["prior_mismatch_flat_hosts_ranked_lumK"]
    lines += ["", f"Prior mismatch (hosts drawn flat, ranked with L_K): top1 {pm['top1']:.2f}, "
                  f"top10 {pm['top10']:.2f}, top100 {pm['top100']:.2f}.", "",
              f"Budget accuracy (catalog halved): predicted P_cat {inj['budget_accuracy']['mean_predicted_p_cat_half']:.3f} "
              f"(full catalog {inj['budget_accuracy']['mean_predicted_p_cat_full']:.3f}; this synthetic volume is "
              f"not exactly at the mean density), halved/full {inj['budget_accuracy']['predicted_relative_completeness']:.3f} "
              f"vs truth {inj['budget_accuracy']['truth_weighted_kept_fraction']:.3f}; empirical host-in-catalog rate "
              f"{inj['budget_accuracy']['empirical_host_in_catalog']:.3f} over {inj['budget_accuracy']['n_empirical']} "
              "injections.", "",
              "## Historical events", ""]
    g = rep["historical"]["GW170817"]
    lines += ["### GW170817", "", "| Strategy | NGC 4993 rank | score fraction |", "| --- | --- | --- |"]
    for s, r in g["ngc4993"]["by_strategy"].items():
        sf = f"{r['score_fraction']:.3f}" if r["score_fraction"] is not None else "-"
        lines.append(f"| {s} | {r['rank']} | {sf} |")
    lines += ["", f"Credible level of NGC 4993: {g['ngc4993']['credible_level']:.3f}. "
                  f"First example-1m visit covering it: {g['example_1m_first_visit_covering_ngc4993']}. "
                  f"Reference: {g['reference']}.", ""]
    for name in ("S190425z", "S200105ae"):
        if name in rep["historical"]:
            e = rep["historical"][name]
            lines += [f"### {name}", "", f"90% area {e['area_90_deg2']:.0f} deg2, wall {e['wall_s']} s.", "",
                      "| Profile | status | visits | cum P_sky | mode |", "| --- | --- | --- | --- | --- |"]
            for k, v in e["schedules"].items():
                lines.append(f"| {k} | {v['status']} | {v['n_visits']} | {v['cum_p_sky']:.3f} | {v['field_mode']} |")
            lines.append("")
    sens = rep.get("sigma_sensitivity", {})
    if sens.get("events"):
        lines += ["## Sensitivity to kilonova diversity sigma_M", "",
                  f"Adopted sigma_M = {sens['adopted_sigma_mag']} mag (Ascenzi et al. 2019). Field overlap is the "
                  "Jaccard index of scheduled fields against the adopted value.", "",
                  "| Event | profile | sigma_M | visits | expected detections | cum P_sky | field overlap |",
                  "| --- | --- | --- | --- | --- | --- | --- |"]
        for ev, e in sens["events"].items():
            for sig, r in e["by_sigma"].items():
                lines.append(f"| {ev} | {e['profile']} | {sig} | {r['n_visits']} | "
                             f"{r['cum_expected_detections']:.3f} | {r['cum_p_sky']:.3f} | "
                             f"{r['field_overlap_with_adopted']:.2f} |")
        lines.append("")
    lines += ["## Greedy vs. provable upper bound (GW170817, example-1m, 1 night)", "",
              "| fields / objective | greedy | upper bound | greedy / bound (>= greedy / optimum) |",
              "| --- | --- | --- | --- |"]
    for k, v in rep["engineering"]["optimality"].items():
        if v.get("upper_bound") is not None:
            lines.append(f"| {k} | {v['greedy']:.4f} | {v['upper_bound']:.4f} | {v['ratio_lower_bound']:.3f} |")
    lines += ["", "## Engineering", "", f"Deterministic: {rep['engineering']['deterministic']}", "",
              "Malformed inputs: " + ", ".join(f"{k}={v}" for k, v in rep["engineering"]["malformed_inputs"].items()),
              "", f"Runtime of this report: {rep['runtime_s']} s", ""]
    return "\n".join(lines)

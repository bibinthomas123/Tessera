"""Command-line interface: ``tessera <command>``.  Run ``tessera --help``."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml


def _cmd_plan(a: argparse.Namespace) -> int:
    from tessera.core.config import PlanConfig
    from tessera.core.plan import plan

    cfg: dict = {}
    if a.config:
        cfg = yaml.safe_load(Path(a.config).read_text()) or {}
    if a.preset:
        cfg["preset"] = a.preset
    if a.strategy:
        cfg["strategy"] = a.strategy
    if a.nights:
        cfg["horizon_nights"] = a.nights
    event = {"superevent_id": a.event_id, "revision": a.revision} if a.event_id else None
    prior = json.loads(Path(a.observations).read_text()) if a.observations else None
    r = plan(a.skymap, a.catalog, a.profile or [], PlanConfig.model_validate(cfg), a.t_start, event=event,
             prior_observations=prior, progress=(lambda st: print(f"  stage: {st}", file=sys.stderr)) if a.verbose else None)
    files = r.write(a.out)
    print(json.dumps({"plan_key": r.plan_key, "out": str(a.out), "files": len(files),
                      "top_galaxies": r.summary["top_galaxies"][:5], "budget": r.summary["budget"],
                      "schedules": {k: {"status": v["status"], "n_visits": v["n_visits"]}
                                    for k, v in r.summary["schedules"].items()},
                      "diagnostics": r.summary["diagnostic_codes"],
                      "timings_s": {k: round(v, 2) for k, v in r.timings.items()}}, indent=2))
    return 0 if all(v["status"] != "infeasible" for v in r.schedules.values()) else 2


def _cmd_ingest(a: argparse.Namespace) -> int:
    from tessera.core.catalog import build_snapshot

    if a.kind == "glade+":
        from tessera.adapters.catalogs import glade_plus

        cols = glade_plus.read(a.path)
        citation = glade_plus.CITATION
        name = a.name or "gladeplus"
    else:
        from tessera.adapters.catalogs import generic

        if not a.mapping:
            print("--mapping is required for generic catalogs", file=sys.stderr)
            return 1
        cols = generic.read(a.path, a.mapping)
        citation = a.citation
        name = a.name or Path(a.path).stem
    ebv = None
    if a.ebv_map:
        import numpy as np
        from astropy.io import fits

        with fits.open(a.ebv_map) as h:
            ebv = np.asarray(h[1].data[h[1].columns.names[0]], float).ravel()
            if str(h[1].header.get("ORDERING", "NESTED")).upper() == "RING":
                from astropy_healpix.core import ring_to_nested

                nside = int(np.sqrt(len(ebv) / 12))
                out = np.empty_like(ebv)
                out[ring_to_nested(np.arange(len(ebv)), nside)] = ebv
                ebv = out
    snap = build_snapshot(cols, a.out, name=name, upstream_version=a.upstream_version, source_url=a.source_url,
                          retrieved=a.retrieved, citation=citation, footprint=a.footprint, ebv_map=ebv)
    print(json.dumps({"snapshot_id": snap.id, **snap.describe(), "duplicates_flagged": snap.meta["duplicates_flagged"],
                      "dist_kind_counts": snap.meta["dist_kind_counts"]}, indent=2))
    return 0


def _cmd_validate(a: argparse.Namespace) -> int:
    from tessera.core.profile import load_profile
    from tessera.core.validate import validate_visits

    doc = json.loads(Path(a.schedule).read_text())
    prof = load_profile(a.profile or doc["profile"])
    rep = validate_visits(doc["visits"], prof)
    print(json.dumps({k: rep[k] for k in ("feasible", "n_visits", "n_infeasible", "violations_by_constraint")}, indent=2))
    if a.verbose:
        for v in rep["visits"]:
            if not v["feasible"]:
                print(json.dumps(v))
    return 0 if rep["feasible"] else 2


def _cmd_export(a: argparse.Namespace) -> int:
    from tessera.adapters.schedulers import get_adapter
    from tessera.core.profile import load_profile

    doc = json.loads(Path(a.schedule).read_text())
    prof = load_profile(a.profile or doc["profile"])
    ad = get_adapter(a.adapter)
    kw = {}
    if a.adapter == "rubin_too":
        from tessera.core.skymap import load_skymap

        if not a.skymap:
            print("--skymap is required for rubin_too", file=sys.stderr)
            return 1
        kw["skymap"] = load_skymap(a.skymap)
    if a.graceid:
        kw["graceid"] = a.graceid
    exp = ad.export(doc, prof, **kw)
    issues = ad.validate(exp)
    data, _, name = exp.as_bytes()
    out = Path(a.out or name)
    out.write_bytes(data)
    print(json.dumps({"written": str(out), "issues": [i.__dict__ for i in issues]}, indent=2))
    return 0 if not any(i.level == "error" for i in issues) else 2


def _cmd_rerun(a: argparse.Namespace) -> int:
    from tessera.core.hashing import sha256_bytes
    from tessera.core.plan import plan

    man = json.loads(Path(a.manifest).read_text())
    from tessera.core.config import PlanConfig

    pc = PlanConfig.model_validate(man["request_config"])
    r = plan(a.skymap, a.catalog, [p["key"].split("@")[0] for p in man["profiles"]], pc, man["t_start"],
             event=man.get("event") or None, prior_observations=man.get("prior_observations") or None)
    new = {k: sha256_bytes(v) for k, v in r.outputs().items() if k != "manifest.json"}
    diff = sorted(k for k in set(new) | set(man["outputs"]) if new.get(k) != man["outputs"].get(k))
    print(json.dumps({"plan_key_then": man["plan_key"], "plan_key_now": r.plan_key, "identical": not diff,
                      "differences": diff}, indent=2))
    return 0 if not diff else 3


def _cmd_compare(a: argparse.Namespace) -> int:
    from tessera.core.compare import compare_plans

    rep = compare_plans(a.dir_a, a.dir_b, top_k=a.top_k)
    text = json.dumps(rep, indent=2, default=str)
    if a.out:
        Path(a.out).write_text(text)
    print(text)
    return {"identical": 0, "within_tolerance": 0}.get(rep["verdict"], 3)


def _cmd_profiles(a: argparse.Namespace) -> int:
    from tessera.core.profile import list_profiles, load_profile

    if a.check:
        p = load_profile(a.check)
        print(json.dumps({"key": p.key, "sha256": p.sha256, "validation": p.validation.model_dump()}, indent=2))
        return 0
    for n in list_profiles():
        p = load_profile(n)
        print(f"{p.key:32s} {p.validation.level:20s} {p.description or ''}")
    return 0


def _cmd_keys(a: argparse.Namespace) -> int:
    from tessera.service.auth import create_key, revoke_key
    from tessera.service.db import session_factory
    from tessera.service.settings import get_settings

    S = session_factory(get_settings().database_url)
    with S() as s:
        if a.action == "create":
            key = create_key(s, a.name, [x.strip() for x in a.scopes.split(",")])
            print(key)
        else:
            print("revoked" if revoke_key(s, a.prefix) else "not found")
    return 0


def _cmd_serve(a: argparse.Namespace) -> int:  # pragma: no cover
    import uvicorn

    uvicorn.run("tessera.service.api:app_from_env", factory=True, host=a.host, port=a.port,
                proxy_headers=True, forwarded_allow_ips="*" if a.behind_proxy else None, log_config=None)
    return 0


def _cmd_worker(_: argparse.Namespace) -> int:  # pragma: no cover
    from tessera.service.worker import run_forever

    run_forever()
    return 0


def _cmd_listen(a: argparse.Namespace) -> int:  # pragma: no cover
    from tessera.service.listener import run

    run(a.source, a.path)
    return 0


def _cmd_replay(a: argparse.Namespace) -> int:
    from tessera.adapters.alerts.sources import directory
    from tessera.service.listener import consume

    outs = consume(directory(a.path, once=True))
    hist: dict[str, int] = {}
    for o in outs:
        hist[o.outcome] = hist.get(o.outcome, 0) + 1
    print(json.dumps({"alerts": len(outs), "outcomes": hist}, indent=2))
    return 0


def _cmd_benchmark(a: argparse.Namespace) -> int:
    from tessera.benchmark import run_benchmark

    rep = run_benchmark(Path(a.data), Path(a.out), n_injections=a.injections, quick=a.quick)
    print(json.dumps(rep["headline"], indent=2))
    return 0 if rep["headline"]["all_schedules_feasible"] else 2


def main(argv: list[str] | None = None) -> int:
    from tessera import __version__

    p = argparse.ArgumentParser(prog="tessera", description="GW follow-up engine")
    p.add_argument("--version", action="version", version=f"tessera {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("plan", help="rank galaxies and build schedules from a local sky map")
    sp.add_argument("--skymap", required=True)
    sp.add_argument("--catalog", required=True, help="catalog snapshot directory")
    sp.add_argument("--profile", action="append", help="telescope profile id or YAML path (repeatable)")
    sp.add_argument("--preset", choices=["bns", "nsbh", "bbh", "custom"])
    sp.add_argument("--strategy")
    sp.add_argument("--nights", type=float)
    sp.add_argument("--config", help="YAML PlanConfig")
    sp.add_argument("--t-start", help="ISO UTC planning start (default: now)")
    sp.add_argument("--event-id")
    sp.add_argument("--revision", type=int)
    sp.add_argument("--observations", help="JSON list of completed observations to credit")
    sp.add_argument("--out", required=True)
    sp.add_argument("-v", "--verbose", action="store_true")
    sp.set_defaults(fn=_cmd_plan)

    sp = sub.add_parser("ingest-catalog", help="build a versioned catalog snapshot")
    sp.add_argument("kind", choices=["glade+", "generic"])
    sp.add_argument("path")
    sp.add_argument("--out", required=True)
    sp.add_argument("--name")
    sp.add_argument("--upstream-version", required=True)
    sp.add_argument("--source-url")
    sp.add_argument("--retrieved")
    sp.add_argument("--citation")
    sp.add_argument("--mapping", help="column-map YAML (generic)")
    sp.add_argument("--footprint", default="allsky", choices=["allsky", "occupied"])
    sp.add_argument("--ebv-map", help="HEALPix FITS E(B-V) map (e.g. SFD)")
    sp.set_defaults(fn=_cmd_ingest)

    sp = sub.add_parser("validate-schedule", help="run the independent feasibility validator")
    sp.add_argument("schedule")
    sp.add_argument("--profile")
    sp.add_argument("-v", "--verbose", action="store_true")
    sp.set_defaults(fn=_cmd_validate)

    sp = sub.add_parser("export", help="export a schedule through a scheduler adapter")
    sp.add_argument("schedule")
    sp.add_argument("--adapter", required=True)
    sp.add_argument("--profile")
    sp.add_argument("--skymap")
    sp.add_argument("--graceid")
    sp.add_argument("--out")
    sp.set_defaults(fn=_cmd_export)

    sp = sub.add_parser("rerun", help="rerun from a manifest and diff outputs")
    sp.add_argument("manifest")
    sp.add_argument("--skymap", required=True)
    sp.add_argument("--catalog", required=True)
    sp.set_defaults(fn=_cmd_rerun)

    sp = sub.add_parser("compare", help="compare two groups' plan outputs (independent validation)")
    sp.add_argument("dir_a")
    sp.add_argument("dir_b")
    sp.add_argument("--top-k", type=int, default=100)
    sp.add_argument("--out", help="write the comparison report as JSON")
    sp.set_defaults(fn=_cmd_compare)

    sp = sub.add_parser("profiles", help="list or check telescope profiles")
    sp.add_argument("--check", help="profile id or YAML path to validate")
    sp.set_defaults(fn=_cmd_profiles)

    sp = sub.add_parser("keys", help="manage API keys")
    sp.add_argument("action", choices=["create", "revoke"])
    sp.add_argument("--name", default="cli")
    sp.add_argument("--scopes", default="read,plan")
    sp.add_argument("--prefix")
    sp.set_defaults(fn=_cmd_keys)

    sp = sub.add_parser("serve", help="run the REST API")
    sp.add_argument("--host", default="0.0.0.0")
    sp.add_argument("--port", type=int, default=8000)
    sp.add_argument("--behind-proxy", action="store_true")
    sp.set_defaults(fn=_cmd_serve)

    sp = sub.add_parser("worker", help="run a compute worker")
    sp.set_defaults(fn=_cmd_worker)

    sp = sub.add_parser("listen", help="run the alert listener")
    sp.add_argument("--source", choices=["gcn", "scimma", "dir"], required=True)
    sp.add_argument("--path")
    sp.set_defaults(fn=_cmd_listen)

    sp = sub.add_parser("replay", help="replay a directory of LVK JSON alerts through the revision logic")
    sp.add_argument("path")
    sp.set_defaults(fn=_cmd_replay)

    sp = sub.add_parser("benchmark", help="generate the benchmark/validation report")
    sp.add_argument("--data", default="tests/data")
    sp.add_argument("--out", default="benchmark-report")
    sp.add_argument("--injections", type=int, default=200)
    sp.add_argument("--quick", action="store_true")
    sp.set_defaults(fn=_cmd_benchmark)

    a = p.parse_args(argv)
    from tessera.core.errors import TesseraError

    try:
        return int(a.fn(a))
    except TesseraError as exc:
        print(json.dumps(exc.to_problem(), indent=2), file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

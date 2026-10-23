"""Compare two plan output directories (independent operational validation, spec 11.3).

Two groups run the same inputs and exchange their output directories; this
reports whether the inputs really were identical, whether outputs agree
exactly or within the documented tolerance, and the likely causes of any
difference (software versions, IERS table, ephemeris, platform, inputs).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from tessera.core.manifest import TOLERANCE


def _load(d: Path) -> dict[str, Any]:
    man = json.loads((d / "manifest.json").read_text())
    gal = json.loads((d / "galaxies.json").read_text())["rows"]
    sch = {p.stem: json.loads(p.read_text()) for p in sorted((d / "schedules").glob("*.json"))} \
        if (d / "schedules").is_dir() else {}
    return {"manifest": man, "galaxies": gal, "schedules": sch}


def _flat_diff(a: dict, b: dict, prefix: str = "") -> list[str]:
    out = []
    for k in sorted(set(a) | set(b)):
        va, vb = a.get(k), b.get(k)
        if isinstance(va, dict) and isinstance(vb, dict):
            out += _flat_diff(va, vb, f"{prefix}{k}.")
        elif va != vb:
            out.append(f"{prefix}{k}")
    return out


def compare_plans(dir_a: str | Path, dir_b: str | Path, *, rel_tol: float = TOLERANCE["scores_rel"],
                  time_tol_s: float = 1.0, top_k: int = 100) -> dict[str, Any]:
    A, B = _load(Path(dir_a)), _load(Path(dir_b))
    ma, mb = A["manifest"], B["manifest"]
    input_diff = _flat_diff(ma["identity"], mb["identity"])
    output_diff = sorted(k for k in set(ma["outputs"]) | set(mb["outputs"])
                         if ma["outputs"].get(k) != mb["outputs"].get(k))
    env_diff = _flat_diff({"versions": ma.get("versions", {}), "environment": ma.get("environment", {})},
                          {"versions": mb.get("versions", {}), "environment": mb.get("environment", {})})

    ga = [g for g in A["galaxies"] if g.get("tier") == 1][:top_k]
    gb = [g for g in B["galaxies"] if g.get("tier") == 1][:top_k]
    ids_a, ids_b = [g["gal_id"] for g in ga], [g["gal_id"] for g in gb]
    sb = {g["gal_id"]: g["S"] for g in B["galaxies"]}
    rel = [abs(g["S"] / sb[g["gal_id"]] - 1) for g in ga if g["gal_id"] in sb and sb[g["gal_id"]]]
    swaps = []
    for i, (x, y) in enumerate(zip(ids_a, ids_b)):
        if x != y:
            sx, sy = next(g["S"] for g in ga if g["gal_id"] == x), sb.get(y)
            tie = sy is not None and abs(sx / sy - 1) <= rel_tol
            swaps.append({"rank": i + 1, "a": x, "b": y, "tie_within_tolerance": tie})
    galaxies = {
        "top_k": top_k,
        "same_ids": set(ids_a) == set(ids_b),
        "same_order": ids_a == ids_b,
        "order_differences_beyond_ties": sum(1 for s in swaps if not s["tie_within_tolerance"]),
        "max_rel_score_diff": float(max(rel)) if rel else 0.0,
        "first_differences": swaps[:10],
    }
    from tessera.core.timeutil import to_time

    schedules = {}
    for key in sorted(set(A["schedules"]) | set(B["schedules"])):
        va = A["schedules"].get(key, {}).get("visits", [])
        vb = B["schedules"].get(key, {}).get("visits", [])
        same_fields = [v["field_id"] for v in va] == [v["field_id"] for v in vb]
        dt = [abs(to_time(x["start"]).jd - to_time(y["start"]).jd) * 86400 for x, y in zip(va, vb)
              if x["field_id"] == y["field_id"]]
        schedules[key] = {"n_visits": [len(va), len(vb)], "same_field_sequence": same_fields,
                          "max_start_diff_s": float(np.max(dt)) if dt else 0.0}

    if not output_diff:
        verdict = "identical"
    elif (not input_diff and galaxies["same_ids"] and galaxies["order_differences_beyond_ties"] == 0
          and galaxies["max_rel_score_diff"] <= rel_tol
          and all(s["same_field_sequence"] and s["max_start_diff_s"] <= time_tol_s for s in schedules.values())):
        verdict = "within_tolerance"
    else:
        verdict = "different"
    causes = []
    if input_diff:
        causes.append(f"inputs differ: {', '.join(input_diff)}")
    for k in env_diff:
        if k.startswith(("versions.", "environment.iers", "environment.ephemeris")):
            causes.append(f"environment differs: {k} ({_get(ma, k)} vs {_get(mb, k)})")
    if verdict == "different" and not causes:
        causes.append("same inputs and software: suspect platform floating-point behaviour; "
                      "rerun both sides in the published container")
    return {"verdict": verdict, "tolerance": {"scores_rel": rel_tol, "times_s": time_tol_s},
            "plan_keys": [ma["plan_key"], mb["plan_key"]], "input_differences": input_diff,
            "output_files_differing": output_diff, "environment_differences": env_diff,
            "galaxies": galaxies, "schedules": schedules, "likely_causes": causes}


def _get(m: dict, dotted: str):
    cur: Any = m
    for part in dotted.split("."):
        cur = cur.get(part) if isinstance(cur, dict) else None
    return cur

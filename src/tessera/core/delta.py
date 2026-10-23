"""Revision delta reports (spec 8.4): how a new localisation changed the plan."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from scipy.stats import kendalltau


def _load(plan_dir: Path) -> tuple[list[str], dict[str, Any], dict[str, set[str]]]:
    gal = json.loads((plan_dir / "galaxies.json").read_text())
    ids = [r["gal_id"] for r in gal["rows"] if r.get("tier") == 1]
    summary = json.loads((plan_dir / "summary.json").read_text())
    fields: dict[str, set[str]] = {}
    fdir = plan_dir / "fields"
    if fdir.is_dir():
        for f in fdir.glob("*.json"):
            doc = json.loads(f.read_text())
            fields[f.stem] = {r["field_id"] for r in doc["rows"] if r.get("status") == "scheduled"}
    return ids, summary, fields


def revision_delta(old_dir: str | Path, new_dir: str | Path, k: int = 100) -> dict[str, Any]:
    a_ids, a_sum, a_f = _load(Path(old_dir))
    b_ids, b_sum, b_f = _load(Path(new_dir))
    top_a, top_b = a_ids[:k], b_ids[:k]
    common = [g for g in top_a if g in set(top_b)]
    tau = None
    if len(common) >= 3:
        ra = {g: i for i, g in enumerate(top_a)}
        rb = {g: i for i, g in enumerate(top_b)}
        tau = float(kendalltau([ra[g] for g in common], [rb[g] for g in common]).statistic)
    fields = {}
    for prof in sorted(set(a_f) | set(b_f)):
        fa, fb = a_f.get(prof, set()), b_f.get(prof, set())
        fields[prof] = {"added": sorted(fb - fa), "dropped": sorted(fa - fb), "kept": len(fa & fb)}

    def pick(sm: dict[str, Any]) -> dict[str, Any]:
        out = {"area_90_deg2": sm["skymap"].get("area_90_deg2"), "distmean": sm["skymap"].get("distmean"),
               "budget": sm.get("budget")}
        out["schedules"] = {p: {k2: v["summary"].get(k2) for k2 in ("cum_p_sky", "cum_expected_detections",
                                                                      "cum_score_fraction_in_field", "cum_p_host_cat")}
                            for p, v in (sm.get("schedules") or {}).items()}
        return out

    return {
        "top_k": k,
        "top_k_overlap": len(common) / max(1, min(len(top_a), len(top_b))) if top_a and top_b else None,
        "kendall_tau_common": tau,
        "fields": fields,
        "before": pick(a_sum),
        "after": pick(b_sum),
        "new_top10": [g for g in b_ids[:10] if g not in set(a_ids[:10])],
    }

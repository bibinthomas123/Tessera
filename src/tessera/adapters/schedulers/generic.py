"""Reference adapter: the core schedule as JSON, CSV or ECSV."""

from __future__ import annotations

import csv
import io
import json
from typing import Any

from tessera.adapters.schedulers.base import BaseAdapter, Export, Issue
from tessera.core.profile import TelescopeProfile

CSV_COLUMNS = ["seq", "field_id", "origin", "ra", "dec", "pa_deg", "filter", "n_exp", "t_exp_s", "start", "end",
               "m_lim_pred", "gain", "gain_rate_per_h", "airmass_max", "moon_sep_deg", "decision_rule", "feasible"]


class GenericAdapter(BaseAdapter):
    name = "generic"
    version = "1"
    media_type = "application/json"

    def export(self, schedule: dict[str, Any], profile: TelescopeProfile, *, fmt: str = "json",
               allow_infeasible: bool = False, **_: Any) -> Export:
        self.check_exportable(schedule, allow_infeasible=allow_infeasible)
        if fmt == "json":
            data = (json.dumps(schedule, indent=2, sort_keys=True) + "\n").encode()
            return Export(self.name, self.version, {"schedule.json": data}, "application/json")
        rows = []
        for v in schedule["visits"]:
            rows.append({**{k: v.get(k) for k in CSV_COLUMNS},
                         "airmass_max": (v.get("airmass") or {}).get("max"),
                         "decision_rule": v.get("decision", {}).get("rule")})
        if fmt == "csv":
            buf = io.StringIO()
            buf.write(f"# {self.header_comment(schedule)}\n")
            w = csv.DictWriter(buf, fieldnames=CSV_COLUMNS, lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
            return Export(self.name, self.version, {"schedule.csv": buf.getvalue().encode()}, "text/csv")
        if fmt == "ecsv":
            from astropy.table import Table

            t = Table(rows=[[r[c] for c in CSV_COLUMNS] for r in rows], names=CSV_COLUMNS) if rows else Table(
                names=CSV_COLUMNS)
            t.meta["tessera"] = self.header_comment(schedule)
            buf = io.StringIO()
            t.write(buf, format="ascii.ecsv")
            return Export(self.name, self.version, {"schedule.ecsv": buf.getvalue().encode()}, "text/x-ecsv")
        raise ValueError(f"unknown format {fmt!r}")

    def validate(self, export: Export) -> list[Issue]:
        issues = []
        if "schedule.json" in export.files:
            doc = json.loads(export.files["schedule.json"])
            if doc.get("schema") != "tessera.schedule/1":
                issues.append(Issue("error", "schema field missing"))
        return issues

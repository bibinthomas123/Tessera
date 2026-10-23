"""Las Cumbres Observatory request-group JSON (Observation Portal API).

Schema checked against the open-source Observation Portal
(github.com/observatorycontrolsystem/observation-portal, requestgroups/models.py and
serializers.py): operator SINGLE|MANY, observation_type NORMAL|RAPID_RESPONSE|TIME_CRITICAL|
DIRECT|REAL_TIME, ipp_value within the portal's configured range (0.5-2.0 by default),
names <= 50 characters, ICRS ra/dec in decimal degrees, max_airmass 1-25,
min_lunar_distance 0-180 deg, acquisition_config and guiding_config always present.
Instrument types, modes and filter codes come from LCO's configdb and are not verifiable
offline; the endpoint path and token header follow LCO's client practice and must be
confirmed before enabling submission.

LCO's own scheduler chooses the exact start inside each request window, so
each visit becomes a request whose window runs from the validated start to
the end of the validated feasibility window; the validator has checked the
whole window.  LCO-side constraints (max airmass, Moon distance) are copied
from the profile so the network enforces them too.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from typing import Any

from tessera.adapters.schedulers.base import (
    Approval,
    BaseAdapter,
    Export,
    Issue,
    SubmitReceipt,
    check_submission,
)
from tessera.core.errors import AdapterError
from tessera.core.profile import TelescopeProfile

DEFAULT_URL = "https://observe.lco.global/api/requestgroups/"


class LCOAdapter(BaseAdapter):
    name = "lco"
    version = "1"
    requires = {"telescope_class", "instrument_type", "proposal"}
    media_type = "application/json"
    can_submit = True

    def export(self, schedule: dict[str, Any], profile: TelescopeProfile, *, allow_infeasible: bool = False,
               **opts: Any) -> Export:
        self.check_exportable(schedule, allow_infeasible=allow_infeasible)
        o = self.options(profile, **opts)
        fmap = o.get("filter_map", {"g": "gp", "r": "rp", "i": "ip", "z": "zs"})
        ms = profile.constraints.moon_sep()
        min_moon = float(ms.constant) if ms.constant else float(max((y for _, y in ms.by_illumination or []),
                                                                     default=0))
        requests = []
        for v in schedule["visits"]:
            requests.append({
                "acceptability_threshold": float(o.get("acceptability_threshold", 90)),
                "configuration_repeats": 1,
                "optimization_type": o.get("optimization_type", "TIME"),
                "configurations": [{
                    "type": "EXPOSE",
                    "instrument_type": o["instrument_type"],
                    "target": {"name": v["field_id"][:50], "type": "ICRS", "ra": v["ra"], "dec": v["dec"]},
                    "constraints": {"max_airmass": profile.constraints.airmass_max, "min_lunar_distance": min_moon},
                    "acquisition_config": {"mode": o.get("acquisition_mode", "OFF")},
                    "guiding_config": ({"mode": o["guiding_mode"]} if o.get("guiding_mode") else {}),
                    "instrument_configs": [{"exposure_time": v["t_exp_s"], "exposure_count": v["n_exp"],
                                            "optical_elements": {"filter": fmap.get(v["filter"], v["filter"])}}],
                }],
                "windows": [{"start": v["start"], "end": v["window"][1]}],
                "location": {"telescope_class": o["telescope_class"]},
            })
        doc = {
            "name": f"tessera-{schedule.get('schedule_id', 'plan')}"[:50],
            "proposal": o["proposal"],
            "ipp_value": float(o.get("ipp_value", 1.0)),
            "operator": "MANY" if len(requests) > 1 else "SINGLE",
            "observation_type": o.get("observation_type", "NORMAL"),
            "requests": requests,
        }
        data = (json.dumps(doc, indent=2, sort_keys=True) + "\n").encode()
        return Export(self.name, self.version, {"requestgroup.json": data}, self.media_type)

    def validate(self, export: Export) -> list[Issue]:
        doc = json.loads(export.files["requestgroup.json"])
        issues = []
        for key in ("name", "proposal", "ipp_value", "operator", "observation_type", "requests"):
            if key not in doc:
                issues.append(Issue("error", f"missing {key}"))
        if doc.get("operator") not in ("SINGLE", "MANY"):
            issues.append(Issue("error", "operator must be SINGLE or MANY"))
        if doc.get("operator") == "SINGLE" and len(doc.get("requests", [])) != 1:
            issues.append(Issue("error", "operator SINGLE requires exactly one request"))
        if doc.get("observation_type") not in ("NORMAL", "RAPID_RESPONSE", "TIME_CRITICAL", "DIRECT", "REAL_TIME"):
            issues.append(Issue("error", "invalid observation_type"))
        if not 0.5 <= float(doc.get("ipp_value", 0)) <= 2.0:
            issues.append(Issue("error", "ipp_value outside the portal default range 0.5-2.0"))
        if len(doc.get("name", "")) > 50:
            issues.append(Issue("error", "name longer than 50 characters"))
        for i, r in enumerate(doc.get("requests", [])):
            if not r.get("windows"):
                issues.append(Issue("error", f"request {i} has no windows"))
            if not 0 <= r.get("acceptability_threshold", 90) <= 100:
                issues.append(Issue("error", f"request {i}: acceptability_threshold must be 0-100"))
            for w in r.get("windows", []):
                if not w.get("end") or (w.get("start") and w["start"] >= w["end"]):
                    issues.append(Issue("error", f"request {i}: window end must be after start"))
            for c in r.get("configurations", []):
                tg = c.get("target", {})
                if tg.get("type") != "ICRS":
                    issues.append(Issue("error", f"request {i}: target type must be ICRS"))
                if not (0 <= tg.get("ra", -1) <= 360 and -90 <= tg.get("dec", 99) <= 90):
                    issues.append(Issue("error", f"request {i}: ra/dec out of range"))
                if len(tg.get("name", "")) > 50:
                    issues.append(Issue("error", f"request {i}: target name longer than 50 characters"))
                cons = c.get("constraints", {})
                if not 1.0 <= cons.get("max_airmass", 1.6) <= 25.0:
                    issues.append(Issue("error", f"request {i}: max_airmass must be 1-25"))
                if not 0.0 <= cons.get("min_lunar_distance", 30) <= 180.0:
                    issues.append(Issue("error", f"request {i}: min_lunar_distance must be 0-180"))
                for key in ("acquisition_config", "guiding_config"):
                    if key not in c:
                        issues.append(Issue("error", f"request {i}: {key} is required"))
                for ic in c.get("instrument_configs", []):
                    if not ic.get("exposure_time") or ic["exposure_time"] <= 0:
                        issues.append(Issue("error", f"request {i}: exposure_time must be > 0"))
        return issues

    def submit(self, export: Export, approval: Approval | None, schedule: dict[str, Any],
               profile: TelescopeProfile) -> SubmitReceipt:
        check_submission(approval, schedule, profile)
        token = os.environ.get("TESSERA_LCO_TOKEN")
        if not token:
            raise AdapterError("TESSERA_LCO_TOKEN is not set", code="ADAPTER_CREDENTIALS_MISSING")
        import httpx

        url = profile.adapter_options.get("lco", {}).get("url", DEFAULT_URL)
        r = httpx.post(url, content=export.files["requestgroup.json"], timeout=30,
                       headers={"Authorization": f"Token {token}", "Content-Type": "application/json"})
        if r.status_code >= 300:
            raise AdapterError(f"LCO rejected the request group: HTTP {r.status_code}", code="ADAPTER_REMOTE_ERROR",
                               details={"body": r.text[:2000]})
        body = r.json()
        return SubmitReceipt(self.name, str(body.get("id")), datetime.now(timezone.utc).isoformat(), {"state": body.get("state")})

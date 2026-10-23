"""Gravitational Wave Treasure Map: planned pointings for cross-group coordination.

Exports the payload for ``POST /api/v1/pointings`` with ``status: planned``.
Request format checked against the official Python client
(github.com/TheTreasureMap/gwtm_api: core/_client_factory.py, pointing.py,
core/apimodels.py): base URL https://treasuremap.space/api/v1, batch body
``{"graceid", "pointings": [...]}``, token sent as the ``api_token`` HTTP header,
per-pointing ra/dec, instrumentid (int), band (bandpass enum), depth,
depth_unit (ab_mag|vega_mag|flux_erg|flux_jy), ISO-8601 time, pos_angle, status
(planned|completed|cancelled).  The token is never written into the export.
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

DEFAULT_URL = "https://treasuremap.space/api/v1/pointings"
# bandpass enum of gwtm_api core/apimodels.py
BANDS = {"U", "B", "V", "R", "I", "J", "H", "K", "u", "g", "r", "i", "z", "UVW1", "UVW2", "UVM2", "XRT", "clear",
         "open", "UHF", "VHF", "L", "S", "C", "X", "other", "TESS", "BAT", "HESS", "WISEL", "q"}


class TreasureMapAdapter(BaseAdapter):
    name = "treasuremap"
    version = "1"
    requires = {"instrumentid"}
    media_type = "application/json"
    can_submit = True

    def export(self, schedule: dict[str, Any], profile: TelescopeProfile, *, graceid: str | None = None,
               allow_infeasible: bool = False, **opts: Any) -> Export:
        self.check_exportable(schedule, allow_infeasible=allow_infeasible)
        o = self.options(profile, **opts)
        bmap = o.get("band_map", {})
        pointings = [{
            "ra": v["ra"],
            "dec": v["dec"],
            "instrumentid": int(o["instrumentid"]),
            "band": bmap.get(v["filter"], v["filter"]),
            "depth": v.get("m_lim_pred"),
            "depth_unit": "ab_mag",
            "time": v["start"].rstrip("Z"),        # ISO-8601 UTC
            "pos_angle": v.get("pa_deg", 0.0),
            "status": "planned",
        } for v in schedule["visits"]]
        doc = {"graceid": graceid or o.get("graceid") or "UNKNOWN", "pointings": pointings}
        data = (json.dumps(doc, indent=2, sort_keys=True) + "\n").encode()
        return Export(self.name, self.version, {"treasuremap_pointings.json": data}, self.media_type)

    def validate(self, export: Export) -> list[Issue]:
        doc = json.loads(export.files["treasuremap_pointings.json"])
        issues = []
        if not doc.get("graceid") or doc["graceid"] == "UNKNOWN":
            issues.append(Issue("warning", "graceid not set; pass graceid= when exporting"))
        for i, p in enumerate(doc.get("pointings", [])):
            if not (0 <= p["ra"] < 360 and -90 <= p["dec"] <= 90):
                issues.append(Issue("error", f"pointing {i}: coordinates out of range"))
            if p["status"] != "planned":
                issues.append(Issue("error", f"pointing {i}: status must be 'planned'"))
            if p["band"] not in BANDS:
                issues.append(Issue("error", f"pointing {i}: band {p['band']!r} is not a Treasure Map bandpass; "
                                             "set adapter_options.treasuremap.band_map"))
            if p["depth_unit"] not in ("ab_mag", "vega_mag", "flux_erg", "flux_jy"):
                issues.append(Issue("error", f"pointing {i}: invalid depth_unit"))
            if not isinstance(p["instrumentid"], int):
                issues.append(Issue("error", f"pointing {i}: instrumentid must be an integer"))
        return issues

    def submit(self, export: Export, approval: Approval | None, schedule: dict[str, Any],
               profile: TelescopeProfile) -> SubmitReceipt:
        check_submission(approval, schedule, profile)
        token = os.environ.get("TESSERA_TREASUREMAP_TOKEN")
        if not token:
            raise AdapterError("TESSERA_TREASUREMAP_TOKEN is not set", code="ADAPTER_CREDENTIALS_MISSING")
        import httpx

        doc = json.loads(export.files["treasuremap_pointings.json"])
        url = profile.adapter_options.get("treasuremap", {}).get("url", DEFAULT_URL)
        r = httpx.post(url, json=doc, headers={"api_token": token}, timeout=30)
        if r.status_code >= 300:
            raise AdapterError(f"Treasure Map rejected the pointings: HTTP {r.status_code}",
                               code="ADAPTER_REMOTE_ERROR", details={"body": r.text[:2000]})
        return SubmitReceipt(self.name, None, datetime.now(timezone.utc).isoformat(), {"response": r.text[:500]})

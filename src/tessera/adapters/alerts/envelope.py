"""Normalised alert envelope for LVK public alerts (spec 8.1, 8.2).

Parses the LVK JSON alert format distributed on the ``igwn.gwalert`` topic
(GCN Kafka and SCiMMA Hopskotch).  Alert types: EARLY_WARNING, PRELIMINARY,
INITIAL, UPDATE, RETRACTION.  The multi-order sky map arrives base64-encoded
in ``event.skymap``.
"""

from __future__ import annotations

import base64
import json
from dataclasses import dataclass, field
from typing import Any

from tessera.core.errors import TesseraError
from tessera.core.hashing import sha256_bytes

ALERT_TYPES = ("EARLY_WARNING", "PRELIMINARY", "INITIAL", "UPDATE", "RETRACTION")


class AlertError(TesseraError):
    code = "ALERT_INVALID"
    http_status = 422


@dataclass
class AlertEnvelope:
    superevent_id: str
    alert_type: str
    time_created: str
    source: str
    skymap_bytes: bytes | None = None
    skymap_filename: str | None = None
    significant: bool | None = None
    far: float | None = None
    group: str | None = None
    pipeline: str | None = None
    search: str | None = None
    instruments: list[str] = field(default_factory=list)
    properties: dict[str, float] = field(default_factory=dict)       # HasNS, HasRemnant, HasMassGap
    classification: dict[str, float] = field(default_factory=dict)   # BNS, NSBH, BBH, Terrestrial
    urls: dict[str, str] = field(default_factory=dict)
    raw: bytes = b""

    @property
    def is_mock(self) -> bool:
        return self.superevent_id.startswith(("MS", "TS"))  # mock (playground) and test superevents

    @property
    def skymap_sha256(self) -> str | None:
        return sha256_bytes(self.skymap_bytes) if self.skymap_bytes else None

    @property
    def alert_key(self) -> str:
        return sha256_bytes("|".join([self.superevent_id, self.alert_type, self.time_created,
                                      self.skymap_sha256 or "-"]).encode())

    def metadata(self) -> dict[str, Any]:
        return {"superevent_id": self.superevent_id, "alert_type": self.alert_type, "time_created": self.time_created,
                "source": self.source, "significant": self.significant, "far": self.far, "group": self.group,
                "pipeline": self.pipeline, "search": self.search, "instruments": self.instruments,
                "properties": self.properties, "classification": self.classification, "urls": self.urls,
                "skymap_sha256": self.skymap_sha256, "skymap_filename": self.skymap_filename}


def parse_lvk_json(payload: bytes | str | dict[str, Any], source: str = "unknown") -> AlertEnvelope:
    if isinstance(payload, dict):
        doc, raw = payload, json.dumps(payload, sort_keys=True).encode()
    else:
        raw = payload.encode() if isinstance(payload, str) else bytes(payload)
        try:
            doc = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise AlertError(f"alert is not valid JSON: {exc}") from exc
    try:
        sid = str(doc["superevent_id"])
        atype = str(doc["alert_type"]).upper()
        tc = str(doc["time_created"])
    except KeyError as exc:
        raise AlertError(f"alert missing required field {exc}") from exc
    if atype not in ALERT_TYPES:
        raise AlertError(f"unknown alert_type {atype!r}")
    env = AlertEnvelope(superevent_id=sid, alert_type=atype, time_created=tc, source=source, raw=raw,
                        urls=dict(doc.get("urls") or {}))
    ev = doc.get("event")
    if ev:
        env.significant = ev.get("significant")
        env.far = ev.get("far")
        env.group = ev.get("group")
        env.pipeline = ev.get("pipeline")
        env.search = ev.get("search")
        env.instruments = list(ev.get("instruments") or [])
        env.properties = {k: float(v) for k, v in (ev.get("properties") or {}).items() if v is not None}
        env.classification = {k: float(v) for k, v in (ev.get("classification") or {}).items() if v is not None}
        sk = ev.get("skymap")
        if sk:
            try:
                env.skymap_bytes = base64.b64decode(sk, validate=False)
            except Exception as exc:
                raise AlertError(f"cannot decode event.skymap: {exc}") from exc
            env.skymap_filename = f"{sid}-{atype.lower()}.multiorder.fits"
    if atype != "RETRACTION" and not ev:
        raise AlertError("non-retraction alert without an event section")
    return env


def make_alert(superevent_id: str, alert_type: str, time_created: str, skymap_bytes: bytes | None = None,
               **event: Any) -> dict[str, Any]:
    """Build an LVK-format alert dict (tests, replay, local injection)."""
    doc: dict[str, Any] = {"superevent_id": superevent_id, "alert_type": alert_type, "time_created": time_created,
                           "urls": {"gracedb": f"https://gracedb.ligo.org/superevents/{superevent_id}/view/"}}
    if alert_type == "RETRACTION":
        doc["event"] = None
    else:
        doc["event"] = {"time": event.get("time", time_created), "far": event.get("far", 1e-10),
                        "significant": event.get("significant", True), "instruments": event.get("instruments", ["H1", "L1"]),
                        "group": event.get("group", "CBC"), "pipeline": event.get("pipeline", "gstlal"),
                        "search": event.get("search", "AllSky"),
                        "properties": event.get("properties", {"HasNS": 1.0, "HasRemnant": 1.0}),
                        "classification": event.get("classification", {"BNS": 0.99}),
                        "skymap": base64.b64encode(skymap_bytes).decode() if skymap_bytes else None}
    return doc

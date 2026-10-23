"""Alert listener process: consumes alert sources and applies revision logic (spec 8)."""

from __future__ import annotations

import os
from collections.abc import Iterable

from tessera.adapters.alerts.envelope import AlertEnvelope
from tessera.service import obs
from tessera.service.db import session_factory
from tessera.service.events import AlertOutcome, TriggerPolicy, handle_alert
from tessera.service.settings import Settings, get_settings
from tessera.service.store import ArtifactStore


def consume(envelopes: Iterable[AlertEnvelope], settings: Settings | None = None) -> list[AlertOutcome]:
    settings = settings or get_settings()
    S = session_factory(settings.database_url)
    store = ArtifactStore(settings.artifacts_dir)
    policy = TriggerPolicy.load(settings.trigger_policy, settings.mode)
    outcomes = []
    for env in envelopes:
        with S() as s:
            try:
                out = handle_alert(s, settings, store, env, policy)
            except Exception as exc:
                obs.log.error("alert_failed", superevent_id=env.superevent_id, alert_type=env.alert_type,
                              error=str(exc))
                obs.ALERTS.labels(source=env.source.split(":")[0], type=env.alert_type, outcome="error").inc()
                continue
        obs.ALERTS.labels(source=env.source.split(":")[0], type=env.alert_type, outcome=out.outcome).inc()
        obs.log.info("alert", superevent_id=env.superevent_id, alert_type=env.alert_type, outcome=out.outcome,
                     revision_id=out.revision_id, jobs=out.jobs, detail=out.detail)
        outcomes.append(out)
    return outcomes


def run(source: str, path: str | None = None) -> None:  # pragma: no cover - long-running
    obs.configure_logging(os.environ.get("TESSERA_LOG_LEVEL", "INFO"))
    from tessera.adapters.alerts import sources

    if source == "gcn":
        it = sources.gcn_kafka()
    elif source == "scimma":
        it = sources.scimma_hop()
    elif source == "dir":
        it = sources.directory(path or "./alerts")
    else:
        raise SystemExit(f"unknown source {source!r}")
    obs.log.info("listener_start", source=source)
    consume(it)

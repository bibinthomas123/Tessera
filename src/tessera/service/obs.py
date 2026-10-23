"""Observability: structured JSON logs and Prometheus metrics (spec 10.4)."""

from __future__ import annotations

import logging
import sys

import structlog
from prometheus_client import Counter, Gauge, Histogram

STAGE_SECONDS = Histogram("tessera_stage_seconds", "Duration of engine stages", ["stage"],
                          buckets=(0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 30, 60, 120, 300))
JOB_SECONDS = Histogram("tessera_job_seconds", "End-to-end job duration",
                        buckets=(1, 2, 5, 10, 20, 30, 60, 120, 300, 600))
ALERT_TO_PROPOSED = Histogram("tessera_alert_to_proposed_seconds", "Alert received to proposed schedule",
                              buckets=(5, 10, 20, 30, 45, 60, 90, 120, 180, 300, 600))
JOBS = Counter("tessera_jobs_total", "Jobs by final state", ["state"])
QUEUE_DEPTH = Gauge("tessera_queue_depth", "Queued and running jobs")
ALERTS = Counter("tessera_alerts_total", "Alerts received", ["source", "type", "outcome"])
VALIDATOR_VIOLATIONS = Counter("tessera_validator_violations_total", "Validator violations", ["constraint"])


def configure_logging(level: str = "INFO") -> None:
    logging.basicConfig(stream=sys.stdout, level=getattr(logging, level.upper(), logging.INFO), format="%(message)s")
    structlog.configure(
        processors=[
            structlog.contextvars.merge_contextvars,
            structlog.processors.add_log_level,
            structlog.processors.TimeStamper(fmt="iso", utc=True),
            structlog.processors.StackInfoRenderer(),
            structlog.processors.format_exc_info,
            structlog.processors.JSONRenderer(sort_keys=True),
        ],
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        cache_logger_on_first_use=True,
    )


log = structlog.get_logger("tessera")

"""Compute worker (spec 10.2, 10.3).

The parent process claims jobs, holds the lease (heartbeat) and supervises a
child process per job.  The child has a wall-time limit and (on Linux) an
address-space limit.  Cancellation and timeouts terminate the child; crashes
are retried with exponential backoff up to ``TESSERA_JOB_MAX_ATTEMPTS``; failures
caused by the inputs are never retried.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import socket
import time
import traceback
from typing import Any

from tessera.core.errors import TesseraError
from tessera.service import obs
from tessera.service.db import Job, Plan, session_factory
from tessera.service.jobs import claim, execute, finish, heartbeat
from tessera.service.settings import Settings, get_settings

_INPUT_ERRORS_PREFIXES = ("SKYMAP_", "CATALOG_", "PROFILE_", "CONFIG_", "HOST_MODEL", "PLAN_KEY_MISMATCH", "WEIGHT_")


def _child(job_id: str, env: dict[str, str], conn) -> None:
    os.environ.update(env)
    try:
        import resource

        lim = int(env.get("TESSERA_JOB_MEMORY_BYTES", "0"))
        if lim and hasattr(resource, "RLIMIT_AS") and os.uname().sysname == "Linux":
            resource.setrlimit(resource.RLIMIT_AS, (lim, lim))
    except Exception:
        pass
    settings = Settings()
    S = session_factory(settings.database_url)
    with S() as s:
        try:
            key = execute(s, settings, job_id)
            conn.send({"ok": True, "plan_key": key})
        except TesseraError as exc:
            conn.send({"ok": False, "error": {"code": exc.code, "message": exc.message, "details": exc.details}})
        except MemoryError:
            conn.send({"ok": False, "error": {"code": "JOB_OUT_OF_MEMORY", "message": "memory limit exceeded"}})
        except Exception as exc:  # unexpected: retryable
            conn.send({"ok": False, "retry": True, "error": {"code": "INTERNAL_ERROR", "message": str(exc),
                                                             "traceback": traceback.format_exc()[-4000:]}})


def run_one(settings: Settings, worker_id: str) -> bool:
    """Claim and run a single job.  Returns False when the queue was empty."""
    S = session_factory(settings.database_url)
    with S() as s:
        job = claim(s, worker_id, settings.lease_s)
        if job is None:
            return False
        job_id, attempts = job.id, job.attempts
        obs.log.info("job_start", job_id=job_id, plan_key=job.plan_key, attempt=attempts)
    t0 = time.monotonic()
    ctx = mp.get_context("spawn")
    parent, child = ctx.Pipe(duplex=False)
    env = {k: v for k, v in os.environ.items() if k.startswith("TESSERA_")}
    env["TESSERA_JOB_MEMORY_BYTES"] = str(settings.job_memory_bytes)
    proc = ctx.Process(target=_child, args=(job_id, env, child), daemon=True)
    proc.start()
    child.close()
    result: dict[str, Any] | None = None
    last_hb = time.monotonic()
    state = None
    while True:
        if parent.poll(1.0):
            try:
                result = parent.recv()
            except EOFError:
                result = None
            proc.join(10)
            break
        if not proc.is_alive():
            break
        now = time.monotonic()
        with S() as s:
            if now - last_hb > settings.lease_s / 3:
                heartbeat(s, job_id, worker_id, settings.lease_s)
                last_hb = now
            j = s.get(Job, job_id)
            if j.cancel_requested:
                proc.terminate()
                state = "superseded" if (j.error or {}).get("code") == "REVISION_SUPERSEDED" else "cancelled"
                break
        if now - t0 > settings.job_timeout_s:
            proc.terminate()
            state = "timeout"
            break
    proc.join(5)
    if proc.is_alive():
        proc.kill()
    elapsed = time.monotonic() - t0
    with S() as s:
        if state in ("cancelled", "superseded"):
            finish(s, job_id, state, error={"code": "JOB_CANCELLED" if state == "cancelled" else "REVISION_SUPERSEDED"})
        elif state == "timeout":
            finish(s, job_id, "failed", error={"code": "JOB_TIMEOUT", "message": f"exceeded {settings.job_timeout_s} s"})
        elif result and result.get("ok"):
            finish(s, job_id, "succeeded")
            plan = s.get(Plan, result["plan_key"])
            for stage, sec in (plan.manifest.get("timings_s") or {}).items():
                obs.STAGE_SECONDS.labels(stage=stage).observe(sec)
            from sqlalchemy import select

            from tessera.service.db import Schedule, ScheduleVersion

            for sv in s.scalars(select(ScheduleVersion).join(Schedule, Schedule.id == ScheduleVersion.schedule_id)
                                .where(Schedule.plan_key == plan.plan_key, ScheduleVersion.version == 0)).all():
                for constraint, n in (sv.doc.get("validation", {}).get("violations_by_constraint") or {}).items():
                    obs.VALIDATOR_VIOLATIONS.labels(constraint=constraint).inc(n)
            obs.JOB_SECONDS.observe(elapsed)
            ev = plan.manifest.get("event") or {}
            if ev.get("time_created"):
                from astropy.time import Time

                from tessera.core.timeutil import to_time

                obs.ALERT_TO_PROPOSED.observe(max(0.0, (Time.now() - to_time(ev["time_created"])).sec))
        else:
            err = (result or {}).get("error") or {"code": "WORKER_CRASH", "message": f"exit code {proc.exitcode}"}
            retryable = (result is None or result.get("retry")) and not str(err.get("code", "")).startswith(
                _INPUT_ERRORS_PREFIXES)
            if retryable and attempts < settings.job_max_attempts:
                finish(s, job_id, "queued", error=err, retry_in_s=5.0 * 2 ** (attempts - 1))
                state = "retry"
            else:
                finish(s, job_id, "failed", error=err)
        j = s.get(Job, job_id)
        obs.JOBS.labels(state=j.state).inc()
        obs.log.info("job_end", job_id=job_id, state=j.state, seconds=round(elapsed, 2),
                     error=(j.error or {}).get("code"))
    return True


def run_forever(poll_s: float = 1.0) -> None:  # pragma: no cover - long-running loop
    settings = get_settings()
    obs.configure_logging(os.environ.get("TESSERA_LOG_LEVEL", "INFO"))
    worker_id = f"{socket.gethostname()}:{os.getpid()}"
    try:
        from prometheus_client import start_http_server

        start_http_server(settings.worker_metrics_port)
    except OSError:
        obs.log.warning("metrics_port_unavailable", port=settings.worker_metrics_port)
    obs.log.info("worker_start", worker_id=worker_id)
    while True:
        try:
            if not run_one(settings, worker_id):
                time.sleep(poll_s)
        except Exception as exc:
            obs.log.error("worker_loop_error", error=str(exc))
            time.sleep(poll_s)

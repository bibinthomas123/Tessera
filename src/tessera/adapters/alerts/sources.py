"""Alert sources: GCN Kafka, SCiMMA Hopskotch, GraceDB, local files (spec 8.1).

Each source yields :class:`AlertEnvelope` objects.  Delivery is at-least-once;
idempotency is handled downstream by the alert key (spec 8.2).
"""

from __future__ import annotations

import json
import os
import time
from collections.abc import Iterator
from pathlib import Path
from urllib.parse import urlparse

from tessera.adapters.alerts.envelope import AlertEnvelope, AlertError, parse_lvk_json
from tessera.core.errors import TesseraError

LVK_TOPIC = "igwn.gwalert"
DEFAULT_ALLOWED_HOSTS = ("gracedb.ligo.org", "gracedb-playground.ligo.org", "gracedb-test.ligo.org")


class FetchError(TesseraError):
    code = "SKYMAP_FETCH_FAILED"
    http_status = 422


def allowed_hosts() -> tuple[str, ...]:
    env = os.environ.get("TESSERA_SKYMAP_URL_ALLOWLIST")
    return tuple(h.strip() for h in env.split(",") if h.strip()) if env else DEFAULT_ALLOWED_HOSTS


def fetch_url(url: str, *, max_bytes: int = 100 * 1024**2, timeout: float = 30.0) -> bytes:
    """Fetch a sky map from an allow-listed https host (SSRF protection, spec 10.3)."""
    import httpx

    def check(u: str) -> None:
        p = urlparse(u)
        if p.scheme != "https":
            raise FetchError("only https URLs are accepted", code="SKYMAP_URL_REFUSED")
        if p.hostname not in allowed_hosts():
            raise FetchError(f"host {p.hostname!r} is not on the allowlist", code="SKYMAP_URL_REFUSED",
                             details={"allowed": list(allowed_hosts())})

    check(url)
    with httpx.Client(timeout=timeout, follow_redirects=False) as client:
        for _ in range(5):
            with client.stream("GET", url) as r:
                if r.is_redirect:
                    url = str(r.headers["location"])
                    if url.startswith("/"):
                        p = urlparse(str(r.url))
                        url = f"{p.scheme}://{p.netloc}{url}"
                    check(url)
                    continue
                if r.status_code != 200:
                    raise FetchError(f"HTTP {r.status_code} fetching {url}")
                buf = bytearray()
                for chunk in r.iter_bytes():
                    buf += chunk
                    if len(buf) > max_bytes:
                        raise FetchError("download exceeds size limit", code="SKYMAP_TOO_LARGE")
                return bytes(buf)
    raise FetchError("too many redirects")


def gracedb_skymap_url(superevent_id: str, filename: str = "bayestar.multiorder.fits",
                       host: str = "gracedb.ligo.org") -> str:
    return f"https://{host}/api/superevents/{superevent_id}/files/{filename}"


def gcn_kafka(client_id: str | None = None, client_secret: str | None = None, *, topics=(LVK_TOPIC,),
              domain: str | None = None) -> Iterator[AlertEnvelope]:
    """Consume LVK alerts from GCN Kafka (requires ``gcn-kafka``)."""
    from gcn_kafka import Consumer  # type: ignore

    consumer = Consumer(client_id=client_id or os.environ["GCN_CLIENT_ID"],
                        client_secret=client_secret or os.environ["GCN_CLIENT_SECRET"],
                        domain=domain or os.environ.get("GCN_DOMAIN", "gcn.nasa.gov"))
    consumer.subscribe(list(topics))
    while True:
        for message in consumer.consume(timeout=1):
            if message.error():
                continue
            try:
                yield parse_lvk_json(message.value(), source="gcn")
            except AlertError:
                continue


def scimma_hop(url: str = f"kafka://kafka.scimma.org/{LVK_TOPIC}") -> Iterator[AlertEnvelope]:
    """Consume LVK alerts from SCiMMA Hopskotch (requires ``hop-client`` and hop credentials)."""
    from hop import stream  # type: ignore

    with stream.open(url, "r") as s:
        for msg in s:
            content = getattr(msg, "content", msg)
            if isinstance(content, (bytes, str)):
                payload = content
            else:
                payload = json.dumps(content)
            try:
                yield parse_lvk_json(payload, source="scimma")
            except AlertError:
                continue


def directory(path: str | Path, *, poll_s: float = 2.0, once: bool = False) -> Iterator[AlertEnvelope]:
    """Watch a directory for ``*.json`` alerts (offline testing, replay, air-gapped sites)."""
    seen: set[str] = set()
    p = Path(path)
    while True:
        for f in sorted(p.glob("*.json")):
            if f.name in seen:
                continue
            seen.add(f.name)
            try:
                yield parse_lvk_json(f.read_bytes(), source=f"file:{f.name}")
            except AlertError:
                continue
        if once:
            return
        time.sleep(poll_s)

"""Deterministic serialisation of tables and documents.

Every output is serialised exactly one way per format so that its sha256 in
the manifest is reproducible.
"""

from __future__ import annotations

import io
import json
import math
from typing import Any

import numpy as np
from astropy.table import Table

from tessera.core.hashing import canonical_json

TABLE_FORMATS = ("csv", "ecsv", "json", "parquet", "votable")
MEDIA_TYPES = {
    "csv": "text/csv",
    "ecsv": "text/x-ecsv",
    "json": "application/json",
    "parquet": "application/vnd.apache.parquet",
    "votable": "application/x-votable+xml",
}


def _py(v: Any) -> Any:
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        f = float(v)
        return None if not math.isfinite(f) else f
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, bytes):
        return v.decode()
    return v if not isinstance(v, np.str_) else str(v)


def table_rows(t: Table) -> list[dict[str, Any]]:
    cols = t.colnames
    data = [t[c].tolist() if hasattr(t[c], "tolist") else list(t[c]) for c in cols]
    return [{c: _py(data[j][i]) for j, c in enumerate(cols)} for i in range(len(t))]


def table_to_bytes(t: Table, fmt: str, *, meta: dict[str, Any] | None = None) -> bytes:
    if fmt == "csv":
        buf = io.StringIO()
        t.write(buf, format="ascii.csv", overwrite=True)
        return buf.getvalue().encode()
    if fmt == "ecsv":
        tt = t.copy(copy_data=False)
        if meta:
            tt.meta.update(json.loads(canonical_json(meta)))
        buf = io.StringIO()
        tt.write(buf, format="ascii.ecsv", overwrite=True)
        return buf.getvalue().encode()
    if fmt == "json":
        doc = {"meta": meta or {}, "columns": {c: {"unit": str(t[c].unit) if t[c].unit else None,
                                                   "description": t[c].description} for c in t.colnames},
               "rows": table_rows(t)}
        return canonical_json(doc).encode()
    if fmt == "parquet":
        import pyarrow as pa
        import pyarrow.parquet as pq

        arrays = {c: pa.array(t[c].tolist()) for c in t.colnames}
        md = {b"tessera": canonical_json(meta or {}).encode()}
        tbl = pa.table(arrays).replace_schema_metadata(md)
        buf = io.BytesIO()
        pq.write_table(tbl, buf, compression="zstd")
        return buf.getvalue()
    if fmt == "votable":
        buf = io.BytesIO()
        t.write(buf, format="votable", overwrite=True)
        return buf.getvalue()
    raise ValueError(f"unknown table format {fmt!r}; use one of {TABLE_FORMATS}")


def doc_to_bytes(doc: dict[str, Any]) -> bytes:
    return (json.dumps(json.loads(canonical_json(doc)), indent=2, sort_keys=True) + "\n").encode()

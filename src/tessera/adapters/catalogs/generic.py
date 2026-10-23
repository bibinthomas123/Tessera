"""Generic CSV/Parquet catalog adapter with a declarative column map.

Use this for DESI-based lists, a group's private catalog, or any table that
is not GLADE+.  Example ``columns.yaml``::

    columns:
      gal_id: TARGETID          # required
      ra: TARGET_RA             # required, deg ICRS
      dec: TARGET_DEC           # required
      z: Z
      z_err: ZERR
      mag_W1: W1MAG
    id_prefix: "DESI-DR1:"
    dist_kind: specz            # constant, or {column: SPECTYPE, map: {GALAXY: specz}}
    compute_distance_from_z: true
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pyarrow.csv as pacsv
import pyarrow.parquet as pq
import yaml

from tessera.core.errors import CatalogError


def _read_table(path: Path) -> dict[str, np.ndarray]:
    if path.suffix == ".parquet":
        tbl = pq.read_table(path)
    else:
        tbl = pacsv.read_csv(path)
    out = {}
    for name in tbl.column_names:
        col = tbl.column(name)
        try:
            out[name] = col.to_numpy(zero_copy_only=False)
        except Exception:
            out[name] = np.asarray(col.to_pylist(), dtype=object)
    return out


def read(path: str | Path, mapping: dict[str, Any] | str | Path) -> dict[str, np.ndarray]:
    if not isinstance(mapping, dict):
        mapping = yaml.safe_load(Path(mapping).read_text())
    raw = _read_table(Path(path))
    cmap: dict[str, str] = mapping.get("columns", {})
    for req in ("gal_id", "ra", "dec"):
        if req not in cmap:
            raise CatalogError(f"column map must define {req!r}", code="CATALOG_MAPPING")
    out: dict[str, Any] = {}
    for canon, src in cmap.items():
        if src not in raw:
            raise CatalogError(f"source column {src!r} (for {canon}) not found", code="CATALOG_MAPPING",
                               details={"available": sorted(raw)[:50]})
        out[canon] = raw[src]
    prefix = mapping.get("id_prefix", "")
    out["gal_id"] = np.array([f"{prefix}{v}" for v in out["gal_id"]], dtype=object)
    n = len(out["gal_id"])
    dk = mapping.get("dist_kind")
    if isinstance(dk, str):
        out["dist_kind"] = np.full(n, dk, dtype=object)
    elif isinstance(dk, dict):
        src = raw[dk["column"]]
        m = dk["map"]
        out["dist_kind"] = np.array([m.get(str(v).strip(), "none") for v in src], dtype=object)
    if mapping.get("compute_distance_from_z") and "z" in out and "dist_mpc" not in out:
        from astropy.cosmology import Planck18

        z = np.asarray(out["z"], float)
        d = np.full(n, np.nan)
        ok = np.isfinite(z) & (z > 0)
        d[ok] = Planck18.luminosity_distance(z[ok]).value
        out["dist_mpc"] = d
        if "dist_kind" not in out:
            out["dist_kind"] = np.where(ok, "specz", "none").astype(object)
    if "dist_kind" not in out:
        has = np.isfinite(np.asarray(out.get("dist_mpc", np.full(n, np.nan)), float))
        out["dist_kind"] = np.where(has, "specz", "none").astype(object)
    return out

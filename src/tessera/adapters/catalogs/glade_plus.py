"""GLADE+ (Dalya et al. 2022, MNRAS 514, 1403) -> canonical schema.

Reads the official whitespace-separated ``GLADE+.txt`` (40 columns, missing
values written as ``null``), a VizieR ``asu-tsv`` export of ``VII/291/gladep``,
or the full fixed-width CDS file ``gladep.dat`` (VizieR FTP, 23.2 million rows,
598-character records).  Flag semantics follow the VizieR ReadMe:

* f_dL (dist flag): 0 no distance, 1 photo-z, 2 measured d_L (redshift
  independent), 3 spectroscopic z
* f_zcmb: 1 = redshift/distance corrected for peculiar velocity
* Type: Q = SDSS-DR16Q quasar, G = galaxy
* M*: stellar mass in units of 1e10 Msun
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow as pa
import pyarrow.csv as pacsv

CITATION = "Dalya G. et al. 2022, MNRAS 514, 1403 (GLADE+); VizieR VII/291"

OFFICIAL_COLUMNS = [
    "GLADE+", "PGC", "GWGC", "HyperLEDA", "2MASS", "WISExSCOS", "SDSS-DR16Q", "Type",
    "RAJ2000", "DEJ2000", "Bmag", "e_Bmag", "f_Bmag", "BMAG", "Jmag", "e_Jmag", "Hmag", "e_Hmag",
    "Kmag", "e_Kmag", "W1mag", "e_W1mag", "W2mag", "e_W2mag", "f_W1mag", "Bjmag", "e_Bjmag",
    "zhelio", "zcmb", "f_zcmb", "e_z", "e_zhelio", "dL", "e_dL", "f_dL", "M*", "e_M*", "f_M*",
    "logRate", "e_logRate",
]
_DIST_KIND = {0: "none", 1: "photoz", 2: "z_indep", 3: "specz"}
_NUMERIC = {"RAJ2000", "DEJ2000", "Bmag", "Kmag", "W1mag", "zcmb", "e_zhelio", "dL", "e_dL",
            "M*", "e_M*", "logRate", "f_dL", "f_zcmb", "PGC", "GLADE+"}


def _f(table: dict[str, Any], key: str, n: int) -> np.ndarray:
    if key not in table:
        return np.full(n, np.nan)
    v = table[key]
    if isinstance(v, np.ndarray) and v.dtype.kind in "fiu":
        return v.astype(float)
    out = np.full(n, np.nan)
    for i, s in enumerate(v):
        if s is None:
            continue
        s = str(s).strip()
        if s and s.lower() not in ("null", "nan", "-"):
            try:
                out[i] = float(s)
            except ValueError:
                pass
    return out


def _s(table: dict[str, Any], key: str, n: int) -> np.ndarray:
    if key not in table:
        return np.full(n, "", dtype=object)
    out = np.empty(n, dtype=object)
    for i, s in enumerate(table[key]):
        s = "" if s is None else str(s).strip()
        out[i] = "" if s.lower() in ("null", "-", "nan") else s
    return out


def to_canonical(t: dict[str, Any]) -> dict[str, np.ndarray]:
    n = len(next(iter(t.values())))
    glade = _f(t, "GLADE+", n)
    if np.any(~np.isfinite(glade)):
        raise ValueError("GLADE+ rows without a GLADE+ number")
    fdl = np.nan_to_num(_f(t, "f_dL", n), nan=0).astype(int)
    fz = np.nan_to_num(_f(t, "f_zcmb", n), nan=0).astype(int)
    dist = _f(t, "dL", n)
    kind = np.array([_DIST_KIND.get(int(k), "none") for k in fdl], dtype=object)
    kind[~np.isfinite(dist)] = "none"
    dist[kind == "none"] = np.nan
    mstar = _f(t, "M*", n) * 1e10
    e_mstar = _f(t, "e_M*", n) * 1e10
    with np.errstate(divide="ignore", invalid="ignore"):
        log_m = np.where(mstar > 0, np.log10(mstar), np.nan)
        log_m_err = np.where((mstar > 0) & np.isfinite(e_mstar), e_mstar / (mstar * np.log(10)), np.nan)
    typ = _s(t, "Type", n)
    alt = []
    pgc = _f(t, "PGC", n)
    names = {k: _s(t, k, n) for k in ("GWGC", "HyperLEDA", "2MASS", "WISExSCOS", "SDSS-DR16Q")}
    for i in range(n):
        parts = []
        if np.isfinite(pgc[i]):
            parts.append(f"PGC:{int(pgc[i])}")
        for k, arr in names.items():
            if arr[i]:
                parts.append(f"{k}:{arr[i]}")
        alt.append("|".join(parts))
    return {
        "gal_id": np.array([f"GLADE+:{int(g)}" for g in glade], dtype=object),
        "ra": _f(t, "RAJ2000", n),
        "dec": _f(t, "DEJ2000", n),
        "dist_mpc": dist,
        "dist_err_mpc": _f(t, "e_dL", n),
        "dist_kind": kind,
        # GLADE+ e_dL for peculiar-velocity-corrected spec-z distances already
        # propagates the peculiar-velocity error estimate (column e_z).
        "dist_err_includes_pv": (fz == 1) & (kind == "specz"),
        "z": _f(t, "zcmb", n),
        "z_err": _f(t, "e_zhelio", n),
        "pv_corrected": fz == 1,
        "mag_B": _f(t, "Bmag", n),
        "mag_K": _f(t, "Kmag", n),
        "mag_W1": _f(t, "W1mag", n),
        "log_mstar": log_m,
        "log_mstar_err": log_m_err,
        "log_bns_rate": _f(t, "logRate", n),
        "obj_type": np.where(typ == "Q", "qso", "galaxy").astype(object),
        "alt_ids": np.array(alt, dtype=object),
    }


def read_vizier_tsv(path: str | Path) -> dict[str, np.ndarray]:
    """Parse a VizieR asu-tsv export (header lines start with '#', then names, units, dashes)."""
    lines = [ln for ln in Path(path).read_text().splitlines() if ln and not ln.startswith("#")]
    header = lines[0].split("\t")
    body = [ln for ln in lines[3:] if ln.strip()]
    cols: dict[str, list[str]] = {h.strip(): [] for h in header}
    keys = [h.strip() for h in header]
    for ln in body:
        parts = ln.split("\t")
        parts += [""] * (len(keys) - len(parts))
        for k, v in zip(keys, parts):
            cols[k].append(v)
    return to_canonical({k: np.array(v, dtype=object) for k, v in cols.items()})


def iter_official(path: str | Path, block_size: int = 1 << 26) -> Iterator[dict[str, np.ndarray]]:
    """Stream the official GLADE+.txt in chunks (pyarrow CSV reader)."""
    read_opts = pacsv.ReadOptions(column_names=OFFICIAL_COLUMNS, block_size=block_size)
    parse_opts = pacsv.ParseOptions(delimiter=" ", ignore_empty_lines=True)
    types = {c: (pa.float64() if c in _NUMERIC else pa.string()) for c in OFFICIAL_COLUMNS}
    conv = pacsv.ConvertOptions(null_values=["null"], strings_can_be_null=True, column_types=types)
    with pacsv.open_csv(str(path), read_options=read_opts, parse_options=parse_opts, convert_options=conv) as rd:
        for batch in rd:
            t: dict[str, Any] = {}
            for name in OFFICIAL_COLUMNS:
                arr = batch.column(name)
                if pa.types.is_floating(arr.type):
                    t[name] = arr.to_numpy(zero_copy_only=False).astype(float)
                else:
                    t[name] = np.asarray(arr.to_pylist(), dtype=object)
            yield to_canonical(t)


def read_official(path: str | Path) -> dict[str, np.ndarray]:
    parts = list(iter_official(path))
    if not parts:
        raise ValueError("empty GLADE+ file")
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


# 1-based inclusive byte ranges from the VizieR ReadMe of VII/291 (file gladep.dat)
CDS_RECORD = 598
CDS_COLUMNS = {
    "GLADE+": (1, 8), "PGC": (10, 16), "GWGC": (18, 45), "HyperLEDA": (47, 75), "2MASS": (77, 92),
    "WISExSCOS": (94, 112), "SDSS-DR16Q": (114, 131), "Type": (133, 133), "RAJ2000": (135, 155),
    "DEJ2000": (157, 179), "Bmag": (181, 198), "Kmag": (310, 327), "W1mag": (348, 365),
    "zcmb": (449, 471), "f_zcmb": (473, 473), "e_zhelio": (498, 519), "dL": (521, 543),
    "e_dL": (545, 564), "f_dL": (566, 566), "M*": (568, 577), "e_M*": (579, 587), "logRate": (591, 593),
}


def _is_cds(head: bytes) -> bool:
    lines = head.split(b"\n")
    return len(lines) > 1 and len(lines[0]) == CDS_RECORD and lines[0][:8].strip().isdigit()


def iter_cds(path: str | Path, rows_per_chunk: int = 500_000) -> Iterator[dict[str, np.ndarray]]:
    """Stream the fixed-width CDS ``gladep.dat`` in chunks, slicing the ReadMe byte columns."""
    import pandas as pd

    rec = CDS_RECORD + 1                       # + newline
    with Path(path).open("rb") as fh:
        while True:
            buf = fh.read(rows_per_chunk * rec)
            if not buf:
                break
            n = len(buf) // rec
            if n == 0:
                raise ValueError("truncated final record in CDS file")
            if len(buf) % rec:
                raise ValueError(f"CDS file is not a whole number of {rec}-byte records (incomplete download?)")
            arr = np.frombuffer(buf, dtype="S1").reshape(n, rec)
            t: dict[str, Any] = {}
            for name, (a, b) in CDS_COLUMNS.items():
                raw = np.ascontiguousarray(arr[:, a - 1:b]).view(f"S{b - a + 1}").ravel()
                s = pd.Series(raw).str.decode("ascii").str.strip()
                if name in _NUMERIC:
                    t[name] = pd.to_numeric(s.where(s != "-"), errors="coerce").to_numpy(float)
                else:
                    t[name] = s.where(~s.isin(["", "-"]), None).to_numpy(object)
            yield to_canonical(t)


def read_cds(path: str | Path) -> dict[str, np.ndarray]:
    parts = list(iter_cds(path))
    if not parts:
        raise ValueError("empty GLADE+ file")
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def read(path: str | Path) -> dict[str, np.ndarray]:
    """Auto-detect official text, VizieR TSV, or the fixed-width CDS gladep.dat."""
    raw = Path(path).open("rb").read(4096)
    if _is_cds(raw):
        return read_cds(path)
    head = raw.decode(errors="replace")
    if head.startswith("#") or "\t" in head.splitlines()[0]:
        return read_vizier_tsv(path)
    return read_official(path)


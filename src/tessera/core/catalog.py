"""Versioned galaxy-catalog snapshots (spec section 4).

A snapshot is an immutable directory::

    snapshot.json          metadata, provenance, sha256 of every file
    galaxies.parquet       canonical schema, sorted by hpx10 (nested, order 10)
    completeness.npz       completeness cube per weight kind (pix6 x distance shell)
    footprint.npz          order-6 footprint mask
    ebv.npz                optional E(B-V) map (nested HEALPix, any order)

The snapshot ID is the sha256 of the canonical ``snapshot.json``.  Nothing in
the planning path touches the network.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from astropy_healpix import lonlat_to_healpix
from scipy.spatial import cKDTree

from tessera.core import science
from tessera.core.errors import CatalogError
from tessera.core.geometry import chord_for_angle, radec_to_xyz
from tessera.core.hashing import canonical_json, sha256_bytes, sha256_file
from tessera.core.weights import tracer

SNAPSHOT_SCHEMA = "tessera.catalog-snapshot/1"
INDEX_ORDER = 10
CUBE_ORDER = 4  # ~13.4 deg2 cells: coarse enough to suppress Poisson noise in C_w
DIST_KINDS = ("none", "photoz", "z_indep", "specz")

#: Canonical column -> numpy dtype ("U" = string)
CANONICAL: dict[str, str] = {
    "gal_id": "U", "ra": "f8", "dec": "f8",
    "dist_mpc": "f8", "dist_err_mpc": "f8", "dist_kind": "U", "dist_err_includes_pv": "?",
    "z": "f8", "z_err": "f8", "pv_corrected": "?",
    "mag_B": "f4", "mag_K": "f4", "mag_W1": "f4",
    "log_mstar": "f4", "log_mstar_err": "f4", "log_bns_rate": "f4",
    "obj_type": "U", "alt_ids": "U", "duplicate_of": "U", "hpx10": "i8",
}
OPTIONAL = {"sfr": "f4"}


def default_shell_edges() -> np.ndarray:
    """10 Mpc shells to 200 Mpc, then logarithmic to 5000 Mpc."""
    return np.unique(np.concatenate([np.arange(0, 200, 10.0), np.geomspace(200, 5000, 34)]))


@dataclass
class Galaxies:
    """Column store in the canonical schema (numpy arrays of equal length)."""

    cols: dict[str, np.ndarray]

    def __len__(self) -> int:
        return len(self.cols["gal_id"])

    def __getitem__(self, key: str) -> np.ndarray:
        return self.cols[key]

    def take(self, idx: np.ndarray) -> Galaxies:
        return Galaxies({k: v[idx] for k, v in self.cols.items()})

    @classmethod
    def empty(cls) -> Galaxies:
        return cls({k: np.array([], dtype=object if t == "U" else t) for k, t in CANONICAL.items()})


def normalise_columns(cols: dict[str, Any]) -> dict[str, np.ndarray]:
    """Coerce adapter output to the canonical schema, filling defaults."""
    n = len(cols["ra"])
    out: dict[str, np.ndarray] = {}
    for name, dt in {**CANONICAL, **OPTIONAL}.items():
        if name == "hpx10":
            continue
        if name not in cols:
            if name in OPTIONAL:
                continue
            if dt == "U":
                default: Any = {"dist_kind": "none", "obj_type": "galaxy"}.get(name, "")
                out[name] = np.full(n, default, dtype=object)
            elif dt == "?":
                out[name] = np.zeros(n, dtype=bool)
            else:
                out[name] = np.full(n, np.nan, dtype=dt)
            continue
        v = np.asarray(cols[name])
        out[name] = v.astype(object) if dt == "U" else v.astype(dt)
    if np.any(~np.isfinite(out["ra"])) or np.any(~np.isfinite(out["dec"])):
        raise CatalogError("catalog rows with non-finite coordinates", code="CATALOG_INVALID")
    if np.any((out["dec"] < -90) | (out["dec"] > 90)):
        raise CatalogError("declination out of range", code="CATALOG_INVALID")
    out["ra"] = np.mod(out["ra"], 360.0)
    bad = ~np.isin(out["dist_kind"].astype(str), DIST_KINDS)
    if bad.any():
        raise CatalogError("invalid dist_kind values", code="CATALOG_INVALID",
                           details={"values": sorted(set(out["dist_kind"][bad].astype(str)))[:10]})
    if len(set(out["gal_id"])) != n:
        raise CatalogError("gal_id values are not unique", code="CATALOG_INVALID")
    return out


# ---------------------------------------------------------------------------
# Duplicates (spec 4.2)
# ---------------------------------------------------------------------------

_KIND_RANK = {"z_indep": 0, "specz": 1, "photoz": 2, "none": 3}


def find_duplicates(cols: dict[str, np.ndarray], radius_arcsec: float = 3.0, nsigma: float = 3.0) -> np.ndarray:
    """Return ``duplicate_of`` (empty string for representatives)."""
    n = len(cols["ra"])
    dup = np.full(n, "", dtype=object)
    if n < 2:
        return dup
    xyz = radec_to_xyz(cols["ra"], cols["dec"])
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(chord_for_angle(radius_arcsec / 3600.0), output_type="ndarray")
    if len(pairs) == 0:
        return dup
    d = cols["dist_mpc"]
    s = np.where(np.isfinite(cols["dist_err_mpc"]), cols["dist_err_mpc"], 0.1 * d)
    i, j = pairs[:, 0], pairs[:, 1]
    both = np.isfinite(d[i]) & np.isfinite(d[j])
    agree = ~both | (np.abs(d[i] - d[j]) <= nsigma * np.sqrt(s[i] ** 2 + s[j] ** 2 + 1e-12))
    pairs = pairs[agree]
    if len(pairs) == 0:
        return dup
    # connected components
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, labels = connected_components(g, directed=False)
    members = np.unique(pairs.ravel())
    kind_rank = np.array([_KIND_RANK.get(str(k), 3) for k in cols["dist_kind"]])
    derr = np.where(np.isfinite(cols["dist_err_mpc"]), cols["dist_err_mpc"], np.inf)
    kmag = np.where(np.isfinite(cols["mag_K"]), cols["mag_K"], np.inf)
    ids = cols["gal_id"].astype(str)
    for lab in np.unique(labels[members]):
        grp = members[labels[members] == lab]
        best = grp[np.lexsort((ids[grp], kmag[grp], derr[grp], kind_rank[grp]))[0]]
        for k in grp:
            if k != best:
                dup[k] = ids[best]
    return dup


# ---------------------------------------------------------------------------
# Building snapshots
# ---------------------------------------------------------------------------

def _hpx(ra: np.ndarray, dec: np.ndarray, order: int) -> np.ndarray:
    return np.asarray(lonlat_to_healpix(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg, 2**order,
                                        order="nested"), dtype=np.int64)


def build_completeness_cube(cols: dict[str, np.ndarray], edges: np.ndarray, weights: list[str],
                            ebv: np.ndarray | None = None) -> dict[str, np.ndarray]:
    """C_w(pix6, shell) = cataloged weight / (rho_bar_w * shell volume per pixel), capped at 1."""
    npix6 = 12 * 4**CUBE_ORDER
    pix6 = _hpx(cols["ra"], cols["dec"], CUBE_ORDER)
    d = cols["dist_mpc"]
    ok = np.isfinite(d) & (d > 0) & (cols["duplicate_of"].astype(str) == "")
    if "obj_type" in cols:
        ok &= cols["obj_type"].astype(str) != "qso"
    shell = np.searchsorted(edges, d, side="right") - 1
    ok &= (shell >= 0) & (shell < len(edges) - 1)
    vol = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3) / npix6
    out: dict[str, np.ndarray] = {}
    for kind in weights:
        ref = science.reference_density(kind)
        if ref is None:
            continue
        raw, _ = tracer(cols, kind, np.where(ok, d, np.nan), ebv)
        good = ok & np.isfinite(raw) & (raw > 0)
        flat_idx = pix6[good] * (len(edges) - 1) + shell[good]
        tot = np.bincount(flat_idx, weights=raw[good], minlength=npix6 * (len(edges) - 1))
        cube = tot.reshape(npix6, len(edges) - 1) / (float(ref["value"]) * vol[None, :])
        out[kind] = np.minimum(_smooth_neighbours(cube), 1.0).astype(np.float32)
    return out


def _smooth_neighbours(cube: np.ndarray) -> np.ndarray:
    """Mean of each cell and its 8 HEALPix neighbours in the same distance shell.

    Smoothing *before* capping at 1 stops single bright galaxies in sparse cells from
    being clipped away, which otherwise biases P_cat low (measured in the benchmark's
    halved-catalog test)."""
    from astropy_healpix import neighbours

    npix = cube.shape[0]
    with np.errstate(invalid="ignore"):
        nb = np.asarray(neighbours(np.arange(npix), 2**CUBE_ORDER, order="nested"))   # (8, npix), -1 = none
    acc = cube.copy()
    cnt = np.ones(npix)
    for row in nb:
        ok = row >= 0
        acc[ok] += cube[row[ok]]
        cnt[ok] += 1
    return acc / cnt[:, None]


def build_snapshot(
    cols: dict[str, Any],
    out_dir: str | Path,
    *,
    name: str,
    upstream_version: str,
    source_url: str | None = None,
    retrieved: str | None = None,
    citation: str | None = None,
    provenance: dict[str, str] | None = None,
    footprint: str | np.ndarray = "allsky",
    ebv_map: np.ndarray | None = None,
    dedupe_radius_arcsec: float = 3.0,
    row_group_size: int = 65536,
    cube_weights: tuple[str, ...] = ("lum_B", "lum_K", "mstar"),
) -> Snapshot:
    """Write an immutable snapshot directory and return the opened snapshot."""
    from tessera import __version__

    out = Path(out_dir)
    if (out / "snapshot.json").exists():
        raise CatalogError(f"snapshot already exists at {out}; snapshots are immutable", code="CATALOG_EXISTS")
    out.mkdir(parents=True, exist_ok=True)
    c = normalise_columns(cols)
    if not np.any(c["duplicate_of"].astype(str) != ""):
        c["duplicate_of"] = find_duplicates(c, dedupe_radius_arcsec)
    c["hpx10"] = _hpx(c["ra"], c["dec"], INDEX_ORDER)
    order = np.lexsort((c["gal_id"].astype(str), c["hpx10"]))
    c = {k: v[order] for k, v in c.items()}

    ebv_gal = None
    if ebv_map is not None:
        ebv_map = np.asarray(ebv_map, dtype=np.float32)
        o = int(round(np.log2(np.sqrt(len(ebv_map) / 12))))
        ebv_gal = ebv_map[_hpx(c["ra"], c["dec"], o)]

    arrays = {}
    for k, v in c.items():
        arrays[k] = pa.array(v.astype(str) if v.dtype == object else v)
    table = pa.table(arrays)
    pq.write_table(table, out / "galaxies.parquet", row_group_size=row_group_size, compression="zstd",
                   write_statistics=True)

    edges = default_shell_edges()
    cube = build_completeness_cube(c, edges, list(cube_weights), ebv_gal)
    np.savez_compressed(out / "completeness.npz", edges=edges, **cube)

    npix6 = 12 * 4**CUBE_ORDER
    if isinstance(footprint, str):
        if footprint == "allsky":
            fp = np.ones(npix6, dtype=bool)
        elif footprint == "occupied":
            fp = np.zeros(npix6, dtype=bool)
            fp[_hpx(c["ra"], c["dec"], CUBE_ORDER)] = True
        else:
            raise CatalogError(f"unknown footprint mode {footprint!r}")
    else:
        fp = np.asarray(footprint, dtype=bool)
        if fp.shape != (npix6,):
            raise CatalogError("footprint mask must have 12*4**6 entries (order 6, nested)")
    np.savez_compressed(out / "footprint.npz", mask=fp)
    files = ["galaxies.parquet", "completeness.npz", "footprint.npz"]
    if ebv_map is not None:
        np.savez_compressed(out / "ebv.npz", ebv=ebv_map)
        files.append("ebv.npz")

    n_dup = int(np.sum(c["duplicate_of"].astype(str) != ""))
    kinds, counts = np.unique(c["dist_kind"].astype(str), return_counts=True)
    meta = {
        "schema": SNAPSHOT_SCHEMA,
        "name": name,
        "upstream_version": upstream_version,
        "source_url": source_url,
        "retrieved": retrieved,
        "citation": citation,
        "ingested_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "ingest_code_version": __version__,
        "science_model": science.defaults()["science_model"],
        "cosmology": science.defaults()["cosmology"],
        "row_count": int(len(c["gal_id"])),
        "duplicates_flagged": n_dup,
        "dist_kind_counts": {str(k): int(v) for k, v in zip(kinds, counts)},
        "index": {"column": "hpx10", "order": INDEX_ORDER, "scheme": "nested"},
        "cube": {"order": CUBE_ORDER, "weights": sorted(cube), "shells_mpc": "see completeness.npz:edges"},
        "footprint": footprint if isinstance(footprint, str) else "mask",
        "has_ebv": ebv_map is not None,
        "provenance": provenance or {},
        "files": {f: sha256_file(out / f) for f in files},
    }
    (out / "snapshot.json").write_text(json.dumps(meta, indent=2, sort_keys=True))
    return open_snapshot(out, verify=False)


# ---------------------------------------------------------------------------
# Reading snapshots
# ---------------------------------------------------------------------------

class Snapshot:
    """Read-only access to a catalog snapshot."""

    def __init__(self, path: Path, meta: dict[str, Any], snapshot_id: str):
        self.path = path
        self.meta = meta
        self.id = snapshot_id
        self._pf = pq.ParquetFile(path / "galaxies.parquet")
        md = self._pf.metadata
        col_idx = self._pf.schema_arrow.get_field_index("hpx10")
        self._rg_min = np.empty(md.num_row_groups, dtype=np.int64)
        self._rg_max = np.empty(md.num_row_groups, dtype=np.int64)
        for i in range(md.num_row_groups):
            st = md.row_group(i).column(col_idx).statistics
            self._rg_min[i] = st.min
            self._rg_max[i] = st.max
        self._cube: dict[str, np.ndarray] | None = None
        self._edges: np.ndarray | None = None
        self._fp: np.ndarray | None = None
        self._ebv: np.ndarray | None = None

    @property
    def name(self) -> str:
        return self.meta["name"]

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "upstream_version": self.meta["upstream_version"], "snapshot_id": self.id,
                "row_count": self.meta["row_count"], "cube_weights": self.meta["cube"]["weights"],
                "footprint": self.meta["footprint"], "has_ebv": self.meta["has_ebv"],
                "citation": self.meta.get("citation")}

    # -- region extraction ---------------------------------------------------
    def extract_ranges(self, ranges: np.ndarray) -> Galaxies:
        """Galaxies whose hpx10 lies in any half-open range [lo, hi) (sorted, merged)."""
        if len(ranges) == 0:
            return Galaxies.empty()
        lo, hi = ranges[:, 0], ranges[:, 1]
        # row groups overlapping any range: rg_min < hi_k and rg_max >= lo_k for some k
        k = np.searchsorted(hi, self._rg_min, side="right")
        k = np.clip(k, 0, len(ranges) - 1)
        sel = np.flatnonzero((lo[k] <= self._rg_max) & (self._rg_min < hi[k]))
        # The vectorised test above checks only the first range ending after rg_min;
        # it overlaps iff that range starts at or before rg_max.
        if len(sel) == 0:
            return Galaxies.empty()
        tbl = self._pf.read_row_groups(sel.tolist())
        hp = tbl.column("hpx10").to_numpy()
        pos = np.searchsorted(lo, hp, side="right") - 1
        keep = (pos >= 0) & (hp < hi[np.clip(pos, 0, None)])
        idx = np.flatnonzero(keep)
        tbl = tbl.take(pa.array(idx))
        cols: dict[str, np.ndarray] = {}
        for name in tbl.column_names:
            arr = tbl.column(name)
            if pa.types.is_string(arr.type) or pa.types.is_large_string(arr.type):
                cols[name] = np.asarray(arr.to_pylist(), dtype=object)
            else:
                cols[name] = arr.to_numpy(zero_copy_only=False)
        return Galaxies(cols)

    def extract_for_map(self, order: np.ndarray, ipix: np.ndarray) -> Galaxies:
        """Galaxies inside the given multi-order pixels."""
        return self.extract_ranges(pixels_to_index_ranges(order, ipix))

    # -- completeness, footprint, extinction ---------------------------------
    def _load_cube(self) -> None:
        if self._cube is None:
            with np.load(self.path / "completeness.npz") as z:
                self._edges = z["edges"]
                self._cube = {k: z[k] for k in z.files if k != "edges"}

    def completeness(self, weight: str) -> tuple[np.ndarray, np.ndarray] | None:
        self._load_cube()
        assert self._cube is not None and self._edges is not None
        if weight not in self._cube:
            return None
        return self._edges, self._cube[weight]

    def footprint_mask(self) -> np.ndarray:
        if self._fp is None:
            with np.load(self.path / "footprint.npz") as z:
                self._fp = z["mask"]
        return self._fp

    def ebv(self, ra: np.ndarray, dec: np.ndarray) -> np.ndarray | None:
        if not self.meta.get("has_ebv"):
            return None
        if self._ebv is None:
            with np.load(self.path / "ebv.npz") as z:
                self._ebv = z["ebv"]
        o = int(round(np.log2(np.sqrt(len(self._ebv) / 12))))
        return self._ebv[_hpx(ra, dec, o)].astype(float)

    def verify(self) -> None:
        for f, digest in self.meta["files"].items():
            if sha256_file(self.path / f) != digest:
                raise CatalogError(f"checksum mismatch for {f}", code="CATALOG_CHECKSUM_MISMATCH")


def pixels_to_index_ranges(order: np.ndarray, ipix: np.ndarray, index_order: int = INDEX_ORDER) -> np.ndarray:
    """Merge multi-order pixels into sorted half-open ranges of nested indices at ``index_order``."""
    order = np.asarray(order, dtype=np.int64)
    ipix = np.asarray(ipix, dtype=np.int64)
    fine = order >= index_order
    lo = np.empty(len(order), dtype=np.int64)
    hi = np.empty(len(order), dtype=np.int64)
    sh = 2 * (order[fine] - index_order)
    lo[fine] = ipix[fine] >> sh
    hi[fine] = lo[fine] + 1
    sh = 2 * (index_order - order[~fine])
    lo[~fine] = ipix[~fine] << sh
    hi[~fine] = (ipix[~fine] + 1) << sh
    srt = np.argsort(lo, kind="stable")
    lo, hi = lo[srt], hi[srt]
    if len(lo) == 0:
        return np.empty((0, 2), dtype=np.int64)
    # merge
    run_max = np.maximum.accumulate(hi)
    starts = np.concatenate([[True], lo[1:] > run_max[:-1]])
    group = np.cumsum(starts) - 1
    out_lo = lo[starts]
    out_hi = np.zeros(len(out_lo), dtype=np.int64)
    np.maximum.at(out_hi, group, hi)
    return np.stack([out_lo, out_hi], axis=1)


def open_snapshot(path: str | Path, *, verify: bool = False) -> Snapshot:
    p = Path(path)
    meta_path = p / "snapshot.json"
    if not meta_path.is_file():
        raise CatalogError(f"no catalog snapshot at {p}", code="CATALOG_SNAPSHOT_MISSING")
    meta = json.loads(meta_path.read_text())
    if meta.get("schema") != SNAPSHOT_SCHEMA:
        raise CatalogError(f"unsupported snapshot schema {meta.get('schema')!r}", code="CATALOG_SCHEMA")
    snap = Snapshot(p, meta, sha256_bytes(canonical_json(meta).encode()))
    if verify:
        snap.verify()
    return snap

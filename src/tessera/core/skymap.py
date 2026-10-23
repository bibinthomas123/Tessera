"""Sky-map ingestion, validation and the internal multi-order representation.

Every accepted input (flat RING/NESTED, optionally gzip-compressed, or
multi-order NUNIQ) is converted to :class:`MultiOrderMap`:

* ``uniq``         HEALPix NUNIQ index of each pixel
* ``probdensity``  probability per steradian (rho)
* ``distmu``, ``distsigma``, ``distnorm``  distance ansatz (optional)

Validation is strict: malformed or unexpectedly large inputs are rejected
with a typed :class:`~tessera.core.errors.SkymapError` before the data section is
parsed.  See science spec section 3.1.
"""

from __future__ import annotations

import io
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy_healpix import healpix_to_lonlat, lonlat_to_healpix
from astropy_healpix.core import ring_to_nested

from tessera.core.distance import valid_distance
from tessera.core.errors import Diagnostics, SkymapError
from tessera.core.hashing import sha256_bytes

FOUR_PI = 4.0 * np.pi
_HEADER_KEYS = ("OBJECT", "DATE-OBS", "MJD-OBS", "CREATOR", "ORIGIN", "DISTMEAN", "DISTSTD", "LOGBCI",
                "LOGBSN", "INSTRUME", "RUNTIME", "MOCORDER", "NSIDE", "ORDERING", "INDXSCHM", "VCSVERS")


@dataclass(frozen=True)
class SkymapLimits:
    """Resource limits applied before and during parsing."""

    max_bytes: int = 100 * 1024**2           # raw upload size
    max_decompressed_bytes: int = 1024**3    # 1 GiB, enforced while streaming gzip
    max_nside: int = 2048
    max_order: int = 29
    normalisation_accept: float = 1e-3
    normalisation_reject: float = 1e-2


@dataclass
class SkymapMeta:
    filename: str
    sha256_raw: str
    sha256_data: str
    format: str                    # "flat" | "multiorder"
    ordering: str                  # "RING" | "NESTED" | "NUNIQ"
    nside: int | None
    max_order: int
    n_pixels: int
    header: dict[str, Any] = field(default_factory=dict)
    distance_prior: str = "unknown"

    @property
    def object(self) -> str | None:
        return self.header.get("OBJECT")

    @property
    def mjd_obs(self) -> float | None:
        v = self.header.get("MJD-OBS")
        return float(v) if v is not None else None

    def to_dict(self) -> dict[str, Any]:
        return {
            "filename": self.filename,
            "sha256_raw": self.sha256_raw,
            "sha256_data": self.sha256_data,
            "format": self.format,
            "ordering": self.ordering,
            "nside": self.nside,
            "max_order": self.max_order,
            "n_pixels": self.n_pixels,
            "header": self.header,
            "distance_prior": self.distance_prior,
        }


class MultiOrderMap:
    """Immutable multi-order HEALPix sky map with optional distance layers."""

    def __init__(
        self,
        uniq: np.ndarray,
        probdensity: np.ndarray,
        distmu: np.ndarray | None = None,
        distsigma: np.ndarray | None = None,
        distnorm: np.ndarray | None = None,
        meta: SkymapMeta | None = None,
        diagnostics: Diagnostics | None = None,
    ):
        self.uniq = np.ascontiguousarray(uniq, dtype=np.int64)
        self.probdensity = np.ascontiguousarray(probdensity, dtype=np.float64)
        self.distmu = None if distmu is None else np.ascontiguousarray(distmu, dtype=np.float64)
        self.distsigma = None if distsigma is None else np.ascontiguousarray(distsigma, dtype=np.float64)
        self.distnorm = None if distnorm is None else np.ascontiguousarray(distnorm, dtype=np.float64)
        self.order, self.ipix = _uniq_to_order_ipix(self.uniq)
        self.area = FOUR_PI / (12.0 * 4.0 ** self.order.astype(np.float64))
        self.meta = meta
        self.diagnostics = diagnostics if diagnostics is not None else Diagnostics()
        self._lookup: dict[int, tuple[np.ndarray, np.ndarray]] | None = None
        self._credible: np.ndarray | None = None
        for arr in (self.uniq, self.probdensity, self.order, self.ipix, self.area):
            arr.setflags(write=False)

    # -- basic properties ---------------------------------------------------
    def __len__(self) -> int:
        return len(self.uniq)

    @property
    def prob(self) -> np.ndarray:
        """Probability contained in each pixel."""
        return self.probdensity * self.area

    @property
    def has_distance_columns(self) -> bool:
        return self.distmu is not None

    @property
    def distance_valid(self) -> np.ndarray:
        if self.distmu is None:
            return np.zeros(len(self), dtype=bool)
        return valid_distance(self.distmu, self.distsigma, self.distnorm)

    @property
    def max_order(self) -> int:
        return int(self.order.max())

    # -- lookups ------------------------------------------------------------
    def _build_lookup(self) -> dict[int, tuple[np.ndarray, np.ndarray]]:
        if self._lookup is None:
            table: dict[int, tuple[np.ndarray, np.ndarray]] = {}
            for o in np.unique(self.order):
                idx = np.flatnonzero(self.order == o)
                ip = self.ipix[idx]
                srt = np.argsort(ip, kind="stable")
                table[int(o)] = (ip[srt], idx[srt])
            self._lookup = table
        return self._lookup

    def lookup(self, ra_deg, dec_deg) -> np.ndarray:
        """Index of the pixel containing each position, or -1 (partial maps).

        No interpolation is performed (spec 3.2): a point belongs to exactly
        one multi-order pixel.
        """
        ra = np.atleast_1d(np.asarray(ra_deg, dtype=float))
        dec = np.atleast_1d(np.asarray(dec_deg, dtype=float))
        out = np.full(ra.shape, -1, dtype=np.int64)
        lon, lat = ra * u.deg, dec * u.deg
        for o, (ip_sorted, idx_sorted) in self._build_lookup().items():
            todo = out < 0
            if not todo.any():
                break
            pix = lonlat_to_healpix(lon[todo], lat[todo], 2**o, order="nested")
            pos = np.searchsorted(ip_sorted, pix)
            pos_c = np.clip(pos, 0, len(ip_sorted) - 1)
            hit = ip_sorted[pos_c] == pix
            sub = np.flatnonzero(todo)
            out[sub[hit]] = idx_sorted[pos_c[hit]]
        return out

    def pixel_centers(self, index: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        """(ra, dec) in degrees of pixel centres."""
        idx = np.arange(len(self)) if index is None else np.asarray(index)
        ra = np.empty(len(idx))
        dec = np.empty(len(idx))
        orders = self.order[idx]
        for o in np.unique(orders):
            m = orders == o
            lon, lat = healpix_to_lonlat(self.ipix[idx[m]], 2**int(o), order="nested")
            ra[m] = lon.to_value(u.deg)
            dec[m] = lat.to_value(u.deg)
        return ra, dec

    # -- credible regions ---------------------------------------------------
    def credible_levels(self) -> np.ndarray:
        """Cumulative probability of the smallest region containing each pixel.

        Pixels are sorted by density (descending, ties broken by UNIQ) and the
        cumulative probability *including* the pixel is returned, matching the
        ligo.skymap ``find_greedy_credible_levels`` convention.
        """
        if self._credible is None:
            order = np.lexsort((self.uniq, -self.probdensity))
            cum = np.cumsum(self.prob[order])
            cl = np.empty(len(self))
            cl[order] = cum
            cl.setflags(write=False)
            self._credible = cl
        return self._credible

    def region_mask(self, level: float) -> np.ndarray:
        return self.credible_levels() <= level

    def region_mask_inclusive(self, level: float) -> np.ndarray:
        """Smallest set of pixels whose cumulative probability reaches ``level``."""
        cl = self.credible_levels()
        mask = cl <= level
        above = cl[~mask]
        if above.size:
            mask |= cl == above.min()
        return mask

    def area_deg2(self, level: float) -> float:
        return float(self.area[self.region_mask(level)].sum() * (180 / np.pi) ** 2)

    # -- event classification (spec 3.1) -----------------------------------
    def distance_class(self, level: float = 0.9) -> str:
        """``"none"``, ``"weak"`` or ``"ok"``."""
        if self.distmu is None:
            return "none"
        mask = self.region_mask_inclusive(level)
        ok = self.distance_valid & mask
        if not ok.any():
            return "none"
        w = self.prob[ok]
        frac = self.distsigma[ok] / np.abs(self.distmu[ok])
        srt = np.argsort(frac)
        cw = np.cumsum(w[srt]) / w.sum()
        median = frac[srt][np.searchsorted(cw, 0.5)]
        return "weak" if median > 0.5 else "ok"

    def summary(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "n_pixels": len(self),
            "max_order": self.max_order,
            "area_50_deg2": self.area_deg2(0.5),
            "area_90_deg2": self.area_deg2(0.9),
            "distance_class": self.distance_class(),
        }
        if self.meta is not None:
            out["distmean"] = self.meta.header.get("DISTMEAN")
            out["diststd"] = self.meta.header.get("DISTSTD")
        return out

    def preview(self, level: float = 0.99, order: int = 7) -> list[list[float]]:
        """Coarse (ra, dec, prob) samples of the credible region, for lightweight plotting."""
        mask = self.region_mask_inclusive(level)
        o, ip, p = self.order[mask], self.ipix[mask], self.prob[mask]
        coarse = np.where(o > order, ip >> (2 * np.maximum(o - order, 0)), ip << (2 * np.maximum(order - o, 0)))
        # pixels coarser than ``order`` are represented by their first child (centre approx.)
        uniq_c, inv = np.unique(coarse, return_inverse=True)
        prob = np.bincount(inv, weights=p)
        lon, lat = healpix_to_lonlat(uniq_c, 2**order, order="nested")
        return [[round(float(a), 4), round(float(b), 4), float(c)]
                for a, b, c in zip(lon.to_value(u.deg), lat.to_value(u.deg), prob)]

    def region_moc(self, level: float, max_order: int | None = None) -> dict[str, list[int]]:
        """Credible region as a MOC JSON dict ``{"order": [ipix, ...]}`` (Aladin format)."""
        mask = self.region_mask_inclusive(level)
        orders = self.order[mask]
        ipix = self.ipix[mask]
        if max_order is not None:
            coarse = orders > max_order
            ipix = np.where(coarse, ipix >> (2 * (orders - max_order)), ipix)
            orders = np.minimum(orders, max_order)
        out: dict[str, list[int]] = {}
        for o in np.unique(orders):
            out[str(int(o))] = sorted(set(int(i) for i in ipix[orders == o]))
        return out


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _uniq_to_order_ipix(uniq: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if np.any(uniq < 4):
        raise SkymapError("UNIQ indices must be >= 4", code="SKYMAP_INVALID_FORMAT")
    order = (np.floor(np.log2(uniq.astype(np.float64)) / 2.0) - 1).astype(np.int64)
    # Correct possible floating rounding at exact powers of four.
    low = (4 * (1 << (2 * order))) > uniq
    order[low] -= 1
    high = (4 * (1 << (2 * (order + 1)))) <= uniq
    order[high] += 1
    ipix = uniq - 4 * (np.int64(1) << (2 * order))
    return order, ipix


def _read_source(source: str | Path | bytes, limits: SkymapLimits) -> tuple[bytes, str]:
    if isinstance(source, (bytes, bytearray, memoryview)):
        raw = bytes(source)
        name = "<bytes>"
    else:
        path = Path(source)
        if not path.is_file():
            raise SkymapError(f"sky map not found: {path}", code="SKYMAP_NOT_FOUND")
        if path.stat().st_size > limits.max_bytes:
            raise SkymapError(f"sky map exceeds {limits.max_bytes} bytes", code="SKYMAP_TOO_LARGE")
        raw = path.read_bytes()
        name = path.name
    if len(raw) > limits.max_bytes:
        raise SkymapError(f"sky map exceeds {limits.max_bytes} bytes", code="SKYMAP_TOO_LARGE")
    return raw, name


def _maybe_decompress(raw: bytes, limits: SkymapLimits) -> bytes:
    if raw[:2] != b"\x1f\x8b":
        return raw
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    out = io.BytesIO()
    total = 0
    data = raw
    try:
        while data:
            chunk = d.decompress(data, 1 << 22)
            total += len(chunk)
            if total > limits.max_decompressed_bytes:
                raise SkymapError("decompressed sky map exceeds limit", code="SKYMAP_TOO_LARGE",
                                  details={"limit_bytes": limits.max_decompressed_bytes})
            out.write(chunk)
            data = d.unconsumed_tail
            if d.eof:
                break
        if not d.eof:
            raise SkymapError("truncated gzip stream", code="SKYMAP_INVALID_FORMAT")
    except zlib.error as exc:
        raise SkymapError(f"corrupt gzip stream: {exc}", code="SKYMAP_INVALID_FORMAT") from exc
    return out.getvalue()


def _column(table, names: tuple[str, ...]) -> np.ndarray | None:
    cols = {c.upper(): c for c in table.columns.names}
    for n in names:
        if n in cols:
            return np.asarray(table.data[cols[n]], dtype=np.float64).ravel()
    return None


def _infer_distance_prior(header: dict[str, Any]) -> str:
    creator = str(header.get("CREATOR", "")).upper()
    return "r2" if "BAYESTAR" in creator else "unknown"


def load_skymap(
    source: str | Path | bytes,
    *,
    filename: str | None = None,
    limits: SkymapLimits | None = None,
) -> MultiOrderMap:
    """Load and validate a LIGO/Virgo/KAGRA HEALPix sky map.

    Accepts flat (RING or NESTED, implicit or explicit indexing) and
    multi-order (NUNIQ) maps, optionally gzip-compressed.
    """
    limits = limits or SkymapLimits()
    raw, name = _read_source(source, limits)
    data = _maybe_decompress(raw, limits)
    sha_raw = sha256_bytes(raw)
    sha_data = sha_raw if data is raw else sha256_bytes(data)
    del raw
    diags = Diagnostics()
    if not data.startswith(b"SIMPLE  ="):
        raise SkymapError("not a FITS file", code="SKYMAP_INVALID_FORMAT")
    try:
        hdul = fits.open(io.BytesIO(data), memmap=False, lazy_load_hdus=True)
        if len(hdul) < 2:
            raise SkymapError("FITS file has no binary-table extension", code="SKYMAP_INVALID_FORMAT")
        hdu = hdul[1]
    except SkymapError:
        raise
    except Exception as exc:  # astropy raises many types for malformed files
        raise SkymapError(f"unreadable FITS: {exc}", code="SKYMAP_INVALID_FORMAT") from exc
    if not isinstance(hdu, fits.BinTableHDU):
        raise SkymapError("HDU 1 is not a binary table", code="SKYMAP_INVALID_FORMAT")

    h = hdu.header
    expected = int(h.get("NAXIS1", 0)) * int(h.get("NAXIS2", 0))
    if expected > len(data):
        raise SkymapError("FITS header declares more data than present (truncated file)",
                          code="SKYMAP_INVALID_FORMAT")
    pixtype = str(h.get("PIXTYPE", "HEALPIX")).upper()
    if pixtype != "HEALPIX":
        raise SkymapError(f"PIXTYPE={pixtype!r} is not HEALPIX", code="SKYMAP_INVALID_FORMAT")
    ordering = str(h.get("ORDERING", "")).upper()
    if ordering not in ("RING", "NESTED", "NUNIQ"):
        raise SkymapError(f"unknown ORDERING {ordering!r}", code="SKYMAP_INVALID_FORMAT")
    header = {k: (h[k] if not isinstance(h[k], bool) else bool(h[k])) for k in _HEADER_KEYS if k in h}
    hdr0 = hdul[0].header
    for k in ("OBJECT", "DATE-OBS", "MJD-OBS", "CREATOR", "ORIGIN"):
        if k not in header and k in hdr0:
            header[k] = hdr0[k]

    try:
        table = hdu
        _ = table.data  # force read now that sizes are checked
    except Exception as exc:
        raise SkymapError(f"cannot read table data: {exc}", code="SKYMAP_INVALID_FORMAT") from exc

    dist_cols = [_column(table, (n,)) for n in ("DISTMU", "DISTSIGMA", "DISTNORM")]
    present = [c is not None for c in dist_cols]
    if any(present) and not all(present):
        raise SkymapError("distance layers must include all of DISTMU, DISTSIGMA, DISTNORM",
                          code="SKYMAP_INVALID_FORMAT")
    distmu, distsigma, distnorm = dist_cols if all(present) else (None, None, None)

    if ordering == "NUNIQ":
        names = {c.upper() for c in table.columns.names}
        if "UNIQ" not in names:
            raise SkymapError("NUNIQ map without UNIQ column", code="SKYMAP_INVALID_FORMAT")
        uniq = np.asarray(table.data["UNIQ"], dtype=np.int64).ravel()
        rho = _column(table, ("PROBDENSITY",))
        if rho is None:
            raise SkymapError("multi-order map without PROBDENSITY column", code="SKYMAP_INVALID_FORMAT")
        order, _ = _uniq_to_order_ipix(uniq)
        if order.max() > limits.max_order:
            raise SkymapError(f"order {order.max()} exceeds limit {limits.max_order}", code="SKYMAP_TOO_LARGE")
        fmt, nside = "multiorder", None
    else:
        nside = int(h.get("NSIDE", 0))
        if nside <= 0 or nside & (nside - 1):
            raise SkymapError(f"invalid NSIDE {nside}", code="SKYMAP_INVALID_FORMAT")
        if nside > limits.max_nside:
            raise SkymapError(f"NSIDE {nside} exceeds limit {limits.max_nside}", code="SKYMAP_TOO_LARGE")
        prob = _column(table, ("PROB", "PROBABILITY", "PROBDENSITY"))
        if prob is None:
            prob = np.asarray(table.data[table.columns.names[0]], dtype=np.float64).ravel()
            diags.add("SKYMAP_PROB_COLUMN_GUESSED", f"using first column {table.columns.names[0]!r} as PROB")
        npix = 12 * nside * nside
        explicit = str(h.get("INDXSCHM", "IMPLICIT")).upper() == "EXPLICIT"
        if explicit:
            pix = np.asarray(table.data["PIXEL"], dtype=np.int64).ravel()
        else:
            if len(prob) != npix:
                raise SkymapError(f"expected {npix} pixels for NSIDE {nside}, got {len(prob)}",
                                  code="SKYMAP_INVALID_FORMAT")
            pix = np.arange(npix, dtype=np.int64)
        if np.any((pix < 0) | (pix >= npix)):
            raise SkymapError("pixel index out of range", code="SKYMAP_INVALID_FORMAT")
        if ordering == "RING":
            pix = np.asarray(ring_to_nested(pix, nside), dtype=np.int64)
        perm = np.argsort(pix, kind="stable")
        pix = pix[perm]
        order_val = int(np.log2(nside))
        uniq = 4 * (np.int64(1) << (2 * order_val)) + pix
        area = FOUR_PI / npix
        rho = prob[perm] / area
        if distmu is not None:
            distmu, distsigma, distnorm = distmu[perm], distsigma[perm], distnorm[perm]
        fmt = "flat"
    hdul.close()
    del data  # release decompressed buffer early

    # -- value checks ---------------------------------------------------------
    if np.any(~np.isfinite(rho)) or np.any(rho < 0):
        raise SkymapError("probability column contains NaN, inf or negative values",
                          code="SKYMAP_INVALID_VALUES")
    if distmu is not None:
        if np.any(np.isnan(distmu)) or np.any(np.isnan(distsigma)) or np.any(np.isnan(distnorm)):
            raise SkymapError("distance layers contain NaN", code="SKYMAP_INVALID_VALUES")
        bad = (distnorm > 0) & np.isfinite(distmu) & ~(distsigma > 0)
        if bad.any():
            raise SkymapError("DISTSIGMA must be > 0 where DISTNORM > 0", code="SKYMAP_INVALID_VALUES",
                              details={"n_bad": int(bad.sum())})
        if np.any(distnorm < 0):
            raise SkymapError("DISTNORM must be >= 0", code="SKYMAP_INVALID_VALUES")

    order, ipix = _uniq_to_order_ipix(uniq)
    _check_no_overlap(order, ipix)
    area = FOUR_PI / (12.0 * 4.0 ** order.astype(np.float64))
    total = float(np.sum(rho * area))
    if not np.isfinite(total) or total <= 0:
        raise SkymapError("sky map has zero total probability", code="SKYMAP_NOT_NORMALISED")
    dev = abs(total - 1.0)
    if dev > limits.normalisation_reject:
        raise SkymapError(f"sky map integrates to {total:.6f}", code="SKYMAP_NOT_NORMALISED",
                          details={"total": total})
    if dev > limits.normalisation_accept:
        rho = rho / total
        diags.add("SKYMAP_RENORMALISED", f"sky map integrated to {total:.6f}; renormalised", total=total)

    prior = _infer_distance_prior(header)
    meta = SkymapMeta(
        filename=filename or name,
        sha256_raw=sha_raw,
        sha256_data=sha_data,
        format=fmt,
        ordering=ordering,
        nside=nside,
        max_order=int(order.max()),
        n_pixels=len(uniq),
        header=header,
        distance_prior=prior,
    )
    if distmu is not None and prior != "r2":
        diags.add("DIST_PRIOR_ASSUMED",
                  "distance prior of this map is not known to be r^2 (CREATOR is not BAYESTAR); "
                  "assuming r^2 unless config.distance_prior is set",
                  creator=header.get("CREATOR"))
    m = MultiOrderMap(uniq, rho, distmu, distsigma, distnorm, meta=meta, diagnostics=diags)
    dclass = m.distance_class()
    if dclass == "none":
        diags.add("EVENT_NO_DIST", "no usable distance information in the 90% region; using 2D scoring")
    elif dclass == "weak":
        diags.add("EVENT_WEAK_DIST", "median sigma/mu in the 90% region exceeds 0.5")
    return m


def _check_no_overlap(order: np.ndarray, ipix: np.ndarray) -> None:
    shift = 2 * (29 - order)
    lo = ipix << shift
    hi = (ipix + 1) << shift
    srt = np.argsort(lo, kind="stable")
    if np.any(lo[srt][1:] < hi[srt][:-1]):
        raise SkymapError("multi-order pixels overlap", code="SKYMAP_INVALID_FORMAT")


def from_arrays(
    uniq: np.ndarray,
    probdensity: np.ndarray,
    distmu: np.ndarray | None = None,
    distsigma: np.ndarray | None = None,
    distnorm: np.ndarray | None = None,
    *,
    header: dict[str, Any] | None = None,
    filename: str = "<arrays>",
) -> MultiOrderMap:
    """Build a map from arrays (tests, synthetic injections).  Writes and
    re-reads through the FITS path so the same validation applies."""
    return load_skymap(to_fits_bytes(uniq, probdensity, distmu, distsigma, distnorm, header=header),
                       filename=filename)


def to_fits_bytes(uniq, probdensity, distmu=None, distsigma=None, distnorm=None, *,
                  header: dict[str, Any] | None = None) -> bytes:
    cols = [fits.Column(name="UNIQ", format="K", array=np.asarray(uniq, dtype=np.int64)),
            fits.Column(name="PROBDENSITY", format="D", unit="sr-1", array=np.asarray(probdensity, float))]
    if distmu is not None:
        cols += [fits.Column(name="DISTMU", format="D", unit="Mpc", array=np.asarray(distmu, float)),
                 fits.Column(name="DISTSIGMA", format="D", unit="Mpc", array=np.asarray(distsigma, float)),
                 fits.Column(name="DISTNORM", format="D", unit="Mpc-2", array=np.asarray(distnorm, float))]
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.header["PIXTYPE"] = "HEALPIX"
    hdu.header["ORDERING"] = "NUNIQ"
    hdu.header["INDXSCHM"] = "EXPLICIT"
    order, _ = _uniq_to_order_ipix(np.asarray(uniq, dtype=np.int64))
    hdu.header["MOCORDER"] = int(order.max())
    for k, v in (header or {}).items():
        hdu.header[k] = v
    buf = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(buf)
    return buf.getvalue()


def write_multiorder(m: MultiOrderMap, path: str | Path) -> None:
    header = dict(m.meta.header) if m.meta else {}
    for k in ("ORDERING", "INDXSCHM", "NSIDE", "MOCORDER"):
        header.pop(k, None)
    Path(path).write_bytes(to_fits_bytes(m.uniq, m.probdensity, m.distmu, m.distsigma, m.distnorm, header=header))

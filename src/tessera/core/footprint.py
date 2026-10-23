"""Instrument footprints: containment tests, rasterisation, sky polygons."""

from __future__ import annotations

from dataclasses import dataclass

import astropy.units as u
import numpy as np
from astropy_healpix import HEALPix, healpix_to_lonlat

from tessera.core.geometry import (
    angsep,
    distance_to_polygon_edges,
    gnomonic,
    inverse_gnomonic,
    points_in_polygon,
    polygon_area,
    rotate,
)
from tessera.core.profile import Footprint

_ROT_SAMPLES = 8


@dataclass
class Raster:
    order: int
    ipix: np.ndarray        # nested indices at ``order``
    weight: np.ndarray      # coverage weight in [0, 1] (fractional when rotation-averaged)
    boundary: np.ndarray    # bool: pixel centre within one pixel of the footprint edge


class InstrumentFootprint:
    def __init__(self, fp: Footprint):
        self.spec = fp
        self.kind = fp.kind
        if fp.kind == "circle":
            self.polygons: list[np.ndarray] = []
            self.radius = float(fp.radius_deg)
        elif fp.kind == "rect":
            w, h = fp.width_deg / 2, fp.height_deg / 2
            self.polygons = [np.array([[-w, -h], [w, -h], [w, h], [-w, h]], float)]
        else:
            self.polygons = [np.asarray(c, float) for c in fp.chips]
        self.rotation_averaged = bool(fp.rotation_averaged) and fp.kind != "circle"

    # -- geometry summaries -------------------------------------------------
    @property
    def bounding_radius_deg(self) -> float:
        if self.kind == "circle":
            return self.radius
        return float(max(np.max(np.hypot(p[:, 0], p[:, 1])) for p in self.polygons))

    @property
    def min_dimension_deg(self) -> float:
        if self.kind == "circle":
            return 2 * self.radius
        return float(min(min(np.ptp(p[:, 0]), np.ptp(p[:, 1])) for p in self.polygons))

    @property
    def area_deg2(self) -> float:
        if self.kind == "circle":
            r = np.radians(self.radius)
            return float(2 * np.pi * (1 - np.cos(r)) * (180 / np.pi) ** 2)
        return float(sum(polygon_area(p) for p in self.polygons))

    def tile_steps(self, overlap: float) -> tuple[float, float]:
        """(width, height) steps for tessellation: the inscribed rectangle shrunk by ``overlap``."""
        if self.kind == "circle":
            side = self.radius * np.sqrt(2)
            return side * (1 - overlap), side * (1 - overlap)
        allp = np.vstack(self.polygons)
        return float(np.ptp(allp[:, 0]) * (1 - overlap)), float(np.ptp(allp[:, 1]) * (1 - overlap))

    # -- containment --------------------------------------------------------
    def _contains_local(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        inside = np.zeros(x.shape, bool)
        for p in self.polygons:
            inside |= points_in_polygon(x, y, p)
        return inside

    def coverage(self, ra0: float, dec0: float, pa: float, ra, dec) -> np.ndarray:
        """Coverage weight in [0, 1] for each point (fraction of rotations if rotation-averaged)."""
        ra = np.asarray(ra, float)
        dec = np.asarray(dec, float)
        if self.kind == "circle":
            return (angsep(ra0, dec0, ra, dec) <= self.radius).astype(float)
        xi, eta, valid = gnomonic(ra0, dec0, ra, dec)
        if not self.rotation_averaged:
            x, y = rotate(xi, eta, -pa)
            return (self._contains_local(x, y) & valid).astype(float)
        acc = np.zeros(ra.shape)
        for k in range(_ROT_SAMPLES):
            x, y = rotate(xi, eta, -(pa + 360.0 * k / _ROT_SAMPLES))
            acc += self._contains_local(x, y) & valid
        return acc / _ROT_SAMPLES

    def contains(self, ra0: float, dec0: float, pa: float, ra, dec) -> np.ndarray:
        return self.coverage(ra0, dec0, pa, ra, dec) > 0

    def edge_distance(self, ra0: float, dec0: float, pa: float, ra, dec) -> np.ndarray:
        ra = np.asarray(ra, float)
        dec = np.asarray(dec, float)
        if self.kind == "circle" or self.rotation_averaged:
            r = self.radius if self.kind == "circle" else self.bounding_radius_deg
            return np.abs(angsep(ra0, dec0, ra, dec) - r)
        xi, eta, _ = gnomonic(ra0, dec0, ra, dec)
        x, y = rotate(xi, eta, -pa)
        return np.min([distance_to_polygon_edges(x, y, p) for p in self.polygons], axis=0)

    # -- rasterisation -------------------------------------------------------
    def raster_order(self, max_order: int = 13, min_order: int = 6) -> int:
        target = self.min_dimension_deg / 4.0
        for o in range(min_order, max_order + 1):
            res = np.degrees(np.sqrt(4 * np.pi / (12 * 4**o)))
            if res <= target:
                return o
        return max_order

    def rasterize(self, ra0: float, dec0: float, pa: float, order: int) -> Raster:
        hp = HEALPix(nside=2**order, order="nested")
        res = np.degrees(np.sqrt(4 * np.pi / (12 * 4**order)))
        cand = np.asarray(hp.cone_search_lonlat(ra0 * u.deg, dec0 * u.deg, (self.bounding_radius_deg + res) * u.deg))
        if cand.size == 0:
            return Raster(order, cand.astype(np.int64), np.zeros(0), np.zeros(0, bool))
        lon, lat = healpix_to_lonlat(cand, 2**order, order="nested")
        ra, dec = lon.to_value(u.deg), lat.to_value(u.deg)
        w = self.coverage(ra0, dec0, pa, ra, dec)
        keep = w > 0
        edge = self.edge_distance(ra0, dec0, pa, ra[keep], dec[keep]) < res
        return Raster(order, cand[keep].astype(np.int64), w[keep], edge)

    # -- sky outlines (UI, adapters) -----------------------------------------
    def sky_polygons(self, ra0: float, dec0: float, pa: float, n_edge: int = 8) -> list[list[tuple[float, float]]]:
        if self.kind == "circle" or self.rotation_averaged:
            r = self.radius if self.kind == "circle" else self.bounding_radius_deg
            t = np.linspace(0, 2 * np.pi, 48, endpoint=False)
            x, y = np.tan(np.radians(r)) * np.sin(t), np.tan(np.radians(r)) * np.cos(t)
            ra, dec = inverse_gnomonic(ra0, dec0, np.degrees(x), np.degrees(y))
            return [list(zip(ra.round(6).tolist(), dec.round(6).tolist()))]
        out = []
        for p in self.polygons:
            pts = []
            for i in range(len(p)):
                a, b = p[i], p[(i + 1) % len(p)]
                for s in np.linspace(0, 1, n_edge, endpoint=False):
                    pts.append(a + s * (b - a))
            pts = np.asarray(pts)
            x, y = rotate(pts[:, 0], pts[:, 1], pa)
            ra, dec = inverse_gnomonic(ra0, dec0, x, y)
            out.append(list(zip(ra.round(6).tolist(), dec.round(6).tolist())))
        return out

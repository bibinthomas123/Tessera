"""Spherical geometry helpers (all angles in degrees unless stated)."""

from __future__ import annotations

import numpy as np

DEG2RAD = np.pi / 180.0
SR_PER_DEG2 = DEG2RAD**2


def radec_to_xyz(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    ra = np.asarray(ra, dtype=float) * DEG2RAD
    dec = np.asarray(dec, dtype=float) * DEG2RAD
    c = np.cos(dec)
    return np.stack([c * np.cos(ra), c * np.sin(ra), np.sin(dec)], axis=-1)


def xyz_to_radec(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xyz = np.asarray(xyz, dtype=float)
    norm = np.linalg.norm(xyz, axis=-1)
    ra = np.mod(np.arctan2(xyz[..., 1], xyz[..., 0]) / DEG2RAD, 360.0)
    dec = np.arcsin(np.clip(xyz[..., 2] / norm, -1, 1)) / DEG2RAD
    return ra, dec


def angsep(ra1, dec1, ra2, dec2) -> np.ndarray:
    """Great-circle separation (Vincenty formula, numerically stable)."""
    ra1, dec1, ra2, dec2 = (np.asarray(x, dtype=float) * DEG2RAD for x in (ra1, dec1, ra2, dec2))
    dra = ra2 - ra1
    s1, c1 = np.sin(dec1), np.cos(dec1)
    s2, c2 = np.sin(dec2), np.cos(dec2)
    num = np.hypot(c2 * np.sin(dra), c1 * s2 - s1 * c2 * np.cos(dra))
    den = s1 * s2 + c1 * c2 * np.cos(dra)
    return np.arctan2(num, den) / DEG2RAD


def chord_for_angle(angle_deg: float) -> float:
    """Chord length on the unit sphere subtending ``angle_deg``."""
    return 2.0 * np.sin(0.5 * min(angle_deg, 180.0) * DEG2RAD)


def gnomonic(ra0: float, dec0: float, ra, dec) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project onto the tangent plane at (ra0, dec0).

    Returns (xi, eta, valid) in degrees; ``xi`` increases towards east (+RA),
    ``eta`` towards north.  Points on the far hemisphere are flagged invalid.
    """
    ra = np.asarray(ra, dtype=float) * DEG2RAD
    dec = np.asarray(dec, dtype=float) * DEG2RAD
    a0, d0 = ra0 * DEG2RAD, dec0 * DEG2RAD
    cosc = np.sin(d0) * np.sin(dec) + np.cos(d0) * np.cos(dec) * np.cos(ra - a0)
    valid = cosc > 1e-6
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = np.cos(dec) * np.sin(ra - a0) / cosc
        eta = (np.cos(d0) * np.sin(dec) - np.sin(d0) * np.cos(dec) * np.cos(ra - a0)) / cosc
    return xi / DEG2RAD, eta / DEG2RAD, valid


def inverse_gnomonic(ra0: float, dec0: float, xi, eta) -> tuple[np.ndarray, np.ndarray]:
    xi = np.asarray(xi, dtype=float) * DEG2RAD
    eta = np.asarray(eta, dtype=float) * DEG2RAD
    a0, d0 = ra0 * DEG2RAD, dec0 * DEG2RAD
    rho = np.hypot(xi, eta)
    c = np.arctan(rho)
    with np.errstate(divide="ignore", invalid="ignore"):
        dec = np.where(
            rho > 0,
            np.arcsin(np.cos(c) * np.sin(d0) + eta * np.sin(c) * np.cos(d0) / np.where(rho > 0, rho, 1)),
            d0,
        )
        ra = a0 + np.arctan2(xi * np.sin(c), rho * np.cos(d0) * np.cos(c) - eta * np.sin(d0) * np.sin(c))
    return np.mod(ra / DEG2RAD, 360.0), dec / DEG2RAD


def points_in_polygon(x: np.ndarray, y: np.ndarray, poly: np.ndarray) -> np.ndarray:
    """Even-odd rule point-in-polygon test, vectorised over points."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    poly = np.asarray(poly, dtype=float)
    inside = np.zeros(x.shape, dtype=bool)
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        cond = (y1 > y) != (y2 > y)
        with np.errstate(divide="ignore", invalid="ignore"):
            xint = (x2 - x1) * (y - y1) / (y2 - y1) + x1
        inside ^= cond & (x < xint)
    return inside


def distance_to_polygon_edges(x: np.ndarray, y: np.ndarray, poly: np.ndarray) -> np.ndarray:
    """Minimum Euclidean distance from points to any polygon edge (tangent-plane degrees)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    poly = np.asarray(poly, dtype=float)
    best = np.full(x.shape, np.inf)
    n = len(poly)
    for i in range(n):
        ax, ay = poly[i]
        bx, by = poly[(i + 1) % n]
        dx, dy = bx - ax, by - ay
        L2 = dx * dx + dy * dy
        t = np.clip(((x - ax) * dx + (y - ay) * dy) / (L2 if L2 > 0 else 1.0), 0.0, 1.0)
        px, py = ax + t * dx, ay + t * dy
        best = np.minimum(best, np.hypot(x - px, y - py))
    return best


def polygon_area(poly: np.ndarray) -> float:
    poly = np.asarray(poly, dtype=float)
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def polygon_is_convex(poly: np.ndarray) -> bool:
    poly = np.asarray(poly, dtype=float)
    n = len(poly)
    if n < 3:
        return False
    signs = []
    for i in range(n):
        a, b, c = poly[i], poly[(i + 1) % n], poly[(i + 2) % n]
        cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
        if abs(cross) > 1e-15:
            signs.append(cross > 0)
    return len(set(signs)) == 1


def rotate(x: np.ndarray, y: np.ndarray, pa_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Rotate tangent-plane coordinates (x east, y north) by a position angle
    measured from north through east: (0, 1) rotated by 90 deg becomes (1, 0)."""
    t = pa_deg * DEG2RAD
    c, s = np.cos(t), np.sin(t)
    return c * x + s * y, -s * x + c * y

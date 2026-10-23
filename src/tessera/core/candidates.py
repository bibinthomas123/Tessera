"""Candidate follow-up (P1): associate a detected transient with ranked galaxies
and the fields that covered it."""

from __future__ import annotations

from typing import Any

import numpy as np

from tessera.core.geometry import angsep


def associate_transient(ra: float, dec: float, galaxies: list[dict[str, Any]], fields: dict[str, list[dict[str, Any]]],
                        *, max_offset_kpc: float = 50.0, max_offset_deg_nodist: float = 0.05) -> dict[str, Any]:
    """Galaxies within a projected offset of the transient, ranked by offset
    (a simple, explicit association; not a probabilistic host model)."""
    out = []
    for g in galaxies:
        sep = float(angsep(ra, dec, g["ra"], g["dec"]))
        d = g.get("dist_mpc")
        if d and np.isfinite(d):
            kpc = np.radians(sep) * d * 1000.0
            if kpc <= max_offset_kpc:
                out.append({**{k: g.get(k) for k in ("rank", "gal_id", "dist_mpc", "S", "score_fraction")},
                            "sep_arcsec": round(sep * 3600, 2), "offset_kpc": round(float(kpc), 2)})
        elif sep <= max_offset_deg_nodist:
            out.append({**{k: g.get(k) for k in ("rank", "gal_id")}, "sep_arcsec": round(sep * 3600, 2),
                        "offset_kpc": None})
    out.sort(key=lambda r: (r["offset_kpc"] if r["offset_kpc"] is not None else 1e9, r["sep_arcsec"]))
    covering = {}
    for prof, rows in fields.items():
        covering[prof] = [r["field_id"] for r in rows if r.get("status") == "scheduled"
                          and float(angsep(ra, dec, r["ra"], r["dec"])) <= float(r.get("_radius_deg", 1.0))]
    return {"ra": ra, "dec": dec, "max_offset_kpc": max_offset_kpc, "galaxies": out[:50], "fields": covering,
            "method": "projected offset at the galaxy distance; fields matched within the footprint bounding radius"}

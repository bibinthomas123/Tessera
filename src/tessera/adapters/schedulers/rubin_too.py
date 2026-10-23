"""Rubin Observatory ToO adapter (spec 9.3).

Rubin's Feature-Based Scheduler builds its own ToO tiling and sequencing and
controls the rotator, so this adapter never sends pointings.  It emits what
``rubin_scheduler.scheduler.utils.TargetoO`` consumes:

* ``too_footprint_ring.fits``  full-sky HEALPix map, RING order, values 1
  (observe) / 0, at ``nside`` (default 32 = ``rubin_scheduler`` DEFAULT_NSIDE).
  Rubin regrids with ``hp.ud_grade`` (mean), so 1/0 masks expand forgivingly.
* ``too_footprint.multiorder.fits``  the same region with per-pixel weights
  (``raw`` sky probability or ``galaxy_boosted``) for humans and other tools.
* ``too_request.json``  ``TargetoO`` arguments: tooid, mjd_start, duration,
  ra_rad_center, dec_rad_center, too_type, posterior_distance (kpc), plus an
  eligibility verdict and the reasons for it.

ToO types follow the labels that ``rubin_scheduler`` (checked against 4.6.0)
uses in ``gen_too_surveys`` and simulates in ``targetofo.gen_events``:
``GW_case_A``–``E`` for events with an expected counterpart, by 90 % area
(<= 30, 50, 100, 250, 500 deg2), and ``BBH_case_A``–``C`` (<= 10, 20, 30 deg2).
Those area cuts are the simulation's definitions; the operational trigger
criteria and the transport into Rubin's EFD are owned by Rubin Observatory's
ToO process, so ``submit`` is deliberately not implemented.
"""

from __future__ import annotations

import io
import json
from typing import Any

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy_healpix import healpix_to_lonlat, level_to_nside, uniq_to_level_ipix

from tessera.adapters.schedulers.base import BaseAdapter, Export, Issue
from tessera.core.errors import AdapterError
from tessera.core.geometry import radec_to_xyz, xyz_to_radec
from tessera.core.profile import TelescopeProfile
from tessera.core.skymap import MultiOrderMap, to_fits_bytes

GW_CASES = (("GW_case_A", 30.0), ("GW_case_B", 50.0), ("GW_case_C", 100.0), ("GW_case_D", 250.0),
            ("GW_case_E", 500.0))
BBH_CASES = (("BBH_case_A", 10.0), ("BBH_case_B", 20.0), ("BBH_case_C", 30.0))
RUBIN_DEFAULT_NSIDE = 32


def classify(area90_deg2: float, classification: dict[str, float] | None,
             properties: dict[str, float] | None) -> tuple[str | None, list[str]]:
    """Return (too_type or None, reasons)."""
    c, p = classification or {}, properties or {}
    em_bright = (p.get("HasRemnant", 0.0) >= 0.1 or p.get("HasNS", 0.0) >= 0.5
                 or c.get("BNS", 0.0) + c.get("NSBH", 0.0) >= 0.5)
    family, cases = ("GW", GW_CASES) if em_bright else ("BBH", BBH_CASES)
    reasons = [f"{family}-type event (HasNS={p.get('HasNS')}, HasRemnant={p.get('HasRemnant')}, "
               f"BNS={c.get('BNS')}, NSBH={c.get('NSBH')}, BBH={c.get('BBH')})",
               f"90% area {area90_deg2:.1f} deg2"]
    for name, limit in cases:
        if area90_deg2 <= limit:
            reasons.append(f"area <= {limit:g} deg2 -> {name}")
            return name, reasons
    reasons.append(f"area exceeds the largest {family} case ({cases[-1][1]:g} deg2): not eligible")
    return None, reasons


def ring_footprint(skymap: MultiOrderMap, level: float, nside: int) -> np.ndarray:
    """1/0 RING-ordered map: pixels whose centre lies in the ``level`` credible region."""
    from astropy_healpix import HEALPix

    hp = HEALPix(nside=nside, order="ring")
    lon, lat = hp.healpix_to_lonlat(np.arange(hp.npix))
    idx = skymap.lookup(lon.to_value(u.deg), lat.to_value(u.deg))
    inside = skymap.region_mask_inclusive(level)
    fp = np.zeros(hp.npix)
    ok = idx >= 0
    fp[ok] = inside[idx[ok]].astype(float)
    if not fp.any():                       # region smaller than one pixel: take the peak pixel
        peak = int(np.argmax(skymap.probdensity))
        ra, dec = skymap.pixel_centers(np.array([peak]))
        fp[int(hp.lonlat_to_healpix(ra * u.deg, dec * u.deg)[0])] = 1.0
    return fp


def ring_fits_bytes(m: np.ndarray, nside: int, header: dict[str, Any]) -> bytes:
    hdu = fits.BinTableHDU.from_columns([fits.Column(name="TOO", format="E", array=m.astype(np.float32))])
    hdu.header.update({"PIXTYPE": "HEALPIX", "ORDERING": "RING", "NSIDE": nside, "INDXSCHM": "IMPLICIT",
                       "FIRSTPIX": 0, "LASTPIX": len(m) - 1, **header})
    buf = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(buf)
    return buf.getvalue()


def rubin_followed_types() -> set[str] | None:
    """ToO types followed by the installed rubin_scheduler's ``gen_too_surveys``.

    Read from the package source rather than by calling it, because building the
    surveys needs Rubin's downloaded site-model data."""
    import inspect
    import re

    try:
        from rubin_scheduler.scheduler.example import too_surveys
    except Exception:
        return None
    src = inspect.getsource(too_surveys)
    out: set[str] = set()
    for lst in re.findall(r"too_types_to_follow=\[([^\]]*)\]", src):
        out.update(re.findall(r"[\"']([A-Za-z0-9_]+)[\"']", lst))
    return out or None


class RubinToOAdapter(BaseAdapter):
    name = "rubin_too"
    version = "2"
    media_type = "application/zip"

    def export(self, schedule: dict[str, Any], profile: TelescopeProfile, *, skymap: MultiOrderMap | None = None,
               galaxies=None, event: dict[str, Any] | None = None, weighting: str = "raw", level: float = 0.9,
               nside: int = RUBIN_DEFAULT_NSIDE, duration_days: float | None = None, allow_infeasible: bool = True,
               **opts: Any) -> Export:
        if skymap is None:
            raise AdapterError("rubin_too export needs the sky map (skymap=...)", code="ADAPTER_OPTIONS_MISSING")
        if weighting not in ("raw", "galaxy_boosted"):
            raise AdapterError(f"unknown weighting {weighting!r}")
        mask = skymap.region_mask_inclusive(level)
        w = skymap.probdensity[mask].astype(float)
        if weighting == "galaxy_boosted" and galaxies is not None and len(galaxies):
            t1 = np.asarray(galaxies["tier"]) == 1
            S = np.asarray(galaxies["S"], float)[t1]
            idx = skymap.lookup(np.asarray(galaxies["ra"], float)[t1], np.asarray(galaxies["dec"], float)[t1])
            dens = np.zeros(len(skymap))
            ok = idx >= 0
            np.add.at(dens, idx[ok], S[ok])
            dens = (dens / skymap.area)[mask]
            if dens.sum() > 0:
                a = skymap.area[mask]
                w = 0.5 * w / np.sum(w * a) + 0.5 * dens / np.sum(dens * a)
        ev = event or {}
        area90 = skymap.area_deg2(0.9)
        too_type, reasons = classify(area90, ev.get("classification"), ev.get("properties"))
        fp = ring_footprint(skymap, level, nside)
        # probability-weighted centre
        r = skymap.region_mask_inclusive(0.9)
        ra, dec = skymap.pixel_centers(np.flatnonzero(r))
        cra, cdec = xyz_to_radec((radec_to_xyz(ra, dec) * skymap.prob[r][:, None]).sum(axis=0))
        from astropy.time import Time

        mjd_start = Time(schedule["horizon"]["start"].rstrip("Z"), scale="utc").mjd
        dur = duration_days or (Time(schedule["horizon"]["end"].rstrip("Z")).mjd - mjd_start)
        distmean = (skymap.meta.header.get("DISTMEAN") if skymap.meta else None)
        req = {
            "schema": "tessera.rubin-too/2",
            "rubin_scheduler_reference": "rubin_scheduler.scheduler.utils.TargetoO (checked against 4.6.0)",
            "targetoo": {
                "tooid": f"{ev.get('superevent_id', 'unknown')}-r{ev.get('revision', 0)}",
                "footprint": {"file": "too_footprint_ring.fits", "nside": nside, "ordering": "RING",
                              "values": "1 observe / 0 do not observe", "credible_level": level,
                              "n_pixels_on": int(fp.sum())},
                "mjd_start": mjd_start,
                "duration": dur,
                "ra_rad_center": float(np.radians(cra)),
                "dec_rad_center": float(np.radians(cdec)),
                "too_type": too_type,
                "posterior_distance": None if distmean is None else float(distmean) * 1000.0,   # kpc
            },
            "eligible": too_type is not None,
            "classification_reasons": reasons,
            "event": {k: ev.get(k) for k in ("superevent_id", "revision", "alert_type", "classification", "properties")},
            "area_90_deg2": area90,
            "weights_file": {"file": "too_footprint.multiorder.fits", "weighting": weighting},
            "advisory": {"filters": [s.filter for s in profile.visit_pattern],
                         "visits_planned_by_tessera": len(schedule.get("visits", [])),
                         "note": "advisory only; Rubin's scheduler builds the tiling"},
            "plan_key": schedule.get("plan_key"),
        }
        hdr = {"CREATOR": "tessera rubin_too", "TESSERALVL": level, "OBJECT": str(ev.get("superevent_id", ""))[:60]}
        return Export(self.name, self.version, {
            "too_footprint_ring.fits": ring_fits_bytes(fp, nside, hdr),
            "too_footprint.multiorder.fits": to_fits_bytes(skymap.uniq[mask], w, header={**hdr, "TESSERAWGT": weighting}),
            "too_request.json": (json.dumps(req, indent=2, sort_keys=True) + "\n").encode(),
        }, self.media_type)

    def validate(self, export: Export) -> list[Issue]:
        issues: list[Issue] = []
        req = json.loads(export.files["too_request.json"])
        t = req["targetoo"]
        with fits.open(io.BytesIO(export.files["too_footprint_ring.fits"])) as h:
            fp = np.asarray(h[1].data["TOO"], float)
            nside = int(h[1].header["NSIDE"])
            if h[1].header["ORDERING"] != "RING" or len(fp) != 12 * nside * nside:
                issues.append(Issue("error", "footprint must be a full-sky RING map"))
        if not set(np.unique(fp)) <= {0.0, 1.0} or fp.sum() == 0:
            issues.append(Issue("error", "footprint must be 1/0 with at least one pixel on"))
        with fits.open(io.BytesIO(export.files["too_footprint.multiorder.fits"])) as h:
            d = h[1].data
            if np.any(~np.isfinite(d["PROBDENSITY"])) or np.any(d["PROBDENSITY"] < 0):
                issues.append(Issue("error", "weights must be finite and >= 0"))
            lvl, ip = uniq_to_level_ipix(np.asarray(d["UNIQ"], np.int64))
            _ = healpix_to_lonlat(ip[:1], level_to_nside(int(lvl[0])), order="nested")
        if t["duration"] <= 0:
            issues.append(Issue("error", "duration must be positive"))
        if not req["eligible"]:
            issues.append(Issue("warning", "event is outside Rubin's GW/BBH ToO cases: "
                                + req["classification_reasons"][-1]))
        try:  # conformance against Rubin's scheduler package, when installed
            import healpy as hp
            from rubin_scheduler.scheduler.utils import TargetoO
        except Exception:
            issues.append(Issue("info", "rubin_scheduler/healpy not installed; scheduler conformance check skipped"))
            return issues
        try:
            too = TargetoO(t["tooid"], fp, t["mjd_start"], t["duration"], ra_rad_center=t["ra_rad_center"],
                           dec_rad_center=t["dec_rad_center"], too_type=t["too_type"] or "",
                           posterior_distance=t["posterior_distance"])
            if not np.any(hp.ud_grade(too.footprint, RUBIN_DEFAULT_NSIDE) > 0):
                issues.append(Issue("error", "footprint vanishes when regridded to Rubin's nside 32"))
            else:
                issues.append(Issue("info", "rubin_scheduler TargetoO built; footprint survives regridding"))
        except Exception as exc:
            issues.append(Issue("error", f"rubin_scheduler rejected the ToO: {exc}"))
            return issues
        followed = rubin_followed_types()
        if followed is None:
            issues.append(Issue("info", "could not read too_types_to_follow from rubin_scheduler"))
        elif t["too_type"] and t["too_type"] not in followed:
            issues.append(Issue("error", f"too_type {t['too_type']!r} is not followed by rubin_scheduler's "
                                         f"gen_too_surveys ({sorted(followed)})"))
        elif t["too_type"]:
            issues.append(Issue("info", f"too_type {t['too_type']} is followed by rubin_scheduler's gen_too_surveys"))
        return issues

    def submit(self, *args: Any, **kwargs: Any):
        raise AdapterError("Rubin ToO submission follows Rubin Observatory's ToO process and is not automated; "
                           "deliver the exported package to your Rubin contact", code="ADAPTER_NO_SUBMIT")

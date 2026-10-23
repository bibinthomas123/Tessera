"""Telescope capability profiles (spec 7.1).

Profiles are immutable once published: an edit bumps ``version``.  Loading a
profile runs JSON-Schema-level validation (Pydantic) plus physical sanity
checks, so that a typo cannot silently produce an infeasible schedule.
"""

from __future__ import annotations

import os
from importlib import resources
from pathlib import Path
from typing import Any, Literal

import numpy as np
import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator, model_validator

from tessera.core.errors import ProfileError
from tessera.core.geometry import polygon_area, polygon_is_convex
from tessera.core.hashing import sha256_json

ValidationLevel = Literal["community", "facility-attested", "measured-overheads"]


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Validation(_Strict):
    level: ValidationLevel = "community"
    by: str | None = None
    date: str | None = None
    overheads: Literal["measured", "approximate"] = "approximate"
    notes: str | None = None


class Site(_Strict):
    name: str | None = None
    lat_deg: float = Field(ge=-90, le=90)
    lon_deg: float = Field(ge=-180, le=360)
    elevation_m: float = Field(ge=-500, le=6000)
    timezone: str = "UTC"


class Axis(_Strict):
    vmax_dps: float = Field(gt=0, le=50)
    acc_dps2: float = Field(gt=0, le=50)


class Mount(_Strict):
    type: Literal["equatorial", "altaz"] = "equatorial"
    ha_limits_h: tuple[float, float] = (-12.0, 12.0)
    dec_limits_deg: tuple[float, float] = (-90.0, 90.0)
    park: dict[str, float] = Field(default_factory=lambda: {"ha_h": 0.0, "dec_deg": 0.0})

    @model_validator(mode="after")
    def _limits(self) -> Mount:
        if not (-12 <= self.ha_limits_h[0] < self.ha_limits_h[1] <= 12):
            raise ValueError("ha_limits_h must satisfy -12 <= lo < hi <= 12")
        if not (-90 <= self.dec_limits_deg[0] < self.dec_limits_deg[1] <= 90):
            raise ValueError("dec_limits_deg must satisfy -90 <= lo < hi <= 90")
        return self


class Slew(_Strict):
    """Either per-axis trapezoidal kinematics, or a simple constant + rate."""

    axes: Literal["simultaneous", "sequential"] = "simultaneous"
    axis1: Axis | None = None          # HA (equatorial) or azimuth (alt-az)
    axis2: Axis | None = None          # Dec (equatorial) or altitude (alt-az)
    rate_dps: float | None = Field(None, gt=0, le=50)
    constant_s: float = Field(0.0, ge=0, le=600)
    dome_dps: float | None = Field(None, gt=0, le=50)
    settle_s: float = Field(0.0, ge=0, le=600)
    acquisition_s: float = Field(0.0, ge=0, le=1800)

    @model_validator(mode="after")
    def _one_model(self) -> Slew:
        if (self.axis1 is None) != (self.axis2 is None):
            raise ValueError("give both axis1 and axis2, or neither")
        if self.axis1 is None and self.rate_dps is None:
            raise ValueError("slew needs either axis1/axis2 kinematics or rate_dps")
        return self


class Footprint(_Strict):
    kind: Literal["circle", "rect", "chips"]
    radius_deg: float | None = Field(None, gt=0, le=10)
    width_deg: float | None = Field(None, gt=0, le=20)
    height_deg: float | None = Field(None, gt=0, le=20)
    chips: list[list[tuple[float, float]]] | None = None
    pa_deg: list[float] = Field(default_factory=lambda: [0.0])
    rotation_averaged: bool = False

    @model_validator(mode="after")
    def _shape(self) -> Footprint:
        if self.kind == "circle" and self.radius_deg is None:
            raise ValueError("circle footprint needs radius_deg")
        if self.kind == "rect" and (self.width_deg is None or self.height_deg is None):
            raise ValueError("rect footprint needs width_deg and height_deg")
        if self.kind == "chips":
            if not self.chips:
                raise ValueError("chips footprint needs at least one chip polygon")
            for i, chip in enumerate(self.chips):
                poly = np.asarray(chip, float)
                if poly.ndim != 2 or poly.shape[0] < 3 or poly.shape[1] != 2:
                    raise ValueError(f"chip {i} must be a list of >= 3 (x, y) vertices")
                if not polygon_is_convex(poly):
                    raise ValueError(f"chip {i} is not convex")
                if np.max(np.abs(poly)) > 10:
                    raise ValueError(f"chip {i} extends more than 10 deg from the boresight")
                if polygon_area(poly) <= 0:
                    raise ValueError(f"chip {i} has zero area")
        if not self.pa_deg:
            raise ValueError("pa_deg must list at least one position angle")
        return self


class MoonModel(_Strict):
    dm_full: float = Field(0.0, ge=0, le=5, description="mag lost at full Moon when the Moon is up")


class FilterSpec(_Strict):
    m_ref: float = Field(ge=10, le=30, description="limiting mag (5 sigma) at t_ref, airmass 1, dark")
    t_ref_s: float = Field(gt=0, le=36000)
    k_ext: float = Field(0.1, ge=0, le=2, description="extinction coefficient, mag per airmass")
    moon: MoonModel = MoonModel()


class Exposure(_Strict):
    t_min_s: float = Field(gt=0)
    t_max_s: float = Field(gt=0)
    t_single_max_s: float = Field(gt=0)
    fixed_s: float | None = Field(None, gt=0)
    n_exp_fixed: int | None = Field(None, ge=1)

    @model_validator(mode="after")
    def _order(self) -> Exposure:
        if self.t_min_s > self.t_max_s:
            raise ValueError("t_min_s must be <= t_max_s")
        return self


class Instrument(_Strict):
    name: str | None = None
    footprint: Footprint
    readout_s: float = Field(ge=0, le=600)
    filter_change_s: float = Field(ge=0, le=600)
    filters: dict[str, FilterSpec]
    exposure: Exposure

    @field_validator("filters")
    @classmethod
    def _filters(cls, v: dict[str, FilterSpec]) -> dict[str, FilterSpec]:
        if not v:
            raise ValueError("at least one filter is required")
        return v


class MoonSep(_Strict):
    by_illumination: list[tuple[float, float]] | None = None
    constant: float | None = None

    def min_sep(self, illum: np.ndarray) -> np.ndarray:
        illum = np.asarray(illum, float)
        if self.by_illumination:
            tab = np.asarray(sorted(self.by_illumination), float)
            return np.interp(illum, tab[:, 0], tab[:, 1])
        return np.full(illum.shape, self.constant or 0.0)


class Constraints(_Strict):
    sun_alt_max_deg: float = Field(-12.0, ge=-30, le=0)
    airmass_max: float = Field(2.0, gt=1, le=10)
    alt_min_deg: float | None = Field(None, ge=0, le=90)
    alt_max_deg: float = Field(90.0, gt=0, le=90)
    horizon_mask: list[tuple[float, float]] | None = None
    sun_sep_min_deg: float = Field(0.0, ge=0, le=180)
    moon_sep_min_deg: MoonSep | float = 0.0
    moon_illum_max: float = Field(1.0, ge=0, le=1)
    moon_illum_only_when_up: bool = True

    @field_validator("horizon_mask")
    @classmethod
    def _mask(cls, v):
        if v is None:
            return v
        az = [p[0] for p in v]
        if az != sorted(az) or any(a < 0 or a > 360 for a in az):
            raise ValueError("horizon_mask azimuths must be increasing within [0, 360]")
        if any(not (0 <= p[1] <= 90) for p in v):
            raise ValueError("horizon_mask altitudes must be within [0, 90]")
        return v

    def moon_sep(self) -> MoonSep:
        if isinstance(self.moon_sep_min_deg, MoonSep):
            return self.moon_sep_min_deg
        return MoonSep(constant=float(self.moon_sep_min_deg))

    def min_alt(self) -> float:
        from_airmass = float(np.degrees(np.arcsin(1.0 / self.airmass_max)))
        return max(from_airmass, self.alt_min_deg or 0.0)


class Window(_Strict):
    start: str
    end: str
    label: str | None = None


class VisitStep(_Strict):
    filter: str
    min_gap_min: float = Field(0.0, ge=0)
    max_gap_min: float | None = Field(None, gt=0)


class TelescopeProfile(_Strict):
    schema_: Literal["tessera.telescope/1"] = Field("tessera.telescope/1", alias="schema")
    id: str = Field(pattern=r"^[a-z0-9][a-z0-9_.-]{1,63}$")
    version: int = Field(ge=1)
    description: str | None = None
    validation: Validation = Validation()
    site: Site
    mount: Mount = Mount()
    slew: Slew
    instrument: Instrument
    constraints: Constraints = Constraints()
    windows: list[Window] | None = None
    visit_pattern: list[VisitStep]
    field_grid: str | None = Field(None, description="CSV with field_id, ra, dec[, pa_deg]; relative to profile")
    adapters: list[str] = Field(default_factory=lambda: ["generic"])
    adapter_options: dict[str, dict[str, Any]] = Field(default_factory=dict)
    submission_enabled: bool = False
    auto_approve: dict[str, Any] | None = None

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    @model_validator(mode="after")
    def _consistency(self) -> TelescopeProfile:
        for step in self.visit_pattern:
            if step.filter not in self.instrument.filters:
                raise ValueError(f"visit_pattern filter {step.filter!r} is not in instrument.filters")
        if not self.visit_pattern:
            raise ValueError("visit_pattern must have at least one step")
        if self.visit_pattern[0].min_gap_min:
            raise ValueError("first visit step cannot have min_gap_min")
        return self

    @property
    def key(self) -> str:
        return f"{self.id}@{self.version}"

    def canonical(self) -> dict[str, Any]:
        return self.model_dump(mode="json", by_alias=True)

    @property
    def sha256(self) -> str:
        return sha256_json(self.canonical())

    # where relative paths (field_grid) are resolved
    _base_dir: Path | None = None


def _profile_dirs() -> list[Path]:
    dirs = [Path(p) for p in os.environ.get("TESSERA_PROFILE_DIRS", "").split(os.pathsep) if p]
    return dirs


def list_profiles() -> list[str]:
    names = set()
    for d in _profile_dirs():
        names.update(p.stem for p in d.glob("*.yaml"))
    pkg = resources.files("tessera.data.profiles")
    names.update(p.name[:-5] for p in pkg.iterdir() if p.name.endswith(".yaml"))
    return sorted(names)


def parse_profile(data: dict[str, Any], base_dir: Path | None = None) -> TelescopeProfile:
    try:
        prof = TelescopeProfile.model_validate(data)
    except ValidationError as exc:
        raise ProfileError("telescope profile failed validation",
                           details={"errors": exc.errors(include_url=False, include_context=False)}) from exc
    prof._base_dir = base_dir
    return prof


def load_profile(name_or_path: str | Path | dict[str, Any]) -> TelescopeProfile:
    """Load a profile from a dict, a YAML path, or a registry name (``id``)."""
    if isinstance(name_or_path, dict):
        return parse_profile(name_or_path)
    p = Path(name_or_path)
    if p.suffix in (".yaml", ".yml", ".json") and p.is_file():
        return parse_profile(yaml.safe_load(p.read_text()), base_dir=p.parent)
    name = str(name_or_path).split("@")[0]
    for d in _profile_dirs():
        cand = d / f"{name}.yaml"
        if cand.is_file():
            return parse_profile(yaml.safe_load(cand.read_text()), base_dir=d)
    pkg = resources.files("tessera.data.profiles").joinpath(f"{name}.yaml")
    if pkg.is_file():
        return parse_profile(yaml.safe_load(pkg.read_text()), base_dir=None)
    raise ProfileError(f"unknown telescope profile {name_or_path!r}", code="PROFILE_NOT_FOUND",
                       details={"known": list_profiles()})

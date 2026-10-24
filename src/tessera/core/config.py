"""Plan configuration.  Presets are expanded into an explicit, fully resolved
config whose canonical JSON is hashed into the plan key (spec 7.5)."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from tessera.core import science

Objective = Literal["sky_prob", "host_prob", "yield"]
FieldMode = Literal["auto", "grid", "tessellate", "galaxy", "hybrid"]


class SchedulerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    algorithm: Literal["greedy-v1"] = "greedy-v1"
    gain_mode: Literal["gain_rate", "gain"] = "gain_rate"
    urgency_tolerance: float = Field(0.10, ge=0, le=1)
    grid_step_s: float = Field(300.0, gt=0, le=3600)
    refine_s: float = Field(10.0, gt=0, le=300)
    require_complete_pattern: bool = False
    max_visits: int = Field(2000, ge=1)
    runners_up: int = Field(3, ge=0, le=20)
    filter_block_min: float | None = Field(None, gt=0, le=720,
                                           description="hold one filter for blocks of this many minutes")
    exposure_choice: Literal["best", "yield_rate", "depth"] = Field(
        "best", description="yield objective only. 'depth': reach the kilonova depth at the field's 90% "
                            "distance (the only behaviour before 1.2). 'yield_rate': pick each visit's "
                            "exposure from a ladder to maximise expected detections per hour. 'best': "
                            "schedule both ways and keep the schedule with more expected detections")


class FieldConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    overlap: float = Field(0.05, ge=0, lt=0.9)
    low_value_eps: float = Field(1e-4, ge=0)
    offset_fill: bool = Field(True, description="add half-tile-shifted grids as lowest-priority fields")
    low_value_rel: float = Field(0.01, ge=0, le=1,
                                 description="also keep fields above this fraction of the best field's P_sky")
    max_fields: int = Field(5000, ge=1)
    galaxy_fields_max: int = Field(1000, ge=1)
    galaxy_fields_cum_score: float = Field(0.99, gt=0, le=1)
    region_level: float = Field(0.99, gt=0, le=1)
    max_raster_order: int = Field(13, ge=6, le=16)
    hybrid_completeness_threshold: float = Field(0.5, ge=0, le=1)


class PlanConfig(BaseModel):
    """User-facing configuration.  ``None`` means "take it from the preset"."""

    model_config = ConfigDict(extra="forbid")

    preset: Literal["bns", "nsbh", "bbh", "custom"] = "bns"
    strategy: str | None = None
    objective: Objective | None = None
    field_mode: FieldMode | None = None
    light_curve: str | None = None
    horizon_nights: float | None = Field(None, gt=0, le=14)
    weight_alpha: float = Field(1.0, ge=0, le=3)

    cr_max: float | None = Field(None, gt=0, le=1)
    zscore_max: float | None = Field(None, gt=0)
    sigma_v_kms: float | None = Field(None, ge=0)
    dist_frac_floor: float | None = Field(None, ge=0)
    photoz_sigma_frac: float | None = Field(None, ge=0)
    exclude_qso: bool = True
    missing_distance: Literal["separate_tier", "exclude"] = "separate_tier"
    distance_prior: Literal["r2", "unknown"] | None = None
    budget_level: float = Field(0.999, gt=0.5, le=1)

    t_merger: str | None = Field(None, description="ISO time of merger; default MJD-OBS of the map")
    coordination: Literal["independent", "sequential"] = "independent"
    compare_strategies: bool = True
    fields: FieldConfig = FieldConfig()
    scheduler: SchedulerConfig = SchedulerConfig()

    @model_validator(mode="after")
    def _check_names(self) -> PlanConfig:
        if self.strategy is not None:
            science.strategy(self.strategy)
        if self.light_curve is not None:
            science.light_curve(self.light_curve)
        return self

    def resolved(self) -> ResolvedConfig:
        sd = science.defaults()
        p = science.preset(self.preset)
        strategy_name = self.strategy or p["strategy"]
        strat = science.strategy(strategy_name)
        return ResolvedConfig(
            preset=self.preset,
            strategy=strategy_name,
            strategy_def=strat,
            objective=self.objective or p["objective"],
            field_mode=self.field_mode or p["field_mode"],
            light_curve=self.light_curve or p["light_curve"],
            horizon_nights=self.horizon_nights or float(p["horizon_nights"]),
            galaxy_targeting=bool(p.get("galaxy_targeting", True)),
            preset_note=p.get("note"),
            weight_alpha=self.weight_alpha,
            cr_max=self.cr_max if self.cr_max is not None else sd["cr_max"],
            zscore_max=self.zscore_max if self.zscore_max is not None else sd["zscore_max"],
            sigma_v_kms=self.sigma_v_kms if self.sigma_v_kms is not None else sd["peculiar_velocity_kms"],
            dist_frac_floor=self.dist_frac_floor if self.dist_frac_floor is not None else sd["dist_frac_floor"],
            photoz_sigma_frac=(self.photoz_sigma_frac if self.photoz_sigma_frac is not None
                               else sd["photoz_sigma_frac"]),
            w_impute_percentile=sd["w_impute_percentile"],
            cosmology=sd["cosmology"],
            exclude_qso=self.exclude_qso,
            missing_distance=self.missing_distance,
            distance_prior=self.distance_prior,
            budget_level=self.budget_level,
            t_merger=self.t_merger,
            coordination=self.coordination,
            compare_strategies=self.compare_strategies,
            fields=self.fields,
            scheduler=self.scheduler,
            science_model=sd["science_model"],
        )


class ResolvedConfig(BaseModel):
    """Fully expanded configuration; this is what the manifest records."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    preset: str
    strategy: str
    strategy_def: dict[str, Any]
    objective: Objective
    field_mode: FieldMode
    light_curve: str
    horizon_nights: float
    galaxy_targeting: bool
    preset_note: str | None
    weight_alpha: float
    cr_max: float
    zscore_max: float
    sigma_v_kms: float
    dist_frac_floor: float
    photoz_sigma_frac: float
    w_impute_percentile: float
    cosmology: str
    exclude_qso: bool
    missing_distance: str
    distance_prior: str | None
    budget_level: float
    t_merger: str | None
    coordination: str
    compare_strategies: bool
    fields: FieldConfig
    scheduler: SchedulerConfig
    science_model: str

    @property
    def weight(self) -> str:
        return self.strategy_def["weight"]

    @property
    def uses_distance(self) -> bool:
        return bool(self.strategy_def["distance"])

    @property
    def host_posterior(self) -> bool:
        return self.strategy_def["output"] == "host_probability"

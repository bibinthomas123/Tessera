# Science model (science-model 1.2.0)

The full specification, with rationale, lives in the design document; this page is the
implementation reference. Every output records the strategy, preset and science-model version.

## Sky map

All inputs become a multi-order map: UNIQ, density ρ (sr⁻¹) and the distance ansatz
(μ, σ, N) of Singer et al. (2016):

    p(r | n) = N r² Gauss(r; μ, σ)

Pixels with μ = +∞ or N = 0 carry no distance information. Normalisation within 1e-3 is
accepted, within 1e-2 renormalised (`SKYMAP_RENORMALISED`), otherwise rejected.

## Galaxy score

    p_sky,i  = ρ of the containing pixel                    [deg⁻²]  (no interpolation)
    p_dist,i = N Gauss(d_i; μ, sqrt(σ² + s_i²))             [Mpc⁻³]
    w_i      = tracer^α
    S_i      = p_sky,i · p_dist,i · w_i                      (ranking score, not a probability)

`p_dist` is the GW posterior density per unit volume divided by the sky density, marginalised over a
Gaussian galaxy-distance error. It assumes the map's distance prior ∝ r² (true for BAYESTAR);
other creators raise `DIST_PRIOR_ASSUMED`.

Galaxy distance error: s² = s_meas² + s_pv², with s_pv = σ_v / H₀ (σ_v = 300 km/s, Planck18)
added only to redshift-derived distances whose catalog error does not already include it.
`credible_level` is the searched probability (agrees with `ligo.skymap` crossmatch to 1e-12).

`dist_zscore = (d − μ)/sqrt(σ² + s²)`; galaxies with |z| > 5 or credible level > 0.99 are pre-cut
for efficiency, and the score mass removed is reported as `p_precut`.

## Weights

| weight | tracer | default for |
| --- | --- | --- |
| `flat` | 1 | BBH |
| `lum_K` | K-band luminosity (W1 fallback) | BNS, NSBH |
| `mstar` | catalog stellar mass | — |
| `lum_B` | B-band luminosity | comparison |
| `bns_rate` | GLADE+ logRate (experimental) | — |

A missing tracer is imputed with the 5th percentile (lower method) of the candidate set, so it can
never outrank a measured galaxy (`W_IMPUTED`).

## Tiers and flags

Tier 1 = galaxies with a distance; tier 2 = no distance, ranked by p_sky only and never interleaved.
Flags: `NO_DIST`, `PHOTOZ`, `LOW_Z_PV`, `DIST_FROM_Z`, `PIX_NO_DIST`, `W_IMPUTED`, `W1_FOR_K`,
`OUTSIDE_CR90`. QSOs and suspected duplicates are excluded (counts in the summary).

## Host probability (`host-posterior`, experimental)

    P_i = S_i / (Σ S_k + U),   U = ρ̄_w (1 − P_cat) · (π/180)²

ρ̄_w is the reference weight density in `data/science/defaults.yaml`. Galaxy probabilities plus
`p_uncat_host_model` sum to one.

| weight | ρ̄_w (h = 0.7) | source | status |
| --- | --- | --- | --- |
| `lum_B` | 1.98×10⁸ L_B☉ Mpc⁻³ | Kopparapu et al. 2008 | final |
| `lum_K` | 4.01×10⁸ L_K☉ Mpc⁻³ | Kochanek et al. 2001, Table 3 full-sample fit (K* = −23.39 + 5 log h, α = −1.09, n* = 1.16×10⁻² h³ Mpc⁻³), integrated as n* L* Γ(α+2) with M☉,K = 3.28 | final (derived) |
| `mstar` | 2.3×10⁸ M☉ Mpc⁻³ | Baldry et al. 2012 (GAMA GSMF, z < 0.06, Chabrier IMF), §3; their double-Schechter fit integrates to 2.24×10⁸ | final: same mass scale as GLADE+ (see below) |

**Why GAMA for `mstar`.** GLADE+ (Dalya et al. 2022, §2.3) derives stellar masses from WISE W1
using the mass-to-light ratios of Kettlety et al. 2018 (M/L_W1 = 0.65 for passive galaxies, a
colour relation for star-forming ones). Kettlety et al. calibrate those ratios against GAMA SED
masses, which use a Chabrier IMF. The GAMA stellar mass density is therefore on the same mass
scale as the catalog masses it is compared with. Bell et al. 2003, used before 1.2.0, is an upper
limit for a bottom-heavy IMF and was about 1.7 times too high on this scale.

## Catalog budget

The completeness cube C_w(pix₄, shell) = cataloged weight / (ρ̄_w · shell volume), on ~13 deg²
cells and 10 Mpc shells (logarithmic beyond 200 Mpc), averaged over each cell and its 8 neighbours
in the same shell, then capped at 1. Order 4 and the smoothing are used because finer or unsmoothed
cells are dominated by Poisson noise from single bright galaxies, which capping turns into a low
bias. With smoothing, GW170817 (≈ 40 Mpc) gets P_cat = 0.99 from GLADE+, consistent with GLADE+
being complete to 47 Mpc in B-band luminosity (Dalya et al. 2022); without it, 0.64.

Completeness is measured against a fixed cosmic mean density, so local over- or under-density moves
P_cat too: an underdense volume looks incomplete even if every galaxy is cataloged.

    P_cat = Σ_pix prob · Σ_shell P(shell|pix) · min(1, C)
    P_uncat = P_in_footprint − P_cat

Probability outside the 0.999 region is `p_unassessed`, never attributed to the catalog.

## Detectability

    m_lim = m_ref + 1.25 log10(t/t_ref) − k(X − 1) − dm_full · illum · [Moon up]
    P_det = Φ((m_lim − M(t) − A − E[DM]) / sqrt(Var[DM] + σ_M²))

The distance-modulus distribution in each pixel is summarised by its mean and variance (Gaussian
approximation).

Default light curve `at2017gfo` (since science model 1.1.0): binned medians of the u, g, r, i, z, y
detections in Villar et al. 2017 (ApJL 851, L21) Table 3, corrected for E(B−V) = 0.105 and placed at
39.5 Mpc as in that paper (`data/science/at2017gfo_absmag.csv`; the extracted rows are in
`at2017gfo_villar2017_ugrizy.csv`). Interpolated in log time; flat before 0.47 d; extrapolated with
the last measured decline.

σ_M = 1.0 mag is the population spread around AT2017gfo. One kilonova cannot measure it; it is taken
from Ascenzi et al. 2019 (MNRAS 486, 672), whose fits to sGRB-associated kilonovae give 95% peak
ranges of [−16.8, −12.3] in g and [−16.2, −13.1] in H, i.e. σ ≈ 1.15 and 0.79 mag. Gompertz et al.
2018 and Kasliwal et al. 2020 (O3 non-detections) also argue against all kilonovae being as bright
as AT2017gfo. `tessera benchmark` re-plans with σ_M = 0.25, 0.5 and 1.0 (section "Sensitivity to
kilonova diversity"): GW170817's schedule is unchanged, and for S190425z (155 Mpc) expected
detections go from 0.53 to 0.46 while about 95% of the scheduled fields stay the same.

## Changes in 1.2.0

* `mstar` reference density 2.3×10⁸ (GAMA, Chabrier; was the Bell et al. 2003 upper limit 3.9×10⁸),
  matched to the mass scale of GLADE+ stellar masses. Affects the host-posterior uncataloged term
  and budgets under the `mstar` weight only.
* `at2017gfo` σ_M 1.0 mag (was an assumed 0.5), sourced from Ascenzi et al. 2019. Lowers expected
  yields for distant events; schedules change little.
* GW170817 golden ranking and schedule are byte-identical apart from the science-model version.

## Changes in 1.1.0

* Tabulated AT2017gfo light curve replaces the parametric `kn-gfo-approx` in the presets.
* Sourced `lum_K` (4.01×10⁸, was a provisional 5.0×10⁸) and `mstar` reference densities.
* Completeness cube smoothed over neighbouring cells before capping (see above). This changes
  budgets and, through the hybrid strategy's low-coverage mask, which areas are tiled.
* Rankings are unchanged (the GW170817 golden ranking is byte-identical). Budgets and yield-driven
  schedules can change, hence a minor bump.

## Scheduling

Stage A builds fields (`grid`, `tessellate`, `galaxy`, `hybrid`), rasterised so a pixel is at most a
quarter of the smallest chip dimension. Stage B is `greedy-v1`: gain on the union of coverage,
divided by (gap + visit time); among candidates within 10% of the best, the one whose window closes
first wins; ties go to higher gain, then field ID. Only visits fitting entirely inside a hard-
constraint window are candidates, and an independent validator rechecks every visit at 10 s.
Optional filter blocking (`scheduler.filter_block_min`) holds one filter for blocks of that many
minutes.

**Exposure times (yield objective).** `scheduler.exposure_choice`:

* `depth` gives each visit the exposure that reaches the kilonova at the field's 90% distance,
  clamped to the profile's limits. This was the only behaviour in 1.1.0 and earlier.
* `yield_rate` tries a ladder (the depth time, then t_min, 2 t_min, … up to t_max) and keeps the
  exposure with the most expected detections per hour.
* `best` (the default) schedules the night both ways from the same starting coverage, and keeps
  the schedule with more expected detections.

The yield-rate choice is myopic. It wins when there is more sky than time, but loses when the
visible region can be covered anyway and depth is what is left to buy. Examples:

| Telescope and event | Expected detections, depth | Expected detections, yield_rate | `best` keeps |
| --- | --- | --- | --- |
| example 1 m, S190425z | 0.013 | 0.017 (60 s exposures instead of 100–600 s) | yield_rate |
| example 1 m, S200105ae | 0.006 | 0.012 | yield_rate |
| example wide-field, S190425z | 0.457 | 0.413 | depth |

`schedule.summary.exposure_choice` records which was kept.

**Shifted pointings and schedule variants.**

* *Offset fields.* With `fields.offset_fill` (the default), grids shifted by half a tile are added
  as lowest-priority `offset` fields. For fields of 2° or more, quarter-tile shifts are added too.
  Fixed tiles whose centres never become observable leave the visible edge of the region
  uncovered; shifted pointings reach it.
* *Scheduling variants.* Taking an offset pointing during a short lull can spoil tiles that become
  visible later, so the scheduler runs every variant from the same starting coverage: with and
  without offset fields, and, under the yield objective, with each exposure policy. It keeps the
  variant with the best value of the plan's own objective, and
  `summary.variants_compared` lists them all.
* *Cost.* Planning takes about 3× longer than 1.2.0 (GW170817, 1 m: 2 s → 6 s). Setting
  `fields.offset_fill: false` and `scheduler.exposure_choice: depth` restores the 1.1 speed.

**Detection accounting.** Several looks at the same pixel are not independent chances to detect: a
counterpart too faint for one exposure is too faint for another of the same depth. So each pixel
keeps the best detection probability of any look. Before 1.3 the miss probabilities were multiplied,
which overstated the yield of overlapping and repeated looks.

**Field limits.** A field is rejected as low value when its P_sky is below `fields.low_value_eps`
(1e-4) **and** below `fields.low_value_rel` (1%) of the best field's P_sky. Without the relative
term, very large maps lost almost all their fields: on S230529ay (24,000 deg²) every 0.25 deg² tile
holds about 1e-5, and 1.1.0 scheduled no tiles at all. When a region has more tiles than
`fields.max_fields`, tiles the site can never observe are dropped before the cap: their
declination keeps them below the minimum altitude, or outside the mount's limits. The drop is
reported as `FIELDS_UNREACHABLE_FROM_SITE` and `REJECT_UNREACHABLE_FROM_SITE`.

In sequential coordination, a later telescope's `cum_p_sky` includes earlier telescopes' coverage
only on pixels inside its own fields. It is not the combined coverage of all telescopes.

`hybrid` fields are galaxy-targeted fields plus tiles where the catalog is incomplete (the primary
fields), and tiles over the rest of the credible region with origin `fill`. A fill tile is a
candidate only when no primary field has positive gain at that moment, so galaxies keep their
priority and the telescope does not sit idle once they are done (decision rule `fill: ...`). In
1.1.0 and earlier there were no fill tiles: on GW170817 the example 1 m telescope stopped after 25 visits and
16% of the sky probability, with an hour of usable time left. It now makes 69 visits and covers
72%, with the same first 25 visits and NGC 4993 still on visit 6.

`tessera.core.optimality` solves a MILP relaxation (time buckets, minimum gaps, capped coverage) whose
optimum bounds every feasible schedule from above, so greedy / bound is a guaranteed lower bound on
greedy / optimum. On GW170817 with the example 1 m profile it is ≥ 0.91.

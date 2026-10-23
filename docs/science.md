# Science model (science-model 1.1.0)

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
| `mstar` | 3.9×10⁸ M☉ Mpc⁻³ | Bell et al. 2003: Ω* h = 0.0020 ± 0.0006, an upper limit for a bottom-heavy IMF | provisional: IMF consistency with GLADE+ masses unverified |

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

Default light curve `at2017gfo` (science model 1.1.0): binned medians of the u, g, r, i, z, y
detections in Villar et al. 2017 (ApJL 851, L21) Table 3, corrected for E(B−V) = 0.105 and placed at
39.5 Mpc as in that paper (`data/science/at2017gfo_absmag.csv`; the extracted rows are in
`at2017gfo_villar2017_ugrizy.csv`). Interpolated in log time; flat before 0.47 d; extrapolated with
the last measured decline. σ_M = 0.5 mag for kilonova diversity is an assumption, not a measurement.

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

`tessera.core.optimality` solves a MILP relaxation (time buckets, minimum gaps, capped coverage) whose
optimum bounds every feasible schedule from above, so greedy / bound is a guaranteed lower bound on
greedy / optimum. On GW170817 with the example 1 m profile it is ≥ 0.91.

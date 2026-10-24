# tessera benchmark report (tessera 1.2.0, science model 1.2.0)

## Headline

| Check | Result |
| --- | --- |
| ngc4993_rank_vol3d_lumK | 1 |
| ngc4993_rank_sky2d | 14 |
| gw170817_plan_wall_s | 6.41 |
| top10_recovery_vol3d_lumK | 0.43 |
| top10_recovery_sky2d | 0.07 |
| 3d_beats_2d_at_all_k | True |
| budget_abs_difference | 0.008939455013649455 |
| all_schedules_feasible | True |
| deterministic | True |
| malformed_all_rejected | True |
| edge_cases_all_handled | True |
| greedy_vs_bound_min_ratio | 0.9072437161190452 |

## Synthetic injections (top-k host recovery)

200 injections, seed 20170817, 10363 synthetic galaxies.

| Strategy | top1 | top10 | top100 |
| --- | --- | --- | --- |
| vol3d-lumK | 0.14 | 0.43 | 0.84 |
| vol3d-flat | 0.15 | 0.33 | 0.81 |
| sky2d-flat | 0.00 | 0.07 | 0.57 |

Prior mismatch (hosts drawn flat, ranked with L_K): top1 0.04, top10 0.16, top100 0.78.

Budget accuracy (catalog halved): predicted P_cat 0.447 (full catalog 0.883; this synthetic volume is not exactly at the mean density), halved/full 0.506 vs truth 0.497; empirical host-in-catalog rate 0.580 over 100 injections.

## Historical events

### GW170817

| Strategy | NGC 4993 rank | score fraction |
| --- | --- | --- |
| vol3d-lumK | 1 | 0.493 |
| vol3d-mstar | 1 | 0.557 |
| vol3d-lumB | 1 | 0.224 |
| vol3d-flat | 2 | 0.077 |
| sky2d-flat | 14 | 0.018 |

Credible level of NGC 4993: 0.383. First example-1m visit covering it: {'start': '2017-08-17T23:21:08Z', 'seq': 6, 'field_id': 'g:GLADE+:1151336', 'hours_after_merger': 10.67}. Reference: SSS17a discovered ~10.9 h after merger by Swope (Coulter et al. 2017).

### S190425z

90% area 10183 deg2, wall 54.2 s.

| Profile | status | visits | cum P_sky | mode |
| --- | --- | --- | --- | --- |
| example-wide-north@1 | proposed | 380 | 0.586 | tessellate |
| rubin-lsstcam-approx@1 | proposed | 949 | 0.691 | tessellate |

### S200105ae

90% area 7719 deg2, wall 47.31 s.

| Profile | status | visits | cum P_sky | mode |
| --- | --- | --- | --- | --- |
| example-wide-north@1 | proposed | 156 | 0.622 | tessellate |
| rubin-lsstcam-approx@1 | proposed | 708 | 0.782 | tessellate |

## Sensitivity to kilonova diversity sigma_M

Adopted sigma_M = 1.0 mag (Ascenzi et al. 2019). Field overlap is the Jaccard index of scheduled fields against the adopted value.

| Event | profile | sigma_M | visits | expected detections | cum P_sky | field overlap |
| --- | --- | --- | --- | --- | --- | --- |
| GW170817 | example-1m | 0.25 | 69 | 0.723 | 0.723 | 1.00 |
| GW170817 | example-1m | 0.5 | 69 | 0.723 | 0.723 | 1.00 |
| GW170817 | example-1m | 1.0 | 69 | 0.722 | 0.723 | 1.00 |
| S190425z | example-wide-north | 0.25 | 386 | 0.525 | 0.587 | 0.96 |
| S190425z | example-wide-north | 0.5 | 387 | 0.499 | 0.588 | 0.96 |
| S190425z | example-wide-north | 1.0 | 380 | 0.443 | 0.586 | 1.00 |

## Greedy vs. provable upper bound (GW170817, example-1m, 1 night)

| fields / objective | greedy | upper bound | greedy / bound (>= greedy / optimum) |
| --- | --- | --- | --- |
| tessellate/sky_prob | 0.8448 | 0.9312 | 0.907 |
| tessellate/host_prob | 0.9998 | 1.0000 | 1.000 |
| galaxy/sky_prob | 0.1551 | 0.1551 | 1.000 |
| galaxy/host_prob | 0.9917 | 0.9917 | 1.000 |

## Engineering

Deterministic: True

Malformed inputs: empty=SKYMAP_INVALID_FORMAT, not_fits=SKYMAP_INVALID_FORMAT, truncated=SKYMAP_INVALID_FORMAT, gzip_garbage=SKYMAP_INVALID_FORMAT, gzip_bomb=SKYMAP_TOO_LARGE

Runtime of this report: 362.9 s

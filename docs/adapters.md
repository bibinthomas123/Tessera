# Scheduler adapters

All adapters read `tessera.schedule/1`, refuse infeasible schedules, and refuse submission without an
approval bound to the exact version and `submission_enabled: true` in the profile. Third-party
adapters register under the `tessera.adapters` entry-point group and implement `export`, `validate`
and optionally `submit`.

## Verification status

| Adapter | Checked against | Still to confirm with the facility | Submission |
| --- | --- | --- | --- |
| `generic` | tessera's own schema | — | — |
| `rubin_too` | installed `rubin-scheduler` 4.6.0: `TargetoO` is built from the export, the RING footprint survives `hp.ud_grade` to nside 32, and `too_type` is one that `gen_too_surveys` follows | the operational trigger criteria (the GW/BBH area cases come from Rubin's event simulator), and how a package enters Rubin's EFD | never (by design) |
| `lco` | open-source Observation Portal models and serializers (operator, observation types, ipp_value range, name lengths, ICRS degrees, constraint ranges, required acquisition/guiding configs) | instrument types, telescope classes, filter codes and guiding/acquisition modes against a snapshot of `https://observe.lco.global/api/instruments/` (`data/lco_instruments.json`, refresh with `tessera.adapters.schedulers.lco.refresh_instruments()`); still to confirm: endpoint path and `Authorization: Token` header, via LCO's validate endpoint with an account | `TESSERA_LCO_TOKEN` |
| `treasuremap` | official client `gwtm_api` (base URL `/api/v1`, batch `pointings` body, `api_token` HTTP header, bandpass and depth-unit enums, integer `instrumentid`) | your `instrumentid`; how planned pointings are cancelled (no client call exists) | `TESSERA_TREASUREMAP_TOKEN` |
| `acp` | **not verified**: the ACP help pages could not be reached; every export carries a warning | the whole directive syntax and date format, against your ACP version; then set `adapter_options.acp.syntax_confirmed: true` | no (deliver the file) |

## Options

Profile options live under `adapter_options.<name>`:

```yaml
adapter_options:
  lco: {telescope_class: 1m0, instrument_type: 1M0-SCICAM-SINISTRO, proposal: MY-PROP,
        filter_map: {r: rp, i: ip}, acquisition_mode: "OFF", guiding_mode: ON}
  treasuremap: {instrumentid: 12, band_map: {w: clear}}
  acp: {filter_map: {r: R}, binning: 1, syntax_confirmed: false}
```

## Rubin

The Rubin adapter never sends pointings; Rubin's scheduler builds its own tiling and controls the
rotator. It exports the `TargetoO` arguments (`tooid`, `mjd_start`, `duration`, `ra_rad_center`,
`dec_rad_center`, `too_type`, `posterior_distance` in kpc), a full-sky RING 1/0 footprint, and a
weighted multi-order copy for people. `too_type` follows `rubin-scheduler`'s labels:

| Event | 90% area | Type |
| --- | --- | --- |
| BNS/NSBH-like (HasNS ≥ 0.5, HasRemnant ≥ 0.1, or P(BNS)+P(NSBH) ≥ 0.5) | ≤ 30 / 50 / 100 / 250 / 500 deg² | `GW_case_A` … `GW_case_E` |
| otherwise (BBH-like) | ≤ 10 / 20 / 30 deg² | `BBH_case_A` … `BBH_case_C` |
| larger | — | not eligible (reported, not exported as a type) |

Install the conformance check with `pip install "tessera-gw[rubin]"`.

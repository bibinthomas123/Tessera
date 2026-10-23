# Contributing

1. `python -m venv .venv && .venv/bin/pip install -e ".[dev,validation]"`
2. `ruff check src tests && pytest` must pass. CI also runs the quick benchmark.
3. Keep `tessera.core` free of network, database and web imports.

## Science-model changes (SMEP)

Any change to formulas, priors or defaults that alters rankings is a Science Model Enhancement
Proposal: open an issue with the SMEP template, get approval from two science maintainers, bump
`SCIENCE_MODEL_VERSION` (major if rankings change), regenerate the golden file with
`TESSERA_UPDATE_GOLDEN=1 pytest tests/regression`, and attach the benchmark diff. No major or minor
science-model change is merged into an observing-run LTS branch.

## Telescope profiles

Submit profiles with the "Telescope profile" issue template. `community` profiles need only pass
`tessera profiles --check`; `facility-attested` requires sign-off from an operator of that facility;
`measured-overheads` additionally requires measured slew/readout/filter times.

## Releases

Semver for software; separate semver for the science model. Deprecations are announced one
release ahead. Each release publishes the benchmark report and a signed container image.

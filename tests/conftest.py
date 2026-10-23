from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

DATA = Path(__file__).parent / "data"
GW170817_MAP = DATA / "GW170817_skymap.fits.gz"
S200105AE_MAP = DATA / "S200105ae_bayestar.multiorder.fits"
S190425Z_MAP = DATA / "S190425z_bayestar.fits.gz"
GLADE_EXCERPT = DATA / "glade_plus_gw170817.tsv"
NGC4993 = "GLADE+:1151336"
T_GW170817 = "2017-08-17T22:00:00"


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA


@pytest.fixture(scope="session")
def gw170817_map():
    from tessera.core.skymap import load_skymap

    return load_skymap(GW170817_MAP)


@pytest.fixture(scope="session")
def s200105ae_map():
    from tessera.core.skymap import load_skymap

    return load_skymap(S200105AE_MAP)


@pytest.fixture(scope="session")
def glade_cols():
    from tessera.adapters.catalogs import glade_plus

    return glade_plus.read(GLADE_EXCERPT)


@pytest.fixture(scope="session")
def glade_snapshot(tmp_path_factory, glade_cols):
    from tessera.core.catalog import build_snapshot

    d = tmp_path_factory.mktemp("catalogs") / "gladeplus-excerpt"
    return build_snapshot(glade_cols, d, name="gladeplus-excerpt", upstream_version="2022-10-31")


@pytest.fixture(scope="session")
def gw170817_plan(gw170817_map, glade_snapshot):
    from tessera import PlanConfig, plan

    return plan(gw170817_map, glade_snapshot, ["example-1m"], PlanConfig(preset="bns", horizon_nights=1),
                T_GW170817, event={"superevent_id": "GW170817", "revision": 1})


@pytest.fixture(scope="session")
def synthetic_map():
    from tessera.testing import gaussian_map

    return gaussian_map([(197.45, -23.38, 2.0, 1.0)], order=8, dist=(40.0, 8.0))


@pytest.fixture()
def rng():
    return np.random.default_rng(1234)

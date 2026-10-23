"""Regression against a stored golden result for GW170817 (spec 11.2).

A science-model change that alters rankings must update the golden file in the
same pull request (and bump the science-model version).  Regenerate with::

    TESSERA_UPDATE_GOLDEN=1 pytest tests/regression
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from tessera import SCIENCE_MODEL_VERSION

GOLDEN = Path(__file__).parent / "golden_gw170817.json"


def _extract(plan) -> dict:
    t = plan.galaxies
    s = plan.schedules["example-1m@1"]
    return {
        "science_model": SCIENCE_MODEL_VERSION,
        "top20": [{"gal_id": str(r["gal_id"]), "S": float(r["S"]), "credible_level": float(r["credible_level"])}
                  for r in t[:20]],
        "n_scored": len(t),
        "budget_p_cat": plan.budget["p_cat"],
        "visits": [{k: v[k] for k in ("field_id", "filter", "n_exp", "t_exp_s", "start", "end")}
                   for v in s["visits"][:15]],
        "n_visits": len(s["visits"]),
    }


def test_gw170817_matches_golden(gw170817_plan):
    got = _extract(gw170817_plan)
    if os.environ.get("TESSERA_UPDATE_GOLDEN") == "1" or not GOLDEN.exists():
        GOLDEN.write_text(json.dumps(got, indent=2) + "\n")
        pytest.skip("golden file written")
    want = json.loads(GOLDEN.read_text())
    assert got["science_model"] == want["science_model"], "science model changed: regenerate the golden file"
    assert [g["gal_id"] for g in got["top20"]] == [g["gal_id"] for g in want["top20"]]
    np.testing.assert_allclose([g["S"] for g in got["top20"]], [g["S"] for g in want["top20"]], rtol=1e-9)
    np.testing.assert_allclose([g["credible_level"] for g in got["top20"]],
                               [g["credible_level"] for g in want["top20"]], rtol=1e-12)
    assert got["n_scored"] == want["n_scored"]
    assert got["budget_p_cat"] == pytest.approx(want["budget_p_cat"], rel=1e-9)
    assert got["visits"] == want["visits"] and got["n_visits"] == want["n_visits"]


def test_ngc4993_observed_on_the_discovery_night(gw170817_plan):
    """The field containing NGC 4993 is scheduled before 00:30 UTC on 18 Aug 2017 (it set ~1.5 h after
    nautical twilight at Las Campanas; SSS17a was found at 23:33 UTC, Coulter et al. 2017)."""
    t = gw170817_plan.galaxies
    row = t[np.asarray(t["gal_id"]).astype(str) == "GLADE+:1151336"][0]
    fids = {f.split(":", 1)[1] for f in str(row["field_ids"]).split("|") if f}
    visits = [v for v in gw170817_plan.schedules["example-1m@1"]["visits"] if v["field_id"] in fids]
    assert visits and visits[0]["start"] < "2017-08-18T00:30:00Z"

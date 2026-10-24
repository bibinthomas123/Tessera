"""How close is greedy-v1 to optimal?  A provable upper bound via MILP (spec 6.2, 11.4).

Exact sequencing with slews is a hard orienteering problem, so instead of
solving it we solve a *relaxation* whose optimum can only be larger than that
of any feasible schedule:

* time is cut into buckets; a visit is assigned to the bucket containing its
  start, and may use it only if the field's feasibility window overlaps it;
* the visits assigned to a bucket must fit in the bucket length plus the
  longest single visit (a visit may run past the bucket edge), each charged
  its exposure + readout + the minimum possible gap (settle + acquisition,
  zero slew, no filter change);
* coverage of a pixel (or galaxy) is at most the sum of the weights of the
  chosen fields covering it, and at most 1;
* each field is counted once (revisits only consume time, never add gain).

Every real schedule maps onto a feasible point of this model, so
greedy / bound is a lower bound on greedy / optimum.  If the solver stops at
its time limit the reported bound is the solver's best bound, still valid.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from tessera.core.fields import FieldSet
from tessera.core.scheduler import DAY, GreedyScheduler


def upper_bound(sched: GreedyScheduler, *, bucket_min: float = 30.0, max_fields: int = 1000,
                time_limit_s: float = 60.0) -> dict[str, Any]:
    try:
        from ortools.linear_solver import pywraplp
    except ImportError:  # pragma: no cover - optional dependency
        return {"status": "skipped", "reason": "ortools not installed (pip install ortools)"}
    fs: FieldSet = sched.fs
    vis = sched.vis
    objective = sched.objective
    if objective not in ("sky_prob", "host_prob"):
        return {"status": "skipped", "reason": f"bound implemented for sky_prob/host_prob, not {objective}"}
    rows = np.array([i for i in fs.active() if vis.windows[i]], dtype=int)
    if len(rows) == 0:
        return {"status": "skipped", "reason": "no visible candidate fields"}
    if len(rows) > max_fields:
        return {"status": "skipped", "reason": f"{len(rows)} fields exceeds max_fields={max_fields}"}

    filt = fs.profile.visit_pattern[0].filter
    starts = np.array([vis.windows[i][0][0] for i in rows])
    n_exp, t_single = sched._exposures(rows, starts, filt)
    dur = n_exp * (t_single + fs.profile.instrument.readout_s) / DAY
    g_min = (fs.profile.slew.settle_s + fs.profile.slew.acquisition_s) / DAY
    cost = dur + g_min
    edges = np.arange(sched.t0, sched.t1 + 1e-9, bucket_min / 1440.0)
    if edges[-1] < sched.t1:
        edges = np.append(edges, sched.t1)

    solver = pywraplp.Solver.CreateSolver("SCIP") or pywraplp.Solver.CreateSolver("CBC")
    if solver is None:  # pragma: no cover
        return {"status": "skipped", "reason": "no MILP backend available in ortools"}
    solver.SetTimeLimit(int(time_limit_s * 1000))
    x: dict[tuple[int, int], Any] = {}
    for j, i in enumerate(rows):
        for b in range(len(edges) - 1):
            lo, hi = edges[b], edges[b + 1]
            if any(wa < hi and wb > lo and wb - max(wa, lo) >= dur[j] for wa, wb in vis.windows[i]):
                x[j, b] = solver.BoolVar(f"x_{j}_{b}")
    if not x:
        return {"status": "skipped", "reason": "no field fits a bucket"}
    for j in range(len(rows)):
        vs = [v for (jj, _), v in x.items() if jj == j]
        if vs:
            solver.Add(solver.Sum(vs) <= 1)
    dmax = float(dur.max())
    for b in range(len(edges) - 1):
        vs = [(v, cost[j]) for (j, bb), v in x.items() if bb == b]
        if vs:
            solver.Add(solver.Sum([v * c for v, c in vs]) <= (edges[b + 1] - edges[b]) + dmax)
    chosen = {j: solver.Sum([v for (jj, _), v in x.items() if jj == j]) for j in range(len(rows))
              if any(jj == j for (jj, _) in x)}

    # coverage terms
    if objective == "sky_prob":
        ptr, idx, w, val = fs.pix_ptr, fs.pix_idx, fs.pix_w, fs.uprob
    else:
        ptr, idx, w = fs.gal_ptr, fs.gal_idx, fs.gal_w
        val_full = sched.S / fs.host_norm if fs.host_norm > 0 else sched.S * 0
        val = val_full
    contrib: dict[int, list[tuple[Any, float]]] = {}
    for j, i in enumerate(rows):
        if j not in chosen:
            continue
        for k in range(ptr[i], ptr[i + 1]):
            contrib.setdefault(int(idx[k]), []).append((chosen[j], float(w[k])))
    obj_terms = []
    for k, terms in contrib.items():
        if val[k] <= 0:
            continue
        y = solver.NumVar(0.0, 1.0, f"y_{k}")
        solver.Add(y <= solver.Sum([v * ww for v, ww in terms]))
        obj_terms.append(y * float(val[k]))
    solver.Maximize(solver.Sum(obj_terms))
    status = solver.Solve()
    names = {pywraplp.Solver.OPTIMAL: "optimal", pywraplp.Solver.FEASIBLE: "time_limit"}
    if status not in names:
        return {"status": "failed", "reason": f"solver status {status}"}
    bound = solver.Objective().BestBound() if status != pywraplp.Solver.OPTIMAL else solver.Objective().Value()
    return {"status": names[status], "upper_bound": float(bound), "n_fields": int(len(rows)),
            "n_buckets": int(len(edges) - 1), "n_assignment_vars": len(x), "bucket_min": bucket_min,
            "objective": objective, "solver": solver.SolverVersion(), "wall_s": solver.wall_time() / 1000.0}


def greedy_vs_bound(fs: FieldSet, vis, cfg, *, t_start, t_end, t_merger, galaxies_S, objective: str = "sky_prob",
                    **kw: Any) -> dict[str, Any]:
    """Run greedy-v1 with ``objective`` and compare with the MILP upper bound."""
    cfg2 = cfg.model_copy(update={"objective": objective})
    g = GreedyScheduler(fs, vis, cfg2, t_start=t_start, t_end=t_end, t_merger=t_merger, galaxies_S=galaxies_S)
    res = g.run()
    key = "cum_p_sky" if objective == "sky_prob" else f"cum_{fs.host_metric}"
    achieved = float(res.summary[key])
    ub = upper_bound(GreedyScheduler(fs, vis, cfg2, t_start=t_start, t_end=t_end, t_merger=t_merger,
                                     galaxies_S=galaxies_S), **kw)
    out = {"greedy": achieved, **ub}
    if ub.get("upper_bound"):
        out["ratio_lower_bound"] = achieved / ub["upper_bound"] if ub["upper_bound"] > 0 else None
    return out

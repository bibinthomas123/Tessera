"""Stage B: deterministic greedy schedule selection (``greedy-v1``, spec 6.2).

Only visits that fit entirely inside a feasible window (computed from hard
constraints) are ever candidates.  Marginal gains are measured on the union of
coverage so far, so overlapping fields are never double-counted.

Loop:
  1. candidates = visits whose start (now + modelled gap) fits a window entirely
  2. none -> jump to the next window opening / revisit release; stop at horizon end
  3. score = gain / (gap + duration) (``gain_rate``) or gain
  4. among candidates within ``urgency_tolerance`` of the best score, pick the
     one whose window closes first (deterministic tie-break, not a penalty)
  5. remaining ties: higher gain, then lower field ID
  6. commit; advance time, pointing and filter; release dependent visits
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from astropy.time import Time

from tessera.core.config import ResolvedConfig
from tessera.core.detectability import R_FILTER, LightCurve, efficiency
from tessera.core.errors import Diagnostics
from tessera.core.fields import FieldSet
from tessera.core.overheads import OverheadModel
from tessera.core.timeutil import iso
from tessera.core.visibility import Visibility

DAY = 86400.0


@dataclass
class Coverage:
    """Coverage state shared across visits (and, in sequential mode, profiles)."""

    sky: np.ndarray        # covered weight per raster pixel of this field set
    miss: np.ndarray       # probability the counterpart is still undetected per raster pixel
    gal: np.ndarray        # covered weight per galaxy (index into the galaxy table)


@dataclass
class ScheduleResult:
    visits: list[dict[str, Any]]
    unscheduled: dict[str, str]          # field_id -> reason
    summary: dict[str, Any]
    coverage: Coverage
    empty_reason: str | None = None
    diagnostics: list[dict[str, Any]] = field(default_factory=list)


class GreedyScheduler:
    def __init__(self, fs: FieldSet, vis: Visibility, cfg: ResolvedConfig, *, t_start: Time, t_end: Time,
                 t_merger: Time | None, galaxies_S: np.ndarray, coverage: Coverage | None = None,
                 diagnostics: Diagnostics | None = None):
        self.fs, self.vis, self.cfg = fs, vis, cfg
        self.prof = fs.profile
        self.sc = cfg.scheduler
        self.ovh = OverheadModel(self.prof)
        self.lc = LightCurve.named(cfg.light_curve)
        self.t0, self.t1 = t_start.jd, t_end.jd
        self.t_merger = t_merger.jd if t_merger is not None else None
        self.S = np.asarray(galaxies_S, float)
        self.diags = diagnostics if diagnostics is not None else Diagnostics()
        self.objective = cfg.objective
        if self.objective == "yield" and (not np.isfinite(fs.u_dm_mean).any() or self.t_merger is None):
            self.diags.add("YIELD_UNAVAILABLE",
                           "expected-yield objective needs distances and a merger time; using sky_prob")
            self.objective = "sky_prob"
        if self.objective == "host_prob" and fs.gal_idx.size == 0:
            self.diags.add("HOST_OBJECTIVE_NO_GALAXIES", "no cataloged galaxies in any field; using sky_prob")
            self.objective = "sky_prob"
        n_up = len(fs.upix)
        self.cov = coverage or Coverage(np.zeros(n_up), np.ones(n_up), np.zeros(len(self.S)))
        # flattened windows for vectorised lookups
        f_idx, a, b = [], [], []
        for i, ws in enumerate(vis.windows):
            for (wa, wb) in ws:
                f_idx.append(i)
                a.append(wa)
                b.append(wb)
        self.w_f = np.asarray(f_idx, np.int64)
        self.w_a = np.asarray(a, float)
        self.w_b = np.asarray(b, float)
        self._key_scale = 1e5
        self.w_key_a = self.w_f * self._key_scale + (self.w_a - self.t0)
        self.w_key_b = self.w_f * self._key_scale + (self.w_b - self.t0)
        self.row_of_pix = np.repeat(np.arange(len(fs)), np.diff(fs.pix_ptr))
        self.row_of_gal = np.repeat(np.arange(len(fs)), np.diff(fs.gal_ptr))

    # -- vectorised helpers -------------------------------------------------
    def _window_for(self, rows: np.ndarray, start: np.ndarray, end: np.ndarray) -> np.ndarray:
        """Index into flattened windows containing [start, end] for each row, or -1."""
        if len(self.w_f) == 0:
            return np.full(len(rows), -1)
        key = rows * self._key_scale + (start - self.t0)
        pos = np.searchsorted(self.w_key_a, key + 1e-9, side="right") - 1
        posc = np.clip(pos, 0, len(self.w_f) - 1)
        ok = (pos >= 0) & (self.w_f[posc] == rows) & (self.w_b[posc] + 1e-9 >= end)
        return np.where(ok, posc, -1)

    def _next_open(self, rows: np.ndarray, t: float) -> np.ndarray:
        if len(self.w_f) == 0:
            return np.full(len(rows), np.inf)
        key = rows * self._key_scale + (t - self.t0)
        pos = np.searchsorted(self.w_key_b, key, side="right")
        posc = np.clip(pos, 0, len(self.w_f) - 1)
        ok = (pos < len(self.w_f)) & (self.w_f[posc] == rows)
        return np.where(ok, np.maximum(self.w_a[posc], t), np.inf)

    def _airmass(self, rows: np.ndarray, jd: np.ndarray) -> np.ndarray:
        k = np.clip(np.round((jd - self.vis.jd[0]) * DAY / self.vis.step_s).astype(int), 0, len(self.vis.jd) - 1)
        alt = self.vis.alt[rows, k].astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            X = np.where(np.isfinite(alt) & (alt > 1), 1.0 / np.sin(np.radians(alt)), 10.0)
        return X

    def _moon(self, jd: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        k = np.clip(np.round((jd - self.vis.jd[0]) * DAY / self.vis.step_s).astype(int), 0, len(self.vis.jd) - 1)
        return self.vis.moon_illum[k], self.vis.moon_alt[k] > 0

    def _exposures(self, rows: np.ndarray, start: np.ndarray, filt: str) -> tuple[np.ndarray, np.ndarray]:
        spec = self.prof.instrument.filters[filt]
        exp = self.prof.instrument.exposure
        X = self._airmass(rows, start)
        illum, up = self._moon(start)
        if exp.fixed_s is not None:
            total = np.full(len(rows), exp.fixed_s * (exp.n_exp_fixed or 1))
        else:
            if self.t_merger is not None:
                M = self.lc.abs_mag((start - self.t_merger), filt)
            else:
                M = np.full(len(rows), np.inf)
            A = R_FILTER.get(filt, 2.5) * self.fs.ebv_mean[rows]
            need = M + self.fs.dm90[rows] + A
            loss = spec.k_ext * (X - 1.0) + spec.moon.dm_full * illum * up
            with np.errstate(over="ignore", invalid="ignore"):
                total = spec.t_ref_s * 10 ** ((need - spec.m_ref + loss) / 1.25)
            total = np.where(np.isfinite(total), total, exp.t_min_s)
            total = np.clip(total, exp.t_min_s, exp.t_max_s)
        n = np.full(len(rows), exp.n_exp_fixed) if exp.n_exp_fixed else np.maximum(
            1, np.ceil(total / exp.t_single_max_s - 1e-9)).astype(int)
        return n.astype(int), np.round(total / n, 1)

    def _m_lim(self, rows, filt, n_exp, t_single, start, end) -> np.ndarray:
        spec = self.prof.instrument.filters[filt]
        X = np.maximum(self._airmass(rows, start), self._airmass(rows, end))
        illum, up = self._moon(start)
        return (spec.m_ref + 1.25 * np.log10(n_exp * t_single / spec.t_ref_s) - spec.k_ext * (X - 1.0)
                - spec.moon.dm_full * illum * up)

    def _gains(self, rows: np.ndarray, filt: str, mlim: np.ndarray, start: np.ndarray) -> np.ndarray:
        fs = self.fs
        sel = np.zeros(len(fs), bool)
        sel[rows] = True
        out = np.zeros(len(fs))
        if self.objective == "sky_prob":
            m = sel[self.row_of_pix]
            p = fs.pix_idx[m]
            g = fs.uprob[p] * np.clip(fs.pix_w[m] - self.cov.sky[p], 0, None)
            np.add.at(out, self.row_of_pix[m], g)
        elif self.objective == "host_prob":
            m = sel[self.row_of_gal]
            gi = fs.gal_idx[m]
            g = self.S[gi] / fs.host_norm * np.clip(fs.gal_w[m] - self.cov.gal[gi], 0, None)
            np.add.at(out, self.row_of_gal[m], g)
        else:  # yield
            m = sel[self.row_of_pix]
            p = fs.pix_idx[m]
            r = self.row_of_pix[m]
            ml = np.zeros(len(fs))
            ml[rows] = mlim
            st = np.zeros(len(fs))
            st[rows] = start
            M = self.lc.abs_mag(st[r] - self.t_merger, filt)
            A = R_FILTER.get(filt, 2.5) * fs.u_ebv[p]
            eff = efficiency(ml[r] - M - A, fs.u_dm_mean[p], fs.u_dm_sd[p], self.lc.sigma_mag)
            pdet = fs.pix_w[m] * np.nan_to_num(eff)
            np.add.at(out, r, fs.uprob[p] * self.cov.miss[p] * pdet)
        return out[rows]

    def _commit(self, row: int, filt: str, mlim: float, start: float) -> None:
        fs = self.fs
        sl = slice(fs.pix_ptr[row], fs.pix_ptr[row + 1])
        p = fs.pix_idx[sl]
        np.maximum.at(self.cov.sky, p, fs.pix_w[sl])
        if self.t_merger is not None and np.isfinite(fs.u_dm_mean[p]).any():
            M = self.lc.abs_mag(np.full(len(p), start - self.t_merger), filt)
            A = R_FILTER.get(filt, 2.5) * fs.u_ebv[p]
            eff = np.nan_to_num(efficiency(mlim - M - A, fs.u_dm_mean[p], fs.u_dm_sd[p], self.lc.sigma_mag))
            np.multiply.at(self.cov.miss, p, 1.0 - fs.pix_w[sl] * eff)
        gs = slice(fs.gal_ptr[row], fs.gal_ptr[row + 1])
        np.maximum.at(self.cov.gal, fs.gal_idx[gs], fs.gal_w[gs])

    def metrics(self) -> dict[str, float]:
        fs = self.fs
        p_sky = float(np.sum(fs.uprob * self.cov.sky))
        p_host = float(np.sum(self.S * self.cov.gal) / fs.host_norm) if fs.host_norm > 0 else 0.0
        exp_det = float(np.sum(fs.uprob * (1.0 - self.cov.miss)))
        return {"p_sky": p_sky, fs.host_metric: p_host, "expected_detections": exp_det}

    # -- main loop ------------------------------------------------------------
    def run(self) -> ScheduleResult:
        fs, vis, sc = self.fs, self.vis, self.sc
        pattern = self.prof.visit_pattern
        nf = len(fs)
        active = np.zeros(nf, bool)
        active[fs.active()] = True
        unscheduled: dict[str, str] = {}
        for i in np.flatnonzero(fs.status == "rejected"):
            unscheduled[fs.field_id[i]] = fs.reason[i]
        for i in np.flatnonzero(active):
            if not vis.windows[i]:
                r = vis.never_reason[i]
                code = "REJECT_DEC_LIMIT" if r == "declination" else "REJECT_NEVER_VISIBLE"
                unscheduled[fs.field_id[i]] = code if r is None else f"{code}:{r}"
                active[i] = False
            else:
                longest = max(b - a for a, b in vis.windows[i]) * DAY
                n, ts = self._exposures(np.array([i]), np.array([vis.windows[i][0][0]]), pattern[0].filter)
                need = self.ovh.visit_duration(int(n[0]), float(ts[0]))
                if longest < need:
                    unscheduled[fs.field_id[i]] = "REJECT_WINDOW_TOO_SHORT"
                    active[i] = False
        step = np.zeros(nf, int)
        last_end = np.full(nf, -np.inf)
        pattern_gain = np.zeros(nf)
        broken = np.zeros(nf, bool)
        first_visit: dict[int, tuple[str, float, float]] = {}
        block_filter: str | None = None
        block_start = 0.0
        n_filter_changes = 0
        visits: list[dict[str, Any]] = []
        t = max(self.t0, vis.jd[0])
        lst = vis.lst_at(t)
        cur_ra, cur_dec = self.ovh.park_radec(lst)
        cur_filter: str | None = None
        ids = fs.field_id.astype(str)
        id_rank = np.argsort(np.argsort(ids, kind="stable"), kind="stable")
        exposure_s = overhead_s = 0.0

        while t < self.t1 and len(visits) < sc.max_visits:
            elig = active & ~broken & (step < len(pattern))
            if not elig.any():
                break
            rows_all = np.flatnonzero(elig)
            # each eligible field's next step decides its filter; evaluate per filter group
            cands: list[tuple] = []
            lst = vis.lst_at(t)
            for fi, stepdef in enumerate(pattern):
                rows = rows_all[step[rows_all] == fi]
                if not len(rows):
                    continue
                filt = stepdef.filter
                gap = self.ovh.gap(cur_ra, cur_dec, fs.ra[rows], fs.dec[rows], lst,
                                   cur_filter is not None and filt != cur_filter) / DAY
                start = t + gap
                if fi > 0:
                    earliest = last_end[rows] + stepdef.min_gap_min / 1440.0
                    latest = (last_end[rows] + stepdef.max_gap_min / 1440.0) if stepdef.max_gap_min else np.inf
                    ready = earliest <= start + 1e-9
                    expired = start > latest
                    for r in rows[expired]:
                        broken[r] = True
                    keep = ready & ~expired
                    rows, gap, start = rows[keep], gap[keep], start[keep]
                    if not len(rows):
                        continue
                n_exp, t_single = self._exposures(rows, start, filt)
                dur = (n_exp * (t_single + self.prof.instrument.readout_s)) / DAY
                end = start + dur
                win = self._window_for(rows, start, end)
                ok = (win >= 0) & (end <= self.t1)
                if not ok.any():
                    continue
                rows, gap, start, end, n_exp, t_single, dur, win = (x[ok] for x in
                                                                     (rows, gap, start, end, n_exp, t_single, dur, win))
                mlim = self._m_lim(rows, filt, n_exp, t_single, start, end)
                if fi == 0:
                    gain = self._gains(rows, filt, mlim, start)
                else:
                    gain = pattern_gain[rows]
                score = gain / ((gap + dur) * 24.0) if sc.gain_mode == "gain_rate" else gain
                for k in range(len(rows)):
                    cands.append((rows[k], fi, filt, gap[k], start[k], end[k], n_exp[k], t_single[k], win[k],
                                  mlim[k], gain[k], score[k]))
            if not cands:
                # jump to the next time something can start
                rows = rows_all[~broken[rows_all]]
                nxt = self._next_open(rows, t + 60.0 / DAY)
                dep = rows[step[rows] > 0]
                if len(dep):
                    rel = np.array([last_end[r] + pattern[step[r]].min_gap_min / 1440.0 for r in dep])
                    nxt = np.concatenate([nxt, np.maximum(rel, t + 60.0 / DAY)])
                tn = float(np.min(nxt)) if len(nxt) else np.inf
                if not np.isfinite(tn) or tn >= self.t1:
                    break
                t = max(tn, t + 60.0 / DAY)
                continue
            # Filter blocking (optional): within an open block, only visits in the block's
            # filter compete, as long as any such visit has positive gain.  This trades a
            # little gain for fewer filter changes; revisits whose max_gap lapses meanwhile
            # are reported as PATTERN_BROKEN_MAX_GAP, never scheduled late.
            blocked = False
            if (sc.filter_block_min and block_filter is not None
                    and (t - block_start) * 1440.0 < sc.filter_block_min):
                same = [c for c in cands if c[2] == block_filter and c[10] > 0]
                if same and len(same) < len(cands):
                    cands = same
                    blocked = True
            C = np.array([(c[10], c[11]) for c in cands])
            gains, scores = C[:, 0], C[:, 1]
            if np.max(gains) <= 0:
                break
            best = float(np.max(scores))
            tol = (scores >= (1.0 - sc.urgency_tolerance) * best) & (gains > 0)
            wclose = np.array([self.w_b[c[8]] for c in cands])
            ranks = np.array([id_rank[c[0]] for c in cands])
            stepi = np.array([c[1] for c in cands])
            idx_tol = np.flatnonzero(tol)
            # dependent (revisit) steps first, then urgency, then gain, then ID
            order = idx_tol[np.lexsort((ranks[idx_tol], -gains[idx_tol], wclose[idx_tol], -stepi[idx_tol]))]
            pick = int(order[0])
            argmax = int(np.flatnonzero(scores == best)[np.argmin(ranks[scores == best])])
            row, fi, filt, gap, start, end, n_exp, t_single, win, mlim, gain, score = cands[pick]
            if stepi[pick] > 0:
                rule = "revisit_pattern"
            elif pick != argmax:
                rule = f"urgency: window closes in {(wclose[pick] - start) * 1440:.0f} min"
            else:
                rule = "max_gain_rate" if sc.gain_mode == "gain_rate" else "max_gain"
            if blocked:
                rule = f"filter_block:{block_filter}; {rule}"
            if filt != block_filter:
                if block_filter is not None:
                    n_filter_changes += 1
                block_filter, block_start = filt, start
            runners = [int(i) for i in np.argsort(-scores, kind="stable") if i != pick][: sc.runners_up]
            if fi == 0:
                pattern_gain[row] = gain
                first_visit[row] = (filt, float(mlim), float(start))
            # Coverage is credited at the first visit, or only once the whole revisit
            # pattern completes when ``require_complete_pattern`` is set.
            credit_now = (fi == 0) if not sc.require_complete_pattern else (fi == len(pattern) - 1)
            if credit_now:
                f0, m0, s0 = first_visit[row]
                self._commit(int(row), f0, m0, s0)
            step[row] += 1
            last_end[row] = end
            exposure_s += float(n_exp * t_single)
            overhead_s += float(gap * DAY + n_exp * self.prof.instrument.readout_s)
            visits.append({
                "seq": len(visits) + 1,
                "field_id": str(ids[row]),
                "origin": str(fs.origin[row]),
                "ra": round(float(fs.ra[row]), 6),
                "dec": round(float(fs.dec[row]), 6),
                "pa_deg": float(fs.pa[row]),
                "filter": filt,
                "pattern_step": int(fi),
                "n_exp": int(n_exp),
                "t_exp_s": float(t_single),
                "start": iso(start),
                "end": iso(end),
                "start_jd": float(start),
                "end_jd": float(end),
                "gap_before_s": round(float(gap * DAY), 1),
                "window": [iso(self.w_a[win]), iso(self.w_b[win])],
                "m_lim_pred": round(float(mlim), 2),
                "gain": float(gain),
                "gain_rate_per_h": float(gain / ((gap + (end - start)) * 24.0)),
                "decision": {"rule": rule, "n_competitors": len(cands),
                             "runners_up": [{"field_id": str(ids[cands[i][0]]), "score": float(scores[i])}
                                            for i in runners]},
                "cum": self.metrics(),
            })
            t = end
            cur_ra, cur_dec, cur_filter = float(fs.ra[row]), float(fs.dec[row]), filt

        done = set(v["field_id"] for v in visits)
        for i in np.flatnonzero(active):
            fid = str(ids[i])
            if fid in done:
                if step[i] < len(pattern) and broken[i]:
                    unscheduled[fid + "#revisit"] = "PATTERN_BROKEN_MAX_GAP"
                elif step[i] < len(pattern):
                    unscheduled[fid + "#revisit"] = "PATTERN_INCOMPLETE_HORIZON"
                continue
            unscheduled.setdefault(fid, "NOT_SELECTED_HORIZON_FULL" if len(visits) >= sc.max_visits
                                   else "LOWER_GAIN_THAN_SCHEDULED")
        summary = {
            "objective": self.objective,
            "scheduler": sc.algorithm,
            "n_visits": len(visits),
            "n_fields_observed": len(done),
            "t_exposure_s": round(exposure_s, 1),
            "t_overhead_s": round(overhead_s, 1),
            "open_shutter_frac": round(exposure_s / (exposure_s + overhead_s), 4) if visits else 0.0,
            "n_filter_changes": n_filter_changes,
            **{f"cum_{k}": v for k, v in self.metrics().items()},
        }
        hours = (exposure_s + overhead_s) / 3600.0
        summary["expected_detections_per_hour"] = summary["cum_expected_detections"] / hours if hours else 0.0
        empty = None
        if not visits:
            hist: dict[str, int] = {}
            for r in unscheduled.values():
                hist[r.split(":")[0]] = hist.get(r.split(":")[0], 0) + 1
            if nf == 0:
                empty = "NO_CANDIDATE_FIELDS"
            elif not active.any():
                empty = "NO_FIELD_VISIBLE_IN_HORIZON"
            else:
                empty = "NO_VISIT_FITS_A_WINDOW"
            summary["rejection_histogram"] = hist
        return ScheduleResult(visits=visits, unscheduled=unscheduled, summary=summary, coverage=self.cov,
                              empty_reason=empty)

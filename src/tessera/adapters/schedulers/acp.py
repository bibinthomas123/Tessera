"""ACP Observatory Control Software plan files (small robotic telescopes).

Emits one target block per visit::

    ; tessera ...
    #waituntil 1, 08/17/2017 23:13:50
    #count 1
    #interval 30
    #filter r
    #binning 1
    g_GLADE+_1151336    13:09:47.71    -23:23:01.8
    ...
    #quitat 08/18/2017 10:00:00

Times are UTC.  ``#waituntil`` holds each visit to its validated start so the
observatory cannot drift ahead of the feasibility checks.

UNVERIFIED: the ACP help pages (acp.dc3.com) could not be reached when this
adapter was written, so the directive syntax and date format here are from
memory.  Every export carries a warning until a site confirms it against its
ACP version (set ``adapter_options.acp.syntax_confirmed: true``).
"""

from __future__ import annotations

import re
from typing import Any

import astropy.units as u
from astropy.coordinates import Angle

from tessera.adapters.schedulers.base import BaseAdapter, Export, Issue
from tessera.core.profile import TelescopeProfile


def _acp_time(iso: str) -> str:
    d, t = iso.rstrip("Z").split("T")
    y, m, dd = d.split("-")
    return f"{m}/{dd}/{y} {t[:8]}"


def _name(fid: str) -> str:
    return re.sub(r"[^A-Za-z0-9_+.-]", "_", fid)[:48]


class ACPAdapter(BaseAdapter):
    name = "acp"
    version = "1"
    media_type = "text/plain"

    def export(self, schedule: dict[str, Any], profile: TelescopeProfile, *, allow_infeasible: bool = False,
               **opts: Any) -> Export:
        self.check_exportable(schedule, allow_infeasible=allow_infeasible)
        o = self.options(profile, **opts)
        fmap = o.get("filter_map", {})
        binning = int(o.get("binning", 1))
        lines = [f"; {self.header_comment(schedule)}", "; times are UTC"]
        if o.get("syntax_confirmed"):
            lines.append("; syntax_confirmed by site configuration")
        lines.append("")
        for v in schedule["visits"]:
            ra = Angle(v["ra"], u.deg).to_string(unit=u.hour, sep=":", precision=2, pad=True)
            dec = Angle(v["dec"], u.deg).to_string(unit=u.deg, sep=":", precision=1, alwayssign=True, pad=True)
            lines += [
                f"; seq {v['seq']} gain {v['gain']:.4g} rule {v.get('decision', {}).get('rule')}",
                f"#waituntil 1, {_acp_time(v['start'])}",
                f"#count {v['n_exp']}",
                f"#interval {v['t_exp_s']:g}",
                f"#filter {fmap.get(v['filter'], v['filter'])}",
                f"#binning {binning}",
                f"{_name(v['field_id'])}\t{ra}\t{dec}",
                "",
            ]
        lines.append(f"#quitat {_acp_time(schedule['horizon']['end'])}")
        return Export(self.name, self.version, {"tessera_plan.txt": ("\n".join(lines) + "\n").encode()}, "text/plain")

    def validate(self, export: Export) -> list[Issue]:
        text = export.files["tessera_plan.txt"].decode()
        issues = []
        targets = [ln for ln in text.splitlines() if ln and not ln.startswith(("#", ";"))]
        for ln in targets:
            parts = ln.split("\t")
            if len(parts) != 3:
                issues.append(Issue("error", f"target line must be name<TAB>RA<TAB>Dec: {ln!r}"))
        if text.count("#waituntil") != len(targets):
            issues.append(Issue("error", "every target needs a #waituntil"))
        if "; syntax_confirmed" not in text:
            issues.append(Issue("warning", "ACP directive syntax not yet confirmed against ACP documentation "
                                           "for this site (adapter_options.acp.syntax_confirmed)"))
        return issues

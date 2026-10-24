/* tessera minimal UI: plain JS, no external dependencies (works offline / air-gapped). */
"use strict";
const API = location.pathname.replace(/\/ui\/?.*$/, "");
const $ = (id) => document.getElementById(id);
const keyInput = $("apikey");
try { keyInput.value = localStorage.getItem("tessera_key") || ""; } catch (e) { /* storage unavailable */ }
keyInput.addEventListener("change", () => { try { localStorage.setItem("tessera_key", keyInput.value); } catch (e) { } });

function headers(extra) {
  const h = Object.assign({}, extra || {});
  if (keyInput.value) h["Authorization"] = "Bearer " + keyInput.value;
  return h;
}
async function api(path, opts) {
  opts = opts || {};
  const r = await fetch(API + path, Object.assign({}, opts, { headers: headers(opts.headers) }));
  const ct = r.headers.get("content-type") || "";
  const body = ct.includes("json") ? await r.json() : await r.blob();
  if (!r.ok) throw new Error((body && (body.detail || body.title)) || r.statusText);
  return body;
}
function esc(s) { return String(s == null ? "" : s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c])); }
function fmt(v, d) {
  if (v == null) return "–";
  if (typeof v !== "number") return esc(v);
  if (Number.isInteger(v) && d == null) return String(v);
  if (v !== 0 && (Math.abs(v) < 1e-3 || Math.abs(v) >= 1e5)) return v.toExponential(2);
  return v.toFixed(d == null ? 3 : d);
}
function setStatus(t) { $("status").textContent = t; }
// Progress of a running plan is shown inside the Plan button; messages and errors go below it.
function setBusy(t) {
  const b = $("planbtn");
  b.disabled = true; b.classList.add("busy"); b.setAttribute("aria-busy", "true");
  $("planlabel").textContent = t;
}
function setIdle() {
  const b = $("planbtn");
  b.disabled = false; b.classList.remove("busy"); b.removeAttribute("aria-busy");
  $("planlabel").textContent = "Plan";
}

async function loadRefs() {
  try {
    const [profiles, catalogs, sci] = await Promise.all([api("/v1/profiles"), api("/v1/catalogs"), api("/v1/science")]);
    $("profiles").innerHTML = profiles.profiles.map((p) =>
      `<label title="${esc(p.description || "")}"><input type="checkbox" value="${esc(p.id)}"> ${esc(p.key)} <small>(${esc(p.validation ? p.validation.level : "?")})</small></label>`).join("");
    $("catalog").innerHTML = catalogs.catalogs.map((c) => `<option value="${esc(c.name)}">${esc(c.name)} — ${esc(c.row_count)} rows</option>`).join("");
    $("strategy").innerHTML = '<option value="">(preset default)</option>' + Object.keys(sci.strategies).map((s) => `<option>${esc(s)}</option>`).join("");
    loadHistory();
  } catch (e) { setStatus("Enter an API key to load profiles and catalogs (" + e.message + ")"); }
}
keyInput.addEventListener("change", loadRefs);

$("profileupload").addEventListener("change", async (ev) => {
  const f = ev.target.files[0]; if (!f) return;
  const text = await f.text();
  let doc;
  try { doc = JSON.parse(text); } catch (e) { setStatus("Upload a JSON profile here; YAML profiles go through `tessera profiles --check` and POST /v1/profiles"); return; }
  try { const r = await api("/v1/profiles", { method: "POST", body: JSON.stringify(doc), headers: { "content-type": "application/json" } }); setStatus("profile published: " + r.key); loadRefs(); }
  catch (e) { setStatus("profile rejected: " + e.message); }
});

$("planform").addEventListener("submit", async (ev) => {
  ev.preventDefault();
  try {
    let skymapId = null;
    const fd = new FormData();
    if ($("skymap").files[0]) fd.append("file", $("skymap").files[0]);
    else if ($("skymapurl").value) fd.append("url", $("skymapurl").value);
    else { setStatus("choose a sky map file or URL"); return; }
    setStatus("");
    setBusy("Uploading…");
    skymapId = (await api("/v1/skymaps", { method: "POST", body: fd })).skymap_id;
    const profiles = [...document.querySelectorAll("#profiles input:checked")].map((i) => i.value);
    const config = { preset: $("preset").value };
    if ($("strategy").value) config.strategy = $("strategy").value;
    if ($("nights").value) config.horizon_nights = parseFloat($("nights").value);
    const body = { skymap_id: skymapId, catalog: $("catalog").value || null, profiles, config, t_start: $("tstart").value || null };
    setBusy("Queued…");
    let job = await api("/v1/plans", { method: "POST", body: JSON.stringify(body), headers: { "content-type": "application/json" } });
    while (job.state === "queued" || job.state === "running") {
      setBusy(job.state === "queued" ? "Queued…" : `Running${job.stage ? " · " + job.stage : "…"}`);
      await new Promise((r) => setTimeout(r, 1500));
      job = await api("/v1/jobs/" + job.job_id);
    }
    setIdle();
    if (job.state !== "succeeded") { setStatus(`Plan ${job.state}: ${job.error ? job.error.code + " " + (job.error.message || "") : ""}`); return; }
    showPlan(job.plan_key);
  } catch (e) { setIdle(); setStatus("Error: " + e.message); }
});

let current = null;
async function showPlan(key) {
  current = key;
  const p = await api("/v1/plans/" + encodeURIComponent(key));
  $("result").classList.remove("hidden");
  $("plankey").textContent = key;
  const b = p.budget;
  $("budget").innerHTML = [["Cataloged (P_cat)", b.p_cat], ["Uncataloged (P_uncat)", b.p_uncat], ["Outside footprint", b.p_out_footprint],
  ["Not assessed", b.p_unassessed], ["Beyond d where C<0.5", b.p_beyond_dmax], ["Weight model", b.weight]]
    .map(([k, v]) => `<tr><th>${k}</th><td>${fmt(v)}</td></tr>`).join("");
  $("comparison").innerHTML = Object.entries(p.strategy_comparison).map(([prof, c]) =>
    `<p><strong>${esc(prof)}</strong> — recommended: <strong>${esc(c.recommendation.strategy)}</strong> <small>${esc(c.recommendation.reasons.join("; "))}</small></p>
     <table><tr><th>strategy</th><th>visits</th><th>P_sky</th><th>host (cat.)</th><th>E[det]</th><th>E[det]/h</th><th>hours</th></tr>` +
    Object.entries(c.strategies).map(([s, v]) => `<tr><td>${esc(s)}${s === c.primary ? " ★" : ""}</td><td>${v.n_visits}</td><td>${fmt(v.cum_p_sky)}</td>
      <td>${fmt(v.cum_score_fraction_in_field != null ? v.cum_score_fraction_in_field : v.cum_p_host_cat)}</td><td>${fmt(v.cum_expected_detections)}</td>
      <td>${fmt(v.expected_detections_per_hour)}</td><td>${fmt(v.telescope_time_h, 2)}</td></tr>`).join("") + "</table>").join("");
  $("diagnostics").innerHTML = p.diagnostics.map((d) => `<li class="lvl-${d.level}"><code>${esc(d.code)}</code> ${esc(d.message)}</li>`).join("") || "<li>none</li>";
  const gal = await api(`/v1/plans/${encodeURIComponent(key)}/galaxies?format=json`);
  const cols = ["rank", "tier", "gal_id", "ra", "dec", "dist_mpc", "credible_level", "p_sky_deg2", "p_dist_mpc3", "w", "score_fraction", "flags"];
  $("galaxies").innerHTML = "<tr>" + cols.map((c) => `<th>${c}</th>`).join("") + "</tr>" +
    gal.rows.slice(0, 200).map((r) => "<tr>" + cols.map((c) => `<td>${fmt(r[c], c === "ra" || c === "dec" ? 4 : undefined)}</td>`).join("") + "</tr>").join("");
  $("schedules").innerHTML = p.schedules.map((s) => `<div class="sched"><p><strong>${esc(s.profile)}</strong> <span class="badge ${esc(s.status)}">${esc(s.status)}</span> v${s.current_version}
      · export: ${["generic", "acp", "lco", "treasuremap", "rubin_too"].map((a) => `<a href="#" data-exp="${esc(s.schedule_id)}|${s.current_version}|${a}">${a}</a>`).join(" · ")}
      ${s.status === "proposed" || s.status === "edited" ? ` <button class="secondary" data-approve="${esc(s.schedule_id)}|${s.current_version}">approve v${s.current_version}</button>` : ""}</p>
      <div class="scroll" id="v-${esc(s.schedule_id)}"></div></div>`).join("") || "<p>No telescope profiles were selected.</p>";
  for (const s of p.schedules) {
    const d = await api("/v1/schedules/" + s.schedule_id);
    const v = d.doc.visits;
    document.getElementById("v-" + s.schedule_id).innerHTML = `<p><small>${v.length} visits · P_sky ${fmt(d.doc.summary.cum_p_sky)} · E[det] ${fmt(d.doc.summary.cum_expected_detections)} · open shutter ${fmt(d.doc.summary.open_shutter_frac, 2)}${d.doc.empty_reason ? " · " + esc(d.doc.empty_reason) : ""}</small></p>` +
      "<table><tr><th>#</th><th>field</th><th>filter</th><th>n×t</th><th>start</th><th>end</th><th>X max</th><th>gain</th><th>why</th><th>ok</th></tr>" +
      v.slice(0, 300).map((x) => `<tr><td>${x.seq}</td><td>${esc(x.field_id)}</td><td>${esc(x.filter)}</td><td>${x.n_exp}×${x.t_exp_s}s</td><td>${esc(x.start)}</td><td>${esc(x.end)}</td>
        <td>${fmt(x.airmass && x.airmass.max, 2)}</td><td>${fmt(x.gain)}</td><td>${esc(x.decision.rule)}</td><td>${x.feasible ? "✓" : "✗"}</td></tr>`).join("") + "</table>";
  }
  drawSky(key, gal.rows);
}

async function drawSky(key, galaxies) {
  const reg = await api(`/v1/plans/${encodeURIComponent(key)}/files/region.json`);
  const cv = $("sky"), ctx = cv.getContext("2d");
  const css = getComputedStyle(document.documentElement);
  ctx.clearRect(0, 0, cv.width, cv.height);
  const pts = reg.preview; if (!pts.length) return;
  // centre on the probability-weighted mean direction; equirectangular around it
  let x = 0, y = 0, z = 0;
  for (const [ra, dec, p] of pts) { const a = ra * Math.PI / 180, d = dec * Math.PI / 180; x += p * Math.cos(d) * Math.cos(a); y += p * Math.cos(d) * Math.sin(a); z += p * Math.sin(d); }
  const ra0 = (Math.atan2(y, x) * 180 / Math.PI + 360) % 360, dec0 = Math.atan2(z, Math.hypot(x, y)) * 180 / Math.PI;
  const dra = (ra) => ((ra - ra0 + 540) % 360) - 180;
  const cosd = Math.cos(dec0 * Math.PI / 180);
  // fit everything drawn (region, fields, top galaxies) with equal scale on both axes
  const polys = Object.values(reg.footprints).flat().flatMap((f) => f.polygons);
  const gals = galaxies.slice(0, 50);
  const xs = [], ys = [];
  for (const [ra, dec] of pts) { xs.push(dra(ra) * cosd); ys.push(dec - dec0); }
  for (const poly of polys) for (const [ra, dec] of poly) { xs.push(dra(ra) * cosd); ys.push(dec - dec0); }
  for (const g of gals) { xs.push(dra(g.ra) * cosd); ys.push(g.dec - dec0); }
  const capH = 26, pad = 14, W = cv.width, H = cv.height - capH;
  const xmin = Math.min(...xs), xmax = Math.max(...xs), ymin = Math.min(...ys), ymax = Math.max(...ys);
  const scale = Math.min((W - 2 * pad) / Math.max(xmax - xmin, 0.5), (H - 2 * pad) / Math.max(ymax - ymin, 0.5));
  const xc = (xmin + xmax) / 2, yc = (ymin + ymax) / 2;
  // east is left on the sky: larger RA plots to the left
  const sx = (ra) => W / 2 - (dra(ra) * cosd - xc) * scale;
  const sy = (dec) => H / 2 - (dec - dec0 - yc) * scale;
  // credible region: each preview pixel (HEALPix order 7, ~0.46 deg) drawn at its true size
  const pix = 0.458 * scale, pmax = Math.max(...pts.map((p) => p[2]));
  ctx.fillStyle = css.getPropertyValue("--region");
  for (const [ra, dec, p] of pts) {
    ctx.globalAlpha = 0.12 + 0.6 * Math.sqrt(p / pmax);
    ctx.fillRect(sx(ra) - pix / 2, sy(dec) - pix / 2, pix + 0.5, pix + 0.5);
  }
  ctx.globalAlpha = 1; ctx.strokeStyle = css.getPropertyValue("--field"); ctx.lineWidth = 1.25;
  for (const poly of polys) {
    ctx.beginPath(); poly.forEach(([ra, dec], i) => (i ? ctx.lineTo(sx(ra), sy(dec)) : ctx.moveTo(sx(ra), sy(dec)))); ctx.closePath(); ctx.stroke();
  }
  ctx.fillStyle = css.getPropertyValue("--gal");
  for (const g of gals) { ctx.beginPath(); ctx.arc(sx(g.ra), sy(g.dec), 3, 0, 2 * Math.PI); ctx.fill(); }
  // caption strip below the plot, never over the data
  ctx.fillStyle = css.getPropertyValue("--bg"); ctx.fillRect(0, H, W, capH);
  ctx.strokeStyle = css.getPropertyValue("--line"); ctx.beginPath(); ctx.moveTo(0, H + 0.5); ctx.lineTo(W, H + 0.5); ctx.stroke();
  ctx.fillStyle = css.getPropertyValue("--muted"); ctx.font = "12px system-ui";
  const deg = 1 / scale * 100;
  ctx.fillText(`centre RA ${ra0.toFixed(1)}°, Dec ${dec0.toFixed(1)}° · 100 px = ${deg.toFixed(1)}° · east left`, 8, H + 17);
}

document.addEventListener("click", async (ev) => {
  const t = ev.target;
  try {
    if (t.dataset.dl && current) { ev.preventDefault(); const blob = await api(`/v1/plans/${encodeURIComponent(current)}/${t.dataset.dl}`); saveBlob(blob, t.dataset.dl.replace(/[?=]/g, ".")); }
    if (t.dataset.exp) { ev.preventDefault(); const [sid, v, a] = t.dataset.exp.split("|"); const blob = await api(`/v1/schedules/${sid}/versions/${v}/export?adapter=${a}`); saveBlob(blob, `${sid}-v${v}-${a}`); }
    if (t.dataset.approve) { const [sid, v] = t.dataset.approve.split("|"); if (!confirm(`Approve ${sid} version ${v}? This is recorded with your identity.`)) return; await api(`/v1/schedules/${sid}/versions/${v}/approve`, { method: "POST" }); showPlan(current); }
    if (t.dataset.plan) { ev.preventDefault(); showPlan(t.dataset.plan); }
  } catch (e) { alert(e.message); }
});
function saveBlob(blob, name) {
  if (!(blob instanceof Blob)) blob = new Blob([JSON.stringify(blob, null, 2)], { type: "application/json" });
  const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
async function loadHistory() {
  try {
    const [plans, events] = await Promise.all([api("/v1/plans?limit=20"), api("/v1/events")]);
    $("plans").innerHTML = plans.plans.map((p) => `<li><a href="#" data-plan="${esc(p.plan_key)}">${esc((p.event && p.event.superevent_id) || "local")} · ${esc(p.created.slice(0, 19))}</a></li>`).join("") || "<li>none yet</li>";
    $("events").innerHTML = events.events.map((e) => `<li>${esc(e.superevent_id)} — ${esc(e.state)}</li>`).join("") || "<li>none yet</li>";
  } catch (e) { /* needs read scope */ }
}
loadRefs();

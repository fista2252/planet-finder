"use strict";

const el = id => document.getElementById(id);
const esc = value => String(value ?? "Unavailable").replace(/[&<>"']/g, c =>
  ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"})[c]);
const num = (value, decimals = 2) => Number.isFinite(value) ? value.toFixed(decimals) : "Unavailable";
const words = value => String(value).replaceAll("_", " ");
const safe = url => typeof url === "string" && url.startsWith("https://") ? esc(url) : "#";
let run, selected, queued;

function table(headers, rows) {
  if (!rows.length) return "<p>No entries in the selected range.</p>";
  return '<div class="table-wrap"><table><thead><tr>' + headers.map(x => '<th scope="col">' + esc(x) + '</th>').join("") +
    '</tr></thead><tbody>' + rows.map(row => "<tr>" + row.map(x => "<td>" + esc(x) + "</td>").join("") + "</tr>").join("") + "</tbody></table></div>";
}

function fields(items) {
  return "<dl>" + items.map(([label, value]) => "<div><dt>" + esc(label) + "</dt><dd>" + esc(value) + "</dd></div>").join("") + "</dl>";
}

function links(target) {
  const refs = target.prior_reports ?? [];
  if (!refs.length) return "<p>No exact report in the checked inputs. Priority remains unresolved.</p>";
  return "<ul>" + refs.map(r => '<li><a href="' + safe(r.source) + '">' + esc(words(r.kind)) + " · " +
    esc(r.label) + '</a>' + (r.period_days == null ? " · single event" : " · " + num(r.period_days, 7) + " d") +
    (r.record?.registered ? " · registered " + esc(r.record.registered) : "") + "</li>").join("") + "</ul>";
}

function eventPlot(events) {
  const measured = events.filter(r => r.measurement);
  if (!measured.length) return "";
  const x = measured.map(r => r.midpoint_btjd);
  const y = measured.map(r => 100 * r.measurement.depth);
  const loX = Math.min(...x), spanX = Math.max(1, Math.max(...x) - loX);
  const loY = Math.min(0, ...y), spanY = Math.max(.05, Math.max(0, ...y) - loY);
  const points = measured.map((r, i) =>
    '<circle cx="' + (55 + 615 * (x[i] - loX) / spanX) + '" cy="' + (178 - 108 * (y[i] - loY) / spanY) +
    '" r="4"><title>Sector ' + esc(r.sector) + ' · ' + num(x[i], 5) + ' BTJD · ' + num(y[i], 3) + '%</title></circle>').join("");
  return '<svg viewBox="0 0 700 225" role="img" aria-label="Measured event depth by midpoint">' +
    '<rect width="700" height="225" fill="#f7f8fc"/><g font-family="system-ui" fill="#18313c" font-size="13">' +
    '<text x="22" y="34" font-size="17">Local event depths (%)</text><path d="M55 64v114h615" fill="none" stroke="#8ba3a6"/>' +
    '<text x="16" y="74">' + num(loY + spanY, 2) + '</text><text x="16" y="180">' + num(loY, 2) +
    '</text><text x="55" y="202">' + num(loX, 1) + '</text><text x="670" y="202" text-anchor="end">' +
    num(loX + spanX, 1) + '</text><text x="355" y="218" text-anchor="middle">BTJD</text><g fill="#155c54">' +
    points + '</g></g></svg>';
}

function archived(t) {
  const obs = t.observed;
  const measurements = fields([
    ["Original sector", obs.original_sector], ["Original midpoint (BTJD)", num(obs.original_midpoint_btjd, 5)],
    ["Original box depth", obs.original_depth_fraction === null ? "Unavailable" : num(100 * obs.original_depth_fraction, 3) + "%"],
    ["Later sectors with measurements", obs.measured_later_sectors],
    ["Indexed products", t.n_products], ["Product errors", obs.product_errors]
  ]);
  const methods = t.solutions.length ? table(["Period (d)", "Conditional σ (d)", "Duration (h)", "Cycle count"],
    t.solutions.map(s => [num(s.period_days, 7), num(s.period_sigma_days, 7), num(s.duration_hours, 2), s.alias_cycles])) : "";
  const image = t.solutions.length ? '<img src="' + esc(t.tic) + '.svg" alt="Fitted event depths versus midpoint for TIC ' + esc(t.tic) + '">' : "";
  const evidence = run.details[t.tic];
  const alias = Array.isArray(evidence?.diagnostics) ? evidence.diagnostics : null;
  const coverage = alias ? table(["Period (d)", "Windows within archived spans", "Measurable timing offsets", "Rejected"],
    alias.map(a => [num(a.period, 5), a.predicted_in_product_span.length,
      a.predicted_in_product_span.reduce((n, p) => n + p.offset_measurements.length, 0), a.rejected ? "Yes" : "No"])) : "";
  return '<p>' + esc(t.caveat) + '</p>' + measurements + methods + image + coverage +
    '<details><summary>Evidence scope</summary><p>TESS-only photometry. Alternate reductions share spacecraft images. Fit and timing errors depend on the adopted shape and linear period. Unavailable cadence, catalog data and spectroscopy leave those tests unresolved.</p></details>';
}

function current(t) {
  const e = t.evidence;
  const review = run.decisions?.[t.tic];
  const reviewState = run.review_queue?.find(row => row.tic === t.tic)?.review_state;
  const evidence = fields([
    ["Measurement status", words(e.status)], ["Mode", e.mode ? words(e.mode) : "Unavailable"],
    [e.period_days == null && e.trial_period_days ? "Leading trial period (d)" : "Period (d)", num(e.period_days ?? e.trial_period_days, 7)],
    ["Trial epoch (BTJD)", num(e.epoch_btjd ?? e.trial_epoch_btjd, 6)],
    ["Supported cycles", e.supported_cycles], ["Sensitive non-detections", e.sensitive_missing_cycles]
  ]);
  const flags = t.review_reasons.length ? "<p>Review flags: " + esc(t.review_reasons.map(words).join(", ")) + ".</p>" : "";
  const events = e.events ?? [];
  const rows = events.map(r => [r.sector, r.cycle, num(r.midpoint_btjd, 6), words(r.event_state ?? r.coverage),
    num(r.measurement?.depth * 100, 3), num(r.measurement?.score, 1)]);
  const productRows = t.provenance.map(r => [r.sequence_number, r.provenance_name, r.cadence_count, r.sha256]);
  const catalogRows = t.catalog_checks.map(r => [r.source_name ?? r.kind, words(r.state), r.rows_checked,
    r.revision?.slice(0, 12), r.revision_date ?? r.source_modified_at, r.retrieved_at]);
  const alternatives = (e.alternatives ?? []).map(a => [num(a.period_days, 7), words(a.method),
    a.supported_cycles, a.sensitive_missing_cycles, a.sensitive_missing_cycles >= 2 ? "Contradicted" :
      a.supported_cycles >= 3 ? "Retained for comparison" : "Insufficient support"]);
  const priorWarning = t.prior_warning ? '<p role="note">' + esc(t.prior_warning) + '</p>' : "";
  return evidence + priorWarning + flags + (review ? "<p>" + (reviewState === "stale" ? "Previous review · evidence changed" : "Review decision") +
    " · " + esc(words(review.decision)) + " · " + esc(review.reason) + "</p>" : "") +
    eventPlot(events) + (e.uncertainty ? "<p>" + esc(e.uncertainty) + "</p>" : "") +
    (alternatives.length ? '<details><summary>Tested periods and harmonics</summary>' +
      table(["Period (d)", "Origin", "Supported cycles", "Sensitive misses", "Comparison"], alternatives) + '</details>' : '') +
    '<details><summary>Observed cadence and local measurements</summary>' + table(["Sector", "Cycle", "BTJD", "Window result", "Depth (%)", "Nominal score"], rows) + '</details>' +
    '<details><summary>Products and SHA-256 digests</summary>' + table(["Sector", "Pipeline", "Cadences", "SHA-256"], productRows) + '</details>' +
    '<details><summary>Stage status and catalog checks</summary>' + table(["Stage", "State"], Object.entries(t.stages).map(([k,v]) => [words(k), words(v)])) +
    table(["Catalog", "State", "Rows", "Revision", "Source date", "Retrieved (UTC)"], catalogRows) +
    (t.failures.length ? "<ul>" + t.failures.map(f => "<li>" + esc(f.stage + " · " + f.error) + "</li>").join("") + "</ul>" : "") + '</details>';
}

function selectTarget(t, updateHash = true) {
  selected = t.tic;
  if (updateHash) history.replaceState(null, "", "#tic-" + t.tic);
  for (const button of el("list").querySelectorAll("button")) button.setAttribute("aria-pressed", button.dataset.tic === selected);
  const coord = t.coordinates_icrs_deg ?? [t.ra_deg, t.dec_deg];
  el("detail").innerHTML = '<span class="pill ' + (t.prior_status.startsWith("recovered") ? "recovered" : "") + '">' +
    esc(t.prior_status === "unresolved" ? "Unresolved priority" : t.prior_status) + '</span><h3>TIC ' + esc(t.tic) + '</h3><p>' +
    esc(t.review || "Run evidence") + ' · ICRS ' + num(coord[0], 5) + '°, ' + num(coord[1], 5) + '°</p>' +
    (run.run_id ? current(t) : archived(t)) + '<div><h3>Prior reports</h3>' + links(t) + '</div>' +
    '<p><a href="https://exofop.ipac.caltech.edu/tess/target.php?id=' + esc(t.tic) + '">ExoFOP target ↗</a></p>';
}

function listTargets() {
  const query = el("search").value.trim().replace(/^TIC\s*/i, "");
  const filter = el("filter").value;
  const visible = [...run.targets].sort((a, b) => Number(Boolean(b.solutions?.length)) - Number(Boolean(a.solutions?.length)))
    .filter(t => t.tic.includes(query) && (filter === "all" ||
    (filter === "flagged" && queued.has(t.tic)) || (filter === "recovered" && t.prior_status.startsWith("recovered")) ||
    (filter === "unresolved" && t.prior_status === "unresolved")));
  el("list").replaceChildren();
  for (const t of visible) {
    const button = document.createElement("button");
    button.className = "target";
    button.dataset.tic = t.tic;
    button.innerHTML = "<strong>TIC " + esc(t.tic) + "</strong><small>" +
      esc(t.prior_status === "unresolved" ? t.review ?? words(t.evidence.status) : t.prior_status) + "</small>";
    button.addEventListener("click", () => selectTarget(t));
    el("list").append(button);
  }
  const initial = visible.find(t => t.tic === selected) ?? visible[0];
  if (initial) selectTarget(initial, false);
  else el("detail").textContent = "No targets match these filters.";
}

function planner() {
  if (!run.run_id) {
    el("schedule-download").hidden = false;
    el("planner-note").textContent = "Frozen planning windows. BTJD is BJD TDB minus 2457000. The 95%-style half-widths are conditional. Future detector coverage was projected for TIC 238187210 only.";
    el("windows").innerHTML = table(["Approx. UTC date", "TIC", "Period (d)", "Midpoint (BTJD)", "Half-width (min)", "Coverage"],
      [...run.predictions].sort((a,b) => a.btjd - b.btjd).map(p => [
        p.utc_date_approx, p.tic, num(p.period_days, 6), num(p.btjd, 5), num(p.planning_halfwidth_minutes, 1),
        p.coverage ? "TESS sector " + p.coverage.sector + " off detector. Sutherland " + num(p.coverage.usable_fraction * 100, 0) +
          "% of full window. Midpoint " + p.coverage.utc_mid.slice(11,16) + " UTC." : "Future coverage not evaluated"
      ]));
    return;
  }
  el("planner-note").textContent = "Enter a BTJD range of up to 366 days. These linear projections have no calibrated timing uncertainty or future footprint check.";
  const form = el("plan-form");
  form.hidden = false;
  form.addEventListener("submit", event => {
    event.preventDefault();
    const lo = Number(el("plan-start").value), hi = Number(el("plan-end").value);
    if (lo < 1000 || hi > 10000 || hi <= lo || hi - lo > 366) {
      el("windows").textContent = "Enter an increasing BTJD interval within 1000–10000 and at most 366 days long.";
      return;
    }
    const rows = [];
    for (const t of run.targets) {
      const e = t.evidence;
      if (e.status !== "measured" || !(e.period_days > 0) || !(e.supported_cycles >= 3)) continue;
      for (let i = Math.ceil((lo-e.epoch_btjd)/e.period_days); i < Math.ceil((hi-e.epoch_btjd)/e.period_days); i++) {
        rows.push([t.tic, e.epoch_btjd + i * e.period_days, 24 * e.duration_days]);
      }
    }
    rows.sort((a,b) => a[1] - b[1]);
    el("windows").innerHTML = table(["TIC", "Midpoint (BTJD)", "Duration (h)"],
      rows.map(([tic, mid, duration]) => [tic, num(mid, 6), num(duration, 2)]));
  });
}

async function initialize() {
  try {
    const response = await fetch("data.json");
    if (!response.ok) throw Error("HTTP " + response.status);
    run = await response.json();
    if (run.schema_version !== 1 || !Array.isArray(run.targets)) throw Error("Unsupported run schema");
    const queueResponse = run.run_id ? null : await fetch("review-queue.json");
    const queue = run.run_id ? run.review_queue : queueResponse.ok ? await queueResponse.json() : [];
    queued = new Set(queue.filter(r => r.review_state !== "reviewed").map(r => r.tic));
    el("run-label").textContent = run.run_id ? "Local run · " + run.run_id : "Imported pilot · " + run.snapshot_date;
    el("summary").innerHTML = [
      [run.targets.length, "Targets"],
      [queued.size, "Flagged"],
      [run.targets.filter(t => t.prior_status.startsWith("recovered")).length, "Prior periods recovered"]
    ].map(([n,label]) => "<div><strong>" + n + "</strong><small>" + label + "</small></div>").join("");
    selected = location.hash.startsWith("#tic-") ? location.hash.slice(5) : run.run_id ? run.targets[0]?.tic : "238187210";
    el("search").addEventListener("input", listTargets);
    el("filter").addEventListener("change", listTargets);
    listTargets();
    planner();
    el("message").hidden = true;
    el("targets").hidden = false;
    el("planner").hidden = false;
  } catch (error) {
    el("message").textContent = "Run data could not be loaded. Serve the export with python -m http.server --directory EXPORT. " + error.message;
  }
}

initialize();

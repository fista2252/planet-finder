"""Validate the pilot snapshot, crossmatch prior reports and build the static site."""
import argparse
import csv
import gzip
import html
import io
import json
import math
import shutil
from pathlib import Path

from import_pilot import ROOT, json_bytes, sha

DATA = ROOT / "data"
CACHE = ROOT / ".audit-cache"
DETAILED = {
    "238187210": ("Repeating eclipse", "Pixel deficits favor the target over catalog neighbors bright enough to cause the dip. An unresolved stellar eclipse remains possible. No checked spectrum settles the mass."),
    "269730302": ("Grazing repeat", "The shared fit hits the impact-parameter bound. Gaia RUWE is 2.83. Stellar companions and unresolved blends remain plausible."),
    "339672790": ("Four surviving aliases", "Two measured dips support four periods. The Sector 90 gap leaves three of five timing offsets unmeasurable for the shortest alias."),
}


def load(path):
    return json.loads(path.read_bytes())


def snapshot():
    result = {}
    for row in load(DATA / "snapshot-manifest.json"):
        path = DATA / "snapshot" / row["path"]
        payload = path.read_bytes()
        if sha(payload) != row["sha256"]:
            raise ValueError(f"Snapshot SHA-256 mismatch: {row['path']}")
        raw = gzip.decompress(payload) if path.suffix == ".gz" else payload
        if sha(raw) != row["uncompressed_sha256"]:
            raise ValueError(f"Uncompressed SHA-256 mismatch: {row['path']}")
        if row["logical_path"] in result:
            raise ValueError(f"Duplicate snapshot path: {row['logical_path']}")
        result[row["logical_path"]] = raw
    return result


def csv_rows(raw):
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))


def atomic(path, raw):
    from tempfile import NamedTemporaryFile

    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        temp = Path(stream.name)
        stream.write(raw)
    try:
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


def checkpoint(path, key, make):
    if path.is_file():
        try:
            saved = load(path)
            if saved["key"] == key and saved["value_sha256"] == sha(json_bytes(saved["value"])):
                return saved["value"]
        except (KeyError, ValueError, OSError):
            pass
    value = make()
    atomic(path, json_bytes(dict(key=key, value=value, value_sha256=sha(json_bytes(value)))))
    return value


def prior_reports(tic, periods, catalogs):
    evidence = []
    for catalog in catalogs:
        usable = catalog["state"] == "checked_snapshot" and catalog.get("rows_checked", 0) > 0
        rows = [r for r in catalog.get("matches", []) if r["tic"] == tic] if usable else []
        for row in rows:
            p = row["period_days"]
            # Tolerance screens period agreement only. Harmonics and epoch
            # consistency can require review, even with an exact TIC match.
            row = dict(row, period_agrees=p is not None and any(abs(p - q) <= .001 * q for q in periods))
            evidence.append(dict(kind=catalog["kind"], source=catalog["url"], **row))
    if any(r["period_agrees"] and r["kind"] in {"confirmed", "eclipsing_binary"} for r in evidence):
        status = "known catalog signal"
    elif any(r["period_agrees"] for r in evidence):
        status = "recovered prior period"
    else:
        status = "unresolved"
    return status, evidence


def number(value, lo, hi):
    value = float(value)
    if not math.isfinite(value) or not lo <= value <= hi:
        raise ValueError(f"Out-of-range measurement: {value}")
    return value


def optional_number(value, lo, hi):
    return number(value, lo, hi) if value != "" else None


def validate_targets(cohort, summary, fits_manifest):
    ids = [x["catalog"]["TICID"] for x in cohort]
    if len(ids) != 24 or len(set(ids)) != 24 or any(not x.isdecimal() for x in ids):
        raise ValueError("Expected the 24 unique frozen TIC identifiers")
    if sorted(ids) != sorted(r["tic"] for r in summary):
        raise ValueError("Summary and cohort disagree")
    products = {p["dataURL"] for row in cohort for p in row["primary"]}
    if products != set(fits_manifest):
        raise ValueError("Primary product manifests disagree")
    for item in cohort:
        for product in item["primary"]:
            if not product["dataURL"].startswith("mast:") or product.get("dataRights") != "PUBLIC":
                raise ValueError("Unexpected or non-public product")
            if product.get("target_name") != item["catalog"]["TICID"]:
                raise ValueError("Product and cohort target disagree")


def prepare(raw, catalogs, fits_manifest):
    read = lambda path: json.loads(raw[path])
    cohort = read("cohort.json")
    summary = csv_rows(raw["pilot_summary.csv"])
    validate_targets(cohort, summary, fits_manifest)
    rows = {r["tic"]: r for r in summary}
    models = read("followup/transit_fits.json")
    diag = read("followup/diagnostics.json")
    timing = read("followup/timing_check.json")
    pixel = read("followup/centroid_check.json")
    visibility = read("followup/tic238_visibility.json")
    target_audit = read("followup/tic238_audit.json")
    targets = []
    for item in cohort:
        catalog = item["catalog"]
        tic = catalog["TICID"]
        row = rows[tic]
        product = item["primary"][0]
        observed = dict(original_sector=int(row["original_sector"]),
                        original_midpoint_btjd=optional_number(row["original_btjd"], 0, 10000),
                        original_depth_fraction=optional_number(row["original_depth"], -1, 1),
                        original_duration_days=optional_number(row["original_duration_days"], 0, 2),
                        original_score=optional_number(row["nominal_original_score"], -1e6, 1e6),
                        measured_later_sectors=int(row["valid_later_sectors"]),
                        product_errors=int(row["product_errors"]),
                        capable_catalog_neighbors=int(row["plausible_neighbor_count"]),
                        matching_later_peaks=int(row["matching_later_peaks"]) if row["matching_later_peaks"] else None)
        solutions = []
        if tic in models:
            fit = models[tic]["nominal"]
            options = diag[tic] if isinstance(diag[tic], list) else [diag[tic]]
            for option in options:
                if option.get("rejected"):
                    continue
                solutions.append(dict(period_days=option["period"], duration_hours=24 * fit["q"][2],
                                      duration_sigma_hours=24 * fit["sigma"][2],
                                      period_sigma_days=option.get("period_sigma_floor", fit["sigma"][1]),
                                      alias_cycles=option.get("integer_cycles")))
        periods = [x["period_days"] for x in solutions]
        status, evidence = prior_reports(tic, periods, catalogs)
        targets.append(dict(tic=tic, gaia_dr3=catalog["GAIAID"],
                            ra_deg=number(product["s_ra"], 0, 360), dec_deg=number(product["s_dec"], -90, 90),
                            coordinate_source="MAST primary product pointing coordinates for target, ICRS",
                            observed=observed, solutions=solutions, prior_status=status, prior_reports=evidence,
                            review=DETAILED.get(tic, ("Pilot only", "Screening and vetting artifacts are available."))[0],
                            caveat=DETAILED.get(tic, ("", "No detailed candidate interpretation."))[1],
                            n_products=len(item["primary"]),
                            missing_product_hashes=sum(fits_manifest[p["dataURL"]] is None for p in item["primary"])))
    by_id = {r["tic"]: r for r in targets}
    predictions = []
    for row in csv_rows(raw["predictions.csv"]):
        tic = row["tic"]
        period = number(row["period_days"], 2, 500)
        model = next(r for r in by_id[tic]["solutions"] if abs(r["period_days"] - period) < 1e-8)
        btjd = number(row["btjd"], 4000, 6000)
        index = len([r for r in predictions if r["tic"] == tic])
        # Sutherland timestamps come from a site-specific barycentric inversion.
        spans = [r for r in visibility["spans"] if r["site"] == "Sutherland"]
        coverage = None
        if tic == "238187210":
            span = spans[index]
            if span["utc_mid"][:10] != row["utc_date_approx"]:
                raise ValueError("UTC site visibility and prediction date disagree")
            coverage = dict(tess="off projected detector", sector=(111, 112, 113)[index], **span)
        predictions.append(dict(tic=tic, solution=row["solution"], period_days=period,
                                duration_hours=model["duration_hours"], btjd=btjd,
                                utc_date_approx=row["utc_date_approx"], coverage=coverage,
                                planning_halfwidth_minutes=number(row["midpoint_halfwidth_minutes"], 0, 10000)))
    details = {tic: dict(model=models[tic], diagnostics=diag[tic], timing=timing.get(tic)) for tic in DETAILED}
    # Use compact event summaries on the public site. Full sampled stacks and
    # pixel inputs remain in the immutable snapshot.
    details["238187210"]["event_ledger"] = target_audit["ledger"]
    details["238187210"]["pixel_models"] = read("followup/tic238_pixels.json")
    return dict(schema_version=1, snapshot_date="2026-11-11", time_scale="BJD_TDB minus 2457000",
                targets=targets, predictions=predictions, details=details, centroids=pixel)


def review_queue(targets):
    queue = []
    for row in targets:
        flags = []
        if row["tic"] in DETAILED:
            flags.append("shortlist")
        if row["observed"]["product_errors"] or row["missing_product_hashes"]:
            flags.append("missing_or_failed_product")
        if row["observed"]["original_midpoint_btjd"] is None:
            flags.append("original_measurement_unavailable")
        if row["observed"]["capable_catalog_neighbors"] and row["observed"]["matching_later_peaks"]:
            flags.append("repeat_like_signal_with_capable_neighbor")
        if not flags and int(sha(row["tic"].encode()), 16) % 5 == 0:
            flags.append("deterministic_quality_sample")
        if flags:
            queue.append(dict(tic=row["tic"], reasons=flags, prior_status=row["prior_status"],
                              review_state="queued", caveat=row["caveat"]))
    return queue


def render_plot(tic, model):
    events = model["nominal"]["events"]
    values = [(float(row["center"]), 100 * float(row["depth"])) for row in events]
    epochs, depths = zip(*values)
    x0, x1 = min(epochs), max(epochs)
    y0, y1 = min(0, min(depths)), max(depths) * 1.12
    xspan, yspan = max(x1 - x0, 1), max(y1 - y0, .05)
    dots = "\n".join(
        f'<circle cx="{69 + (x-x0)/xspan*690:.2f}" cy="{212 - (y-y0)/yspan*140:.2f}" r="4"><title>{x:.4f} BTJD, {y:.3f}%</title></circle>'
        for x, y in values)
    ticks = "".join(f'<text x="{69+i*172.5}" y="238" text-anchor="middle">{x0 + i*xspan/4:.1f}</text>' for i in range(5))
    levels = "".join(f'<text x="58" y="{216-i*35}" text-anchor="end">{y0 + i*yspan/4:.2f}</text>' for i in range(5))
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 820 280" role="img">
<title>Fitted event depths for TIC {tic}</title>
<rect width="820" height="280" fill="#f7f8fc"/>
<g font-family="system-ui,sans-serif" fill="#253346" font-size="13">
<text x="24" y="32" font-size="19">Fitted event depths · TIC {tic}</text>
<text x="24" y="57">Depth (%)</text><path d="M69 72v140h690" fill="none" stroke="#91a0ae"/>
{levels}{ticks}<text x="414" y="267" text-anchor="middle">Event midpoint (BTJD)</text>
<g fill="#236a6c">{dots}</g></g></svg>"""


def build(destination):
    raw = snapshot()
    catalogs = load(DATA / "catalogs.json")
    fits_manifest = json.loads(gzip.decompress((DATA / "fits-manifest.json.gz").read_bytes()))
    key = sha(json_bytes([sha(path.read_bytes()) for path in
                          (Path(__file__), ROOT / "import_pilot.py", DATA / "snapshot-manifest.json",
                           DATA / "catalogs.json", DATA / "fits-manifest.json.gz")]))
    result = checkpoint(CACHE / "publication.json", key, lambda: prepare(raw, catalogs, fits_manifest))
    destination.mkdir(parents=True, exist_ok=True)
    for path in (ROOT / "site").iterdir():
        if path.is_file():
            shutil.copyfile(path, destination / path.name)
    atomic(destination / "data.json", json_bytes(result))
    atomic(destination / "review-queue.json", json_bytes(review_queue(result["targets"])))
    atomic(destination / "predictions.csv", raw["predictions.csv"])
    for tic, details in result["details"].items():
        atomic(destination / f"{tic}.svg", render_plot(html.escape(tic), details["model"]).encode())
    print(f"{len(result['targets'])} targets, {len(result['predictions'])} prediction rows -> {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "docs")
    args = parser.parse_args()
    build(args.out)

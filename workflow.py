"""Create, resume, review and export independent Planet Finder runs."""
import argparse
import csv
import io
import json
import math
import re
import shutil
import sys
from contextlib import contextmanager
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import quote, urlsplit
from uuid import uuid4

from audit import atomic, checkpoint, load, prior_reports
from import_pilot import COMMIT, ROOT, json_bytes, sha

DEFAULT_CATALOG = dict(name="eumemic", kind="community", repository="eumemic/exoplanets", ref=COMMIT,
                       path="candidates/ephemerides.csv", tic_column="TIC", period_column="period_d",
                       label_column="candidate_id")
DEFAULT_CATALOG_DATE = "2026-11-04"
PINNED_PRIOR_REASON = "pinned_historical_prior_catalog"
PINNED_PRIOR_WARNING = (
    f"The default Eumemic catalog is pinned to {DEFAULT_CATALOG_DATE}. Newer entries in that catalog "
    "and other prior work may be unchecked. Start a new run with examples/refresh-priors.json to check a newer revision."
)
DEFAULTS = dict(max_sectors=6, download_mib=72, min_period_days=1.0, max_period_days=35.0,
                durations_days=[.08, .16, .24], grid_points=8000, seeds={}, prior_sources=[DEFAULT_CATALOG])
RUNS = ROOT / "runs"
SHARED = ROOT / ".audit-cache" / "shared"
MIB = 1024 * 1024
MAX_PRODUCT_BYTES = 12 * MIB
ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,62}")
TIC = re.compile(r"[1-9][0-9]{0,15}")
PIPELINES = ["SPOC", "TESS-SPOC", "QLP", "T16"]
MAST = "https://mast.stsci.edu/api/v0/invoke"
NASA = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
REVISION = re.compile(r"[0-9a-f]{40}")


def now():
    return datetime.now(UTC).isoformat(timespec="seconds")


def write(path, value):
    atomic(path, json_bytes(value))


def validate_sources(sources):
    if not isinstance(sources, list) or len(sources) > 6:
        raise ValueError("Specify at most six prior_sources")
    fields = {"name", "kind", "tic_column", "period_column", "label_column"}
    names = {"toi", "confirmed"}
    for source in sources:
        if not isinstance(source, dict) or set(source) not in (fields | {"url"}, fields | {"repository", "ref", "path"}):
            raise ValueError("Prior source requires CSV columns, name and kind, plus url or repository/ref/path")
        if any(not isinstance(value, str) or not value.strip() for value in source.values()):
            raise ValueError("Prior source values must be nonempty strings")
        name = source["name"]
        if not ID.fullmatch(name) or name in names or source["kind"] not in {"community", "ctoi", "eclipsing_binary", "tce"}:
            raise ValueError("Invalid, duplicate or reserved prior source name or kind")
        names.add(name)
        if "url" in source:
            url = urlsplit(source["url"])
            if url.scheme != "https" or not url.hostname or url.username is not None or url.password is not None:
                raise ValueError("Prior CSV URL must use HTTPS with no credentials")
        else:
            repo = source["repository"]
            ref = source["ref"]
            path = source["path"]
            if (not re.fullmatch(r"[\w.-]+/[\w.-]+", repo, re.ASCII)
                    or not re.fullmatch(r"[\w./-]{1,160}", ref, re.ASCII)
                    or not re.fullmatch(r"[\w./-]{1,240}\.csv", path, re.ASCII)
                    or any(part in {"", ".", ".."} for value in (repo, ref, path) for part in value.split("/"))):
                raise ValueError("Invalid GitHub repository, ref or CSV path")


def validate_config(config):
    if not isinstance(config, dict) or set(config) - {"tics", *DEFAULTS}:
        raise ValueError("Unknown settings. See the example JSON configs.")
    config = DEFAULTS | config
    tics = config.get("tics")
    if not isinstance(tics, list) or not 1 <= len(tics) <= 8 or any(type(tic) is not str or not TIC.fullmatch(tic) for tic in tics):
        raise ValueError("Supply one to eight TIC identifiers as strings")
    if len(tics) != len(set(tics)):
        raise ValueError("TIC identifiers must be unique")
    for name, low, high in (("max_sectors", 1, 12), ("download_mib", 1, 128), ("grid_points", 50, 12000)):
        value = config[name]
        if type(value) is not int or not low <= value <= high:
            raise ValueError(f"{name} must be an integer in [{low}, {high}]")
    lo, hi = config["min_period_days"], config["max_period_days"]
    if any(type(x) not in (int, float) for x in (lo, hi)) or not .5 <= lo < hi <= 100:
        raise ValueError("Require .5 <= min_period_days < max_period_days <= 100")
    ds = config["durations_days"]
    if not isinstance(ds, list) or not 1 <= len(ds) <= 5 or any(
        type(d) not in (int, float) or not .04 <= d <= min(.5, .4 * lo) for d in ds
    ):
        raise ValueError("Set one to five durations in days, each less than 40% of minimum period")
    seeds = config["seeds"]
    if not isinstance(seeds, dict) or set(seeds) - set(tics):
        raise ValueError("Seed TIC must be in tics")
    for seed in seeds.values():
        if not isinstance(seed, dict) or set(seed) != {"period_days", "epoch_btjd", "duration_days"}:
            raise ValueError("Seed needs period_days, epoch_btjd and duration_days")
        period, epoch, dur = (seed[k] for k in ("period_days", "epoch_btjd", "duration_days"))
        if any(type(x) not in (int, float) for x in (period, epoch, dur)) or not (
            .5 <= period <= 500 and 1000 <= epoch <= 10000 and .04 <= dur <= min(1, .4 * period)
        ):
            raise ValueError("Seed outside supported period, duration or BTJD range")
    validate_sources(config["prior_sources"])
    return config


def run_path(run_id):
    if not ID.fullmatch(run_id):
        raise ValueError("Invalid run ID")
    return RUNS / run_id


def create_run(config, run_id):
    config = validate_config(config)
    root = run_path(run_id)
    root.mkdir(parents=True)
    try:
        write(root / "run.json", dict(schema_version=1, run_id=run_id, created_at=now(),
                                     config_sha256=sha(json_bytes(config)), config=config))
    except BaseException:
        root.rmdir()
        raise
    return root


def discover(tic):
    import requests

    filters = [dict(paramName="obs_collection", values=["TESS", "HLSP"]),
               dict(paramName="target_name", values=[tic, tic.zfill(16)]),
               dict(paramName="provenance_name", values=PIPELINES),
               dict(paramName="dataproduct_type", values=["timeseries"])]
    query = dict(service="Mast.Caom.Filtered", params=dict(columns="*", filters=filters),
                 format="json", pagesize=2000, page=1)
    # A target with more than one response page needs explicit scope review.
    response = requests.post(MAST, data={"request": json.dumps(query)}, timeout=(8, 35))
    response.raise_for_status()
    result = response.json()
    data = result.get("data")
    if result.get("status") != "COMPLETE" or not isinstance(data, list) or len(data) >= 2000:
        raise ValueError("MAST returned incomplete inventory")
    products = []
    for row in data:
        uri = row.get("dataURL") or ""
        if (str(row.get("target_name", "")).lstrip("0") == tic and row.get("dataRights") == "PUBLIC"
                and uri.startswith("mast:") and uri.endswith(("_lc.fits", "_llc.fits"))
                and row.get("sequence_number") is not None and row.get("provenance_name") in PIPELINES):
            products.append({k: row.get(k) for k in ("dataURL", "obs_id", "provenance_name", "sequence_number",
                                                    "t_exptime", "s_ra", "s_dec", "dataRights")})
    return dict(retrieved_at=now(), source=MAST, products=products, raw_row_count=len(data),
                response_sha256=sha(response.content))


def choose_products(inventory, limit):
    best = {}
    ordered = sorted(inventory["products"], key=lambda p: (PIPELINES.index(p["provenance_name"]),
                     p["t_exptime"] or 1e8, p["dataURL"]))
    for product in ordered:
        best.setdefault(int(product["sequence_number"]), product)
    sectors = sorted(best)
    if len(sectors) > limit:
        indexes = [round(i * (len(sectors) - 1) / (limit - 1)) for i in range(limit)] if limit > 1 else [-1]
        sectors = [sectors[i] for i in indexes]
    return [best[s] for s in sectors]


def product_cache(product, cache):
    uri = product["dataURL"]
    key = sha(uri.encode())
    folder = cache / "fits" / key
    return folder / uri.rsplit("/", 1)[-1], folder / "receipt.json"


def fetch_product(product, cache, budget, local_cache):
    from astropy.io import fits
    from download_curves import open_product

    uri = product["dataURL"]
    path, receipt = product_cache(product, cache)
    if path.is_file() and receipt.is_file():
        try:
            meta = load(receipt)
            if meta["uri"] == uri and path.stat().st_size == meta["bytes"] and sha(path.read_bytes()) == meta["sha256"]:
                return path, meta
        except (OSError, ValueError, KeyError):
            pass
    local = local_cache / path.name if local_cache else None
    path.parent.mkdir(parents=True, exist_ok=True)
    if local and local.is_file():
        if local.stat().st_size > MAX_PRODUCT_BYTES:
            raise ValueError("Local FITS exceeds the per-product limit")
        raw = local.read_bytes()
        method = "user supplied archive cache"
    else:
        allowance = min(MAX_PRODUCT_BYTES, budget[0])
        if allowance < 2880:
            raise ValueError("Download budget exhausted")
        with open_product(uri) as response:
            response.raise_for_status()
            advertised = response.headers.get("Content-Length", "")
            if advertised.isdecimal() and int(advertised) > allowance:
                raise ValueError("Response exceeds the per-product or remaining download limit")
            chunks = []
            read_bytes = 0
            while allowance:
                chunk = response.raw.read(min(65536, allowance), decode_content=True)
                if not chunk:
                    break
                budget[0] -= len(chunk)
                allowance -= len(chunk)
                read_bytes += len(chunk)
                chunks.append(chunk)
            if not allowance and (not advertised.isdecimal() or read_bytes != int(advertised)):
                raise ValueError("Response reached transfer limit without a complete length")
            raw = b"".join(chunks)
        method = "public archive"
    if not raw.startswith(b"SIMPLE  =") or len(raw) % 2880:
        raise ValueError("Response is not a complete block-aligned FITS file")
    with fits.open(io.BytesIO(raw)) as hdus:
        if len(hdus) < 2:
            raise ValueError("FITS has no light-curve extension")
        for hdu in hdus:
            info = hdu.fileinfo()
            if info and info["datLoc"] + info["datSpan"] > len(raw):
                raise ValueError("Truncated FITS extension")
        target = hdus[0].header.get("TICID")
        if target is not None and str(target).strip().removeprefix("TIC ") != product["tic"]:
            raise ValueError("FITS target differs from archive inventory")
    meta = dict(uri=uri, bytes=len(raw), sha256=sha(raw), retrieved_at=now(), method=method)
    atomic(path, raw)
    write(receipt, meta)
    return path, meta


def fetch_catalog(url, args):
    import requests

    response = requests.get(url, params=args, timeout=(8, 20), stream=True)
    with response:
        response.raise_for_status()
        buffer = io.BytesIO()
        for chunk in response.iter_content(65536):
            buffer.write(chunk)
            if buffer.tell() > 4 * MIB:
                raise ValueError("Catalog exceeded transfer limit")
    raw = buffer.getvalue()
    return raw, response.url, response.headers.get("Last-Modified")


def parse_catalog(raw, columns, tic, allow_empty):
    source = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    if not set(columns) <= set(source.fieldnames or []):
        raise ValueError("Unexpected catalog schema")
    id_column, period_column, label_column = columns
    rows = list(source)
    if not rows and not allow_empty:
        raise ValueError("Catalog returned no rows")
    matches = []
    for row in rows:
        if row[id_column] is None or row[id_column].strip().removeprefix("TIC ").strip().lstrip("0") != tic:
            continue
        period = float(row[period_column]) if row[period_column] else None
        if period is not None and (not math.isfinite(period) or period <= 0):
            raise ValueError(f"Invalid period in matching catalog row for TIC {tic}")
        matches.append(dict(tic=tic, period_days=period, label=row[label_column]))
    return rows, matches


def catalog_rows(url, args, kind, columns, tic, allow_empty, *, source=None):
    raw, final_url, modified = fetch_catalog(url, args)
    rows, matches = parse_catalog(raw, columns, tic, allow_empty)
    return dict(kind=kind, state="checked_snapshot" if rows else "checked_no_rows", rows_checked=len(rows),
                url=final_url, sha256=sha(raw), source_name=source or kind,
                retrieved_at=now(), source_modified_at=modified, scope="exact TIC", matches=matches)


def resolve_github(source):
    import requests

    repo, ref = source["repository"], source["ref"]
    if REVISION.fullmatch(ref):
        date = DEFAULT_CATALOG_DATE if source == DEFAULT_CATALOG else None
        return dict(revision=ref, revision_date=date)
    response = requests.get(f"https://api.github.com/repos/{repo}/commits/{quote(ref, safe='')}",
                            timeout=(8, 20))
    response.raise_for_status()
    data = response.json()
    if not REVISION.fullmatch(data.get("sha", "")):
        raise ValueError("GitHub did not return an immutable revision")
    return dict(revision=data["sha"], revision_date=data["commit"]["committer"]["date"])


def source_snapshot(root, source):
    """Freeze an optional catalog and receipt before matching any TIC."""
    folder = root / "prior-snapshots" / source["name"]
    metadata = folder / "receipt.json"
    data = folder / "catalog.csv"
    if metadata.is_file():
        receipt = load(metadata)
        raw = data.read_bytes()
        if sha(raw) != receipt["sha256"]:
            raise ValueError("Frozen catalog SHA-256 mismatch")
        return raw, receipt
    if "repository" in source:
        resolved = checkpoint(folder / "revision.json", sha(json_bytes(source)), lambda: resolve_github(source))
        repo, path, rev = source["repository"], source["path"], resolved["revision"]
        download = f"https://raw.githubusercontent.com/{repo}/{rev}/{path}"
        link = f"https://github.com/{repo}/blob/{rev}/{path}"
    else:
        resolved = dict(revision=None, revision_date=None)
        download = link = source["url"]
    raw, _, modified = fetch_catalog(download, None)
    receipt = dict(url=link, sha256=sha(raw), source_modified_at=modified, retrieved_at=now(),
                   requested_ref=source.get("ref"), **resolved)
    atomic(data, raw)
    write(metadata, receipt)
    return raw, receipt


def catalogs_for(tic, cache, failures, root, sources):
    specs = (
        ("toi", NASA, {"query": f"select tid,pl_orbper,toi from toi where tid={tic}", "format": "csv"}, ("tid", "pl_orbper", "toi"), True),
        ("confirmed", NASA, {"query": f"select tic_id,pl_orbper,pl_name from pscomppars where tic_id='TIC {tic}'", "format": "csv"}, ("tic_id", "pl_orbper", "pl_name"), True),
    )
    records = []
    for kind, url, args, cols, allow_empty in specs:
        key = sha(json_bytes([kind, url, args, cols]))
        try:
            if kind != "community":
                key = stage_key(key, now()[:10])
            records.append(checkpoint(cache / "catalogs" / f"{tic}-{kind}.json", key,
                                      lambda: catalog_rows(url, args, kind, cols, tic, allow_empty)))
        except (ValueError, OSError) as exc:
            failures.append(dict(stage="catalog", source=kind, error=str(exc)))
            records.append(dict(kind=kind, state="failed", url=url, rows_checked=0, matches=[]))
    for source in sources:
        name, kind = source["name"], source["kind"]
        try:
            raw, receipt = source_snapshot(root, source)
            columns = [source[k] for k in ("tic_column", "period_column", "label_column")]
            rows, matches = parse_catalog(raw, columns, tic, allow_empty=False)
            records.append(dict(kind=kind, source_name=name, state="checked_snapshot",
                                rows_checked=len(rows), scope="exact TIC", matches=matches, **receipt))
        except (ValueError, OSError, KeyError) as exc:
            failures.append(dict(stage="catalog", source=name, error=str(exc)))
            records.append(dict(kind=kind, source_name=name, state="failed", rows_checked=0, matches=[],
                                requested_ref=source.get("ref"), url=source.get("url")))
    return records


def stage_key(*values):
    return sha(json_bytes(values))


def run_target(root, config, tic, cache, budget, local_cache, code):
    from screen import read_curve
    from search import search
    from urllib3.exceptions import HTTPError

    target_dir = root / "targets" / tic
    path = target_dir / "status.json"
    failures = []
    stages = {}

    def stage(name, state):
        stages[name] = state
        write(path, dict(tic=tic, stages=stages, failures=failures, updated_at=now()))
        print(f"  {name}: {state}", flush=True)

    stage("discovery", "running")
    try:
        # Keep the same frozen inventory when resuming this run. Newly created
        # runs reuse the per-TIC archive response on the same UTC day.
        def get_inventory():
            return checkpoint(cache / "inventory" / f"{tic}.json", stage_key(code, tic, now()[:10]), lambda: discover(tic))

        inventory = checkpoint(target_dir / "inventory.json", stage_key(1, tic, "MAST"), get_inventory)
        products = [dict(tic=tic, **p) for p in choose_products(inventory, config["max_sectors"])]
        stage("discovery", "complete" if products else "no_products")
    except (OSError, ValueError) as exc:
        failures.append(dict(stage="discovery", error=str(exc)))
        stage("discovery", "failed")
        products = []
    stage("ingestion", "running")
    curves, inputs = [], []
    input_lock = target_dir / "fits-lock.json"
    locked = load(input_lock) if input_lock.is_file() else {}
    for product in products:
        try:
            file, receipt = fetch_product(product, cache, budget, local_cache)
            uri = receipt["uri"]
            if uri in locked and locked[uri]["sha256"] != receipt["sha256"]:
                raise ValueError("Archive bytes differ from this run's frozen input")
            time, flux, _ = read_curve(product, path=file)
            if uri not in locked:
                locked[uri] = receipt
                write(input_lock, locked)
            curves.append(dict(sector=product["sequence_number"], time=time, flux=flux))
            inputs.append(dict(**product, **locked[uri], cadence_count=len(time), first_btjd=float(time[0]),
                               last_btjd=float(time[-1])))
        except (OSError, HTTPError, ValueError, KeyError, IndexError) as exc:
            failures.append(dict(stage="ingestion", uri=product["dataURL"], error=str(exc)))
    stage("ingestion", "complete" if curves and len(curves) == len(products) else "partial" if curves else "unavailable")
    stage("screening_vetting", "running")
    try:
        result = checkpoint(target_dir / "measurements.json", stage_key(code, config, inputs),
                            lambda: search(curves, config, config["seeds"].get(tic)))
        stage("screening_vetting", "complete" if result["status"] == "measured" else result["status"])
    except (OSError, ValueError, ArithmeticError) as exc:
        failures.append(dict(stage="screening_vetting", error=str(exc)))
        result = {"status": "failed"}
        stage("screening_vetting", "failed")
    stage("catalog", "running")
    catalog_root = target_dir / "catalogs.json"
    if catalog_root.is_file():
        catalogs = load(catalog_root)
    else:
        catalogs = catalogs_for(tic, cache, failures, root, config["prior_sources"])
        if all(r["state"] != "failed" for r in catalogs):
            write(catalog_root, catalogs)
    stage("catalog", "incomplete" if any(r["state"] == "failed" for r in catalogs) else "complete")
    period = result.get("period_days")
    prior, records = prior_reports(tic, [period] if period else [], catalogs)
    reasons = list(result.get("flags", []))
    if failures or result["status"] != "measured":
        reasons.append("incomplete_evidence")
    if result["status"] == "measured" and result.get("supported_cycles", 0) >= 3:
        reasons.append("repeat_like_signal")
    if records:
        reasons.append("prior_report")
    if not reasons and int(sha(tic.encode()), 16) % 5 == 0:
        reasons.append("deterministic_quality_sample")
    prior_warning = None
    if DEFAULT_CATALOG in config["prior_sources"]:
        reasons.append(PINNED_PRIOR_REASON)
        prior_warning = PINNED_PRIOR_WARNING
    return dict(tic=tic, coordinates_icrs_deg=[inputs[0]["s_ra"], inputs[0]["s_dec"]] if inputs else None,
                provenance=inputs, stages=stages, failures=failures, evidence=result, prior_status=prior,
                prior_reports=records, catalog_checks=catalogs, review_reasons=reasons,
                prior_warning=prior_warning)


@contextmanager
def run_lock(root):
    lock = root / "running.lock"
    with lock.open("x", encoding="utf-8") as f:
        f.write(now())
    try:
        yield
    finally:
        lock.unlink()


def evidence_fingerprint(row):
    products = [(p["uri"], p["sha256"], p["cadence_count"]) for p in row["provenance"]]
    checks = [(c.get("source_name", c["kind"]), c["state"]) for c in row["catalog_checks"]]
    # The historical-catalog notice explains coverage already present in the
    # evidence. Adding the notice must not invalidate a prior review.
    reasons = [reason for reason in row["review_reasons"] if reason != PINNED_PRIOR_REASON]
    return stage_key(row["tic"], row["evidence"], products, checks, row["stages"],
                     row["prior_status"], row["prior_reports"], reasons)


def saved_decisions(root):
    return {path.stem: load(path) for path in sorted((root / "decisions").glob("*.json"))}


def review_queue(report, decisions):
    queue = []
    for row in report["targets"]:
        tic = row["tic"]
        decision = decisions.get(tic)
        if not row["review_reasons"] and decision is None:
            continue
        fingerprint = evidence_fingerprint(row)
        state = "queued"
        if decision is not None:
            state = "reviewed" if decision.get("evidence_sha256") == fingerprint else "stale"
        queue.append(dict(tic=tic, reasons=row["review_reasons"], prior_status=row["prior_status"],
                          review_state=state, evidence_sha256=fingerprint, prior_warning=row.get("prior_warning")))
    return queue


def record_review(root, tic, decision, reason):
    with run_lock(root):
        report = load(root / "report.json")
        target = next((row for row in report["targets"] if row["tic"] == tic), None)
        if target is None:
            raise ValueError("TIC has no evidence in this run")
        write(root / "decisions" / f"{tic}.json",
              dict(tic=tic, decision=decision, reason=reason, reviewed_at=now(),
                   evidence_sha256=evidence_fingerprint(target)))
        write(root / "review-queue.json", review_queue(report, saved_decisions(root)))


def process(root, cache, local_cache=None):
    import numpy
    import astropy

    sys.path.insert(0, str(ROOT / "legacy"))
    record = load(root / "run.json")
    if record["config_sha256"] != sha(json_bytes(record["config"])):
        raise ValueError("Run config changed. Create a new run to change settings.")
    config = validate_config(record["config"])
    code = [sys.version.split()[0], numpy.__version__, astropy.__version__, version("requests")]
    code.extend(sha(path.read_bytes()) for path in (Path(__file__), ROOT / "search.py", ROOT / "audit.py",
                                                   ROOT / "legacy/screen.py", ROOT / "legacy/download_curves.py"))
    budget = [config["download_mib"] * MIB]
    with run_lock(root):
        rows = []
        for tic in config["tics"]:
            print(f"{record['run_id']} TIC {tic}", flush=True)
            rows.append(run_target(root, config, tic, cache, budget, local_cache, code))
        report = dict(**record, processed_at=now(), code_fingerprint=stage_key(code),
                      time_scale="BJD_TDB minus 2457000", new_fits_bytes=config["download_mib"] * MIB - budget[0],
                      targets=rows)
        write(root / "report.json", report)
        write(root / "review-queue.json", review_queue(report, saved_decisions(root)))
    return report


def plan(report, start, end):
    if not 1000 <= start < end <= 10000 or end - start > 366:
        raise ValueError("Enter a positive BTJD interval of at most 366 days")
    rows = []
    for row in report["targets"]:
        m = row["evidence"]
        if m["status"] != "measured" or m.get("supported_cycles", 0) < 3:
            continue
        epoch, period = m["epoch_btjd"], m["period_days"]
        for cycle in range(math.ceil((start - epoch) / period), math.ceil((end - epoch) / period)):
            rows.append(dict(tic=row["tic"], cycle=cycle, btjd=epoch + period * cycle,
                             duration_hours=24 * m["duration_days"], timing_uncertainty_minutes=None,
                             future_coverage="not evaluated", mode=m["mode"]))
    return sorted(rows, key=lambda r: r["btjd"])


def export(root, out):
    report = load(root / "report.json")
    decisions = saved_decisions(root)
    queue = review_queue(report, decisions)
    if out.exists():
        raise FileExistsError(f"Export already exists: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=out.parent) as temporary:
        staging = Path(temporary)
        for path in (ROOT / "site").iterdir():
            if path.is_file():
                shutil.copyfile(path, staging / path.name)
        write(staging / "data.json", dict(**report, review_queue=queue, decisions=decisions))
        write(staging / "manifest.json", {path.name: sha(path.read_bytes()) for path in sorted(staging.iterdir())})
        staging.rename(out)


def target_list(path):
    raw = path.read_text(encoding="utf-8-sig")
    if path.suffix.lower() != ".csv":
        return raw.split()
    rows = csv.DictReader(io.StringIO(raw))
    columns = {name.lower().replace(" ", "").replace("_", ""): name for name in rows.fieldnames or []}
    column = columns.get("tic") or columns.get("ticid")
    if column is None:
        raise ValueError("CSV requires a TIC or TIC ID column")
    return [row[column].strip().removeprefix("TIC ").strip() for row in rows]


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    new = commands.add_parser("new")
    inputs = new.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--config", type=Path)
    inputs.add_argument("--tic", action="append")
    inputs.add_argument("--tic-file", type=Path, help="Text TIC list or CSV with a TIC column")
    new.add_argument("--id", default=datetime.now(UTC).strftime("%Y%m%d") + "-" + uuid4().hex[:8])
    for command in ("run", "status", "review", "export", "plan"):
        sub = commands.add_parser(command)
        sub.add_argument("id")
        if command == "run":
            sub.add_argument("--local-fits", type=Path, help="Optional archive cache with original FITS filenames")
        if command == "review":
            sub.add_argument("--tic", required=True)
            sub.add_argument("--decision", choices=("needs_followup", "likely_false_positive", "insufficient_data", "prior_recovery"), required=True)
            sub.add_argument("--reason", required=True)
        if command == "export":
            sub.add_argument("--out", required=True, type=Path)
        if command == "plan":
            sub.add_argument("--start-btjd", required=True, type=float)
            sub.add_argument("--end-btjd", required=True, type=float)
    args = parser.parse_args()
    root = run_path(args.id)
    if args.command == "new":
        config = load(args.config) if args.config else {"tics": args.tic or target_list(args.tic_file)}
        print(create_run(config, args.id))
    elif args.command == "run":
        process(root, SHARED, args.local_fits)
    elif args.command == "status":
        record = load(root / "run.json")
        print(f"{record['run_id']} · {len(record['config']['tics'])} targets · {record['created_at']}")
        for path in sorted((root / "targets").glob("*/status.json")):
            print(path.read_text(encoding="utf-8").strip())
    elif args.command == "review":
        record_review(root, args.tic, args.decision, args.reason)
    elif args.command == "export":
        export(root, args.out)
        print(args.out)
    else:
        print(json_bytes(plan(load(root / "report.json"), args.start_btjd, args.end_btjd)).decode(), end="")


if __name__ == "__main__":
    try:
        cli()
    except (ValueError, OSError) as error:
        raise SystemExit(str(error))

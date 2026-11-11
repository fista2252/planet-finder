"""Create, resume, review and export independent Planet Finder runs."""
import argparse
import csv
import io
import json
import math
import re
import shutil
import sys
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory
from uuid import uuid4

from audit import atomic, checkpoint, load, prior_reports
from import_pilot import COMMIT, COMMUNITY, ROOT, json_bytes, sha

DEFAULTS = dict(max_sectors=6, download_mib=72, min_period_days=1.0, max_period_days=35.0,
                durations_days=[.08, .16, .24], grid_points=8000, seeds={})
RUNS = ROOT / "runs"
SHARED = ROOT / ".audit-cache" / "shared"
MIB = 1024 * 1024
MAX_PRODUCT_BYTES = 12 * MIB
ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,62}")
TIC = re.compile(r"[1-9][0-9]{0,15}")
PIPELINES = ["SPOC", "TESS-SPOC", "QLP", "T16"]
MAST = "https://mast.stsci.edu/api/v0/invoke"
NASA = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
RAW_COMMUNITY = f"https://raw.githubusercontent.com/eumemic/exoplanets/{COMMIT}/candidates/ephemerides.csv"


def now():
    return datetime.now(UTC).isoformat(timespec="seconds")


def write(path, value):
    atomic(path, json_bytes(value))


def validate_config(config):
    if not isinstance(config, dict) or set(config) - {"tics", *DEFAULTS}:
        raise ValueError("Unknown settings. See examples/new-run.json.")
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


def catalog_rows(url, args, kind, columns, tic, allow_empty):
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
    source = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    if not set(columns) <= set(source.fieldnames or []):
        raise ValueError("Unexpected catalog schema")
    id_column, period, label = columns
    rows = list(source)
    if not rows and not allow_empty:
        raise ValueError("Catalog returned no rows")
    matches = [dict(tic=tic, period_days=float(row[period]) if row[period] else None, label=row[label],
                    **({"record": row} if kind == "community" else {}))
               for row in rows if row[id_column].strip().removeprefix("TIC ").strip() == tic]
    return dict(kind=kind, state="checked_snapshot" if rows else "checked_no_rows", rows_checked=len(rows),
                url=COMMUNITY if kind == "community" else response.url, sha256=sha(raw),
                retrieved_at=now(), matches=matches)


def catalogs_for(tic, cache, failures):
    specs = (
        ("community", RAW_COMMUNITY, None, ("TIC", "period_d", "candidate_id"), False),
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
        catalogs = catalogs_for(tic, cache, failures)
        if all(r["state"] != "failed" for r in catalogs):
            write(catalog_root, catalogs)
    stage("catalog", "incomplete" if any(r["state"] == "failed" for r in catalogs) else "complete")
    period = result.get("period_days")
    prior, records = prior_reports(tic, [period] if period else [], catalogs)
    reasons = list(result.get("flags", []))
    if failures or result["status"] != "measured":
        reasons.append("incomplete_evidence")
    if result.get("supported_cycles", 0) >= 3:
        reasons.append("repeat_like_signal")
    if records:
        reasons.append("prior_report")
    if not reasons and int(sha(tic.encode()), 16) % 5 == 0:
        reasons.append("deterministic_quality_sample")
    return dict(tic=tic, coordinates_icrs_deg=[inputs[0]["s_ra"], inputs[0]["s_dec"]] if inputs else None,
                provenance=inputs, stages=stages, failures=failures, evidence=result, prior_status=prior,
                prior_reports=records, catalog_checks=catalogs, review_reasons=reasons)


def process(root, cache, local_cache=None):
    import numpy
    import astropy

    sys.path.insert(0, str(ROOT / "legacy"))
    record = load(root / "run.json")
    config = validate_config(record["config"])
    if record["config_sha256"] != sha(json_bytes(config)):
        raise ValueError("Run config changed. Create a new run to change settings.")
    code = [sys.version.split()[0], numpy.__version__, astropy.__version__, version("requests")]
    code.extend(sha(path.read_bytes()) for path in (Path(__file__), ROOT / "search.py", ROOT / "audit.py",
                                                   ROOT / "legacy/screen.py", ROOT / "legacy/download_curves.py"))
    budget = [config["download_mib"] * MIB]
    lock = root / "running.lock"
    with lock.open("x", encoding="utf-8") as f:
        f.write(now())
    try:
        rows = []
        for tic in config["tics"]:
            print(f"{record['run_id']} TIC {tic}", flush=True)
            rows.append(run_target(root, config, tic, cache, budget, local_cache, code))
        report = dict(**record, processed_at=now(), code_fingerprint=stage_key(code),
                      time_scale="BJD_TDB minus 2457000", new_fits_bytes=config["download_mib"] * MIB - budget[0],
                      targets=rows)
        write(root / "report.json", report)
        queue = [dict(tic=r["tic"], reasons=r["review_reasons"], prior_status=r["prior_status"],
                      review_state="queued") for r in rows if r["review_reasons"]]
        write(root / "review-queue.json", queue)
    finally:
        lock.unlink()
    return report


def plan(report, start, end):
    if not 1000 <= start < end <= 10000 or end - start > 366:
        raise ValueError("Enter a positive BTJD interval of at most 366 days")
    rows = []
    for row in report["targets"]:
        m = row["evidence"]
        if m["status"] != "measured":
            continue
        epoch, period = m["epoch_btjd"], m["period_days"]
        for cycle in range(math.ceil((start - epoch) / period), math.ceil((end - epoch) / period)):
            rows.append(dict(tic=row["tic"], cycle=cycle, btjd=epoch + period * cycle,
                             duration_hours=24 * m["duration_days"], timing_uncertainty_minutes=None,
                             future_coverage="not evaluated", mode=m["mode"]))
    return sorted(rows, key=lambda r: r["btjd"])


def export(root, out):
    report = load(root / "report.json")
    queue = load(root / "review-queue.json")
    decisions = {path.stem: load(path) for path in sorted((root / "decisions").glob("*.json"))}
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
        if args.tic not in load(root / "run.json")["config"]["tics"]:
            raise ValueError("TIC is outside this run")
        write(root / "decisions" / f"{args.tic}.json",
              dict(tic=args.tic, decision=args.decision, reason=args.reason, reviewed_at=now()))
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

"""Run bounded offline maintenance checks and persist a weekly ledger."""
import argparse
import json
import re
import subprocess
import sys
import time
from contextlib import contextmanager
from datetime import UTC, datetime, timedelta, timezone
from pathlib import Path
from tempfile import TemporaryDirectory
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

from audit import atomic, load
from import_pilot import ROOT, json_bytes, sha

SPECS = ROOT / "maintenance_specs"
LEDGER = ROOT / ".audit-cache/maintenance"
MINUTES = 30
MAX_ATTEMPTS = 2
FINGERPRINT_FILES = (
    "workflow.py", "search.py", "audit.py", "import_pilot.py", "legacy/screen.py",
    "legacy/download_curves.py", "maintenance.py", "maintenance_benchmark.py",
    "requirements.txt", "maintenance_specs/benchmark-v1.json", "maintenance_specs/registry.json",
    "maintenance_specs/backlog.json", "data/catalogs.json", "data/snapshot-manifest.json",
    "data/fits-manifest.json.gz", "examples/results/example-recovery-v2.json",
    "examples/results/example-unseeded-v2.json",
)


def now():
    return datetime.now(UTC)


def week_start(stamp):
    try:
        karachi = ZoneInfo("Asia/Karachi")
    except ZoneInfoNotFoundError:
        karachi = timezone(timedelta(hours=5))
    local_date = stamp.astimezone(karachi).date()
    return (local_date - timedelta(days=local_date.weekday())).isoformat()


@contextmanager
def ledger_lock(folder):
    folder.mkdir(parents=True, exist_ok=True)
    lock = folder / ".maintenance.lock"
    with lock.open("x", encoding="utf-8") as f:
        f.write(now().isoformat())
    try:
        yield
    finally:
        lock.unlink()


def run_command(args, deadline, timeout):
    remaining = deadline - time.monotonic()
    start = time.monotonic()
    shown = ["python" if arg == sys.executable else str(arg) for arg in args]
    shown = [arg.replace(str(ROOT), ".") for arg in shown]
    if remaining <= 0:
        return dict(command=shown, status="blocked_environment", reason="Weekly time budget exhausted")
    try:
        result = subprocess.run(args, cwd=ROOT, capture_output=True, encoding="utf-8", errors="replace",
                                timeout=min(remaining, timeout))
    except (OSError, subprocess.TimeoutExpired) as exc:
        return dict(command=shown, status="blocked_environment", reason=type(exc).__name__,
                    elapsed_seconds=round(time.monotonic() - start, 3))
    def clean(value):
        return value.replace(str(ROOT), ".").replace(sys.executable, "python")
    return dict(command=shown, status="passed" if result.returncode == 0 else "failed",
                exit_code=result.returncode, stdout=clean(result.stdout), stderr=clean(result.stderr),
                elapsed_seconds=round(time.monotonic() - start, 3))


def repository(deadline, remote, dry_run):
    commands = [run_command(["git", *args], deadline, 5) for args in (
        ("rev-parse", "HEAD"), ("symbolic-ref", "--short", "HEAD"),
        ("status", "--porcelain"), ("rev-parse", "HEAD^{tree}"),
    )]
    if any(x["status"] != "passed" for x in commands):
        return dict(status="blocked_environment", commands=commands, reason="Cannot reconcile local Git state")
    head, branch, dirty, tree = (x["stdout"].strip() for x in commands)
    status, reason = "passed", None
    if dirty and not dry_run:
        status, reason = "blocked_environment", "Shared checkout has changes"
    elif not dry_run and not verified_remote(remote, head):
        status, reason = "blocked_environment", "Current remote revision and successful CI need verification"
    return dict(status=status, reason=reason, head=head, branch=branch, dirty=bool(dirty),
                tree=tree, remote=remote, commands=commands)


def verified_remote(remote, head):
    if not remote or remote.get("commit") != head or remote.get("ci_status") != "success" or not remote.get("ci_url"):
        return False
    try:
        age = now() - datetime.fromisoformat(remote["checked_at"])
    except (KeyError, TypeError, ValueError):
        return False
    return timedelta(0) <= age <= timedelta(hours=1)


def registry():
    spec = load(SPECS / "registry.json")
    backlog = load(SPECS / "backlog.json")
    entries = {row["id"]: row for row in spec["entries"]}
    keys = {"id", "problem", "reference", "license", "expected_benefit", "state", "acceptance",
            "benchmark", "estimated_cost", "dependencies", "rejection_rationale"}
    allowed = {"proposed", "experiment", "validated", "enabled", "rejected"}
    if spec["version"] != 1 or backlog["version"] != 1 or len(entries) != len(spec["entries"]):
        raise ValueError("Unexpected registry version or duplicate IDs")
    for item in entries.values():
        if (not keys <= item.keys() or item["state"] not in allowed or not item["acceptance"]
                or (item["state"] == "rejected" and not item["rejection_rationale"])):
            raise ValueError(f"Invalid registry entry {item['id']}")
    ordered = sorted(backlog["items"], key=lambda item: item["priority"])
    if any(row["id"] not in entries for row in ordered) or len({row["priority"] for row in ordered}) != len(ordered):
        raise ValueError("Invalid backlog links or priorities")
    return entries, ordered


def catalog_check():
    from workflow import DEFAULT_CATALOG, DEFAULT_CATALOG_DATE

    revision_age = (now().date() - datetime.fromisoformat(DEFAULT_CATALOG_DATE).date()).days
    return dict(status="warning", revision=DEFAULT_CATALOG["ref"], revision_date=DEFAULT_CATALOG_DATE,
                age_days=revision_age, scope="Pinned Eumemic CSV plus the configured exact-TIC archive queries",
                coverage_limit="Newer prior work remains unchecked in the default pinned source.",
                refresh_config="examples/refresh-priors.json", live_refresh="not_run")


def examples_check():
    records = []
    for name, status, cycles in (("example-unseeded", "ambiguous_period", 4), ("example-recovery", "measured", 3)):
        path = ROOT / f"examples/results/{name}-v2.json"
        report = load(path)
        target, = report["targets"]
        evidence = target["evidence"]
        ok = (evidence["status"] == status and evidence["supported_cycles"] == cycles
              and report["new_fits_bytes"] == 0 and target["prior_warning"]
              and report["review_queue"][0]["prior_warning"] == target["prior_warning"])
        records.append(dict(path=path.relative_to(ROOT).as_posix(), sha256=sha(path.read_bytes()),
                            status=evidence["status"], supported_cycles=evidence["supported_cycles"], passed=bool(ok)))
    return dict(status="passed" if all(r["passed"] for r in records) else "failed", cached_examples=records)


def compare_benchmarks(baseline, new):
    if baseline["definition_sha256"] != new["definition_sha256"]:
        return ["Benchmark definition differs from baseline"]
    regressions = list(new["violations"])
    for key in ("recovered", "wrong_period", "false_alarm"):
        old, current = baseline["metrics"][key], new["metrics"][key]
        if old["n"] != current["n"] or (current["count"] < old["count"] if key == "recovered"
                                        else current["count"] > old["count"]):
            regressions.append(key)
    return regressions


def checks(deadline):
    results = {}
    results["dependencies"] = run_command([sys.executable, "-c", "import numpy, astropy, requests"], deadline, 30)
    if results["dependencies"]["status"] == "failed":
        results["dependencies"]["status"] = "blocked_environment"
    results["javascript"] = run_command(["node", "--check", "site/app.js"], deadline, 20)
    if results["dependencies"]["status"] != "passed":
        return results, None
    results["regression_suite"] = run_command(
        [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"], deadline, 120)
    with TemporaryDirectory(prefix="planet-finder-health-") as temp:
        results["publication"] = run_command([sys.executable, "audit.py", "--out", temp], deadline, 60)
        publication = results["publication"]
        publication["command"][-1] = "<temporary-publication>"
        for stream in ("stdout", "stderr"):
            if stream in publication:
                publication[stream] = publication[stream].replace(temp, "<temporary-publication>")
        if results["publication"]["status"] == "passed":
            files = sorted(Path(temp).iterdir())
            differing = [path.name for path in files if not (ROOT / "docs" / path.name).is_file()
                         or path.read_bytes() != (ROOT / "docs" / path.name).read_bytes()]
            results["publication"].update(compared_files=len(files), differing_files=differing,
                                          status="failed" if differing else "passed")
    try:
        results["cached_examples"] = examples_check()
    except (OSError, ValueError, KeyError, TypeError) as exc:
        results["cached_examples"] = dict(status="failed", reason=str(exc).replace(str(ROOT), "."))
    result = run_command([sys.executable, "maintenance_benchmark.py"], deadline, 120)
    results["benchmark"] = result
    benchmark = None
    if result["status"] == "passed":
        try:
            benchmark = json.loads(result["stdout"])
            baseline_file = SPECS / "baseline-v1.json"
            baseline = load(baseline_file)
            regressions = compare_benchmarks(baseline, benchmark)
            result.pop("stdout")
            result.update(status="failed" if regressions else "passed", regressions=regressions,
                          baseline_sha256=sha(baseline_file.read_bytes()),
                          baseline_metrics=baseline["metrics"], metrics=benchmark["metrics"])
        except (OSError, ValueError, KeyError, TypeError) as exc:
            result.update(status="failed", reason=str(exc).replace(str(ROOT), "."))
    return results, benchmark


def previous_item(folder, week):
    earlier = sorted(path for path in folder.glob("????-??-??.json") if path.stem < week)
    if earlier:
        prior = load(earlier[-1])
        if prior.get("disposition") in {"active", "retry", "checkpoint", "blocked", "expired", "ready",
                                       "ready_for_review", "awaiting_review"}:
            return prior.get("item"), "awaiting_review" if prior.get("review_pending") else prior["disposition"]
    return None, None


def health(folder, path, *, dry_run, remote):
    if path.is_file():
        return load(path)
    started = now()
    deadline = time.monotonic() + MINUTES * 60
    week = week_start(started)
    registry_check = dict(status="passed")
    try:
        _, backlog = registry()
    except (OSError, ValueError, KeyError, TypeError) as exc:
        backlog = []
        registry_check = dict(status="failed", reason=str(exc).replace(str(ROOT), "."))
    repo = repository(deadline, remote, dry_run)
    results, benchmark = ({}, None) if repo.get("dirty") and not dry_run else checks(deadline)
    results["registry"] = registry_check
    catalog = catalog_check()
    env_block = repo["status"] == "blocked_environment" or any(
        item["status"] == "blocked_environment" for item in results.values())
    defects = [name for name, item in results.items() if item["status"] == "failed"]
    item, previous = previous_item(folder, week) if not dry_run else (None, None)
    disposition = "blocked" if env_block else "needs_investigation" if defects else "ready" if item else "no_change"
    review_pending = bool(item and previous in {"ready_for_review", "awaiting_review"})
    if review_pending:
        if item["id"] in (remote or {}).get("resolved_items", []) and not env_block:
            item = None
            review_pending = False
            disposition = "needs_investigation" if defects else "no_change"
        elif not env_block:
            disposition = "awaiting_review"
    report = dict(schema_version=1, week=week, dry_run=dry_run, started_at=started.isoformat(),
                  deadline_at=(started + timedelta(minutes=MINUTES)).isoformat(), completed_at=now().isoformat(),
                  repository=repo, catalog=catalog, checks=results, benchmark=benchmark, defects=defects,
                  input_hashes={name: sha((ROOT / name).read_bytes()) if (ROOT / name).is_file() else None
                                for name in FINGERPRINT_FILES},
                  disposition=disposition, item=item, attempts=[], review_pending=review_pending,
                  next_item=item or (backlog[0] if backlog else None),
                  limits=dict(window_minutes=MINUTES, attempts=MAX_ATTEMPTS, network="offline checks"),
                  limitations=["Remote CI is supplied by the scheduler for live work.",
                               "Cached reports and synthetic fixtures are checked without fetching FITS or live catalogs."])
    atomic(path, json_bytes(report))
    return report


def claim(record, item, evidence):
    if record["dry_run"] or record["disposition"] not in {"no_change", "needs_investigation", "ready", "retry"}:
        raise ValueError("This ledger cannot start an attempt")
    if now() >= datetime.fromisoformat(record["deadline_at"]):
        raise ValueError("Weekly window has expired")
    if record["item"] and record["item"]["id"] != item:
        raise ValueError("Continue the previously selected item")
    if len(record["attempts"]) >= MAX_ATTEMPTS:
        raise ValueError("Weekly attempt limit reached")
    bug = re.fullmatch(r"BUG-[A-Z0-9-]+", item) and evidence.get("reproduction") and evidence.get("observed_failure")
    if not bug:
        entries, backlog = registry()
        eligible = next((row["id"] for row in backlog if entries[row["id"]]["state"] == "validated"), None)
        if item != eligible or record["defects"]:
            raise ValueError("Select a reproducible bug or the highest-priority validated backlog item")
    remaining = (datetime.fromisoformat(record["deadline_at"]) - now()).total_seconds()
    git = repository(time.monotonic() + min(remaining, 20), None, True)
    if (git["status"] != "passed" or git["dirty"] or not git["branch"].startswith("maintenance/")
            or git["head"] != record["repository"]["head"] and not record["attempts"]):
        raise ValueError("Use a clean maintenance/* branch at the evaluated revision")
    if now() >= datetime.fromisoformat(record["deadline_at"]):
        raise ValueError("Weekly window has expired")
    record["item"] = record["item"] or dict(id=item, kind="bug" if bug else "technique", evidence=evidence)
    record["attempts"].append(dict(number=len(record["attempts"]) + 1, started_at=now().isoformat(),
                                   branch=git["branch"], commit=git["head"], disposition="active"))
    record["disposition"] = "active"


def finish(record, evidence):
    if record["disposition"] != "active":
        raise ValueError("No active attempt to finish")
    requested = evidence.get("disposition")
    if requested not in {"retry", "checkpoint", "blocked", "ready_for_review", "no_change"}:
        raise ValueError("Unknown attempt disposition")
    if not all(isinstance(evidence.get(key), list) for key in ("changed_files", "tests", "limitations", "links")):
        raise ValueError("Record changed_files, tests, limitations and links as lists")
    disposition = "expired" if now() >= datetime.fromisoformat(record["deadline_at"]) else requested
    if disposition == "ready_for_review":
        baseline = record["benchmark"]
        new = evidence.get("benchmark")
        if not baseline or not new or compare_benchmarks(baseline, new) or not evidence["tests"] or not evidence["links"]:
            raise ValueError("Review needs a passing benchmark comparison, tests and review links")
        if record["item"]["kind"] == "technique" and evidence.get("experimental_option") is not True:
            raise ValueError("A new technique needs an experimental option")
    if disposition == "retry" and len(record["attempts"]) >= MAX_ATTEMPTS:
        disposition = "checkpoint"
    record["attempts"][-1].update(disposition=disposition, finished_at=now().isoformat(), evidence=evidence)
    record["disposition"] = disposition
    record["review_pending"] = disposition == "ready_for_review"
    if disposition in {"retry", "checkpoint", "blocked", "expired", "ready_for_review"}:
        record["next_item"] = record["item"]


def add_parser(commands):
    parser = commands.add_parser("maintenance", help=__doc__)
    parser.add_argument("--action", choices=("health", "claim", "finish"), default="health")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--ledger", type=Path, default=LEDGER)
    parser.add_argument("--reconcile", type=Path, help="Remote HEAD and CI verification JSON")
    parser.add_argument("--item", help="Stable item ID for claim")
    parser.add_argument("--evidence", type=Path, help="Reproduction or attempt evidence JSON")


def cli(args):
    if args.dry_run and args.action != "health":
        raise ValueError("Dry runs only execute health checks")
    week = week_start(now())
    folder = args.ledger.resolve()
    path = folder / (f"{week}-dry.json" if args.dry_run else f"{week}.json")
    with ledger_lock(folder):
        if args.action == "health":
            record = health(folder, path, dry_run=args.dry_run, remote=load(args.reconcile) if args.reconcile else None)
        else:
            if not args.evidence or not path.is_file() or args.action == "claim" and not args.item:
                raise ValueError("Start with a live health ledger, then supply evidence and the item for a claim")
            record = load(path)
            evidence = load(args.evidence)
            if args.action == "claim":
                claim(record, args.item, evidence)
            else:
                finish(record, evidence)
            atomic(path, json_bytes(record))
    issues = record["defects"]
    blocked = [name for name, check in record["checks"].items() if check["status"] == "blocked_environment"]
    if record["repository"]["status"] == "blocked_environment":
        blocked.append("repository")
    print(json.dumps(dict(week=week, disposition=record["disposition"], defects=issues,
                          blocked_checks=blocked, catalog_warning=record["catalog"]["coverage_limit"],
                          metrics=record["benchmark"]["metrics"] if record["benchmark"] else None,
                          ledger=str(path))))
    return 1 if issues or record["disposition"] == "blocked" else 0

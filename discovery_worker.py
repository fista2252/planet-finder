"""Run Planet Finder without a model and queue evidence for bounded review."""
import argparse
import json
import math
import re
import sqlite3
import subprocess
import sys
import time
from contextlib import closing, contextmanager
from datetime import UTC, datetime, timedelta, timezone
from pathlib import Path

from audit import atomic, load
from import_pilot import ROOT, json_bytes, sha

STATE = ROOT / ".audit-cache/discovery"
REVIEW_LIMIT = 3
COOLDOWN = 86400
LEASE = 1800
PAKISTAN = timezone(timedelta(hours=5))


@contextmanager
def database(state):
    state.mkdir(parents=True, exist_ok=True)
    with closing(sqlite3.connect(state / "events.sqlite", timeout=5)) as db:
        db.executescript("""
        CREATE TABLE IF NOT EXISTS events (
            id TEXT PRIMARY KEY, subject TEXT NOT NULL, kind TEXT NOT NULL,
            packet TEXT NOT NULL, state TEXT NOT NULL DEFAULT 'pending',
            created INTEGER NOT NULL, lease INTEGER, attempts INTEGER NOT NULL DEFAULT 0,
            outcome TEXT
        );
        CREATE TABLE IF NOT EXISTS claim_budget (day INTEGER NOT NULL, kind TEXT NOT NULL, count INTEGER NOT NULL,
            PRIMARY KEY(day,kind));
        CREATE TABLE IF NOT EXISTS failures (
            fingerprint TEXT PRIMARY KEY, run_id TEXT NOT NULL, count INTEGER NOT NULL
        );
        """)
        yield db


def packet(row, run_id, artifact):
    evidence = row["evidence"]
    period = evidence.get("period_days")
    support = evidence.get("supported_cycles", 0)
    flags = evidence.get("flags", [])
    kind = "uncertain"
    if (evidence["status"] == "measured" and period and math.isfinite(period) and support >= 3
            and evidence.get("mode") == "blind" and evidence.get("sensitive_missing_cycles") == 0
            and row["stages"]["ingestion"] == "complete" and row["stages"]["catalog"] == "complete"
            and row["prior_status"] == "unresolved" and not row.get("prior_reports")
            and not set(flags) - {"period_grid_truncated"}):
        kind = "candidate"
    elif evidence["status"] not in {"ambiguous_period", "inconclusive_search", "measured"} or support < 3:
        return None
    elif row["prior_status"] != "unresolved":
        return None
    # Only signal-relevant results participate in identity. Transfer timestamps,
    # retrieval retries and archive inventory churn cannot trigger another call.
    signal = dict(tic=row["tic"], status=evidence["status"], period_days=period,
                  leading_trial_days=evidence.get("trial_period_days"), supported_cycles=support,
                  sensitive_missing_cycles=evidence.get("sensitive_missing_cycles"), flags=sorted(flags),
                  inputs=sorted((p["uri"], p["sha256"]) for p in row["provenance"]))
    event_id = sha(json_bytes([kind, signal]))[:24]
    return dict(id=event_id, subject=row["tic"], kind=kind, **signal, run_id=run_id,
                depth_fraction=evidence.get("depth_fraction"),
                alternatives=[dict(period_days=a["period_days"], supported_cycles=a["supported_cycles"])
                              for a in evidence.get("alternatives", [])[:3]],
                prior_coverage="Configured catalog checks only. A candidate requires review before any priority claim.",
                artifact=str(artifact))


def ingest(state, report, artifact, stamp=None):
    stamp = int(time.time()) if stamp is None else stamp
    created = 0
    with database(state) as db, db:
        db.execute("BEGIN IMMEDIATE")
        for row in report["targets"]:
            event = packet(row, report["run_id"], artifact)
            if event is None and row["evidence"]["status"] == "failed":
                errors = [x["error"][:240] for x in row["failures"] if x["stage"] == "screening_vetting"
                          and x.get("failure_class") == "computation"]
                if errors:
                    signature = sha(json_bytes([row["tic"], errors]))[:24]
                    db.execute("""INSERT INTO failures VALUES (?, ?, 1)
                        ON CONFLICT(fingerprint) DO UPDATE SET count=count+1,run_id=excluded.run_id
                        WHERE run_id != excluded.run_id""", (signature, report["run_id"]))
                    count, = db.execute("SELECT count FROM failures WHERE fingerprint=?", (signature,)).fetchone()
                    if count >= 2:
                        event = dict(id=signature, subject=row["tic"], kind="repeated_failure",
                                     error=errors, run_id=report["run_id"], artifact=str(artifact))
            if event is None:
                continue
            payload = json_bytes(event)
            created += db.execute("INSERT OR IGNORE INTO events(id,subject,kind,packet,created) VALUES(?,?,?,?,?)",
                                  (event["id"], event["subject"], event["kind"], payload.decode(), stamp)).rowcount
    return created


def claim(state, *, uncertainty=False, stamp=None):
    stamp = int(time.time()) if stamp is None else stamp
    kinds = ("uncertain",) if uncertainty else ("candidate", "repeated_failure")
    day = (stamp + 5 * 3600) // 86400
    lane = "uncertain" if uncertainty else "actionable"
    with database(state) as db, db:
        db.execute("BEGIN IMMEDIATE")
        used = db.execute("SELECT count FROM claim_budget WHERE day=? AND kind=?", (day, lane)).fetchone()
        remaining = REVIEW_LIMIT - (used[0] if used else 0)
        # A crash after handing off a packet may already have incurred a call.
        # Quarantine expired claims until the operator reconciles their outcome.
        db.execute("UPDATE events SET state='uncertain_delivery' WHERE state='leased' AND lease<=?", (stamp,))
        if remaining <= 0:
            return []
        params = ",".join("?" for _ in kinds)
        available = db.execute(f"""SELECT id,subject,packet FROM events AS e
            WHERE kind IN ({params}) AND state='pending' AND attempts<2
            AND (lease IS NULL OR lease<=?)
            AND NOT EXISTS (SELECT 1 FROM events AS sibling WHERE sibling.subject=e.subject
                AND sibling.id!=e.id AND (sibling.state IN ('leased','uncertain_delivery')
                    OR sibling.state='acknowledged' AND sibling.lease>?))
            ORDER BY created,id""", (*kinds, stamp, stamp - COOLDOWN))
        packets, subjects = [], set()
        for event_id, subject, raw in available.fetchall():
            if subject in subjects:
                continue
            db.execute("UPDATE events SET state='leased',lease=?,attempts=attempts+1 WHERE id=?",
                       (stamp + LEASE, event_id))
            packets.append(json.loads(raw))
            subjects.add(subject)
            if len(packets) == remaining:
                break
        db.execute("""INSERT INTO claim_budget VALUES (?,?,?)
                   ON CONFLICT(day,kind) DO UPDATE SET count=count+excluded.count""", (day, lane, len(packets)))
    return packets


def acknowledge(state, event_id, outcome, *, never_delivered=False, stamp=None):
    stamp = int(time.time()) if stamp is None else stamp
    with database(state) as db, db:
        db.execute("BEGIN IMMEDIATE")
        row = db.execute("SELECT state,attempts FROM events WHERE id=?", (event_id,)).fetchone()
        if row is None or row[0] not in {"leased", "uncertain_delivery"}:
            raise ValueError("Event has no claim to reconcile")
        status = "pending" if never_delivered and row[1] < 2 else "exhausted" if never_delivered else "acknowledged"
        deadline = stamp + 300 if status == "pending" else stamp
        db.execute("UPDATE events SET state=?,lease=?,outcome=? WHERE id=?", (status, deadline, outcome, event_id))


def scheduled_tick(state, config_path, slot=None):
    import workflow

    config = workflow.validate_config(load(config_path))
    if config["download_mib"] > 48 or len(config["tics"]) > 4 or config["max_sectors"] > 3:
        raise ValueError("Worker limit is four TICs, three sectors and 48 MiB per run")
    slot = slot or datetime.now(UTC).astimezone(PAKISTAN).date().isoformat()
    if not re.fullmatch(r"\d{4}-\d\d-\d\d|order-[a-z0-9-]{1,32}", slot):
        raise ValueError("Invalid slot or work order ID")
    digest = sha(json_bytes(config))[:8]
    run_id = f"worker-{slot}-{digest}"
    root = workflow.run_path(run_id)
    # The OS releases this database lock if the worker crashes. A separate
    # database keeps queue reads and acknowledgments responsive during search.
    state.mkdir(parents=True, exist_ok=True)
    with closing(sqlite3.connect(state / "worker.sqlite", timeout=1)) as lock:
        lock.execute("BEGIN EXCLUSIVE")
        summary_path = state / f"{run_id}.json"
        if summary_path.is_file():
            return load(summary_path)
        if not (root / "run.json").is_file():
            workflow.create_run(config, run_id)
        path = root / "report.json"
        started = time.monotonic()
        if not path.is_file():
            started_path = state / f"{run_id}.started"
            if started_path.is_file():
                return dict(run_id=run_id, status="interrupted", model_calls=0,
                            reason="Slot transfer budget already reserved. Inspect the run before issuing a new work order.")
            atomic(started_path, json_bytes(dict(started_at=datetime.now(UTC).isoformat())))
            # Invoked without a shell and limited to one archive pass per slot.
            # This supervisor owns the subprocess, and can safely remove the
            # run lock after subprocess.run has killed and reaped a timeout.
            try:
                subprocess.run([sys.executable, "workflow.py", "run", run_id], cwd=ROOT, check=True, timeout=420)
            except subprocess.SubprocessError as error:
                if isinstance(error, subprocess.TimeoutExpired):
                    (root / "running.lock").unlink(missing_ok=True)
                failure = dict(run_id=run_id, status="interrupted", error=type(error).__name__, model_calls=0)
                atomic(summary_path, json_bytes(failure))
                return failure
        count = ingest(state, load(path), path.relative_to(ROOT).as_posix())
        with database(state) as db:
            pending = db.execute("SELECT kind,count(*),sum(length(CAST(packet AS BLOB))) FROM events WHERE state='pending' GROUP BY kind").fetchall()
        summary = dict(run_id=run_id, artifact=path.relative_to(ROOT).as_posix(),
                       new_events=count, pending=[dict(kind=k, count=n, packet_bytes=b) for k,n,b in pending],
                       model_calls=0, wall_seconds=round(time.monotonic() - started, 3))
        atomic(summary_path, json_bytes(summary))
        return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", type=Path, default=STATE)
    commands = parser.add_subparsers(dest="command", required=True)
    tick = commands.add_parser("tick")
    tick.add_argument("--config", type=Path, required=True)
    tick.add_argument("--order", help="Explicit work order ID with a bounded config")
    collect = commands.add_parser("collect")
    collect.add_argument("report", type=Path)
    review = commands.add_parser("claim")
    review.add_argument("--monday-batch", action="store_true")
    ack = commands.add_parser("ack")
    ack.add_argument("event")
    ack.add_argument("outcome")
    ack.add_argument("--never-delivered", action="store_true")
    args = parser.parse_args()
    if args.command == "tick":
        if not args.order and datetime.now(UTC).astimezone(PAKISTAN).weekday() == 0:
            result = dict(status="Monday reserved for maintenance", model_calls=0)
        else:
            result = scheduled_tick(args.state, args.config, f"order-{args.order}" if args.order else None)
    elif args.command == "collect":
        result = dict(new_events=ingest(args.state, load(args.report), args.report.as_posix()), model_calls=0)
    elif args.command == "claim":
        if args.monday_batch and datetime.now(UTC).astimezone(PAKISTAN).weekday() != 0:
            raise ValueError("Uncertainty batch is available on Monday Pakistan time")
        result = claim(args.state, uncertainty=args.monday_batch)
    else:
        acknowledge(args.state, args.event, args.outcome, never_delivered=args.never_delivered)
        result = dict(event=args.event, recorded=True)
    print(json.dumps(result, ensure_ascii=True))
    if isinstance(result, dict) and result.get("status") == "interrupted":
        sys.exit(1)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, sqlite3.Error, subprocess.SubprocessError) as error:
        sys.exit(str(error))

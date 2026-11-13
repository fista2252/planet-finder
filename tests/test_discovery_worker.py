import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import discovery_worker as worker
import workflow
from audit import atomic
from import_pilot import json_bytes


def report(run_id="sample", digest="one", status="measured", **changes):
    row = dict(tic="1234", stages={"ingestion": "complete", "catalog": "complete"},
               prior_status="unresolved", prior_reports=[], failures=[],
               evidence=dict(mode="blind", status=status, period_days=3.5, supported_cycles=4,
                             sensitive_missing_cycles=0, flags=[]),
               provenance=[dict(uri="mast:test", sha256=digest)])
    return dict(run_id=run_id, targets=[dict(row, **changes)])


class WorkerTests(unittest.TestCase):
    def test_claim_budget_applies_across_commands(self):
        with tempfile.TemporaryDirectory() as tmp:
            state = Path(tmp)
            for i in range(5):
                worker.ingest(state, report(tic=str(1000+i)), "artifact", 0)
            self.assertEqual(len(worker.claim(state, stamp=1)), 3)
            self.assertEqual(worker.claim(state, stamp=2), [])
            self.assertEqual(len(worker.claim(state, stamp=86401)), 2)

    def test_deduplication_change_cooldown_and_ack(self):
        with tempfile.TemporaryDirectory() as tmp:
            state = Path(tmp)
            self.assertEqual(worker.ingest(state, report(), "artifact", 10), 1)
            self.assertEqual(worker.ingest(state, report("new-slot"), "elsewhere", 11), 0)
            packet, = worker.claim(state, stamp=12)
            self.assertEqual(packet["kind"], "candidate")
            worker.acknowledge(state, packet["id"], "reviewed", stamp=13)
            self.assertEqual(worker.ingest(state, report(digest="two"), "artifact", 14), 1)
            self.assertEqual(worker.claim(state, stamp=15), [])
            new, = worker.claim(state, stamp=13 + worker.COOLDOWN)
            self.assertNotEqual(new["id"], packet["id"])

    def test_crashed_delivery_requires_reconciliation_and_backoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            state = Path(tmp)
            worker.ingest(state, report(), "artifact", 0)
            packet, = worker.claim(state, stamp=1)
            expired = worker.LEASE + 2
            self.assertEqual(worker.claim(state, stamp=expired), [])
            worker.acknowledge(state, packet["id"], "confirmed no delivery", never_delivered=True, stamp=expired)
            self.assertEqual(worker.claim(state, stamp=expired + 299), [])
            self.assertEqual(worker.claim(state, stamp=expired + 300), [packet])
            worker.acknowledge(state, packet["id"], "confirmed no delivery", never_delivered=True, stamp=expired+301)
            self.assertEqual(worker.claim(state, stamp=expired + 1000), [])

    def test_uncertainty_known_and_transient_failures(self):
        with tempfile.TemporaryDirectory() as tmp:
            state = Path(tmp)
            self.assertEqual(worker.ingest(state, report(status="no_usable_observations"), "a"), 0)
            self.assertEqual(worker.ingest(state, report(prior_status="known catalog signal"), "a"), 0)
            self.assertEqual(worker.ingest(state, report(status="ambiguous_period"), "a"), 1)
            self.assertEqual(worker.claim(state), [])
            self.assertEqual(worker.claim(state, uncertainty=True)[0]["kind"], "uncertain")
            error = [{"stage": "screening_vetting", "failure_class": "computation", "error": "invalid box shape"}]
            self.assertEqual(worker.ingest(state, report(status="failed", failures=error), "a"), 0)
            self.assertEqual(worker.ingest(state, report(status="failed", failures=error), "a"), 0)
            self.assertEqual(worker.ingest(state, report("different", status="failed", failures=error), "a"), 1)

    def test_slot_budget_and_restart_do_not_repeat_transfer(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            cfg = folder / "config.json"
            with patch.object(workflow, "RUNS", folder / "runs"), patch.object(worker, "ROOT", folder), \
                    patch.object(worker.subprocess, "run") as run:
                atomic(cfg, json_bytes(dict(tics=["123"], download_mib=49)))
                with self.assertRaisesRegex(ValueError, "Worker limit"):
                    worker.scheduled_tick(folder / "state", cfg, "order-test")
                atomic(cfg, json_bytes(dict(tics=["123"], download_mib=48, max_sectors=2)))
                run.side_effect = worker.subprocess.TimeoutExpired("fixture", 420)
                result = worker.scheduled_tick(folder / "state", cfg, "order-test")
                self.assertEqual(result["status"], "interrupted")
                self.assertEqual(result, worker.scheduled_tick(folder / "state", cfg, "order-test"))
                run.assert_called_once()

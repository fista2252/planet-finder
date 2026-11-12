import tempfile
import unittest
from datetime import timedelta
from pathlib import Path
from unittest.mock import patch

import maintenance


class MaintenanceTests(unittest.TestCase):
    def test_weekly_health_replay_skips_checks_and_preserves_ledger(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            path = folder / "dry.json"
            with patch.object(maintenance, "repository", return_value={"status": "passed"}), \
                    patch.object(maintenance, "checks", return_value=({}, None)) as checks:
                report = maintenance.health(folder, path, dry_run=True, remote=None)
                original = path.read_bytes()
                self.assertEqual(maintenance.health(folder, path, dry_run=True, remote=None), report)
                self.assertEqual(path.read_bytes(), original)
                checks.assert_called_once()

    def test_weekly_attempts_reserve_one_item_and_stop_after_two(self):
        record = dict(dry_run=False, disposition="needs_investigation", item=None, attempts=[],
                      defects=["regression_suite"], repository={"head": "base"},
                      deadline_at=(maintenance.now() + timedelta(minutes=5)).isoformat())
        evidence = dict(reproduction="python -m unittest tests", observed_failure="Example failed")
        attempt = dict(disposition="retry", changed_files=[], tests=[], limitations=[], links=[])
        with patch.object(maintenance, "repository", return_value=dict(
            status="passed", branch="maintenance/test", dirty=False, head="base"
        )):
            maintenance.claim(record, "BUG-1", evidence)
            maintenance.finish(record, attempt)
            with self.assertRaisesRegex(ValueError, "previously selected"):
                maintenance.claim(record, "BUG-2", evidence)
            maintenance.claim(record, "BUG-1", evidence)
            maintenance.finish(record, attempt)
        self.assertEqual(record["disposition"], "checkpoint")
        self.assertEqual(len(record["attempts"]), 2)
        self.assertEqual(record["next_item"]["id"], "BUG-1")

    def test_expired_attempt_cannot_be_marked_ready_for_review(self):
        record = dict(disposition="active", item={"id": "BUG-1"}, attempts=[{}],
                      deadline_at=(maintenance.now() - timedelta(seconds=1)).isoformat())
        maintenance.finish(record, dict(disposition="ready_for_review", changed_files=[], tests=[],
                                        limitations=[], links=[]))
        self.assertEqual(record["disposition"], "expired")

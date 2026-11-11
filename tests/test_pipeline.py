import argparse
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from astropy.io import fits

import audit
import workflow
from import_pilot import ROOT, json_bytes
from search import search

sys.path.insert(0, str(ROOT / "legacy"))
from screen import read_curve


def series(t0, period=3.25, secondary=False):
    time = t0 + np.arange(0, 22, 1 / 48)
    phase = (time - 1501 + period / 2) % period - period / 2
    flux = np.random.default_rng(733).normal(0, .001, size=len(time))
    flux[abs(phase) < .09] -= .015
    if secondary:
        flux[abs(abs(phase) - period / 2) < .09] -= .012
    return time, flux


def fixture_fits(path, tic, time, flux):
    columns = [fits.Column(name=name, format=code, array=values) for name, code, values in (
        ("TIME", "D", time), ("PDCSAP_FLUX", "D", 1 + flux), ("QUALITY", "J", np.zeros(len(time), dtype=np.int32)),
    )]
    lightcurve = fits.BinTableHDU.from_columns(columns)
    lightcurve.header["TIMESYS"] = "TDB"
    lightcurve.header["BJDREFI"] = 2457000
    primary = fits.PrimaryHDU()
    primary.header["TICID"] = int(tic)
    fits.HDUList([primary, lightcurve]).writeto(path)


class PilotTests(unittest.TestCase):
    def test_frozen_pilot_and_priority_correction(self):
        result = audit.prepare(audit.snapshot(), audit.load(audit.DATA / "catalogs.json"),
                               json.loads(gzip.decompress((audit.DATA / "fits-manifest.json.gz").read_bytes())))
        self.assertEqual(len(result["targets"]), 24)
        targets = {x["tic"]: x for x in result["targets"]}
        self.assertEqual(targets["238187210"]["prior_status"], "recovered prior period")
        prior = next(r for r in targets["238187210"]["prior_reports"] if r["kind"] == "community")
        self.assertEqual(prior["period_days"], 26.0268916)
        self.assertEqual(prior["record"]["registered"], "2026-10-10")
        self.assertEqual(len(targets["339672790"]["solutions"]), 4)
        self.assertEqual(targets["141423461"]["observed"]["original_midpoint_btjd"], None)
        self.assertEqual(len(result["predictions"]), 18)
        self.assertTrue(all(p["coverage"] is None for p in result["predictions"] if p["tic"] != "238187210"))

    def test_alias_with_gap_retains_shortest_period(self):
        raw = audit.snapshot()
        alias = json.loads(raw["followup/diagnostics.json"])["339672790"][0]
        window = alias["predicted_in_product_span"][0]
        self.assertEqual(len(window["offset_measurements"]), 2)
        self.assertFalse(window["rejected_here"])
        self.assertFalse(alias["rejected"])

    def test_failed_catalog_cannot_credit_a_recovery(self):
        catalog = dict(kind="community", state="failed", rows_checked=5, url="https://example.com",
                       matches=[dict(tic="9", period_days=2.0)])
        self.assertEqual(audit.prior_reports("9", [2.0], [catalog]), ("unresolved", []))

    def test_checkpoint_rebuilds_after_mutation(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "stage.json"
            self.assertEqual(audit.checkpoint(path, "hash", lambda: {"score": 1}), {"score": 1})
            saved = audit.load(path)
            saved["value"]["score"] = 9
            audit.atomic(path, json_bytes(saved))
            self.assertEqual(audit.checkpoint(path, "hash", lambda: {"score": 2}), {"score": 2})


class SearchTests(unittest.TestCase):
    def test_seeded_binary_warning_and_local_injection(self):
        time, flux = series(1500, secondary=True)
        result = search([dict(time=time, flux=flux, sector=1)], workflow.DEFAULTS,
                        dict(period_days=3.25, epoch_btjd=1501, duration_days=.18))
        self.assertGreaterEqual(result["supported_cycles"], 5)
        self.assertIn("possible_half_phase_secondary", result["flags"])
        self.assertGreater(result["local_controls"]["recovered_at_or_above_5"], 0)

    def test_observed_cadence_gap_is_unmeasurable(self):
        time, flux = series(1500)
        gap = abs(time - 1504.25) > .4
        result = search([dict(time=time[gap], flux=flux[gap], sector=1)], workflow.DEFAULTS,
                        dict(period_days=3.25, epoch_btjd=1501, duration_days=.18))
        event = next(x for x in result["events"] if x["cycle"] == 1)
        self.assertEqual(event["coverage"], "insufficient_cadence")
        self.assertIsNone(event["measurement"])

    def test_unseeded_period_is_reproducible(self):
        time, flux = series(1500)
        config = workflow.DEFAULTS | dict(min_period_days=2.6, max_period_days=3.8)
        curve = [dict(time=time, flux=flux, sector=1)]
        result = search(curve, config)
        self.assertEqual(result, search(curve, config))
        self.assertAlmostEqual(result["period_days"], 3.25, delta=.02)
        self.assertEqual(result["mode"], "period_grid")
        self.assertGreater(result["search"]["search_snr"], result["search"]["inverted_max_snr"])

    def test_no_data_is_explicit(self):
        self.assertEqual(search([], workflow.DEFAULTS), {"status": "no_usable_observations"})


class RunTests(unittest.TestCase):
    def test_fits_validation_rejects_unexpected_timescale(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "test.fits"
            fixture_fits(path, "123456", *series(1500))
            with fits.open(path, mode="update") as hdus:
                hdus[1].header["TIMESYS"] = "UTC"
            with self.assertRaisesRegex(ValueError, "Unexpected time convention"):
                read_curve({"provenance_name": "SPOC"}, path=path)

    def test_invalid_config_or_duplicate_run_rejected(self):
        for config in ({"tics": [123]}, {"tics": ["999", "999"]}, {"tics": ["999"], "download_mib": 500},
                       {"tics": ["999"], "durations_days": 0}, {"tics": ["999"], "seeds": {"998": {}}}):
            with self.subTest(config=config), self.assertRaises(ValueError):
                workflow.validate_config(config)
        with tempfile.TemporaryDirectory() as temp, patch.object(workflow, "RUNS", Path(temp)):
            workflow.create_run({"tics": ["123"]}, "run")
            with self.assertRaises(FileExistsError):
                workflow.create_run({"tics": ["456"]}, "run")

    def test_fits_download_budget_refuses_before_network(self):
        product = {"dataURL": "mast:fixture/foo_lc.fits", "tic": "123"}
        with tempfile.TemporaryDirectory() as temp, patch("download_curves.open_product") as remote:
            with self.assertRaisesRegex(ValueError, "budget exhausted"):
                workflow.fetch_product(product, Path(temp), [0], None)
            remote.assert_not_called()

    def test_process_resume_export_and_fresh_run(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            local = root / "fits"
            local.mkdir()
            tic = "123456"
            uri = "mast:fixture/valid_lc.fits"
            fixture_fits(local / "valid_lc.fits", tic, *series(1500))
            catalog = dict(source="fixture", products=[
                dict(dataURL=uri, provenance_name="SPOC", sequence_number=2, t_exptime=1800, s_ra=20, s_dec=-60)
            ])
            with patch.object(workflow, "RUNS", root / "runs"), patch.object(workflow, "discover", return_value=catalog) as discover, \
                    patch.object(workflow, "catalogs_for", return_value=[]):
                config = dict(tics=[tic], min_period_days=2.6, max_period_days=3.8)
                run = workflow.create_run(config, "fixture-a")
                first = workflow.process(run, root / "cache", local)
                second = workflow.process(run, root / "cache")
                self.assertEqual(first["targets"][0]["evidence"], second["targets"][0]["evidence"])
                fresh = workflow.create_run(config | {"max_sectors": 1}, "fixture-b")
                workflow.process(fresh, root / "cache")
                discover.assert_called_once_with(tic)
                self.assertNotEqual(audit.load(run / "run.json")["config_sha256"],
                                    audit.load(fresh / "run.json")["config_sha256"])
                out = root / "publication"
                workflow.export(run, out)
                self.assertEqual(audit.load(out / "data.json")["run_id"], "fixture-a")
                self.assertTrue(workflow.plan(first, 1700, 1720))
                with self.assertRaises(FileExistsError):
                    workflow.export(run, out)

    def test_failed_discovery_is_exported_for_review(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(workflow, "RUNS", Path(temp)), \
                patch.object(workflow, "discover", side_effect=ValueError("bad response")), \
                patch.object(workflow, "catalogs_for", return_value=[]):
            run = workflow.create_run({"tics": ["123"]}, "failure")
            report = workflow.process(run, Path(temp) / "cache")
            target = report["targets"][0]
            self.assertEqual(target["stages"]["discovery"], "failed")
            self.assertEqual(target["evidence"]["status"], "no_usable_observations")
            self.assertIn("incomplete_evidence", audit.load(run / "review-queue.json")[0]["reasons"])


if __name__ == "__main__":
    unittest.main()

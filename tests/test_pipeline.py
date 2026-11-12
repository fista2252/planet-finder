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
        self.assertEqual(result["mode"], "blind")
        self.assertGreater(result["search"]["search_snr"], result["search"]["inverted_max_snr"])

    def test_no_data_is_explicit(self):
        self.assertEqual(search([], workflow.DEFAULTS), {"status": "no_usable_observations"})

    def test_sparse_multiyear_search_rejects_coarse_harmonic(self):
        period = 26.026849
        rng = np.random.default_rng(19)
        curves = []
        for sector, start in enumerate([1500, 1865, 2230, 2595, 2960, 3325]):
            time = start + np.arange(0, 27, 1 / 48)
            flux = rng.normal(0, .001, len(time))
            phase = (time - 1501 + period / 2) % period - period / 2
            flux[abs(phase) < .12] -= .011
            curves.append(dict(time=time, flux=flux, sector=sector))
        result = search(curves, workflow.DEFAULTS)
        coarse = next(r for r in result["alternatives"] if r["method"] == "bls_grid")
        self.assertAlmostEqual(coarse["period_days"], period / 6, delta=.002)
        self.assertGreaterEqual(coarse["sensitive_missing_cycles"], 20)
        self.assertAlmostEqual(result["trial_period_days"], period, delta=.003)
        self.assertEqual(result["supported_cycles"], 7)
        self.assertEqual(result["status"], "inconclusive_search")
        self.assertIsNone(result["period_days"])
        self.assertEqual(workflow.plan({"targets": [{"tic": "123", "evidence": result}]}, 4300, 4400), [])

    def test_short_harmonic_contradicted_only_with_sensitive_coverage(self):
        time, flux = series(1500)
        trial = dict(period_days=1.625, epoch_btjd=1501, duration_days=.18)
        data = dict(time=time, flux=flux, sector=1)
        with_coverage = search([data], workflow.DEFAULTS, trial)
        self.assertEqual(with_coverage["status"], "contradicted_ephemeris")
        self.assertGreaterEqual(with_coverage["sensitive_missing_cycles"], 4)
        phase = (time - 1501 + 1.625) % 3.25 - 1.625
        keep = abs(abs(phase) - 1.625) > .4
        gap = search([dict(data, time=time[keep], flux=flux[keep])], workflow.DEFAULTS, trial)
        self.assertEqual(gap["sensitive_missing_cycles"], 0)
        self.assertTrue(any(e["event_state"] == "insufficient_cadence" for e in gap["events"]))

    def test_alternating_eclipse_depths_do_not_promote_a_short_harmonic(self):
        time, flux = series(1500, secondary=True)
        config = workflow.DEFAULTS | dict(min_period_days=1.4, max_period_days=3.8)
        result = search([dict(time=time, flux=flux, sector=1)], config)
        self.assertIn("possible_odd_even_difference", result["flags"])
        self.assertEqual(result["status"], "ambiguous_period")
        self.assertTrue(any(abs(a["period_days"] - 3.25) < .03 for a in result["alternatives"]))


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
                workflow.record_review(run, tic, "needs_followup", "Test review")
                second = workflow.process(run, root / "cache")
                self.assertEqual(first["targets"][0]["evidence"], second["targets"][0]["evidence"])
                self.assertEqual(audit.load(run / "review-queue.json")[0]["review_state"], "reviewed")
                stale = json.loads(json.dumps(second))
                stale["targets"][0]["evidence"]["supported_cycles"] -= 1
                queue = workflow.review_queue(stale, workflow.saved_decisions(run))
                self.assertEqual(queue[0]["review_state"], "stale")
                fresh = workflow.create_run(config | {"max_sectors": 1}, "fixture-b")
                workflow.process(fresh, root / "cache")
                discover.assert_called_once_with(tic)
                self.assertNotEqual(audit.load(run / "run.json")["config_sha256"],
                                    audit.load(fresh / "run.json")["config_sha256"])
                out = root / "publication"
                workflow.export(run, out)
                self.assertEqual(audit.load(out / "data.json")["run_id"], "fixture-a")
                self.assertEqual(audit.load(out / "data.json")["review_queue"][0]["review_state"], "reviewed")
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

    def test_catalog_branch_resolved_once_per_run(self):
        source = dict(workflow.DEFAULT_CATALOG, ref="main")
        raw = b"TIC,period_d,candidate_id\n123,26.02,old\n"
        updated = b"TIC,period_d,candidate_id\n123,26.02,new\n"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with patch.object(workflow, "resolve_github", side_effect=[
                dict(revision="a" * 40, revision_date="2026-11-04"),
                dict(revision="b" * 40, revision_date="2026-11-12"),
            ]) as resolve, patch.object(workflow, "fetch_catalog", side_effect=[
                (raw, "https://example.com", None), (updated, "https://example.com", None),
            ]) as fetch:
                first, receipt = workflow.source_snapshot(root / "first", source)
                self.assertEqual(workflow.source_snapshot(root / "first", source), (first, receipt))
                second, new_receipt = workflow.source_snapshot(root / "second", source)
                self.assertEqual((resolve.call_count, fetch.call_count), (2, 2))
                self.assertIn("/" + "a" * 40 + "/", fetch.call_args_list[0].args[0])
                self.assertEqual(receipt["requested_ref"], "main")
                self.assertNotEqual(new_receipt["revision"], receipt["revision"])
                _, matches = workflow.parse_catalog(second, ["TIC", "period_d", "candidate_id"], "123", False)
                self.assertEqual(matches[0]["label"], "new")

    def test_legacy_config_hash_still_resumes(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(workflow, "RUNS", Path(temp)), \
                patch.object(workflow, "discover", return_value={"products": []}), \
                patch.object(workflow, "catalogs_for", return_value=[]):
            run = workflow.create_run({"tics": ["123"]}, "legacy")
            record = audit.load(run / "run.json")
            record["config"].pop("prior_sources")
            record["config_sha256"] = workflow.sha(json_bytes(record["config"]))
            workflow.write(run / "run.json", record)
            self.assertEqual(workflow.process(run, Path(temp) / "cache")["run_id"], "legacy")


if __name__ == "__main__":
    unittest.main()

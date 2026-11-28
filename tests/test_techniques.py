import unittest

import numpy as np

from techniques import compare_windows, measure_window


class ExperimentalTechniquesTests(unittest.TestCase):
    def test_morphology_distinguishes_injected_shapes_and_preserves_gap(self):
        rng = np.random.default_rng(907)
        time = np.arange(1498, 1502, .002)
        distance = abs(time - 1500) / .20
        noise = rng.normal(0, .0004, len(time))
        for injected, template in (
            ("trapezoid", np.clip((.5 - distance) / .25, 0, 1)),
            ("box", (distance <= .5).astype(float)),
        ):
            with self.subTest(template=injected):
                out = measure_window(time, noise - .008 * template, 1500, .20, .001)
                fits = out["fits"]
                self.assertEqual(min(fits, key=lambda k: fits[k]["residual_ss"]), injected)
                self.assertAlmostEqual(fits[injected]["depth_fraction"], .008, places=3)
        gap = abs(time - 1500) > .09
        self.assertIsNone(measure_window(time[gap], noise[gap], 1500, .20, .001))

    def test_correlated_control_windows_downweight_nominal_score(self):
        rng = np.random.default_rng(14)
        time = np.arange(1490, 1518, .005)
        noise = rng.normal(0, .001, len(time))
        correlated = noise.copy()
        for i in range(1, len(correlated)):
            correlated[i] = .95 * correlated[i - 1] + np.sqrt(1 - .95**2) * noise[i]
        eclipse = abs(time - 1504) < .10
        white = compare_windows(time, noise - .007 * eclipse, [1504], .20, .005)
        red = compare_windows(time, correlated - .007 * eclipse, [1504], .20, .005)
        white_fit = white["events"][0]["measurement"]["fits"]["box"]
        red_fit = red["events"][0]["measurement"]["fits"]["box"]
        self.assertGreaterEqual(len(red["controls"]), 8)
        self.assertGreater(red_fit["control_noise_scale"], white_fit["control_noise_scale"])
        self.assertLess(red_fit["adjusted_score"], red_fit["nominal_score"])
        self.assertTrue(all(abs(r["btjd"] - 1504) > .80 for r in red["controls"]))

    def test_reports_insufficient_controls_without_fabricating_adjustment(self):
        rng = np.random.default_rng(867)
        time = np.arange(1499, 1501, .005)
        flux = rng.normal(0, .001, len(time)) - .01 * (abs(time - 1500) < .08)
        out = compare_windows(time, flux, [1500], .16, .005)
        self.assertIsNone(out["events"][0]["measurement"]["fits"]["box"]["adjusted_score"])


if __name__ == "__main__":
    unittest.main()

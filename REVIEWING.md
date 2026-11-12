# Reviewing the evidence

## Pilot targets

| TIC | What the photometry supports | What remains open |
| --- | --- | --- |
| 238187210 | The independent trapezoid fit to 18 retained event timings gives 26.026849 d, with a conditional 95% bootstrap interval of 26.026794–26.026915 d. Median depth 1.12%. Full-cutout deficit comparisons favor the target over the bright neighbor. | The period was already registered by Eumemic at 26.0268916 ± 0.000035 d on 2026-10-10 as TIC238187210_P26.027. An unresolved stellar system remains possible. The processed ESO search returned no records and the raw ESO query failed. |
| 269730302 | Repeating 28.430790 d signal with a 3.66 h conditional duration. | Grazing fit reaches the bound and Gaia RUWE is 2.83. An eclipsing stellar system or blend remains plausible. Priority for the fitted period remains unresolved. |
| 339672790 | Two dips allow 119.78418, 143.74102, 239.56837 and 359.35255 d. | The short alias had two measurable Sector 90 timing offsets out of five. Local injections were detectable at those two offsets. The other offsets lack sufficient cadence, so all four aliases survive. |

The community catalog is pinned to [commit e2eee92](https://github.com/eumemic/exoplanets/blob/e2eee928fd621f75c1621f6b4206d34bb222ae59/candidates/ephemerides.csv). The pilot recovered the published TIC 238187210 candidate and period. Its 26.02683864 d limb-darkened result and 26.026849 d trapezoid result share observations. The original T16 single-event entries for all three are in [Roth et al.](https://arxiv.org/abs/2604.18579) and their [data deposit](https://doi.org/10.7910/DVN/CWEUGW).

## Data dictionary

| Field or artifact | Interpretation |
| --- | --- |
| `tic`, `gaia_dr3` | String identifiers. Gaia DR3 IDs exceed JavaScript's safe integer range. |
| `ra_deg`, `dec_deg`, `coordinates_icrs_deg` | Archive target coordinates, ICRS degrees. Missing coordinates are null. |
| `*_btjd` | BJD TDB minus 2457000, in days. Site-specific UTC appears only where computed. Approximate UTC dates label the archived pilot predictions. |
| `*_depth_fraction`, `measurement.depth` | Positive is a dip, dimensionless relative flux. Box fits can be negative. Multiply by 100 for percent, 1,000,000 for ppm. |
| `*_duration_days`, `duration_hours` | The field name declares the units. Pilot durations come from fit assumptions. Generic durations come from a grid or user seed. |
| `original_score`, `measurement.score` | Local sideband-based nominal score. Red noise and selection effects limit interpretation. |
| `period_sigma_days` | Conditional pilot fit or timing-floor uncertainty, not a statement about alias probability. Generic runs have no calibrated error bar. |
| `planning_halfwidth_minutes` | Pilot planning half-width using fit covariance and timing floors. Excludes transit duration and baseline. |
| `prior_status` | Known catalog signal, recovered prior period, or unresolved. Exact TIC and period agreement are required for the first two labels. |
| `status`, `stages`, `failures` | Explicit measurement and input states. Null, missing, unavailable, incomplete and failed values do not count as passed vetting. |
| `data/snapshot-manifest.json` | SHA-256 hashes of stored and decompressed scientific snapshot files. |
| `data/fits-manifest.json.gz` | Public MAST URI, expected SHA-256 and byte size. A null receipt means the original bytes were unavailable at packaging. |
| `data/catalogs.json` | Rows checked and source hashes for the pilot catalog snapshots, plus only the pilot-matching rows. |
| `run.json`, `report.json`, `review-queue.json` | Live-run configuration, evidence and review handoff. The export is versioned by run ID and an explicit destination name. |

## Assumptions and unavailable tests

The pilot’s alternate light-curve reductions share spacecraft images. They test pipeline sensitivity. No blind holdout was reserved. Difference-image source checks use a limited number of sectors and approximate centroid models. The legacy timing and transit fits use conditional assumptions about limb darkening, duration and period. There is no planet confirmation, survey completeness estimate, calibrated astrophysical false-positive probability, generic PRF fit, or spectroscopic mass inference.

The generic workflow uses a narrower vetter than the imported pilot. The blind search combines a bounded BLS grid with periods derived from pairs of strong sector dips. It compares a shortlist of refined trials and integer harmonics. A trial needs at least three supported cycles to enter the planner. A rival with comparable support, evidence of alternating eclipse depths, or untested comparable trials can leave the period unresolved. This is a bounded heuristic search, with no guarantee of recovering a true period.

For a sensitive non-detection, each of three trial centers spanning ±¼ of the duration must pass the cadence gate. At each center a sampled box at half the median supported depth must reach a nominal score of five, and the measured depth plus five local errors must be shallower than that injection. Two distinct sensitive misses contradict a stable trial. Timing gaps and low sensitivity leave the trial open. The test is conditional on the depth, shape, timing bracket and local noise. Inversion and fixed-grid injection counts do not calibrate a survey false-positive rate.

Generic runs have no automated Gaia or imaging checks, orbital-shape inference, or blend localization. A seeded run applies the input ephemeris. The planner's predictions are conditional projections with future coverage unevaluated.

The frozen source code is under `legacy/`. `import_pilot.py` copied it from the local pilot and adapted the screen reader and vetter for replay. To rebuild deeper historical fit stages, unpack snapshot inputs in their original relative paths beneath `legacy/` and use `requirements-followup.txt`. Some scripts also require alternate archive products, public catalog snapshots or TESS-Point. Those dependencies are separate from the offline site build and the bounded new-run path. Original FITS and external code are retrieved from their cited sources. [MAST T16](https://archive.stsci.edu/hlsp/t16) · [QLP](https://archive.stsci.edu/hlsp/qlp) · [TESS-Point](https://github.com/tessgi/tess-point)

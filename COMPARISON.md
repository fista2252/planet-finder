# Comparison with eumemic/exoplanets

Reviewed the reference at [e2eee92](https://github.com/eumemic/exoplanets/tree/e2eee928fd621f75c1621f6b4206d34bb222ae59), including the README, reviewing guide, MIT license, DAG and lock file, catalog dictionary, validation notes, tests, candidate explorer, transit planner, and GitHub Pages workflow. Its code and frontend assets were not copied. The project credits its published ephemeris for TIC 238187210.

| Capability | Reference repository | Planet Finder |
| --- | --- | --- |
| Input and search | Large TESS search using indexed public light curves, BLS, sharding and several target selections. | Bounded local runs from TIC lists or configs. Public MAST discovery, seeded checks and a generic BLS grid. The 24-target single-event pilot is imported separately. |
| Checkpoints and provenance | Dependency graph tracks inputs, outputs, parameters, code and environment. CI checks the build lock. | Per-run configs, inventories, FITS digests, versioned evidence and per-target checkpoint hashes. Offline CI checks deterministic site artifacts. Shared archive cache across runs. |
| Prior reports | TOIs, CTOIs, confirmed planets, SPOC TCEs, literature and dated community listings. | Generic checks cover two NASA Archive tables plus configurable community, CTOI, EB and TCE CSVs. Each run freezes source bytes and GitHub revisions. Automated literature searches are outside this release. |
| Vetting | LEO-Vetter, multi-sector pixels, catalog neighbors, binary checks, and TRICERATOPS with companion priors. | Generic cadence, sensitive non-detections, harmonic comparisons, parity, half-phase, local injection and inversion. The three imported targets have additional hand-designed pixel, shape, timing and coverage checks. |
| Calibration | Large inverted control population, injection tests and sensitivity and catalog studies. | Real cadence controls and seeded tests on small samples. They do not measure completeness or a false-positive rate. |
| Transit planner | Published candidate ephemerides, uncertainty propagation and predicted TESS windows. | Imported conditional timing and one target's projected TESS footprint. New runs project BTJD windows, with future visibility and timing errors explicitly unavailable. |
| Publication | Candidate downloads, vetting folders, separate explorer and planner, automated GitHub Pages build. | An original explorer for the pilot and any local export, SHA-256 manifest, documented review queue, CI and reproducible static build. New computation runs locally. |

The current deliverable is a reusable local research workflow. The reference has broader survey validation and a more mature vetting stack. The available evidence supports comparisons at the feature level and no ranking of scientific performance.

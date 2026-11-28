# TIC 99769876 independent evidence

**Disposition: recurring eclipse supported by TESS and ZTF. Companion nature unresolved.** The [earlier campaign](README.md) localized deficits near the target in two sectors, with one event in a sector omitted from the original search. Three additional reductions of the same Sector 7 TESS exposures recover the recurring deep dip. T16 TFA2, TGLC calibrated PSF and Eleanor corrected flux give 3.39%, 5.48% and 1.89% dimming respectively. The separately catalogued T16 epoch at BTJD 1496.11022 lacks a comparable dip in the tested reductions.

## Method and result

Fit an exposure-averaged common trapezoid, independent event amplitudes and local baselines to the frozen TESS windows. Fit timing using Sectors 7, 34 and 87 only and preserve original cycle IDs. The TESS-only model gives epoch 1506.1704128152 BTJD, period 19.0858737816 days and 0.24 d box duration. The separately evaluated Sector 33 center falls 0.90 min after the prediction. This conditional fit estimates timing inside previously selected windows and does not establish a formal period uncertainty or exclude all aliases.

Queried the public ZTF DR24 collection with a 3 arcsec cone at (108.520256973°, −18.7895365481°). Kept 735 usable catalog observations after matching within 2 arcsec, verifying zero catflags, dropping duplicates within each OID and retaining series with at least 20 rows. The median OID positions are 0.12–0.14 arcsec from the TIC target. Match each measurement to the frozen clock using barycentric TDB, computing exposure midpoint from MJD UTC plus half of EXPTIME. The published HJD minus reconstructed heliocentric UTC at the listed MJD equaled half of EXPTIME to within 2 ms for this payload. Use per-OID out-of-eclipse medians outside ±0.03 phase and max(reported error, baseline MAD × 1.4826) as the weight denominator's standard deviation. Collapse exposures within each cycle before combining.

Nine ZTF exposures in seven distinct cycles fall inside the 0.24-day window. Mean dimming is 24.9 mmag (2.27% flux), with a nominal score of 4.26. The six displaced controls with at least seven sampled cycles score lower. Leaving out any one cycle yields 19.6–28.5 mmag. This is independent observational support for the TESS prediction and a descriptive statistic with sparse phase coverage. The strongest non-primary phase from the initial QLP sectors (0.56) and the prespecified phase 0.5 produced no reproducible secondary in another TESS sector or in ZTF. Sparse coverage and eccentric binaries limit interpretation. [ZTF API](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html) · [ZTF timestamp definitions](https://irsa.ipac.caltech.edu/onlinehelp/irsaviewer/irsaviewer/overviewTSV.html) · [TGLC method](https://arxiv.org/abs/2301.03704)

Gaia DR3 exact source 2931196677085790208 has combined RV 28.03 ± 3.39 km/s across 30 transits, with no published epoch velocities, robust amplitude, constancy p-value or non-single-star flag. These data cannot constrain the companion mass. A separate public multi-epoch spectrum or a reproducible secondary would help distinguish a stellar companion from other explanations. [Gaia DR3 table definitions](https://gaia.aip.de/metadata/gaiadr3/gaia_source/)

## Provenance

Source checkout 25512680d1e4a9ac2fda85eb6836367f88bf8e61. Raw inputs, query receipts, timing and morphology grids, controls, reproduction code and figure are retained in the isolated local research run `planet-discovery-20261128-tic99769876`. Archive inputs include 426,240 new FITS bytes and 187,941 bytes of ZTF CSV. No change to search thresholds or production code.

| Artifact | SHA-256 |
| --- | --- |
| ZTF DR24 CSV response | `9145a2b4c3c93fc640d49522b72979efe9373fdabddf4fee509045a4a4682952` |
| TESS alternate and secondary results | `f72766901cd880f4d2d7a58f47c5ead98da40450abdafda2d614828e81374c24` |
| TESS-only morphology and ephemeris | `50461e303af6440bba20e0fdd2dffc121e441d1c722478c4df60309337ee165d` |
| ZTF result (mid-exposure convention) | `679399ff2feaa9895ebf6b76e309416bbc2d7fdd91ec54f744764e5dacb84fcd` |
| Gaia DR3 VOTable | `7c047abbc1aa86e5e7759b53d4883c429ab26b713c584f77d851dd7216586975` |

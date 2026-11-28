# Planet Finder campaign, 2026-11-28

The expanded survey completed **93 new T16 TICs and 279 public TESS light-curve products**, with 92 inconclusive outcomes and one ambiguous period. The existing period gate accepted zero. Twenty-six trials showed at least three nominally supported cycles in at least two sectors and were retained for vetting. Original and alternate trial outcomes were retained independently. This is a targeted sample of previously eligible T16 single events, so these counts do not measure survey completeness.

## Dispositions

| TIC | Outcome | Evidence and next check |
| --- | --- | --- |
| 140920261 | Reject as source of the examined dips | Six measured pixel deficit maps across Sectors 8, 9, 61 and 62 localize 2.68–2.87 pixels away from the target. Their centroids are 0.05–0.21 pixel from TIC 140851769, 55.15 arcsec away, and farther from TIC 140851787. The six independent time windows include the original T16 event in Sector 8 and the disputed Sector 9 window. A calibrated image-response fit is still needed to distinguish the two close neighbors conclusively. This supersedes the on-target disposition in the [earlier period checks](../tic140920261-unused-sectors/README.md). |
| 418817415 | Reject as source of the examined dips | Three Sector 88 events and the Sector 8 event at the published T16 time localize 2.24–2.32 pixels from the target, within 0.06–0.10 pixel of TIC 418817427 (46.15 arcsec). The ∼10.006 d trial is a lead for source-specific photometry of that neighbor. |
| 99769876 | Follow up, period unresolved | Nominal 19.08605 d trial adds a Sector 33 event. Sector 33 and 87 pixel centroids are 0.20 and 0.19 pixel from the target. There are comparable or larger residuals in some S87 controls. Selected raw QLP depths span 3.78–4.45%, giving an illustrative uniform-disk occulting size of 0.247 solar radii from the TIC radius. Check binary explanations before considering a planet interpretation. |
| 112847512 | Unresolved, known variable | Two Sector 34 events add nominal support to the 1.04157 d trial. Difference-image locations disagree and a bright neighbor is 11.6 arcsec away. Gaia DR3 marks the target variable. The exact-ID DR3 eclipsing-binary lookup returned no entry. |
| 307815853 | Unresolved | Original T16 Sector 6 FITS contains a 1.125 d gap around the reported dip with only one in-box cadence. Five unused QLP sectors supply 55 tested trial windows: 40 have insufficient sensitivity and 15 lack sufficient cadence. Sector 87 pixel maps have large background changes and inconsistent source locations. This extends the [anchor test](../discovery-20261113/tic307815853-anchor-timing.md). |

## Scope and reproduction

The campaign used source `25512680d1e4a9ac2fda85eb6836367f88bf8e61`. It froze the prior eligibility audit and MAST product inventories, excluded all previous cohort and report IDs, and selected three widely spaced sectors per TIC. The search used unchanged settings of 1–35 d, box durations 0.08, 0.16 and 0.24 d, and an 8,000-point grid. The held-out exercise selected five highest-ranked cross-sector trials without known or nearby prior matches, sensitive misses or a comparable inverted search, plus the top remaining trial over 10 d. Five products were available across four of these six TICs. Three TICs picked up nominal support, including the blended TIC 418817415. These are additional TESS sectors and not independent instruments.

For the pixel checks, align TESScut and QLP timestamps one-to-one within 0.001 d, require quality zero in both, and use the fixed 0.16 d QLP-time box. Interpolate median images from the preceding and following 0.75–2 duration sidebands, subtract the event median image, and remove a plane fit to the rim pixels. The centroid is weighted by positive deficit in the 3-by-3 patch around the peak. Fit a local image-gradient translation proxy and compare windows displaced by ±3, ±6, and ±9 durations when they have sufficient coverage. Actual TESScut position correction columns were constant, so the translation proxy is not independent pointing telemetry. Centroid distances are descriptive and no calibrated host posterior was computed. Published astrometry was used to check catalog identity, not to infer a period. [MAST TESScut](https://mast.stsci.edu/tesscut/docs/getting_started.html) · [TESS source localization methods](https://doi.org/10.1093/mnras/stad3286) · [T16 monotransit paper](https://arxiv.org/html/2604.18579v1)

The expanded research ledger, source scripts, query responses and raw archive inputs are retained in the campaign research directory `planet-discovery-20261128-scale`. New archive-file subtotal: 287,448,536 bytes, excluding smaller metadata. Screen runtime: 506.2 aggregate target seconds. The refreshed catalog files include NASA TOIs and confirmed planets, ExoFOP CTOIs and the pinned Eumemic revision. The Gaia DR3 positions and EB tables and previous TESS EB and T16 periodic releases also informed follow-up. Catalogue non-matches do not establish novelty. The previously blocked manual maintenance dry run passed all 27 tests and matched the fixed benchmark and 10 publication artifacts. No production code or acceptance criteria changed.

| Campaign artifact | SHA-256 |
| --- | --- |
| scan-plan.json | `e2d6000b8474d4cc31ac382ced8b64b32054e2efa71bb1bbe809dd38aa441673` |
| survey-ledger.json | `1ef6b7894c704cad1d983d3aea41f6d1aa60e1abb702c30f1ed3ae36ca1dcb0c` |
| holdout-plan.json | `1ca18dda1baf7df169fc2d717e64ae41bf2c52ccb3b72d294f0ddcfe396a263b` |
| pixel-early-results.json | `ce2f10e3c71a04891913b6d9f374178767290040cb9c52a0f929296bb26d301b` |
| pixel-wide-results.json | `678d5a6313d4a6c112dacb92fa7d86b55d37a5bda53336f9e3ecb2208761571d` |
| pixel-shortlist-results.json | `f55e48d90699eccb36f170dce6ea40349047b9bdb47f1e0cbaa59eb8802ca4c1` |
| pixel-shortlist-heldout.json | `f9884a3a5743763cd902edcccd1534f457779273fbc68d7f49a49790b31412e5` |
| tic307815853-followup.json | `e7d2350a57d6b034d3b4fc92072a68323b40e9c6eb0646d116738fd454eb59e8` |

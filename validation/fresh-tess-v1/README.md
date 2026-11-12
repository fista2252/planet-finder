# Five-target TESS validation

The initial blind run selected no period for any of the five targets. It ranked the catalog period first for WASP-18, retained a catalog-consistent TOI-270 period behind its harmonic, and abstained on the eclipsing stellar system and the nominal negative control. It also exposed a fixable ingestion failure on TOI-1231. The repair restored the missing sector while leaving that period inconclusive.

## Frozen inputs

[Selection and expected outcomes](selection.json) were saved at 12:11:03 UTC on 2026-11-12, before the blind run started at 12:11:21 UTC. None of the five TICs appeared in the 24-target pilot or the development examples. The [configuration](config.json) had no seed periods. It allowed three sectors per target, 96 MiB per run, and a 0.75–35 day period range. The search used the existing grid, thresholds and candidate limits. Period comparisons used a prespecified 0.1% relative tolerance for planet periods and 1% for the rounded stellar binary period.

The [NASA TAP result](nasa-reference.csv) retains the reference rows with LF line endings. The selection records the hash of the original response, and the [summary](summary.json) records the saved CSV hash. Tau Ceti has catalogued RV signals and serves as a nominal negative for transit-like detections. Its entry supplies neither a guaranteed quiet light curve nor evidence that it has no planets. [SIMBAD target](https://simbad.u-strasbg.fr/simbad/sim-id?Ident=tau+Ceti&protocol=html)

## Blind results

| Target | Reference | Leading trial | Outcome |
| --- | --- | --- | --- |
| WASP-18 · 100100827 | 0.94145223 d | 0.94145298 d, 63 supported cycles | Ambiguous. The runner retained a nearby alternative with the same support count and a possible half-phase secondary flag. [NASA Archive](https://exoplanetarchive.ipac.caltech.edu/overview/WASP-18) |
| TOI-270 · 259377017 | 3.35992, 5.66051 and 11.38194 d | 1.88664560 d, 12 supported cycles | Inconclusive. The leading trial is close to one third of TOI-270 c's period. The 5.65993679 d alternative had 11 supported cycles. [NASA Archive](https://exoplanetarchive.ipac.caltech.edu/overview/TOI-270) |
| TOI-1231 · 447061717 | 24.245586 d | 32.23917550 d, one supported cycle | Inconclusive. The runner initially analyzed two sectors. The ingestion repair restored a third sector and left the blind outcome inconclusive. [NASA Archive](https://exoplanetarchive.ipac.caltech.edu/overview/TOI-1231) |
| TIC 120362137 | About 3.28 d for the inner eclipsing binary | 1.64396794 d, 33 supported cycles | Inconclusive. The leading trial is near half the stellar binary period. A 3.28795235 d alternative had 18 supported cycles. The longer outer eclipse cycle lies outside this search. [Borkovits et al.](https://doi.org/10.1038/s41467-026-69223-4) |
| Tau Ceti · 419015728 | Nominal negative for transits | 0.84928703 d, nine supported cycles | Inconclusive. The number of supported cycles shows why the search also needs tests of competing trials and cadence coverage. |

All five targets reached the search and appeared in the review queue. The first run verified 15 FITS digests, ingested 14 files, and transferred 29,053,440 FITS bytes in 92.861 seconds of processing. Nominal exposures ranged from 120 to 1,800 seconds across SPOC, TESS-SPOC and QLP. No catalog or transport failures occurred in this run. The export passed its file-hash manifest check and matched the report and review queue.

### Ingestion repair

TOI-1231 sector 37 had 16,175 finite, quality-zero SPOC samples but a median PDCSAP flux of −1,057,754. The normalizer correctly rejected that input. The product selector had already discarded alternatives for that sector, which left a coverage gap even though MAST listed other products.

The repair tries up to two alternate products per sector within the original transfer limit. It retains rejected attempts and pins the successful product when resuming. The focused rerun reused the original inventory and catalog results. The TESS-SPOC alternative also failed the flux check. The QLP file supplied 2,807 accepted cadences. Ingestion went from partial to complete and the period remained inconclusive. The rerun used 599,040 additional FITS bytes and took 8.799 seconds. The regression test checks that a valid alternate replaces an unusable preferred file and stays pinned on resume.

### What the run exposed

This is an ingestion and search stress test of five targets. The three planet hosts yielded zero accepted periods under the fixed settings, even though two had a matching trial somewhere in the search results. Sparse coverage, multiple signals and close aliases remain material recovery limits. The binary and nominal negative each lacked an accepted period, which gives no calibrated estimate of false-alarm performance. These cases are now public diagnostic examples and should be excluded from claims of independent validation for any tuning they inform.

## Reproduce and inspect

- [Original evidence bundle](original-evidence.json.gz) holds the unmodified initial report, review queue, frozen MAST inventories, FITS receipts and export manifest. It includes the rejected SPOC file's URI and SHA-256. The FITS remain in MAST.
- [Execution record](execution.json) includes the actual CLI steps and timing. [Repair evidence](ingestion-fix.json) holds the separate one-target report, queue, verified FITS receipts and export manifest. [Comparison summary](summary.json) links trial periods and cadence coverage to the independent references.
- To recreate the original code path, use [commit a440c3e](https://github.com/fista2252/planet-finder/commit/a440c3ef8f3512db91f2cf2fc402345c4c44c4ad) with the saved config. Use a fresh run ID and verify retrieved bytes against the recorded receipts before treating a replay as equivalent. Current code applies the alternate-product repair.

```sh
python workflow.py new --id fresh-tess-replay --config validation/fresh-tess-v1/config.json
python workflow.py run fresh-tess-replay
python workflow.py export fresh-tess-replay --out exports/fresh-tess-replay
```

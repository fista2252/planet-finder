# Planet Finder

Planet Finder runs bounded transit searches and follow-up checks on public TESS light curves. Start with a TIC list, review the evidence and archive a versioned report. The [hosted explorer](https://fista2252.github.io/planet-finder/) shows the first imported pilot. The local Python runner creates new runs and can export the same explorer for their results.

**The pilot recovered a previously reported period.** Eumemic listed TIC 238187210 at **26.0268916 d** on **2026-10-10**. The pilot's independent fit near 26.02685 d recovers that period. TIC 269730302 has a grazing fit and Gaia astrometric concerns. TIC 339672790 retains four aliases. The three signals remain unconfirmed in this project. Roth et al. originally published the three single events. [Pinned prior ephemerides](https://github.com/eumemic/exoplanets/blob/e2eee928fd621f75c1621f6b4206d34bb222ae59/candidates/ephemerides.csv) · [Roth et al.](https://arxiv.org/abs/2604.18579)

## Start a run

Use Python 3.13. Create a virtual environment and activate it before installing dependencies.

```sh
python -m venv .venv
python -m pip install -r requirements.txt
python workflow.py new --id example-a --tic-file examples/unseeded.txt
python workflow.py run example-a
python workflow.py status example-a
python workflow.py export example-a --out exports/example-a-v1
python -m http.server 8000 --directory exports/example-a-v1
```

Open [the local explorer](http://127.0.0.1:8000). The run command prints stage progress. The status command shows completed, unavailable, partial and failed stages. `runs/example-a/review-queue.json` contains the selected cases for review. The report gives each target's evidence, catalog queries and FITS receipts.

Use repeated `--tic` flags, a text list with one TIC per line, a CSV with a `TIC` or `TIC ID` column, or a JSON config. The second example applies an explicit prior ephemeris to four selected sectors.

```sh
python workflow.py new --id example-b --config examples/recover-tic238187210.json
python workflow.py run example-b
python workflow.py plan example-b --start-btjd 4370 --end-btjd 4450
python workflow.py review example-b --tic 238187210 --decision prior_recovery --reason "Matches the pinned community ephemeris"
python workflow.py export example-b --out exports/example-b-v1
```

Each run has its own config digest and ID. Use a new ID to change settings. To retry a failed step, run the same command again. Complete inputs and measurements are checked and reused. Changed upstream FITS bytes trigger a discrepancy for the existing run. New reports use their own export directories and SHA-256 manifests. A killed process may leave `running.lock` in its run directory. Remove that file after verifying the process has stopped, then resume.

### How the runner works

1. **Discover.** Query MAST for public SPOC, TESS-SPOC, QLP and T16 light curves. Select one product per sector with documented priority. If the available sectors exceed the limit, spread the selection across the observed baseline. Freeze the inventory and hash its source response.
2. **Ingest.** Cache FITS by MAST URI, verify local bytes and time convention, and record the product, sector, retained cadences, hash and retrieval timestamp. `--local-fits` can supply an existing archive cache with the original filenames. Check that the cache came from the stated archive before using it.
3. **Screen and vet.** Combine Box Least Squares with a bounded search over separations between individual dips, refine trial periods, and compare integer harmonics. Check expected events on the retained cadences for positive support and sensitive non-detections. Ambiguous or inconclusive periods stay out of the planner. A seeded run tests the supplied period and epoch. The report keeps odd/even, half-phase and inversion controls visible.
4. **Crossmatch and review.** Query NASA Exoplanet Archive TOI and confirmed-planet tables, plus configured CSV catalogs. Match exact TIC IDs and flag period agreement within 0.1%. Failed queries remain visible and unresolved means priority is open. The review queue records queued, reviewed and stale cases. An unchanged evidence fingerprint preserves a saved decision when a run resumes. The runner uses no AI service or cloud compute.
5. **Plan and export.** Use BTJD for linear midpoint projections. Generic runs leave future TESS and ground visibility, and timing uncertainty, uncomputed. The imported pilot supplies its previously checked conditional timing and footprint for TIC 238187210. Export only the report files you intend to publish.

Default limits are eight targets per run, six sectors per target, 72 MiB of new FITS per invocation and 12 MiB per FITS. Hard caps are 12 sectors and 128 MiB per invocation. Metadata queries have separate fixed timeouts and bounds. The config also limits the period grid to at most 12,000 points. The runner leaves cases outside its cadence, event or sample limits in the review queue. An unseeded period may be a harmonic or an artifact. Local scores and inversion counts have no calibrated survey false-alarm interpretation.

### Refresh or extend the prior search

The default community source is pinned to the November 4 Eumemic snapshot. The report and review queue flag that cutoff because newer entries may be unchecked. To check a newer revision, create a new run with [the refresh config](examples/refresh-priors.json), which resolves the public GitHub `main` branch once and also checks ExoFOP CTOIs.

```sh
python workflow.py new --id current-priors --config examples/refresh-priors.json
python workflow.py run current-priors
```

Set `prior_sources` to a list of up to six CSV sources. Each entry declares `name`, `kind`, `tic_column`, `period_column` and `label_column`, plus either `repository`, `ref` and `path`, or a public HTTPS `url`. Kinds include `community`, `ctoi`, `eclipsing_binary` and `tce`. Use the EB or TCE kind with an archive CSV that provides a TIC column and periods in days. Built-in NASA checks run in addition to these sources. Each run freezes CSV bytes and their digest, and records the resolved GitHub commit and date when available. An endpoint without a version exposes the content hash, retrieval time and any server modification date. Matching only covers the selected sources and their exact TIC records.

## Imported pilot and verification

```sh
python audit.py
python -m unittest discover -s tests -v
python -m http.server 8000 --directory docs
```

The repository contains the original 24-target pilot's selection ledger, MAST metadata, screening, alias and vetting outputs, FITS hash manifest, and the three-target follow-up. The saved FITS are retrieved separately. `replay.py` remeasures one to three selected pilot targets through the original screen, alias and vetting stages.

```sh
python replay.py --tic 238187210 --curves /path/to/original/curves --controls
# For a fresh cache, omit --curves and request bounded public retrieval with --fetch.
```

See [methods, units and data limitations](REVIEWING.md), [the reference comparison](COMPARISON.md), [attribution](ATTRIBUTION.md), and [the changelog](CHANGELOG.md). Saved live-run checks are in [examples/results](examples/results).

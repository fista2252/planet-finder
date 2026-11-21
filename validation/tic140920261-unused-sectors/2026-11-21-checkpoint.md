# TIC 140920261 investigation checkpoint — 2026-11-21

Disposition: unresolved. No planet or new candidate is established.

## Objective and selection

Resume the saved investigation, without a new period search. The most discriminating next observation is a Sector 62 pixel cutout at ICRS RA 128.659046515 deg, Dec -47.1748911729 deg. The [Sector 62 aperture and position check](sector62-aperture-check.md) shows that both selected dips deepen with aperture size and coincide with reported position changes. Extracted light curves cannot locate the missing light. Pixel deficit maps can distinguish an on-target deficit from a neighbor or an image-motion residual, subject to registration and noise checks.

Preserve the two frozen event centers BTJD 3006.81430 and 3010.66565 and duration 0.16 d. Test events separately using local sidebands and pointing-matched controls; inspect positive/negative image residual structure before interpreting a position. Do not treat different reductions of the same exposures as independent confirmation. The earlier Sector 9 non-detection and reduction discrepancy remain relevant. Keep catalog identification, proper-motion/epoch uncertainty and source-location uncertainty separate from signal detection.

## Earlier computation recovered from the chat handoff

The November 19 cached Sector 61 computation had not been published successfully. The following rounded values are transcribed from that saved computation, not newly measured on November 21. A complete executable export and machine-readable output still need to be saved and independently checked before stronger use.

Input: `mast:HLSP/qlp/s0061/0000/0001/4092/0261/hlsp_qlp_tess_ffi_s0061-0000000140920261_tess_v01_llc.fits`.
SHA-256: `4b2462b61fd3eb7212339eb54c3c0d1b9da9f718ecf48b223e79bdb61dcf0afd`.
Local code checkpoint: `25512680d1e4a9ac2fda85eb6836367f88bf8e61`.
Frozen centers: BTJD 2981.780478007012 and 2985.631834841218; duration 0.16 d. The existing local-sideband scorer used 9,326 common finite quality-zero timestamps. Flux was divided by its median and decreased values denote dimming. Positions retain pixel units.

| Measurement | First event | Second event |
| --- | ---: | ---: |
| SAP depth | 0.903% | 1.023% |
| Systematics-removed depth | 0.891% | 1.021% |
| Small detrended aperture depth | 0.300% | 0.378% |
| Default detrended aperture depth | 0.534% | 0.641% |
| Large detrended aperture depth | 0.905% | 1.007% |
| SAP_X decrease, px | 0.0090213 | 0.0080171 |
| SAP_Y decrease, px | 0.0024008 | 0.0029013 |

Controls used centers `arange(t[0]+0.32,t[-1]-0.32,0.64)`, excluding centers within 0.64 d of either event. Twenty-eight control windows met the scorer requirements. Event position-shift vector magnitudes were 0.0093353 and 0.0085260 px; no retained control exceeded either. This is descriptive, not a calibrated significance or false-alarm probability. The repeating aperture/position pattern across Sectors 61 and 62 leaves blends and pointing effects viable. These sectors use the same instrument; this is not independent-instrument confirmation.

## November 21 execution outcome

Read the current local investigation ledger, maintenance protocol and clean isolated-checkout status. Process enumeration was denied, so no absence of other researchers/processes is claimed. No local files were modified. Inspected the existing remote Sector 62 note and verified that this checkpoint path did not already exist.

The shell initially ran read-only commands. A narrow permission request then returned a turn-scoped grant of network access and write access only to a new isolated research directory, `research_runs/planet-discovery-20261121-tic140920261`. The first directory-creation command failed before process launch with:

`Failed to create unified exec process: workspace-write sandbox has no writable root capability SIDs`

The grant therefore did not produce a usable execution environment. No cutout request was issued, no new archive data were downloaded, no fresh measurement was made, and no shell retry was attempted. The original dirty checkout and clean isolated checkout remain untouched. This checkpoint is published via the already available GitHub connector to preserve the interrupted state; it does not establish that local write/network execution works.

## Resume contract

After the shell sandbox can launch with a writable root, reserve a fresh isolated directory and verify writer state. Download at most one 7 by 7 pixel Sector 62 cutout, streaming with an 8 MiB byte ceiling and bounded timeout; retain HTTP metadata, retrieval time, full archive URL, raw bytes and SHA-256. If larger or unavailable, record that outcome instead of increasing limits silently. Preserve FITS calibration/time/WCS metadata and inspect them before aligning with QLP times. Do not search/refit event timing on these images. Save the frozen protocol before inspecting deficit maps, separate event maps and controls, exact code/dependency versions, failures and runtime. Stop within the 30-minute research budget.

Official API documentation was re-opened November 21: [STScI TESScut getting started](https://mast.stsci.edu/tesscut/docs/getting_started.html). It documents a sector-bounded cutout request and ZIP response. Planned URL: `https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=128.659046515&dec=-47.1748911729&x=7&y=7&units=px&sector=62`. This URL was not fetched in this session.

No novelty claim, catalog promotion, software fix or acceptance-threshold change was made. New astronomical bytes: 0. Numerical-analysis runtime this session: none; earlier Sector 61 runtime was not recorded in the available handoff.

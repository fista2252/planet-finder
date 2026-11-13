# TIC 140920261 follow-up

**Disposition: unresolved, with evidence against a stable on-target 1.92568-day eclipse.** The 1.92568-day trial picked up three nominally supported events in previously unused QLP sectors. Its Sector 9 event disappears in TESS-SPOC SAP and PDCSAP even when both pipelines use the same cadences. That leaves aperture contamination, reduction systematics and an ephemeris error as plausible explanations. Neither reduction establishes the location of the source.

The earlier search used Sectors 8, 61 and 89. This test froze three trial periods, two reference epochs per period and a 0.16-day box before loading QLP Sectors 9, 36 and 62. The injection depth was half the median supported depth measured in the original sectors. A sensitive miss required enough cadences and a detectable injection at all three timing offsets. Two distinct sensitive misses conditionally contradict the trial at the tested depth and timing. [Frozen protocol](protocol.json) · [MAST QLP documentation](https://archive.stsci.edu/hlsp/qlp)

| Period (d) | Epoch (BTJD) | New support | Sensitive misses | Result |
| ---: | ---: | ---: | ---: | --- |
| 3.851357 | 2985.631835 | 0 | 1 | Unresolved |
| 3.851357 | 1522.170900 | 0 | 1 | Unresolved |
| 1.925678 | 2985.631835 | 3 | 1 | QLP favors this trial, with a discrepant Sector 9 reduction |
| 1.925678 | 1522.170900 | 1 | 1 | Unresolved |
| 1.283786 | 2985.631835 | 0 | 1 | Unresolved |
| 1.283786 | 1522.170900 | 0 | 3 | Fixed trial contradicted under the tested assumptions |

All three periods at the search epoch share the sensitive miss near BTJD 1564.48116. The 1.925678-day QLP detections are near BTJD 1547.15006, 3006.81430 and 3010.66565. In Sector 9, 952 timestamps match within 0.002 d. The first dip is +0.501% in QLP and −0.064% in PDCSAP and −0.107% in SAP on the matched cadences, with positive values denoting dimming. The SAP check also gives negative depths at the two timing offsets around that center. TESS-SPOC and QLP use the same TESS images, so this is a reduction and aperture comparison. [Window ledger](result.json) · [Sector 9 comparison](sector9-reduction-check.json)

The exact-TIC TOI and confirmed-planet checks returned no rows. The pinned Eumemic table contained no matching row. The original event is in the T16 catalog, and the coverage of the prior search leaves priority unresolved. [Prior check and revision](prior-check.json) · [T16 source data](https://doi.org/10.7910/DVN/CWEUGW)

The four FITS transfers totaled 1,091,520 bytes. The scripts recorded 5.6 seconds for the archive and measurement passes, excluding the later matched-cadence analysis and catalog checks. Input URIs and SHA-256 digests are recorded with the results. The fixed-window analysis uses the existing scorer and leaves search thresholds unchanged. To rerun it, place the three original FITS in the Planet Finder shared cache at their URI-derived paths, verify their digests, and run `python validation/tic140920261-unused-sectors/check.py` from a replay checkout with the protocol and an empty result path.

The next discriminating test is a Sector 62 pixel difference image at the two fixed QLP event times. A deficit centered on the target would favor a real event on TIC 140920261. A deficit centered elsewhere would favor contamination. The shared Sector 9 miss remains relevant to the assumption of a stable repeating eclipse.

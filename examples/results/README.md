# Live workflow verification

Both runs used public MAST products and live NASA Archive queries on 2026-11-11. Each also checked the pinned Eumemic catalog. The JSON files preserve the initial-release results. Each resume transferred zero new FITS bytes.

| Run | Input | Initial report | Updated cached report |
| --- | --- | --- | --- |
| Unseeded | TIC 183374187, outside the imported pilot. One public T16 light curve. Initial FITS transfer 313,920 bytes. | [Initial JSON](example-unseeded.json) selects 5.058648 d with four supported cycles. | [Updated JSON](example-unseeded-v2.json) retains 5.058648 d as the leading trial and records `ambiguous_period` because the light curve cannot reject shorter harmonics. |
| Seeded recovery | TIC 238187210, seeded with the published ephemeris. Four MAST sectors. Initial FITS transfer 2,416,320 bytes. | [Initial JSON](example-recovery.json) records three supported cycles at the input ephemeris and matches Eumemic's candidate. | [Updated JSON](example-recovery-v2.json) records the same three supported cycles with the new evidence fields. |

Both updated reports include the review queue and the November 4 catalog cutoff notice. They use cached FITS and the original prior inputs. The configurable live catalog refresh was tested with simulated responses because the local runner could not reach the remote endpoint.

[Pilot replay](pilot-replay.json) checks the original three computation stages for TIC 238187210 against the frozen outputs. The extra fixed-phase controls tested 1,049 windows on 21 sectors and recovered 351 half-depth box injections at the nominal threshold of 8. The inverted controls produced zero scores at or above 8 in those windows. These fixed-grid tests describe the chosen light curves, depth, template and threshold.

# Live workflow verification

Both runs used public MAST products and live NASA Archive queries on 2026-11-11. Each also checked the pinned Eumemic catalog. The JSON files preserve the initial-release results. Each resume transferred zero new FITS bytes.

| Run | Input | Result |
| --- | --- | --- |
| [example-unseeded](example-unseeded.json) | TIC 183374187, outside the imported pilot. One available public T16 light curve. Initial FITS transfer 313,920 bytes. | BLS measured 5.058648 d with four supported cycles. The period agrees with the known TIC 183374187 b catalog entry. The search grid hit its resolution limit. |
| [example-recovery](example-recovery.json) | TIC 238187210, seeded with the published ephemeris. Four MAST sectors. Initial FITS transfer 2,416,320 bytes. | Three supported cycles at the input ephemeris. The prior-report check recovered Eumemic's exact candidate ID and period. |

On the same cached FITS, the updated search retained a leading 5.058648-day trial for TIC 183374187 but marked the period ambiguous. The single light curve provided insufficient sensitivity to reject shorter harmonics. The seeded TIC 238187210 result still had three supported cycles. Both reruns used zero new FITS bytes.

[Pilot replay](pilot-replay.json) checks the original three computation stages for TIC 238187210 against the frozen outputs. The extra fixed-phase controls tested 1,049 windows on 21 sectors and recovered 351 half-depth box injections at the nominal threshold of 8. The inverted controls produced zero scores at or above 8 in those windows. These fixed-grid tests describe the chosen light curves, depth, template and threshold.

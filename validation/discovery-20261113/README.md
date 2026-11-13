# Bounded TESS discovery search, November 13, 2026

The run tested four known systems and two additional T16 single-event targets outside the original 24-target pilot. All six searches were blind and inconclusive. The additional T16 targets were selected by the number of later sectors in the frozen ledger, then depth and TIC. The configurations and source expectations were saved before running each batch.

| Target | Leading trial | Supported cycles | Interpretation |
| --- | ---: | ---: | --- |
| WASP-121 | 1.274187 d | 36 | Within 0.1% of the NASA catalog period. Competing trials remained. |
| TOI-700 | 16.055801 d | 2 | Near TOI-700 c's 16.051104-day period. Two supported cycles leave limited evidence. |
| TOI-178 | 0.603816 d | 4 | Leading trial did not match a catalog period within the fixed tolerance. |
| LHS 3844 | 27.318989 d | 0 | Its known 0.462930-day period lies below the configured search floor. |
| TIC 307815853 | 13.926528 d | 1 | An additional cycle is needed before periodic follow-up. |
| TIC 140920261 | 3.851357 d | 4 | Retained for bounded uncertainty review. The trial did not pass the blind period selection. |

The known-system batch transferred 17,910,720 FITS bytes in 66.321 seconds including setup and export. A worker replay transferred zero new FITS bytes. The additional T16 batch transferred 2,373,120 FITS bytes in 26.353 seconds. The worker preserved three uncertainty packets totaling 2,997 UTF-8 bytes across the two worker receipts. Duplicate collection added zero events, actionable claim returned an empty list, and both worker runs made zero model calls. Token usage was unavailable. The byte count measures packet size only.

The configured minimum period also constrains duration. A first configuration with the default 0.24-day duration and a 0.5-day minimum failed validation before any search or transfer. The frozen configuration supplied compatible durations. The search thresholds and benchmark outcomes were unchanged. The trial summaries are diagnostics from a small sample and do not estimate survey performance.

- [Known-system selection](selection.json), [NASA catalog query receipt](reference-receipt.json), [catalog rows](reference.csv), [configuration](config.json) and [report](original-report.json.gz)
- [Additional T16 selection](pilot-extension-selection.json), [configuration](pilot-extension-config.json) and [report with inventory](pilot-extension-report.json.gz)
- [Execution timing](execution.json) and [worker verification](worker-verification.json)

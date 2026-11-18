# TIC 140920261 Sector 62 aperture and position checks

The source remains unresolved. Both dips appear in raw QLP flux, become deeper in larger detrended apertures, and coincide with shifts in the reported image position. A blended variable source and pointing-related systematics remain possible. The earlier Sector 9 reduction discrepancy and sensitive miss still apply. [Sector 9 comparison](sector9-reduction-check.json)

The frozen centers were BTJD 3006.81430 and 3010.66565. The 0.16-day local sideband box used 10,038 common finite, quality-zero timestamps in the verified original QLP FITS. Flux columns were median normalized. Positive values indicate decreases. [Fixed windows and archive provenance](result.json) · [QLP data conventions](https://tess.mit.edu/qlp-update/)

| Measurement | First dip | Second dip |
| --- | ---: | ---: |
| Raw SAP depth | 1.019% | 1.062% |
| Systematics-removed depth | 1.014% | 1.066% |
| Small detrended aperture depth | 0.290% | 0.354% |
| Default detrended aperture depth | 0.596% | 0.588% |
| Large detrended aperture depth | 1.026% | 0.994% |
| SAP_X decrease (px) | 0.00719 | 0.00696 |
| SAP_Y decrease (px) | 0.00650 | 0.00630 |

On the November 18 check, the same scorer measured raw SAP flux and SAP_X/Y at comparison centers spaced by 0.64 d starting 0.32 d after the first retained timestamp. Centers within 0.64 d of either event were excluded. Thirty-one windows passed the cadence requirements. The median and maximum control position shift magnitudes were 0.00063 and 0.00400 px. The event shifts were 0.00969 and 0.00939 px, both larger than every control. Neither the two selected events nor the control grid supplies a calibrated false-alarm probability. This check measures position changes in an extracted light curve and leaves the source of the missing flux unidentified.

This analysis used the cached Sector 62 FITS (SHA-256 `5d7b63c09667fab29a8bf92ef9bf34fd93abff36a5f2948f6e3411e60754904b`) and transferred zero FITS bytes. The scheduled session lacked direct network access to TESScut. The next check is a Sector 62 pixel cutout at RA 128.659046515°, Dec −47.1748911729°. Compare the two fixed-time deficit maps separately against local and pointing-matched controls, then check their locations against the target and catalog neighbors. [STScI TESScut API](https://mast.stsci.edu/tesscut/docs/getting_started.html)

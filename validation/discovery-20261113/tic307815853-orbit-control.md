# TIC 307815853: Sector 87 spacecraft-phase comparison

**Disposition: unresolved.** The two QLP dimmings are about one spacecraft orbit apart, but the flux feature did not recur at either fixed delay from the observation segment start in the three other segments. A large local background change did recur without a comparably deep flux dip. A repeatable artifact at a fixed delay from segment start is therefore an incomplete explanation. A phase-dependent background pattern with changes in geometry, a blend, or other extraction effects remain plausible.

The [Sector 87 release notes, DRN 120 v01, Table 1](https://tasoc.dk/docs/release_notes/tess_sector_87_drn120_v01.pdf) supply timing boundaries for four observation segments. Momentum dumps took place before science collection resumed. The candidate windows fall roughly 2.95 and 3.27 days after starts of 181a and 182a, respectively, well away from the start and end. The released flags and table are context, and the phase test does not rule out disturbances at other times. This joins published spacecraft operations timing to the cached target measurements without trying a new orbital period. The prior [fixed-window aperture result](tic307815853-fixed-windows.md) and [local background regression](tic307815853-background-control.md) motivated the comparison.

## Fixed spacecraft-phase test

Convert DRN Table 1 UTC starts to TDB, then project each frozen center's delay from its own segment start onto all three other starts. QLP times also include a target barycentric correction not applied to the DRN times. The projection is an approximate schedule comparison. Evaluate the center and offsets ±0.04 d, without optimizing a phase, and retain the 0.16-day quality-zero box and existing sideband scorer.

| Frozen reference | Center BTJD | Delay after start | Raw dimming | Maximum dimming in its three alternate segments, including offsets |
| --- | ---: | ---: | ---: | ---: |
| 181a | 3665.989244 | 2.9540 d | 0.460% | 0.136% |
| 182a | 3679.915773 | 3.2653 d | 0.739% | 0.102% |

All six alternate center windows contained 69 or 70 quality-zero samples. At the 182a delay projected onto 181a, the estimated background decreased by 407.6 native column units, larger than at either selected dip, while raw SAP dimmed by only 0.0785%. Because both delays were identified from the same searched sector, these are targeted descriptive comparisons rather than a calibrated significance test. The target is on camera 1, CCD 4. No other cached Sector 87 target or pixel series was available for a contemporaneous spatial comparison.

## Correction to the November 26 background control

The original script used `ORBITID` as its same-orbit guard. This FITS column is zero for 11,348 of 11,350 rows, and does not reliably partition this table. Re-ran the original model using the four DRN observation intervals converted to TDB and padded by ±0.01 d to allow for the missing target barycentric term. Require evaluation and training samples within the same mapped segment. Both event fits use the same 329 and 439 training cadences and their depths are unchanged. One fewer control qualifies. For the remaining 23, raw versus background-adjusted box-depth RMS is 0.0916% versus 0.1013%, superseding the earlier 24-window figures. The qualitative result is unchanged. This is a correction to the research diagnostic and does not establish a product or production-code bug. Any future orbit-matched analysis of this FITS product should check that the grouping metadata actually vary or use the archive segment boundaries.

## Provenance and next check

The single input is the previous `mast:HLSP/qlp/s0087/0000/0003/0781/5853/hlsp_qlp_tess_ffi_s0087-0000000307815853_tess_v01_llc.fits` cache file, SHA-256 `d0a763d16aa02b427ec6dc973c51925aec07586d86d096b4f458491dd0933a8f`. No new FITS bytes, no local writes or search threshold changes. Source checkout `25512680d1e4a9ac2fda85eb6836367f88bf8e61`. In-process work took 0.659 seconds across these two checks. The published release notes were accessed 2026-11-27. The web tool could not render the release note figures, and no plotted telemetry comparison is claimed. Continue with pixel deficit maps at the two frozen centers and use the strongest background-only window as an additional control. Pixel acquisition remains blocked under the recorded local sandbox limits.

### Reproduce spacecraft-phase check

```python
import hashlib,json,sys,time
from pathlib import Path
import astropy,numpy as np
from astropy.io import fits
from astropy.time import Time
sys.path.insert(0,'legacy')
from screen import scores

start=time.perf_counter()
p,=Path('.audit-cache/shared/fits').glob('*/hlsp_qlp_tess_ffi_s0087-0000000307815853_tess_v01_llc.fits')
segments={
 '181a':('2024-12-18T12:49:33','2024-12-25T00:19:33'),
 '181b':('2024-12-25T05:19:33','2024-12-31T22:33:33'),
 '182a':('2025-01-01T03:35:33','2025-01-07T23:23:32'),
 '182b':('2025-01-08T04:25:32','2025-01-14T10:29:32')}
starts={s:Time(a,scale='utc').tdb.jd-2457000 for s,(a,_) in segments.items()}
ends={s:Time(b,scale='utc').tdb.jd-2457000 for s,(_,b) in segments.items()}
events={'181a':3665.9892441863212,'182a':3679.9157725651703}
d=.16
with fits.open(p,memmap=False) as f:
 h=f[1].header
 assert h['TIMESYS']=='TDB' and h['BJDREFI']==2457000 and h.get('BJDREFF',0)==0
 a=f[1].data
 cols=['SAP_FLUX','DET_FLUX','SAP_BKG']
 keep=(a['QUALITY']==0)&np.isfinite(a['TIME'])
 for k in cols:
  keep &= np.isfinite(a[k])
 b=a[keep]
 t,ix=np.unique(b['TIME'],return_index=True)
 b=b[ix]
 def meas(center):
  c=np.array([center-d/4,center,center+d/4])
  out={}
  for k in cols:
   v=b[k].astype(float)
   if k.endswith('FLUX'):
    v=v/np.median(v)-1
   lookup={r['tc']:r for r in scores(t,v,c,d)}
   out[k]=[lookup.get(u) for u in c]
  return out
 out=[]
 for source,c in events.items():
  aoff=c-starts[source]
  for s in segments:
   center=starts[s]+aoff
   q=a['QUALITY'][abs(a['TIME']-center)<d/2]
   nearest=a[np.nanargmin(abs(a['TIME']-center))]
   out.append(dict(phase_source=source,segment=s,is_original=s==source,btjd=center,
    offset_from_start_days=aoff,days_to_end=ends[s]-center,qlp_cadence=int(nearest['CADENCENO']),
    qlp_orbit_id=int(nearest['ORBITID']),quality={str(int(k)):int(n) for k,n in zip(*np.unique(q,return_counts=True))},
    measurements=meas(center)))
 orbit_ids={str(int(k)):int(n) for k,n in zip(*np.unique(a['ORBITID'],return_counts=True))}
print(json.dumps(dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),n=len(t),
   python=sys.version.split()[0],numpy=np.__version__,astropy=astropy.__version__,seconds=time.perf_counter()-start,
   starts_tdb_without_barycentric=starts,orbit_id_counts=orbit_ids,windows=out),indent=2))
```

```json
{
  "file": "hlsp_qlp_tess_ffi_s0087-0000000307815853_tess_v01_llc.fits",
  "sha256": "d0a763d16aa02b427ec6dc973c51925aec07586d86d096b4f458491dd0933a8f",
  "n": 8437,
  "python": "3.13.15",
  "numpy": "2.5.4",
  "astropy": "8.0.2",
  "seconds": 0.28734040004201233,
  "starts_tdb_without_barycentric": {
    "181a": 3663.0352104576305,
    "181b": 3669.722710459493,
    "182a": 3676.650488239713,
    "182b": 3683.685198890511
  },
  "orbit_id_counts": {
    "0": 11348,
    "181": 1,
    "182": 1
  },
  "windows": [
    {
      "phase_source": "181a",
      "segment": "181a",
      "is_original": true,
      "btjd": 3665.9892441863212,
      "offset_from_start_days": 2.9540337286907743,
      "days_to_end": 3.5251329401489784,
      "qlp_cadence": 1053442,
      "qlp_orbit_id": 0,
      "quality": {
        "0": 69
      },
      "measurements": {
        "SAP_FLUX": [
          {
            "tc": 3665.9492441863213,
            "duration": 0.16,
            "depth": 0.004534016617634322,
            "error": 0.001422315905366568,
            "score": 3.187770452771382,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3665.9892441863212,
            "duration": 0.16,
            "depth": 0.004596368935568644,
            "error": 0.0015012291373165774,
            "score": 3.0617370935023143,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3666.029244186321,
            "duration": 0.16,
            "depth": 0.002960130241406069,
            "error": 0.0014775173430404548,
            "score": 2.0034487279280753,
            "points": 70,
            "cadence": 0.0023148116069933167
          }
        ],
        "DET_FLUX": [
          {
            "tc": 3665.9492441863213,
            "duration": 0.16,
            "depth": 0.004346952298988391,
            "error": 0.0014320533799561269,
            "score": 3.0354680627349007,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3665.9892441863212,
            "duration": 0.16,
            "depth": 0.004453110949046289,
            "error": 0.0015152714403616862,
            "score": 2.938820616841665,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3666.029244186321,
            "duration": 0.16,
            "depth": 0.0027150384635381275,
            "error": 0.001480365756467612,
            "score": 1.8340321989186246,
            "points": 70,
            "cadence": 0.0023148116069933167
          }
        ],
        "SAP_BKG": [
          {
            "tc": 3665.9492441863213,
            "duration": 0.16,
            "depth": 298.02828511422103,
            "error": 102.03710189384478,
            "score": 2.920783514846173,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3665.9892441863212,
            "duration": 0.16,
            "depth": 286.78723627102136,
            "error": 93.53169955059496,
            "score": 3.066203625604888,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3666.029244186321,
            "duration": 0.16,
            "depth": 280.9081442252118,
            "error": 88.32305286499894,
            "score": 3.180462349445479,
            "points": 70,
            "cadence": 0.0023148116069933167
          }
        ]
      }
    },
    {
      "phase_source": "181a",
      "segment": "181b",
      "is_original": false,
      "btjd": 3672.676744188184,
      "offset_from_start_days": 2.9540337286907743,
      "days_to_end": 3.76402182923448,
      "qlp_cadence": 1056331,
      "qlp_orbit_id": 0,
      "quality": {
        "0": 69
      },
      "measurements": {
        "SAP_FLUX": [
          {
            "tc": 3672.636744188184,
            "duration": 0.16,
            "depth": 0.00027873535377822677,
            "error": 0.0012871662215451641,
            "score": 0.21654961815547186,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3672.676744188184,
            "duration": 0.16,
            "depth": 0.00009794166865873498,
            "error": 0.0012272379384324823,
            "score": 0.07980658484517943,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3672.716744188184,
            "duration": 0.16,
            "depth": 0.00003225320923344946,
            "error": 0.001236353602540023,
            "score": 0.026087366241491876,
            "points": 70,
            "cadence": 0.0023148116069933167
          }
        ],
        "DET_FLUX": [
          {
            "tc": 3672.636744188184,
            "duration": 0.16,
            "depth": 0.00043759415229717566,
            "error": 0.0012982052380437789,
            "score": 0.3370762491734907,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3672.676744188184,
            "duration": 0.16,
            "depth": 0.0002658378757988951,
            "error": 0.001231582011732196,
            "score": 0.21585072960346288,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3672.716744188184,
            "duration": 0.16,
            "depth": 0.00023439868348595292,
            "error": 0.001239281760219287,
            "score": 0.1891407515305291,
            "points": 70,
            "cadence": 0.0023148116069933167
          }
        ],
        "SAP_BKG": [
          {
            "tc": 3672.636744188184,
            "duration": 0.16,
            "depth": 13.302324759381548,
            "error": 11.554100125361831,
            "score": 1.1513077275643713,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3672.676744188184,
            "duration": 0.16,
            "depth": 17.472358784286172,
            "error": 10.254018112015281,
            "score": 1.7039524012360292,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3672.716744188184,
            "duration": 0.16,
            "depth": 22.06635521577573,
            "error": 9.92516305741482,
            "score": 2.22327382312279,
            "points": 70,
            "cadence": 0.0023148116069933167
          }
        ]
      }
    },
    {
      "phase_source": "181a",
      "segment": "182a",
      "is_original": false,
      "btjd": 3679.604521968404,
      "offset_from_start_days": 2.9540337286907743,
      "days_to_end": 3.8709546998125006,
      "qlp_cadence": 1059324,
      "qlp_orbit_id": 0,
      "quality": {
        "0": 69
      },
      "measurements": {
        "SAP_FLUX": [
          {
            "tc": 3679.564521968404,
            "duration": 0.16,
            "depth": 0.0013618919056503226,
            "error": 0.001500428260226115,
            "score": 0.9076687914722994,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.604521968404,
            "duration": 0.16,
            "depth": -0.0015731261019054955,
            "error": 0.001469122739330699,
            "score": -1.0707928342475852,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.6445219684038,
            "duration": 0.16,
            "depth": -0.00417017364152082,
            "error": 0.0014354644843702958,
            "score": -2.9051040181953205,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ],
        "DET_FLUX": [
          {
            "tc": 3679.564521968404,
            "duration": 0.16,
            "depth": 0.0006955376389025003,
            "error": 0.0014765848370848905,
            "score": 0.471044820069836,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.604521968404,
            "duration": 0.16,
            "depth": -0.002294245092786911,
            "error": 0.001451081659446275,
            "score": -1.5810585695517525,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.6445219684038,
            "duration": 0.16,
            "depth": -0.004866867332870379,
            "error": 0.0014212381495574378,
            "score": -3.4243855151129194,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ],
        "SAP_BKG": [
          {
            "tc": 3679.564521968404,
            "duration": 0.16,
            "depth": 321.6382261386616,
            "error": 167.20570057883725,
            "score": 1.9236080171023213,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.604521968404,
            "duration": 0.16,
            "depth": 96.02258388830342,
            "error": 138.78113818315555,
            "score": 0.6918993830528916,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.6445219684038,
            "duration": 0.16,
            "depth": -3.85949018570318,
            "error": 102.90806299795197,
            "score": -0.0375042545089979,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ]
      }
    },
    {
      "phase_source": "181a",
      "segment": "182b",
      "is_original": false,
      "btjd": 3686.6392326192017,
      "offset_from_start_days": 2.9540337286907743,
      "days_to_end": 3.29874405078408,
      "qlp_cadence": 1062363,
      "qlp_orbit_id": 0,
      "quality": {
        "0": 69
      },
      "measurements": {
        "SAP_FLUX": [
          {
            "tc": 3686.5992326192018,
            "duration": 0.16,
            "depth": 0.0008744716024661342,
            "error": 0.0010821076212824989,
            "score": 0.8081188832490827,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.6392326192017,
            "duration": 0.16,
            "depth": 0.00021334105557839662,
            "error": 0.0011153270686149028,
            "score": 0.1912811601025156,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.6792326192017,
            "duration": 0.16,
            "depth": -0.0002806290931987813,
            "error": 0.001114099765050055,
            "score": -0.251888656655603,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ],
        "DET_FLUX": [
          {
            "tc": 3686.5992326192018,
            "duration": 0.16,
            "depth": 0.0008831261751080155,
            "error": 0.0010483015398506305,
            "score": 0.8424352550639672,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.6392326192017,
            "duration": 0.16,
            "depth": 0.000254841135586412,
            "error": 0.0010947348078691951,
            "score": 0.23278800834189042,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.6792326192017,
            "duration": 0.16,
            "depth": -0.0001950766309127488,
            "error": 0.0011067373543910543,
            "score": -0.1762628053880799,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ],
        "SAP_BKG": [
          {
            "tc": 3686.5992326192018,
            "duration": 0.16,
            "depth": -22.59474954336406,
            "error": 13.020873287776276,
            "score": -1.735271440247832,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.6392326192017,
            "duration": 0.16,
            "depth": -21.356885044726482,
            "error": 12.151615901107379,
            "score": -1.7575345714128625,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.6792326192017,
            "duration": 0.16,
            "depth": -28.767312983825263,
            "error": 12.264406634541809,
            "score": -2.3455935408080992,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ]
      }
    },
    {
      "phase_source": "182a",
      "segment": "181a",
      "is_original": false,
      "btjd": 3666.3004947830877,
      "offset_from_start_days": 3.265284325457287,
      "days_to_end": 3.2138823433824655,
      "qlp_cadence": 1053577,
      "qlp_orbit_id": 0,
      "quality": {
        "0": 69
      },
      "measurements": {
        "SAP_FLUX": [
          {
            "tc": 3666.260494783088,
            "duration": 0.16,
            "depth": -0.0007665967823830971,
            "error": 0.0015394006085594793,
            "score": -0.49798394136043206,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3666.3004947830877,
            "duration": 0.16,
            "depth": 0.0007846928383609467,
            "error": 0.001440530225816722,
            "score": 0.544725007707532,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3666.3404947830877,
            "duration": 0.16,
            "depth": 0.0006573275532573584,
            "error": 0.0013987096665650889,
            "score": 0.4699528207820316,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ],
        "DET_FLUX": [
          {
            "tc": 3666.260494783088,
            "duration": 0.16,
            "depth": -0.000759360960580873,
            "error": 0.001550721439439311,
            "score": -0.4896823770331263,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3666.3004947830877,
            "duration": 0.16,
            "depth": 0.0008228388615143124,
            "error": 0.001452227709670626,
            "score": 0.5666045731223082,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3666.3404947830877,
            "duration": 0.16,
            "depth": 0.0005543234320632843,
            "error": 0.0014087685588478871,
            "score": 0.3934808372757968,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ],
        "SAP_BKG": [
          {
            "tc": 3666.260494783088,
            "duration": 0.16,
            "depth": 395.34525233747036,
            "error": 65.56453363207204,
            "score": 6.029864477585185,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3666.3004947830877,
            "duration": 0.16,
            "depth": 407.5853998412583,
            "error": 61.415022440231866,
            "score": 6.6365749558736065,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3666.3404947830877,
            "duration": 0.16,
            "depth": 399.7491098783868,
            "error": 59.99362455566773,
            "score": 6.663193178259499,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ]
      }
    },
    {
      "phase_source": "182a",
      "segment": "181b",
      "is_original": false,
      "btjd": 3672.9879947849504,
      "offset_from_start_days": 3.265284325457287,
      "days_to_end": 3.452771232467967,
      "qlp_cadence": 1056466,
      "qlp_orbit_id": 0,
      "quality": {
        "0": 69
      },
      "measurements": {
        "SAP_FLUX": [
          {
            "tc": 3672.9479947849504,
            "duration": 0.16,
            "depth": -0.0008753307879443388,
            "error": 0.0011740143034142958,
            "score": -0.7455878394314971,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3672.9879947849504,
            "duration": 0.16,
            "depth": -0.0023795938000402365,
            "error": 0.0011506420838768704,
            "score": -2.068057333712883,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3673.0279947849504,
            "duration": 0.16,
            "depth": -0.002308025069195891,
            "error": 0.0010677601954858509,
            "score": -2.1615575097793345,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ],
        "DET_FLUX": [
          {
            "tc": 3672.9479947849504,
            "duration": 0.16,
            "depth": -0.0011216441698883568,
            "error": 0.0011705385460499724,
            "score": -0.9582291618446812,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3672.9879947849504,
            "duration": 0.16,
            "depth": -0.0024561006210095473,
            "error": 0.001138458901682127,
            "score": -2.1573906773275193,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3673.0279947849504,
            "duration": 0.16,
            "depth": -0.0023127514960763077,
            "error": 0.0010479495358569387,
            "score": -2.2069302165252678,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ],
        "SAP_BKG": [
          {
            "tc": 3672.9479947849504,
            "duration": 0.16,
            "depth": -3.0940727581851206,
            "error": 9.05057048362671,
            "score": -0.3418649425229685,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3672.9879947849504,
            "duration": 0.16,
            "depth": 6.389496566687484,
            "error": 9.430276491954414,
            "score": 0.6775513498611394,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3673.0279947849504,
            "duration": 0.16,
            "depth": 0.442956616427324,
            "error": 9.310811344741103,
            "score": 0.047574437933114506,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ]
      }
    },
    {
      "phase_source": "182a",
      "segment": "182a",
      "is_original": true,
      "btjd": 3679.9157725651703,
      "offset_from_start_days": 3.265284325457287,
      "days_to_end": 3.5597041030459877,
      "qlp_cadence": 1059458,
      "qlp_orbit_id": 0,
      "quality": {
        "0": 70
      },
      "measurements": {
        "SAP_FLUX": [
          {
            "tc": 3679.8757725651703,
            "duration": 0.16,
            "depth": 0.005666668349279871,
            "error": 0.0012608666361834597,
            "score": 4.4942646483472775,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.9157725651703,
            "duration": 0.16,
            "depth": 0.007385813889099514,
            "error": 0.0012778018646507757,
            "score": 5.780093215874327,
            "points": 70,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.9557725651703,
            "duration": 0.16,
            "depth": 0.006177100661063738,
            "error": 0.001291788649474712,
            "score": 4.781819892577296,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ],
        "DET_FLUX": [
          {
            "tc": 3679.8757725651703,
            "duration": 0.16,
            "depth": 0.00508345122290972,
            "error": 0.001289007440282897,
            "score": 3.943694244149639,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.9157725651703,
            "duration": 0.16,
            "depth": 0.0069020900684554345,
            "error": 0.0012903401726936813,
            "score": 5.349046875016537,
            "points": 70,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.9557725651703,
            "duration": 0.16,
            "depth": 0.005783196938114946,
            "error": 0.0013071797477531666,
            "score": 4.4241788078918285,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ],
        "SAP_BKG": [
          {
            "tc": 3679.8757725651703,
            "duration": 0.16,
            "depth": 253.46129926003323,
            "error": 73.626453997664,
            "score": 3.442530306675845,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.9157725651703,
            "duration": 0.16,
            "depth": 230.95847445772188,
            "error": 71.78893157977491,
            "score": 3.2171877944871072,
            "points": 70,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3679.9557725651703,
            "duration": 0.16,
            "depth": 218.2554713856921,
            "error": 73.71009171261971,
            "score": 2.960998505287779,
            "points": 69,
            "cadence": 0.0023148116069933167
          }
        ]
      }
    },
    {
      "phase_source": "182a",
      "segment": "182b",
      "is_original": false,
      "btjd": 3686.9504832159682,
      "offset_from_start_days": 3.265284325457287,
      "days_to_end": 2.9874934540175673,
      "qlp_cadence": 1062498,
      "qlp_orbit_id": 0,
      "quality": {
        "0": 70
      },
      "measurements": {
        "SAP_FLUX": [
          {
            "tc": 3686.9104832159683,
            "duration": 0.16,
            "depth": -0.00019174012410530863,
            "error": 0.0012661821857379227,
            "score": -0.15143170253462676,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.9504832159682,
            "duration": 0.16,
            "depth": -0.0003177464935004079,
            "error": 0.0012239067005424967,
            "score": -0.25961659770272255,
            "points": 70,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.990483215968,
            "duration": 0.16,
            "depth": 0.0010184408707035481,
            "error": 0.0012723058310799252,
            "score": 0.8004686026150661,
            "points": 68,
            "cadence": 0.0023148116069933167
          }
        ],
        "DET_FLUX": [
          {
            "tc": 3686.9104832159683,
            "duration": 0.16,
            "depth": -0.000043467689275918456,
            "error": 0.001253139564568564,
            "score": -0.034687029685224004,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.9504832159682,
            "duration": 0.16,
            "depth": -0.00014200567320002893,
            "error": 0.001209933443573419,
            "score": -0.11736651627764681,
            "points": 70,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.990483215968,
            "duration": 0.16,
            "depth": 0.0011859014464704431,
            "error": 0.001235924725827217,
            "score": 0.9595256261878792,
            "points": 68,
            "cadence": 0.0023148116069933167
          }
        ],
        "SAP_BKG": [
          {
            "tc": 3686.9104832159683,
            "duration": 0.16,
            "depth": 37.63894899281553,
            "error": 13.720961115827738,
            "score": 2.743171464089154,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.9504832159682,
            "duration": 0.16,
            "depth": 19.908330262748677,
            "error": 13.703673396235065,
            "score": 1.4527732591918219,
            "points": 70,
            "cadence": 0.0023148116069933167
          },
          {
            "tc": 3686.990483215968,
            "duration": 0.16,
            "depth": -11.156092212628664,
            "error": 13.566196358974727,
            "score": -0.8223448870580694,
            "points": 68,
            "cadence": 0.0023148116069933167
          }
        ]
      }
    }
  ]
}
```

### Reproduce corrected background comparison

```python
import hashlib,json,sys,time
from pathlib import Path
import astropy,numpy as np
from astropy.io import fits
sys.path.insert(0,'legacy')
from screen import scores

start=time.perf_counter()
path,=Path('.audit-cache/shared/fits').glob('*/hlsp_qlp_tess_ffi_s0087-0000000307815853_tess_v01_llc.fits')
with fits.open(path,memmap=False) as hdus:
    h=hdus[1].header
    assert h['TIMESYS']=='TDB' and h['BJDREFI']==2457000 and h.get('BJDREFF',0)==0
    raw=hdus[1].data
    cols=['TIME','SAP_FLUX','SAP_BKG','ORBITID']
    keep=(raw['QUALITY']==0)
    for name in cols:
        keep &= np.isfinite(raw[name])
    a=raw[keep]
t,idx=np.unique(a['TIME'],return_index=True)
a=a[idx]
y=a['SAP_FLUX'].astype(float)
y=y/np.median(y)-1
bg=a['SAP_BKG'].astype(float)
from astropy.time import Time
utc_pairs=[('2024-12-18T12:49:33','2024-12-25T00:19:33'),
 ('2024-12-25T05:19:33','2024-12-31T22:33:33'),
 ('2025-01-01T03:35:33','2025-01-07T23:23:32'),
 ('2025-01-08T04:25:32','2025-01-14T10:29:32')]
orb=np.full(t.size,-1)
for i,(lo,hi) in enumerate(utc_pairs):
    start_tdb,end_tdb=Time([lo,hi],scale='utc').tdb.jd-2457000
    orb[(t>=start_tdb-.01)&(t<=end_tdb+.01)]=i
events=np.array([3665.9892441863212,3679.9157725651703])
d=.16
grid=np.arange(t[0]+2*d,t[-1]-2*d,4*d)
grid=grid[np.min(abs(grid[:,None]-events),axis=1)>4*d]
event_free=np.min(abs(t[:,None]-events),axis=1)>2*d
results=[]
for c in np.r_[events,grid]:
    event=bool(c in events)
    baseline=scores(t,y,np.array([c]),d)
    result=dict(btjd=float(c),event=event,status='invalid_sidebands')
    results.append(result)
    if not baseline:
        continue
    delta=t-c
    inside=(abs(delta)<d/2)
    eval_orbits=np.unique(orb[abs(delta)<2*d])
    result['status']='insufficient_local_training'
    if len(eval_orbits)!=1 or eval_orbits[0]<0:
        continue
    train=(abs(delta)<=6*d)&(abs(delta)>2*d)&(orb==eval_orbits[0])&event_free
    if min(np.count_nonzero(train&(delta<0)),np.count_nonzero(train&(delta>0)))<20:
        continue
    bscale=bg[train].std()
    if not bscale>0:
        continue
    bs=(bg-bg[train].mean())/bscale
    x=np.column_stack([np.ones(len(t)),delta,bs])
    coef=np.linalg.lstsq(x[train],y[train],rcond=None)[0]
    b=scores(t,bs,np.array([c]),d)[0]['depth']
    depth=baseline[0]['depth']
    depth_adjusted=scores(t,y-coef[2]*bs,np.array([c]),d)[0]['depth']
    result.update(status='measured',train=int(train.sum()),box=int(inside.sum()),
        raw_depth=depth,adjusted_depth=depth_adjusted,predicted_background_depth=float(coef[2]*b),
        normalized_flux_per_native_bkg=float(coef[2]/bscale),
        bkg_z_max_abs_in_box=float(max(abs(bs[inside]))),
        in_box_bkg_outside_training_range=int(np.count_nonzero(inside&((bs<bs[train].min())|(bs>bs[train].max())))))
controls=[r for r in results if not r['event'] and r['status']=='measured']
print(json.dumps(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    n=len(t),windows=results,controls=len(controls),
    control_rms_raw=float(np.sqrt(np.mean([r['raw_depth']**2 for r in controls]))),
    control_rms_adjusted=float(np.sqrt(np.mean([r['adjusted_depth']**2 for r in controls]))),
    python=sys.version.split()[0],numpy=np.__version__,astropy=astropy.__version__,
    seconds=time.perf_counter()-start),indent=2))
```

```json
{
  "file": "hlsp_qlp_tess_ffi_s0087-0000000307815853_tess_v01_llc.fits",
  "sha256": "d0a763d16aa02b427ec6dc973c51925aec07586d86d096b4f458491dd0933a8f",
  "n": 8437,
  "windows": [
    {
      "btjd": 3665.9892441863212,
      "event": true,
      "status": "measured",
      "train": 329,
      "box": 69,
      "raw_depth": 0.004596368935568644,
      "adjusted_depth": 0.004574322272399422,
      "predicted_background_depth": 0.000022046663169221864,
      "normalized_flux_per_native_bkg": 7.687463171613002e-8,
      "bkg_z_max_abs_in_box": 0.8905665282505033,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3679.9157725651703,
      "event": true,
      "status": "measured",
      "train": 439,
      "box": 70,
      "raw_depth": 0.007385813889099514,
      "adjusted_depth": 0.007087596297897123,
      "predicted_background_depth": 0.0002982175912022844,
      "normalized_flux_per_native_bkg": 0.000001291217358023536,
      "bkg_z_max_abs_in_box": 0.21937081003792933,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3667.0492858755238,
      "event": false,
      "status": "measured",
      "train": 457,
      "box": 69,
      "raw_depth": -0.0006655374250797664,
      "adjusted_depth": -0.00025563193288724237,
      "predicted_background_depth": -0.0004099054921925472,
      "normalized_flux_per_native_bkg": -0.000018991697027364844,
      "bkg_z_max_abs_in_box": 1.2270938512230412,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3667.6892858755236,
      "event": false,
      "status": "measured",
      "train": 554,
      "box": 69,
      "raw_depth": -0.0011922341194355538,
      "adjusted_depth": -0.0014901434849752702,
      "predicted_background_depth": 0.0002979093655409296,
      "normalized_flux_per_native_bkg": -0.00005064940637940046,
      "bkg_z_max_abs_in_box": 2.8751725634982233,
      "in_box_bkg_outside_training_range": 1
    },
    {
      "btjd": 3668.3292858755235,
      "event": false,
      "status": "measured",
      "train": 552,
      "box": 69,
      "raw_depth": -0.00005115762125199215,
      "adjusted_depth": -0.0007915023022894508,
      "predicted_background_depth": 0.0007403446810371143,
      "normalized_flux_per_native_bkg": -0.000059115572907857384,
      "bkg_z_max_abs_in_box": 3.0599153821879974,
      "in_box_bkg_outside_training_range": 2
    },
    {
      "btjd": 3668.9692858755234,
      "event": false,
      "status": "measured",
      "train": 377,
      "box": 69,
      "raw_depth": -0.0007980532147605773,
      "adjusted_depth": -0.001291485205651634,
      "predicted_background_depth": 0.0004934319908909757,
      "normalized_flux_per_native_bkg": -0.00005770623426733741,
      "bkg_z_max_abs_in_box": 2.6681499209297495,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3669.6092858755233,
      "event": false,
      "status": "invalid_sidebands"
    },
    {
      "btjd": 3670.249285875523,
      "event": false,
      "status": "measured",
      "train": 361,
      "box": 69,
      "raw_depth": -0.0010933084302461576,
      "adjusted_depth": -0.0007724325934905943,
      "predicted_background_depth": -0.0003208758367550008,
      "normalized_flux_per_native_bkg": -0.0000426082136475119,
      "bkg_z_max_abs_in_box": 2.0965851690796753,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3670.889285875523,
      "event": false,
      "status": "measured",
      "train": 551,
      "box": 69,
      "raw_depth": -0.0007889840091456633,
      "adjusted_depth": -0.0018669323993235742,
      "predicted_background_depth": 0.0010779483901788287,
      "normalized_flux_per_native_bkg": -0.000043563630360421226,
      "bkg_z_max_abs_in_box": 2.6012918676370034,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3671.529285875523,
      "event": false,
      "status": "measured",
      "train": 554,
      "box": 70,
      "raw_depth": -0.0007919404090318344,
      "adjusted_depth": 0.00004964128995673181,
      "predicted_background_depth": -0.0008415816989881815,
      "normalized_flux_per_native_bkg": -0.00004413794348596567,
      "bkg_z_max_abs_in_box": 2.4787027470544003,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3672.1692858755227,
      "event": false,
      "status": "measured",
      "train": 552,
      "box": 69,
      "raw_depth": -0.0003938756561572222,
      "adjusted_depth": 0.00007444394624207015,
      "predicted_background_depth": -0.0004683196023994726,
      "normalized_flux_per_native_bkg": -0.000041515088509692364,
      "bkg_z_max_abs_in_box": 2.4866075567244983,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3672.8092858755226,
      "event": false,
      "status": "measured",
      "train": 554,
      "box": 70,
      "raw_depth": 0.000676819381454195,
      "adjusted_depth": 0.0010343220390356192,
      "predicted_background_depth": -0.0003575026575811766,
      "normalized_flux_per_native_bkg": -0.000042279466916610417,
      "bkg_z_max_abs_in_box": 1.7932384977787028,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3673.4492858755225,
      "event": false,
      "status": "measured",
      "train": 552,
      "box": 69,
      "raw_depth": -0.0020166299788425676,
      "adjusted_depth": -0.001412393336781277,
      "predicted_background_depth": -0.0006042366420618884,
      "normalized_flux_per_native_bkg": -0.00004354005866020814,
      "bkg_z_max_abs_in_box": 2.3552194275953378,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3674.0892858755224,
      "event": false,
      "status": "measured",
      "train": 554,
      "box": 70,
      "raw_depth": 0.00022042005268731873,
      "adjusted_depth": 0.001045885537915541,
      "predicted_background_depth": -0.0008254654852287511,
      "normalized_flux_per_native_bkg": -0.000046846428816299884,
      "bkg_z_max_abs_in_box": 2.790746189945513,
      "in_box_bkg_outside_training_range": 1
    },
    {
      "btjd": 3674.7292858755222,
      "event": false,
      "status": "measured",
      "train": 537,
      "box": 69,
      "raw_depth": 0.00019408477608081134,
      "adjusted_depth": 0.00046254412922489517,
      "predicted_background_depth": -0.0002684593531444914,
      "normalized_flux_per_native_bkg": -0.0000476371074944933,
      "bkg_z_max_abs_in_box": 2.5116675781249116,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3675.369285875522,
      "event": false,
      "status": "measured",
      "train": 485,
      "box": 67,
      "raw_depth": 0.0005656371197355087,
      "adjusted_depth": 0.0005580704327047409,
      "predicted_background_depth": 0.000007566687030806984,
      "normalized_flux_per_native_bkg": 6.314143029525638e-7,
      "bkg_z_max_abs_in_box": 0.5601278949814459,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3676.009285875522,
      "event": false,
      "status": "insufficient_local_training"
    },
    {
      "btjd": 3676.649285875522,
      "event": false,
      "status": "invalid_sidebands"
    },
    {
      "btjd": 3677.2892858755217,
      "event": false,
      "status": "invalid_sidebands"
    },
    {
      "btjd": 3677.9292858755216,
      "event": false,
      "status": "invalid_sidebands"
    },
    {
      "btjd": 3678.5692858755215,
      "event": false,
      "status": "invalid_sidebands"
    },
    {
      "btjd": 3679.2092858755213,
      "event": false,
      "status": "invalid_sidebands"
    },
    {
      "btjd": 3681.129285875521,
      "event": false,
      "status": "measured",
      "train": 501,
      "box": 69,
      "raw_depth": 0.0006463468340721264,
      "adjusted_depth": 0.0006878230359573211,
      "predicted_background_depth": -0.00004147620188506419,
      "normalized_flux_per_native_bkg": -0.000005180844745437442,
      "bkg_z_max_abs_in_box": 1.3489231397929453,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3681.769285875521,
      "event": false,
      "status": "measured",
      "train": 514,
      "box": 65,
      "raw_depth": 0.0010127616098617696,
      "adjusted_depth": 0.0003724696680542641,
      "predicted_background_depth": 0.00064029194180633,
      "normalized_flux_per_native_bkg": -0.000028428307697298708,
      "bkg_z_max_abs_in_box": 2.5754580755376923,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3682.4092858755207,
      "event": false,
      "status": "measured",
      "train": 482,
      "box": 60,
      "raw_depth": 0.0010836334039453963,
      "adjusted_depth": 0.0010191281178705365,
      "predicted_background_depth": 0.00006450528607358434,
      "normalized_flux_per_native_bkg": -0.000026261609061694015,
      "bkg_z_max_abs_in_box": 2.1393570308515444,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3683.0492858755206,
      "event": false,
      "status": "measured",
      "train": 267,
      "box": 55,
      "raw_depth": -0.00069333549876186,
      "adjusted_depth": -0.0019851729399095242,
      "predicted_background_depth": 0.001291837441141895,
      "normalized_flux_per_native_bkg": -0.000037953982812664396,
      "bkg_z_max_abs_in_box": 2.7754862970655028,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3683.6892858755205,
      "event": false,
      "status": "invalid_sidebands"
    },
    {
      "btjd": 3684.3292858755203,
      "event": false,
      "status": "measured",
      "train": 410,
      "box": 69,
      "raw_depth": 0.000009804762052194327,
      "adjusted_depth": -0.00015981872002524856,
      "predicted_background_depth": 0.00016962348207460892,
      "normalized_flux_per_native_bkg": -0.000037388300564279244,
      "bkg_z_max_abs_in_box": 2.0357867619251913,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3684.96928587552,
      "event": false,
      "status": "measured",
      "train": 554,
      "box": 69,
      "raw_depth": 0.0007961414280171893,
      "adjusted_depth": 0.00024946435723980644,
      "predicted_background_depth": 0.0005466770707798172,
      "normalized_flux_per_native_bkg": -0.00004499637549149866,
      "bkg_z_max_abs_in_box": 3.078995334958061,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3685.60928587552,
      "event": false,
      "status": "measured",
      "train": 551,
      "box": 69,
      "raw_depth": -0.000630917713328895,
      "adjusted_depth": 0.000015983324132263794,
      "predicted_background_depth": -0.0006469010374610723,
      "normalized_flux_per_native_bkg": -0.00003409841928992503,
      "bkg_z_max_abs_in_box": 2.2694176151379493,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3686.24928587552,
      "event": false,
      "status": "measured",
      "train": 553,
      "box": 69,
      "raw_depth": -0.0020377171518735663,
      "adjusted_depth": -0.0019041359805536864,
      "predicted_background_depth": -0.00013358117131675578,
      "normalized_flux_per_native_bkg": -0.0000368858033096793,
      "bkg_z_max_abs_in_box": 2.628853777743507,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3686.88928587552,
      "event": false,
      "status": "measured",
      "train": 534,
      "box": 69,
      "raw_depth": -0.00034967393784987935,
      "adjusted_depth": 0.0008500486983912623,
      "predicted_background_depth": -0.001199722636240471,
      "normalized_flux_per_native_bkg": -0.000030120482845449875,
      "bkg_z_max_abs_in_box": 3.3211129930820458,
      "in_box_bkg_outside_training_range": 1
    },
    {
      "btjd": 3687.5292858755197,
      "event": false,
      "status": "measured",
      "train": 463,
      "box": 67,
      "raw_depth": -0.0009178863370613508,
      "adjusted_depth": -0.0004893704936653458,
      "predicted_background_depth": -0.00042851584339343953,
      "normalized_flux_per_native_bkg": -0.0000179373101965142,
      "bkg_z_max_abs_in_box": 1.5968578861573992,
      "in_box_bkg_outside_training_range": 0
    },
    {
      "btjd": 3688.1692858755196,
      "event": false,
      "status": "insufficient_local_training"
    }
  ],
  "controls": 23,
  "control_rms_raw": 0.0009156269853127205,
  "control_rms_adjusted": 0.0010134212672659802,
  "python": "3.13.15",
  "numpy": "2.5.4",
  "astropy": "8.0.2",
  "seconds": 0.37190860003465787
}
```

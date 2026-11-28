# TIC 307815853: conditional T16 anchor timing and coverage

**Disposition: unresolved.** Anchoring the period to the catalog event and the strongest S87 trial window moves the 13.92653 d trial to 13.93571 d under the nearest-cycle assumption. None of the three resulting trials has a new supported cycle or a sensitive non-detection in the cached data. The November 13 search did not independently establish either timing uncertainty or a period.

## Fixed conditional test

Compare the T16 catalog center BTJD 1478.07283 and the strongest earlier QLP trial center BTJD 3679.9157725651703. For the existing 13.926528 d trial and its one-half and one-third harmonics, take the nearest integer number of cycles spanning those two centers and divide their gap by that count. Those periods condition on both adopted anchors and are not period measurements. The S87 center was itself inherited from the earlier grid. Use a 0.16 d box and freeze the injection sensitivity depth to half the original 0.69021% event measurement. Use the established score threshold 5 and the same three timing samples and conservative non-detection rule as the [previous fixed-trial check](../tic140920261-unused-sectors/check.py). Preserve all insufficient cadence or sensitivity outcomes.

| Trial period | Anchor gap in cycles | Checked windows | Nominal support | Sensitive misses |
| ---: | ---: | ---: | --- | ---: |
| 13.935714826 d | 158 | 4 | S87 defining anchor only | 0 |
| 6.967857413 d | 316 | 10 | S87 defining anchor only | 0 |
| 4.645238275 d | 474 | 15 | S87 defining anchor only | 0 |

The S6 catalog time still has insufficient QLP coverage. The earlier S45 trial center at 2537.94045 was unmeasurable. All three T16-anchored trials move the shared S45 center to 2537.18716, where a measurement exists but it is insufficiently sensitive at the frozen 0.34510% injection depth. The earlier S87 interval shifts to 3665.98006 and has nominal score 2.66, still below the support threshold. None of these windows provides new independent evidence. The conditional construction also shows why a 1.45 d phase discrepancy over 158 cycles alone cannot rule out an approximate 13.9-day period. The required period change is about 13.23 minutes. Alternative counts, timing variation, blends and instrumental effects remain open. See the [spacecraft-phase comparison](tic307815853-orbit-control.md) for the evidence against one simple instrumental explanation.

## Next data

The saved MAST inventory identifies one public **T16 Sector 6 FITS** by archive URI, camera/CCD and Gaia identifier, plus unused QLP Sectors 33, 43, 44, 71 and 72. Its MAST position is approximately 0.56 arcsec from the QLP position and its filename's numerical Gaia identifier matches the frozen selection. The [MAST T16 description](https://archive.stsci.edu/hlsp/t16) uses differing DR2 and DR3 labels on the page. Preserve the actual FITS identity and coordinates before assigning a catalog release. Fetch only the T16 S6 product first when permitted. Verify its FITS time, magnitude, quality and aperture conventions before measuring the original event. T16's difference-image extraction is an alternate reduction of the same exposures. The original [T16 publication](https://arxiv.org/html/2604.18579v1) established the catalog signal. There is no novelty claim. Preserve the two S87 pixel windows as the primary source-location test.

## Reproduction

Run from the cache checkout at `25512680d1e4a9ac2fda85eb6836367f88bf8e61`. Three original QLP FITS verified, no new astronomical bytes or local writes. In-process runtime 0.072 s. The archive execution limitation from the [November 21 checkpoint](../tic140920261-unused-sectors/2026-11-21-checkpoint.md) remained in force. The script emits the full result; the run ledger below retains the per-window states, measurements and injection scores.

```python
import hashlib,json,sys,time
from pathlib import Path
import astropy,numpy as np
sys.path.insert(0,'legacy')
from screen import read_curve
from search import cadence_trials,MIN_SCORE

start=time.perf_counter()
root=Path('.audit-cache/shared')
inv=root/'inventory/307815853.json'
report=Path('runs/worker-order-pilot-extension-3b8bef46/report.json')
old=next(t for t in json.loads(report.read_bytes())['targets'] if t['tic']=='307815853')
inventory=json.loads(inv.read_bytes())['value']['products']
curves=[]
receipts=[]
for receipt in old['provenance']:
    s=receipt['sequence_number']
    path,=root.glob(f'fits/*/hlsp_qlp_tess_ffi_s{s:04d}-0000000307815853_tess_v01_llc.fits')
    assert hashlib.sha256(path.read_bytes()).hexdigest()==receipt['sha256']
    t,y,_=read_curve(receipt,path=path)
    curves.append(dict(sector=s,time=t,flux=y))
    receipts.append({k:receipt[k] for k in ('sequence_number','uri','sha256','bytes')})
earlier=1478.07283
later=3679.9157725651703
original=old['evidence']['trial_period_days']
depth=old['evidence']['depth_fraction']/2
out=[]
for harmonic in (1,2,3):
    cycles=round((later-earlier)/(original/harmonic))
    period=(later-earlier)/cycles
    model=dict(period_days=period,epoch_btjd=earlier,duration_days=.16)
    events=cadence_trials(curves,model)['events']
    for e in events:
        ms=e['timing_samples']
        m=e['measurement']
        injection=[(k['depth']+depth)/k['error'] if k else None for k in ms]
        sensitive=all(k and j>=MIN_SCORE and max(0,k['depth'])+MIN_SCORE*k['error']<depth for k,j in zip(ms,injection))
        e['test_state']=('supported' if m and m['score']>=MIN_SCORE else
            'sensitive_non_detection' if sensitive else
            'insufficient_sensitivity' if m else 'insufficient_cadence')
        e['test_scores']=injection
        e['is_anchor']=e['cycle'] in (0,cycles)
    out.append(dict(harmonic=harmonic,base_trial=original/harmonic,cycles_between_anchors=cycles,**model,events=events))
print(json.dumps(dict(tic='307815853',input_receipts=receipts,
   inventory_sha256=hashlib.sha256(inv.read_bytes()).hexdigest(),
   source_report_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),
   unprocessed_t16=[v for v in inventory if v['provenance_name']=='T16'],
   unused_qlp_sectors=[v['sequence_number'] for v in inventory if v['provenance_name']=='QLP' and v['sequence_number'] not in (6,45,87)],
   assumed_depth=depth,results=out,python=sys.version.split()[0],numpy=np.__version__,
   astropy=astropy.__version__,seconds=time.perf_counter()-start),indent=2))
```

```json
{
  "tic": "307815853",
  "input_receipts": [
    {
      "sequence_number": 6,
      "uri": "mast:HLSP/qlp/s0006/0000/0003/0781/5853/hlsp_qlp_tess_ffi_s0006-0000000307815853_tess_v01_llc.fits",
      "sha256": "981a506772088dfa007ebb00dc39912ce716ecc8193de9273a9f726bc32c6eca",
      "bytes": 66240
    },
    {
      "sequence_number": 45,
      "uri": "mast:HLSP/qlp/s0045/0000/0003/0781/5853/hlsp_qlp_tess_ffi_s0045-0000000307815853_tess_v01_llc.fits",
      "sha256": "e0cfade54f0eae4a2d5155c69d14985a64649ff5ccc5964df3dd952ae5f0442e",
      "bytes": 204480
    },
    {
      "sequence_number": 87,
      "uri": "mast:HLSP/qlp/s0087/0000/0003/0781/5853/hlsp_qlp_tess_ffi_s0087-0000000307815853_tess_v01_llc.fits",
      "sha256": "d0a763d16aa02b427ec6dc973c51925aec07586d86d096b4f458491dd0933a8f",
      "bytes": 691200
    }
  ],
  "inventory_sha256": "81021be609b9ffc565f76aafb40c11210b1ea0baff67d34fe4f1dc8dd5c381ce",
  "source_report_sha256": "83add6622715fb26885bb74680d2254cbf68d640b8cb2c3971a2640d071b5586",
  "unprocessed_t16": [
    {
      "dataRights": "PUBLIC",
      "dataURL": "mast:HLSP/t16/s0006/cam1-ccd3/0033/5614/8674/1498/hlsp_t16_tess_ffi_s0006-cam1-ccd3-03356148674149834752_tess_v01_llc.fits",
      "obs_id": "hlsp_t16_tess_ffi_s0006-cam1-ccd3-03356148674149834752_tess_v01",
      "provenance_name": "T16",
      "s_dec": 14.8888786326,
      "s_ra": 99.2583686621,
      "sequence_number": 6,
      "t_exptime": 1800
    }
  ],
  "unused_qlp_sectors": [
    33,
    43,
    44,
    71,
    72
  ],
  "assumed_depth": 0.0034510450342277173,
  "results": [
    {
      "harmonic": 1,
      "base_trial": 13.926528378848928,
      "cycles_between_anchors": 158,
      "period_days": 13.935714826361837,
      "epoch_btjd": 1478.07283,
      "duration_days": 0.16,
      "events": [
        {
          "sector": 6,
          "cycle": 0,
          "btjd": 1478.07283,
          "anchor": true,
          "state": "insufficient_cadence",
          "baseline": null,
          "trial_injected_scores": [
            null,
            null,
            null
          ]
        },
        {
          "sector": 45,
          "cycle": 76,
          "btjd": 2537.1871568034994,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 2537.1871568034994,
            "duration": 0.16,
            "depth": -0.0011480569285179083,
            "error": 0.0008220845861888249,
            "score": -1.396519224183837,
            "points": 22,
            "cadence": 0.006944895575998089
          },
          "trial_injected_scores": [
            1.2203938674446926,
            2.8014004208331365,
            2.7121077858609937
          ]
        },
        {
          "sector": 87,
          "cycle": 157,
          "btjd": 3665.9800577388082,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 3665.9800577388082,
            "duration": 0.16,
            "depth": 0.00398486239584227,
            "error": 0.0014983232551800865,
            "score": 2.6595478526183065,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          "trial_injected_scores": [
            4.838767668937847,
            4.962819207645716,
            4.590331138961522
          ]
        },
        {
          "sector": 87,
          "cycle": 158,
          "btjd": 3679.9157725651703,
          "anchor": true,
          "state": "supported",
          "baseline": {
            "tc": 3679.9157725651703,
            "duration": 0.16,
            "depth": 0.0069020900684554345,
            "error": 0.0012903401726936813,
            "score": 5.349046875016537,
            "points": 70,
            "cadence": 0.0023148116069933167
          },
          "trial_injected_scores": [
            6.620982928744291,
            8.023570312524805,
            7.064248041032499
          ]
        }
      ]
    },
    {
      "harmonic": 2,
      "base_trial": 6.963264189424464,
      "cycles_between_anchors": 316,
      "period_days": 6.967857413180918,
      "epoch_btjd": 1478.07283,
      "duration_days": 0.16,
      "events": [
        {
          "sector": 6,
          "cycle": -1,
          "btjd": 1471.1049725868193,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 1471.1049725868193,
            "duration": 0.16,
            "depth": 0.0006508334500136414,
            "error": 0.0007944901961086893,
            "score": 0.8191837397130133,
            "points": 8,
            "cadence": 0.02083352332647337
          },
          "trial_injected_scores": [
            5.740359963244502,
            5.162906357223577,
            4.573855238803208
          ]
        },
        {
          "sector": 6,
          "cycle": 0,
          "btjd": 1478.07283,
          "anchor": true,
          "state": "insufficient_cadence",
          "baseline": null,
          "trial_injected_scores": [
            null,
            null,
            null
          ]
        },
        {
          "sector": 6,
          "cycle": 1,
          "btjd": 1485.0406874131809,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 1485.0406874131809,
            "duration": 0.16,
            "depth": 0.00025852405262208973,
            "error": 0.0009260725503686539,
            "score": 0.27916177033772965,
            "points": 7,
            "cadence": 0.02083352332647337
          },
          "trial_injected_scores": [
            3.2584866579742777,
            4.005700293538655,
            3.9335592127309558
          ]
        },
        {
          "sector": 45,
          "cycle": 151,
          "btjd": 2530.219299390319,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 2530.219299390319,
            "duration": 0.16,
            "depth": 0.0007167022711862704,
            "error": 0.0010377819044458752,
            "score": 0.6906097207090486,
            "points": 20,
            "cadence": 0.006944895575998089
          },
          "trial_injected_scores": [
            2.7876102888746557,
            4.016014624613599,
            2.5095210074956658
          ]
        },
        {
          "sector": 45,
          "cycle": 152,
          "btjd": 2537.1871568034994,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 2537.1871568034994,
            "duration": 0.16,
            "depth": -0.0011480569285179083,
            "error": 0.0008220845861888249,
            "score": -1.396519224183837,
            "points": 22,
            "cadence": 0.006944895575998089
          },
          "trial_injected_scores": [
            1.2203938674446926,
            2.8014004208331365,
            2.7121077858609937
          ]
        },
        {
          "sector": 45,
          "cycle": 153,
          "btjd": 2544.1550142166807,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 2544.1550142166807,
            "duration": 0.16,
            "depth": 0.0008703684552591987,
            "error": 0.0012023141591682998,
            "score": 0.7239110083019197,
            "points": 20,
            "cadence": 0.006944895575998089
          },
          "trial_injected_scores": [
            3.8903518895376252,
            3.594246525780127,
            3.9313237262191274
          ]
        },
        {
          "sector": 87,
          "cycle": 314,
          "btjd": 3665.9800577388082,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 3665.9800577388082,
            "duration": 0.16,
            "depth": 0.00398486239584227,
            "error": 0.0014983232551800865,
            "score": 2.6595478526183065,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          "trial_injected_scores": [
            4.838767668937847,
            4.962819207645716,
            4.590331138961522
          ]
        },
        {
          "sector": 87,
          "cycle": 315,
          "btjd": 3672.9479151519895,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 3672.9479151519895,
            "duration": 0.16,
            "depth": -0.0011216441698883568,
            "error": 0.0011705385460499724,
            "score": -0.9582291618446812,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          "trial_injected_scores": [
            3.518856570914584,
            1.9900249096452343,
            0.8739396843821876
          ]
        },
        {
          "sector": 87,
          "cycle": 316,
          "btjd": 3679.9157725651703,
          "anchor": true,
          "state": "supported",
          "baseline": {
            "tc": 3679.9157725651703,
            "duration": 0.16,
            "depth": 0.0069020900684554345,
            "error": 0.0012903401726936813,
            "score": 5.349046875016537,
            "points": 70,
            "cadence": 0.0023148116069933167
          },
          "trial_injected_scores": [
            6.620982928744291,
            8.023570312524805,
            7.064248041032499
          ]
        },
        {
          "sector": 87,
          "cycle": 317,
          "btjd": 3686.883629978351,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 3686.883629978351,
            "duration": 0.16,
            "depth": -0.00004401756427985267,
            "error": 0.001275582470807924,
            "score": -0.034507815282200434,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          "trial_injected_scores": [
            2.126097089852474,
            2.6709582076570344,
            2.881192415954858
          ]
        }
      ]
    },
    {
      "harmonic": 3,
      "base_trial": 4.642176126282976,
      "cycles_between_anchors": 474,
      "period_days": 4.645238275453946,
      "epoch_btjd": 1478.07283,
      "duration_days": 0.16,
      "events": [
        {
          "sector": 6,
          "cycle": -2,
          "btjd": 1468.7823534490922,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 1468.7823534490922,
            "duration": 0.16,
            "depth": 0.0014961793771717867,
            "error": 0.0015295995282505017,
            "score": 0.9781510451190191,
            "points": 8,
            "cadence": 0.02083352332647337
          },
          "trial_injected_scores": [
            2.0267535220129207,
            3.234326580276834,
            5.878489743034389
          ]
        },
        {
          "sector": 6,
          "cycle": -1,
          "btjd": 1473.4275917245461,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 1473.4275917245461,
            "duration": 0.16,
            "depth": -0.0009653682256261897,
            "error": 0.0007980392732043863,
            "score": -1.209675084974106,
            "points": 8,
            "cadence": 0.02083352332647337
          },
          "trial_injected_scores": [
            2.395317758115975,
            3.114729928792514,
            3.7320094749652677
          ]
        },
        {
          "sector": 6,
          "cycle": 0,
          "btjd": 1478.07283,
          "anchor": true,
          "state": "insufficient_cadence",
          "baseline": null,
          "trial_injected_scores": [
            null,
            null,
            null
          ]
        },
        {
          "sector": 6,
          "cycle": 1,
          "btjd": 1482.718068275454,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 1482.718068275454,
            "duration": 0.16,
            "depth": 0.0011290300086629799,
            "error": 0.0008998313890921986,
            "score": 1.2547128521511235,
            "points": 8,
            "cadence": 0.02083352332647337
          },
          "trial_injected_scores": [
            3.9820634881615904,
            5.08992584434217,
            4.74409148869063
          ]
        },
        {
          "sector": 6,
          "cycle": 2,
          "btjd": 1487.363306550908,
          "anchor": false,
          "state": "insufficient_cadence",
          "baseline": null,
          "trial_injected_scores": [
            null,
            null,
            1.8532916722357948
          ]
        },
        {
          "sector": 45,
          "cycle": 226,
          "btjd": 2527.896680252592,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 2527.896680252592,
            "duration": 0.16,
            "depth": -0.0003153775918407969,
            "error": 0.0012257707788388625,
            "score": -0.25728920715465664,
            "points": 21,
            "cadence": 0.006944895575998089
          },
          "trial_injected_scores": [
            4.249465084088247,
            2.558118937504162,
            2.248861653391307
          ]
        },
        {
          "sector": 45,
          "cycle": 227,
          "btjd": 2532.5419185280457,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 2532.5419185280457,
            "duration": 0.16,
            "depth": -0.0002513260160280843,
            "error": 0.0012228469826057996,
            "score": -0.20552531886902686,
            "points": 21,
            "cadence": 0.006944895575998089
          },
          "trial_injected_scores": [
            1.4587369020297192,
            2.6166143955159953,
            2.224944360313396
          ]
        },
        {
          "sector": 45,
          "cycle": 228,
          "btjd": 2537.1871568034994,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 2537.1871568034994,
            "duration": 0.16,
            "depth": -0.0011480569285179083,
            "error": 0.0008220845861888249,
            "score": -1.396519224183837,
            "points": 22,
            "cadence": 0.006944895575998089
          },
          "trial_injected_scores": [
            1.2203938674446926,
            2.8014004208331365,
            2.7121077858609937
          ]
        },
        {
          "sector": 45,
          "cycle": 229,
          "btjd": 2541.8323950789536,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 2541.8323950789536,
            "duration": 0.16,
            "depth": -0.0002731632736049646,
            "error": 0.0010470392036445438,
            "score": -0.26089116114672245,
            "points": 23,
            "cadence": 0.006944895575998089
          },
          "trial_injected_scores": [
            4.0548276820695754,
            3.0351124862957874,
            2.8805052862469487
          ]
        },
        {
          "sector": 45,
          "cycle": 230,
          "btjd": 2546.4776333544078,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 2546.4776333544078,
            "duration": 0.16,
            "depth": 0.0012591347156271245,
            "error": 0.0010357467975920932,
            "score": 1.2156781160746661,
            "points": 22,
            "cadence": 0.006944895575998089
          },
          "trial_injected_scores": [
            4.251875019005028,
            4.547617005242091,
            4.523037375399982
          ]
        },
        {
          "sector": 87,
          "cycle": 471,
          "btjd": 3665.9800577388087,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 3665.9800577388087,
            "duration": 0.16,
            "depth": 0.00398486239584227,
            "error": 0.0014983232551800865,
            "score": 2.6595478526183065,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          "trial_injected_scores": [
            4.838767668937847,
            4.962819207645716,
            4.590331138961522
          ]
        },
        {
          "sector": 87,
          "cycle": 472,
          "btjd": 3670.6252960142624,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 3670.6252960142624,
            "duration": 0.16,
            "depth": 0.0006395658626079802,
            "error": 0.0010457228440124402,
            "score": 0.6116016937661656,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          "trial_injected_scores": [
            4.195138610472769,
            3.9117543623126916,
            3.5546508999599724
          ]
        },
        {
          "sector": 87,
          "cycle": 473,
          "btjd": 3675.2705342897166,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 3675.2705342897166,
            "duration": 0.16,
            "depth": -0.0011413879784090081,
            "error": 0.001235320401432885,
            "score": -0.9239610849825504,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          "trial_injected_scores": [
            1.7833506171831608,
            1.8696825966281048,
            1.580910277987094
          ]
        },
        {
          "sector": 87,
          "cycle": 474,
          "btjd": 3679.9157725651703,
          "anchor": true,
          "state": "supported",
          "baseline": {
            "tc": 3679.9157725651703,
            "duration": 0.16,
            "depth": 0.0069020900684554345,
            "error": 0.0012903401726936813,
            "score": 5.349046875016537,
            "points": 70,
            "cadence": 0.0023148116069933167
          },
          "trial_injected_scores": [
            6.620982928744291,
            8.023570312524805,
            7.064248041032499
          ]
        },
        {
          "sector": 87,
          "cycle": 475,
          "btjd": 3684.5610108406245,
          "anchor": false,
          "state": "insufficient_sensitivity",
          "baseline": {
            "tc": 3684.5610108406245,
            "duration": 0.16,
            "depth": -0.000580956577412698,
            "error": 0.0011106993641915998,
            "score": -0.5230547492349882,
            "points": 69,
            "cadence": 0.0023148116069933167
          },
          "trial_injected_scores": [
            2.2249917712431055,
            2.5840371835487232,
            2.556979836738705
          ]
        }
      ]
    }
  ],
  "python": "3.13.15",
  "numpy": "2.5.4",
  "astropy": "8.0.2",
  "seconds": 0.07213420001789927
}
```

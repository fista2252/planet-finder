# TIC 307815853: fixed-window flux and background check

**Disposition: unresolved.** The original T16 single event is a known catalog entry and cannot be checked with the quality-zero QLP cadences at its stated time. Two later QLP dimmings near the 13.92653 d trial have aperture-dependent depths and simultaneous background changes. This warrants an alternate reduction or pixel check before any periodic interpretation.

## Objective and input selection

The Sector 62 pixel test for TIC 140920261 remains queued behind the [archive execution blocker](../tic140920261-unused-sectors/2026-11-21-checkpoint.md). This session used the other non-pilot T16 target in the [frozen November 13 selection](pilot-extension-selection.json), TIC 307815853, Gaia DR3 3356148674149834752, ICRS 99.2582080078°, +14.8888976002°. The [Roth et al. T16 catalog](https://www.ipac.caltech.edu/publication/2026ApJS..284...19R) reports an approximate S6 event time of BTJD 1478.07283 and duration 0.1425 d in the preserved selection. The 13.926528378848928 d blind-search trial instead used epoch 1479.5242887070394, 1.45146 days later. Do not count the catalog event as support for that ephemeris. The earlier run found one nominally supported cycle, no accepted period, and a truncated period grid.

Reused exactly the three cached QLP v01 FITS from the original worker report, which have 824, 2,949, and 8,437 common finite quality-zero cadences. This adds aperture, local background and reported position comparisons to the existing time-window analysis, with no period or timing optimization. FITS time is TDB with BJD reference 2457000. The box duration is 0.16 d for the frozen blind-trial windows and 0.1425 d for the published T16 time. Normalize each flux column by its sector median and apply the existing sideband scorer. For each duration, controls are spaced by four durations and exclude centers within four durations of any listed center in that sector. These controls are descriptive: the windows were chosen from an earlier search of the same data, and controls are neither independent observations nor a calibrated false-alarm estimate. [MAST QLP product documentation](https://archive.stsci.edu/hlsp/qlp)

## Results

| Window | Retained in-box cadences | Raw SAP depth | Small / default / large detrended-aperture depth |
| --- | ---: | ---: | ---: |
| S6 catalog, 1478.07283 | 0 | Unmeasurable | Unmeasurable |
| S6 trial, 1479.52429 | 8 | 0.150% | 0.091% / 0.152% / 0.182% |
| S45 trial, 2537.94045 | 8 | Unmeasurable with sidebands | Unmeasurable with sidebands |
| S87 cycle 157, 3665.98924 | 69 | 0.460% | 0.204% / 0.445% / 0.649% |
| S87 cycle 158, 3679.91577 | 70 | 0.739% | 0.326% / 0.690% / 0.813% |

The cached S6 QLP table has two rows in the T16 event box, both flagged QUALITY=4096. S45 has 18 rows in-box, eight pass the mask, and the existing coverage rule rejects the window. These are coverage limits. They do not count as missing transits.

Both S87 windows have a decreased QLP background estimate, by 286.79 and 230.96 in the native column units. For each event only one of the 25 usable controls has a background change of at least that absolute size. No control has raw or default-aperture dimming as deep as either selected window. Position-change magnitudes are 0.00171 and 0.00271 px, with 10 and six control windows at least as large. Column changes may reflect extraction residuals or a blend and do not locate the dimming source. Cycle 157 remains below the earlier nominal support threshold and does not become a new supported cycle. A true eclipse could coexist with a background feature, so the measured correlation alone does not establish cause.

November 13 exact-TIC TOI and confirmed-planet queries returned zero rows, and the 740-row pinned Eumemic table at `e2eee928fd621f75c1621f6b4206d34bb222ae59` had no exact match. Exact TIC and Gaia web searches today surfaced no usable additional match. This does not cover all CTOIs, eclipsing binaries or sky-coordinate duplicates, and the T16 entry has published priority for its original event. No novelty or new candidate is claimed.

## Reproduction and next discriminating check

Run the script below from isolated checkout `25512680d1e4a9ac2fda85eb6836367f88bf8e61`. Prior report `runs/worker-order-pilot-extension-3b8bef46/report.json` SHA-256 `83add6622715fb26885bb74680d2254cbf68d640b8cb2c3971a2640d071b5586`; frozen selection SHA-256 `c00896cfd40f9ac9416072c6feba9c94a8d3292abeb44075a3708a6f99999581`. For each FITS the archive path is `mast:HLSP/qlp/s{sector:04d}/0000/0003/0781/5853/{file}`, and hashes are in the output below. New archive bytes: zero. Computation: 0.065 s in-process, no code or threshold changes. The local checkout was left untouched.

First check public archive availability for an S87 alternate light curve or 7-by-7-pixel cutout, then compare the two frozen S87 windows against nearby background controls and examine whether the deficit is spatially concentrated. If S87 pixel access is unavailable, check the original T16 vetting artifact or a reduction covering the flagged S6 interval before using that epoch in another search.

```python
import hashlib, json, sys, time
from pathlib import Path
import astropy, numpy as np
from astropy.io import fits
sys.path.insert(0,'legacy')
from screen import scores

start = time.perf_counter()
frozen = {6:[('catalog',1478.07283,.1425),('trial',1479.5242887070394,.16)],
          45:[('trial',2537.940445499558,.16)],
          87:[('trial_157',3665.9892441863212,.16),('trial_158',3679.9157725651703,.16)]}
results = []
for sector, windows in frozen.items():
    path, = Path('.audit-cache/shared/fits').glob(f'*/hlsp_qlp_tess_ffi_s{sector:04d}-0000000307815853_tess_v01_llc.fits')
    with fits.open(path, memmap=False) as f:
        hp, h = f[0].header, f[1].header
        assert h['TIMESYS']=='TDB' and h['BJDREFI']==2457000 and h.get('BJDREFF',0)==0
        a = f[1].data
        values = {'raw':'SAP_FLUX','x':'SAP_X','y':'SAP_Y','background':'SAP_BKG'}
        for label, suffix in [('default',''),('small','_SML'),('large','_LAG')]:
            values[label] = next(k+suffix for k in ('DET_FLUX','KSPSAP_FLUX') if k+suffix in a.names)
        mask = a['QUALITY']==0
        for col in ['TIME',*values.values()]:
            mask &= np.isfinite(a[col])
        b = a[mask]
        t, ix = np.unique(b['TIME'],return_index=True)
        b = b[ix]
        out = dict(sector=sector,camera=hp.get('CAMERA'),ccd=hp.get('CCD'),n=len(t),
                   file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),columns=values,windows=[])
        for label, center, duration in windows:
            grid = np.arange(t[0]+2*duration,t[-1]-2*duration,4*duration)
            for _, event_center, _ in windows:
                grid = grid[abs(grid-event_center)>4*duration]
            event, controls = {}, {}
            for name, key in values.items():
                v = b[key].astype(float)
                if name in ('raw','default','small','large'):
                    v = v/np.median(v)-1
                e = scores(t,v,np.array([center]),duration)
                event[name] = e[0] if e else None
                controls[name] = {r['tc']:r['depth'] for r in scores(t,v,grid,duration)}
            base = dict(label=label,center=center,duration=duration,
                        in_total=int(np.count_nonzero(abs(a['TIME']-center)<duration/2)),
                        in_retained=int(np.count_nonzero(abs(t-center)<duration/2)),
                        in_quality={str(int(q)):int(k) for q,k in zip(*np.unique(a['QUALITY'][abs(a['TIME']-center)<duration/2],return_counts=True))},
                        measurements=event)
            if event['raw'] is not None:
                vec = float(np.hypot(event['x']['depth'],event['y']['depth']))
                xy = np.array([np.hypot(x,controls['y'][c]) for c,x in controls['x'].items()])
                base.update(xy_norm=vec,n_controls=len(xy),
                            controls_with_xy_ge=int(sum(xy>=vec)),
                            controls_with_raw_depth_ge=int(sum(d>=event['raw']['depth'] for d in controls['raw'].values())),
                            controls_with_default_depth_ge=int(sum(d>=event['default']['depth'] for d in controls['default'].values())),
                            controls_with_bkg_change_abs_ge=int(sum(abs(d)>=abs(event['background']['depth']) for d in controls['background'].values())))
            out['windows'].append(base)
        results.append(out)
print(json.dumps(dict(results=results,seconds=time.perf_counter()-start,
    versions=dict(python=sys.version.split()[0],numpy=np.__version__,astropy=astropy.__version__)),indent=2))
```

```json
{
  "results": [
    {
      "sector": 6,
      "camera": 1,
      "ccd": 3,
      "n": 824,
      "file": "hlsp_qlp_tess_ffi_s0006-0000000307815853_tess_v01_llc.fits",
      "sha256": "981a506772088dfa007ebb00dc39912ce716ecc8193de9273a9f726bc32c6eca",
      "columns": {
        "raw": "SAP_FLUX",
        "x": "SAP_X",
        "y": "SAP_Y",
        "background": "SAP_BKG",
        "default": "KSPSAP_FLUX",
        "small": "KSPSAP_FLUX_SML",
        "large": "KSPSAP_FLUX_LAG"
      },
      "windows": [
        {
          "label": "catalog",
          "center": 1478.07283,
          "duration": 0.1425,
          "in_total": 2,
          "in_retained": 0,
          "in_quality": {
            "4096": 2
          },
          "measurements": {
            "raw": null,
            "x": null,
            "y": null,
            "background": null,
            "default": null,
            "small": null,
            "large": null
          }
        },
        {
          "label": "trial",
          "center": 1479.5242887070394,
          "duration": 0.16,
          "in_total": 8,
          "in_retained": 8,
          "in_quality": {
            "0": 8
          },
          "measurements": {
            "raw": {
              "tc": 1479.5242887070394,
              "duration": 0.16,
              "depth": 0.001504088873518705,
              "error": 0.0016532198704428595,
              "score": 0.9097936096762461,
              "points": 8,
              "cadence": 0.02083352332647337
            },
            "x": {
              "tc": 1479.5242887070394,
              "duration": 0.16,
              "depth": -0.0001284208524339192,
              "error": 0.0007226020417068718,
              "score": -0.17772002433119882,
              "points": 8,
              "cadence": 0.02083352332647337
            },
            "y": {
              "tc": 1479.5242887070394,
              "duration": 0.16,
              "depth": 0.0008609395171674805,
              "error": 0.0012527779298852884,
              "score": 0.6872243648531653,
              "points": 8,
              "cadence": 0.02083352332647337
            },
            "background": {
              "tc": 1479.5242887070394,
              "duration": 0.16,
              "depth": -10.671516640543473,
              "error": 93.22028202021403,
              "score": -0.11447633936818004,
              "points": 8,
              "cadence": 0.02083352332647337
            },
            "default": {
              "tc": 1479.5242887070394,
              "duration": 0.16,
              "depth": 0.001524813737158113,
              "error": 0.00169539026406479,
              "score": 0.8993880462084816,
              "points": 8,
              "cadence": 0.02083352332647337
            },
            "small": {
              "tc": 1479.5242887070394,
              "duration": 0.16,
              "depth": 0.0009108614835042956,
              "error": 0.001035491347356125,
              "score": 0.8796418104602792,
              "points": 8,
              "cadence": 0.02083352332647337
            },
            "large": {
              "tc": 1479.5242887070394,
              "duration": 0.16,
              "depth": 0.0018241297112665841,
              "error": 0.0022043705747726416,
              "score": 0.8275059248850319,
              "points": 8,
              "cadence": 0.02083352332647337
            }
          },
          "xy_norm": 0.0008704646848439223,
          "n_controls": 24,
          "controls_with_xy_ge": 10,
          "controls_with_raw_depth_ge": 2,
          "controls_with_default_depth_ge": 2,
          "controls_with_bkg_change_abs_ge": 21
        }
      ]
    },
    {
      "sector": 45,
      "camera": 1,
      "ccd": 1,
      "n": 2949,
      "file": "hlsp_qlp_tess_ffi_s0045-0000000307815853_tess_v01_llc.fits",
      "sha256": "e0cfade54f0eae4a2d5155c69d14985a64649ff5ccc5964df3dd952ae5f0442e",
      "columns": {
        "raw": "SAP_FLUX",
        "x": "SAP_X",
        "y": "SAP_Y",
        "background": "SAP_BKG",
        "default": "KSPSAP_FLUX",
        "small": "KSPSAP_FLUX_SML",
        "large": "KSPSAP_FLUX_LAG"
      },
      "windows": [
        {
          "label": "trial",
          "center": 2537.940445499558,
          "duration": 0.16,
          "in_total": 18,
          "in_retained": 8,
          "in_quality": {
            "0": 8,
            "2048": 1,
            "4096": 8,
            "6144": 1
          },
          "measurements": {
            "raw": null,
            "x": null,
            "y": null,
            "background": null,
            "default": null,
            "small": null,
            "large": null
          }
        }
      ]
    },
    {
      "sector": 87,
      "camera": 1,
      "ccd": 4,
      "n": 8437,
      "file": "hlsp_qlp_tess_ffi_s0087-0000000307815853_tess_v01_llc.fits",
      "sha256": "d0a763d16aa02b427ec6dc973c51925aec07586d86d096b4f458491dd0933a8f",
      "columns": {
        "raw": "SAP_FLUX",
        "x": "SAP_X",
        "y": "SAP_Y",
        "background": "SAP_BKG",
        "default": "DET_FLUX",
        "small": "DET_FLUX_SML",
        "large": "DET_FLUX_LAG"
      },
      "windows": [
        {
          "label": "trial_157",
          "center": 3665.9892441863212,
          "duration": 0.16,
          "in_total": 69,
          "in_retained": 69,
          "in_quality": {
            "0": 69
          },
          "measurements": {
            "raw": {
              "tc": 3665.9892441863212,
              "duration": 0.16,
              "depth": 0.004596368935568644,
              "error": 0.0015012291373165774,
              "score": 3.0617370935023143,
              "points": 69,
              "cadence": 0.0023148116069933167
            },
            "x": {
              "tc": 3665.9892441863212,
              "duration": 0.16,
              "depth": -0.0016270452825324355,
              "error": 0.00101127450404479,
              "score": -1.608905669058945,
              "points": 69,
              "cadence": 0.0023148116069933167
            },
            "y": {
              "tc": 3665.9892441863212,
              "duration": 0.16,
              "depth": -0.0005265323551100209,
              "error": 0.0007456962176561366,
              "score": -0.7060949789513632,
              "points": 69,
              "cadence": 0.0023148116069933167
            },
            "background": {
              "tc": 3665.9892441863212,
              "duration": 0.16,
              "depth": 286.78723627102136,
              "error": 93.53169955059496,
              "score": 3.066203625604888,
              "points": 69,
              "cadence": 0.0023148116069933167
            },
            "default": {
              "tc": 3665.9892441863212,
              "duration": 0.16,
              "depth": 0.004453110949046289,
              "error": 0.0015152714403616862,
              "score": 2.938820616841665,
              "points": 69,
              "cadence": 0.0023148116069933167
            },
            "small": {
              "tc": 3665.9892441863212,
              "duration": 0.16,
              "depth": 0.002035885420053831,
              "error": 0.0008790635523604676,
              "score": 2.3159706878837794,
              "points": 69,
              "cadence": 0.0023148116069933167
            },
            "large": {
              "tc": 3665.9892441863212,
              "duration": 0.16,
              "depth": 0.006494137670040106,
              "error": 0.002022892802298376,
              "score": 3.2103221993086235,
              "points": 69,
              "cadence": 0.0023148116069933167
            }
          },
          "xy_norm": 0.0017101206601841748,
          "n_controls": 25,
          "controls_with_xy_ge": 10,
          "controls_with_raw_depth_ge": 0,
          "controls_with_default_depth_ge": 0,
          "controls_with_bkg_change_abs_ge": 1
        },
        {
          "label": "trial_158",
          "center": 3679.9157725651703,
          "duration": 0.16,
          "in_total": 70,
          "in_retained": 70,
          "in_quality": {
            "0": 70
          },
          "measurements": {
            "raw": {
              "tc": 3679.9157725651703,
              "duration": 0.16,
              "depth": 0.007385813889099514,
              "error": 0.0012778018646507757,
              "score": 5.780093215874327,
              "points": 70,
              "cadence": 0.0023148116069933167
            },
            "x": {
              "tc": 3679.9157725651703,
              "duration": 0.16,
              "depth": -0.002692049777692773,
              "error": 0.0009203568960047006,
              "score": -2.9250063637041777,
              "points": 70,
              "cadence": 0.0023148116069933167
            },
            "y": {
              "tc": 3679.9157725651703,
              "duration": 0.16,
              "depth": -0.0003440418598188444,
              "error": 0.0006762757366257328,
              "score": -0.5087301542643476,
              "points": 70,
              "cadence": 0.0023148116069933167
            },
            "background": {
              "tc": 3679.9157725651703,
              "duration": 0.16,
              "depth": 230.95847445772188,
              "error": 71.78893157977491,
              "score": 3.2171877944871072,
              "points": 70,
              "cadence": 0.0023148116069933167
            },
            "default": {
              "tc": 3679.9157725651703,
              "duration": 0.16,
              "depth": 0.0069020900684554345,
              "error": 0.0012903401726936813,
              "score": 5.349046875016537,
              "points": 70,
              "cadence": 0.0023148116069933167
            },
            "small": {
              "tc": 3679.9157725651703,
              "duration": 0.16,
              "depth": 0.0032577812797690113,
              "error": 0.0007718525787038695,
              "score": 4.220729929074834,
              "points": 70,
              "cadence": 0.0023148116069933167
            },
            "large": {
              "tc": 3679.9157725651703,
              "duration": 0.16,
              "depth": 0.00812809934811831,
              "error": 0.001743574156029224,
              "score": 4.661745713545708,
              "points": 70,
              "cadence": 0.0023148116069933167
            }
          },
          "xy_norm": 0.0027139448791166187,
          "n_controls": 25,
          "controls_with_xy_ge": 6,
          "controls_with_raw_depth_ge": 0,
          "controls_with_default_depth_ge": 0,
          "controls_with_bkg_change_abs_ge": 1
        }
      ]
    }
  ],
  "seconds": 0.06469819997437298,
  "versions": {
    "python": "3.13.15",
    "numpy": "2.5.4",
    "astropy": "8.0.2"
  }
}
```

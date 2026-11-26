# TIC 307815853: Sector 87 local background control

**Disposition: unresolved.** A local linear model of residual flux versus the estimated background explains only 0.5% and 4.0% of the two dips' raw depths. It also worsens the depth scatter in comparison windows. The simple model supplies no basis to subtract these events or to adopt it as a search correction. Spatially varying background, a blend, and other reduction effects remain plausible.

The [November 25 fixed-window check](tic307815853-fixed-windows.md) found background changes at both Sector 87 dips. With the [archive access limitation](../tic140920261-unused-sectors/2026-11-21-checkpoint.md) still in place, this session tested whether the association predicts comparable flux changes in nearby data outside the events. This is an experimental diagnostic using the same raw QLP SAP flux as before. QLP reports a background estimate that has already been used in its extraction, so this checks for a residual association and does not treat SAP_BKG as unapplied sky flux. MAST confirms that this sector predates the Sector 94 extraction change. [QLP archive and product history](https://archive.stsci.edu/hlsp/qlp) · [QLP data product description](https://archive.stsci.edu/hlsps/qlp/hlsp_qlp_tess_ffi_all_tess_v1_data-prod-desc.pdf)

## Fixed comparison

Freeze both previous centers, BTJD 3665.9892441863212 and 3679.9157725651703, and the 0.16-day box. For each event or control, fit ordinary least squares on quality-zero, finite observations 0.32–0.96 d away on either side, in the same spacecraft orbit, omitting the central neighborhoods of both events and requiring at least 20 training points per side. Fit an intercept, local time trend, and standardized background. Subtract just the fitted background term, then remeasure with the existing local sideband scorer. The sideband scoring removes the local linear time trend itself. Use the earlier 0.64-day control spacing with centers more than 0.64 d from either event. This reproduces held-out prediction for control windows on the same light curve, not independent validation.

| Frozen cycle | Training cadences | Raw depth | Background-adjusted depth | Estimated background contribution |
| --- | ---: | ---: | ---: | ---: |
| 157 | 329 | 0.4596% | 0.4574% | 0.0022 percentage points |
| 158 | 439 | 0.7386% | 0.7088% | 0.0298 percentage points |

Neither event's in-box background values exceeded the local training range. The control box-depth RMS rose from 0.0989% to 0.1359% across 24 usable control windows. Eight proposed controls failed baseline cadence or local training requirements. This RMS is descriptive and does not estimate a detection probability. The scalar archive background estimate does not measure a spatial residual in the image, and the fit explores one functional form only. The trial retains one nominally supported event under the original pipeline and no accepted period. No alternative model was tuned or applied to production.

## Archive and field definitions

MAST documents the [TGLC archive](https://stdatu.stsci.edu/hlsp/tglc) and [TESS-SPOC archive](https://archive.stsci.edu/hlsp/tess-spoc). This session did not verify a target-specific alternate Sector 87 product or obtain pixels. TGLC aperture and PSF light curves would be distinct reductions of the same TESS exposures if a matching product is available. The next discriminating step is a sector-bounded pixel cutout, comparing each frozen box with local sidebands and checking background and residual structure before any location claim. For TESScut timestamp alignment, MAST states that FFI barycentric corrections use the CCD center. Verify the target timing convention rather than assuming it matches QLP. [MAST TESS product overview](https://outerspace.stsci.edu/spaces/TESS/pages/14563420/2.0%2B-%2BData%2BProduct%2BOverview)

The QLP linked format guide explicitly documents SAP_X/Y as aperture flux-weighted centroid coordinates for primary-mission products. This clarifies the [earlier Sector 9 note](../tic140920261-unused-sectors/sector9-position-check.md), which could not determine that from the FITS headers. It still does not make QLP and SPOC centroids interchangeable across different apertures or establish that every extended-mission release uses the identical implementation.

Cache input `mast:HLSP/qlp/s0087/0000/0003/0781/5853/hlsp_qlp_tess_ffi_s0087-0000000307815853_tess_v01_llc.fits`. SHA-256 and dependency versions below. Source checkout `25512680d1e4a9ac2fda85eb6836367f88bf8e61`; clean at start. New FITS bytes: zero. In-process runtime 0.149 s. The existing source files and thresholds were preserved.

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
orb=a['ORBITID']
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
    if len(eval_orbits)!=1:
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
      "train": 461,
      "box": 69,
      "raw_depth": -0.0007980532147605773,
      "adjusted_depth": -0.0012435337394671236,
      "predicted_background_depth": 0.00044548052470623793,
      "normalized_flux_per_native_bkg": -0.00005209837220699644,
      "bkg_z_max_abs_in_box": 2.40239507011276,
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
      "train": 461,
      "box": 69,
      "raw_depth": -0.0010933084302461576,
      "adjusted_depth": -0.0007645586398622346,
      "predicted_background_depth": -0.0003287497903842551,
      "normalized_flux_per_native_bkg": -0.000043653774141654294,
      "bkg_z_max_abs_in_box": 2.0225630700874637,
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
      "status": "measured",
      "train": 374,
      "box": 50,
      "raw_depth": 0.002049324108650611,
      "adjusted_depth": 0.004566628292312071,
      "predicted_background_depth": -0.0025173041836613168,
      "normalized_flux_per_native_bkg": -0.000002030121143582957,
      "bkg_z_max_abs_in_box": 0.6588205966890972,
      "in_box_bkg_outside_training_range": 0
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
      "train": 401,
      "box": 55,
      "raw_depth": -0.00069333549876186,
      "adjusted_depth": -0.001979219808723013,
      "predicted_background_depth": 0.0012858843099600333,
      "normalized_flux_per_native_bkg": -0.000037779080745700963,
      "bkg_z_max_abs_in_box": 2.8330531164700217,
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
      "train": 437,
      "box": 69,
      "raw_depth": 0.000009804762052194327,
      "adjusted_depth": -0.0001229869143314329,
      "predicted_background_depth": 0.00013279167638418027,
      "normalized_flux_per_native_bkg": -0.000029269857264824622,
      "bkg_z_max_abs_in_box": 1.7428675429647482,
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
  "controls": 24,
  "control_rms_raw": 0.0009891558512955298,
  "control_rms_adjusted": 0.0013587362141085667,
  "python": "3.13.15",
  "numpy": "2.5.4",
  "astropy": "8.0.2",
  "seconds": 0.1491712999995798
}
```

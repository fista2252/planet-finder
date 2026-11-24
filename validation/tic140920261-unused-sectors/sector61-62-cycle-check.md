# TIC 140920261: fixed cycle comparison in Sectors 61 and 62

**Disposition: unresolved.** The deep QLP dips occupy opposite cycle parities in the two sectors. Within each sector, the usable intervening cycle lacks a comparable flux dip or QLP position excursion. This weakens a stable single-period eclipse at either frozen trial period.

## Why this check

The earlier [Sector 61](2026-11-21-checkpoint.md) and [Sector 62](sector62-aperture-check.md) aperture work examined the strongest events against generic control windows. This check compares **every** in-span cycle at the frozen 1.9256784171030665 d ephemeris, anchored at BTJD 2985.631834841218. It tests whether the pair spacing reflects a common 3.851356834206133 d phase, and whether a comparable signal occurs on intervening 1.9256784171030665 d cycles. A fixed phase at twice the trial period would retain the same parity across sectors. Phase, duration and search thresholds were not refit. Depth and position comparisons here are descriptive. The former search treated many individual windows as insufficiently sensitive, and this check does not change those formal dispositions.

| Sector | Deep cycles | Intervening cycle | Raw SAP depth at intervening cycle | QLP position change at intervening cycle |
| --- | --- | ---: | ---: | ---: |
| 61 | −2, 0 (even) | −1 | +0.050% | 0.00048 px |
| 62 | 11, 13 (odd) | 12 | −0.087% | 0.00055 px |

The four previously identified windows give raw flux depths of 0.903–1.062% and position-vector magnitudes of 0.00853–0.00969 px. The intervening windows contain 70 and 69 quality-zero samples, respectively. Across both sectors, 23 of 26 in-span windows satisfied the existing cadence rules. The 19 usable windows outside the four selected events reach at most 0.297% raw dimming and 0.00381 px position change. Positive flux depths indicate dimming. The exact values and missing windows are below.

The two sector pairs alone can each resemble a 3.85136 d sequence, but they cannot both lie on one parity of the frozen clock. A double-eclipse interpretation with variable depths, timing errors, a blend, or an instrumental or extraction effect remain possibilities. The four strong events were previously selected using these data. These comparisons therefore supply neither a survey false-alarm estimate nor a blind validation. The earlier [Sector 9 comparison](sector9-position-check.md) also leaves the source unresolved. QLP position columns are retained as diagnostic measurements and are not assumed to be target flux-weighted centroids. The [QLP extraction paper](https://research.usq.edu.au/download/7f686e690658d98f289ed50a75204b6365664979c03301e05ec31c705897de7c/919161/Huang_2020_Res._Notes_AAS_4_204.pdf) describes how QLP estimates per-frame positions for aperture extraction.

## Provenance and next test

Used read-only cached QLP FITS, no new bytes or paid services. See [original baseline](protocol.json) and [unused-sector receipts](result.json) for archive URIs. Both hashes match the receipts. The clean isolated source checkout was `25512680d1e4a9ac2fda85eb6836367f88bf8e61`. Measurement took 0.098 s in-process. No code change or acceptance threshold change. The cutoff at the prior archive execution blocker was preserved. Next compare separately the two Sector 62 pixel deficit maps and their local and pointing-matched controls, with the intervening cycle 12 as an additional fixed comparison. Check the target-versus-neighbor location with the cutout WCS.

Run from the repository root in the shared cache checkout:

```python
import hashlib, json, sys, time
from pathlib import Path
import astropy, numpy as np
from astropy.io import fits
sys.path.insert(0,'legacy')
from screen import scores

start = time.perf_counter()
epoch = 2985.631834841218
half = 1.9256784171030665
duration = .16
fields = ['SAP_FLUX','DET_FLUX','DET_FLUX_SML','DET_FLUX_LAG','SAP_X','SAP_Y','SAP_BKG']
results = {}
for sector, selected in [(61,[2981.780478007012,2985.631834841218]), (62,[3006.81430,3010.66565])]:
    p, = Path('.audit-cache/shared/fits').glob(f'*/hlsp_qlp_tess_ffi_s{sector:04d}-0000000140920261_tess_v01_llc.fits')
    with fits.open(p,memmap=False) as h:
        assert h[1].header['TIMESYS']=='TDB' and h[1].header['BJDREFI']==2457000
        a = h[1].data
        mask = a['QUALITY']==0
        for key in ['TIME',*fields]:
            mask &= np.isfinite(a[key])
        a = a[mask]
    t, ix = np.unique(a['TIME'],return_index=True)
    a = a[ix]
    n = np.arange(int(np.ceil((t[0]-epoch)/half)),int(np.floor((t[-1]-epoch)/half))+1)
    cs = epoch + n*half
    # Use the same exact ephemeris for every row, including the selected events.
    selected_n = np.rint((np.array(selected)-epoch)/half).astype(int)
    cols = {}
    for key in fields:
        y = a[key].astype(float)
        if 'FLUX' in key:
            y = y/np.median(y)-1
        cols[key] = {v['tc']:v for v in scores(t,y,cs,duration)}
    rows = []
    for k,c in zip(n,cs):
        row = dict(cycle=int(k),btjd=float(c),selected=bool(k in selected_n),valid=bool(c in cols['SAP_FLUX']))
        if row['valid']:
            row['depths'] = {key:cols[key][c]['depth'] for key in fields}
            row['sap_nominal_score'] = cols['SAP_FLUX'][c]['score']
            row['points'] = cols['SAP_FLUX'][c]['points']
            row['xy_norm_px'] = float(np.hypot(row['depths']['SAP_X'],row['depths']['SAP_Y']))
        rows.append(row)
    results[str(sector)] = dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),cadences=len(t),windows=rows)
print(json.dumps(dict(epoch=epoch,half_period=half,duration=duration,python=sys.version.split()[0],
    numpy=np.__version__,astropy=astropy.__version__,sectors=results,seconds=time.perf_counter()-start),indent=2))
```

```json
{
  "epoch": 2985.631834841218,
  "half_period": 1.9256784171030665,
  "duration": 0.16,
  "python": "3.13.15",
  "numpy": "2.5.4",
  "astropy": "8.0.2",
  "sectors": {
    "61": {
      "file": "hlsp_qlp_tess_ffi_s0061-0000000140920261_tess_v01_llc.fits",
      "sha256": "4b2462b61fd3eb7212339eb54c3c0d1b9da9f718ecf48b223e79bdb61dcf0afd",
      "cadences": 9326,
      "windows": [
        {
          "cycle": -11,
          "btjd": 2964.449372253084,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": 0.0004435278033714386,
            "DET_FLUX": 0.0006741911680522903,
            "DET_FLUX_SML": 0.00040366167316617757,
            "DET_FLUX_LAG": 0.000826417538241101,
            "SAP_X": 0.00006763685655641893,
            "SAP_Y": 0.0001788690879038768,
            "SAP_BKG": -6.411733779938221
          },
          "sap_nominal_score": 0.5816867241804344,
          "points": 69,
          "xy_norm_px": 0.00019122995312554565
        },
        {
          "cycle": -10,
          "btjd": 2966.375050670187,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.00009751515602433565,
            "DET_FLUX": -0.00022400244221101537,
            "DET_FLUX_SML": -0.000055092112635534695,
            "DET_FLUX_LAG": -0.00033975703703457137,
            "SAP_X": -0.00048808978829129046,
            "SAP_Y": -0.0003574244633455237,
            "SAP_BKG": -7.0695177633163695
          },
          "sap_nominal_score": -0.14613373511470237,
          "points": 69,
          "xy_norm_px": 0.000604966022543475
        },
        {
          "cycle": -9,
          "btjd": 2968.30072908729,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.0009970484764518267,
            "DET_FLUX": -0.0011196036726206307,
            "DET_FLUX_SML": -0.0007537589742059787,
            "DET_FLUX_LAG": -0.0016738562191847012,
            "SAP_X": 0.0005223624778523117,
            "SAP_Y": 0.0005092481414976646,
            "SAP_BKG": 17.228600500833693
          },
          "sap_nominal_score": -1.5287243248465123,
          "points": 70,
          "xy_norm_px": 0.0007295178050512765
        },
        {
          "cycle": -8,
          "btjd": 2970.226407504393,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": 0.00013140893327884392,
            "DET_FLUX": 0.00032665698555708964,
            "DET_FLUX_SML": 0.00008496637797929298,
            "DET_FLUX_LAG": 0.0003946830801366451,
            "SAP_X": -0.00008870307914321529,
            "SAP_Y": 0.00032233347792498535,
            "SAP_BKG": -5.461552879481744
          },
          "sap_nominal_score": 0.1865985546860531,
          "points": 69,
          "xy_norm_px": 0.00033431587943246805
        },
        {
          "cycle": -7,
          "btjd": 2972.1520859214966,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": 0.002967137847286308,
            "DET_FLUX": 0.0021197187210549684,
            "DET_FLUX_SML": 0.0013185714188790997,
            "DET_FLUX_LAG": 0.002633888124440019,
            "SAP_X": -0.00047960242272893083,
            "SAP_Y": 0.00013320683888196072,
            "SAP_BKG": -24.96437815209702
          },
          "sap_nominal_score": 4.2703322141991675,
          "points": 69,
          "xy_norm_px": 0.0004977575170827506
        },
        {
          "cycle": -6,
          "btjd": 2974.0777643385995,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.001461622105528432,
            "DET_FLUX": -0.0015545712862929583,
            "DET_FLUX_SML": -0.0008905360935235391,
            "DET_FLUX_LAG": -0.0021670041324487266,
            "SAP_X": 0.0002887501983650509,
            "SAP_Y": -0.00018185452381658251,
            "SAP_BKG": 12.139933339870296
          },
          "sap_nominal_score": -1.989392195713762,
          "points": 69,
          "xy_norm_px": 0.000341244406384064
        },
        {
          "cycle": -5,
          "btjd": 2976.0034427557025,
          "selected": false,
          "valid": false
        },
        {
          "cycle": -4,
          "btjd": 2977.9291211728055,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": 0.0009748813744972309,
            "DET_FLUX": 0.0008725814804980918,
            "DET_FLUX_SML": 0.0005189265804907909,
            "DET_FLUX_LAG": 0.0012214990342375286,
            "SAP_X": 0.0008459133716769429,
            "SAP_Y": 0.0006839536047209549,
            "SAP_BKG": -11.562242207692137
          },
          "sap_nominal_score": 1.3441699948447043,
          "points": 69,
          "xy_norm_px": 0.0010878244186414653
        },
        {
          "cycle": -3,
          "btjd": 2979.8547995899085,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.00001186319314349987,
            "DET_FLUX": -0.000050744294244166766,
            "DET_FLUX_SML": -0.00004459442611930895,
            "DET_FLUX_LAG": -0.00007621748156463241,
            "SAP_X": -0.0003064461582198419,
            "SAP_Y": 0.001058329056149887,
            "SAP_BKG": 5.3738108249444565
          },
          "sap_nominal_score": -0.0177787231795839,
          "points": 69,
          "xy_norm_px": 0.0011018029038711104
        },
        {
          "cycle": -2,
          "btjd": 2981.780478007012,
          "selected": true,
          "valid": true,
          "depths": {
            "SAP_FLUX": 0.009026228457885768,
            "DET_FLUX": 0.005336961410594736,
            "DET_FLUX_SML": 0.003002265530198074,
            "DET_FLUX_LAG": 0.00904604777765626,
            "SAP_X": 0.009021343382812574,
            "SAP_Y": 0.0024007599752167152,
            "SAP_BKG": 7.990709473813865
          },
          "sap_nominal_score": 11.025657678743718,
          "points": 51,
          "xy_norm_px": 0.00933532457331928
        },
        {
          "cycle": -1,
          "btjd": 2983.706156424115,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": 0.0004975828700557407,
            "DET_FLUX": 0.0002286151111269148,
            "DET_FLUX_SML": 0.00031623225495208576,
            "DET_FLUX_LAG": 0.0003379920918458784,
            "SAP_X": 0.0004388091923601678,
            "SAP_Y": -0.00018616717193253862,
            "SAP_BKG": 11.907863367956267
          },
          "sap_nominal_score": 0.7234996170377996,
          "points": 70,
          "xy_norm_px": 0.00047666730872291015
        },
        {
          "cycle": 0,
          "btjd": 2985.631834841218,
          "selected": true,
          "valid": true,
          "depths": {
            "SAP_FLUX": 0.010229010770575052,
            "DET_FLUX": 0.006408824866572115,
            "DET_FLUX_SML": 0.0037813788898780508,
            "DET_FLUX_LAG": 0.010073126537595522,
            "SAP_X": 0.008017144285986433,
            "SAP_Y": 0.0029013237312938145,
            "SAP_BKG": -23.490317655323196
          },
          "sap_nominal_score": 16.9558150266007,
          "points": 69,
          "xy_norm_px": 0.008525976888081129
        },
        {
          "cycle": 1,
          "btjd": 2987.557513258321,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.0003652678973767899,
            "DET_FLUX": 0.00044732856884548014,
            "DET_FLUX_SML": -0.000027817120479510613,
            "DET_FLUX_LAG": 0.00043643868356630657,
            "SAP_X": -0.002387389852344768,
            "SAP_Y": -0.0003273447916853911,
            "SAP_BKG": 240.95070973287324
          },
          "sap_nominal_score": -0.4037004119374332,
          "points": 68,
          "xy_norm_px": 0.00240972714632224
        }
      ]
    },
    "62": {
      "file": "hlsp_qlp_tess_ffi_s0062-0000000140920261_tess_v01_llc.fits",
      "sha256": "5d7b63c09667fab29a8bf92ef9bf34fd93abff36a5f2948f6e3411e60754904b",
      "cadences": 10038,
      "windows": [
        {
          "cycle": 2,
          "btjd": 2989.4831916754238,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": 0.0007752640505160332,
            "DET_FLUX": -0.00015840914622071835,
            "DET_FLUX_SML": 0.00039943618400982526,
            "DET_FLUX_LAG": -0.00022031959298927142,
            "SAP_X": -0.003287174438042939,
            "SAP_Y": -0.0019156679577463365,
            "SAP_BKG": 398.74749879694355
          },
          "sap_nominal_score": 1.0687436743847873,
          "points": 69,
          "xy_norm_px": 0.003804641837342765
        },
        {
          "cycle": 3,
          "btjd": 2991.408870092527,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.0007021730217173679,
            "DET_FLUX": -0.0004160074381771309,
            "DET_FLUX_SML": -0.0001045585263392205,
            "DET_FLUX_LAG": -0.0008512270431357525,
            "SAP_X": -0.0006981411462163578,
            "SAP_Y": 0.00032073203442450904,
            "SAP_BKG": 30.880497242398803
          },
          "sap_nominal_score": -1.096078942348867,
          "points": 69,
          "xy_norm_px": 0.0007682903734567904
        },
        {
          "cycle": 4,
          "btjd": 2993.33454850963,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.0002745806845459416,
            "DET_FLUX": -0.0006433306154048953,
            "DET_FLUX_SML": 0.00001286581155030323,
            "DET_FLUX_LAG": -0.00115483415667097,
            "SAP_X": -0.000280362221502628,
            "SAP_Y": 0.00018286155182067887,
            "SAP_BKG": -4.426339490321496
          },
          "sap_nominal_score": -0.3748998943268139,
          "points": 69,
          "xy_norm_px": 0.00033472574203391567
        },
        {
          "cycle": 5,
          "btjd": 2995.260226926733,
          "selected": false,
          "valid": false
        },
        {
          "cycle": 6,
          "btjd": 2997.185905343836,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.00017462816053123308,
            "DET_FLUX": 0.0007466042434680774,
            "DET_FLUX_SML": -0.0001019954120968675,
            "DET_FLUX_LAG": 0.0011251891900976156,
            "SAP_X": -0.00013226224274376364,
            "SAP_Y": 0.000563973905855164,
            "SAP_BKG": 12.909023504047454
          },
          "sap_nominal_score": -0.25553971989646485,
          "points": 69,
          "xy_norm_px": 0.0005792752949514934
        },
        {
          "cycle": 7,
          "btjd": 2999.111583760939,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.0004626194628888666,
            "DET_FLUX": 0.00036106197734054127,
            "DET_FLUX_SML": -0.00026755080131993605,
            "DET_FLUX_LAG": 0.0005067331543448706,
            "SAP_X": 0.00007397918784590729,
            "SAP_Y": -0.0003989832780462166,
            "SAP_BKG": -7.454602850316931
          },
          "sap_nominal_score": -0.646949535664435,
          "points": 69,
          "xy_norm_px": 0.0004057839035679516
        },
        {
          "cycle": 8,
          "btjd": 3001.0372621780425,
          "selected": false,
          "valid": false
        },
        {
          "cycle": 9,
          "btjd": 3002.9629405951455,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.00012833957111607792,
            "DET_FLUX": -0.00034978103957168096,
            "DET_FLUX_SML": 0.000020955155968552084,
            "DET_FLUX_LAG": -0.0004253531125895911,
            "SAP_X": -0.00039010532324823544,
            "SAP_Y": 0.0005342704175745894,
            "SAP_BKG": -0.9418917118945274
          },
          "sap_nominal_score": -0.16861544409194618,
          "points": 69,
          "xy_norm_px": 0.0006615338557639634
        },
        {
          "cycle": 10,
          "btjd": 3004.8886190122485,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.0007340822262463,
            "DET_FLUX": -0.001090379821898066,
            "DET_FLUX_SML": -0.0005462774120149872,
            "DET_FLUX_LAG": -0.001643235411071385,
            "SAP_X": 0.00043044283773951975,
            "SAP_Y": 0.0011291406578948227,
            "SAP_BKG": -8.012310661612588
          },
          "sap_nominal_score": -1.2657413952772105,
          "points": 69,
          "xy_norm_px": 0.0012084037660783763
        },
        {
          "cycle": 11,
          "btjd": 3006.8142974293514,
          "selected": true,
          "valid": true,
          "depths": {
            "SAP_FLUX": 0.010188862439119659,
            "DET_FLUX": 0.0059608555923919854,
            "DET_FLUX_SML": 0.002897012657950565,
            "DET_FLUX_LAG": 0.010259521628235194,
            "SAP_X": 0.0071883036877693485,
            "SAP_Y": 0.0064976518408457196,
            "SAP_BKG": -23.01437726664031
          },
          "sap_nominal_score": 14.65412965044718,
          "points": 61,
          "xy_norm_px": 0.009689746609300172
        },
        {
          "cycle": 12,
          "btjd": 3008.7399758464544,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": -0.0008672866678698201,
            "DET_FLUX": -0.0007567829582896372,
            "DET_FLUX_SML": -0.0005655564689755028,
            "DET_FLUX_LAG": -0.0010480013206708252,
            "SAP_X": -0.0005109063799295654,
            "SAP_Y": -0.00020653109118029533,
            "SAP_BKG": 20.83381682515389
          },
          "sap_nominal_score": -1.249704110317449,
          "points": 69,
          "xy_norm_px": 0.0005510720648670707
        },
        {
          "cycle": 13,
          "btjd": 3010.665654263558,
          "selected": true,
          "valid": true,
          "depths": {
            "SAP_FLUX": 0.010622912467952005,
            "DET_FLUX": 0.005878063797200067,
            "DET_FLUX_SML": 0.0035401134110987806,
            "DET_FLUX_LAG": 0.009943471558080829,
            "SAP_X": 0.006955560009544115,
            "SAP_Y": 0.0062988832355586055,
            "SAP_BKG": -71.37449756824411
          },
          "sap_nominal_score": 13.797387255379528,
          "points": 69,
          "xy_norm_px": 0.009383802271018427
        },
        {
          "cycle": 14,
          "btjd": 3012.591332680661,
          "selected": false,
          "valid": true,
          "depths": {
            "SAP_FLUX": 0.00026686855283467297,
            "DET_FLUX": -0.00015705408141440436,
            "DET_FLUX_SML": -0.00032415012514070015,
            "DET_FLUX_LAG": -0.00006101418526899022,
            "SAP_X": 0.00003412058450180666,
            "SAP_Y": 0.00001357377050226205,
            "SAP_BKG": -2.567383419155476
          },
          "sap_nominal_score": 0.31483542488598093,
          "points": 69,
          "xy_norm_px": 0.00003672140428133173
        }
      ]
    }
  },
  "seconds": 0.09800360001099762
}
```

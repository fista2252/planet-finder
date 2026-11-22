# TIC 140920261: Sector 9 position and pointing comparison

Disposition remains **unresolved**. The QLP-only Sector 9 dip coincides with an unusually large change in the QLP position columns. The TESS-SPOC flux-weighted centroid also shifts in the same detector quadrant but by a much smaller amount. The TESS-SPOC position correction itself is ordinary relative to local controls. This comparison does not localize the signal or rule out image motion or pipeline systematics.

## Choice and protocol

The Sector 62 pixels remain the primary next observation. The [previous checkpoint](2026-11-21-checkpoint.md) records why acquisition was blocked. Today we asked a narrower question using the already cached Sector 9 FITS: does the position change at the disputed QLP dip appear in the alternate reduction, and does it track the alternate pipeline's pointing correction? These are new position measurements, using the same original exposures and [the frozen Sector 9 centers](sector9-reduction-check.json), not an independent observation. The cached TESS-SPOC observation identity also explicitly gives a Sector 9 target-pixel product ID, `hlsp_tess-spoc_tess_phot_0000000140920261-s0009_tess_v1_tp`. Its pixel file was not fetched. A web-only Gaia DR3 45-arcsecond CDS query was inaccessible, so no neighbor crossmatch is claimed. The cache contained no second QLP target in Sectors 9, 61 or 62 for a contemporaneous field-star control.

Use quality-zero finite timestamps matched one-to-one within 0.002 days, with TDB BJD reference 2457000. Compare signed sideband-minus-event offsets with a fixed 0.16-day duration. Use the existing linear sideband scorer at the two frozen times and on the comparison grid used for the previous aperture checks: centers every 0.64 d from the first retained time plus 0.32 d, excluding centers within 0.64 d of either frozen event. The second frozen center is the sensitive-miss comparison, not an independent random control. Keep each pair of columns on common finite cadences and report all available pairs. A magnitude outside the sampled control range is descriptive only. See the [original QLP extraction method](https://research.usq.edu.au/download/7f686e690658d98f289ed50a75204b6365664979c03301e05ec31c705897de7c/919161/Huang_2020_Res._Notes_AAS_4_204.pdf) and [TESS products description](https://ntrs.nasa.gov/api/citations/20205008729/downloads/EXP-TESS-ARC-ICD-TM-0014-Rev-F_v2.pdf?attachment=true). QLP computes per-frame astrometric positions for extraction. These files did not establish that SAP_X and SAP_Y are target flux-weighted centroids, so their amplitudes should not be interpreted as equivalent to the SPOC moment centroids or used to triangulate a blend.

## Results

952 matched cadences, maximum time separation 1.97 seconds. Position changes are magnitudes in pixels. Each usable pair had 28 valid control windows.

| Series | QLP dip, BTJD 1547.15006 | Sensitive miss, BTJD 1564.48116 | Controls at least as large as dip |
| --- | ---: | ---: | ---: |
| QLP SAP_X/Y | 0.011036 | 0.000254 | 0/28 |
| SPOC moment centroid | 0.000793 | 0.000136 | 0/28 |
| SPOC position correction | 0.000190 | 0.000274 | 20/28 |
| SPOC moment minus position correction | 0.000716 | 0.000234 | 2/28 |

The PSF-derived SPOC centroids had no finite pairs. The simple subtraction is an exploratory residual, not a fitted correction or a false-alarm calculation. The QLP flux feature is absent from matched-cadence SPOC flux as previously reported. The difference between position measures is consistent with their different estimators or apertures and is insufficient to distinguish a blend from local reduction or image effects. FITS primary headers verify that Sectors 61 and 62 use camera 3, CCDs 1 and 2 respectively. This weakens explanations specific to a single CCD without excluding instrument-wide pointing or shared reduction behavior.

Next, compare the two frozen Sector 62 event deficit maps to local and pointing-matched controls after checking WCS and time metadata. The Sector 9 archive pixel product is a secondary route for comparing the two apertures at the disputed event.

## Reproduction

Local checkout `25512680d1e4a9ac2fda85eb6836367f88bf8e61`. Run the code from that checkout with the shared cache and compare SHA-256 values to the archived [test receipts](result.json) and [SPOC receipt](sector9-reduction-check.json). No newly downloaded FITS bytes, local edits, threshold changes, or paid services. The measured computation took 0.17 s in-process.

```python
import hashlib, json, sys, time
from pathlib import Path
import astropy, numpy as np
from astropy.io import fits
sys.path.insert(0, 'legacy')
from screen import scores

start = time.perf_counter()
cache = Path('.audit-cache/shared/fits')
patterns = ['*qlp*s0009*140920261*.fits', '*tess-spoc*140920261*s0009*.fits']
paths = [next(cache.glob('*/' + pattern)) for pattern in patterns]
centers = np.array([1547.15006, 1564.48116])
d = 0.16
def read(p):
    with fits.open(p, memmap=False) as h:
        hd = h[1].header
        assert hd['TIMESYS'] == 'TDB' and hd['BJDREFI'] == 2457000 and hd.get('BJDREFF', 0) == 0
        a = h[1].data
        a = a[(a['QUALITY'] == 0) & np.isfinite(a['TIME'])]
        return a[np.argsort(a['TIME'])]
q,s = map(read, paths)
qt, st = q['TIME'], s['TIME']
j = np.clip(np.searchsorted(st, qt), 1, len(st)-1)
j -= (abs(st[j-1]-qt) < abs(st[j]-qt))
ok = abs(st[j]-qt) <= .002
assert np.unique(j[ok]).size == ok.sum()
q,s = q[ok],s[j[ok]]
t = q['TIME'].astype(float)
grid = np.arange(t[0]+2*d, t[-1]-2*d, 4*d)
grid = grid[np.min(abs(grid[:,None]-centers), axis=1)>4*d]
def matrix(a, keys):
    return np.column_stack([a[k].astype(float) for k in keys])
qlp = matrix(q, ['SAP_X','SAP_Y'])
mom = matrix(s, ['MOM_CENTR1','MOM_CENTR2'])
pos = matrix(s, ['POS_CORR1','POS_CORR2'])
psf = matrix(s, ['PSF_CENTR1','PSF_CENTR2'])
def measure(a):
    valid = np.isfinite(a).all(axis=1)
    if valid.sum() < 3:
        return dict(finite=int(valid.sum()), events=[])
    def vectors(w):
        x,y = [scores(t[valid],z[valid],w,d) for z in a.T]
        return {m['tc']:[m['depth'],n['depth']] for m,n in zip(x,y)}
    ctl = vectors(grid)
    norm = np.linalg.norm(list(ctl.values()),axis=1) if ctl else np.array([])
    return dict(finite=int(valid.sum()), controls=len(ctl),
        control_median=float(np.median(norm)) if ctl else None,
        control_max=float(max(norm)) if ctl else None,
        events=[dict(btjd=c, xy=v, norm=float(np.linalg.norm(v)),
                      controls_ge=int(np.count_nonzero(norm>=np.linalg.norm(v))))
                for c,v in vectors(centers).items()])
out = dict(python=sys.version.split()[0], numpy=np.__version__, astropy=astropy.__version__,
           files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
           matched=int(ok.sum()), max_time_diff_d=float(max(abs(t-s['TIME']))),
           camera_ccd={str(sector):[fits.getval(next(cache.glob(f'*/hlsp_qlp_tess_ffi_s{sector:04d}-0000000140920261_tess_v01_llc.fits')),key) for key in ('CAMERA','CCD')] for sector in (61,62)},
           measured={key:measure(a) for key,a in
                     dict(qlp=qlp, spoc_moment=mom, spoc_pointing=pos,
                          spoc_moment_minus_pointing=mom-pos, spoc_psf=psf).items()})
out['seconds'] = time.perf_counter()-start
print(json.dumps(out, indent=2))
```

```json
{
  "python": "3.13.15",
  "numpy": "2.5.4",
  "astropy": "8.0.2",
  "files": {
    "hlsp_qlp_tess_ffi_s0009-0000000140920261_tess_v01_llc.fits": "eb044137333e1efae5b7bab02f36c408c99e2da07a3020d773f933de1b22d6b3",
    "hlsp_tess-spoc_tess_phot_0000000140920261-s0009_tess_v1_lc.fits": "fbb0f8f7b708fbbebcd54df733ee7de582e243fdc6ed27cbc509615298c40a0b"
  },
  "matched": 952,
  "max_time_diff_d": 0.000022701654643242364,
  "camera_ccd": {
    "61": [
      3,
      1
    ],
    "62": [
      3,
      2
    ]
  },
  "measured": {
    "qlp": {
      "finite": 952,
      "controls": 28,
      "control_median": 0.00046103152249844033,
      "control_max": 0.004153418370545099,
      "events": [
        {
          "btjd": 1547.15006,
          "xy": [
            0.005834791241341009,
            0.009367201187842511
          ],
          "norm": 0.011035816549922707,
          "controls_ge": 0
        },
        {
          "btjd": 1564.48116,
          "xy": [
            0.00024244513360827114,
            0.00007629418018950673
          ],
          "norm": 0.00025416617544654,
          "controls_ge": 22
        }
      ]
    },
    "spoc_moment": {
      "finite": 952,
      "controls": 28,
      "control_median": 0.0002676126287568955,
      "control_max": 0.0006508880934267766,
      "events": [
        {
          "btjd": 1547.15006,
          "xy": [
            0.0004267320216513326,
            0.0006685135772386275
          ],
          "norm": 0.0007931018984058857,
          "controls_ge": 0
        },
        {
          "btjd": 1564.48116,
          "xy": [
            0.00011534427767401212,
            0.00007216812537080841
          ],
          "norm": 0.0001360607978503593,
          "controls_ge": 24
        }
      ]
    },
    "spoc_pointing": {
      "finite": 952,
      "controls": 28,
      "control_median": 0.00029639162701352917,
      "control_max": 0.00093833256316422,
      "events": [
        {
          "btjd": 1547.15006,
          "xy": [
            0.00019027308160872447,
            -0.0000071232058316827995
          ],
          "norm": 0.00019040636976267586,
          "controls_ge": 20
        },
        {
          "btjd": 1564.48116,
          "xy": [
            0.0002449749138955237,
            -0.0001229666245534458
          ],
          "norm": 0.00027410490545079155,
          "controls_ge": 15
        }
      ]
    },
    "spoc_moment_minus_pointing": {
      "finite": 952,
      "controls": 28,
      "control_median": 0.0002918811134186691,
      "control_max": 0.0008506769719827237,
      "events": [
        {
          "btjd": 1547.15006,
          "xy": [
            0.00023645893998036627,
            0.0006756367831712851
          ],
          "norm": 0.0007158197350385643,
          "controls_ge": 2
        },
        {
          "btjd": 1564.48116,
          "xy": [
            -0.00012963063636561856,
            0.0001951347499016265
          ],
          "norm": 0.000234268377088598,
          "controls_ge": 19
        }
      ]
    },
    "spoc_psf": {
      "finite": 0,
      "events": []
    }
  },
  "seconds": 0.1696049000020139
}
```

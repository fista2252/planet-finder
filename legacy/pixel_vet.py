"""Compute local difference images at the fixed epochs of three sparse repeat cases."""
import io
import json
import zipfile

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from screen import ROOT

EVENTS = (("339672790", 7, 1492.3037435589358, .3044),
          ("339672790", 33, 2211.004655372317, .3044),
          ("273607159", 27, 2056.8050497034546, .18555),
          ("271693283", 64, 3054.3236873269766, .160425))
EPHEMERIS_SECTORS = {"389667951": (1,96), "269730302": (8,68), "238187210": (1,66),
                     "150108429": (63,), "340562993": (61,)}


def difference(t, cube, epoch, duration, shift=0):
    rel = (t - epoch) / duration
    a = (rel >= shift-.5) & (rel < shift+.5)
    before = (rel >= shift-2) & (rel < shift-.75)
    after = (rel >= shift+.75) & (rel < shift+2)
    if min(a.sum(), before.sum(), after.sum()) < 3:
        return None
    weight = (np.mean(t[a]) - np.mean(t[before])) / (np.mean(t[after])-np.mean(t[before]))
    baseline = (1-weight)*np.median(cube[before], axis=0) + weight*np.median(cube[after],axis=0)
    missing = baseline - np.median(cube[a], axis=0)
    missing -= np.median(np.r_[missing[0], missing[-1], missing[1:-1, 0], missing[1:-1, -1]])
    return missing, [int(x.sum()) for x in (before,a,after)], baseline


def analyze(tic, sector, epoch, duration):
    path = ROOT / "data/cutouts" / f"{tic}-s{sector:04}.zip"
    out = dict(tic=tic, sector=sector, epoch=epoch, duration=duration)
    if not path.exists():
        return dict(**out, error="cutout unavailable")
    nearby = json.loads((ROOT / "data/nearby" / f"{tic}.json").read_text())
    stars = sorted((s for s in nearby if s["Tmag"] is not None and s["dstArcSec"] <= 42),
                   key=lambda s: (str(s["ID"]) != tic, s["Tmag"]))
    with zipfile.ZipFile(path) as z:
        member = next(n for n in z.namelist() if n.endswith(".fits"))
        with fits.open(io.BytesIO(z.read(member))) as f:
            assert f[1].header["BJDREFI"] == 2457000 and f[1].header["TIMESYS"] == "TDB"
            table = f[1].data
            t = table["TIME"].astype(float) + float(f[1].header.get("BJDREFF",0))
            flux = table["FLUX"].astype(float)
            valid = np.isfinite(t) & np.isfinite(flux).all(axis=(1,2)) & (table["QUALITY"] == 0)
            wcs = WCS(f[2].header).celestial
            x,y = wcs.world_to_pixel_values([s["ra"] for s in stars], [s["dec"] for s in stars])
            out["catalog_pixels"] = [dict(tic=str(s["ID"]), Tmag=s["Tmag"], x=float(xx), y=float(yy))
                                     for s,xx,yy in zip(stars,x,y)]
    t, flux = t[valid], flux[valid]
    event = difference(t, flux, epoch, duration)
    if event is None:
        return dict(**out,error="inadequate local quality-zero pixel coverage")
    missing, counts, baseline = event
    peak_y, peak_x = np.unravel_index(np.argmax(missing), missing.shape)
    yy, xx = np.indices(missing.shape)
    patch = (abs(xx-peak_x) <= 1) & (abs(yy-peak_y) <= 1)
    weights = np.maximum(0, missing) * patch
    centroid = [float((xx*weights).sum()/weights.sum()), float((yy*weights).sum()/weights.sum())] if weights.sum() else None
    control = [difference(t, flux, epoch, duration, shift) for shift in (-9,-6,-3,3,6,9)]
    null_peaks = [float(abs(a[0]).max()) for a in control if a is not None]
    return dict(**out, n_before_in_after=counts, missing=missing.tolist(), peak_pixel_xy=[int(peak_x),int(peak_y)],
                peak_flux=float(missing.max()), centroid_xy=centroid, null_max_abs_pixel=null_peaks,
                local_aperture_depth=float(missing[patch].sum()/baseline[patch].sum()))


if __name__ == "__main__":
    events = list(EVENTS)
    for target in json.loads((ROOT / "ephemeris_vetting.json").read_text()):
        for sector in EPHEMERIS_SECTORS[target["tic"]]:
            candidates = [m for m in target["measurements"] if m["sector"] == sector and m["phase"] == 0]
            if candidates:
                match = max(candidates, key=lambda m:m["score"])
                events.append((target["tic"], sector, match["tc"], match["duration"]))
    results = [analyze(*args) for args in events]
    (ROOT / "pixel_vetting.json").write_text(json.dumps(results,indent=2),encoding="utf-8")
    for row in results:
        print(row["tic"], row["sector"], row.get("centroid_xy"), row.get("local_aperture_depth"), row.get("error"))

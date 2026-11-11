"""Fixed-ephemeris measurements on unused sectors and TGLC PSF photometry."""
import io
import json
import zipfile

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from fit_transits import HERE, ROOT, REPEAT, PAIR, event, template
from diagnostics import centers


def depth(t,y,epoch,d,period,k,b):
    e=event(t,y,epoch,d,0,0,"extra")
    if e is None:
        return None
    x,y=e["x"],e["y"]
    shape=template(x,period,0,d,k,b,e["cadence"])
    matrix=np.column_stack((np.ones(len(x)),x,shape))
    a=np.linalg.lstsq(matrix,y,rcond=None)[0]
    resid=y-matrix@a
    sigma=max(e["sigma"],np.std(resid)*np.sqrt(max(1,(20/1440)/e["cadence"])))
    dep=-shape.min()*a[-1]
    error=-shape.min()*sigma*np.sqrt(np.linalg.pinv(matrix.T@matrix)[-1,-1])
    return dict(depth=float(dep),error=float(error),nominal_score=float(dep/error),
                in_template=int(np.count_nonzero(shape<0)),points=len(x))


def pixel_curve(path,ra,dec):
    with zipfile.ZipFile(path) as archive:
        member=next(n for n in archive.namelist() if n.endswith(".fits"))
        with fits.open(io.BytesIO(archive.read(member))) as h:
            hdr=h[1].header
            assert hdr["BJDREFI"]==2457000 and hdr["TIMESYS"]=="TDB"
            tab=h[1].data
            t=tab["TIME"]+float(hdr.get("BJDREFF",0))
            c=tab["FLUX"].astype(float)
            x,y=WCS(h[2].header).celestial.world_to_pixel_values(ra,dec)
            yy,xx=np.indices(c.shape[1:])
            aperture=(abs(xx-round(float(x)))<=1)&(abs(yy-round(float(y)))<=1)
            perimeter=(yy==0)|(xx==0)|(yy==yy.max())|(xx==xx.max())
            a=c[:,aperture].sum(axis=1)-aperture.sum()*np.nanmedian(c[:,perimeter],axis=1)
            return t.copy(),a,tab["QUALITY"].copy()


def prepared(t,y,q):
    take=np.isfinite(t)&np.isfinite(y)&(q==0)
    t,ix=np.unique(t[take],return_index=True)
    flux=y[take][ix]
    return t,flux/np.median(flux)-1


if __name__=="__main__":
    f=json.loads((HERE/"transit_fits.json").read_text())
    archive=HERE/"archive"
    manifest=json.loads((archive/"extra_manifest.json").read_text())
    cohort={a["catalog"]["TICID"]:a for a in json.loads((ROOT/"cohort.json").read_text())}
    out=[]
    for label in manifest["extra"]:
        if not label["file"].endswith(".zip"):
            continue
        path=archive/label["file"]
        tic=path.stem.split("_")[1]
        sector=int(path.stem.split("_")[-1][1:])
        p=cohort[tic]["primary"][0]
        t,flux,flags=pixel_curve(path,p["s_ra"],p["s_dec"])
        valid=np.isfinite(t)&np.isfinite(flux)&(flags==0)
        net=float(np.median(flux[valid]))
        # Compare against a known source-bearing cutout for the same TIC.
        refsector={"269730302":8,"238187210":1}[tic]
        refpath=ROOT/"data/cutouts"/f"{tic}-s{refsector:04}.zip"
        rt,rf,rq=pixel_curve(refpath,p["s_ra"],p["s_dec"])
        reference=float(np.median(rf[np.isfinite(rt)&np.isfinite(rf)&(rq==0)]))
        t=t[valid]
        q=f[tic]["nominal"]["q"]
        t0,p=REPEAT[tic][0]+q[0],REPEAT[tic][1]+q[1]
        _,pred=centers(t,t0,p)
        before=t0+np.floor((t[0]-t0)/p)*p
        after=t0+np.ceil((t[-1]-t0)/p)*p
        row=dict(tic=tic,sector=sector,product="unused TESScut",span=[float(t[0]),float(t[-1])],
                 nearest_outside=[float(before),float(after)],events=[],
                 median_net_aperture=net,reference_source_aperture=reference)
        if not reference>0 or net<.1*reference:
            row["error"]="No useful target flux in aperture compared with verified source-bearing cutout"
            row["predicted_epochs"]=pred.tolist()
            out.append(row)
            continue
        t,y=prepared(*pixel_curve(path,p["s_ra"],p["s_dec"]))
        for epoch in pred:
            row["events"].append(dict(epoch=float(epoch),result=depth(t,y,epoch,q[2],p,q[3],q[4])))
        out.append(row)
    for prod in manifest["products"]:
        tic=prod["target_name"].lstrip("0")
        if tic not in f:
            tic=next(a for a in f if cohort[a]["catalog"]["GAIAID"] in prod["dataURL"])
        path=ROOT/"data/curves"/prod["dataURL"].rsplit("/",1)[1]
        q=f[tic]["nominal"]["q"]
        if tic in REPEAT:
            t0,p=REPEAT[tic][0]+q[0], REPEAT[tic][1]+q[1]
        else:
            t0,p=PAIR[0]+q[0],239.56837
        with fits.open(path) as h:
            a=h[1].data
            for name in ("cal_psf_flux","cal_aper_flux"):
                t,y=prepared(a["time"].astype(float),a[name].astype(float),a["TESS_flags"]|a["TGLC_flags"])
                out.append(dict(tic=tic,sector=prod["sequence_number"],product="TGLC "+name,
                                result=depth(t,y,t0,q[2],p,q[3],q[4])))
    (HERE/"unused_observations.json").write_text(json.dumps(out,indent=2),encoding="utf-8")

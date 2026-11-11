"""Block-resample pixel difference centroids against plausible TIC neighbors."""
import io
import json
import zipfile
from pathlib import Path

import numpy as np
from astropy.io import fits

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

RNG=np.random.default_rng(31110)


def sample(a,n):
    block=max(1,min(n,len(a)))
    starts=RNG.integers(0,len(a)-block+1,size=int(np.ceil(len(a)/block)))
    return a[(starts[:,None]+np.arange(block)).ravel()[:len(a)]]


def calculate(row):
    path=ROOT/"data/cutouts"/f'{row["tic"]}-s{row["sector"]:04}.zip'
    with zipfile.ZipFile(path) as z:
        name=next(n for n in z.namelist() if n.endswith(".fits"))
        with fits.open(io.BytesIO(z.read(name))) as h:
            a=h[1].data
            t=a["TIME"].astype(float)
            y=a["FLUX"].astype(float)
            ok=(a["QUALITY"]==0)&np.isfinite(t)&np.isfinite(y).all(axis=(1,2))
    t,y=t[ok],y[ok]
    dt=(t-row["epoch"])/row["duration"]
    groups=[y[(dt>=lo)&(dt<hi)] for lo,hi in ((-2,-.75),(-.5,.5),(.75,2))]
    block=max(1,round((30/1440)/np.median(np.diff(t))))
    grid_y,grid_x=np.indices(y.shape[1:])
    target=next(x for x in row["catalog_pixels"] if x["tic"]==row["tic"])
    tx,ty=target["x"],target["y"]
    center_mask=(abs(grid_x-round(tx))<=1)&(abs(grid_y-round(ty))<=1)
    perimeter=(grid_x==0)|(grid_y==0)|(grid_x==grid_x.max())|(grid_y==grid_y.max())
    median_depth=row["local_aperture_depth"]
    plausible=[s for s in row["catalog_pixels"] if s["tic"]!=row["tic"]
               and 10**(-.4*(s["Tmag"]-target["Tmag"]))/(1+10**(-.4*(s["Tmag"]-target["Tmag"])))>=median_depth]
    xy=[]
    for _ in range(500):
        left,inside,right=(np.median(sample(g,block),axis=0) for g in groups)
        z=(left+right)/2-inside
        z=np.maximum(0,z-np.median(z[perimeter]))*center_mask
        if z.sum()>0:
            xy.append([(grid_x*z).sum()/z.sum(),(grid_y*z).sum()/z.sum()])
    xy=np.asarray(xy)
    dist=np.linalg.norm(xy-[tx,ty],axis=1)
    comparisons=[]
    for star in plausible:
        other=np.linalg.norm(xy-[star["x"],star["y"]],axis=1)
        comparisons.append(dict(tic=star["tic"],separation_pixels=float(np.hypot(star["x"]-tx,star["y"]-ty)),
                                target_closer_fraction=float(np.mean(dist<other))))
    return dict(tic=row["tic"],sector=row["sector"],draws=len(xy),block=block,
                target_distance_pixels_95=np.quantile(dist,[.025,.5,.975]).tolist(),
                plausible_neighbors=comparisons)


if __name__=="__main__":
    rows=json.loads((ROOT/"pixel_vetting.json").read_text())
    out=[calculate(r) for r in rows if r["tic"] in {"269730302","238187210","339672790"} and r.get("centroid_xy")]
    (HERE/"centroid_check.json").write_text(json.dumps(out,indent=2))

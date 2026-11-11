"""Audit TIC 238187210 using only the saved light curves."""
import json
import sys
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
from screen import read_curve

TIC="238187210"
RNG=np.random.default_rng(238187210)


def residualize(matrix,x):
    A=np.column_stack([np.ones(len(x)),x])
    return matrix-A@np.linalg.lstsq(A,matrix,rcond=None)[0]


def fit_window(t,y,epoch,duration):
    x=t-epoch
    take=abs(x)<.55
    x,y=x[take],y[take]
    if len(x)<12:
        return dict(state="insufficient local cadences")
    cad=float(np.median(np.diff(x)))
    counts=[int(mask.sum()) for mask in (x<-.75*duration,abs(x)<.5*duration,x>.75*duration)]
    if min(counts)<3:
        return dict(state="insufficient in-window or baseline cadences",counts=counts)
    offsets=np.linspace(-.05,.05,101)
    ex=cad*np.linspace(-.5,.5,7)
    trials=[]
    for d in np.linspace(.16,.30,36):
        # Fixed fractional ingress; independent of the original limb-darkened fit.
        shape=-np.clip((d/2-abs(x[:,None,None]-offsets[None,:,None]+ex))/(.22*d),0,1).mean(2)
        z=residualize(shape,x)
        den=np.sum(z*z,axis=0)
        yy=residualize(y,x)
        amp=(z.T@yy)/np.maximum(den,1e-12)
        rss=np.sum(yy**2)-np.maximum(amp,0)**2*den
        trials.append((rss,d,amp,z,den))
    rss,d,amp,z,den,k=min(((*v,int(np.argmin(v[0]))) for v in trials),key=lambda v:v[0][v[-1]])
    a=float(amp[k])
    e=yy-a*z[:,k]
    block=max(1,round((30/1440)/cad))
    scores=np.array([np.dot(z[i:i+block,k],e[i:i+block]) for i in range(0,len(e),block)])
    err=max(np.sqrt(np.sum(scores**2))/den[k],np.std(e)/np.sqrt(den[k]))
    return dict(state="measured",epoch=float(epoch+offsets[k]),duration=d,depth=a,error=float(err),
                score=float(a/err),counts=counts,offset_edge=bool(k in (0,len(offsets)-1)),
                shape_edge=bool(d in (.16,.30)))


def harmonic(curves,t0,p,raw):
    matrices=[]
    rhs=[]
    for prod in curves:
        t,y,_=read_curve(prod,raw=raw)
        phase=(t-t0+p/2)%p-p/2
        out=abs(phase)>.20
        t,y,phase=t[out],y[out],phase[out]
        x=(t-t.mean())/max(np.ptp(t),.01)
        drift=np.column_stack((np.ones(len(t)),x,x*x))
        y=y-drift@np.linalg.lstsq(drift,y,rcond=None)[0]
        scale=max(1.4826*np.median(abs(y-np.median(y))),1e-6)
        keep=abs(y)<5*scale
        # Uniform 1-hour bins prevent high-cadence sectors from dominating.
        bins=np.floor(24*t[keep]).astype(int)
        _,ix,n=np.unique(bins,return_inverse=True,return_counts=True)
        def bin(v):
            return np.bincount(ix,weights=v[keep])/n
        phase=2*np.pi*phase/p
        X=np.column_stack([bin(g) for g in (np.sin(phase),np.cos(phase),np.sin(2*phase),np.cos(2*phase))])
        D=np.column_stack((np.ones(len(n)),bin(x),bin(x*x)))
        yy=bin(y)
        X-=D@np.linalg.lstsq(D,X,rcond=None)[0]
        yy-=D@np.linalg.lstsq(D,yy,rcond=None)[0]
        matrices.append(X.T@X/scale**2)
        rhs.append(X.T@yy/scale**2)
    A,B=np.array(matrices),np.array(rhs)
    coef=np.linalg.pinv(A.sum(0))@B.sum(0)
    boot=[]
    for _ in range(1000):
        j=RNG.integers(0,len(A),len(A))
        boot.append(np.linalg.pinv(A[j].sum(0))@B[j].sum(0))
    boot=np.array(boot)
    return dict(n_sectors=len(A),coefficients=coef.tolist(),
                amplitude_ppm=[float(1e6*np.linalg.norm(coef[i:i+2])) for i in (0,2)],
                bootstrap_amplitude_95_percentile_ppm=[float(1e6*np.quantile(np.linalg.norm(boot[:,i:i+2],axis=1),.95)) for i in (0,2)])


if __name__=="__main__":
    item=next(r for r in json.loads((ROOT/"cohort.json").read_text()) if r["catalog"]["TICID"]==TIC)
    fit=json.loads((HERE/"diagnostics.json").read_text())[TIC]
    p,t0=fit["period"],fit["epoch"]
    products=[v for v in item["primary"] if v["sequence_number"]!=1 or v["provenance_name"]=="T16"]
    ledger=[]
    stacks=[]
    for product in products:
        t,y,_=read_curve(product)
        for cycle in range(int(np.ceil((t[0]-t0)/p)),int(np.floor((t[-1]-t0)/p))+1):
            ep=t0+cycle*p
            m=fit_window(t,y,ep,.236)
            ledger.append(dict(sector=product["sequence_number"],cycle=cycle,predicted=float(ep),**m))
            if m["state"]=="measured":
                x=t-ep
                take=abs(x)<.34
                side=take&(abs(x)>.19)
                baseline=np.polyfit(x[side],y[side],1)
                stacks.append(dict(sector=product["sequence_number"],cycle=cycle,
                                   hour=(24*x[take]).tolist(),relative=(y[take]-np.polyval(baseline,x[take])).tolist()))
    usable=[r for r in ledger if r.get("score",0)>=5 and not r.get("offset_edge",True)]
    n=np.array([r["cycle"] for r in usable])
    t=np.array([r["epoch"] for r in usable])
    periods=[]
    for _ in range(1000):
        ix=RNG.integers(0,len(n),len(n))
        if np.ptp(n[ix])>0:
            periods.append(np.polyfit(n[ix],t[ix]-t0,1)[0])
    core=dict(period=float(np.polyfit(n,t-t0,1)[0]),period_95=np.quantile(periods,[.025,.975]).tolist(),
              n_timing=len(n),depth=float(np.median([r["depth"] for r in usable])),
              duration_hours=float(24*np.median([r["duration"] for r in usable])),
              ledger=ledger,stacks=stacks,
              harmonics={k:harmonic([v for v in products if v["sequence_number"]!=1],t0,p,k=="raw")
                         for k in ("raw","detrended")})
    (HERE/"tic238_audit.json").write_text(json.dumps(core,indent=2),encoding="utf-8")

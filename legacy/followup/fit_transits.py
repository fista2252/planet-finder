"""Exposure-integrated, circular transit fits and event-cluster uncertainty."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(HERE/"deps"), str(ROOT)]
import setuptools  # Provides the distutils compatibility import used by batman.
import batman
import numpy as np
from scipy.optimize import least_squares
from screen import read_curve

REPEAT = {"269730302": (1520.40917,28.430676,.1116),
          "238187210": (1334.25317,26.027341,.16075)}
PAIR = (1492.30374, 2211.00466, .3044)


def template(t, period, epoch, duration, k, b, exposure, limb=(.35,.25)):
    p = batman.TransitParams()
    p.t0, p.per, p.rp = epoch, period, k
    p.a = np.sqrt(b*b+((1+k)**2-b*b)/np.sin(np.pi*duration/period)**2)
    p.inc, p.ecc, p.w = np.degrees(np.arccos(b/p.a)), 0, 90
    p.limb_dark, p.u = "quadratic", list(limb)
    return batman.TransitModel(p, t, supersample_factor=7, exp_time=exposure).light_curve(p)-1


def event(t,y,center,duration,cycle,sector,pipeline):
    x=t-center
    mask=abs(x)<max(.42,2.7*duration)
    x,y=x[mask],y[mask]
    if min(np.sum((x<-.85*duration)),np.sum((x>.85*duration)),np.sum(abs(x)<.47*duration))<3:
        return None
    side=abs(x)>.85*duration
    trend=np.polyfit(x[side],y[side],1)
    err=y[side]-np.polyval(trend,x[side])
    sigma=max(1.4826*np.median(abs(err-np.median(err))),1e-6)
    # One 20-minute block supplies no more weight than one independent datum.
    cadence=float(np.median(np.diff(x)))
    inflation=np.sqrt(max(1, (20/1440)/cadence))
    return dict(x=x,y=y,center=center,cycle=cycle,sector=sector,pipeline=pipeline,
                cadence=cadence,sigma=sigma*inflation)


def gather(item, repeated, raw=False, extra=()):
    tic=item["catalog"]["TICID"]
    orig=int(item["catalog"]["Sector"])
    if repeated:
        t0,p,d=REPEAT[tic]
    else:
        t0,t1,d=PAIR
        p=239.57
    products=[a for a in item["primary"] if a["sequence_number"]!=orig or a["provenance_name"]=="T16"]
    if not repeated:
        products=[a for a in products if a["sequence_number"] in (7,33)]
    windows=[]
    for a in products:
        t,y,_=read_curve(a,raw=raw)
        s=a["sequence_number"]
        if repeated:
            centers=[(n,t0+n*p) for n in range(int(np.ceil((t[0]-t0)/p)),int(np.floor((t[-1]-t0)/p))+1)]
        else:
            centers=[(int(s==33),t1 if s==33 else t0)]
        for n,center in centers:
            e=event(t,y,center,d,n,s,a["provenance_name"])
            if e is not None:
                windows.append(e)
    return windows,p,d


def baseline(x,y,shape):
    mat=np.column_stack((np.ones(len(x)), x))
    trend=np.linalg.lstsq(mat,y-shape,rcond=None)[0]
    return y-shape-mat@trend


def evaluate(windows,p,d,repeated,limb=(.35,.25)):
    def residual(q):
        dt,dr,duration,k,b=q
        return np.concatenate([
            baseline(e["x"],e["y"],template(e["x"],p+dr if repeated else p,
                        dt+e["cycle"]*dr if repeated else (dr if e["cycle"] else dt),
                        duration,k,b,e["cadence"],limb))/e["sigma"] for e in windows])
    low=[-.18,-.002 if repeated else -.18,.45*d,.015,0]
    high=[.18,.002 if repeated else .18,1.9*d,.45,.97]
    starts=([0,0,d,.085,b] for b in (.15,.5,.8))
    result=min((least_squares(residual,x0,bounds=(low,high),loss="soft_l1",f_scale=2,
                               x_scale=[.03,.0002 if repeated else .03,d,.04,.3],max_nfev=160)
                for x0 in starts), key=lambda m:m.cost)
    dt,dr,duration,k,b=result.x
    # Sandwich covariance clustered by transit epoch. This accounts for within-
    # event correlation, but not model misspecification or target blends.
    J=result.jac
    inv=np.linalg.pinv(J.T@J)
    offset=0
    outer=np.zeros((5,5))
    diagnostics=[]
    for e in windows:
        n=len(e["x"])
        q=result.fun[offset:offset+n]
        jac=J[offset:offset+n]
        v=jac.T@np.clip(q,-3,3)
        outer+=np.outer(v,v)
        offset+=n
        ep=dt+e["cycle"]*dr if repeated else (dr if e["cycle"] else dt)
        shape=template(e["x"],p+dr if repeated else p,ep,duration,k,b,e["cadence"],limb)
        design=np.column_stack((np.ones(n),e["x"],shape))
        beta=np.linalg.lstsq(design,e["y"],rcond=None)[0]
        resid=e["y"]-design@beta
        std=max(e["sigma"],np.std(resid)*np.sqrt(max(1,(20/1440)/e["cadence"])))
        amp_err=std*np.sqrt(np.linalg.pinv(design.T@design)[2,2])
        diagnostics.append(dict(sector=e["sector"],cycle=e["cycle"],center=e["center"]+ep,
                                depth=float(beta[2]*(-shape.min())), error=float(amp_err*(-shape.min())),
                                n=n,pipeline=e["pipeline"]))
    model_cov=inv*max(1,np.sum(np.minimum(result.fun**2,9))/(len(result.fun)-5))
    cov=inv@outer@inv*len(windows)/max(len(windows)-1,1)
    if repeated:
        # Add the model covariance as a floor for small-cluster fluctuations.
        cov=cov+model_cov
    else:
        cov=model_cov
    return dict(q=result.x.tolist(),covariance=cov.tolist(),sigma=np.sqrt(np.diag(cov)).tolist(),
                success=bool(result.success),cost=float(result.cost),events=diagnostics,
                limb=list(limb),duration=float(duration),radius_ratio=float(k),impact=float(b))


if __name__=="__main__":
    cohort={x["catalog"]["TICID"]:x for x in json.loads((ROOT/"cohort.json").read_text())}
    out={}
    for tic in (*REPEAT,"339672790"):
        repeated=tic in REPEAT
        windows,p,d=gather(cohort[tic],repeated)
        nominal=evaluate(windows,p,d,repeated)
        sensitivity=evaluate(windows,p,d,repeated,limb=(.55,.10))
        out[tic]=dict(nominal=nominal,limb_sensitivity=sensitivity,
                      reference_period=p,reference_epoch=REPEAT[tic][0] if repeated else PAIR[0],
                      reference_second_epoch=None if repeated else PAIR[1])
        print(tic,len(windows),nominal["q"],nominal["sigma"],flush=True)
    (HERE/"transit_fits.json").write_text(json.dumps(out,indent=2),encoding="utf-8")

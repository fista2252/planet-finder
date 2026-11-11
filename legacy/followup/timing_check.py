"""Estimate event-to-event timing scatter and conservative planning windows."""
import json
import sys
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
from screen import read_curve

REPEAT={"269730302": (1520.40917,28.430676,.1116),
        "238187210": (1334.25317,26.027341,.16075)}


def estimate(e,q,p):
    dt,dp,d,k,b=q
    ep=dt+e["cycle"]*dp
    x,y=e["x"],e["y"]
    def fit(offset):
        times=x[:,None]-(ep+offset)+e["cadence"]*np.linspace(-.5,.5,7)
        s=-np.clip((d/2-np.abs(times))/(.22*d),0,1).mean(axis=1)*k*k
        X=np.column_stack((np.ones(len(x)),x,s))
        a=np.linalg.lstsq(X,y,rcond=None)[0]
        resid=y-X@a
        loss=np.sum(np.hypot(1,resid/(2*e["sigma"]))-1)
        return loss, a[2]
    trials=np.linspace(-d*.32,d*.32,121)
    offset,(_,scale)=min(((v,fit(v)) for v in trials),key=lambda item:item[1][0])
    return dict(sector=e["sector"],cycle=e["cycle"],offset_days=float(offset),amplitude_scale=float(scale),
                edge=bool(abs(offset)>.30*d))


def gather(item,t0,p,d):
    original=int(item["catalog"]["Sector"])
    for prod in item["primary"]:
        sector=prod["sequence_number"]
        if sector==original and prod["provenance_name"]!="T16":
            continue
        t,y,_=read_curve(prod)
        for n in range(int(np.ceil((t[0]-t0)/p)),int(np.floor((t[-1]-t0)/p))+1):
            x=t-(t0+n*p)
            take=abs(x)<max(.42,2.7*d)
            a,b=x[take],y[take]
            if min(np.sum(a<-.85*d),np.sum(a>.85*d),np.sum(abs(a)<.47*d))<3:
                continue
            side=abs(a)>.85*d
            err=b[side]-np.polyval(np.polyfit(a[side],b[side],1),a[side])
            cadence=float(np.median(np.diff(a)))
            sigma=max(1.4826*np.median(abs(err-np.median(err))),1e-6)
            yield dict(x=a,y=b,cycle=n,sector=sector,cadence=cadence,
                       sigma=sigma*np.sqrt(max(1,(20/1440)/cadence)))


if __name__=="__main__":
    f=json.loads((HERE/"transit_fits.json").read_text())
    diag=json.loads((HERE/"diagnostics.json").read_text())
    cohort={a["catalog"]["TICID"]:a for a in json.loads((ROOT/"cohort.json").read_text())}
    out={}
    for tic, (t0,p,d) in REPEAT.items():
        windows=gather(cohort[tic],t0,p,d)
        q=f[tic]["nominal"]["q"]
        measures=[estimate(e,q,p+q[1]) for e in windows]
        usable=[m["offset_days"] for m in measures if .18<m["amplitude_scale"]<3 and not m["edge"]]
        jitter=float(np.sqrt(np.mean(np.square(usable)))) if usable else 20/1440
        floor=max(jitter,10/1440)
        next_events=[dict(**n,planning_halfwidth_minutes=float(np.hypot(n["timing_95_minutes"],1.96*1440*floor)))
                     for n in diag[tic]["next"]]
        out[tic]=dict(timing_rms_minutes=jitter*1440,valid_events=len(usable),events=measures,next=next_events)
    (HERE/"timing_check.json").write_text(json.dumps(out,indent=2),encoding="utf-8")

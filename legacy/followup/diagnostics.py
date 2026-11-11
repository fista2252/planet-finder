"""Odd/even, secondary-search, and cadence-conditioned alias injection tests."""
import json
from datetime import datetime

import numpy as np
from astropy.time import Time
from fit_transits import HERE,ROOT,REPEAT,PAIR,template
from screen import read_curve,scores

RNG=np.random.default_rng(91109)
TICS=(*REPEAT,"339672790")


def products(item):
    original=int(item["catalog"]["Sector"])
    return [p for p in item["primary"] if p["sequence_number"]!=original or p["provenance_name"]=="T16"]


def centers(t,t0,p):
    n=np.arange(np.ceil((t[0]-t0)/p),np.floor((t[-1]-t0)/p)+1,dtype=int)
    return n,t0+n*p


def parity(events):
    a,b=[np.array([e["depth"] for e in events if e["cycle"]%2==j]) for j in (0,1)]
    diff=np.median(RNG.choice(a,(4000,len(a))),axis=1)-np.median(RNG.choice(b,(4000,len(b))),axis=1)
    return dict(n=[len(a),len(b)],median_depth=[float(np.median(a)),float(np.median(b))],
                even_minus_odd_95=np.quantile(diff,[.025,.975]).tolist())


def secondary(curves,t0,p,d):
    phases=np.linspace(.08,.92,169)
    layers=[]
    half_depths=[]
    for t,y,_ in curves:
        _,ep=centers(np.array([t[0]-p,t[-1]]),t0,p)
        trial=(ep[:,None]+phases[None,:]*p).ravel()
        phase_scores=np.zeros(len(phases))
        n=np.zeros(len(phases),dtype=int)
        for q in scores(t,y,trial,d):
            ph=(q["tc"]-t0)%p/p
            j=np.argmin(abs(phases-ph))
            phase_scores[j]+=np.clip(q["score"],-6,6)
            n[j]+=1
        if n.any():
            layers.append((phase_scores,n))
        half_depths += scores(t,y,ep+.5*p,d)
    a,n=map(np.asarray,zip(*layers))
    total=n.sum(0)
    scale=np.sqrt(np.maximum(total,1))
    stats=a.sum(0)/scale
    valid=total>=4
    rank=np.where(valid,stats,-np.inf)
    best=int(np.argmax(rank))
    shuffled=(RNG.choice([-1,1],(4000,len(a)))@a)/scale
    fap=(1+np.sum(np.max(shuffled[:,valid],axis=1)>=rank[best]))/(1+len(shuffled))
    half=int(np.argmin(abs(phases-.5)))
    # Descriptive stacked bound, using robust per-epoch scatter rather than independence of cadences.
    d50=np.array([q["depth"] for q in half_depths])
    bounds=np.quantile(np.median(RNG.choice(d50,(4000,len(d50))),axis=1),[.025,.975]) if len(d50) else [None,None]
    return dict(best_phase=float(phases[best]),best_score=float(rank[best]),best_n=int(total[best]),
                phase_scan_signflip_p=float(fap),half_n=len(d50),half_median_depth=float(np.median(d50)) if len(d50) else None,
                half_median_95=np.asarray(bounds).tolist(),half_score=float(stats[half]))


def alias_coverage(curves, fit, ncycle):
    q=fit["nominal"]["q"]
    s=fit["nominal"]["sigma"]
    t0,t1=PAIR[0]+q[0],PAIR[1]+q[1]
    period=(t1-t0)/ncycle
    duration,k,b=q[2:]
    # Half the smaller observed depth, expressed as an injected central depth.
    depth_floor=.5*min(x["depth"] for x in fit["nominal"]["events"])
    rows=[]
    for t,y,sector in curves:
        cycles,ep=centers(t,t0,period)
        cadence=float(np.median(np.diff(t)))
        for n,center in zip(cycles,ep):
            if n in (0,ncycle):
                continue
            alpha=n/ncycle
            # Conservative fit errors and a half-cadence timing floor on each measured event.
            dt=1.96*np.hypot(abs(1-alpha)*np.hypot(s[0],.0104),abs(alpha)*np.hypot(s[1],.0035))
            phase_times=center+dt*np.linspace(-1,1,5)
            results=[]
            for point in phase_times:
                actual=scores(t,y,np.array([point]),duration)
                if not actual:
                    continue
                where=abs(t-point)<2.1*duration
                local_time=t[where]-point
                curve=template(local_time,period,0,duration,k,b,cadence)
                peak=-template(np.array([0.]),period,0,duration,k,b,cadence)[0]
                sim=scores(t[where],y[where]+curve*(depth_floor/peak),np.array([point]),duration)
                if sim:
                    z=actual[0]
                    gain=sim[0]["depth"]-z["depth"]
                    results.append(dict(center=float(point),actual_depth=z["depth"],error=z["error"],points=z["points"],
                                        injected_response=float(gain),injected_score=sim[0]["score"]))
            rejects=len(results)==5 and all(r["injected_score"]>=7 and r["injected_response"]>=7*r["error"]
                                      and max(0,r["actual_depth"])+5*r["error"]<r["injected_response"] for r in results)
            rows.append(dict(sector=sector,cycle=int(n),center=float(center),center_halfwidth_95=float(dt),
                             rejected_here=rejects,offset_measurements=results))
    epoch_errors=np.hypot(s[:2],[.0104,.0035])
    covariance=np.array([[epoch_errors[0]**2,-epoch_errors[0]**2/ncycle],
                         [-epoch_errors[0]**2/ncycle,(epoch_errors**2).sum()/ncycle**2]])
    return dict(integer_cycles=ncycle,period=float(period),period_sigma_floor=float(np.sqrt(covariance[1,1])),
                next=dates(t0,period,covariance),
                central_injection_depth=depth_floor,rejected=any(r["rejected_here"] for r in rows),predicted_in_product_span=rows)


def dates(t0,p,cov):
    # UTC date labels are approximate. Barycentric epochs remain authoritative.
    now=Time("2026-11-09T00:00:00",scale="utc").tdb.jd-2457000
    n=int(np.ceil((now-t0)/p))
    return [dict(cycle=j,btjd=float(t0+j*p),date_approx=Time(2457000+t0+j*p,format="jd",scale="tdb").utc.isot[:10],
                 timing_95_minutes=float(1.96*1440*np.sqrt(np.array([1,j])@cov@np.array([1,j]))))
            for j in range(n,n+3)]


if __name__=="__main__":
    f=json.loads((HERE/"transit_fits.json").read_text())
    catalog={x["catalog"]["TICID"]:x for x in json.loads((ROOT/"cohort.json").read_text())}
    out={}
    for tic in TICS:
        curves=[(*read_curve(p)[:2],p["sequence_number"]) for p in products(catalog[tic])]
        m=f[tic]["nominal"]
        q=m["q"]
        if tic in REPEAT:
            t0,p=REPEAT[tic][0]+q[0],REPEAT[tic][1]+q[1]
            cov=np.array(m["covariance"])[:2,:2]
            out[tic]=dict(period=p,epoch=t0,odd_even=parity(m["events"]),secondary=secondary(curves,t0,p,q[2]),
                          next=dates(t0,p,cov))
        else:
            out[tic]=[alias_coverage(curves,f[tic],n) for n in (6,5,3,2)]
    (HERE/"diagnostics.json").write_text(json.dumps(out,indent=2),encoding="utf-8")

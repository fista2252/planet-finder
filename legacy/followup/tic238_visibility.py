"""Project current NASA pointings and compute approximate site visibility."""
import io
import json
import re
import sys
from html import unescape
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
from astropy.utils import iers

HERE=Path(__file__).resolve().parent
DATA=HERE/"tic238_public"
sys.path[:0]=[str(DATA),str(HERE/"visibility_deps")]
from tess_stars2px import tess_stars2px_function_entry

TARGET=SkyCoord(54.9519878477464,-77.5497475569085,unit="deg")
SITES={"CTIO":(-70.8065,-30.1690,2200),"Sutherland":(20.8107,-32.3794,1798),
       "Siding Spring":(149.0661,-31.2733,1165)}
iers.conf.auto_download=False


def pointings():
    source=(DATA/"heasarc_sector.html").read_text(encoding="utf-8")
    table={}
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>",source,re.S):
        cells=re.findall(r"<td[^>]*>(.*?)</td>",row,re.S)
        if len(cells)<6:
            continue
        fields=[unescape(re.sub(r"<[^>]*>","",x)).strip() for x in cells]
        if re.fullmatch(r"s\d{3}",fields[0]) and 98<=int(fields[0][1:])<=121:
            s=int(fields[0][1:])
            ra,dec=map(float,fields[2].split(","))
            table[s]=[s,ra,dec,float(fields[3])]
    # NASA's HTML lacks a starting tr tag on this row.
    if "sector106_summary.html" in source:
        table[106]=[106,25.0148,-40.1221,279.1045]
    return table


if __name__=="__main__":
    table=pointings()
    assert all(i in table for i in (111,112,113))
    text="\n".join(" ".join(map(str,table[i])) for i in sorted(table))
    (DATA/"nasa_pointings.txt").write_text(text+"\n")
    _,_,_,sector,cam,ccd,col,row,_=tess_stars2px_function_entry(
        238187210,TARGET.ra.deg,TARGET.dec.deg,sectorOverrideFile=io.StringIO(text))
    footprint=[dict(sector=int(s),camera=int(c),ccd=int(d),col=float(x),row=float(y))
               for s,c,d,x,y in zip(sector,cam,ccd,col,row) if s>0]
    fit=json.loads((HERE/"timing_check.json").read_text())["238187210"]
    spans=[]
    for event in fit["next"]:
        bary=Time(event["btjd"]+2457000,format="jd",scale="tdb")
        half=2.83+1+event["planning_halfwidth_minutes"]/60
        for name,location in SITES.items():
            loc=EarthLocation.from_geodetic(location[0]*u.deg,location[1]*u.deg,location[2]*u.m)
            mid=Time(bary.jd,format="jd",scale="tdb",location=loc)
            mid=mid-mid.light_travel_time(TARGET,kind="barycentric")
            grid=mid+np.arange(-half,half+.01,.1)*u.h
            frame=AltAz(obstime=grid,location=loc)
            alt=TARGET.transform_to(frame).alt.deg
            sun=get_sun(grid).transform_to(frame).alt.deg
            good=(sun<-12)&(alt>=25)
            spans.append(dict(site=name,utc_mid=mid.utc.isot,required_start=grid[0].utc.isot,
                              required_end=(mid+half*u.h).utc.isot,usable_fraction=float(good.mean()),
                              sun_max_deg=float(sun.max()),target_min_alt_deg=float(alt.min())))
    (HERE/"tic238_visibility.json").write_text(json.dumps(dict(footprint=footprint,sectors_checked=sorted(table),spans=spans),indent=2))

"""Fetch narrow public metadata for the TIC 238187210 review."""
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import archive_check
from inventory import invoke

SAVE=HERE/"tic238_public"
archive_check.DATA=SAVE
RA,DEC=54.9519878477464,-77.5497475569085
G=4626968478435805568


if __name__=="__main__":
    SAVE.mkdir(exist_ok=True)
    gaia="https://gaia.aip.de/tap/sync"
    eso="https://archive.eso.org/tap_obs/sync"
    q1=f"select source_id,radius_gspphot,radius_gspphot_lower,radius_gspphot_upper,distance_gspphot,teff_gspphot,logg_gspphot from gaiadr3.astrophysical_parameters where source_id={G}"
    q2=f"select top 100 dataproduct_type,target_name,instrument_name,s_ra,s_dec,em_res_power,s_resolution,obs_publisher_did from ivoa.ObsCore where 1=contains(point('ICRS',s_ra,s_dec),circle('ICRS',{RA},{DEC},0.008333333))"
    q3=f"select top 100 instrument,target,ra,dec,dp_cat,dp_type,dp_id from dbo.raw where 1=contains(point('ICRS',ra,dec),circle('ICRS',{RA},{DEC},0.008333333)) and dp_cat='SCIENCE'"
    jobs=[("gaia_ap.xml",gaia,dict(REQUEST="doQuery",LANG="ADQL",QUERY=q1)),
          ("eso_reduced.csv",eso,dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",QUERY=q2)),
          ("eso_raw.csv",eso,dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",QUERY=q3)),
          ("mit_year9.html","https://tess.mit.edu/tess-year-9-observations/",None),
          ("heasarc_sector.html","https://heasarc.gsfc.nasa.gov/docs/tess/sector.html",None),
          ("heasarc_ddt.html","https://heasarc.gsfc.nasa.gov/docs/tess/ddt.html",None)]
    with ThreadPoolExecutor(max_workers=5) as ex:
        result=list(ex.map(lambda args:archive_check.fetch(*args),jobs))
    try:
        mast=invoke("Mast.Caom.Cone",dict(ra=RA,dec=DEC,radius=1/60))
        (SAVE/"mast_cone.json").write_text(json.dumps(mast["data"],indent=2))
    except Exception as error:
        result.append(dict(file="mast_cone.json",error=str(error)))
    (SAVE/"manifest.json").write_text(json.dumps(result,indent=2))
    print(result,flush=True)

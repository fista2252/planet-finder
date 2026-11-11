"""Fetch small pixel cutouts around the leading sparse repeat-like cases."""
from concurrent.futures import ThreadPoolExecutor
import json

import requests

from inventory import ROOT

REQUESTS = (("339672790", 7), ("339672790", 33), ("273607159", 27), ("271693283", 64),
            ("389667951", 1), ("389667951", 96), ("269730302", 8), ("269730302", 68),
            ("238187210", 1), ("238187210", 66), ("150108429", 63), ("340562993", 61))
URL = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut"


def fetch(tic, sector):
    nearby = json.loads((ROOT / "data/nearby" / f"{tic}.json").read_text())
    star = next(row for row in nearby if str(row["ID"]) == tic)
    params = dict(ra=star["ra"], dec=star["dec"], x=7, y=7, units="px", sector=sector)
    path = ROOT / "data/cutouts" / f"{tic}-s{sector:04}.zip"
    try:
        if not path.exists():
            r = requests.get(URL, params=params, timeout=(15, 110))
            r.raise_for_status()
            if not r.content.startswith(b"PK") or len(r.content) > 65_000_000:
                raise ValueError("Unexpected cutout response")
            path.write_bytes(r.content)
        return dict(tic=tic, sector=sector, params=params, bytes=path.stat().st_size)
    except (ValueError, OSError, requests.RequestException) as error:
        return dict(tic=tic, sector=sector, params=params, error=str(error))


if __name__ == "__main__":
    (ROOT / "data/cutouts").mkdir(exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        out = list(pool.map(lambda pair: fetch(*pair), REQUESTS))
    (ROOT / "cutout_manifest.json").write_text(json.dumps(out, indent=2))
    print(out)

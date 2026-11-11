"""Record TIC neighbors within one arcminute of each selected target."""
from concurrent.futures import ThreadPoolExecutor
import json

from inventory import ROOT, invoke


def query(target):
    tic = target["catalog"]["TICID"]
    product = target["primary"][0]
    path = ROOT / "data/nearby" / f"{tic}.json"
    if path.exists():
        return tic, "cached"
    try:
        response = invoke("Mast.Catalogs.Tic.Cone",
                          {"ra": product["s_ra"], "dec": product["s_dec"], "radius": 1/60})
        fields = ("ID", "ra", "dec", "Tmag", "dstArcSec", "GAIA", "objType")
        rows = [{name: row.get(name) for name in fields} for row in response["data"]]
        path.write_text(json.dumps(rows, indent=2), encoding="utf-8")
        return tic, len(rows)
    except Exception as error:
        return tic, str(error)


def main():
    (ROOT / "data/nearby").mkdir(exist_ok=True)
    cohort = json.loads((ROOT / "cohort.json").read_text())
    with ThreadPoolExecutor(max_workers=4) as pool:
        for row in pool.map(query, cohort):
            print(*row, flush=True)


if __name__ == "__main__":
    main()

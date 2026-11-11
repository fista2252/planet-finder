"""Retrieve the archived comparison apertures used by the neighbor diagnostic."""
import json

from download_curves import CACHE, ROOT, download

if __name__ == "__main__":
    CACHE.mkdir(exist_ok=True)
    products = json.loads((ROOT/"data/neighbor_220429131_inventory.json").read_text())
    for p in products:
        sector = p["sequence_number"]
        if sector in (1,2,5,7,31) and (p["provenance_name"] == "QLP" or sector == 1 and p["provenance_name"] == "T16"):
            print(download(p["dataURL"]))

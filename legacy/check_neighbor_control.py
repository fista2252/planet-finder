"""Compare TIC 220429133 and its brighter 22-arcsecond neighbor at shared epochs."""
import json

import numpy as np

from screen import ROOT, read_curve, scores

# Fixed epoch and duration from the initial target-only scan.
EVENTS = {1: (1328.4901418532245, .2286), 2: (1362.8241029328656, .28575),
          5: (1462.3035477232686, .17145), 7: (1496.6779542290205, .28575),
          31: (2165.7409531103312, .17145)}


def main():
    targets = json.loads((ROOT / "cohort.json").read_text())
    target = next(row for row in targets if row["catalog"]["TICID"] == "220429133")
    neighbor = json.loads((ROOT / "data/neighbor_220429131_inventory.json").read_text())
    comparisons = []
    for sector, (epoch, duration) in EVENTS.items():
        for tic, products in (("220429133", target["primary"]), ("220429131", neighbor)):
            for p in products:
                if p["sequence_number"] != sector or p["provenance_name"] not in {"QLP", "T16"}:
                    continue
                if sector != 1 and p["provenance_name"] == "T16":
                    continue
                try:
                    t, y, _ = read_curve(p)
                    result = scores(t, y, np.asarray([epoch]), duration)
                except (OSError, KeyError, ValueError) as error:
                    comparisons.append(dict(tic=tic, sector=sector, uri=p["dataURL"], error=str(error)))
                    continue
                comparisons.append(dict(tic=tic, sector=sector, uri=p["dataURL"], measurement=result[0] if result else None))
    (ROOT / "neighbor_220429131_check.json").write_text(json.dumps(comparisons, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

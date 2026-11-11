"""Reproduce the aggregate ephemeris and difference-image diagnostics."""
import csv
import json
from math import dist
from statistics import median

from screen import ROOT


def timing_metrics():
    rows = []
    for result in json.loads((ROOT/"ephemeris_vetting.json").read_text()):
        measurements = result["measurements"]
        primary = [m for m in measurements if m["phase"] == 0 and m["score"] >= 7]
        later = [m for m in primary if m["sector"] >= 55]
        half = [m for m in measurements if m["phase"] == .5]
        depth = median(m["depth"] for m in primary)
        rows.append(dict(tic=result["tic"], **result["fit"], fit_used_late=result["fit_used_late"],
                         primary_events=len(primary), median_depth=depth, later_events=len(later),
                         later_sectors=len({m["sector"] for m in later}),
                         halfway_usable=len(half), halfway_over_5=sum(m["score"] >= 5 for m in half),
                         halfway_quiet_sensitive=sum(max(0,m["depth"])+5*m["error"] < depth/2 for m in half),
                         odd_even=[dict(count=len(d), median_depth=median(d) if d else None)
                                   for parity in (0,1) for d in [[m["depth"] for m in primary if m["cycle"] % 2 == parity]]]))
    return rows


def image_metrics():
    rows = []
    for result in json.loads((ROOT/"pixel_vetting.json").read_text()):
        row = {k:result[k] for k in ("tic","sector","epoch","duration")}
        c = result.get("centroid_xy")
        if c:
            ranked = sorted((dist(c,(s["x"],s["y"])),s["tic"]) for s in result["catalog_pixels"])
            row.update(target_distance_pixels=next(d for d,tic in ranked if tic==result["tic"]),
                       nearest_tic=ranked[0][1], nearest_distance_pixels=ranked[0][0],
                       local_aperture_depth=result["local_aperture_depth"],
                       peak_relative_to_largest_local_control=result["peak_flux"]/max(result["null_max_abs_pixel"])
                       if result["null_max_abs_pixel"] and max(result["null_max_abs_pixel"]) > 0 else None)
        else:
            row["error"] = result.get("error")
        rows.append(row)
    return rows


if __name__ == "__main__":
    (ROOT/"ephemeris_metrics.json").write_text(json.dumps(timing_metrics(),indent=2),encoding="utf-8")
    rows = image_metrics()
    with (ROOT/"pixel_summary.csv").open("w",encoding="utf-8",newline="") as f:
        fields = list(dict.fromkeys(k for row in rows for k in row))
        writer = csv.DictWriter(f,fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

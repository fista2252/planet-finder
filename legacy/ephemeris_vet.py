"""Follow up coherent multi-sector peaks with timing fits and fixed predictions."""
import json

import numpy as np

from aliases import compatible_peaks
from screen import ROOT, read_curve, scores

PERIOD_HINTS = {"389667951": 15.796873, "269730302": 28.430766, "238187210": 26.026897,
                "150108429": 49.559558, "340562993": 20.235459}


def fit_epochs(result, period, late=False):
    first = result["original"]
    epoch = first["tc"]
    candidates = compatible_peaks(result)
    selected = {0: first}
    for peak in candidates:
        if not late and peak["sector"] >= 55:
            continue
        number = round((peak["tc"]-epoch)/period)
        if number not in selected and abs(peak["tc"]-epoch-number*period) <= .085:
            selected[number] = peak
    cycles = np.array(sorted(selected))
    times = np.array([selected[k]["tc"] for k in cycles])
    if len(times) < 3:
        return None
    p, t0 = np.polyfit(cycles, times, 1)
    err = times - (t0 + p * cycles)
    return dict(period=float(p), epoch=float(t0), rms_minutes=float(np.sqrt(np.mean(err**2))*1440),
                events=len(times), last_training_epoch=float(max(times)))


def fixed_measurements(item, first, fit):
    period, epoch = fit["period"], fit["epoch"]
    d = first["duration"]
    target_sector = int(item["catalog"]["Sector"])
    results = []
    for product in item["primary"]:
        sector = product["sequence_number"]
        if sector == target_sector and product["provenance_name"] != "T16":
            continue
        t, y, _ = read_curve(product)
        cycle = np.arange(np.ceil((t[0]-epoch)/period)-1, np.floor((t[-1]-epoch)/period)+1)
        for phase in (0, .5):
            centers = epoch + (cycle + phase)*period
            for m in scores(t, y, centers, d):
                number = round((m["tc"]-epoch)/period - phase)
                results.append(dict(**m, sector=sector, cycle=number, phase=phase, uri=product["dataURL"]))
    return results


def analyze(item, result):
    tic = item["catalog"]["TICID"]
    hint = PERIOD_HINTS[tic]
    fit = fit_epochs(result, hint)
    used_late = fit is None
    fit = fit if fit else fit_epochs(result, hint, late=True)
    if fit is None:
        return dict(tic=tic, error="fewer than three compatible epoch estimates")
    points = fixed_measurements(item, result["original"], fit)
    return dict(tic=tic, fit=fit, fit_used_late=used_late, measurements=points)


if __name__ == "__main__":
    cohort = json.loads((ROOT / "cohort.json").read_text())
    screen = json.loads((ROOT / "screen.json").read_text())
    result = [analyze(i,r) for i,r in zip(cohort,screen,strict=True) if i["catalog"]["TICID"] in PERIOD_HINTS]
    (ROOT / "ephemeris_vetting.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    for row in result:
        print(row["tic"], row.get("fit"), row.get("fit_used_late"), len(row.get("measurements",[])))

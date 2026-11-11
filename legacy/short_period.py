"""Check for sub-two-day repetition using one search sector and separate holdouts."""
import json

import numpy as np
from astropy.timeseries import BoxLeastSquares

from screen import ROOT, read_curve


def folded_depth(t, y, period, t0, duration):
    phase = (t - t0 + period/2) % period - period/2
    inside = abs(phase) < duration/2
    outside = (abs(phase) > duration) & (abs(phase) < .45*period)
    cycles = np.unique(np.rint((t[inside]-t0)/period))
    if inside.sum() < 6 or outside.sum() < 10 or len(cycles) < 3:
        return None
    # Use individual orbit contrasts so a single deep feature cannot drive the check.
    differences = []
    for cycle in cycles:
        same = abs(t-(t0+cycle*period)) < period/2
        a, b = y[same & inside], y[same & outside]
        if len(a) >= 2 and len(b) >= 4:
            differences.append(float(np.median(b)-np.median(a)))
    if len(differences) < 3:
        return None
    return dict(depth=float(np.median(differences)), cycles=len(differences),
                positive_fraction=float(np.mean(np.asarray(differences) > 0)))


def evaluate(item, result):
    later = sorted((s for s in result.get("sectors", []) if s["sector"] > int(item["catalog"]["Sector"])
                    and s.get("peaks")), key=lambda s: -s["peaks"][0]["score"])
    if not later:
        return None
    product = {p["dataURL"]: p for p in item["primary"]}
    # The strongest later sector supplies the trial period. Holdouts are the
    # three nearest other sectors, with times kept fixed when scoring them.
    training = later[0]
    t, y, _ = read_curve(product[training["uri"]])
    cadence = np.median(np.diff(t))
    d = float(item["catalog"]["Duration"])
    durations = np.unique(np.clip(np.asarray([.3, .55, .8])*d, 2*cadence, .20))
    low = max(.25, 3*max(durations))
    if low >= 2:
        return None
    model = BoxLeastSquares(t, y)
    fit = model.autopower(durations, minimum_period=low, maximum_period=2, frequency_factor=2)
    k = int(np.argmax(fit.power))
    p, ep, dur = map(float, (fit.period[k], fit.transit_time[k], fit.duration[k]))
    training_result = folded_depth(t,y,p,ep,dur)
    held = []
    for s in sorted(later[1:], key=lambda s: abs(s["t0"]-training["t0"]))[:3]:
        t, y, _ = read_curve(product[s["uri"]])
        held.append(dict(sector=s["sector"], result=folded_depth(t,y,p,ep,dur)))
    return dict(period=p, epoch=ep, duration=dur, training_sector=training["sector"],
                training_result=training_result, holdout=held)


if __name__ == "__main__":
    cohort = json.loads((ROOT / "cohort.json").read_text())
    screen = json.loads((ROOT / "screen.json").read_text())
    out = [dict(tic=item["catalog"]["TICID"], short_period=evaluate(item, result))
           for item,result in zip(cohort,screen,strict=True)]
    (ROOT / "short_period.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

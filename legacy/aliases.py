"""Bound nominal period aliases using measured cadence coverage and local sensitivity."""
import json

import numpy as np

from screen import ROOT, read_curve, scores


def compatible_peaks(result):
    event = result.get("original")
    if not event or event["depth"] <= 0:
        return []
    sector = int(result["catalog"]["Sector"])
    compatible = []
    for s in result["sectors"]:
        if s["sector"] <= sector:
            continue
        for p in s["peaks"] if "peaks" in s else ():
            if p["score"] >= max(8, s["inverse_p995"] or 0) and .35 <= p["depth"] / event["depth"] <= 2.5:
                compatible.append(dict(**p, sector=s["sector"], uri=s["uri"], inverse_max=s["inverse_max"]))
    return sorted(compatible, key=lambda p: -p["score"])


def predictions(start, end, t0, periods):
    owners, epochs = [], []
    for i, p in enumerate(periods):
        cycle = np.arange(np.ceil((start - t0)/p), np.floor((end - t0)/p)+1, dtype=int)
        cycle = cycle[cycle != 0]
        owners.extend([i] * len(cycle))
        epochs.extend(t0 + cycle * p)
    return np.asarray(owners, dtype=int), np.asarray(epochs, dtype=float)


def exclude(curves, first, periods, second=None):
    excluded = np.zeros(len(periods), dtype=bool)
    coverage = excluded.copy()
    t0, duration = first["tc"], first["duration"]
    # Require evidence against a transit even at half the measured depth.
    threshold = .5 * min(first["depth"], second["depth"] if second else first["depth"])
    base_tol = max(.1 * duration, first["cadence"])
    for t, y in curves:
        owners, epochs = predictions(t[0], t[-1], t0, periods)
        tol = np.full(len(epochs), base_tol)
        if second:
            interval = second["tc"] - t0
            pos = (epochs - t0) / interval
            tol = np.abs(1-pos) * base_tol + np.abs(pos) * max(.1 * duration, second["cadence"])
            other = np.abs(epochs-second["tc"]) > duration
            owners, epochs, tol = owners[other], epochs[other], tol[other]
        # Widely uncertain predicted centers are left unresolved.
        fixed = tol <= duration/3
        owners, epochs, tol = owners[fixed], epochs[fixed], tol[fixed]
        if not len(epochs):
            continue
        requested = np.unique((epochs[:, None] + tol[:, None] * np.array([-1, -.5, 0, .5, 1])).ravel())
        lookup = {m["tc"]: m for m in scores(t, y, requested, duration)}
        for owner, epoch, width in zip(owners, epochs, tol):
            samples = [lookup.get(epoch + width * shift) for shift in (-1, -.5, 0, .5, 1)]
            if any(sample is None for sample in samples):
                continue
            coverage[owner] = True
            if all(max(0, sample["depth"]) + 5 * sample["error"] < threshold for sample in samples):
                excluded[owner] = True
    return dict(tested=len(periods), sufficient_cadence=int(coverage.sum()), excluded=int(excluded.sum()),
                surviving_periods_days=[float(x) for x in periods[~excluded]])


def summarize(item, result):
    original = result.get("original")
    if not original or original["score"] < 7:
        return dict(tic=item["catalog"]["TICID"], status="original_below_screen_threshold")
    curves = []
    original_sector = int(item["catalog"]["Sector"])
    for p in item["primary"]:
        if p["sequence_number"] == original_sector and p["provenance_name"] != "T16":
            continue
        try:
            t, y, _ = read_curve(p)
        except (OSError, KeyError, ValueError):
            continue
        curves.append((t, y))
    grid = np.geomspace(2, 500, 1000)
    single = exclude(curves, original, grid)
    # In the report, single-event exclusions are discrete diagnostics, not intervals or posteriors.
    single["unexcluded"] = len(single["surviving_periods_days"])
    candidates = compatible_peaks(result)
    pairs = []
    for second in candidates[:4]:
        dt = second["tc"] - original["tc"]
        n = np.arange(max(1, int(np.ceil(dt/500))), int(dt/2)+1)
        periods = dt / n
        pairs.append(dict(second=second, **exclude(curves, original, periods, second)))
    return dict(tic=item["catalog"]["TICID"], status="screened", curves=len(curves),
                single_grid=single, compatible_peaks=len(candidates), pairs=pairs)


def main():
    cohort = json.loads((ROOT / "cohort.json").read_text(encoding="utf-8"))
    screen = json.loads((ROOT / "screen.json").read_text(encoding="utf-8"))
    rows = []
    for item, result in zip(cohort, screen, strict=True):
        assert item["catalog"]["TICID"] == result["catalog"]["TICID"]
        row = summarize(item, result)
        rows.append(row)
        print(row["tic"], row["status"], row.get("single_grid"), row.get("compatible_peaks"), flush=True)
    (ROOT / "aliases.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

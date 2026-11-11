"""Bounded periodic box search and per-event cadence-aware diagnostics."""
import math

import numpy as np
from astropy.timeseries import BoxLeastSquares


def search(curves, config, seed=None):
    from screen import scores

    if not curves:
        return {"status": "no_usable_observations"}
    times = np.concatenate([x["time"] for x in curves])
    flux = np.concatenate([x["flux"] - np.median(x["flux"]) for x in curves])
    baseline = float(np.ptp(times))
    if len(times) > 90000:
        return {"status": "sample_limit", "cadences": len(times)}
    if baseline < 2 * config["min_period_days"] and seed is None:
        return {"status": "insufficient_baseline", "baseline_days": baseline}
    settings = None
    if seed is None:
        durations = config["durations_days"]
        max_period = min(config["max_period_days"], baseline / 2)
        fmin, fmax = 1 / max_period, 1 / config["min_period_days"]
        desired = max(2, math.ceil((fmax - fmin) * baseline**2 / min(durations)) + 1)
        frequencies = np.linspace(fmin, fmax, min(config["grid_points"], desired))
        periods = np.unique(1 / frequencies)
        trial = BoxLeastSquares(times, flux).power(periods, durations, objective="snr")
        finite = np.flatnonzero(np.isfinite(trial.power))
        if not len(finite):
            return {"status": "no_finite_search_scores"}
        best = finite[np.argmax(trial.power[finite])]
        period, epoch, duration = map(float, (trial.period[best], trial.transit_time[best], trial.duration[best]))
        depth, score = map(float, (trial.depth[best], trial.depth_snr[best]))
        inverted = BoxLeastSquares(times, -flux).power(periods, durations, objective="snr")
        finite_null = inverted.power[np.isfinite(inverted.power)]
        settings = dict(grid_points=len(periods), ideal_grid_points=desired,
                        search_snr=score, inverted_max_snr=float(max(finite_null)) if len(finite_null) else None,
                        grid_limited=len(periods) < desired)
    else:
        period, epoch, duration = (seed[k] for k in ("period_days", "epoch_btjd", "duration_days"))
        depth = None
    estimated_events = sum(max(0, math.floor((x["time"][-1] - epoch) / period) -
                               math.ceil((x["time"][0] - epoch) / period) + 1) for x in curves)
    if estimated_events > 1000:
        return {"status": "event_limit", "predicted_events": estimated_events}
    events, secondary, controls = [], [], []
    for curve in curves:
        t, y = curve["time"], curve["flux"]
        sector = curve["sector"]
        cycles = np.arange(math.ceil((t[0] - epoch) / period), math.floor((t[-1] - epoch) / period) + 1)
        centers = epoch + period * cycles
        measurements = {row["tc"]: row for row in scores(t, y, centers, duration)}
        half_epoch = epoch + period / 2
        half_cycles = np.arange(math.ceil((t[0] - half_epoch) / period), math.floor((t[-1] - half_epoch) / period) + 1)
        halves = half_epoch + period * half_cycles
        secondary.extend(dict(sector=sector, **r) for r in scores(t, y, halves, duration))
        for cycle, center in zip(cycles, centers):
            m = measurements.get(center)
            events.append(dict(sector=sector, cycle=int(cycle), midpoint_btjd=float(center),
                               coverage="measured" if m else "insufficient_cadence", measurement=m))
        trial_times = np.linspace(t[0] + 2 * duration, t[-1] - 2 * duration, 32)
        phase_dist = np.abs((trial_times - epoch + period / 2) % period - period / 2)
        controls.extend((t, y, x) for x in trial_times[phase_dist > 2 * duration])
    measured = [e for e in events if e["measurement"]]
    support = [e for e in measured if e["measurement"]["score"] >= 5]
    if measured and depth is None:
        depth = float(np.median([e["measurement"]["depth"] for e in measured]))
    null, injected = 0, 0
    testable = 0
    injection_depth = max(.0005, .5 * depth) if depth and depth > 0 else None
    for t, y, center in controls:
        trial = np.array([center])
        inv = scores(t, -y, trial, duration)
        if not inv:
            continue
        testable += 1
        null += inv[0]["score"] >= 5
        if injection_depth is not None:
            model = (abs(t - center) < duration / 2) * injection_depth
            recovered = scores(t, y - model, trial, duration)
            injected += bool(recovered and recovered[0]["score"] >= 5)
    parity = [[e["measurement"] for e in measured if e["cycle"] % 2 == i] for i in range(2)]
    odd_even = None
    if all(len(group) >= 2 for group in parity):
        medians = [float(np.median([r["depth"] for r in group])) for group in parity]
        error = math.hypot(*(float(np.median([r["error"] for r in group])) for group in parity))
        odd_even = dict(even_depth=medians[0], odd_depth=medians[1], nominal_difference_score=(medians[0] - medians[1]) / error)
    max_secondary = max((r["score"] for r in secondary), default=None)
    flags = []
    if len({e["cycle"] for e in support}) < 3:
        flags.append("fewer_than_three_supported_cycles")
    if max_secondary is not None and max_secondary >= 5:
        flags.append("possible_half_phase_secondary")
    if odd_even and abs(odd_even["nominal_difference_score"]) >= 5:
        flags.append("possible_odd_even_difference")
    if settings and settings["inverted_max_snr"] is not None and settings["inverted_max_snr"] >= settings["search_snr"]:
        flags.append("inverted_search_as_strong")
    if not testable:
        flags.append("local_controls_unavailable")
    if settings and settings["grid_limited"]:
        flags.append("period_grid_truncated")
    return dict(status="measured" if measured else "no_measurable_events", mode="seeded" if seed else "period_grid",
                period_days=period, epoch_btjd=epoch, duration_days=duration, depth_fraction=depth,
                search=settings, events=events, supported_cycles=len({e["cycle"] for e in support}),
                odd_even=odd_even, half_phase_measurements=secondary, flags=flags,
                local_controls=dict(testable=testable, inverted_at_or_above_5=int(null), box_depth_fraction=injection_depth,
                                    recovered_at_or_above_5=injected if injection_depth is not None else None),
                uncertainty="No calibrated period or timing uncertainty. Local errors use sideband scatter and are nominal.")

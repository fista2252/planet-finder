"""Bounded periodic box search and per-event cadence-aware diagnostics."""
import math

import numpy as np
from astropy.timeseries import BoxLeastSquares

MIN_SCORE = 5
SHORTLIST = 24
MIN_CYCLES = 3


def event_peaks(curves, durations):
    from screen import scores, separated_peaks

    peaks = []
    for curve in curves:
        t, y = curve["time"], curve["flux"]
        rows = [row for duration in durations for row in scores(t, y, t, duration)
                if row["score"] >= 7]
        peaks.extend(separated_peaks(rows, max(durations), count=3))
    return sorted(peaks, key=lambda row: -row["score"])


def refine_events(peaks, minimum, maximum):
    """Generate aliases from separated dip pairs, then fit matching event times."""
    times = np.array([row["tc"] for row in peaks])
    durations = np.array([row["duration"] for row in peaks])
    if len(times) < MIN_CYCLES:
        return [], False
    proposals = []
    dropped_support = 0
    # Try each of the six strongest dips paired with every other retained dip.
    # This preserves short and long gaps without a quadratic scan of all cadences.
    for i in range(min(6, len(times))):
        for j in range(i + 1, len(times)):
            span = abs(times[j] - times[i])
            cycles = np.arange(max(1, math.ceil(span / maximum)), math.floor(span / minimum) + 1)
            if not len(cycles):
                continue
            periods = span / cycles
            phase = (times - times[i]) / periods[:, None]
            offset = abs(phase - phase.round()) * periods[:, None]
            matches = offset <= .65 * durations
            counts = matches.sum(axis=1)
            valid = np.flatnonzero(counts >= MIN_CYCLES)
            rms = np.sqrt(np.sum(np.square(offset[valid]) * matches[valid], axis=1) / counts[valid])
            order = np.lexsort((rms, -counts[valid]))
            if len(order) > 2 * SHORTLIST:
                dropped_support = max(dropped_support, int(counts[valid[order[2 * SHORTLIST]]]))
            for idx in order[:2 * SHORTLIST]:
                k = valid[idx]
                proposals.append((int(counts[k]), float(rms[idx]), float(periods[k]), float(times[i]), matches[k]))
    proposals.sort(key=lambda item: (-item[0], item[1], -item[2]))
    refined = []
    baseline = float(np.ptp(times))
    for count, _, period, epoch, matches in proposals:
        indices = np.flatnonzero(matches)
        cycle = np.rint((times[indices] - epoch) / period).astype(int)
        if len(set(cycle)) < MIN_CYCLES:
            continue
        # The nearest, strongest anchors define the trial only. Full cadence
        # measurements decide whether these events actually support the trial.
        period, offset = np.polyfit(cycle, times[indices] - epoch, 1)
        duration = float(np.median(durations[indices]))
        if not minimum <= period <= maximum:
            continue
        model = dict(period_days=float(period), epoch_btjd=float(epoch + offset),
                     duration_days=duration, method="event_pairs", anchors=count)
        if any(same_solution(model, item, baseline) for item in refined):
            continue
        if len(refined) >= SHORTLIST:
            return refined, max(count, dropped_support) >= refined[0]["anchors"] - 1
        refined.append(model)
    return refined, bool(refined and dropped_support >= refined[0]["anchors"] - 1)


def same_solution(a, b, baseline):
    p, q = a["period_days"], b["period_days"]
    tolerance = .35 * min(a["duration_days"], b["duration_days"])
    drift = abs(p - q) * max(1, baseline / min(p, q))
    epoch_delta = abs((a["epoch_btjd"] - b["epoch_btjd"] + p / 2) % p - p / 2)
    return drift + epoch_delta <= tolerance


def window_centers(t, epoch, period):
    cycles = np.arange(math.ceil((t[0] - epoch) / period), math.floor((t[-1] - epoch) / period) + 1)
    return cycles, epoch + period * cycles


def cadence_trials(curves, model):
    """Check a trial with box injections at three offsets on actual cadences."""
    from screen import scores

    p, t0, duration = (model[k] for k in ("period_days", "epoch_btjd", "duration_days"))
    windows = [(curve, *window_centers(curve["time"], t0, p)) for curve in curves]
    if sum(len(cycles) for _, cycles, _ in windows) > 1000:
        return dict(events=[], supported_cycles=0, sensitive_missing_cycles=0, event_limit=True)
    events = []
    for curve, cycles, centers in windows:
        offsets = centers[:, None] + duration * np.array([-.25, 0, .25])
        observations = {row["tc"]: row for row in scores(curve["time"], curve["flux"], offsets.ravel(), duration)}
        for cycle, center, trials in zip(cycles, centers, offsets):
            m = observations.get(center)
            events.append(dict(sector=curve["sector"], cycle=int(cycle), midpoint_btjd=float(center),
                               coverage="measured" if m else "insufficient_cadence", measurement=m,
                               timing_samples=[observations.get(offset) for offset in trials]))
    support = [e for e in events if e["measurement"] and e["measurement"]["score"] >= MIN_SCORE]
    # Use a deliberately shallower box than the median observed dip.
    reference_depth = .5 * float(np.median([e["measurement"]["depth"] for e in support])) if support else None
    for event in events:
        measurement = event["measurement"]
        samples = event["timing_samples"]
        # The local statistic is linear: subtracting D from the sampled in-box
        # cadences adds exactly D to the measured depth, with unchanged sidebands.
        # The cadence gate must pass independently at every tested offset.
        injected = [(m["depth"] + reference_depth) / m["error"] if m and reference_depth else None for m in samples]
        sensitive = bool(reference_depth and all(
            m and score >= MIN_SCORE and max(0, m["depth"]) + MIN_SCORE * m["error"] < reference_depth
            for m, score in zip(samples, injected)
        ))
        if measurement and measurement["score"] >= MIN_SCORE:
            state = "supported"
        elif sensitive:
            state = "sensitive_non_detection"
        else:
            state = "insufficient_sensitivity" if measurement else "insufficient_cadence"
        event.update(event_state=state, injected_depth_fraction=reference_depth, injected_scores=injected)
    supported = {event["cycle"] for event in support}
    missing = {event["cycle"] for event in events if event["event_state"] == "sensitive_non_detection"} - supported
    return dict(events=events, supported_cycles=len(supported), sensitive_missing_cycles=len(missing),
                median_support_score=float(np.median([e["measurement"]["score"] for e in support])) if support else 0)


def choose_period(curves, coarse, config, baseline, grid_step):
    peaks = event_peaks(curves, config["durations_days"])
    lo, hi = config["min_period_days"], min(config["max_period_days"], baseline / 2)
    candidates, capped = refine_events(peaks, lo, hi)
    candidates.insert(0, coarse)
    # Explicitly test integer harmonics of the coarse and strongest event trials.
    for trial in candidates[:3]:
        for factor in range(2, 9):
            for period in (trial["period_days"] * factor, trial["period_days"] / factor):
                model = dict(trial, period_days=period, method="integer_harmonic")
                if lo <= period <= hi and not any(same_solution(model, item, baseline) for item in candidates):
                    candidates.append(model)
    evaluated = []
    for candidate in candidates:
        if any(same_solution(candidate, item, baseline) for item in evaluated):
            continue
        result = cadence_trials(curves, candidate)
        # Retain the detail ledger only for the candidate shown to the reader.
        evaluated.append(dict(candidate, **result))
    evaluated.sort(key=lambda r: (r["sensitive_missing_cycles"] >= 2, -r["supported_cycles"],
                                  r["sensitive_missing_cycles"], -r.get("median_support_score", 0)))
    best = evaluated[0]
    viable = [r for r in evaluated if r["supported_cycles"] >= MIN_CYCLES and r["sensitive_missing_cycles"] < 2]
    comparable = [r for r in viable if r["supported_cycles"] >= best["supported_cycles"] - 1]
    alternatives = []
    for row in evaluated:
        misses = [dict(sector=e["sector"], cycle=e["cycle"], midpoint_btjd=e["midpoint_btjd"],
                       depth_fraction=e["measurement"]["depth"], injected_scores=e["injected_scores"])
                  for e in row["events"] if e["event_state"] == "sensitive_non_detection"]
        alternatives.append(dict({k: v for k, v in row.items() if k != "events"},
                                 sensitive_non_detection_examples=misses[:6]))
    coarse_phase_drift = .5 * grid_step * best["period_days"] * baseline
    has_refined_support = any(item["method"] == "event_pairs" and same_solution(best, item, baseline)
                              for item in candidates)
    if not viable or capped or (not has_refined_support and coarse_phase_drift > .25 * best["duration_days"]):
        status = "inconclusive_search"
    elif len(comparable) > 1:
        status = "ambiguous_period"
    else:
        status = "measured"
    return best, status, alternatives, len(peaks), capped


def search(curves, config, seed=None):
    from screen import scores

    if not curves:
        return {"status": "no_usable_observations"}
    times = np.concatenate([x["time"] for x in curves])
    flux = np.concatenate([x["flux"] - np.median(x["flux"]) for x in curves])
    baseline = float(np.ptp(times))
    if len(times) > 90000:
        return {"status": "sample_limit", "cadences": len(times)}
    if baseline > 10000:
        return {"status": "baseline_limit", "baseline_days": baseline}
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
        coarse = dict(period_days=float(trial.period[best]), epoch_btjd=float(trial.transit_time[best]),
                      duration_days=float(trial.duration[best]), method="bls_grid")
        score = float(trial.depth_snr[best])
        inverted = BoxLeastSquares(times, -flux).power(periods, durations, objective="snr")
        finite_null = inverted.power[np.isfinite(inverted.power)]
        grid_step = (fmax - fmin) / (len(frequencies) - 1)
        model, status, alternatives, anchors, capped = choose_period(curves, coarse, config, baseline, grid_step)
        settings = dict(grid_points=len(periods), ideal_grid_points=desired,
                        search_snr=score, inverted_max_snr=float(max(finite_null)) if len(finite_null) else None,
                        grid_limited=len(periods) < desired, event_anchors=anchors, event_shortlist_limited=capped,
                        coarse_period_days=coarse["period_days"], trial_method=model["method"])
    else:
        model = dict(seed, **cadence_trials(curves, seed))
        if model.get("event_limit"):
            return {"status": "event_limit"}
        status = "contradicted_ephemeris" if model["sensitive_missing_cycles"] >= 2 else "measured"
        alternatives = []
    period, epoch, duration = (model[k] for k in ("period_days", "epoch_btjd", "duration_days"))
    events = model["events"]
    secondary, controls = [], []
    for curve in curves:
        t, y = curve["time"], curve["flux"]
        sector = curve["sector"]
        _, halves = window_centers(t, epoch + period / 2, period)
        secondary.extend(dict(sector=sector, **r) for r in scores(t, y, halves, duration))
        trial_times = np.linspace(t[0] + 2 * duration, t[-1] - 2 * duration, 32)
        phase_dist = np.abs((trial_times - epoch + period / 2) % period - period / 2)
        controls.extend((t, y, x) for x in trial_times[phase_dist > 2 * duration])
    measured = [e for e in events if e["measurement"]]
    support = [e for e in measured if e["event_state"] == "supported"]
    depth = float(np.median([e["measurement"]["depth"] for e in support])) if support else None
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
            injection = (abs(t - center) < duration / 2) * injection_depth
            recovered = scores(t, y - injection, trial, duration)
            injected += bool(recovered and recovered[0]["score"] >= 5)
    # Use one measurement per cycle so overlapping coverage cannot inflate
    # the nominal precision of the alternating-depth comparison.
    per_cycle = {}
    for event in measured:
        measurement = event["measurement"]
        current = per_cycle.get(event["cycle"])
        if current is None or measurement["error"] < current["error"]:
            per_cycle[event["cycle"]] = measurement
    parity = [[m for cycle, m in per_cycle.items() if cycle % 2 == i] for i in range(2)]
    odd_even = None
    if all(len(group) >= 2 for group in parity):
        medians = [float(np.median([r["depth"] for r in group])) for group in parity]
        errors = []
        for group, center in zip(parity, medians):
            local = np.median([r["error"] for r in group])
            scatter = 1.4826 * np.median([abs(r["depth"] - center) for r in group])
            errors.append(float(max(local, scatter) / np.sqrt(len(group))))
        error = math.hypot(*errors)
        odd_even = dict(even_depth=medians[0], odd_depth=medians[1], nominal_difference_score=(medians[0] - medians[1]) / error)
    max_secondary = max((r["score"] for r in secondary), default=None)
    flags = []
    if model["supported_cycles"] < MIN_CYCLES:
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
    if model["sensitive_missing_cycles"]:
        flags.append("sensitive_non_detections")
    if not seed and odd_even and abs(odd_even["nominal_difference_score"]) >= 5 and status == "measured":
        status = "ambiguous_period"
    if seed and not measured:
        status = "no_measurable_events"
    chosen = status == "measured" or seed is not None
    return dict(status=status, mode="seeded" if seed else "blind",
                period_days=period if chosen else None, epoch_btjd=epoch if chosen else None,
                trial_period_days=period, trial_epoch_btjd=epoch, duration_days=duration, depth_fraction=depth,
                search=settings, alternatives=alternatives, events=events, supported_cycles=model["supported_cycles"],
                sensitive_missing_cycles=model["sensitive_missing_cycles"],
                odd_even=odd_even, half_phase_measurements=secondary, flags=flags,
                local_controls=dict(testable=testable, inverted_at_or_above_5=int(null), box_depth_fraction=injection_depth,
                                    recovered_at_or_above_5=injected if injection_depth is not None else None),
                uncertainty="No calibrated period or timing uncertainty. Local errors use sideband scatter and are nominal.")

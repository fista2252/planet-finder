"""Small deterministic benchmark for Planet Finder maintenance."""
import math
import sys
import time

import numpy as np

from audit import load
from import_pilot import ROOT, json_bytes, sha
from search import search
from workflow import DEFAULTS

DEFINITION = ROOT / "maintenance_specs/benchmark-v1.json"


def curves_for(case):
    rng = np.random.default_rng(case["rng"])
    for sector, start in enumerate(case["starts"]):
        t = start + np.arange(0, case["days"], 1 / 48)
        y = rng.normal(0, .001, len(t))
        if rho := case.get("rho"):
            for i in range(1, len(y)):
                y[i] = rho * y[i-1] + math.sqrt(1 - rho**2) * y[i]
        if p := case.get("period"):
            phase = (t - 1501 + p / 2) % p - p / 2
            y[abs(phase) < case["duration"] / 2] -= case["depth"] * case.get("flux_fraction", 1)
            half = abs(abs(phase) - p / 2)
            if case.get("secondary_depth"):
                y[half < case["duration"] / 2] -= case["secondary_depth"]
            if case.get("omit_half_phase"):
                keep = half > .4
                t, y = t[keep], y[keep]
        yield dict(sector=sector, time=t, flux=y)


def evaluate():
    sys.path.insert(0, str(ROOT / "legacy"))
    definition = load(DEFINITION)
    outcomes, violations = [], []
    for case in definition["cases"]:
        t0, cpu0 = time.perf_counter(), time.process_time()
        curves = list(curves_for(case))
        cfg = DEFAULTS.copy()
        if bounds := case.get("bounds"):
            cfg.update(min_period_days=bounds[0], max_period_days=bounds[1])
        seed = (dict(period_days=case["seed_period"], epoch_btjd=1501, duration_days=case["duration"])
                if "seed_period" in case else None)
        result = search(curves, cfg, seed)
        accepted = result["status"] == "measured" and result.get("supported_cycles", 0) >= 3
        period = result.get("period_days") if accepted else None
        role = case["role"]
        if role == "injected":
            outcome = ("abstained" if not accepted else
                       "recovered" if abs(period - case["period"]) <= definition["period_tolerance_days"] else "wrong_period")
        elif role == "null":
            outcome = "false_alarm" if accepted else "no_alarm"
        elif role == "gap":
            gaps = sum(event["event_state"] == "insufficient_cadence" for event in result["events"])
            outcome = "gap_preserved" if gaps and not result["sensitive_missing_cycles"] else "gap_misclassified"
        else:
            outcome = "eclipse_held" if not accepted and "possible_odd_even_difference" in result["flags"] else "eclipse_gate_failed"
        if outcome in {"wrong_period", "false_alarm", "gap_misclassified", "eclipse_gate_failed"}:
            violations.append(case["id"])
        outcomes.append(dict(id=case["id"], role=role, status=result["status"], outcome=outcome,
                             accepted_period_days=period if seed is None else None,
                             leading_trial_days=result.get("trial_period_days"),
                             sensitive_missing_cycles=result.get("sensitive_missing_cycles"),
                             cadence_samples=sum(len(x["time"]) for x in curves),
                             wall_seconds=round(time.perf_counter()-t0, 6), cpu_seconds=round(time.process_time()-cpu0, 6)))
    injections = [r for r in outcomes if r["role"] == "injected"]
    nulls = [r for r in outcomes if r["role"] == "null"]
    def rate(rows, outcome):
        return dict(count=sum(row["outcome"] == outcome for row in rows), n=len(rows))
    metrics = {name: rate(rows, name) for rows, name in ((injections, "recovered"), (injections, "wrong_period"),
                                                        (injections, "abstained"), (nulls, "false_alarm"))}
    metrics["ambiguous"] = dict(count=sum(row["status"] == "ambiguous_period" for row in injections), n=len(injections))
    try:
        import resource
    except ImportError:
        peak = None
    else:
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    return dict(id=definition["id"], definition_sha256=sha(DEFINITION.read_bytes()), cases=outcomes,
                violations=violations, metrics=metrics, peak_process_rss_bytes=peak,
                limitations=["Small synthetic fixture set; rates are counts over the listed cases.",
                             "Validation fixtures are public and fixed; there is no blind holdout.",
                             "Dilution scales depth only. This does not simulate source localization."])


if __name__ == "__main__":
    sys.stdout.buffer.write(json_bytes(evaluate()))

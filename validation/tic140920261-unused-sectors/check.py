"""Replay six frozen period/epoch trials using original and unused QLP sectors."""
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "legacy")]

from audit import atomic, load
from import_pilot import json_bytes, sha
from screen import read_curve
from search import MIN_SCORE, cadence_trials
from workflow import SHARED, fetch_product, now, product_cache


def read(product, path):
    time_days, flux, _ = read_curve(product, path=path)
    return dict(sector=product["sequence_number"], time=time_days, flux=flux)


def evaluate():
    from requests.exceptions import RequestException
    from urllib3.exceptions import HTTPError

    protocol_file = HERE / "protocol.json"
    spec = load(protocol_file)
    output = HERE / "result.json"
    if output.exists():
        raise SystemExit("Existing result is frozen. Use a separate checkout for a replay.")
    old, new, receipts, failures = [], [], [], []
    start = time.monotonic()
    budget = [spec["fits_budget_bytes"]]
    for product in spec["old_products"]:
        file, _ = product_cache(product, SHARED)
        if sha(file.read_bytes()) != product["sha256"]:
            raise ValueError("Original FITS digest changed")
        old.append(read(product, file))
    for product in spec["new_products"]:
        try:
            file, receipt = fetch_product(dict(tic=spec["tic"], **product), SHARED, budget, None)
            curve = read(product, file)
            new.append(curve)
            receipts.append(dict(**product, **receipt, cadence_count=len(curve["time"]),
                                 btjd_start=float(curve["time"][0]), btjd_end=float(curve["time"][-1])))
        except (OSError, ValueError, KeyError, HTTPError, RequestException) as error:
            failures.append(dict(uri=product["dataURL"], error=str(error)))
    results = []
    for model in spec["models"]:
        trial = dict(model, duration_days=spec["duration_days"])
        baseline = cadence_trials(old, trial)["events"]
        depth = next((row["injected_depth_fraction"] for row in baseline if row["injected_depth_fraction"]), None)
        events = cadence_trials(new, trial)["events"]
        for event in events:
            samples = event["timing_samples"]
            injected = [(m["depth"] + depth) / m["error"] if m and depth else None for m in samples]
            sensitive = bool(depth and all(
                m and score >= MIN_SCORE and max(0, m["depth"]) + MIN_SCORE * m["error"] < depth
                for m, score in zip(samples, injected)))
            measurement = event["measurement"]
            if measurement and measurement["score"] >= MIN_SCORE:
                state = "supported"
            elif sensitive:
                state = "sensitive_non_detection"
            else:
                state = "insufficient_sensitivity" if measurement else "insufficient_cadence"
            event.update(event_state=state, injected_scores=injected, injected_depth_fraction=depth)
        supported = {e["cycle"] for e in events if e["event_state"] == "supported"}
        missed = {e["cycle"] for e in events if e["event_state"] == "sensitive_non_detection"} - supported
        results.append(dict(**trial, depth_fraction=depth, supported_cycles=len(supported),
                            sensitive_missing_cycles=len(missed), conditionally_contradicted=len(missed) >= 2,
                            events=events))
    atomic(output, json_bytes(dict(tic=spec["tic"], completed_at=now(), protocol_sha256=sha(protocol_file.read_bytes()),
                                   results=results, products=receipts, failures=failures,
                                   fits_bytes=spec["fits_budget_bytes"] - budget[0], seconds=time.monotonic() - start)))


if __name__ == "__main__":
    evaluate()

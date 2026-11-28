"""Experimental transit morphology and empirical noise diagnostics for frozen windows."""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np


VERSION = 1
MAX_EVENTS = 96
TEMPLATES = ("box", "trapezoid")


def measure_window(t, flux, center, duration, exposure):
    """Fit equal-parameter box and trapezoid models on the same retained cadences."""
    left, right = np.searchsorted(t, [center - 2 * duration, center + 2 * duration])
    x = (t[left:right] - center) / duration
    y = flux[left:right]
    before = x < -.75
    after = x > .75
    in_box = abs(x) <= .5
    if min(before.sum(), after.sum(), in_box.sum()) < 3:
        return None
    # Require coverage on both sides of the midpoint and exclude partially covered gaps.
    cadence = float(np.median(np.diff(t)))
    if (not np.any(in_box & (x < 0)) or not np.any(in_box & (x >= 0))
            or max(np.diff(x[in_box])) > max(.45, 2.5 * cadence / duration)):
        return None

    # Integrate both templates over the exposure. Ingress and egress each occupy one
    # quarter of the trapezoid's duration. Center and duration are fixed by the report.
    z = np.abs(x[:, None] + np.linspace(-.5, .5, 15) * exposure / duration)
    templates = (np.mean(z <= .5, axis=1), np.mean(np.clip((.5 - z) / .25, 0, 1), axis=1))
    baseline = np.column_stack((np.ones(len(x)), x))
    side = before | after
    fit = np.linalg.lstsq(baseline[side], y[side], rcond=None)[0]
    scatter = np.std(y[side] - baseline[side] @ fit, ddof=2)
    if not math.isfinite(scatter) or scatter <= 1e-12:
        return None
    orthogonal, _ = np.linalg.qr(baseline)
    residual_flux = y - orthogonal @ (orthogonal.T @ y)
    result = {}
    for name, template in zip(TEMPLATES, templates):
        residual_template = template - orthogonal @ (orthogonal.T @ template)
        norm = float(residual_template @ residual_template)
        if norm <= 1e-12:
            return None
        depth = -float(residual_template @ residual_flux) / norm
        error = float(scatter / math.sqrt(norm))
        residual = residual_flux + depth * residual_template
        result[name] = dict(depth_fraction=depth, iid_error=error, nominal_score=depth / error,
                            residual_ss=float(residual @ residual))
    return dict(in_box=int(in_box.sum()), fits=result)


def compare_windows(t, flux, centers, duration, exposure):
    """Compare event windows with nonoverlapping control windows from the same product."""
    centers = np.asarray(centers, dtype=float)
    start, stop = t[0] + 2 * duration, t[-1] - 2 * duration
    n_steps = max(0, math.floor((stop - start) / (4 * duration)))
    # Sample at most 256 candidate controls across the product. Controls and event
    # neighborhoods never share sideband or in-event cadences.
    indices = np.unique(np.linspace(0, n_steps, min(n_steps + 1, 256)).astype(int)) if start <= stop else []
    grid = start + np.asarray(indices) * (4 * duration)
    allowed = [float(c) for c in grid if np.all(abs(c - centers) > 4 * duration)]
    control_centers = np.linspace(0, len(allowed) - 1, min(32, len(allowed))).astype(int)
    controls = []
    for index in control_centers:
        c = allowed[index]
        if result := measure_window(t, flux, c, duration, exposure):
            controls.append(dict(btjd=c, **result))

    reference = {}
    for template in TEMPLATES:
        scores = np.asarray([r["fits"][template]["nominal_score"] for r in controls])
        if len(scores) >= 8:
            mad = float(1.4826 * np.median(abs(scores - np.median(scores))))
            reference[template] = (scores, max(1.0, mad))

    events = []
    for center in centers:
        event = dict(btjd=float(center), measurement=measure_window(t, flux, center, duration, exposure))
        if event["measurement"]:
            for template, fit in event["measurement"]["fits"].items():
                if template in reference:
                    scores, scale = reference[template]
                    fit.update(control_noise_scale=scale, adjusted_score=fit["nominal_score"] / scale,
                               controls_at_or_above_nominal_score=int(np.sum(scores >= fit["nominal_score"])))
                else:
                    fit.update(control_noise_scale=None, adjusted_score=None, controls_at_or_above_nominal_score=None)
        events.append(event)
    return dict(events=events, controls=controls)


def diagnose(report, tic, cache):
    from workflow import ROOT, product_cache

    sys.path.insert(0, str(ROOT / "legacy"))
    from screen import read_curve

    target = next((t for t in report["targets"] if t["tic"] == tic), None)
    if target is None:
        raise ValueError(f"TIC {tic} is not in the report")
    events = target["evidence"].get("events", [])
    if not events or len(events) > MAX_EVENTS:
        raise ValueError(f"Expected 1–{MAX_EVENTS} frozen event windows in the report")
    by_sector = {}
    for event in events:
        m = event.get("measurement")
        duration = (m or {}).get("duration", target["evidence"].get("duration_days"))
        center = event["midpoint_btjd"]
        if not isinstance(duration, (int, float)) or not math.isfinite(duration) or not .01 <= duration <= 2:
            raise ValueError("The report must provide a finite duration of 0.01–2 days")
        if not math.isfinite(center):
            raise ValueError("Expected a finite BTJD midpoint")
        by_sector.setdefault(int(event["sector"]), {}).setdefault(duration, set()).add(center)

    results = []
    for product in target["provenance"]:
        sector = int(product["sequence_number"])
        if sector not in by_sector:
            continue
        path, _ = product_cache({"dataURL": product["uri"]}, cache)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != product["sha256"]:
            raise ValueError(f"FITS digest differs from report: {path.name}")
        t, flux, _ = read_curve(product, path=path)
        if len(t) > 90_000:
            raise ValueError("Experimental sample limit exceeded")
        exposure_seconds = product.get("t_exptime")
        measured = isinstance(exposure_seconds, (int, float)) and math.isfinite(exposure_seconds) and exposure_seconds > 0
        exposure = exposure_seconds / 86400 if measured else float(np.median(np.diff(t)))
        for duration, centers in sorted(by_sector[sector].items()):
            all_centers = sorted({c for group in by_sector[sector].values() for c in group})
            result = compare_windows(t, flux, all_centers, duration, exposure)
            result["events"] = [r for r in result["events"] if r["btjd"] in centers]
            results.append(dict(sector=sector, uri=product["uri"], sha256=digest,
                                duration_days=duration, exposure_days=exposure,
                                exposure_source="archive" if measured else "median cadence", **result))
    if not results:
        raise ValueError("No frozen event windows map to the reported FITS products")
    return dict(schema_version=VERSION, experimental=True, tic=tic, baseline_status=target["evidence"]["status"],
                notes=["Trapezoids describe morphology. These fits do not determine the companion's nature.",
                       "Scores are diagnostic and do not account for all noise correlations or the search trials.",
                       "Control windows can contain uncatalogued variability and are not a calibrated false-alarm test."],
                results=results)


def main():
    from workflow import SHARED

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tic", required=True)
    parser.add_argument("--cache", type=Path, default=SHARED, help="Shared cache containing fits/ subdirectory")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    raw = args.report.read_bytes()
    result = diagnose(json.loads(raw), args.tic, args.cache)
    result.update(report_sha256=hashlib.sha256(raw).hexdigest(),
                  technique_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  numpy_version=np.__version__)
    # Reject non-finite output before opening the destination.
    body = json.dumps(result, indent=2, allow_nan=False) + "\n"
    with args.out.open("x", encoding="utf-8") as out:
        out.write(body)
    print(args.out)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as error:
        raise SystemExit(str(error))

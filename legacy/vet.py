"""Vetting diagnostics for matched events, including same-time controls."""
import json

import numpy as np
from astropy.io import fits

from aliases import compatible_peaks
from screen import ROOT, read_curve, scores


def measure_at(t, y, event):
    result = scores(t, y, np.array([event["tc"]]), event["duration"])
    return result[0] if result else None


def product_diagnostics(product, event):
    path = ROOT / "data/curves" / product["dataURL"].rsplit("/", 1)[-1]
    raw = measure_at(*read_curve(product, raw=True)[:2], event)
    with fits.open(path) as f:
        table, header = f[1].data, f[1].header
        t = table["TMID_BJD"] if "TMID_BJD" in table.names else table["TIME"] + float(header.get("BJDREFF", 0))
        d = np.abs(t - event["tc"]) / event["duration"]
        inside, side = d < .5, (d >= .75) & (d < 2)
        diagnostics = {}
        for name in ("SAP_BKG", "SAP_X", "SAP_Y", "MOM_CENTR1", "MOM_CENTR2", "XPOS", "YPOS", "BGV"):
            if name not in table.names:
                continue
            values = table[name].astype(float)
            a, b = values[inside & np.isfinite(values)], values[side & np.isfinite(values)]
            if len(a) < 3 or len(b) < 5:
                continue
            scatter = 1.4826 * np.median(np.abs(b-np.median(b)))
            diagnostics[name] = dict(delta=float(np.median(a)-np.median(b)),
                                     robust_sigma=float(scatter), n_in=len(a))
        frac_bad = None
        if "QUALITY" in table.names and np.any(inside):
            frac_bad = float(np.mean(table["QUALITY"][inside] != 0))
        return dict(raw=raw, diagnostics=diagnostics, fraction_flagged_in_window=frac_bad,
                    crowding=header.get("CROWDSAP"))


def phase_support(p, first, second, later):
    t0 = first["tc"]
    matched = []
    for peak in later:
        if abs(peak["tc"]-second["tc"]) <= first["duration"]:
            continue
        phase_error = abs((peak["tc"]-t0) - round((peak["tc"]-t0)/p)*p)
        if phase_error <= max(.15*first["duration"], first["cadence"], peak["cadence"]):
            matched.append(peak)
    return matched


def run_target(target, result, bounds, screen):
    tic = target["catalog"]["TICID"]
    first = result.get("original")
    peaks = compatible_peaks(result)
    rows = []
    for pair in bounds.get("pairs", []):
        second = pair["second"]
        product = next(p for p in target["primary"] if p["dataURL"] == second["uri"])
        try:
            diagnostics = product_diagnostics(product, second)
        except (OSError, ValueError, KeyError) as exc:
            diagnostics = dict(error=str(exc))
        coincident = []
        for other in screen:
            if other["catalog"]["TICID"] == tic:
                continue
            for segment in other.get("sectors", []):
                if segment["sector"] != second["sector"]:
                    continue
                for peak in segment.get("peaks", []):
                    if peak["score"] >= 6 and abs(peak["tc"]-second["tc"]) < max(.035, second["duration"]/3):
                        coincident.append(dict(tic=other["catalog"]["TICID"], **peak))
        support = sorted(((p, phase_support(p, first, second, peaks)) for p in pair["surviving_periods_days"]),
                         key=lambda v: (-len({k["sector"] for k in v[1]}), -len(v[1]), v[0]))
        rows.append(dict(second=second, diagnostics=diagnostics, coincident=coincident,
                         remaining_aliases=len(support),
                         supported_aliases=[dict(period=p, additional=others) for p, others in support[:6] if others]))
    return dict(tic=tic, vetting=rows)


def main():
    cohort = json.loads((ROOT / "cohort.json").read_text())
    screen = json.loads((ROOT / "screen.json").read_text())
    alias = json.loads((ROOT / "aliases.json").read_text())
    output = []
    for target, result, bounds in zip(cohort, screen, alias, strict=True):
        output.append(run_target(target, result, bounds, screen))
    (ROOT / "vetting.json").write_text(json.dumps(output, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

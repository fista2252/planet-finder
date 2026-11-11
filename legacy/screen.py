"""Measure local negative and inverted box events in the fixed T16 cohort."""
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent


def read_curve(product, raw=False, *, path=None):
    if path is None:
        path = ROOT / "data/curves" / product["dataURL"].rsplit("/", 1)[-1]
    with fits.open(path) as hdus:
        table, header = hdus[1].data, hdus[1].header
        if product["provenance_name"] == "T16":
            t = table["TMID_BJD"].astype(float)
            mag = table["IRM2" if raw else "TFA2"].astype(float)
            mag[(mag < -5) | (mag > 35)] = np.nan
            y = 10 ** (-0.4 * (mag - np.nanmedian(mag)))
            q = np.zeros(len(t), dtype=int)
        else:
            if header.get("TIMESYS", "").upper() != "TDB" or header.get("BJDREFI") != 2457000:
                raise ValueError(f"Unexpected time convention in {path.name}")
            t = table["TIME"].astype(float) + float(header.get("BJDREFF", 0))
            if raw:
                name = "SAP_FLUX"
            elif product["provenance_name"] == "QLP":
                name = "DET_FLUX" if "DET_FLUX" in table.names else "KSPSAP_FLUX"
            else:
                name = "PDCSAP_FLUX"
            y = table[name].astype(float)
            q = table["QUALITY"]
        keep = np.isfinite(t) & np.isfinite(y) & (q == 0)
        t, index = np.unique(t[keep], return_index=True)
        y = y[keep][index]
        if len(t) < 10 or np.median(y) <= 0:
            raise ValueError(f"Insufficient valid flux in {path.name}")
        return t, y / np.median(y) - 1, header.get("CROWDSAP")


def scores(t, y, centers, duration):
    """Subtract a local linear baseline interpolated from two adjacent sidebands."""
    cadence = float(np.median(np.diff(t)))
    if duration < 2 * cadence:
        return []
    bounds = centers[:, None] + duration * np.array([-2, -.75, -.5, .5, .75, 2])
    ix = np.searchsorted(t, bounds)
    left, right = ix[:, [0, 2, 4]], ix[:, [1, 3, 5]]
    n = right - left
    valid = np.all(n >= np.maximum(3, np.floor(duration / cadence * np.array([.5, .6, .5]))), axis=1)
    if not valid.any():
        return []
    n, left, right = n[valid], left[valid], right[valid]
    def mean(v):
        cumulative = np.r_[0, np.cumsum(v)]
        return (cumulative[right] - cumulative[left]) / n
    mean_y, mean_t = mean(y), mean(t - t[0])
    weight = (mean_t[:, 1] - mean_t[:, 0]) / (mean_t[:, 2] - mean_t[:, 0])
    depth = (1 - weight) * mean_y[:, 0] + weight * mean_y[:, 2] - mean_y[:, 1]
    # Conservative local scatter proxy. Significance is nominal under correlated noise.
    variance = np.maximum(0, mean(y * y) - mean_y ** 2)
    sigma = np.maximum(np.sqrt(np.maximum(variance[:, 0], variance[:, 2])), 1e-6)
    error = sigma * np.sqrt(1 / n[:, 1] + (1 - weight) ** 2 / n[:, 0] + weight ** 2 / n[:, 2])
    return [dict(tc=float(c), duration=float(duration), depth=float(d), error=float(e), score=float(d/e),
                 points=int(k), cadence=cadence)
            for c, d, e, k in zip(centers[valid], depth, error, n[:, 1])]


def scan(t, y, duration, tc=None, radius=None):
    limit = max(.20, duration) if radius is None else radius
    trial = t if tc is None else np.unique(np.r_[tc, t[np.abs(t-tc) <= limit]])
    result = []
    for scale in (.75, 1.0, 1.25):
        d = duration * scale
        stride = max(1, int(d / (20 * np.median(np.diff(t)))))
        result.extend(scores(t, y, trial[::stride] if tc is None else trial, d))
    return result


def separated_peaks(results, duration, count=6):
    kept = []
    for row in sorted(results, key=lambda r: -r["score"]):
        if row["score"] <= 0:
            break
        if all(abs(r["tc"] - row["tc"]) > duration for r in kept):
            kept.append(row)
        if len(kept) >= count:
            break
    return kept


def run_target(item):
    c = item["catalog"]
    tc, duration, sector = float(c["Tc"]), float(c["Duration"]), int(c["Sector"])
    orig = next(p for p in item["primary"] if p["provenance_name"] == "T16" and p["sequence_number"] == sector)
    t, y, _ = read_curve(orig)
    original = max(scan(t, y, duration, tc), key=lambda r: r["score"], default=None)
    output = dict(catalog=c, original=original, sectors=[])
    if original:
        t, y, _ = read_curve(orig, raw=True)
        raw = scores(t, y, np.array([original["tc"]]), original["duration"])
        output["original_raw"] = raw[0] if raw else None
    for p in item["primary"]:
        if p is orig:
            continue
        desc = dict(sector=p["sequence_number"], pipeline=p["provenance_name"], uri=p["dataURL"])
        try:
            t, y, crowding = read_curve(p)
            d = original["duration"] if original else duration
            if p["sequence_number"] == sector:
                trial = scan(t, y, d, original["tc"] if original else tc, radius=d/2)
            else:
                trial = scan(t, y, d)
            desc.update(n=len(t), t0=float(t[0]), t1=float(t[-1]), crowding=crowding, trials=len(trial),
                        peaks=separated_peaks(trial, d), inverse_max=-min((v["score"] for v in trial), default=0),
                        inverse_p995=float(np.quantile([-v["score"] for v in trial], .995)) if trial else None)
        except (OSError, KeyError, ValueError) as error:
            desc["error"] = str(error)
        output["sectors"].append(desc)
    return output


def main():
    cohort = json.loads((ROOT / "cohort.json").read_text(encoding="utf-8"))
    output = []
    for item in cohort:
        try:
            row = run_target(item)
        except (OSError, KeyError, ValueError) as error:
            row = dict(catalog=item["catalog"], error=str(error))
        output.append(row)
        print(row["catalog"]["TICID"], row.get("original"), len(row.get("sectors", [])), flush=True)
    (ROOT / "screen.json").write_text(json.dumps(output, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

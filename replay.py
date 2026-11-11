"""Replay screening and vetting with a bounded, verified FITS cache."""
import argparse
import gzip
import json
import platform
import shutil
import sys
from importlib.metadata import version
from pathlib import Path

from audit import CACHE, DATA, atomic, checkpoint, snapshot
from import_pilot import ROOT, json_bytes, sha

LEGACY = ROOT / "legacy"
MAX_BYTES = 300_000_000


def verify_curve(path, expected):
    if not path.is_file():
        return "missing"
    if path.stat().st_size != expected["bytes"] or sha(path.read_bytes()) != expected["sha256"]:
        return "hash_mismatch"
    return "verified"


def resolve_products(item, manifest, source, work, fetch):
    from download_curves import download

    records = []
    for product in item["primary"]:
        uri = product["dataURL"]
        filename = uri.rsplit("/", 1)[-1]
        expected = manifest[uri]
        path = source / filename
        if expected is None:
            state = "no_frozen_hash"
        else:
            state = verify_curve(path, expected)
            if state == "missing" and fetch:
                attempt = download(uri)
                state = verify_curve(path, expected)
                if state == "missing":
                    records.append(dict(uri=uri, status="download_failed", message=attempt["status"]))
                    continue
            if state == "verified":
                destination = work / filename
                if path.resolve() != destination.resolve() and verify_curve(destination, expected) != "verified":
                    shutil.copyfile(path, destination)
        records.append(dict(uri=uri, status=state))
    return records


def local_controls(product, ephemeris, original_duration, limit=64):
    """Run fixed-duration, offset-phase controls on the actual retained cadence."""
    import numpy as np
    from screen import read_curve, scores

    t, y, _ = read_curve(product)
    period, epoch, depth = ephemeris
    centers = np.linspace(t[0] + 2 * original_duration, t[-1] - 2 * original_duration, limit)
    phase = abs((centers - epoch + period / 2) % period - period / 2)
    centers = centers[phase > 4 * original_duration]
    null_scores, recovered = [], 0
    for center in centers:
        trial = np.array([center])
        null = scores(t, -y, trial, original_duration)
        if not null:
            continue
        null_scores.append(null[0]["score"])
        # Fixed-depth box injection on a fixed grid with original noise retained.
        injection = (abs(t - center) < original_duration / 2) * depth
        measured = scores(t, y - injection, trial, original_duration)
        recovered += bool(measured and measured[0]["score"] >= 8)
    return dict(sector=product["sequence_number"], pipeline=product["provenance_name"],
                duration_days=original_duration, injected_depth_fraction=depth,
                testable=len(null_scores), recovered=recovered,
                inverted_at_or_above_8=sum(value >= 8 for value in null_scores),
                inverted_scores=null_scores)


def replay(args):
    sys.path.insert(0, str(LEGACY))
    import aliases
    import download_curves
    import screen
    import vet

    raw = snapshot()
    cohort = json.loads(raw["cohort.json"])
    frozen_screen = json.loads(raw["screen.json"])
    detailed = json.loads(raw["followup/diagnostics.json"])
    manifest = json.loads(gzip.decompress((DATA / "fits-manifest.json.gz").read_bytes()))
    selected = [item for item in cohort if item["catalog"]["TICID"] in args.tic]
    if len(selected) != len(args.tic):
        raise ValueError("Unknown or duplicate TIC")
    projected = sum((manifest[p["dataURL"]] or {}).get("bytes", 0) for item in selected for p in item["primary"])
    if projected > MAX_BYTES:
        raise ValueError("Exceeds the 300 MB source budget")
    run_root = CACHE / "replay"
    source = args.curves.resolve()
    source.mkdir(parents=True, exist_ok=True)
    download_curves.CACHE = source
    run_root.mkdir(parents=True, exist_ok=True)
    if args.fetch:
        # The legacy downloader can accept up to 32 MiB per response. Bound
        # the worst-case transfer, independent of saved expected file sizes.
        unfilled = sum(manifest[p["dataURL"]] is not None and not (source / p["dataURL"].rsplit("/", 1)[-1]).is_file()
                       for item in selected for p in item["primary"])
        if unfilled > 48:
            raise ValueError("At most 48 uncached responses per run. Select fewer TICs or provide a populated cache.")
    work = run_root / "data/curves"
    work.mkdir(parents=True, exist_ok=True)
    screen.ROOT = run_root
    vet.ROOT = run_root
    code = [sha(p.read_bytes()) for p in (Path(__file__), ROOT / "audit.py", LEGACY / "screen.py",
                                          LEGACY / "aliases.py", LEGACY / "vet.py")]
    environment = [platform.python_version(), version("numpy"), version("astropy")]
    for item in selected:
        tic = item["catalog"]["TICID"]
        records = resolve_products(item, manifest, source, work, args.fetch)
        target_dir = run_root / tic
        atomic(target_dir / "ingestion.json", json_bytes(records))
        if any(r["status"] != "verified" for r in records):
            print(f"TIC {tic}: incomplete input, see {target_dir / 'ingestion.json'}")
            continue
        key = sha(json_bytes([item, {r["uri"]: manifest[r["uri"]] for r in records}, code, environment]))
        scan = checkpoint(target_dir / "screen.json", key, lambda: screen.run_target(item))
        bounds = checkpoint(target_dir / "aliases.json", key, lambda: aliases.summarize(item, scan))
        # Retain the original coincidence comparison against the other targets.
        others = [r for r in frozen_screen if r["catalog"]["TICID"] != tic]
        vet_key = sha(json_bytes([key, sha(raw["screen.json"])]))

        checkpoint(target_dir / "vet.json", vet_key, lambda: vet.run_target(item, scan, bounds, [scan, *others]))
        if args.controls and tic in {"238187210", "269730302"}:
            diag = detailed[tic]
            original = scan["original"]
            if original:
                ephemeris = (diag["period"], diag["epoch"], original["depth"] / 2)
                checkpoint(target_dir / "controls.json", key, lambda: [
                    local_controls(product, ephemeris, original["duration"])
                    for product in item["primary"] if product["sequence_number"] != item["catalog"]["Sector"]
                ])
        print(f"TIC {tic}: replay stages in {target_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tic", action="append", required=True)
    parser.add_argument("--curves", type=Path, default=CACHE / "downloads")
    parser.add_argument("--fetch", action="store_true", help="Fetch missing public FITS within the transfer limit")
    parser.add_argument("--controls", action="store_true", help="Run descriptive local injection and inverted controls")
    args = parser.parse_args()
    if not 1 <= len(args.tic) <= 3:
        parser.error("Select one to three TICs per run")
    replay(args)

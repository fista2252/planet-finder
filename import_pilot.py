"""Freeze selected public scientific artifacts from the original pilot."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import shutil
from pathlib import Path
from textwrap import dedent, indent

ROOT = Path(__file__).resolve().parent
COMMIT = "e2eee928fd621f75c1621f6b4206d34bb222ae59"
COMMUNITY = f"https://github.com/eumemic/exoplanets/blob/{COMMIT}/candidates/ephemerides.csv"
T16 = "https://doi.org/10.7910/DVN/CWEUGW"
SOURCES = (
    ("toi_archive.csv", "tid", "pl_orbper", "toi", "toi", "https://exoplanetarchive.ipac.caltech.edu/TAP/"),
    ("confirmed_planets.csv", "tic_id", "pl_orbper", "pl_name", "confirmed", "https://exoplanetarchive.ipac.caltech.edu/TAP/"),
    ("exofop_ctoi.csv", "TIC ID", "Period (days)", "CTOI", "ctoi", "https://exofop.ipac.caltech.edu/tess/"),
    ("tess_ebs.csv", "tess_id", "period", "signal_id", "eclipsing_binary", "https://tessebs.villanova.edu/"),
    ("t16_periodic.csv", "TICID", "P_med", "TICID", "t16_periodic", T16),
    ("singleTransit_candidate_table.csv", "TICID", None, "TICID", "single_event", T16),
)
ARTIFACTS = (
    "cohort.json", "catalog_sources.json", "selection_audit.csv", "download_manifest.json",
    "screen.json", "aliases.json", "vetting.json", "short_period.json", "raven_matches.json",
    "ephemeris_metrics.json", "ephemeris_vetting.json", "neighbor_220429131_check.json",
    "pixel_vetting.json", "pixel_summary.csv", "cutout_manifest.json", "pilot_summary.csv",
    "followup/transit_fits.json", "followup/diagnostics.json", "followup/timing_check.json",
    "followup/centroid_check.json", "followup/unused_observations.json",
    "followup/tic238_audit.json", "followup/tic238_pixel_inputs.json",
    "followup/tic238_pixels.json", "followup/tic238_visibility.json",
    "followup/tic238_public/gaia_ap.xml", "followup/tic238_public/nasa_pointings.txt",
    "followup/tic238_public/eso_cone_retry.csv", "followup/tic238_public/manifest.json",
)
CODE = (
    "screen.py", "aliases.py", "vet.py", "download_curves.py", "cohort.py", "inventory.py",
    "cutouts.py", "neighbors.py", "ephemeris_vet.py", "short_period.py", "report_metrics.py",
    "check_neighbor_control.py", "fetch_neighbor_control.py", "summarize.py", "pixel_vet.py",
    "followup/fit_transits.py", "followup/diagnostics.py", "followup/timing_check.py",
    "followup/centroid_check.py", "followup/unused_observations.py",
    "followup/tic238_audit.py", "followup/tic238_pixels.py", "followup/tic238_plot.py",
    "followup/tic238_visibility.py", "followup/tic238_public.py",
)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def json_bytes(value):
    return (json.dumps(value, ensure_ascii=True, allow_nan=False, sort_keys=True, indent=2) + "\n").encode()


def read_catalog(path, source, ids):
    _, id_col, period_col, label_col, kind, url = source
    raw = path.read_bytes()
    rows = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    if not {id_col, label_col, *([period_col] if period_col else [])} <= set(rows.fieldnames or []):
        raise ValueError(f"Unexpected catalog schema: {path.name}")
    matches = []
    count = 0
    for row in rows:
        count += 1
        tic = row[id_col].strip().removeprefix("TIC ").strip()
        if tic in ids:
            matches.append(dict(tic=tic, period_days=float(row[period_col]) if period_col and row[period_col] else None,
                                label=row[label_col], **({"record": row} if kind == "community" else {})))
    if not count:
        raise ValueError(f"Catalog is empty: {path.name}")
    return dict(kind=kind, url=url, sha256=sha(raw), rows_checked=count, state="checked_snapshot", matches=matches)


def adapt_vetter(source):
    """Expose the original per-target vetter for checkpointed replay."""
    start = source.index('        tic = target["catalog"]')
    end = source.index("        output.append", start)
    block = indent(dedent(source[start:end]), "    ")
    function = "def run_target(target, result, bounds, screen):\n" + block + "    return dict(tic=tic, vetting=rows)\n\n\n"
    output_start = source.index('    (ROOT / "vetting.json")', end)
    source = source[:start] + "        output.append(run_target(target, result, bounds, screen))\n" + source[output_start:]
    return source.replace("def main():", function + "def main():")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pilot", type=Path)
    parser.add_argument("predictions", type=Path)
    parser.add_argument("community_ephemerides", type=Path, help=f"CSV from commit {COMMIT}")
    args = parser.parse_args()
    pilot = args.pilot
    snapshot = ROOT / "data/snapshot"
    snapshot.mkdir(parents=True, exist_ok=True)
    files = []
    for name in (*ARTIFACTS, "predictions.csv"):
        raw = (args.predictions if name == "predictions.csv" else pilot / name).read_bytes()
        filename = name + ".gz" if len(raw) > 12000 else name
        target = snapshot / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = gzip.compress(raw, mtime=0) if filename.endswith(".gz") else raw
        target.write_bytes(payload)
        files.append(dict(logical_path=name, path=filename, sha256=sha(payload), uncompressed_sha256=sha(raw)))
    (ROOT / "data/snapshot-manifest.json").write_bytes(json_bytes(files))
    for name in CODE:
        output = ROOT / "legacy" / name
        output.parent.mkdir(parents=True, exist_ok=True)
        if name == "vet.py":
            output.write_text(adapt_vetter((pilot / name).read_text(encoding="utf-8")), encoding="utf-8")
        elif name == "screen.py":
            original = (pilot / name).read_text(encoding="utf-8")
            original = original.replace(
                'def read_curve(product, raw=False):\n    path = ROOT / "data/curves" / product["dataURL"].rsplit("/", 1)[-1]',
                'def read_curve(product, raw=False, *, path=None):\n    if path is None:\n        path = ROOT / "data/curves" / product["dataURL"].rsplit("/", 1)[-1]',
            )
            output.write_text(original, encoding="utf-8")
        else:
            shutil.copyfile(pilot / name, output)
    cohort = json.loads((pilot / "cohort.json").read_bytes())
    ids = {x["catalog"]["TICID"] for x in cohort}
    catalogs = [read_catalog(pilot / "data" / source[0], source, ids) for source in SOURCES]
    catalogs.append(read_catalog(args.community_ephemerides,
                                 ("ephemerides.csv", "TIC", "period_d", "candidate_id", "community", COMMUNITY), ids))
    (ROOT / "data/catalogs.json").write_bytes(json_bytes(catalogs))
    uris = sorted({p["dataURL"] for row in cohort for p in row["primary"]})
    hashes = {}
    for uri in uris:
        file = pilot / "data/curves" / uri.rsplit("/", 1)[-1]
        if file.is_file():
            hashes[uri] = dict(sha256=sha(file.read_bytes()), bytes=file.stat().st_size)
        else:
            hashes[uri] = None
    (ROOT / "data/fits-manifest.json.gz").write_bytes(gzip.compress(json_bytes(hashes), mtime=0))
    print(f"{len(ids)} targets, {len(uris)} primary files, {len(files)} frozen artifacts")


if __name__ == "__main__":
    main()

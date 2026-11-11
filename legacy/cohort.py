"""Select a fixed pilot cohort from catalog matches and MAST observations."""
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
PRIORITY = {"SPOC": 0, "TESS-SPOC": 1, "QLP": 2, "T16": 3}


def read(name):
    with (DATA / name).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def tic(value):
    return str(int(re.search(r"\d+", str(value)).group()))


def main():
    mono = read("singleTransit_candidate_table.csv")
    counts = Counter(row["TICID"] for row in mono)
    catalogs = {
        "TOI": ("toi_archive.csv", "tid"),
        "CTOI": ("exofop_ctoi.csv", "TIC ID"),
        "confirmed": ("confirmed_planets.csv", "tic_id"),
        "T16_periodic": ("t16_periodic.csv", "TICID"),
    }
    match = defaultdict(set)
    for label, (file, column) in catalogs.items():
        for row in read(file):
            if row.get(column, "").strip():
                match[tic(row[column])].add(label)
    observations = defaultdict(list)
    covered_ids = set()
    for path in sorted((DATA / "mast_batches").glob("*.json")):
        batch = json.loads(path.read_text(encoding="utf-8"))
        covered_ids.update(batch["ids"])
        for obs in batch["observations"]:
            if obs.get("dataRights") == "PUBLIC" and (obs.get("dataURL") or "").endswith(("_lc.fits", "_llc.fits")):
                observations[tic(obs["target_name"])].append(obs)
    if covered_ids != counts.keys():
        raise RuntimeError(f"Catalog inventory incomplete: {len(covered_ids)}/{len(counts)} TICs")

    eligible = []
    audit = []
    for row in mono:
        target = row["TICID"]
        sector = int(row["Sector"])
        products = defaultdict(list)
        for obs in observations[target]:
            if obs.get("sequence_number") is not None:
                products[int(obs["sequence_number"])].append(obs)
        originals = [p for p in products[sector] if p["provenance_name"] == "T16"
                     and f'cam{row["Camera"]}-ccd{row["CCD"]}' in p["obs_id"]]
        later = [s for s, pp in products.items() if s > sector and any(p["provenance_name"] != "T16" for p in pp)]
        reason = (
            "multiple_mono_rows" if counts[target] > 1 else
            "catalog_match" if match[target] else
            "no_original_t16" if not originals else
            "insufficient_later_products" if len(later) < 2 or max(later) < 27 else "eligible"
        )
        audit.append(dict(**row, exclusion=reason, matches=";".join(sorted(match[target])),
                          later_sectors=";".join(map(str, sorted(later)))))
        if reason != "eligible":
            continue
        ordered = {s: sorted(pp, key=lambda p: (PRIORITY[p["provenance_name"]], p.get("t_exptime") or 0, p["obs_id"]))
                   for s, pp in sorted(products.items())}
        primary = [pp[0] for s, pp in ordered.items() if s != sector and pp]
        primary.append(sorted(originals, key=lambda p: p["dataURL"])[0])
        alternate = [p for p in ordered[sector] if p["provenance_name"] != "T16"]
        if alternate:
            primary.append(alternate[0])
        eligible.append(dict(catalog=row, later_count=len(later), primary=primary, all_products=ordered))
    eligible.sort(key=lambda r: (-r["later_count"], -float(r["catalog"]["Duration"]), int(r["catalog"]["TICID"])))
    if len(eligible) < 24:
        raise RuntimeError(f"Only {len(eligible)} eligible targets, review the feasibility gate")
    picked = eligible[:24]
    with (ROOT / "selection_audit.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(audit[0]))
        writer.writeheader()
        writer.writerows(audit)
    (ROOT / "cohort.json").write_text(json.dumps(picked, indent=2), encoding="utf-8")
    print(Counter(row["exclusion"] for row in audit))
    print([(r["catalog"]["TICID"], r["later_count"]) for r in picked])


if __name__ == "__main__":
    main()

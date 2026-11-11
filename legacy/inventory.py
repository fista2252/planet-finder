"""Snapshot MAST light-curve observations for the T16 single-transit catalog."""
import csv
import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent
MAST = "https://mast.stsci.edu/api/v0/invoke"


def invoke(service, params, page=1):
    query = dict(service=service, params=params, format="json", pagesize=2000, page=page)
    response = requests.post(MAST, data={"request": json.dumps(query)}, timeout=65)
    response.raise_for_status()
    result = response.json()
    if result.get("status") != "COMPLETE":
        raise RuntimeError(f"{service}: {result.get('status')} {result.get('msg')}")
    return result


def inventory(tics):
    filters = [
        {"paramName": "obs_collection", "values": ["TESS", "HLSP"]},
        {"paramName": "target_name", "values": tics + [tic.zfill(16) for tic in tics]},
        {"paramName": "provenance_name", "values": ["QLP", "SPOC", "TESS-SPOC", "T16"]},
        {"paramName": "dataproduct_type", "values": ["timeseries"]},
    ]
    page = 1
    while True:
        result = invoke("Mast.Caom.Filtered", {"columns": "*", "filters": filters}, page)
        yield from result["data"]
        if len(result["data"]) < 2000:
            break
        page += 1


def main():
    with (ROOT / "data/singleTransit_candidate_table.csv").open(newline="") as f:
        tics = sorted({row["TICID"] for row in csv.DictReader(f)}, key=int)
    assert len(tics) >= 20
    cache = ROOT / "data/mast_batches"
    cache.mkdir(exist_ok=True)
    for i in range(0, len(tics), 40):
        ids = tics[i:i + 40]
        path = cache / f"{i:03}.json"
        if path.exists():
            continue
        observations = list(inventory(ids))
        path.write_text(json.dumps(dict(ids=ids, observations=observations)), encoding="utf-8")
        print(i, len(ids), len(observations), flush=True)


if __name__ == "__main__":
    main()

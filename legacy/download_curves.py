"""Cache only the selected public light-curve products."""
import concurrent.futures
import json
import threading
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "data/curves"
DOWNLOAD = "https://mast.stsci.edu/api/v0.1/Download/file"
FILE_LIMIT = 32 * 1024 * 1024
PILOT_LIMIT = 1_500_000_000
LOCAL = threading.local()


def open_product(uri):
    if not hasattr(LOCAL, "session"):
        LOCAL.session = requests.Session()
    if uri.startswith("mast:HLSP/qlp/"):
        key = uri.replace("mast:HLSP/", "mast/hlsp/", 1)
        try:
            mirror = LOCAL.session.get(f"https://stpubdata.s3.amazonaws.com/{key}", stream=True, timeout=(10, 18))
            if mirror.ok:
                return mirror
            mirror.close()
        except requests.RequestException:
            pass
    return LOCAL.session.get(DOWNLOAD, params={"uri": uri}, stream=True, timeout=(15, 35))


def download(uri):
    target = CACHE / uri.rsplit("/", 1)[-1]
    if target.exists():
        return dict(uri=uri, bytes=target.stat().st_size, status="cached")
    temporary = target.with_suffix(".part")
    try:
        with open_product(uri) as response:
            response.raise_for_status()
            total = 0
            with temporary.open("wb") as f:
                for block in response.iter_content(128 * 1024):
                    if total == 0 and not block.startswith(b"SIMPLE  ="):
                        raise ValueError("Response is not a FITS file")
                    total += len(block)
                    if total > FILE_LIMIT:
                        raise ValueError("Exceeds the per-product size bound")
                    f.write(block)
            if total == 0 or total % 2880:
                raise ValueError("Empty or incomplete FITS blocks")
        temporary.replace(target)
        return dict(uri=uri, bytes=total, status="downloaded")
    except (OSError, ValueError, requests.RequestException) as error:
        temporary.unlink(missing_ok=True)
        return dict(uri=uri, bytes=0, status=str(error))


def main():
    selected = json.loads((ROOT / "cohort.json").read_text(encoding="utf-8"))
    uris = list(dict.fromkeys(p["dataURL"] for target in selected for p in target["primary"]))
    if len(uris) > 650:
        raise RuntimeError("Cohort exceeds the 650 light-curve transfer limit")
    CACHE.mkdir(exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
        results = []
        for i in range(0, len(uris), 16):
            total = sum(path.stat().st_size for path in CACHE.glob("*.fits"))
            batch = uris[i:i+16]
            new = sum(not (CACHE / uri.rsplit("/", 1)[-1]).exists() for uri in batch)
            if total + new * FILE_LIMIT > PILOT_LIMIT:
                print("Stopped at conservative archive transfer bound", flush=True)
                break
            results.extend(pool.map(download, batch))
            (ROOT / "download_manifest.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
            print(f"{len(results)}/{len(uris)} products attempted; cache {total:,} bytes before batch", flush=True)


if __name__ == "__main__":
    main()

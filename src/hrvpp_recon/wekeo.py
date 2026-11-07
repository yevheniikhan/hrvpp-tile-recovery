"""WEkEO HDA download, faithful to MONIFUN `download_copernicus_pheno.py`.

Differences to the original script (all additive, behaviour is unchanged):
* credentials can come from Colab Secrets as well as environment variables;
* `probe_bbox()` builds a tiny bbox around a point so that only one neighbouring tile is
  returned (used for validation, where only the tiles overlapping the gap are needed);
* SHA-256 checksums are written so repeated downloads can be compared (reproducibility).
"""
from __future__ import annotations

import csv
import hashlib
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, List

from .config import HDA_DATASET_ID, HDA_PRODUCT_GROUP

MAX_PASSES_DEFAULT = 5
SLEEP_BETWEEN_PASSES_DEFAULT = 10


def load_credentials() -> bool:
    """Populate WEKEO_USER / WEKEO_PASSWORD from Colab Secrets if not already set."""
    if os.getenv("WEKEO_USER") and os.getenv("WEKEO_PASSWORD"):
        return True
    try:
        from google.colab import userdata  # type: ignore
        for k in ("WEKEO_USER", "WEKEO_PASSWORD"):
            try:
                v = userdata.get(k)
            except Exception:
                v = None
            if v:
                os.environ[k] = v
    except ImportError:
        pass
    return bool(os.getenv("WEKEO_USER") and os.getenv("WEKEO_PASSWORD"))


def make_client(workers: int = 1):
    from hda import Client, Configuration  # imported lazily

    user, password = os.getenv("WEKEO_USER", ""), os.getenv("WEKEO_PASSWORD", "")
    if user and password:
        return Client(config=Configuration(user=user, password=password), max_workers=workers)
    return Client(max_workers=workers)


def query(year: int, product: str, bbox: List[float]) -> dict:
    return {
        "dataset_id": HDA_DATASET_ID, "productType": product, "productGroupId": HDA_PRODUCT_GROUP,
        "start": f"{year}-01-01T00:00:00.000Z", "end": f"{year}-12-31T00:00:00.000Z", "bbox": bbox,
    }


def _result_id(item) -> str:
    if isinstance(item, dict):
        return str(item.get("id") or item.get("productIdentifier") or item)
    return str(item)


def search_ids(year: int, product: str, bbox: List[float], workers: int = 1) -> List[str]:
    client = make_client(workers)
    return [_result_id(r) for r in client.search(query(year, product, bbox)).results]


@dataclass
class ResultRow:
    year: int
    productType: str
    expected_n: int
    downloaded_n: int
    passes_used: int
    seconds: float
    path: str
    status: str


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def probe_bbox(lon: float, lat: float, d: float = 0.005) -> List[float]:
    return [lon - d, lat - d, lon + d, lat + d]


def download_one(*, year: int, product: str, outdir: Path, bbox: List[float], workers: int = 1,
                 max_passes: int = MAX_PASSES_DEFAULT,
                 sleep_between_passes: int = SLEEP_BETWEEN_PASSES_DEFAULT) -> ResultRow:
    """Same contract as MONIFUN download_one(): files land in OUTDIR/PRODUCT/YEAR."""
    t0 = time.time()
    client = make_client(workers)
    matches = client.search(query(year, product, bbox))
    expected_n = len(matches.results)
    path = outdir / product / str(year)
    path.mkdir(parents=True, exist_ok=True)
    count = lambda: sum(1 for p in path.iterdir() if p.is_file() and not p.name.startswith("."))  # noqa: E731
    downloaded_n, passes_used = count(), 0
    if expected_n == 0:
        status = "no matches"
    else:
        for p in range(1, max_passes + 1):
            passes_used = p
            matches.download(download_dir=str(path))
            downloaded_n = count()
            print(f"[{year} {product}] pass {p}/{max_passes}: {downloaded_n} of {expected_n} files present")
            if downloaded_n >= expected_n:
                break
            if p < max_passes:
                time.sleep(sleep_between_passes)
        status = "ok" if downloaded_n >= expected_n else f"incomplete: {downloaded_n}/{expected_n}"
    return ResultRow(year, product, expected_n, downloaded_n, passes_used, round(time.time() - t0, 2),
                     str(path), status)


def download_all(outdir: Path, years: Iterable[int], products: Iterable[str], bbox: List[float],
                 workers: int = 1) -> List[ResultRow]:
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for y in years:
        for prod in products:
            try:
                row = download_one(year=y, product=prod, outdir=outdir, bbox=bbox, workers=workers)
            except Exception as e:  # keep going, as in the MONIFUN script
                row = ResultRow(y, prod, -1, -1, 0, 0.0, str(outdir / prod / str(y)), f"error: {e}")
            print(row)
            rows.append(row)
    with (outdir / "download_summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))
    return rows


def write_checksums(root: Path, out_csv: Path) -> int:
    files = sorted(p for p in root.rglob("*.tif"))
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "bytes", "sha256"])
        for p in files:
            w.writerow([str(p.relative_to(root)), p.stat().st_size, sha256(p)])
    return len(files)

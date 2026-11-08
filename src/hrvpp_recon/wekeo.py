"""WEkEO HDA download, faithful to MONIFUN `download_copernicus_pheno.py`.

Differences to the original script (all additive, behaviour is unchanged):
* credentials can come from Colab Secrets as well as environment variables;
* `probe_bbox()` builds a tiny bbox around a point so that only one neighbouring tile is
  returned (used for validation, where only the tiles overlapping the gap are needed);
* SHA-256 checksums are written so repeated downloads can be compared (reproducibility);
* every downloaded GeoTIFF is opened and its last row is read; truncated files (an interrupted
  HTTP stream leaves a partial file that the original script counts as "present") are deleted
  and fetched again in the next pass. Observed on 08/11/2026: T34UFV TPROD 2021 stopped at
  78.6 of 211 MB and was reported as complete by the file-count logic.
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


def is_valid_tif(path: Path) -> bool:
    """True if the GeoTIFF opens and its last row can be decoded (detects truncated downloads)."""
    try:
        import rasterio
        from rasterio.windows import Window
        with rasterio.open(path) as src:
            src.read(1, window=Window(0, src.height - 1, src.width, 1))
        return True
    except Exception:
        return False


def _valid_count(path: Path, ids=None) -> int:
    """Count valid files; with `ids`, count only files belonging to this query's results."""
    n = 0
    for p in path.iterdir():
        if not p.is_file() or p.name.startswith("."):
            continue
        if p.suffix.lower() in (".tif", ".tiff") and not is_valid_tif(p):
            print(f"  removing truncated/corrupt file: {p.name} ({p.stat().st_size} bytes)")
            p.unlink()
            continue
        if ids and not any(p.stem.startswith(i) or i.startswith(p.stem) for i in ids):
            continue
        n += 1
    return n


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
    ids = {_result_id(r) for r in matches.results}
    ids = {i for i in ids if i.startswith("VPP_") or i.startswith("ST_")} or None
    downloaded_n, passes_used = _valid_count(path, ids), 0
    if expected_n == 0:
        status = "no matches"
    else:
        for p in range(1, max_passes + 1):
            passes_used = p
            try:
                matches.download(download_dir=str(path))
            except Exception as e:  # interrupted stream: keep going, the check below decides
                print(f"[{year} {product}] pass {p}: download error: {str(e)[:150]}")
            downloaded_n = _valid_count(path, ids)
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

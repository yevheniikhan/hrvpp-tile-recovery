"""Anonymous availability checks against the VITO OpenSearch catalogue of HR-VPP.

No credentials are needed. Used to (re)confirm that a tile is absent from the
operational product before deciding to reconstruct it.
"""
from __future__ import annotations

import re
from typing import Iterable, List, Optional

import pandas as pd
import requests

from .config import VITO_COLLECTION_VPP_UTM, VITO_OPENSEARCH

TILE_RE = re.compile(r"_(T\d{2}[A-Z]{3})-")


def search(bbox: List[float], year: int, product_type: str = "AMPL", group: str = "s1",
           collection: str = VITO_COLLECTION_VPP_UTM, count: int = 200, timeout: int = 60) -> pd.DataFrame:
    """Return one row per catalogue product intersecting `bbox` (lon/lat) in `year`."""
    params = {
        "collection": collection, "bbox": ",".join(map(str, bbox)),
        "start": f"{year}-01-01", "end": f"{year}-12-31", "count": count,
    }
    if product_type:
        params["productType"] = product_type
    r = requests.get(VITO_OPENSEARCH, params=params, timeout=timeout)
    r.raise_for_status()
    rows = []
    for f in r.json().get("features", []):
        pid = f["id"]
        if group and f"_{group}_" not in pid and "-st_" not in collection:
            continue
        m = TILE_RE.search(pid)
        data = (f["properties"].get("links", {}).get("data") or [{}])[0]
        rows.append({"year": year, "id": pid, "tile": m.group(1) if m else None,
                     "bbox": f.get("bbox"), "size_bytes": data.get("length"), "href": data.get("href")})
    return pd.DataFrame(rows)


def tile_availability(tile_id: str, bbox: List[float], years: Iterable[int],
                      product_type: str = "AMPL", collection: str = VITO_COLLECTION_VPP_UTM) -> pd.DataFrame:
    """Matrix year x {available, neighbours} for one tile."""
    out = []
    for y in years:
        df = search(bbox, y, product_type, collection=collection)
        tiles = sorted(set(df["tile"].dropna())) if len(df) else []
        out.append({"year": y, "tile": tile_id, "available": tile_id in tiles,
                    "tiles_in_bbox": ", ".join(tiles)})
    return pd.DataFrame(out)


def point_has_product(lon: float, lat: float, year: int, collection: str = VITO_COLLECTION_VPP_UTM,
                      product_type: Optional[str] = None) -> int:
    """Number of catalogue products covering a single point (0 => gap)."""
    d = 1e-4
    return len(search([lon - d, lat - d, lon + d, lat + d], year, product_type or "", group="",
                      collection=collection))

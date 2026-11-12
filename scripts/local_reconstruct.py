#!/usr/bin/env python3
"""Step 3b. Reconstruct without Earth Engine: Sentinel-2 L2A from Microsoft Planetary Computer, NumPy or CuPy (GPU).

    pip install -r requirements-local.txt          # pystac-client, planetary-computer (+ cupy for a GPU)
    python scripts/local_reconstruct.py --region T34UGVE --years 2021 --out data/raw_local
    python scripts/local_reconstruct.py --region T34UGV --years 2024 --cpu --block-rows 128

Same algorithm and the same output bands as the Earth Engine export (see src/hrvpp_recon/s2local.py for the
known differences: SCL cloud mask instead of Cloud Score+, no L1C fallback). A whole region and year reads
~400-650 scenes per block; on a Colab L4 GPU the 7000 x 16980 px extension takes about 1.5 h per year,
dominated by reading. If the GPU runs out of memory, lower --gpu-chunk.
"""
import argparse
from pathlib import Path

from _common import log, parse_years
from hrvpp_recon import config, s2local

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
g = ap.add_mutually_exclusive_group(required=True)
g.add_argument("--region", choices=sorted(config.REGIONS))
g.add_argument("--grid", help="EPSG,x0,y1,width,height")
ap.add_argument("--name", default="AREA")
ap.add_argument("--years", required=True)
ap.add_argument("--out", default="data/raw_local")
ap.add_argument("--mdvi-bounds", help="x0,y0,x1,y1 of the area for the tile-constant MDVI term (default: the region)")
ap.add_argument("--block-rows", type=int, default=192)
ap.add_argument("--gpu-chunk", type=int, default=100_000)
ap.add_argument("--workers", type=int, default=32, help="parallel scene reads")
ap.add_argument("--cpu", action="store_true", help="NumPy instead of CuPy")
a = ap.parse_args()

if a.region:
    grid = config.REGIONS[a.region]
else:
    e, x0, y1, w, h = a.grid.split(",")
    grid = config.Grid(a.name, int(e), float(x0), float(y1), int(w), int(h))
db = tuple(float(v) for v in a.mdvi_bounds.split(",")) if a.mdvi_bounds else grid.bounds
out_dir = Path(a.out)
out_dir.mkdir(parents=True, exist_ok=True)
for y in parse_years(a.years):
    out = out_dir / f"LOCAL_VPP_{y}_S2_{grid.name}-010m_s1.tif"
    if out.exists():
        log(out, "exists, skipped")
        continue
    tmp = out.with_name("_" + out.name)
    s2local.run_region(y, grid.crs, grid.x0, grid.y1, grid.width, grid.height, str(tmp), delta_bounds=db,
                       block_rows=a.block_rows, workers=a.workers, use_gpu=not a.cpu, gpu_chunk=a.gpu_chunk, log=log)
    tmp.replace(out)
    log(y, "->", out)

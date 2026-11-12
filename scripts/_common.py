"""Shared helpers for the command-line scripts: import path, year parsing, logging."""
import os
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))


def parse_years(s: str):
    """'2017-2020,2023' -> [2017, 2018, 2019, 2020, 2023]"""
    out = []
    for part in str(s).replace(" ", "").split(","):
        if "-" in part:
            a, b = map(int, part.split("-"))
            out += list(range(a, b + 1))
        elif part:
            out.append(int(part))
    return out


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def data_root(arg=None) -> Path:
    """Working folder for downloads and outputs: --data, else $HRVPP_DATA, else ./data."""
    return Path(arg or os.getenv("HRVPP_DATA", "data")).expanduser()

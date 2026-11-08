"""Fetch Earth Engine Drive exports when the Earth Engine account differs from the Colab Drive account.

Earth Engine writes `Export.image.toDrive` results to the Drive of the account that authenticated
Earth Engine. If Colab's mounted Drive belongs to another account, the files are not visible under
/content/drive. This helper downloads them through the Drive API with the Earth Engine credentials.
The Drive API must be enabled in the quota project (console: APIs & Services > Google Drive API).
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import List

import requests


def _token() -> str:
    import google.auth.transport.requests
    from google.oauth2.credentials import Credentials

    c = json.load(open(os.path.expanduser("~/.config/earthengine/credentials")))
    creds = Credentials(None, refresh_token=c["refresh_token"], client_id=c["client_id"],
                        client_secret=c["client_secret"], token_uri="https://oauth2.googleapis.com/token")
    creds.refresh(google.auth.transport.requests.Request())
    return creds.token


def fetch(name_prefix: str, dest: Path, quota_project: str) -> List[Path]:
    """Download every non-trashed Drive file whose name starts with `name_prefix` into `dest`."""
    h = {"Authorization": "Bearer " + _token(), "X-Goog-User-Project": quota_project}
    q = f"name contains '{name_prefix}' and trashed = false"
    r = requests.get("https://www.googleapis.com/drive/v3/files", headers=h,
                     params={"q": q, "fields": "files(id,name,size)"}, timeout=60)
    r.raise_for_status()
    dest.mkdir(parents=True, exist_ok=True)
    out = []
    for f in r.json().get("files", []):
        if not f["name"].startswith(name_prefix):
            continue
        p = dest / f["name"]
        if p.exists() and p.stat().st_size == int(f.get("size", -1)):
            out.append(p)
            continue
        with requests.get(f"https://www.googleapis.com/drive/v3/files/{f['id']}", headers=h,
                          params={"alt": "media"}, stream=True, timeout=600) as resp:
            resp.raise_for_status()
            with open(p, "wb") as o:
                for chunk in resp.iter_content(1 << 22):
                    o.write(chunk)
        out.append(p)
    return out

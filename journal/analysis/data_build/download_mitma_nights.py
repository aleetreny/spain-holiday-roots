"""Download the MITMA municipal overnight-stay daily files for 2022-2024 that are not yet local.

Source: the official MITMA open-data RSS index (https://movilidad-opendata.mitma.es/RSS.xml).
Only standard daily municipal files are taken (YYYYMMDD_Pernoctaciones_municipios.csv.gz);
2025 is excluded because of the methodological break of 1 July 2025.
Files already in the main checkout (_local/external/mitma_pernoct) are skipped, never re-downloaded.
Output: <journal worktree>/_local/external_journal/mitma_pernoct/ (not in git) + manifest.csv.

Usage: python -I journal/analysis/data_build/download_mitma_nights.py
"""
from __future__ import annotations

import concurrent.futures as cf
import csv
import hashlib
import re
import time
import urllib.request
import os
from pathlib import Path

JOURNAL = Path(__file__).resolve().parents[3]
MAIN = Path(os.environ.get("SHR_MAIN_ROOT", JOURNAL)).resolve()
LOCAL_RAW = MAIN / "_local" / "external" / "mitma_pernoct"
OUT = JOURNAL / "_local" / "external_journal" / "mitma_pernoct"
RSS = "https://movilidad-opendata.mitma.es/RSS.xml"
YEARS = ("2022", "2023", "2024")
PAT = re.compile(r"https://[^<\"\s]+por-municipios/pernoctaciones/ficheros-diarios/[^<\"\s]+/(\d{8})_Pernoctaciones_municipios\.csv\.gz")
UA = {"User-Agent": "research-client/1.0"}


def fetch(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rss = fetch(RSS).decode("utf-8", "replace")
    urls = {}
    for m in PAT.finditer(rss):
        d = m.group(1)
        if d[:4] in YEARS:
            urls[d] = m.group(0)
    have_main = {p.name[:8] for p in LOCAL_RAW.glob("*_Pernoctaciones_municipios.csv.gz")}
    todo = [(d, u) for d, u in sorted(urls.items())
            if d not in have_main and not (OUT / f"{d}_Pernoctaciones_municipios.csv.gz").exists()]
    print(f"listed {len(urls)}; in main checkout {len(have_main & set(urls))}; to download {len(todo)}", flush=True)

    def one(item):
        d, u = item
        dest = OUT / f"{d}_Pernoctaciones_municipios.csv.gz"
        tmp = dest.with_suffix(".part")
        for attempt in range(1, 6):
            try:
                data = fetch(u)
                if len(data) < 10_000:
                    raise RuntimeError(f"too small ({len(data)} B)")
                tmp.write_bytes(data)
                tmp.replace(dest)
                return d, u, "OK", len(data), hashlib.sha256(data).hexdigest()
            except Exception as e:  # noqa: BLE001
                if attempt == 5:
                    return d, u, f"FAIL {e}", 0, ""
                time.sleep(2 ** attempt)

    rows = []
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        for k, res in enumerate(ex.map(one, todo), 1):
            rows.append(res)
            if k % 25 == 0 or k == len(todo):
                print(f"{k}/{len(todo)} done, {time.time() - t0:.0f}s", flush=True)
    man = OUT / "manifest.csv"
    new = not man.exists()
    with man.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["file_date", "url", "status", "bytes", "sha256", "downloaded"])
        stamp = time.strftime("%Y-%m-%d")
        for r in rows:
            w.writerow([*r, stamp])
    fails = [r for r in rows if r[2] != "OK"]
    print(f"OK {len(rows) - len(fails)}, FAIL {len(fails)}")
    for r in fails:
        print("  ", r[0], r[2])


if __name__ == "__main__":
    main()

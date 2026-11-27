"""Download MITMA municipal overnight-stay files for 2025 (Confirmation II, pre-registration addendum).

Takes every standard daily municipal file dated 2025 and the six method-change recalculations
(folder cambio_metodológico_2025/por_municipios/pernoctaciones). Output (not in git):
<journal>/_local/external_journal/mitma_pernoct_2025/{standard,recalc}/ + manifest.csv.
"""
from __future__ import annotations

import concurrent.futures as cf
import csv
import hashlib
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

JOURNAL = Path(__file__).resolve().parents[3]
OUT = JOURNAL / "_local" / "external_journal" / "mitma_pernoct_2025"
RSS = "https://movilidad-opendata.mitma.es/RSS.xml"
UA = {"User-Agent": "research-client/1.0"}
STD = re.compile(r"https://[^<\"\s]+por-municipios/pernoctaciones/ficheros-diarios/2025-\d{2}/(2025\d{4})_Pernoctaciones_municipios\.csv\.gz")
REC = re.compile(r"https://[^<\"\s]+cambio_metodol[^/]+/por_municipios/pernoctaciones/(\d{8})_Pernoctaciones_municipios\.csv\.gz")


def fetch(url: str, timeout: int = 120) -> bytes:
    parts = urllib.parse.urlsplit(url)
    safe = urllib.parse.urlunsplit(parts._replace(path=urllib.parse.quote(parts.path)))
    with urllib.request.urlopen(urllib.request.Request(safe, headers=UA), timeout=timeout) as r:
        return r.read()


def main() -> None:
    rss = fetch(RSS).decode("utf-8", "replace")
    items = {("standard", m.group(1)): m.group(0) for m in STD.finditer(rss)}
    items |= {("recalc", m.group(1)): m.group(0) for m in REC.finditer(rss)}
    todo = []
    for (kind, d), u in sorted(items.items()):
        dest = OUT / kind / f"{d}_Pernoctaciones_municipios.csv.gz"
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            todo.append((kind, d, u, dest))
    print(f"listed standard {sum(k == 'standard' for k, _ in items)}, recalc {sum(k == 'recalc' for k, _ in items)}; to download {len(todo)}", flush=True)

    def one(item):
        kind, d, u, dest = item
        for attempt in range(1, 6):
            try:
                data = fetch(u)
                if len(data) < 10_000:
                    raise RuntimeError(f"too small ({len(data)} B)")
                tmp = dest.with_suffix(".part")
                tmp.write_bytes(data)
                tmp.replace(dest)
                return kind, d, u, "OK", len(data), hashlib.sha256(data).hexdigest()
            except Exception as e:  # noqa: BLE001
                if attempt == 5:
                    return kind, d, u, f"FAIL {e}", 0, ""
                time.sleep(2 ** attempt)

    rows = []
    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        for k, res in enumerate(ex.map(one, todo), 1):
            rows.append(res)
            if k % 25 == 0 or k == len(todo):
                print(f"{k}/{len(todo)}", flush=True)
    man = OUT / "manifest.csv"
    new = not man.exists()
    with man.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["kind", "file_date", "url", "status", "bytes", "sha256", "downloaded"])
        for r in rows:
            w.writerow([*r, time.strftime("%Y-%m-%d")])
    fails = [r for r in rows if r[3] != "OK"]
    print(f"OK {len(rows) - len(fails)}, FAIL {len(fails)}")


if __name__ == "__main__":
    main()

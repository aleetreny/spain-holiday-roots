"""Fetch only reviewed MITMA manifest URLs, with at most three workers."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
import subprocess

folder = Path('_local/summer_territories_downloads')
folder.mkdir(parents=True, exist_ok=True)
manifest = Path('metadata_cache/data_links_v2_2026-11-23.xml').read_text()
urls = sorted(set(re.findall(
    r'https://movilidad-opendata\.mitma\.es/estudios_basicos/por-municipios/pernoctaciones/'
    r'ficheros-diarios/(?:2024-07|202[23]-08)/\d{8}_Pernoctaciones_municipios\.csv\.gz',
    manifest,
)))
# Full July and a matched seven-night sensitivity window in earlier Augusts.
urls = [u for u in urls if '2024-07' in u or 12 <= int(u.rsplit('/',1)[1][6:8]) <= 18]

def download(url):
    destination = folder / url.rsplit('/', 1)[1]
    if not destination.exists():
        temporary = destination.with_suffix('.part')
        subprocess.run(['curl', '-fLsS', '--retry', '3', '--retry-all-errors', '--http1.1', '--max-time', '90',
                        '-o', str(temporary), url], check=True)
        subprocess.run(['gzip', '-t', str(temporary)], check=True)
        temporary.replace(destination)

with ThreadPoolExecutor(max_workers=3) as pool:
    jobs = [pool.submit(download, url) for url in urls]
    failed = [url for url, job in zip(urls, jobs) if job.exception()]
    if failed:
        raise SystemExit("Missing files:\n" + "\n".join(failed))

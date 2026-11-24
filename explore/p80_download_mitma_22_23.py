# Download MITMA municipal overnight-stay daily files for 2022 and 2023
# (Feb, Aug, Nov 1 / Oct 31, Dec 24) from metadata_cache/data_links_v2_2026-11-23.xml
# Max 3 parallel connections with retries, saving into _local/external/mitma_pernoct/
import os, re, time, concurrent.futures, requests

XML_PATH = 'metadata_cache/data_links_v2_2026-11-23.xml'
OUT_DIR = '_local/external/mitma_pernoct'
os.makedirs(OUT_DIR, exist_ok=True)

txt = open(XML_PATH).read()
all_urls = sorted(set(re.findall(
    r'https://[^\s<]+por-municipios/pernoctaciones/ficheros-diarios/[^\s<]+csv\.gz', txt
)))

# Select: all of Feb & Aug for 2022 and 2023; All Saints window (Oct 29 - Nov 4) for 2022 & 2023; Dec 24 for 2022 & 2023
selected = []
for u in all_urls:
    m = re.search(r'/(\d{8})_Pernoctaciones_municipios\.csv\.gz$', u)
    if not m:
        continue
    d = m.group(1)
    yr, ym, md = d[:4], d[:6], d[4:]
    if yr in ('2022', '2023'):
        if ym in ('202202', '202208', '202302', '202308') or md in ('1029', '1030', '1031', '1101', '1102', '1224'):
            selected.append((d, u))

print(f"Matched {len(selected)} target files for 2022-2023.")
to_dl = []
for d, u in selected:
    fn = os.path.join(OUT_DIR, os.path.basename(u))
    if os.path.exists(fn) and os.path.getsize(fn) > 10000:
        continue
    to_dl.append((d, u, fn))
print(f"Already on disk: {len(selected) - len(to_dl)}; to download: {len(to_dl)}")

def download_one(item):
    d, u, fn = item
    tmp = fn + '.part'
    for attempt in range(1, 5):
        try:
            with requests.get(u, stream=True, timeout=(15, 30), headers={'User-Agent': 'Spain-Mobility-Research/1.0'}) as r:
                r.raise_for_status()
                with open(tmp, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1 << 16):
                        if chunk:
                            f.write(chunk)
            if os.path.getsize(tmp) < 1000:
                raise RuntimeError("Downloaded file suspiciously small")
            os.replace(tmp, fn)
            return (d, 'OK', os.path.getsize(fn))
        except Exception as e:
            if os.path.exists(tmp):
                try: os.remove(tmp)
                except OSError: pass
            if attempt < 4:
                time.sleep(2 ** attempt)
            else:
                return (d, f'FAIL: {e}', 0)

t0 = time.time()
ok = 0
fail = []
total_bytes = 0
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
    futs = [ex.submit(download_one, item) for item in to_dl]
    for fut in concurrent.futures.as_completed(futs):
        d, st, sz = fut.result()
        if st == 'OK':
            ok += 1
            total_bytes += sz
            if ok % 15 == 0 or ok == len(to_dl):
                print(f"  Progress: {ok}/{len(to_dl)} OK ({total_bytes/1e6:.1f} MB in {time.time()-t0:.1f}s)", flush=True)
        else:
            fail.append((d, st))
            print(f"  ERROR {d}: {st}", flush=True)

print(f"Done in {time.time()-t0:.1f}s. OK={ok}, FAIL={len(fail)}, Total={total_bytes/1e6:.1f} MB")
if fail:
    print("Failed:", fail)

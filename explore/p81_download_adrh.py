# Download INE Atlas de Distribución de Renta de los Hogares (ADRH, Op 353)
# All 52 Spanish provinces: municipal and district average/median income (2015-2024).
# Saves raw gzipped CSVs (municipal+district rows) in _local/external/ine_adrh/raw/
# and clean tidy CSVs + README.md in _local/external/ine_adrh/.
import csv, gzip, io, os, re, time, requests

OUT_DIR = '_local/external/ine_adrh'
RAW_DIR = os.path.join(OUT_DIR, 'raw')
os.makedirs(RAW_DIR, exist_ok=True)

# Fetch exact table list from INE Tempus3 API
sess = requests.Session()
sess.headers.update({
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) Spain-Mobility-Research/1.0',
    'Accept-Encoding': 'gzip, deflate',
})

tabs = sess.get('https://servicios.ine.es/wstempus/js/ES/TABLAS_OPERACION/353', timeout=30).json()
# Provincial tables have 'Indicadores de renta media' and exclude national summary 53689
cand_ids = sorted(set(
    t['Id'] for t in tabs
    if 'indicadores de renta media' in t.get('Nombre', '').lower() and t['Id'] != 53689
))
print(f"Found {len(cand_ids)} provincial ADRH tables: {cand_ids[:5]} ... {cand_ids[-5:]}", flush=True)

mun_rows = {}   # (ine_code, indicator, year) -> (name, value)
dist_rows = {}  # (dist_code, indicator, year) -> (mun_code, value)
provs_seen = {}

def parse_es_num(s):
    s = s.strip().strip('"')
    if not s or s in ('.', '..', '-', '""'):
        return None
    # Spanish formatting: thousands '.', decimal ','
    s = s.replace('.', '').replace(',', '.')
    try:
        return float(s)
    except ValueError:
        return None

t0 = time.time()
for i, tid in enumerate(cand_ids, 1):
    raw_path = os.path.join(RAW_DIR, f'adrh_{tid}_mundist.csv.gz')
    lines_kept = []
    if os.path.exists(raw_path) and os.path.getsize(raw_path) > 200:
        with gzip.open(raw_path, 'rt', encoding='utf-8-sig') as f:
            lines_kept = f.read().splitlines()
    else:
        url = f'https://www.ine.es/jaxiT3/files/t/csv_bdsc/{tid}.csv'
        ok = False
        for attempt in range(1, 5):
            try:
                with sess.get(url, stream=True, timeout=(20, 90)) as r:
                    r.raise_for_status()
                    # Stream line by line so we keep full municipal + district rows (skipping Secciones != '')
                    header = None
                    lines_kept = []
                    for raw_line in r.iter_lines(decode_unicode=False):
                        if not raw_line:
                            continue
                        line = raw_line.decode('utf-8-sig', errors='replace')
                        if header is None:
                            header = line
                            lines_kept.append(header)
                            continue
                        parts = line.split(';')
                        if len(parts) >= 6 and parts[2].strip() == '':
                            lines_kept.append(line)
                with gzip.open(raw_path, 'wt', encoding='utf-8') as f:
                    f.write('\n'.join(lines_kept) + '\n')
                ok = True
                break
            except Exception as e:
                print(f"  [retry {attempt}] tid={tid}: {e}", flush=True)
                time.sleep(3 * attempt)
        if not ok:
            print(f"  FAILED tid={tid}!", flush=True)
            continue
        time.sleep(0.25)

    # Parse kept lines
    if len(lines_kept) < 2:
        continue
    reader = csv.reader(lines_kept, delimiter=';')
    hdr = next(reader)
    prov_here = set()
    for row in reader:
        if len(row) < 6:
            continue
        mun_col, dist_col, sec_col, ind, yr, val_s = [x.strip() for x in row[:6]]
        if sec_col != '':
            continue
        val = parse_es_num(val_s)
        if val is None:
            continue
        if dist_col == '' and mun_col != '':
            m = re.match(r'^(\d{5})\s+(.*)$', mun_col)
            if m:
                code, name = m.group(1), m.group(2).strip()
                prov_here.add(code[:2])
                mun_rows[(code, ind, yr)] = (name, val)
        elif dist_col != '':
            # District column format: e.g. '0801901 Barcelona distrito 01'
            m = re.match(r'^(\d{7})\b', dist_col)
            if m:
                dcode = m.group(1)
                prov_here.add(dcode[:2])
                dist_rows[(dcode, ind, yr)] = (dcode[:5], dist_col, val)
    for p in prov_here:
        provs_seen.setdefault(p, []).append(tid)
    if i % 10 == 0 or i == len(cand_ids):
        print(f"  [{i}/{len(cand_ids)}] elapsed={time.time()-t0:.1f}s provs={len(provs_seen)} mun_obs={len(mun_rows)} dist_obs={len(dist_rows)}", flush=True)

print(f"\nTotal provinces covered: {len(provs_seen)}/52")
missing_provs = [f'{p:02d}' for p in range(1, 53) if f'{p:02d}' not in provs_seen]
print("Missing provinces:", missing_provs)

# Write clean CSVs
mun_out = os.path.join(OUT_DIR, 'adrh_municipios.csv')
with open(mun_out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['ine_code', 'municipality_name', 'indicator', 'year', 'value'])
    for (code, ind, yr), (name, val) in sorted(mun_rows.items()):
        w.writerow([code, name, ind, yr, val])

dist_out = os.path.join(OUT_DIR, 'adrh_distritos.csv')
with open(dist_out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['dist_code', 'ine_code', 'dist_label', 'indicator', 'year', 'value'])
    for (dcode, ind, yr), (mcode, lbl, val) in sorted(dist_rows.items()):
        w.writerow([dcode, mcode, lbl, ind, yr, val])

# Summary diagnostics
import pandas as pd
dm = pd.read_csv(mun_out, dtype={'ine_code': str})
dd = pd.read_csv(dist_out, dtype={'dist_code': str, 'ine_code': str})
print("\n=== ADRH Municipal Summary ===")
print("Rows:", len(dm), "Unique munis:", dm.ine_code.nunique(), "Years:", sorted(dm.year.unique()))
print("Indicators:", dm.indicator.unique().tolist())
for yr in sorted(dm.year.unique())[-3:]:
    sub = dm[(dm.year == yr) & (dm.indicator == 'Renta neta media por persona')]
    print(f"  Year {yr} (Renta neta media por persona): {len(sub)} munis, median={sub.value.median():.0f}, Madrid={sub[sub.ine_code=='28079'].value.values}, BCN={sub[sub.ine_code=='08019'].value.values}, Pozuelo={sub[sub.ine_code=='28115'].value.values}")

print("\n=== ADRH District Summary ===")
print("Rows:", len(dd), "Unique districts:", dd.dist_code.nunique(), "Unique munis:", dd.ine_code.nunique())
for city, cname in [('08019', 'Barcelona'), ('28079', 'Madrid')]:
    sub = dd[(dd.ine_code == city) & (dd.indicator == 'Renta neta media por persona') & (dd.year == dm.year.max())]
    print(f"  {cname} ({city}) {dm.year.max()} districts: n={len(sub)}, min={sub.value.min():.0f}, max={sub.value.max():.0f}")

# Write README.md
readme = f"""# INE Atlas de Distribución de Renta de los Hogares (ADRH)

- **Source:** Instituto Nacional de Estadística (INE), Operation `353` (*Atlas de distribución de renta de los hogares*, tables `Indicadores de renta media y mediana` for all 52 provinces: `{cand_ids[0]}`..`{cand_ids[-1]}`).
- **URL pattern:** `https://www.ine.es/jaxiT3/files/t/csv_bdsc/<table_id>.csv`
- **Downloaded:** {time.strftime('%Y-%m-%d')}
- **Processing script:** `explore/p81_download_adrh.py`
- **Raw files:** `raw/adrh_<table_id>_mundist.csv.gz` (original semicolon-delimited INE rows restricted to municipality and district aggregates, excluding census-tract `Secciones` rows to save disk space).
- **Clean outputs:**
  - `adrh_municipios.csv`: {len(dm):,} rows covering {dm.ine_code.nunique():,} Spanish municipalities across {len(provs_seen)} provinces, years {dm.year.min()}–{dm.year.max()}.
  - `adrh_distritos.csv`: {len(dd):,} rows covering {dd.dist_code.nunique():,} INE districts across {dd.ine_code.nunique():,} municipalities, years {dd.year.min()}–{dd.year.max()}.
- **Indicators included:** {', '.join(dm.indicator.unique().tolist())}.
- **License:** INE Open Data (free reuse with attribution).
"""
with open(os.path.join(OUT_DIR, 'README.md'), 'w', encoding='utf-8') as f:
    f.write(readme)
print("Wrote README.md and clean CSVs successfully.")

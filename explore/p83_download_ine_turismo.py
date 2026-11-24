# Download INE Experimental Statistics: "Medición del turismo a partir de teléfonos móviles" (Turismo interno)
# Tables:
#   52997: Province of origin x Province of destination (tourists, overnights, avg stay, monthly 2019-2026)
#   53000: Municipality of origin x Province of destination (tourists, monthly 2019-2026)
#   53001: Province of origin x Municipality of destination (tourists, monthly 2019-2026)
# Saves raw/filtered gzipped files and clean CSVs + README in _local/external/ine_turismo_moviles/.
import csv, gzip, os, re, time, requests

OUT_DIR = '_local/external/ine_turismo_moviles'
RAW_DIR = os.path.join(OUT_DIR, 'raw')
os.makedirs(RAW_DIR, exist_ok=True)

sess = requests.Session()
sess.headers.update({
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) Spain-Mobility-Research/1.0',
    'Accept-Encoding': 'gzip, deflate',
})

# Load 151 cities list from ine_birth_residence so we know all key origin municipalities
import pandas as pd
B = pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv',
                dtype={'ine_code': str})
CITIES_151 = set(B.ine_code.unique())

def parse_es_num(s):
    s = s.strip().strip('"')
    if not s or s in ('.', '..', '-', '""'):
        return None
    s = s.replace('.', '').replace(',', '.')
    try:
        return float(s)
    except ValueError:
        return None

def stream_download_table(tid, keep_fn, out_gz_path, desc):
    if os.path.exists(out_gz_path) and os.path.getsize(out_gz_path) > 5000:
        print(f"  [{tid}] Using cached {out_gz_path} ({os.path.getsize(out_gz_path)/1e6:.2f} MB)", flush=True)
        return
    url = f'https://www.ine.es/jaxiT3/files/t/csv_bdsc/{tid}.csv'
    print(f"  [{tid}] Downloading {desc} from {url} ...", flush=True)
    t0 = time.time()
    tmp = out_gz_path + '.part'
    for attempt in range(1, 5):
        try:
            n_total = 0
            n_kept = 0
            with sess.get(url, stream=True, timeout=(30, 180)) as r:
                r.raise_for_status()
                with gzip.open(tmp, 'wt', encoding='utf-8') as f:
                    header = None
                    for raw_line in r.iter_lines(decode_unicode=False):
                        if not raw_line:
                            continue
                        line = raw_line.decode('utf-8-sig', errors='replace')
                        if header is None:
                            header = line
                            f.write(header + '\n')
                            continue
                        n_total += 1
                        if keep_fn(line):
                            f.write(line + '\n')
                            n_kept += 1
                        if n_total % 2_000_000 == 0:
                            print(f"    [{tid}] scanned {n_total/1e6:.1f}M lines, kept {n_kept/1e6:.2f}M ({time.time()-t0:.1f}s)", flush=True)
            os.replace(tmp, out_gz_path)
            print(f"  [{tid}] DONE: scanned {n_total:,} lines, kept {n_kept:,} lines, {os.path.getsize(out_gz_path)/1e6:.2f} MB in {time.time()-t0:.1f}s", flush=True)
            return
        except Exception as e:
            if os.path.exists(tmp):
                try: os.remove(tmp)
                except OSError: pass
            print(f"  [{tid}] attempt {attempt} failed: {e}", flush=True)
            time.sleep(5 * attempt)
    raise RuntimeError(f"Failed to download table {tid}")

# 1. Table 52997: Province x Province (keep all non-empty rows!)
stream_download_table(
    52997,
    lambda line: not line.endswith(';') and not line.endswith(';""') and not line.endswith(';.'),
    os.path.join(RAW_DIR, 'turismo_52997_prov_x_prov.csv.gz'),
    'Province x Province domestic tourists & overnights',
)
time.sleep(1.0)

# 2. Table 53000: Municipality of origin x Province of destination
# Keep ALL non-empty rows for August (M08), November (M11), February (M02), December (M12), April (M04), March (M03)
# AND ALL months for the 151 Census cities!
TARGET_MONTHS = ('M02', 'M03', 'M04', 'M08', 'M11', 'M12')
def keep_53000(line):
    # Format: Municipio de origen;CCAA_dest;Prov_dest;...;Concepto;Periodo;Total
    if line.endswith(';') or line.endswith(';""') or line.endswith(';.'):
        return False
    parts = line.split(';')
    if len(parts) < 7:
        return False
    val = parts[6].strip().strip('"')
    if not val or val in ('.', '..'):
        return False
    per = parts[5].strip()
    mun5 = parts[0][:5]
    return (mun5 in CITIES_151) or any(per.endswith(m) for m in TARGET_MONTHS)

stream_download_table(
    53000,
    keep_53000,
    os.path.join(RAW_DIR, 'turismo_53000_munorig_x_provdest.csv.gz'),
    'Municipality of origin x Province of destination domestic tourists',
)
time.sleep(1.0)

# 3. Table 53001: Province of origin x Municipality of destination
# Keep ALL Total Nacional rows (for municipal village-refill across ALL 8,131 munis, ALL months 2019-2026!)
# PLUS provincial-origin rows for M02, M08, M11 of 2022-2025.
def keep_53001(line):
    if line.endswith(';') or line.endswith(';""') or line.endswith(';.'):
        return False
    parts = line.split(';')
    if len(parts) < 7:
        return False
    val = parts[6].strip().strip('"')
    if not val or val in ('.', '..'):
        return False
    # Total Nacional origin has parts[0] == 'Total Nacional' and parts[1] == '' and parts[2] == ''
    if parts[0].strip() == 'Total Nacional' and parts[1].strip() == '' and parts[2].strip() == '':
        return True
    per = parts[5].strip()
    if per[:4] in ('2022', '2023', '2024', '2025') and per[4:] in ('M02', 'M08', 'M11'):
        # Keep province-level origin rows (parts[2] non-empty or uniprovincial CCAA)
        return True
    return False

stream_download_table(
    53001,
    keep_53001,
    os.path.join(RAW_DIR, 'turismo_53001_provorig_x_mundest.csv.gz'),
    'Province of origin x Municipality of destination domestic tourists',
)

# Write README.md
readme = f"""# INE Experimental Statistics: Medición del turismo a partir de teléfonos móviles (Turismo interno)

- **Source:** Instituto Nacional de Estadística (INE), *Estadística experimental — Medición del turismo interno interprovincial a partir de la posición de los teléfonos móviles* (3 major Spanish mobile operators: Movistar, Vodafone, Orange).
- **Portal:** `https://www.ine.es/dynt3/inebase/es/index.htm?padre=8578&capsel=8581`
- **Downloaded:** {time.strftime('%Y-%m-%d')}
- **Download script:** `explore/p83_download_ine_turismo.py`
- **Files in `raw/`:**
  - `turismo_52997_prov_x_prov.csv.gz`: Table `52997` (Province of origin × Province of destination, tourists / overnights / avg duration, all non-empty rows).
  - `turismo_53000_munorig_x_provdest.csv.gz`: Table `53000` (Municipality of origin × Province of destination, non-empty rows for all 151 Census cities across all months plus all municipalities for Feb/Mar/Apr/Aug/Nov/Dec).
  - `turismo_53001_provorig_x_mundest.csv.gz`: Table `53001` (Province of origin × Municipality of destination, all non-empty Total Nacional arrivals across all municipalities & months, plus province-origin disaggregation for Feb/Aug/Nov 2022–2025).
- **Important methodological note:** INE defines a domestic tourist as a resident observed sleeping in a **different province** (interprovincial) or at least a threshold distance away for ≥1 night; cells below anonymity disclosure threshold (<15 or <30 tourists) are suppressed (`""`).
- **License:** INE Open Data (free reuse with attribution).
"""
with open(os.path.join(OUT_DIR, 'README.md'), 'w', encoding='utf-8') as f:
    f.write(readme)
print("Finished downloading INE Turismo Moviles tables & README.")

# Check 3: Independent 3-operator data source — INE Experimental Statistics
# "Medición del turismo a partir de teléfonos móviles" (Turismo interno, 2019-2026)
# Replicates:
#   (A) 151-city neighbouring-city roots elasticity (Table 53000: mun of origin x prov of dest)
#       across 2021-2025 (Aug, Nov, Feb, Dec) WITH & WITHOUT ADRH income control & inland-provinces-only!
#   (B) Municipal village-refill / tourist-intensity gradient (Table 53001: Total Nacional domestic tourists
#       arriving in each inland municipality / resident population, vs historical population loss since peak)
#       across 2019-2025 (Aug, Nov, Feb), WITH & WITHOUT ADRH municipal income control!
#   (C) Province x Province tourists & overnights (Table 52997) cross-check against MITMA 2024.
import gzip, os, re, sys, unicodedata, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, 'explore')
import numpy as np, pandas as pd, pyfixest as pf
from common import nodes

COASTAL_PROVS = {'03','04','08','11','12','15','17','18','20','21','27','29','30','33','36','39','43','46','48'}

def norm(s):
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode("ascii")
    return re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()

# Exact 52 INE province name mapping (matching INE JAXI labels)
PROV_MAP = {
    "araba alava": "01", "alava": "01", "araba": "01", "albacete": "02",
    "alicante alacant": "03", "alicante": "03", "almeria": "04", "avila": "05",
    "badajoz": "06", "balears illes": "07", "illes balears": "07", "barcelona": "08",
    "burgos": "09", "caceres": "10", "cadiz": "11", "castellon castello": "12",
    "castellon": "12", "ciudad real": "13", "cordoba": "14", "coruna a": "15", "a coruna": "15",
    "cuenca": "16", "girona": "17", "granada": "18", "guadalajara": "19",
    "gipuzkoa": "20", "guipuzcoa": "20", "huelva": "21", "huesca": "22", "jaen": "23",
    "leon": "24", "lleida": "25", "rioja la": "26", "la rioja": "26", "lugo": "27",
    "madrid": "28", "madrid comunidad de": "28", "comunidad de madrid": "28",
    "malaga": "29", "murcia": "30", "murcia region de": "30", "region de murcia": "30",
    "navarra": "31", "navarra comunidad foral de": "31", "comunidad foral de navarra": "31",
    "ourense": "32", "asturias": "33", "asturias principado de": "33", "principado de asturias": "33",
    "palencia": "34", "palmas las": "35", "las palmas": "35", "pontevedra": "36",
    "salamanca": "37", "santa cruz de tenerife": "38", "cantabria": "39", "segovia": "40",
    "sevilla": "41", "soria": "42", "tarragona": "43", "teruel": "44", "toledo": "45",
    "valencia valencia": "46", "valencia": "46", "valladolid": "47", "bizkaia": "48", "vizcaya": "48",
    "zamora": "49", "zaragoza": "50", "ceuta": "51", "melilla": "52",
}
UNIPROV_CCAA = {
    "asturias principado de": "33", "principado de asturias": "33",
    "balears illes": "07", "illes balears": "07", "cantabria": "39",
    "madrid comunidad de": "28", "comunidad de madrid": "28",
    "murcia region de": "30", "region de murcia": "30",
    "navarra comunidad foral de": "31", "comunidad foral de navarra": "31",
    "rioja la": "26", "la rioja": "26", "ceuta": "51", "melilla": "52",
}

def parse_es_num(s):
    s = s.strip().strip('"')
    if not s or s in ('.', '..', '-'):
        return None
    try:
        return float(s.replace('.', '').replace(',', '.'))
    except ValueError:
        return None

def extract_prov(ccaa_s, prov_s):
    ps = norm(prov_s)
    if ps:
        return PROV_MAP.get(ps)
    cs = norm(ccaa_s)
    if cs in UNIPROV_CCAA:
        return UNIPROV_CCAA[cs]
    return None

# Load shared reference data
B_all = pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv',
                    dtype={'ine_code': str, 'prov_birth': str})
adrh_m = pd.read_csv('_local/external/ine_adrh/adrh_municipios.csv', dtype={'ine_code': str})
N = nodes().set_index('id')
M = N[~N.foreign].copy()
M['w'] = M.population.fillna(0) + 1
pc = M.groupby('prov').apply(lambda g: pd.Series({'x': np.average(g.x, weights=g.w), 'y': np.average(g.y, weights=g.w)}))
cities = sorted(set(B_all.ine_code) & set(N.index))
city_set = set(cities)

# ==============================================================================
# PART A: TABLE 53000 — 151 CITIES OF ORIGIN x DESTINATION PROVINCE (2019-2025)
# ==============================================================================
print("=" * 110)
print("PART A: INE TURISMO MOVILES TABLE 53000 — 151 CITIES x DEST PROVINCE ROOTS ELASTICITY")
print("=" * 110)

p53000 = '_local/external/ine_turismo_moviles/raw/turismo_53000_munorig_x_provdest.csv.gz'
recs = []
unmatched_p = set()
with gzip.open(p53000, 'rt', encoding='utf-8') as f:
    hdr = f.readline()
    for line in f:
        parts = line.rstrip('\n').split(';')
        if len(parts) < 7:
            continue
        mun5 = parts[0][:5]
        if mun5 not in city_set:
            continue
        # parts[1] = Total Nacional, parts[2] = CCAA, parts[3] = Prov
        pcode = extract_prov(parts[2], parts[3])
        if pcode is None:
            if parts[3].strip():
                unmatched_p.add(parts[3].strip())
            continue
        per = parts[5].strip()
        if not re.match(r'^20\d\dM(02|08|11|12)$', per):
            continue
        val = parse_es_num(parts[6])
        if val is None:
            continue
        recs.append((mun5, pcode, per, val))

print(f"Parsed {len(recs):,} (city, dest_prov, period) rows for 151 cities; unmatched prov labels: {unmatched_p}")
df53000 = pd.DataFrame(recs, columns=['c', 'p', 'per', 'tourists'])
# Deduplicate if any (e.g. uniprovincial CCAA also has prov row — take max/first; verify!)
dup_cnt = df53000.duplicated(['c', 'p', 'per']).sum()
if dup_cnt:
    print(f"  Note: {dup_cnt} duplicate uniprovincial rows collapsed (max).")
    df53000 = df53000.groupby(['c', 'p', 'per'], as_index=False).tourists.max()

print("Unique cities in 53000:", df53000.c.nunique(), "Unique dest provs:", df53000.p.nunique(),
      "Periods:", sorted(df53000.per.unique())[:8], "...", sorted(df53000.per.unique())[-8:])

# Build base grid
base_grid = pd.MultiIndex.from_product([cities, [f'{i:02d}' for i in range(1, 51)]], names=['c', 'p']).to_frame(index=False)
base_grid['cprov'] = base_grid.c.str[:2]
base_grid = base_grid[base_grid.cprov != base_grid.p].copy()
base_grid['km'] = np.hypot(N.loc[base_grid.c, 'x'].values - pc.loc[base_grid.p, 'x'].values,
                           N.loc[base_grid.c, 'y'].values - pc.loc[base_grid.p, 'y'].values) / 1000.0
base_grid['dbin'] = np.digitize(base_grid.km, np.r_[0:100:10, 100:400:20, 400:1200:50])
base_grid['pairFE'] = base_grid.cprov + '_' + base_grid.p

# Run for 2024 (Feb, Aug, Nov, Dec) + all Augusts 2019-2025 & Novembers 2019-2025!
inc24 = adrh_m[(adrh_m.year == 2024) & (adrh_m.indicator == 'Renta neta media por persona')].set_index('ine_code').value.reindex(cities)
linc24 = np.log(inc24) - np.log(inc24).mean()

target_pers = ['2019M08', '2019M11', '2021M08', '2022M08', '2022M11', '2023M08', '2023M11',
               '2024M02', '2024M08', '2024M11', '2024M12', '2025M08', '2025M11']
rows_a = []
for per in target_pers:
    sub_f = df53000[df53000.per == per]
    if len(sub_f) == 0:
        continue
    yr = int(per[:4])
    byr = min(max(yr, 2021), 2025)
    By = B_all[(B_all.year == byr) & B_all.prov_birth.str.match(r'^\d\d$')].set_index(['ine_code', 'prov_birth']).persons
    inc_yr = min(max(yr, 2019), 2024)
    inc_y = adrh_m[(adrh_m.year == inc_yr) & (adrh_m.indicator == 'Renta neta media por persona')].set_index('ine_code').value.reindex(cities)
    linc_y = np.log(inc_y) - np.log(inc_y).mean()

    f_s = sub_f.set_index(['c', 'p']).tourists
    g = base_grid.copy()
    g['w'] = f_s.reindex(list(zip(g.c, g.p))).fillna(0).values
    g['lborn'] = np.log1p(By.reindex(list(zip(g.c, g.p))).fillna(0).values)
    g['linc'] = g.c.map(linc_y).values

    # All 50 dest provs & Inland-only 28 dest provs
    t_d = pf.fepois("w ~ lborn | c + p + dbin", data=g, vcov={'CRV1': 'c'}).tidy()
    t_p = pf.fepois("w ~ lborn | c + pairFE", data=g, vcov={'CRV1': 'c'}).tidy()
    t_pi = pf.fepois("w ~ lborn | c + pairFE + p[linc]", data=g, vcov={'CRV1': 'c'}).tidy()
    g_inl = g[~g.p.isin(COASTAL_PROVS)].copy()
    t_inl = pf.fepois("w ~ lborn | c + pairFE + p[linc]", data=g_inl, vcov={'CRV1': 'c'}).tidy()
    rows_a.append(dict(
        period=per, nz_cells=int((g.w > 0).sum()), total_M=g.w.sum()/1e6,
        distFE=t_d.loc['lborn','Estimate'],
        pairFE=t_p.loc['lborn','Estimate'], pair_se=t_p.loc['lborn','Std. Error'],
        pair_inc=t_pi.loc['lborn','Estimate'], pair_inc_se=t_pi.loc['lborn','Std. Error'],
        inland_pair_inc=t_inl.loc['lborn','Estimate'], inland_se=t_inl.loc['lborn','Std. Error'],
    ))

df_ra = pd.DataFrame(rows_a)
df_ra.to_csv('explore/out_p89_ine_turismo_roots.csv', index=False)
pd.set_option('display.width', 240)
print(df_ra.round(3).to_string(index=False))

# Check Barcelona & Madrid within-metro correlations in INE Turismo 2024M08!
AND_ALL = ['04','11','14','18','21','23','29','41']
EXT = ['06','10']
B24 = B_all[(B_all.year == 2024) & B_all.prov_birth.str.match(r'^\d\d$')]
tot_b = B24.groupby('ine_code').persons.sum()
born_and = 100 * B24[B24.prov_birth.isin(AND_ALL)].groupby('ine_code').persons.sum() / tot_b
born_ext = 100 * B24[B24.prov_birth.isin(EXT)].groupby('ine_code').persons.sum() / tot_b
f24a = df53000[(df53000.per == '2024M08') & (df53000.c.str[:2] != df53000.p)]
tot_t = f24a.groupby('c').tourists.sum()
sh_and = 100 * f24a[f24a.p.isin(AND_ALL)].groupby('c').tourists.sum() / tot_t
sh_ext = 100 * f24a[f24a.p.isin(EXT)].groupby('c').tourists.sum() / tot_t
bcn_c = [c for c in cities if c.startswith('08')]
mad_c = [c for c in cities if c.startswith('28')]
print(f"\nWithin-metro correlations in INE Turismo 2024M08 (3 operators):")
print(f"  Barcelona (n={len(bcn_c)}): corr(% born AND, % Aug tourists to AND) = {born_and.reindex(bcn_c).corr(sh_and.reindex(bcn_c)):+.3f}")
print(f"  Barcelona (n={len(bcn_c)}): corr(% born EXT, % Aug tourists to EXT) = {born_ext.reindex(bcn_c).corr(sh_ext.reindex(bcn_c)):+.3f}")
print(f"  Madrid    (n={len(mad_c)}): corr(% born AND, % Aug tourists to AND) = {born_and.reindex(mad_c).corr(sh_and.reindex(mad_c)):+.3f}")
print(f"  Madrid    (n={len(mad_c)}): corr(% born EXT, % Aug tourists to EXT) = {born_ext.reindex(mad_c).corr(sh_ext.reindex(mad_c)):+.3f}")

# ==============================================================================
# PART B: TABLE 53001 — MUNICIPAL VILLAGE-REFILL GRADIENT (ALL ~8,131 MUNIS!)
# ==============================================================================
print("\n" + "=" * 110)
print("PART B: INE TURISMO MOVILES TABLE 53001 — MUNICIPAL TOURIST ARRIVALS / RESIDENTS vs EXODUS LOSS")
print("=" * 110)
p53001 = '_local/external/ine_turismo_moviles/raw/turismo_53001_provorig_x_mundest.csv.gz'
m_recs = []
with gzip.open(p53001, 'rt', encoding='utf-8') as f:
    hdr = f.readline()
    for line in f:
        parts = line.rstrip('\n').split(';')
        if len(parts) < 7:
            continue
        # Total Nacional origin only!
        if parts[0].strip() != 'Total Nacional' or parts[1].strip() != '' or parts[2].strip() != '':
            continue
        m_ine = parts[3][:5]
        if not re.match(r'^\d{5}$', m_ine):
            continue
        per = parts[5].strip()
        if not re.match(r'^20\d\dM(02|08|11)$', per):
            continue
        val = parse_es_num(parts[6])
        if val is None:
            continue
        m_recs.append((m_ine, per, val))

df53001 = pd.DataFrame(m_recs, columns=['ine_code', 'per', 'tourists'])
print(f"Parsed {len(df53001):,} (municipality, period) Total Nacional arrival rows across {df53001.ine_code.nunique():,} unique municipalities.")

# Load historical pop (8,132 munis), coastal flag (mapped via relacion_ine_zonificacionMitma), and ADRH 2024 income
H = pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv',
                dtype={'ine_code': str, 'province_code': str})
H['peak'] = H[['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960']].max(axis=1)
H['loss_peak_2021'] = 1 - H.pop_2021 / H.peak
H['pop'] = H.padron_2024.fillna(H.pop_2021)
rel = pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv', sep='|', dtype=str)[
    ['municipio_ine', 'municipio_mitma']
].drop_duplicates()
cf = pd.read_csv('explore/out_p86_coastal_flags.csv', dtype={'id': str}).rename(columns={'id': 'municipio_mitma'})
H = H.merge(rel, left_on='ine_code', right_on='municipio_ine', how='left').merge(cf, on='municipio_mitma', how='left')
H['inc24'] = H.ine_code.map(adrh_m[(adrh_m.year == 2024) & (adrh_m.indicator == 'Renta neta media por persona')].set_index('ine_code').value)
# Strictly mainland inland municipalities with pop >= 100 and non-missing historical loss
H_inl = H[(~H.province_code.isin(['07','35','38','51','52'])) & (H.coastal_true == 0) & (H['pop'] >= 100) & H.loss_peak_2021.notna()].copy()
H_inl['lpop'] = np.log(H_inl['pop'])
H_inl['linc'] = np.log(H_inl['inc24'])
H_inl['loss_rank'] = H_inl.loss_peak_2021.rank(pct=True)
H_inl['loss_q'] = pd.qcut(H_inl.loss_peak_2021, 5, labels=['Q1 grew','Q2','Q3','Q4','Q5 emptied'])

print(f"Mainland inland municipalities (pop>=100): {len(H_inl):,}")
rows_b = []
for per in ['2019M08', '2022M08', '2023M08', '2024M02', '2024M08', '2024M11', '2025M08', '2025M11']:
    t_sub = df53001[df53001.per == per].set_index('ine_code').tourists
    D = H_inl.copy()
    D['tourists'] = D.ine_code.map(t_sub)
    n_obs = int(D.tourists.notna().sum())
    if n_obs < 1000:
        continue
    # Tourists per resident (observed non-suppressed; also check filling suppressed with 0!)
    D_obs = D[D.tourists.notna() & (D.tourists > 0)].copy()
    D_obs['ratio'] = D_obs.tourists / D_obs['pop']
    D_obs['lratio'] = np.log(D_obs['ratio'])
    D_all = D.copy()
    D_all['ratio0'] = D_all.tourists.fillna(0) / D_all['pop']
    qt = D_obs.groupby('loss_q', observed=True)['ratio'].median()
    qt0 = D_all.groupby('loss_q', observed=True)['ratio0'].median()
    m0 = pf.feols("lratio ~ loss_rank + lpop | province_code", data=D_obs, vcov={'CRV1': 'province_code'}).tidy()
    m1 = pf.feols("lratio ~ loss_rank + lpop + linc | province_code", data=D_obs.dropna(subset=['linc']), vcov={'CRV1': 'province_code'}).tidy()
    rows_b.append(dict(
        period=per, n_obs=n_obs, cov_pct=100*n_obs/len(D),
        q1_ratio=qt.loc['Q1 grew'], q5_ratio=qt.loc['Q5 emptied'],
        q1_r0=qt0.loc['Q1 grew'], q5_r0=qt0.loc['Q5 emptied'],
        q5_over_q1=qt.loc['Q5 emptied'] / qt.loc['Q1 grew'],
        ols_loss=m0.loc['loss_rank','Estimate'], se=m0.loc['loss_rank','Std. Error'],
        ols_loss_inc=m1.loc['loss_rank','Estimate'], se_inc=m1.loc['loss_rank','Std. Error'],
    ))

df_rb = pd.DataFrame(rows_b)
df_rb.to_csv('explore/out_p89_ine_turismo_refill.csv', index=False)
print(df_rb.round(3).to_string(index=False))

# Size-matched check in 2024M08 (terciles of population among inland municipalities)
D24 = H_inl.copy()
D24['tourists'] = D24.ine_code.map(df53001[df53001.per == '2024M08'].set_index('ine_code').tourists)
D24['ratio0'] = D24.tourists.fillna(0) / D24['pop']
D24['size_t'] = pd.qcut(D24['pop'], 3, labels=['small (100-350)', 'mid (350-1.5k)', 'large (>1.5k)'])
print("\n=== INE Turismo 2024M08 Domestic Tourists / Resident by Exodus Quintile x Size Tercile (ALL inland munis, suppressed=0) ===")
print(D24.pivot_table(index='loss_q', columns='size_t', values='ratio0', aggfunc='median', observed=True).round(3).to_string())

# Bonus: within-municipality seasonal ratio (2024M08 / 2024M02) to remove any permanent tourism baseline!
t_aug24 = df53001[df53001.per == '2024M08'].set_index('ine_code').tourists
t_nov24 = df53001[df53001.per == '2024M11'].set_index('ine_code').tourists
t_feb24 = df53001[df53001.per == '2024M02'].set_index('ine_code').tourists
Ds = H_inl.copy()
Ds['aug'] = Ds.ine_code.map(t_aug24)
Ds['nov'] = Ds.ine_code.map(t_nov24)
Ds['feb'] = Ds.ine_code.map(t_feb24)
Ds_af = Ds.dropna(subset=['aug', 'feb'])
Ds_af = Ds_af[(Ds_af.aug > 0) & (Ds_af.feb > 0)].copy()
Ds_af['aug_over_feb'] = Ds_af.aug / Ds_af.feb
Ds_af['l_af'] = np.log(Ds_af.aug_over_feb)
qt_af = Ds_af.groupby('loss_q', observed=True).agg(
    n=('pop','size'), med_loss=('loss_peak_2021','median'),
    aug_per_res=('aug', lambda s: (s / Ds_af.loc[s.index, 'pop']).median()),
    feb_per_res=('feb', lambda s: (s / Ds_af.loc[s.index, 'pop']).median()),
    aug_over_feb=('aug_over_feb', 'median'),
)
m_af = pf.feols("l_af ~ loss_rank + lpop + linc | province_code", data=Ds_af.dropna(subset=['linc']), vcov={'CRV1': 'province_code'}).tidy()
print(f"\n=== INE Turismo 2024 Within-Municipality Seasonal Surge (Aug / Feb, n={len(Ds_af):,} inland munis) ===")
print(qt_af.round(3).to_string())
print(f"OLS log(Aug/Feb) ~ loss_rank + lpop + linc | prov: loss_rank={m_af.loc['loss_rank','Estimate']:+.3f} (se {m_af.loc['loss_rank','Std. Error']:.3f})")

# Same for Nov / Feb!
Ds_nf = Ds.dropna(subset=['nov', 'feb'])
Ds_nf = Ds_nf[(Ds_nf.nov > 0) & (Ds_nf.feb > 0)].copy()
Ds_nf['nov_over_feb'] = Ds_nf.nov / Ds_nf.feb
Ds_nf['l_nf'] = np.log(Ds_nf.nov_over_feb)
m_nf = pf.feols("l_nf ~ loss_rank + lpop + linc | province_code", data=Ds_nf.dropna(subset=['linc']), vcov={'CRV1': 'province_code'}).tidy()
print(f"OLS log(Nov/Feb) ~ loss_rank + lpop + linc | prov (n={len(Ds_nf):,}): loss_rank={m_nf.loc['loss_rank','Estimate']:+.3f} (se {m_nf.loc['loss_rank','Std. Error']:.3f})")

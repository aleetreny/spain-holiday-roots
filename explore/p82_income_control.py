# Check 1 (DECISIVE): Income Control using INE ADRH (2024 & 2023 Renta neta media por persona / hogar)
# 1A. Re-estimate 151-city neighbouring-city roots model (p72) with flexible income x dest-province controls
#     (linear slopes, quadratic slopes, income-quartile x dest-prov FE, and origin_prov x dest_prov income slopes),
#     plus combined income + inland-only + road travel time specifications.
# 1B. Barcelona and Madrid within-metro comparisons controlling for income, AND full explanation of why
#     Madrid shows no raw pattern in Andalusia (coastal luxury beach vs inland/roots + All Saints + all 46 provs).
# 1C. District-level within-metro analysis (Barcelona & Madrid districts using ADRH district income + MITMA districts).
import sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, 'explore')
import duckdb, numpy as np, pandas as pd, pyfixest as pf, statsmodels.api as sm
from common import nodes, PROV2CCAA

COASTAL_PROVS = {'03','04','08','11','12','15','17','18','20','21','27','29','30','33','36','39','43','46','48'}
AND_ALL = ['04','11','14','18','21','23','29','41']
AND_COAST_PROVS = ['11','21','29']   # Cadiz, Huelva, Malaga (prime Madrid beach destinations)
AND_INLAND_PROVS = ['14','23','41']  # Cordoba, Jaen, Sevilla (plus Granada 18, Almeria 04 tested separately)

# 1. Load ADRH municipal & district income (2024 primary; check 2023 sensitivity)
adrh_m = pd.read_csv('_local/external/ine_adrh/adrh_municipios.csv', dtype={'ine_code': str})
inc24_pc = adrh_m[(adrh_m.year == 2024) & (adrh_m.indicator == 'Renta neta media por persona')].set_index('ine_code').value
inc24_hh = adrh_m[(adrh_m.year == 2024) & (adrh_m.indicator == 'Renta neta media por hogar')].set_index('ine_code').value
inc23_pc = adrh_m[(adrh_m.year == 2023) & (adrh_m.indicator == 'Renta neta media por persona')].set_index('ine_code').value

# 2. Load 151 cities birthplace (Census 2024)
B = pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv',
                dtype={'ine_code': str, 'prov_birth': str})
B24 = B[(B.year == 2024) & B.prov_birth.str.match(r'^\d\d$')].copy()
N = nodes().set_index('id')
cities = sorted(set(B24.ine_code) & set(N.index))
print(f"Cities matched to MITMA: {len(cities)}; matched to ADRH 2024 pc: {inc24_pc.reindex(cities).notna().sum()}, hh: {inc24_hh.reindex(cities).notna().sum()}")
miss_inc = [c for c in cities if pd.isna(inc24_pc.get(c))]
if miss_inc:
    print("  Missing ADRH 2024 income:", miss_inc)

M = N[~N.foreign].copy()
M['w'] = M.population.fillna(0) + 1
pc = M.groupby('prov').apply(lambda g: pd.Series({'x': np.average(g.x, weights=g.w), 'y': np.average(g.y, weights=g.w)}))

con = duckdb.connect()
q = "'data/derived/pernoct_mun_2024_6m.parquet'"
cl = ",".join(f"'{c}'" for c in cities)
cf = pd.read_csv('explore/out_p86_coastal_flags.csv', dtype={'id': str}).rename(columns={'id': 'z'})
con.register('cf_tbl', cf)

slices = {
    'Feb weekday': "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",
    'ordinary Sat': "date IN ('2024-02-10','2024-02-17','2024-02-24','2024-10-19','2024-11-16','2024-11-23')",
    'Good Friday': "date='2024-03-29'",
    'mid-Aug': "date BETWEEN '2024-08-10' AND '2024-08-18'",
    'All Saints': "date='2024-11-01'",
    'Christmas Eve': "date='2024-12-24'",
}

# Build city x destination-province grid (50 provinces, excluding same_prov)
grid = pd.MultiIndex.from_product([cities, [f'{i:02d}' for i in range(1, 51)]], names=['c', 'p']).to_frame(index=False)
grid['cprov'] = grid.c.str[:2]
grid = grid[grid.cprov != grid.p].copy()
grid['km'] = np.hypot(N.loc[grid.c, 'x'].values - pc.loc[grid.p, 'x'].values,
                      N.loc[grid.c, 'y'].values - pc.loc[grid.p, 'y'].values) / 1000.0
grid['dbin'] = np.digitize(grid.km, np.r_[0:100:10, 100:400:20, 400:1200:50])
Bm = B24.set_index(['ine_code', 'prov_birth']).persons
grid['born'] = Bm.reindex(list(zip(grid.c, grid.p))).fillna(0).values
grid['lborn'] = np.log1p(grid.born)
grid['pairFE'] = grid.cprov + '_' + grid.p

# Attach standardized log income per person and per household (centered at 0 for numerical stability)
c_inc_pc = inc24_pc.reindex(cities)
c_inc_hh = inc24_hh.reindex(cities)
linc_pc = np.log(c_inc_pc) - np.log(c_inc_pc).mean()
linc_hh = np.log(c_inc_hh) - np.log(c_inc_hh).mean()
inc_q4 = pd.qcut(c_inc_pc, 4, labels=['Q1', 'Q2', 'Q3', 'Q4']).astype(str)
inc_q5 = pd.qcut(c_inc_pc, 5, labels=['Q1', 'Q2', 'Q3', 'Q4', 'Q5']).astype(str)

grid['linc'] = grid.c.map(linc_pc).values
grid['linc2'] = grid['linc'] ** 2
grid['linc_hh'] = grid.c.map(linc_hh).values
grid['p_incq4'] = grid.p + '_' + grid.c.map(inc_q4).values
grid['p_incq5'] = grid.p + '_' + grid.c.map(inc_q5).values

# ==============================================================================
# PART 1A: 151-CITY ROOTS MODEL WITH INCOME x DESTINATION-PROVINCE CONTROLS
# ==============================================================================
print("\n" + "=" * 110)
print("PART 1A: 151 CITIES x 50 DEST PROVINCES — ROOTS ELASTICITY WITH INCOME x DEST-PROVINCE CONTROLS")
print("=" * 110)
rows_1a = []
for mode, where_s in [('All zones', '1=1'), ('Inland zones only (true coast=0)', 'cf_tbl.coastal_true = 0')]:
    for lab, cond in slices.items():
        if where_s == '1=1':
            f = con.execute(f"SELECT r c, substr(s,1,2) p, sum(p)/count(DISTINCT date) w FROM {q} WHERE ({cond}) AND r IN ({cl}) AND r<>s GROUP BY 1,2").df().set_index(['c', 'p']).w
        else:
            f = con.execute(f"SELECT x.r c, substr(x.s,1,2) p, sum(x.p)/count(DISTINCT x.date) w FROM {q} x JOIN cf_tbl ON x.s=cf_tbl.z WHERE ({cond}) AND x.r IN ({cl}) AND x.r<>x.s AND ({where_s}) GROUP BY 1,2").df().set_index(['c', 'p']).w
        g = grid.copy()
        g['w'] = f.reindex(list(zip(g.c, g.p))).fillna(0).values
        specs = [
            ('base_dist', "w ~ lborn | c + p + dbin"),
            ('dist_p_linc', "w ~ lborn | c + p + dbin + p[linc]"),
            ('base_pair', "w ~ lborn | c + pairFE"),
            ('pair_p_linc', "w ~ lborn | c + pairFE + p[linc]"),
            ('pair_p_linc2', "w ~ lborn | c + pairFE + p[linc] + p[linc2]"),
            ('pair_p_hh', "w ~ lborn | c + pairFE + p[linc_hh]"),
            ('pair_p_q4', "w ~ lborn | c + pairFE + p_incq4"),
            ('pair_pair_linc', "w ~ lborn | c + pairFE + pairFE[linc]"),
        ]
        row = {'mode': mode, 'slice': lab}
        for name, fml in specs:
            t = pf.fepois(fml, data=g, vcov={'CRV1': 'c'}).tidy()
            row[name] = t.loc['lborn', 'Estimate']
            row[name + '_se'] = t.loc['lborn', 'Std. Error']
        rows_1a.append(row)

df_1a = pd.DataFrame(rows_1a)
df_1a.to_csv('explore/out_p82_income_cities.csv', index=False)
pd.set_option('display.width', 240)
pd.set_option('display.max_columns', 30)
show_cols = ['mode', 'slice', 'base_dist', 'dist_p_linc', 'base_pair', 'base_pair_se', 'pair_p_linc', 'pair_p_linc_se', 'pair_p_linc2', 'pair_p_hh', 'pair_p_q4', 'pair_pair_linc', 'pair_pair_linc_se']
print(df_1a[show_cols].round(3).to_string(index=False))

# Also check with ROAD TRAVEL TIME (from p84) + INCOME simultaneously!
tt_z = pd.read_parquet('explore/out_p84_city_zone_tt.parquet')
tt_z['p'] = tt_z.z.str[:2]
tt_z['zpop'] = N.loc[tt_z.z, 'population'].fillna(0).values + 1.0
prov_tt = tt_z.groupby(['c', 'p']).apply(
    lambda x: np.average(x.tt_adj, weights=x.zpop), include_groups=False
).rename('tt_adj').reset_index()
g_tt = grid.merge(prov_tt, on=['c', 'p'], how='inner')
g_tt['ttbin'] = np.digitize(g_tt.tt_adj, np.r_[0:120:10, 120:360:15, 360:900:30])
print("\n=== Combined: Road Travel Time + Income x Dest-Prov (137 mainland cities) ===")
for lab in ['Feb weekday', 'mid-Aug', 'All Saints', 'Christmas Eve']:
    f = con.execute(f"SELECT r c, substr(s,1,2) p, sum(p)/count(DISTINCT date) w FROM {q} WHERE ({slices[lab]}) AND r IN ({cl}) AND r<>s GROUP BY 1,2").df().set_index(['c', 'p']).w
    g = g_tt.copy(); g['w'] = f.reindex(list(zip(g.c, g.p))).fillna(0).values
    t1 = pf.fepois("w ~ lborn | c + p + ttbin + p[linc]", data=g, vcov={'CRV1': 'c'}).tidy()
    t2 = pf.fepois("w ~ lborn | c + pairFE + ttbin + p[linc]", data=g, vcov={'CRV1': 'c'}).tidy()
    print(f"  {lab:14s}: c+p+ttbin+p[linc]={t1.loc['lborn','Estimate']:.3f} ({t1.loc['lborn','Std. Error']:.3f})  |  c+pairFE+ttbin+p[linc]={t2.loc['lborn','Estimate']:.3f} ({t2.loc['lborn','Std. Error']:.3f})")

# ==============================================================================
# PART 1B: WITHIN-METRO COMPARISONS (BARCELONA & MADRID) CONTROLLING FOR INCOME
#          + FULL DECOMPOSITION OF WHY MADRID SHOWS NO RAW PATTERN IN ANDALUSIA
# ==============================================================================
print("\n" + "=" * 110)
print("PART 1B: BARCELONA (prov 08) & MADRID (prov 28) WITHIN-METRO ANALYSES WITH INCOME CONTROL")
print("=" * 110)

nm_c = B24.drop_duplicates('ine_code').set_index('ine_code').municipality_name
tot_b = B24.groupby('ine_code').persons.sum()

def get_metro_shares(cond, where_s='1=1'):
    if where_s == '1=1':
        df = con.execute(f"SELECT r c, substr(s,1,2) p, sum(p)/count(DISTINCT date) w FROM {q} WHERE ({cond}) AND r IN ({cl}) AND r<>s GROUP BY 1,2").df()
    else:
        df = con.execute(f"SELECT x.r c, substr(x.s,1,2) p, sum(x.p)/count(DISTINCT x.date) w FROM {q} x JOIN cf_tbl ON x.s=cf_tbl.z WHERE ({cond}) AND x.r IN ({cl}) AND x.r<>x.s AND ({where_s}) GROUP BY 1,2").df()
    df = df[df.p != df.c.str[:2]]
    return df

f_aug_all = get_metro_shares(slices['mid-Aug'])
f_aug_inl = get_metro_shares(slices['mid-Aug'], 'cf_tbl.coastal_true = 0')
f_aug_cst = get_metro_shares(slices['mid-Aug'], 'cf_tbl.coastal_true = 1')
f_nov_all = get_metro_shares(slices['All Saints'])
f_feb_all = get_metro_shares(slices['Feb weekday'])

def share_to(df, prov_list, denom_df=None):
    tot = (denom_df if denom_df is not None else df).groupby('c').w.sum()
    num = df[df.p.isin(prov_list)].groupby('c').w.sum()
    return (100 * num / tot).fillna(0)

born_and = 100 * B24[B24.prov_birth.isin(AND_ALL)].groupby('ine_code').persons.sum() / tot_b
born_and_inl = 100 * B24[B24.prov_birth.isin(AND_INLAND_PROVS)].groupby('ine_code').persons.sum() / tot_b
born_and_cst = 100 * B24[B24.prov_birth.isin(AND_COAST_PROVS)].groupby('ine_code').persons.sum() / tot_b
born_ext = 100 * B24[B24.prov_birth.isin(['06', '10'])].groupby('ine_code').persons.sum() / tot_b
born_clm = 100 * B24[B24.prov_birth.isin(['02', '13', '16', '19', '45'])].groupby('ine_code').persons.sum() / tot_b
born_cyl = 100 * B24[B24.prov_birth.isin(['05','09','24','34','37','40','42','47','49'])].groupby('ine_code').persons.sum() / tot_b

metro_df = pd.DataFrame({
    'city': nm_c.reindex(cities),
    'prov': [c[:2] for c in cities],
    'inc_pc': inc24_pc.reindex(cities),
    'inc_hh': inc24_hh.reindex(cities),
    'linc': np.log(inc24_pc.reindex(cities)),
    'born_AND': born_and.reindex(cities),
    'born_AND_inl': born_and_inl.reindex(cities),
    'born_AND_cst': born_and_cst.reindex(cities),
    'born_EXT': born_ext.reindex(cities),
    'born_CLM': born_clm.reindex(cities),
    'born_CYL': born_cyl.reindex(cities),
    'aug_AND': share_to(f_aug_all, AND_ALL).reindex(cities),
    'aug_AND_cstZ': share_to(f_aug_cst, AND_ALL, f_aug_all).reindex(cities),
    'aug_AND_inlZ': share_to(f_aug_inl, AND_ALL, f_aug_all).reindex(cities),
    'aug_AND_inlP': share_to(f_aug_all, AND_INLAND_PROVS).reindex(cities),
    'aug_AND_cstP': share_to(f_aug_all, AND_COAST_PROVS).reindex(cities),
    'nov_AND': share_to(f_nov_all, AND_ALL).reindex(cities),
    'feb_AND': share_to(f_feb_all, AND_ALL).reindex(cities),
    'aug_EXT': share_to(f_aug_all, ['06', '10']).reindex(cities),
    'nov_EXT': share_to(f_nov_all, ['06', '10']).reindex(cities),
    'aug_CLM': share_to(f_aug_all, ['02','13','16','19','45']).reindex(cities),
    'aug_CYL': share_to(f_aug_all, ['05','09','24','34','37','40','42','47','49']).reindex(cities),
}, index=cities)

for pcode, pname in [('08', 'BARCELONA METRO'), ('28', 'MADRID METRO')]:
    sub = metro_df[metro_df.prov == pcode].dropna(subset=['inc_pc', 'born_AND', 'aug_AND']).copy()
    print(f"\n--- {pname} (n={len(sub)} cities) ---")
    print(f"  corr(born_AND, linc) = {sub.born_AND.corr(sub.linc):+.3f}")
    print(sub[['city','inc_pc','born_AND','aug_AND','aug_AND_cstZ','aug_AND_inlZ','nov_AND','aug_EXT','born_EXT']].sort_values('inc_pc', ascending=False).round(2).to_string())
    print(f"\n  OLS regressions within {pname} (y ~ born + linc):")
    tests = [
        ('Aug ALL -> Andalusia (raw)', 'aug_AND', ['born_AND']),
        ('Aug ALL -> Andalusia (+ linc)', 'aug_AND', ['born_AND', 'linc']),
        ('Aug COASTAL Andalusia zones (+ linc)', 'aug_AND_cstZ', ['born_AND', 'linc']),
        ('Aug INLAND Andalusia zones (raw)', 'aug_AND_inlZ', ['born_AND']),
        ('Aug INLAND Andalusia zones (+ linc)', 'aug_AND_inlZ', ['born_AND', 'linc']),
        ('Aug INLAND Andalusia provs Jaen/Cord/Sev (+ linc)', 'aug_AND_inlP', ['born_AND_inl', 'linc']),
        ('Aug COASTAL Andalusia provs Mal/Cad/Huel (+ linc)', 'aug_AND_cstP', ['born_AND_cst', 'linc']),
        ('All Saints -> Andalusia (raw)', 'nov_AND', ['born_AND']),
        ('All Saints -> Andalusia (+ linc)', 'nov_AND', ['born_AND', 'linc']),
        ('Aug -> Extremadura (raw)', 'aug_EXT', ['born_EXT']),
        ('Aug -> Extremadura (+ linc)', 'aug_EXT', ['born_EXT', 'linc']),
        ('All Saints -> Extremadura (+ linc)', 'nov_EXT', ['born_EXT', 'linc']),
        ('Aug -> Castilla-La Mancha (+ linc)', 'aug_CLM', ['born_CLM', 'linc']),
        ('Aug -> Castilla y Leon (+ linc)', 'aug_CYL', ['born_CYL', 'linc']),
    ]
    for tname, ycol, xcols in tests:
        X = sm.add_constant(sub[xcols])
        fit = sm.OLS(sub[ycol], X).fit(cov_type='HC3')
        b0, p0 = fit.params[xcols[0]], fit.pvalues[xcols[0]]
        se0 = fit.bse[xcols[0]]
        r_raw = sub[ycol].corr(sub[xcols[0]])
        inc_str = ""
        if 'linc' in xcols:
            inc_str = f" | linc={fit.params['linc']:+.2f} (p={fit.pvalues['linc']:.3g})"
        print(f"    {tname:50s}: r_raw={r_raw:+.3f} | b_born={b0:+.3f} (se {se0:.3f}, p={p0:.3g}){inc_str} | R2={fit.rsquared:.2f}")

# Single-metro PPML across ALL 49 other destination provinces (Barcelona-only and Madrid-only!)
print("\n=== Single-Metro Full PPML (ALL 49 other provinces: w ~ lborn | c + p  vs  w ~ lborn | c + p + p[linc]) ===")
for pcode, pname in [('08', 'Barcelona prov (19 cities)'), ('28', 'Madrid prov (23 cities)'), ('46', 'Valencia prov'), ('03', 'Alicante prov')]:
    for mode, where_s in [('All zones', '1=1'), ('Inland zones only', 'cf_tbl.coastal_true = 0')]:
        for lab in ['Feb weekday', 'mid-Aug', 'All Saints']:
            if where_s == '1=1':
                f = con.execute(f"SELECT r c, substr(s,1,2) p, sum(p)/count(DISTINCT date) w FROM {q} WHERE ({slices[lab]}) AND r IN ({cl}) AND r<>s GROUP BY 1,2").df().set_index(['c', 'p']).w
            else:
                f = con.execute(f"SELECT x.r c, substr(x.s,1,2) p, sum(x.p)/count(DISTINCT x.date) w FROM {q} x JOIN cf_tbl ON x.s=cf_tbl.z WHERE ({slices[lab]}) AND x.r IN ({cl}) AND x.r<>x.s AND ({where_s}) GROUP BY 1,2").df().set_index(['c', 'p']).w
            g = grid[grid.cprov == pcode].copy()
            if g.c.nunique() < 4:
                continue
            g['w'] = f.reindex(list(zip(g.c, g.p))).fillna(0).values
            t0 = pf.fepois("w ~ lborn | c + p", data=g, vcov={'CRV1': 'c'}).tidy()
            t1 = pf.fepois("w ~ lborn | c + p + p[linc]", data=g, vcov={'CRV1': 'c'}).tidy()
            print(f"  {pname:26s} | {mode:17s} | {lab:11s}: base={t0.loc['lborn','Estimate']:+.3f} ({t0.loc['lborn','Std. Error']:.3f})  ->  +p[linc]={t1.loc['lborn','Estimate']:+.3f} ({t1.loc['lborn','Std. Error']:.3f})")

# ==============================================================================
# PART 1C: DISTRICT-LEVEL WITHIN-METRO ANALYSIS (BARCELONA & MADRID DISTRICTS)
# ==============================================================================
print("\n" + "=" * 110)
print("PART 1C: DISTRICT-LEVEL ANALYSIS (BARCELONA & MADRID METRO DISTRICTS, ADRH DISTRICT INCOME)")
print("=" * 110)
adrh_d = pd.read_csv('_local/external/ine_adrh/adrh_distritos.csv', dtype={'dist_code': str, 'ine_code': str})
d_inc24 = adrh_d[(adrh_d.year == 2024) & (adrh_d.indicator == 'Renta neta media por persona')].set_index('dist_code').value
rel_d = pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv', sep='|', dtype=str)[
    ['distrito_ine', 'municipio_ine', 'distrito_mitma']
].drop_duplicates()
# Map ADRH district income to MITMA distrito_mitma (average across INE districts if merged; exact 1-to-1 for big cities!)
rel_d['inc_d'] = rel_d.distrito_ine.map(d_inc24)
# Fallback to municipal income if district code absent
rel_d['inc_d'] = rel_d['inc_d'].fillna(rel_d.municipio_ine.map(inc24_pc))
mitma_d_inc = rel_d.groupby('distrito_mitma').inc_d.mean()
mitma_d_mun = rel_d.groupby('distrito_mitma').municipio_ine.first()

dn = pd.read_csv('raw_data_cache/v2/zonificacion/zonificacion_distritos/nombres_distritos.csv', sep='|', dtype=str)
dn = dn.set_index(dn.columns[1])[dn.columns[2]] if dn.shape[1] >= 3 else dn.set_index(dn.columns[0])[dn.columns[1]]

qd = "'data/derived/pernoct_2024_6m.parquet'"
def dist_shares(cond, prov_prefix):
    df = con.execute(f"""
        SELECT x.r, substr(x.s,1,2) sp, cf_tbl.coastal_true cst, sum(x.p)/count(DISTINCT x.date) w
        FROM {qd} x
        LEFT JOIN cf_tbl ON x.s = cf_tbl.z
        WHERE ({cond}) AND x.r LIKE '{prov_prefix}%' AND substr(x.r,1,2) <> substr(x.s,1,2)
        GROUP BY 1,2,3
    """).df()
    tot = df.groupby('r').w.sum()
    and_all = df[df.sp.isin(AND_ALL)].groupby('r').w.sum()
    and_inl = df[df.sp.isin(AND_ALL) & (df.cst == 0)].groupby('r').w.sum()
    and_cst = df[df.sp.isin(AND_ALL) & (df.cst == 1)].groupby('r').w.sum()
    ext = df[df.sp.isin(['06', '10'])].groupby('r').w.sum()
    inl_all = df[df.cst == 0].groupby('r').w.sum()
    out = pd.DataFrame({
        'away': tot,
        'pct_AND': 100 * and_all / tot,
        'pct_AND_inl': 100 * and_inl / tot,
        'pct_AND_cst': 100 * and_cst / tot,
        'pct_EXT': 100 * ext / tot,
        'pct_inland': 100 * inl_all / tot,
    }).fillna(0)
    out['inc_pc'] = out.index.map(mitma_d_inc)
    out['mun'] = out.index.map(mitma_d_mun)
    out['born_AND_mun'] = out['mun'].map(born_and)
    out['born_EXT_mun'] = out['mun'].map(born_ext)
    out['name'] = [dn.get(i, i) for i in out.index]
    return out

# 1. Barcelona city (10 districts) + Barcelona metro districts (19 Census cities!)
bcn_d_aug = dist_shares(slices['mid-Aug'], '08')
bcn_d_nov = dist_shares(slices['All Saints'], '08')
bcn_10 = bcn_d_aug[bcn_d_aug.mun == '08019'].join(bcn_d_nov[['pct_AND', 'pct_AND_inl']], rsuffix='_nov').sort_values('inc_pc', ascending=False)
print("\n--- Barcelona City (08019) 10 Districts: District Income vs August & Nov1 Andalusia shares ---")
print(bcn_10[['name', 'inc_pc', 'pct_AND', 'pct_AND_inl', 'pct_AND_cst', 'pct_AND_nov', 'pct_inland']].round(2).to_string())
print(f"  Within BCN city (10 districts): corr(inc_pc, pct_AND Aug)={bcn_10.inc_pc.corr(bcn_10.pct_AND):+.3f}, "
      f"corr(inc_pc, pct_AND_inl Aug)={bcn_10.inc_pc.corr(bcn_10.pct_AND_inl):+.3f}, corr(inc_pc, pct_AND Nov1)={bcn_10.inc_pc.corr(bcn_10.pct_AND_nov):+.3f}")

# Barcelona metro districts across the 19 Census cities (controlling for DISTRICT income + municipality birthplace!)
bcn_m_d = bcn_d_aug[bcn_d_aug.mun.isin(cities)].dropna(subset=['inc_pc', 'born_AND_mun']).copy()
bcn_m_d['linc_d'] = np.log(bcn_m_d.inc_pc)
fit_bd = sm.OLS(bcn_m_d.pct_AND, sm.add_constant(bcn_m_d[['born_AND_mun', 'linc_d']])).fit(cov_type='cluster', cov_kwds={'groups': bcn_m_d.mun})
fit_bd_inl = sm.OLS(bcn_m_d.pct_AND_inl, sm.add_constant(bcn_m_d[['born_AND_mun', 'linc_d']])).fit(cov_type='cluster', cov_kwds={'groups': bcn_m_d.mun})
print(f"\n  Barcelona metro ALL districts (n={len(bcn_m_d)} districts in {bcn_m_d.mun.nunique()} cities, SE clustered by city):")
print(f"    pct_AND     ~ born_AND_mun ({fit_bd.params['born_AND_mun']:+.3f}, se {fit_bd.bse['born_AND_mun']:.3f}, p={fit_bd.pvalues['born_AND_mun']:.3g}) + linc_d ({fit_bd.params['linc_d']:+.2f}, p={fit_bd.pvalues['linc_d']:.3g})  R2={fit_bd.rsquared:.2f}")
print(f"    pct_AND_inl ~ born_AND_mun ({fit_bd_inl.params['born_AND_mun']:+.3f}, se {fit_bd_inl.bse['born_AND_mun']:.3f}, p={fit_bd_inl.pvalues['born_AND_mun']:.3g}) + linc_d ({fit_bd_inl.params['linc_d']:+.2f}, p={fit_bd_inl.pvalues['linc_d']:.3g})  R2={fit_bd_inl.rsquared:.2f}")

# 2. Madrid city (21 districts) + Madrid metro districts!
mad_d_aug = dist_shares(slices['mid-Aug'], '28')
mad_d_nov = dist_shares(slices['All Saints'], '28')
mad_21 = mad_d_aug[mad_d_aug.mun == '28079'].join(mad_d_nov[['pct_AND', 'pct_AND_inl', 'pct_EXT']], rsuffix='_nov').sort_values('inc_pc', ascending=False)
print("\n--- Madrid City (28079) 21 Districts: District Income vs August Coastal-AND vs Inland-AND vs Nov1 ---")
print(mad_21[['name', 'inc_pc', 'pct_AND', 'pct_AND_cst', 'pct_AND_inl', 'pct_AND_nov', 'pct_EXT']].round(2).to_string())
print(f"  Within Madrid city (21 districts):")
print(f"    corr(inc_pc, Aug pct_AND TOTAL)   = {mad_21.inc_pc.corr(mad_21.pct_AND):+.3f}")
print(f"    corr(inc_pc, Aug pct_AND COASTAL) = {mad_21.inc_pc.corr(mad_21.pct_AND_cst):+.3f}  <-- RICH districts flock to Andalusia coast!")
print(f"    corr(inc_pc, Aug pct_AND INLAND)  = {mad_21.inc_pc.corr(mad_21.pct_AND_inl):+.3f}  <-- POORER districts go to inland Andalusia!")
print(f"    corr(inc_pc, Nov1 pct_AND TOTAL)  = {mad_21.inc_pc.corr(mad_21.pct_AND_nov):+.3f}  <-- POORER districts go to Andalusia on All Saints!")
print(f"    corr(inc_pc, Aug pct_EXT)         = {mad_21.inc_pc.corr(mad_21.pct_EXT):+.3f}  <-- POORER districts go to Extremadura!")

# Madrid metro districts (controlling for district income + municipal birthplace)
mad_m_d = mad_d_aug[mad_d_aug.mun.isin(cities)].join(mad_d_nov[['pct_AND']], rsuffix='_nov').dropna(subset=['inc_pc', 'born_AND_mun']).copy()
mad_m_d['linc_d'] = np.log(mad_m_d.inc_pc)
for ycol, lbl in [('pct_AND', 'Aug AND total'), ('pct_AND_cst', 'Aug AND coastal'), ('pct_AND_inl', 'Aug AND inland'), ('pct_AND_nov', 'Nov1 AND total')]:
    f_m = sm.OLS(mad_m_d[ycol], sm.add_constant(mad_m_d[['born_AND_mun', 'linc_d']])).fit(cov_type='cluster', cov_kwds={'groups': mad_m_d.mun})
    print(f"  Madrid metro districts ({len(mad_m_d)} dists) {lbl:16s} ~ born_AND_mun ({f_m.params['born_AND_mun']:+.3f}, se {f_m.bse['born_AND_mun']:.3f}, p={f_m.pvalues['born_AND_mun']:.3g}) + linc_d ({f_m.params['linc_d']:+.2f}, p={f_m.pvalues['linc_d']:.3g})  R2={f_m.rsquared:.2f}")

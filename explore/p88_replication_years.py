# Check 2: Replication over years — 2022, 2023, and 2024 MITMA municipal overnight stays
# Re-runs BOTH:
#   (A) Village-refill gradient (p73): inland zones Q1..Q5 + OLS (prov FE + lres + linc) for 2022, 2023, 2024
#   (B) City roots elasticity (p72): 151 cities (matched to year's Censo 2022/2023/2024 birthplace + year's ADRH income)
#       for Feb weekday, mid-Aug, All Saints (2022, 2024; note MITMA outage 2023-10-30..2023-11-03), and Christmas Eve!
#   (C) Residence stability check across 2022 & 2023 as well (for Check 7 verification).
import glob, os, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, 'explore')
import duckdb, numpy as np, pandas as pd, pyfixest as pf
from common import nodes

con = duckdb.connect()
con.execute("""
    CREATE TABLE rel AS
    SELECT DISTINCT distrito_mitma AS dz, municipio_mitma AS mz
    FROM read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv',
                  delim='|', header=true, all_varchar=true)
""")

# Check downloaded 2022 & 2023 files
for yr in ['2022', '2023']:
    files = sorted(glob.glob(f'_local/external/mitma_pernoct/{yr}*_Pernoctaciones_municipios.csv.gz'))
    print(f"Year {yr}: {len(files)} daily csv.gz files on disk.")

# Build compact derived parquet for 2022-2023 municipal-level overnight stays
out_pq = 'data/derived/pernoct_mun_2022_2023.parquet'
files_22_23 = sorted(
    glob.glob('_local/external/mitma_pernoct/2022*_Pernoctaciones_municipios.csv.gz') +
    glob.glob('_local/external/mitma_pernoct/2023*_Pernoctaciones_municipios.csv.gz')
)
assert len(files_22_23) >= 115, f"Expected ~126 files for 2022-2023, found {len(files_22_23)}"

# Verify schema & unmapped residence codes
con.execute(f"""
    CREATE TEMP VIEW raw_22_23 AS
    SELECT strptime(CAST(fecha AS VARCHAR), '%Y%m%d')::DATE AS date,
           CAST(zona_residencia AS VARCHAR) AS r,
           CAST(zona_pernoctacion AS VARCHAR) AS s,
           CAST(personas AS DOUBLE) AS p
    FROM read_csv({files_22_23}, delim='|', header=true,
                  columns={{'fecha':'VARCHAR','zona_residencia':'VARCHAR','zona_pernoctacion':'VARCHAR','personas':'DOUBLE'}})
""")
miss = con.execute("SELECT count(DISTINCT r) FROM raw_22_23 x LEFT JOIN rel ON x.r=rel.dz WHERE rel.mz IS NULL").fetchone()[0]
print("Unmapped residence district codes in 2022-2023:", miss)

con.execute(f"""
    COPY (
        SELECT x.date, rel.mz AS r, x.s, sum(x.p) AS p
        FROM raw_22_23 x JOIN rel ON x.r = rel.dz
        GROUP BY 1, 2, 3
    ) TO '{out_pq}' (FORMAT PARQUET)
""")
print("Built", out_pq, con.execute(f"SELECT year(date), count(DISTINCT date), min(date), max(date), round(avg(day_tot)/1e6,3) FROM (SELECT date, sum(p) day_tot FROM '{out_pq}' GROUP BY 1) GROUP BY 1 ORDER BY 1").fetchall())

# Load common metadata
H = pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv',
                dtype={'ine_code': str, 'province_code': str})
rel_df = pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv', sep='|', dtype=str)[
    ['municipio_ine', 'municipio_mitma']
].drop_duplicates()
H = H.merge(rel_df, left_on='ine_code', right_on='municipio_ine', how='inner')
cols = ['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960','pop_1970','pop_1981','pop_2001','pop_2021']
Z = H.groupby('municipio_mitma')[cols].sum(min_count=1)
Z['peak'] = Z[['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960']].max(axis=1)
Z['loss_peak_2021'] = 1 - Z.pop_2021 / Z.peak

cf = pd.read_csv('explore/out_p86_coastal_flags.csv', dtype={'id': str}).set_index('id')
con.register('cf_tbl', cf.reset_index().rename(columns={'id': 'z'}))

# ==============================================================================
# PART A: VILLAGE-REFILL GRADIENT BY YEAR (2022, 2023, 2024)
# ==============================================================================
print("\n" + "=" * 110)
print("PART A: VILLAGE-REFILL GRADIENT REPLICATION (2022, 2023, 2024)")
print("=" * 110)
year_cfg = [
    (2022, f"'{out_pq}'", "date BETWEEN '2022-02-01' AND '2022-02-28' AND dayofweek(date) IN (2,3,4)",
     "date BETWEEN '2022-08-10' AND '2022-08-18'", "date='2022-11-01'", "date='2022-12-24'"),
    (2023, f"'{out_pq}'", "date BETWEEN '2023-02-01' AND '2023-02-28' AND dayofweek(date) IN (2,3,4)",
     "date BETWEEN '2023-08-10' AND '2023-08-18'", "date='2023-10-29'", "date='2023-12-24'"),
    (2024, "'data/derived/pernoct_mun_2024_6m.parquet'", "date BETWEEN '2024-02-01' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",
     "date BETWEEN '2024-08-10' AND '2024-08-18'", "date='2024-11-01'", "date='2024-12-24'"),
]
refill_rows = []
for yr, pq, c_feb, c_aug, c_nov, c_dec in year_cfg:
    R = con.execute(f"SELECT r z, sum(p)/count(DISTINCT date) resid FROM {pq} WHERE {c_feb} GROUP BY 1").df().set_index('z').resid
    R_aug = con.execute(f"SELECT r z, sum(p)/count(DISTINCT date) resid_aug FROM {pq} WHERE {c_aug} GROUP BY 1").df().set_index('z').resid_aug
    P_feb = con.execute(f"SELECT s z, sum(p)/count(DISTINCT date) v FROM {pq} WHERE {c_feb} GROUP BY 1").df().set_index('z').v
    P_aug = con.execute(f"SELECT s z, sum(p)/count(DISTINCT date) v FROM {pq} WHERE {c_aug} GROUP BY 1").df().set_index('z').v
    P_nov = con.execute(f"SELECT s z, sum(p)/count(DISTINCT date) v FROM {pq} WHERE {c_nov} GROUP BY 1").df().set_index('z').v
    P_dec = con.execute(f"SELECT s z, sum(p)/count(DISTINCT date) v FROM {pq} WHERE {c_dec} GROUP BY 1").df().set_index('z').v
    D = pd.DataFrame({'resid': R, 'resid_aug': R_aug, 'p_feb': P_feb, 'p_aug': P_aug, 'p_nov': P_nov, 'p_dec': P_dec})
    D = D[D.index.str.match(r'^\d')].join(Z, how='inner').join(cf, how='inner')
    D['prov'] = D.index.str[:2]
    D = D[(D.resid >= 200) & ~D.prov.isin(['07','35','38','51','52'])].copy()
    D['mult_aug'] = D.p_aug / D.resid
    D['mult_nov'] = D.p_nov / D.resid
    D['mult_dec'] = D.p_dec / D.resid
    D['mult_feb'] = D.p_feb / D.resid
    D['lmult_aug'] = np.log(D.mult_aug)
    D['lmult_nov'] = np.log(D.mult_nov)
    D['lmult_dec'] = np.log(D.mult_dec)
    D['lres'] = np.log(D.resid)
    D['loss_rank'] = D.loss_peak_2021.rank(pct=True)

    # Check residence ratio Aug/Feb for Check 7!
    print(f"  [Check 7 yr={yr}] national resid Feb={D.resid.sum()/1e6:.3f}M, Aug={D.resid_aug.sum()/1e6:.3f}M (ratio={D.resid_aug.sum()/D.resid.sum():.4f}, median zone ratio={(D.resid_aug/D.resid).median():.4f})")

    for cflag, clab in [('coastal_p73', 'p73 inland'), ('coastal_true', 'true inland')]:
        inl = D[D[cflag] == 0].copy()
        inl['loss_q'] = pd.qcut(inl.loss_peak_2021, 5, labels=['Q1','Q2','Q3','Q4','Q5'])
        qt = inl.groupby('loss_q', observed=True)[['mult_aug','mult_nov','mult_dec','mult_feb']].median()
        m_a = pf.feols("lmult_aug ~ loss_rank + lres | prov", data=inl, vcov={'CRV1':'prov'}).tidy()
        m_n = pf.feols("lmult_nov ~ loss_rank + lres | prov", data=inl, vcov={'CRV1':'prov'}).tidy()
        m_d = pf.feols("lmult_dec ~ loss_rank + lres | prov", data=inl, vcov={'CRV1':'prov'}).tidy()
        refill_rows.append(dict(
            year=yr, flag=clab, n=len(inl),
            q1_aug=qt.loc['Q1','mult_aug'], q5_aug=qt.loc['Q5','mult_aug'],
            q1_nov=qt.loc['Q1','mult_nov'], q5_nov=qt.loc['Q5','mult_nov'],
            q5_dec=qt.loc['Q5','mult_dec'], q5_feb=qt.loc['Q5','mult_feb'],
            ols_aug=m_a.loc['loss_rank','Estimate'], se_aug=m_a.loc['loss_rank','Std. Error'],
            ols_nov=m_n.loc['loss_rank','Estimate'], se_nov=m_n.loc['loss_rank','Std. Error'],
            ols_dec=m_d.loc['loss_rank','Estimate'], se_dec=m_d.loc['loss_rank','Std. Error'],
        ))
df_ref = pd.DataFrame(refill_rows)
df_ref.to_csv('explore/out_p88_refill_years.csv', index=False)
pd.set_option('display.width', 240)
print(df_ref.round(3).to_string(index=False))

# ==============================================================================
# PART B: CITY ROOTS ELASTICITY BY YEAR (2022, 2023, 2024) — WITH & WITHOUT INCOME
# ==============================================================================
print("\n" + "=" * 110)
print("PART B: 151-CITY ROOTS ELASTICITY REPLICATION (2022, 2023, 2024; MATCHED-YEAR CENSO + ADRH)")
print("=" * 110)
B_all = pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv',
                    dtype={'ine_code': str, 'prov_birth': str})
adrh_m = pd.read_csv('_local/external/ine_adrh/adrh_municipios.csv', dtype={'ine_code': str})
N = nodes().set_index('id')
M = N[~N.foreign].copy()
M['w'] = M.population.fillna(0) + 1
pc = M.groupby('prov').apply(lambda g: pd.Series({'x': np.average(g.x, weights=g.w), 'y': np.average(g.y, weights=g.w)}))

cities = sorted(set(B_all.ine_code) & set(N.index))
cl = ",".join(f"'{c}'" for c in cities)
base_grid = pd.MultiIndex.from_product([cities, [f'{i:02d}' for i in range(1, 51)]], names=['c', 'p']).to_frame(index=False)
base_grid['cprov'] = base_grid.c.str[:2]
base_grid = base_grid[base_grid.cprov != base_grid.p].copy()
base_grid['km'] = np.hypot(N.loc[base_grid.c, 'x'].values - pc.loc[base_grid.p, 'x'].values,
                           N.loc[base_grid.c, 'y'].values - pc.loc[base_grid.p, 'y'].values) / 1000.0
base_grid['dbin'] = np.digitize(base_grid.km, np.r_[0:100:10, 100:400:20, 400:1200:50])
base_grid['pairFE'] = base_grid.cprov + '_' + base_grid.p

city_rows = []
for yr, pq, c_feb, c_aug, c_nov, c_dec in year_cfg:
    By = B_all[(B_all.year == yr) & B_all.prov_birth.str.match(r'^\d\d$')]
    Bm = By.set_index(['ine_code', 'prov_birth']).persons
    inc_y = adrh_m[(adrh_m.year == yr) & (adrh_m.indicator == 'Renta neta media por persona')].set_index('ine_code').value.reindex(cities)
    linc_y = np.log(inc_y) - np.log(inc_y).mean()

    g0 = base_grid.copy()
    g0['lborn'] = np.log1p(Bm.reindex(list(zip(g0.c, g0.p))).fillna(0).values)
    g0['linc'] = g0.c.map(linc_y).values

    yr_slices = [('Feb weekday', c_feb), ('mid-Aug', c_aug), ('All Saints' if yr != 2023 else 'Oct 29 (Nov1 missing)', c_nov), ('Christmas Eve', c_dec)]
    for slab, cond in yr_slices:
        for mode, where_s in [('All zones', '1=1'), ('Inland only', 'cf_tbl.coastal_true = 0')]:
            if where_s == '1=1':
                f = con.execute(f"SELECT r c, substr(s,1,2) p, sum(p)/count(DISTINCT date) w FROM {pq} WHERE ({cond}) AND r IN ({cl}) AND r<>s GROUP BY 1,2").df().set_index(['c', 'p']).w
            else:
                f = con.execute(f"SELECT x.r c, substr(x.s,1,2) p, sum(x.p)/count(DISTINCT x.date) w FROM {pq} x JOIN cf_tbl ON x.s=cf_tbl.z WHERE ({cond}) AND x.r IN ({cl}) AND x.r<>x.s AND ({where_s}) GROUP BY 1,2").df().set_index(['c', 'p']).w
            g = g0.copy(); g['w'] = f.reindex(list(zip(g.c, g.p))).fillna(0).values
            t1 = pf.fepois("w ~ lborn | c + p + dbin", data=g, vcov={'CRV1': 'c'}).tidy()
            t2 = pf.fepois("w ~ lborn | c + pairFE", data=g, vcov={'CRV1': 'c'}).tidy()
            t3 = pf.fepois("w ~ lborn | c + pairFE + p[linc]", data=g, vcov={'CRV1': 'c'}).tidy()
            city_rows.append(dict(
                year=yr, mode=mode, slice=slab, flow_M=g.w.sum()/1e6,
                distFE=t1.loc['lborn','Estimate'], dist_se=t1.loc['lborn','Std. Error'],
                pairFE=t2.loc['lborn','Estimate'], pair_se=t2.loc['lborn','Std. Error'],
                pair_inc=t3.loc['lborn','Estimate'], pair_inc_se=t3.loc['lborn','Std. Error'],
            ))

df_c = pd.DataFrame(city_rows)
df_c.to_csv('explore/out_p88_roots_years.csv', index=False)
print(df_c.round(3).to_string(index=False))

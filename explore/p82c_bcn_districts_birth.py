# Check 1E: Within-Barcelona City (10 districts) — exact district-level birthplace (Resto de España) + ADRH district income
import os, time, requests, duckdb, numpy as np, pandas as pd, statsmodels.api as sm

raw_csv = '_local/external/bcn_padron_birth/raw/2025_pad_mdbas_lloc-naix_sexe.csv'
url = 'https://opendata-ajuntament.barcelona.cat/data/dataset/7dcf78e6-8444-401d-a744-257ffc0e1a9b/resource/3da24718-fd60-44f8-b5c1-c305be7debe7/download'
if not os.path.exists(raw_csv) or os.path.getsize(raw_csv) < 10000:
    for attempt in range(1, 5):
        try:
            r = requests.get(url, timeout=30, headers={'User-Agent': 'Mozilla/5.0'})
            r.raise_for_status()
            open(raw_csv, 'wb').write(r.content)
            break
        except Exception as e:
            print(f"retry {attempt}: {e}")
            time.sleep(2 * attempt)

b = pd.read_csv(raw_csv)
b['Valor'] = pd.to_numeric(b['Valor'], errors='coerce').fillna(0)
g = b.groupby(['Codi_Districte', 'Nom_Districte', 'LLOC_NAIX'])['Valor'].sum().unstack(fill_value=0)
g['tot'] = g.sum(axis=1)
g['pct_born_restES'] = 100 * g[3] / g['tot']
g['pct_born_CAT'] = 100 * (g[1] + g[2]) / g['tot']
g['pct_born_abroad'] = 100 * (g[4] + g[5]) / g['tot']
g = g.reset_index()
g['r'] = g['Codi_Districte'].apply(lambda x: f'08019{int(x):02d}')
g = g.set_index('r')

adrh_d = pd.read_csv('_local/external/ine_adrh/adrh_distritos.csv', dtype={'dist_code': str})
inc = adrh_d[(adrh_d.year == 2024) & (adrh_d.indicator == 'Renta neta media por persona')].set_index('dist_code').value
g['inc_pc'] = inc.reindex(g.index)
g['linc'] = np.log(g['inc_pc'])

con = duckdb.connect()
cf = pd.read_csv('explore/out_p86_coastal_flags.csv', dtype={'id': str}).rename(columns={'id': 'z'})
con.register('cf_tbl', cf)
AND = ['04','11','14','18','21','23','29','41']
CAT = ['08','17','25','43']
for lab, cond in [('aug', "date BETWEEN '2024-08-10' AND '2024-08-18'"), ('nov1', "date='2024-11-01'")]:
    df = con.execute(f"""
        SELECT x.r, substr(x.s,1,2) sp, cf_tbl.coastal_true cst, sum(x.p)/count(DISTINCT x.date) w
        FROM 'data/derived/pernoct_2024_6m.parquet' x LEFT JOIN cf_tbl ON x.s=cf_tbl.z
        WHERE ({cond}) AND x.r LIKE '08019%' AND substr(x.s,1,2)<>'08'
        GROUP BY 1,2,3
    """).df()
    tot = df.groupby('r').w.sum()
    out_cat = df[~df.sp.isin(CAT)].groupby('r').w.sum()
    out_cat_inl = df[~df.sp.isin(CAT) & (df.cst == 0)].groupby('r').w.sum()
    and_all = df[df.sp.isin(AND)].groupby('r').w.sum()
    and_inl = df[df.sp.isin(AND) & (df.cst == 0)].groupby('r').w.sum()
    g[f'pct_outCAT_{lab}'] = 100 * out_cat / tot
    g[f'pct_outCAT_inl_{lab}'] = 100 * out_cat_inl / tot
    g[f'pct_AND_{lab}'] = 100 * and_all / tot
    g[f'pct_AND_inl_{lab}'] = 100 * and_inl / tot

pd.set_option('display.width', 220)
print("=== WITHIN BARCELONA CITY (10 DISTRICTS): EXACT DISTRICT BIRTHPLACE (Resto de España) + ADRH INCOME ===")
print(g[['Nom_Districte','inc_pc','pct_born_restES','pct_born_abroad','pct_outCAT_aug','pct_outCAT_inl_aug','pct_AND_aug','pct_AND_inl_aug','pct_AND_nov1']].sort_values('inc_pc', ascending=False).round(2).to_string())
print(f"\ncorr(pct_born_restES, linc) = {g.pct_born_restES.corr(g.linc):+.3f}")
for y in ['pct_outCAT_aug','pct_outCAT_inl_aug','pct_AND_aug','pct_AND_inl_aug','pct_outCAT_nov1','pct_AND_nov1']:
    r0 = g[y].corr(g.pct_born_restES)
    fit = sm.OLS(g[y], sm.add_constant(g[['pct_born_restES','linc']])).fit()
    print(f"  {y:20s}: r_raw={r0:+.3f} | b_restES={fit.params['pct_born_restES']:+.3f} (se {fit.bse['pct_born_restES']:.3f}, p={fit.pvalues['pct_born_restES']:.3g}) | b_linc={fit.params['linc']:+.2f} (p={fit.pvalues['linc']:.3g}) | R2={fit.rsquared:.2f}")

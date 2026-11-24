# Check 7: MITMA residence-definition caveat
# Do total residents per zone in August drop in cities and jump in villages (reassignment)?
# Quantify magnitude and direction of bias for both village refill and city roots elasticity.
import sys
sys.path.insert(0, 'explore')
import duckdb, numpy as np, pandas as pd, geopandas as gpd

con = duckdb.connect()
q = "'data/derived/pernoct_mun_2024_6m.parquet'"

# 1. National total residents by date
daily_nat = con.execute(f"""
    SELECT date, strftime(date, '%Y-%m') ym, sum(p) tot_p,
           sum(CASE WHEN r<>s THEN p ELSE 0 END) away_p
    FROM {q}
    GROUP BY 1,2 ORDER BY 1
""").df()
print("=== National daily totals by month (millions) ===")
print(daily_nat.groupby('ym').agg(
    days=('tot_p', 'size'),
    tot_min=('tot_p', lambda s: s.min()/1e6),
    tot_mean=('tot_p', lambda s: s.mean()/1e6),
    tot_max=('tot_p', lambda s: s.max()/1e6),
    away_pct_mean=('away_p', lambda s: 100*(s/daily_nat.loc[s.index,'tot_p']).mean()),
).round(3).to_string())

# 2. Zone-level resident count (sum(p) grouped by r) across months
slices = {
    'feb': "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",
    'mar': "date BETWEEN '2024-03-04' AND '2024-03-21' AND dayofweek(date) IN (2,3,4)",
    'aug_early': "date BETWEEN '2024-08-01' AND '2024-08-07'",
    'aug_mid': "date BETWEEN '2024-08-10' AND '2024-08-18'",
    'aug_late': "date BETWEEN '2024-08-24' AND '2024-08-31'",
    'oct': "date BETWEEN '2024-10-01' AND '2024-10-25' AND dayofweek(date) IN (2,3,4)",
    'nov1': "date='2024-11-01'",
    'nov_mid': "date BETWEEN '2024-11-12' AND '2024-11-28' AND dayofweek(date) IN (2,3,4)",
    'dec24': "date='2024-12-24'",
}
res_df = {}
pres_df = {}
for k, cond in slices.items():
    res_df[k] = con.execute(f"SELECT r z, sum(p)/count(DISTINCT date) v FROM {q} WHERE {cond} GROUP BY 1").df().set_index('z').v
    pres_df[k] = con.execute(f"SELECT s z, sum(p)/count(DISTINCT date) v FROM {q} WHERE {cond} GROUP BY 1").df().set_index('z').v

R = pd.DataFrame(res_df)
P = pd.DataFrame(pres_df)
R = R[R.index.str.match(r'^\d')].copy()

# Load historical population & coastal flag (exact p73 definition)
H = pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv',
                dtype={'ine_code': str, 'province_code': str})
rel = pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv', sep='|', dtype=str)[
    ['municipio_ine', 'municipio_mitma']
].drop_duplicates()
H = H.merge(rel, left_on='ine_code', right_on='municipio_ine', how='inner')
cols = ['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960','pop_1970','pop_1981','pop_2001','pop_2021','padron_2024']
Z = H.groupby('municipio_mitma')[cols].sum(min_count=1)
Z['peak'] = Z[['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960']].max(axis=1)
Z['loss_peak_2021'] = 1 - Z.pop_2021 / Z.peak

g = gpd.read_file('clean_data/v2/zones/municipios_mitma.gpkg')
outer = g.geometry.union_all().boundary
g['coastal'] = (g.geometry.buffer(50).intersects(outer) & g.id.str.match(r'^\d')).astype(int)

D = R.join(Z[['padron_2024', 'loss_peak_2021']], how='inner').join(g.set_index('id')[['coastal']])
D['pres_aug'] = P['aug_mid'].reindex(D.index)
D['pres_nov1'] = P['nov1'].reindex(D.index)
D['prov'] = D.index.str[:2]
D = D[(D['feb'] >= 200) & ~D.prov.isin(['07', '35', '38', '51', '52'])].copy()

# Ratios of assigned residents in each period vs Feb baseline
for k in ['mar', 'aug_early', 'aug_mid', 'aug_late', 'oct', 'nov1', 'nov_mid', 'dec24']:
    D[f'r_ratio_{k}'] = D[k] / D['feb']
D['r_vs_padron_feb'] = D['feb'] / D['padron_2024']
D['r_vs_padron_aug'] = D['aug_mid'] / D['padron_2024']

print("\n=== Zone-level resident count ratio (assigned residents in period / Feb residents), mainland N =", len(D), "===")
for k in ['mar', 'aug_early', 'aug_mid', 'aug_late', 'oct', 'nov1', 'nov_mid', 'dec24']:
    s = D[f'r_ratio_{k}']
    print(f"  {k:10s}: mean={s.mean():.4f}  p05={s.quantile(.05):.4f}  p25={s.quantile(.25):.4f}  median={s.median():.4f}  p75={s.quantile(.75):.4f}  p95={s.quantile(.95):.4f}  sum_ratio={D[k].sum()/D['feb'].sum():.4f}")

print("\n=== Assigned residents ratio (Aug_mid / Feb) by coastal vs inland exodus quintiles ===")
inl = D[D.coastal == 0].copy()
inl['loss_q'] = pd.qcut(inl.loss_peak_2021, 5, labels=['Q1 grew', 'Q2', 'Q3', 'Q4', 'Q5 emptied'])
print(inl.groupby('loss_q', observed=True).agg(
    n=('feb', 'size'),
    res_feb_M=('feb', lambda s: s.sum()/1e6),
    res_aug_M=('aug_mid', lambda s: s.sum()/1e6),
    ratio_aug_feb_med=('r_ratio_aug_mid', 'median'),
    ratio_aug_feb_mean=('r_ratio_aug_mid', 'mean'),
    ratio_nov1_feb_med=('r_ratio_nov1', 'median'),
    ratio_feb_padron_med=('r_vs_padron_feb', 'median'),
    ratio_aug_padron_med=('r_vs_padron_aug', 'median'),
    mult_aug_febres=('pres_aug', lambda s: (s / inl.loc[s.index, 'feb']).median()),
    mult_aug_augres=('pres_aug', lambda s: (s / inl.loc[s.index, 'aug_mid']).median()),
    mult_aug_padron=('pres_aug', lambda s: (s / inl.loc[s.index, 'padron_2024']).median()),
).round(4).to_string())

cst = D[D.coastal == 1]
print(f"\nCoastal zones (n={len(cst)}): median r_ratio_aug_mid={cst.r_ratio_aug_mid.median():.4f} (p95={cst.r_ratio_aug_mid.quantile(.95):.4f}), "
      f"ratio_aug_padron_med={cst.r_vs_padron_aug.median():.4f}, sum_aug/sum_feb={cst.aug_mid.sum()/cst.feb.sum():.4f}")

# Top 15 biggest cities (like Madrid, Barcelona, Valencia, Bilbao, etc.)
nm = pd.read_csv('raw_data_cache/v2/zonificacion/zonificacion_municipios/nombres_municipios.csv', sep='|', dtype=str).set_index('ID')['name']
big = D.sort_values('feb', ascending=False).head(15).copy()
big['name'] = nm.reindex(big.index).values
print("\n=== Top 15 biggest cities: assigned MITMA residents Feb vs Aug vs Oct vs Padron 2024 ===")
print(big[['name', 'padron_2024', 'feb', 'aug_mid', 'oct', 'r_ratio_aug_mid', 'r_ratio_oct', 'r_ratio_nov1', 'r_vs_padron_feb', 'r_vs_padron_aug']].round(4).to_string())

# Biggest gainers and losers in assigned residents (Aug_mid / Feb)
D['name'] = nm.reindex(D.index.str.replace('_AM', '')).values
print("\n=== Top 15 zones gaining assigned residents in Aug vs Feb ===")
print(D.sort_values('r_ratio_aug_mid', ascending=False).head(15)[['name', 'prov', 'coastal', 'feb', 'aug_mid', 'r_ratio_aug_mid', 'loss_peak_2021']].round(3).to_string())
print("\n=== Top 15 zones losing assigned residents in Aug vs Feb ===")
print(D.sort_values('r_ratio_aug_mid', ascending=True).head(15)[['name', 'prov', 'coastal', 'feb', 'aug_mid', 'r_ratio_aug_mid', 'loss_peak_2021']].round(3).to_string())

# Check daily resident totals in Madrid (28079), Barcelona (08019), and inland Q5 sum across all 180 days!
q5_ids = inl[inl.loss_q == 'Q5 emptied'].index.tolist()
q5_sql = ",".join(f"'{z}'" for z in q5_ids)
ts = con.execute(f"""
    SELECT date,
           sum(CASE WHEN r='28079' THEN p ELSE 0 END) mad_res,
           sum(CASE WHEN r='08019' THEN p ELSE 0 END) bcn_res,
           sum(CASE WHEN r IN ({q5_sql}) THEN p ELSE 0 END) q5_res,
           sum(CASE WHEN s IN ({q5_sql}) THEN p ELSE 0 END) q5_pres,
           sum(CASE WHEN r IN ({q5_sql}) AND s IN ({q5_sql}) THEN p ELSE 0 END) q5_self
    FROM {q}
    GROUP BY 1 ORDER BY 1
""").df()
ts['ym'] = ts['date'].astype(str).str[:7]
print("\n=== Monthly averages of resident counts (Madrid, BCN, inland Q5) and Q5 present/self ===")
print(ts.groupby('ym')[['mad_res', 'bcn_res', 'q5_res', 'q5_pres', 'q5_self']].mean().round(0).to_string())
ts.to_csv('explore/out_p87_residence_ts.csv', index=False)
D.to_csv('explore/out_p87_residence_zones.csv')

import pandas as pd, numpy as np
n=pd.read_parquet('data/derived/nodes_municipios.parquet'); m=n[~n.foreign&~n.island&~n.prov.isin(['51','52'])]
print('NaN-pop mainland nodes:', m[m.population.isna()][['id','name']].head(10).values.tolist())
rel = pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv', sep='|', dtype=str)
r2=rel[['municipio_ine','municipio_mitma']].drop_duplicates(); d=r2.groupby('municipio_ine').municipio_mitma.nunique(); print('split INE muni:', r2[r2.municipio_ine.isin(d[d>1].index)].values.tolist())
H = pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv', dtype={'ine_code': str})
cols=['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960','pop_1970','pop_1981','pop_2001','pop_2021']
M=H.merge(r2,left_on='ine_code',right_on='municipio_ine'); Z=M.groupby('municipio_mitma')[cols].sum(min_count=1)
loss=1-Z.pop_2021/Z[cols[:7]].max(axis=1)
nol=m[m.id.map(loss).isna()]
print('no-loss mainland zones:', len(nol), 'examples', nol[['id','name','population']].head(8).values.tolist())
print('no-loss in rel?', nol.id.isin(r2.municipio_mitma).sum())
# compare pop_2021 aggregated vs nodes population
c=m.set_index('id').population.to_frame().join(Z.pop_2021)
print('ratio nodes pop / H 2021 quantiles', (c.population/c.pop_2021).quantile([.01,.5,.99]).round(3).to_dict())
print('nodes population total mainland', m.population.sum()/1e6)

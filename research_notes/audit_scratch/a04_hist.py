import pandas as pd, numpy as np
H = pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv', dtype={'ine_code': str, 'province_code': str})
print(H.shape, H.columns.tolist()[:30])
rel = pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv', sep='|', dtype=str)[['municipio_ine','municipio_mitma']].drop_duplicates()
print('rel ine codes mapping to >1 mitma zone:', (rel.groupby('municipio_ine').municipio_mitma.nunique()>1).sum())
print('ine codes in rel not in H:', len(set(rel.municipio_ine)-set(H.ine_code)), 'H not in rel:', len(set(H.ine_code)-set(rel.municipio_ine)))
cols=['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960','pop_1970','pop_1981','pop_2001','pop_2021']
M=H.merge(rel,left_on='ine_code',right_on='municipio_ine')
print('missing share per year', M[cols].isna().mean().round(3).to_dict())
part = M.groupby('municipio_mitma')[cols].apply(lambda d: ((d.isna().any()) & (d.notna().any())).any())
print('zones where some but not all constituent municipalities miss a census year:', int(part.sum()), 'of', M.municipio_mitma.nunique())
n=pd.read_parquet('data/derived/nodes_municipios.parquet'); m=n[~n.foreign&~n.island&~n.prov.isin(['51','52'])]
print('mainland zones', len(m), 'with partial-missing history', int(part.reindex(m.id).fillna(False).sum()))
Z=M.groupby('municipio_mitma')[cols].sum(min_count=1)
print('mainland zones without loss:', m.id.map(1-Z.pop_2021/Z[cols[:7]].max(axis=1)).isna().sum())
print('population NaN in nodes (mainland):', m.population.isna().sum(), 'pop<200:', (m.population<200).sum())
print(n.columns.tolist())

import pandas as pd, duckdb
n=pd.read_parquet('data/derived/nodes_municipios.parquet'); m=n[~n.foreign&~n.island&~n.prov.isin(['51','52'])]
m=m.assign(L=m.id.str.len()); print(m.L.value_counts().to_dict())
ph=m[m.L>8]; print('phantom ids', len(ph), 'prov assigned', ph.prov.value_counts().to_dict(), 'pop', ph.population.describe()[['count','max']].to_dict())
con=duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"
zs=set(con.execute(f"select distinct s from {q}").df().s)|set(con.execute(f"select distinct r from {q}").df().r)
print('phantom ids appearing in pernoct:', len(set(ph.id)&zs), '; mainland ids never in pernoct:', len(set(m.id)-zs))
print('mainland ids in pernoct:', len(set(m.id)&zs))
c=pd.read_csv('explore/out_p86_coastal_flags.csv',dtype={'id':str}).set_index('id')
print('phantoms flagged coastal by p73:', int(c.loc[ph.id,'coastal_p73'].sum()))

import pandas as pd, geopandas as gpd
c=pd.read_csv('explore/out_p86_coastal_flags.csv',dtype={'id':str}).set_index('id')
n=pd.read_parquet('data/derived/nodes_municipios.parquet').set_index('id'); m=n[~n.foreign&~n.island&~n.prov.isin(['51','52'])]
z=pd.read_csv('data/derived/summer_territories/zones.csv',dtype={'id':str})
z=z.set_index(z.columns[0])
cc=c.loc[m.index].join(z['coastal'].rename('coastal_p100'))
print(pd.crosstab(cc.coastal_true,cc.coastal_p100))
g=gpd.read_file('clean_data/v2/zones/municipios_mitma.gpkg').set_index('id')
fn=cc[(cc.coastal_true==0)&(cc.coastal_p100==True)].index; fp=cc[(cc.coastal_true==1)&(cc.coastal_p100==False)].index
print('p100 coastal but not p86:', [(i,g.loc[i,'name']) for i in fn][:40])
print('p86 coastal but not p100:', [(i,g.loc[i,'name']) for i in fp][:40])
fpb=cc[(cc.coastal_p73==1)&(cc.coastal_true==0)].index
print('p73 false positives by prov', pd.Series([i[:2] for i in fpb]).value_counts().head(15).to_dict())
print('foreign zones in gpkg', (~g.index.str.match(r'^\d')).sum(), g[~g.index.str.match(r'^\d')].index[:10].tolist())

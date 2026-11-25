import duckdb, pandas as pd, numpy as np
nd=pd.read_parquet('data/derived/nodes_municipios.parquet'); nd=nd[~nd.foreign&~nd.island&~nd.prov.isin(['51','52'])]
bm=pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv',dtype={'ine_code':str})
cities=sorted(set(bm[bm.year==2024].ine_code)&set(nd.id)); cl=",".join(f"'{c}'" for c in cities); ml=",".join(f"'{z}'" for z in nd.id)
con=duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"
for lab,cond,k in [('Feb_wd',"date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",12),('Mid_Aug',"date BETWEEN '2024-08-10' AND '2024-08-18'",9)]:
    d=con.execute(f"SELECT r,s,sum(p) sp,count(DISTINCT date) nd FROM {q} WHERE {cond} AND r IN ({cl}) AND s IN ({ml}) AND r<>s GROUP BY 1,2").df()
    b=d.sp/d.nd; c=d.sp/k
    print(lab,'BiCM links (w>=5) buggy',int((b>=5).sum()),'correct',int((c>=5).sum()),'| total flow buggy %.3fM correct %.3fM'%(b.sum()/1e6,c.sum()/1e6))

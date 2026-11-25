import pandas as pd, numpy as np, duckdb
adrh=pd.read_csv('_local/external/ine_adrh/adrh_municipios.csv',dtype={'ine_code':str})
inc=adrh[(adrh.year==2024)&(adrh.indicator=='Renta neta media por persona')].set_index('ine_code').value
B=pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv',dtype={'ine_code':str,'prov_birth':str})
n=pd.read_parquet('data/derived/nodes_municipios.parquet').set_index('id')
cities=sorted(set(B[B.year==2024].ine_code)&set(n.index)); print('cities',len(cities),'missing ADRH 2024:',[c for c in cities if pd.isna(inc.get(c))])
print('island/CeutaMelilla cities in 151:', [c for c in cities if c[:2] in ('07','35','38','51','52')].__len__())
ml=[c for c in cities if c[:2] not in ('07','35','38','51','52')]; s=pd.Series([c[:2] for c in ml]).value_counts()
print('mainland cities',len(ml),'provinces with >=2 cities',(s>=2).sum(),'cities in them',s[s>=2].sum(),'single-city provs',(s==1).sum())
con=duckdb.connect(); q="'data/derived/pernoct_mun_2022_2023.parquet'"
d=con.execute(f"select strftime(date,'%Y-%m') m, count(distinct date) n, min(date), max(date) from {q} group by 1 order by 1").df(); print(d.to_string())
d2=set(pd.to_datetime(con.execute(f"select distinct date from {q}").df().date).dt.date)
for y in (2022,2023):
    for rng in [pd.date_range(f'{y}-02-01',f'{y}-02-28'),pd.date_range(f'{y}-08-01',f'{y}-08-31'),pd.date_range(f'{y}-10-25',f'{y}-11-05')]:
        miss=sorted(set(rng.date)-d2); print(y, rng[0].date(), 'missing:', miss)

import runpy
from pathlib import Path
import numpy as np
import pandas as pd
import duckdb
from scipy.spatial.distance import cdist

audit=Path(__file__).resolve().parent
root=runpy.run_path(str(audit/"roots.py"))
nodes=root["nodes"]
f=root["controls"](root["r2_data"](),"pop")
p=root["panel"](f,["all"])["all"]
zones=set(nodes.prov)
coastal_provs=set(nodes.loc[nodes.coast,"prov"])
print("INLAND provinces",len(zones-coastal_provs),flush=True)
early=["2022-02-16","2022-08-15","2022-11-01","2023-02-15","2023-08-15"]
db=duckdb.connect(config={"threads":4})
db.register("nodes",nodes)
city_indices=tuple(int(i) for i in f.idx.unique())
year_data=db.sql(f"""SELECT CAST(date AS VARCHAR) date,i,prov,SUM(w) as w
                    FROM '{audit/"earlier.parquet"}' JOIN nodes ON j=idx
                    WHERE date IN {tuple(early)} AND i in {city_indices} GROUP BY ALL""").df()
year=year_data.pivot(index=["i","prov"],columns="date",values="w").reindex(p.index).fillna(0)
rows=[]
for scope,table,allowed in (("earlier",year,np.ones(len(f),bool)),("inland_provinces",p,~f.prov.isin(coastal_provs))):
    for d in early if scope=="earlier" else [*root["days"],"oct_mean"]:
        y=table[d].to_numpy()
        keep=allowed & pd.Series(y,index=f.index).groupby(f.origin+f.prov).transform("sum").gt(0)
        ff=f[keep]
        b,_,_=root["solve"](root["design"](ff,["ln_roots"]),y[keep])
        rows.append(dict(scope=scope,date=d,beta=b[0],N=len(ff),G=ff.origin.nunique()))
pd.DataFrame(rows).to_csv(audit/"r2_extra.csv",index=False)

# Same-city roots index across destination provinces with city and distance removed.
all_cities=pd.read_csv("results/r2_city_set.csv",dtype={"ine_code":str,"prov":str}).rename(columns={"ine_code":"id","prov":"origin","income_2024":"income"})
all_roots=pd.read_csv("results/r2_roots_city_prov.csv",dtype={"ine_code":str,"prov_birth":str}).rename(columns={"ine_code":"id","prov_birth":"prov"})
design=root["controls"](all_cities.merge(nodes[["id","idx","x","y"]],on="id").merge(
    all_roots[all_roots.prov.isin(zones)&~all_roots.own_province],on="id").reset_index(drop=True),"pop")
daily=root["panel"](design,["all"])["all"]
city=pd.get_dummies(design.id).to_numpy(dtype=float)
X=np.column_stack([city,np.log1p(np.exp(design.ld))])
q,_=np.linalg.qr(X,mode="reduced")
res=lambda v:v-q@(q.T@v)
y=np.log1p(daily.to_numpy())
born=res(np.log1p(design.born.to_numpy()))
fitted=res(y)
correlation=(born@fitted)/np.linalg.norm(born)/np.linalg.norm(fitted,axis=0)
pd.DataFrame({"date":daily.columns,"corr_within_city_over_province":correlation}).to_csv(audit/"roots_index.csv",index=False)

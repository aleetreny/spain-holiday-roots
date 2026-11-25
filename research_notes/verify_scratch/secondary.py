from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

audit=Path(__file__).resolve().parent
nodes=pd.read_parquet(audit/"nodes.parquet")
prov_groups=["04 11 14 18 21 23 29 41","22 44 50","33","07","35 38","39","05 09 24 34 37 40 42 47 49",
             "02 13 16 19 45","08 17 25 43","03 12 46","06 10","15 27 32 36","28","30","31","01 20 48",
             "26","51","52"]
prov_to_ccaa={code:f"{idx:02d}" for idx,group in enumerate(prov_groups,1) for code in group.split()}

ece=Path("_local/external/ine_ecepov_2021")
homes=pd.read_parquet(ece/"ecepov_vivienda_2021.parquet")
people=pd.read_parquet(ece/"ecepov_hogar_2021.parquet")
home=homes[(homes.SEGUNRESI=="1")&(homes.LUGSEGUNRESI=="4")].copy()
home["origin"]=home.IDQ_PV.map(prov_to_ccaa)
p=people[people.IDEN.isin(home.IDEN)].merge(home[["IDEN","OTRACCAA"]],on="IDEN")
births=p[["PROVN","PROVNP1","PROVNP2"]].map(prov_to_ccaa.get)
p["g1"]=births.PROVN.eq(p.OTRACCAA)
p["g2"]=births[["PROVNP1","PROVNP2"]].eq(p.OTRACCAA,axis=0).any(axis=1)
nac=p[["NACIMP1","NACIMP2"]]
p["missing_parent"]=(nac.isin(["0"," ","","9","99"]) | (nac.eq("1") & births[["PROVNP1","PROVNP2"]].isna().set_axis(nac.columns,axis=1))).any(axis=1)
home=home.join(p.groupby("IDEN")[["g1","g2","missing_parent"]].max(),on="IDEN")
home["any_match"]=home.g1|home.g2
inland={"02","07","08","11","13","15","17"}
for label,sub in [("all",home),("inland",home[home.OTRACCAA.isin(inland)]),
                  ("inland excl dest Madrid",home[home.OTRACCAA.isin(inland-{"13"})]),
                  ("inland excl origin Madrid",home[home.OTRACCAA.isin(inland)&home.origin.ne("13")])]:
    print("ECEPOV",label,len(sub), "weight", sub.FACTOR.sum(),
          "g1 g2 any missing, madrid",
          np.average(sub[["g1","g2","any_match","missing_parent"]].astype(float),axis=0,weights=sub.FACTOR),
          np.average(sub.origin=="13",weights=sub.FACTOR))

birth_flags=(births.assign(IDEN=p.IDEN).melt(id_vars="IDEN").dropna()
             .drop_duplicates(["IDEN","value"]))
for label,sub in [("all",home),("inland",home[home.OTRACCAA.isin(inland)])]:
    observed=np.average(sub.any_match,weights=sub.FACTOR)
    for ref in (sub,home) if label=="inland" else (sub,):
        global_perm=ref is home and label=="inland"
        pop_ref=ref if global_perm else sub
        chance=(ref.groupby("IDQ_PV").OTRACCAA.value_counts(normalize=True).rename("prob")
                .reset_index().rename(columns={"OTRACCAA":"value"}))
        if global_perm:
            chance=chance[chance["value"].isin(inland)]
        exposure=(chance.groupby("IDQ_PV").prob.sum())
        per_home=birth_flags.merge(pop_ref[["IDEN","FACTOR","IDQ_PV"]],on="IDEN")
        linked=per_home.merge(chance,on=["IDQ_PV","value"])
        expected=(linked.FACTOR*linked.prob).sum()/(pop_ref.FACTOR*pop_ref.IDQ_PV.map(exposure).fillna(0)).sum()
        print("ECEPOV null",label,"global" if global_perm else "within",expected,observed/expected)

raw="_local/external/ine_censo_viviendas_2021/raw/censo2021_59525.csv"
housing=pd.read_csv(raw,sep=";",thousands=".").dropna(subset=["Municipios"]).rename(columns={"Tipo de vivienda (principal o no)":"type"})
housing["ine_code"]=housing.Municipios.str[:5]
housing=housing.pivot(index="ine_code",columns="type",values="Total")
print("housing fields",housing.columns.tolist())
local=pd.read_parquet("results/n1_mun_nodes.parquet",columns=["ine_code","prov","padron_2024","zone_id","area_km2"])
raw_geo=Path("raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv")
crosswalk=pd.read_csv(raw_geo,sep="|",dtype=str)[["municipio_ine","municipio_mitma"]].drop_duplicates()
hist=pd.read_csv("_local/external/ine_historical_pop/municipal_population_1900_2025.csv",dtype={"ine_code":str}).set_index("ine_code")
historic_cols=[f"pop_{yr}" for yr in range(1900,1970,10)]
local["loss"]=local.ine_code.map(1-hist.pop_2021/hist[historic_cols].max(axis=1))
xy=nodes[["x","y"]]
nodes["citydist"]=cdist(xy,xy[nodes.city]).min(axis=1)/1000
nodes["coastdist"]=cdist(xy,xy[nodes.coast]).min(axis=1)/1000
local=local.merge(nodes[["id","citydist","coastdist","coast","n_mun"]],left_on="zone_id",right_on="id")
local=local.join(housing,on="ine_code")
local["housing_share"]=local["Vivienda no principal"]/local["Total"]
local["rank"]=local.loss.rank(pct=True)
print("municipal controls", local[["loss","housing_share","area_km2"]].count())
c=Path("_local/external/ine_turismo_moviles/raw")
target_periods=[f"{year}M{month:02d}" for year in [2022,2023,2024,2025] for month in [2,3,8,11,12]]
for fname,outcol in [("turismo_53001_provorig_x_mundest.csv.gz","in"),("turismo_53000_munorig_x_provdest.csv.gz","out")]:
    frame=pd.read_csv(c/fname,sep=";",dtype=str)
    print("INE raw",list(frame),frame.shape)
    if outcol=="in":
        m=frame.iloc[:,0].eq("Total Nacional") & frame.iloc[:,1].isna() & frame.iloc[:,2].isna() & frame.iloc[:,4].eq("Turistas")
        part=frame.loc[m,[frame.columns[3],"Periodo","Total"]].copy()
    else:
        m=frame.iloc[:,1].eq("Total Nacional") & frame.iloc[:,2].isna() & frame.iloc[:,3].isna() & frame.iloc[:,4].eq("Turistas")
        part=frame.loc[m,[frame.columns[0],"Periodo","Total"]].copy()
    part=part[part.Periodo.isin(target_periods)]
    part["code"]=part.iloc[:,0].str[:5]
    part["Total"]=pd.to_numeric(part.Total.str.replace(".","",regex=False),errors="coerce")
    print(part.iloc[:2],part.shape)
    mat=part.pivot(index="code",columns="Periodo",values="Total").add_prefix(outcol)
    local=local.join(mat,on="ine_code")

def partial_share_interval(f, rural_rank, y, reduced, mediator, seed):
    values=np.column_stack([rural_rank,reduced,mediator,y])
    sufficient=[]
    for indices in f.reset_index(drop=True).groupby("prov").indices.values():
        centered=values[indices]-values[indices].mean(axis=0)
        sufficient.append(centered.T@centered)
    sufficient=np.array(sufficient)
    N=len(sufficient)
    rng=np.random.default_rng(seed)
    n_draws=int(sys.argv[1]) if len(sys.argv)>1 else 1999
    selected=rng.integers(N,size=(n_draws,N))
    counts=np.array([np.bincount(i,minlength=N) for i in selected])
    increments=np.tensordot(counts,sufficient,axes=1)
    base=np.linalg.solve(increments[:,:-2,:-2],increments[:,:-2,-1,None])[:,0,0]
    added=np.linalg.solve(increments[:,:-1,:-1],increments[:,:-1,-1,None])[:,0,0]
    return np.quantile(1-added/base,[.025,.975])

def linear(f,y,drops=None):
    variables=pd.DataFrame({"rank":f["rank"],"logpop":np.log(f.padron_2024),
                            "coast":f.coast,"citydist":np.log1p(f.citydist),
                            "coastdist":np.log1p(f.coastdist),
                            "area":np.log(f.area_km2),"housing":f.housing_share})
    if drops:
        variables=variables.drop(columns=drops)
    d=pd.concat([variables,pd.get_dummies(f.prov,drop_first=True)],axis=1)
    X=np.column_stack((np.ones(len(d)),d.to_numpy(float)))
    return np.linalg.lstsq(X,y,rcond=None)[0][1]

for term in ["2024M03","2024M11","2024M12"]:
    a,b="in"+term,"in2024M02"
    f=local[(local[a]>0)&(local[b]>0)&local["rank"].notna() & local.area_km2.notna()]
    y=np.log(f[a]/f[b])
    p0,p1=linear(f,y,["housing"]),linear(f,y)
    print("INE additional housing",term,len(f),p0,p1,1-p1/p0)
    if term=="2024M12":
        inputs=np.column_stack([np.log(f.padron_2024),f.coast,np.log1p(f.citydist),
                                np.log1p(f.coastdist),np.log(f.area_km2)])
        print("INE Dec interval",partial_share_interval(f,f["rank"],y,inputs,f.housing_share,38974))

for year in (2022,2023,2024,2025):
    a,b,oa,ob=(f"in{year}M08",f"in{year}M02",f"out{year}M08",f"out{year}M02")
    f=local.dropna(subset=[a,b,oa,ob,"rank"])
    for typ,z in [("all",f),("single",f[f.n_mun.eq(1)]),("AM",f[f.n_mun.gt(1)])]:
        print("INE net",year,typ,len(z),linear(z,(z[a]-z[oa]-z[b]+z[ob])/z.padron_2024,["area","housing"]))
    a,b=f"in{year}M08",f"in{year}M02"
    f=local[(local[a]>0)&(local[b]>0)&local["rank"].notna()].copy()
    for typ,data in [("all",f),("single",f[f.n_mun.eq(1)]),("AM",f[f.n_mun.gt(1)]),("pop<1k",f[f.padron_2024.lt(1000)])]:
        y=np.log(data[a]/data[b])
        print("INE log slope",year,typ,len(data),linear(data,y,["area","housing"]))
    for area in [True,False]:
        sub=f[f.housing_share.notna()&f.area_km2.notna()]
        drop=[] if area else ["area"]
        y=np.log(sub[a]/sub[b])
        p0=linear(sub,y,drop+["housing"])
        p1=linear(sub,y,drop)
        print("INE housing",year,"area",area,len(sub),p0,p1,1-p1/p0)
        if year==2024 and area:
            inputs=np.column_stack([np.log(sub.padron_2024),sub.coast,np.log1p(sub.citydist),
                                    np.log1p(sub.coastdist),np.log(sub.area_km2)])
            print("INE Aug interval",partial_share_interval(sub,sub["rank"],y,inputs,sub.housing_share,417))

# Rebuild zone housing shares from nationally tabulated counts and the raw crosswalk.
zone_homes=(crosswalk.merge(housing,left_on="municipio_ine",right_index=True)
            .groupby("municipio_mitma")[["Total","Vivienda no principal"]].sum())
m=nodes.loc[nodes.loss.notna() & ~nodes.id.isin(["31905_AM","19103_AM","40149_AM"])].copy()
m["rank"]=m.loss.rank(pct=True)
m["housing_share"]=m.id.map(zone_homes["Vivienda no principal"]/zone_homes["Total"])
m=m.rename(columns={"pop":"padron_2024","area":"area_km2"})
daily=pd.read_parquet(audit/"node_daily.parquet")
daily.date=pd.to_datetime(daily.date).dt.strftime("%Y-%m-%d")
daily["log_in"]=np.log1p(daily.visitors)
mat=daily.pivot(index="idx",columns="date",values="log_in")
m=m.join(mat,on="idx")
for h in ["2024-08-15","2024-11-01","2024-12-25","2024-03-29"]:
    y=m[h]-m["2024-02-14"]
    x=pd.concat([m[["rank","coast","housing_share"]],np.log(m[["padron_2024","area_km2","n_mun"]]),
                 np.log1p(m[["citydist","coastdist"]]),pd.get_dummies(m.prov)],axis=1).to_numpy(float)
    reduced=np.linalg.lstsq(np.delete(x,2,axis=1),y,rcond=None)[0][0]
    full=np.linalg.lstsq(x,y,rcond=None)[0][0]
    print("MITMA housing",h,len(m),reduced,full,1-full/reduced)

    print("MITMA housing interval",h,partial_share_interval(m,x[:,0],y,x[:,[1,3,4,5,6,7]],x[:,2],4571))

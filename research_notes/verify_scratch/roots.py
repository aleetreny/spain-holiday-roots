import os
os.environ["VECLIB_MAXIMUM_THREADS"]="1"
os.environ["OPENBLAS_NUM_THREADS"]="1"
import sys
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import brentq

audit=Path(__file__).resolve().parent
nodes=pd.read_parquet(audit/"nodes.parquet")
dates_feb=[d.strftime("%Y-%m-%d") for d in pd.date_range("2024-02-01","2024-02-29") if d.dayofweek in (1,2,3)]
dates_aug=[f"2024-08-{d:02}" for d in range(10,19)]
dates_nov=[f"2024-11-{d:02}" for d in range(1,4)]
dates_oct=[d.strftime("%Y-%m-%d") for d in pd.date_range("2024-10-01","2024-10-31") if d.dayofweek in (1,2,3)]
days=["2024-02-14","2024-08-15","2024-11-01","2024-12-25","2024-03-29"]
all_dates=sorted(set(dates_feb+dates_aug+dates_nov+dates_oct+days))
flows=pd.read_parquet(audit/"city_prov.parquet")
flows.date=flows.date.dt.strftime("%Y-%m-%d")
flows=flows[flows.date.isin(all_dates)]
income=pd.read_csv("_local/external/ine_adrh/adrh_municipios.csv",dtype={"ine_code":str})
income=income.loc[income.indicator.eq("Renta neta media por persona")].sort_values("year").drop_duplicates("ine_code",keep="last").set_index("ine_code")
print("income source",income.year.value_counts().to_dict(),flush=True)

def indicator(codes, values=None):
    codes=np.asarray(codes)
    active=codes>=0
    vals=np.ones(sum(active)) if values is None else np.asarray(values)[active]
    return sparse.csr_matrix((vals,(np.flatnonzero(active),codes[active])),
                             shape=(len(codes),int(codes.max())+1))

def design(f, terms):
    omitted=set(f.groupby("origin").id.min())
    city=pd.factorize(f.id.where(~f.id.isin(omitted)),sort=True)[0]
    pair=pd.factorize(f.origin+f.prov,sort=True)[0]
    dest=pd.factorize(f.prov.where(f.prov!=min(f.prov)),sort=True)[0]
    features=f[terms+["ld","ld2"]].to_numpy()
    return sparse.hstack([sparse.csr_matrix(features),indicator(city),indicator(pair),indicator(dest,f.li)],format="csr")

def solve(X,y,start=None,offset=0.):
    b=np.zeros(X.shape[1]) if start is None else start.copy()
    intercept=np.log(y.mean()) if start is None else 0.0
    obj_old=np.inf
    for step in range(80):
        eta=X@b+intercept+offset
        mu=np.exp(eta)
        err=mu-y
        delta=spsolve((X.T@X.multiply(mu[:,None])).tocsc(),X.T@err,permc_spec="MMD_AT_PLUS_A")
        update=X@delta
        if max(abs(update))<1e-7:
            break
        alpha=1.
        obj_old=(mu-y*eta).sum()
        while alpha>1e-9:
            trial=eta-alpha*update
            if ((np.exp(np.clip(trial,-700,700))-y*trial).sum() <= obj_old+1e-8):
                break
            alpha/=2
        b-=alpha*delta
    else:
        raise RuntimeError("no convergence")
    return b,mu,intercept

def controls(f, method):
    target=nodes[["prov","x","y","pop"]].copy()
    if method=="pop":
        target[["x","y"]]=target[["x","y"]].mul(target["pop"],axis=0)
        center=target.groupby("prov")[["x","y","pop"]].sum()
        center=center[["x","y"]].div(center["pop"],axis=0)
    else:
        center=target.groupby("prov")[["x","y"]].mean()
    f=f.merge(center,how="left",left_on="prov",right_index=True,suffixes=("","_p"))
    f["ld"]=np.log(np.hypot(f.x-f.x_p,f.y-f.y_p)/1000)
    f["ld2"]=f.ld**2
    f["li"]=np.log(f.income)-np.log(f.income).mean()
    return f

def census_data(transform=None):
    db=duckdb.connect(config={"threads":4,"temp_directory":str(audit/"duckdb_tmp")})
    blocks=db.sql("""SELECT CPRO||CMUN id, CPRO_NAC prov, count(*) born,
        count_if(VARANOP='0001') pre1981,count_if(VARANOP between '2016' and '2023') recent,
        count_if(CMUN_NAC='991') village, count_if(CMUN_NAC not in ('991','992','993','999','000','')) urban,
        count_if(VARANOP='0001' and VAREDAD>=65) pre65,
        count_if(VARANOP='0001' and VAREDAD<65) pre_u65,
        count_if(VARANOP!='0001' and VAREDAD>=65) post65,
        count_if(VARANOP between '2016' and '2023' and VAREDAD<65) recent_u65,
        count_if(VAREDAD>=65) born65,
        count_if(CMUN_NAC='991' and VAREDAD<65) village_u65,
        count_if(CMUN_NAC='991' and VAREDAD>=65) village65,
        count_if(CMUN_NAC not in ('991','992','993','999','000','') and VAREDAD<65) urban_u65,
        count_if(CMUN_NAC not in ('991','992','993','999','000','') and VAREDAD>=65) urban65
        FROM '_local/external/ine_censo_microdatos_2024/md_CensoAnual_muestra_2024.parquet'
        WHERE CMUN not in ('991','992','993','999','000','')
        GROUP BY ALL""").df()
    ids=blocks.id.unique()
    cities=nodes.loc[nodes.id.isin(ids)&nodes.n_mun.eq(1),["id","idx","x","y","prov"]].rename(columns={"prov":"origin"})
    cities["income"]=cities.id.map(income.value)
    cities=cities.dropna(subset=["income"])
    cities=cities[cities.groupby("origin").id.transform("size")>1]
    f=cities.merge(nodes[["prov"]].drop_duplicates(),how="cross")
    f=f[f.origin!=f.prov].merge(blocks,how="left",on=["id","prov"]).reset_index(drop=True)
    if transform is None:
        transform=lambda counts,total: np.log(counts+0.5)
    totals=db.sql("""SELECT CPRO||CMUN id,count(*) n
                     FROM '_local/external/ine_censo_microdatos_2024/md_CensoAnual_muestra_2024.parquet'
                     GROUP BY ALL""").df().set_index("id").n
    for col in blocks.columns[2:]:
        f[col]=transform(f[col].fillna(0),f.id.map(totals))
    print("V1 towns, provinces, cells",len(cities),cities.origin.nunique(),len(f),flush=True)
    return f

def r2_data():
    cities=(pd.read_csv("results/r2_city_set.csv",dtype={"ine_code":str,"prov":str})
            .rename(columns={"ine_code":"id","prov":"origin","income_2024":"income"}))
    cities=cities[cities.identifying].merge(nodes[["id","idx","x","y"]],on="id")
    roots=(pd.read_csv("results/r2_roots_city_prov.csv",dtype={"ine_code":str,"prov_birth":str})
           .rename(columns={"ine_code":"id","prov_birth":"prov"}))
    f=cities.merge(roots[roots.prov.isin(nodes.prov)&~roots.own_province],on="id")
    f["ln_roots"]=np.log(f.share)
    return f[np.isfinite(f.ln_roots)].reset_index(drop=True)

def panel(f, scopes):
    ff=flows[flows.i.isin(f.idx)].copy()
    out={}
    for scope in scopes:
        sub=ff if scope=="all" else ff[ff.coast.eq(scope=="coast")]
        grouped=sub.groupby(["i","prov","date"]).w.sum().unstack(fill_value=0)
        out[scope]=grouped.reindex(pd.MultiIndex.from_frame(f[["idx","prov"]]),fill_value=0)
        out[scope]["feb_mean"]=out[scope][dates_feb].mean(axis=1)
        out[scope]["oct_mean"]=out[scope][dates_oct].mean(axis=1)
    return out

def run(name, mode, extensive):
    spec=({"R2":["ln_roots"]} if name=="R2" else {
          "cohort":["pre1981","recent"],"birth":["village","urban"],
          "age":["pre65","post65","pre_u65","recent_u65"],
          "horserace":["pre1981","recent","born65"],"youngvillage":["village_u65","village65","urban_u65","urban65"]})
    f=controls(r2_data() if name=="R2" else census_data(),mode)
    scopes=["all","inland","coast"] if extensive else ["all"]
    panel_by_scope=panel(f,scopes)
    out=[]
    for scope in scopes:
        p=panel_by_scope[scope]
        if name=="R2" and extensive:
            nights=all_dates+["feb_mean","oct_mean"]
        else:
            nights=all_dates if extensive else days[:2]
        for label,terms in spec.items():
            if label in ("age","horserace","youngvillage"):
                use=days[:2]
            else:
                use=nights
            for night in use:
                y=p[night].to_numpy()
                active=pd.Series(y,index=f.index).groupby(f.origin+f.prov).transform("sum").gt(0).to_numpy()
                ff=f[active]
                b,mu,_=solve(design(ff,terms),y[active])
                out.append(dict(model=label,scope=scope,date=night,n=len(ff),G=ff.origin.nunique(),**dict(zip(terms,b))))
                if extensive and scope=="all" and label in ("cohort","birth","R2") and night in days[:2]:
                    for province in ff.origin.unique():
                        keep=active&f.origin.ne(province)
                        j=f[keep]
                        b,mu,_=solve(design(j,terms),y[keep])
                        out.append(dict(model=label,scope=scope,date=night,drop=province,n=len(j),G=j.origin.nunique(),**dict(zip(terms,b))))
            print(name,mode,scope,label,[(r["date"],[r[t] for t in terms],r["n"]) for r in out
                  if r["model"]==label and r["scope"]==scope and r["date"] in days[:2] and "drop" not in r],flush=True)
            pd.DataFrame(out).to_csv(audit/f"{name}_{mode}.csv",index=False)

if __name__=="__main__":
    run(sys.argv[1],sys.argv[2],"full" in sys.argv)

from pathlib import Path
from itertools import product
import numpy as np
import pandas as pd
import igraph as ig
from scipy.spatial.distance import cdist

audit=Path(__file__).resolve().parent
nodes=pd.read_parquet(audit/"nodes.parquet")
raw=pd.read_parquet(audit/"sister_w10.parquet")
daily=pd.read_parquet(audit/"node_daily.parquet")
main={"2024-02-14","2024-08-15","2024-11-01"}
sister_dates=sorted(raw.date.unique())
sister_dates=[x.strftime("%Y-%m-%d") for x in sister_dates]
for key in ("2024-12-25","2024-03-29"):
    extra=pd.read_parquet(audit/f"edges_{key}.parquet")
    raw=pd.concat([raw,extra.loc[extra.w.ge(10)].assign(date=pd.Timestamp(key))])
raw.date=pd.to_datetime(raw.date).dt.strftime("%Y-%m-%d")
daily.date=pd.to_datetime(daily.date).dt.strftime("%Y-%m-%d")
earlier=pd.read_parquet(audit/"earlier.parquet",filters=[("date","in",[
    pd.Timestamp(x).date() for x in ["2022-02-16","2022-08-15","2022-11-01","2023-02-15","2023-08-15"]])])
earlier.date=pd.to_datetime(earlier.date).dt.strftime("%Y-%m-%d")
early_daily=[]
for night,df in earlier.groupby("date"):
    early_daily.append(pd.DataFrame({"date":night,"idx":nodes.idx,
                                    "visitors":np.bincount(df.j,df.w,minlength=len(nodes)),
                                    "travelers":np.bincount(df.i,df.w,minlength=len(nodes)),
                                    "in_degree":np.bincount(df.j,df.w.ge(10).astype(float),minlength=len(nodes)),
                                    "out_degree":np.bincount(df.i,df.w.ge(10).astype(float),minlength=len(nodes))}))
daily=pd.concat([daily,*early_daily],ignore_index=True)
raw=pd.concat([raw,earlier[earlier.w.ge(10)]],ignore_index=True)
equal=pd.concat([pd.read_parquet(audit/f"edges_{date}.parquet").nlargest(35812,columns="w")
                 .assign(date=date) for date in main])
lookup={}
nr=[]
for mode,edges in (("w10",raw),("equal_edges",equal)):
    for day,b in edges.groupby("date"):
        n=len(nodes)
        g=ig.Graph(n,b[["i","j"]].to_numpy(),directed=True)
        u=g.as_undirected()
        k=np.asarray(u.degree())
        order=np.lexsort((nodes.idx,-k))
        prefix=0
        for end in range(1,len(order)):
            if not set(order[:end]).issubset(u.neighbors(int(order[end]))):
                break
            prefix=end+1
        nr.append(dict(date=day,mode=mode,prefix=prefix,core=max(u.coreness()),
                       C=u.transitivity_avglocal_undirected(mode="zero"),phi50=u.induced_subgraph(order[:50]).density(),
                       diam=u.connected_components().giant().diameter(),away=daily.loc[daily.date==day,"visitors"].sum()))
        if mode=="w10":
            lookup[day]=(u.coreness(),g.coreness(mode="in"))
pd.DataFrame(nr).to_csv(audit/"sister_network.csv",index=False)
for day in set(sister_dates)|set(lookup):
    select=daily.date==day
    for col,val in zip(["core","incore"],lookup[day]):
        daily.loc[select,col]=np.asarray(val)[daily.loc[select,"idx"]]

db=nodes.copy()
xy=db[["x","y"]]
db["log_city"]=np.log1p(cdist(xy,xy[db.city]).min(axis=1)/1000)
db["log_coast"]=np.log1p(cdist(xy,xy[db.coast]).min(axis=1)/1000)
for col in ["pop","n_mun","area"]:
    db["log_"+col]=np.log(db[col])
db=db[db.loss.notna()&~db.id.isin(["31905_AM","19103_AM","40149_AM"])]
db["rank"]=db.loss.rank(pct=True)
s=daily.merge(db,on="idx")
s["net"]=(s.visitors-s.travelers)/(s.visitors+s.travelers)
s["net_rate"]=(s.visitors-s.travelers)/s["pop"]
s["importer"]=1.*s.visitors.gt(s.travelers)
s["log_vis"]=np.log1p(s.visitors)
s["log_in"]=np.log1p(s.in_degree)
s["log_out"]=np.log1p(s.out_degree)
s["log_away_rate"]=np.log(np.maximum(s.travelers,1)/s["pop"])
s["log_vis_rate"]=np.log(s.visitors.where(s.visitors.gt(0))/s["pop"])
s["log_share"]=np.log(s.visitors.where(s.visitors.gt(0))/s.groupby("date").visitors.transform("sum"))
covs=["log_pop","log_n_mun","log_area","coast","log_city","log_coast"]
dates={"Feb":[x for x in sister_dates if x[5:7]=="02"],
       "Aug":[x for x in sister_dates if x[5:7]=="08"],
       "Nov":[x for x in sister_dates if x[5:7]=="11"]}
level_rows=[]
anc_rows=[]

for scope,pop in [("all",s),("single",s[s.n_mun.eq(1)])]:
    f=pop[pop.date=="2024-02-14"].sort_values("idx")
    rank=f["rank"].to_numpy()
    raw_x=pd.concat([f[covs],pd.get_dummies(f.prov,drop_first=True)],axis=1).to_numpy(float)
    if scope=="single":
        raw_x=raw_x[:,np.arange(raw_x.shape[1])!=1]
    X=np.column_stack((np.ones(len(f)),raw_x))
    rank_error=rank-X@np.linalg.lstsq(X,rank,rcond=None)[0]
    den=rank_error@rank
    outcome_cols=["core","incore","log_vis","log_in","log_out","in_degree","net","net_rate","importer",
                  "log_away_rate","log_vis_rate","log_share"]
    for day, data in pop.groupby("date"):
        block=data.sort_values("idx")
        row=dict(scope=scope,date=day,volume=daily.loc[daily.date==day,"visitors"].sum())
        for col in outcome_cols:
            y=block[col].to_numpy()
            good=np.isfinite(y)
            if good.all():
                value=(rank_error@y)/den
            elif good.any():
                value=np.linalg.lstsq(np.column_stack((rank[good],X[good])),y[good],rcond=None)[0][0]
            else:
                value=np.nan
            row[col]=value
        level_rows.append(row)
        if day=="2024-02-14":
            continue
        baseline=f.log_vis.to_numpy()
        y=block.log_vis.to_numpy()
        anc_rows.append(dict(scope=scope,date=day,beta=np.linalg.lstsq(np.column_stack((rank,X,baseline)),y,rcond=None)[0][0]))
levels=pd.DataFrame(level_rows)
levels.to_csv(audit/"r1_level.csv",index=False)
pd.DataFrame(anc_rows).to_csv(audit/"r1_ancova.csv",index=False)
results=[]
for scope,block in levels.groupby("scope"):
    b=block.set_index("date")
    for period in ("Aug","Nov"):
        pairs=list(product(dates[period],dates["Feb"]))
        for col in outcome_cols:
            vals=np.array([b.loc[h,col]-b.loc[o,col] for h,o in pairs])
            results.append(dict(scope=scope,period=period,outcome=col,n=len(vals),
                                positive=(vals>0).sum(),negative=(vals<0).sum(),
                                lo=np.nanmin(vals),hi=np.nanmax(vals),median=np.nanmedian(vals)))
pd.DataFrame(results).to_csv(audit/"r1_sisters.csv",index=False)

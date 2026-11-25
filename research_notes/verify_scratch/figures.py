from pathlib import Path
import numpy as np
import pandas as pd
import igraph as ig
from scipy.spatial.distance import cdist
from scipy.stats import t
from scipy.linalg import qr

audit=Path(__file__).resolve().parent
nodes=pd.read_parquet(audit/"nodes.parquet")
calendar=pd.read_parquet(audit/"calendar.parquet")
days=["2024-02-14","2024-03-29","2024-08-15","2024-11-01","2024-12-25"]
key=calendar[calendar.date.isin(pd.to_datetime(days))].pivot(index="date",columns="kind")
print("Calendar main")
print(key.to_string())
for dist in ["dist50","dist100","dist200","centroid100"]:
    tab=calendar.pivot(index="date",columns="kind",values=dist)
    chosen=tab.loc[days]
    share=chosen.inland/chosen.sum(axis=1)
    print(dist,"ratios Aug/Feb all, coast, inland:",
          chosen.sum(axis=1).iloc[2]/chosen.sum(axis=1).iloc[0], (chosen.iloc[2]/chosen.iloc[0]).to_dict(),
          "shares ",share.tolist(),"inland/coast:",(chosen.inland/chosen.coast).tolist())
    if dist=="dist100":
        reference_share=share
        ref=pd.read_csv("reports/figures/fig1_calendar_data.csv",parse_dates=["date"]).set_index("date")
        print("calendar ref dates",ref.index.to_series().groupby(ref.index.month).size().to_dict())
        print("rebuild by month")
        for month, block in tab.groupby(tab.index.month):
            print(month,len(block),int((block.inland>block.coast).sum()),int((block.coast>block.inland).sum()))
        cmp=tab.join(ref,how="outer",rsuffix="_y")
        print("calendar max common abs discrepancies",
              (cmp.inland-cmp.towns).abs().max(),(cmp.coast-cmp.coast_y).abs().max() if "coast_y" in cmp else None)
    if dist in ["dist50","dist200"]:
        base=key["dist100"].inland/key["dist100"].sum(axis=1)
        print("vs 100km percentage point range", ((share-base)*100).to_dict())

d=pd.read_parquet(audit/"node_daily.parquet")
d["date"]=d.date.dt.strftime("%Y-%m-%d")
d=d.merge(nodes,on="idx")
d["net"]=(d.visitors-d.travelers)/(d.visitors+d.travelers)
d["importer"]=d.visitors>d.travelers
d["visitor_rate"]=d.visitors/d["pop"]
d["net_rate"]=(d.visitors-d.travelers)/d["pop"]
d["away_log"]=np.log1p(d.travelers/d["pop"])
xy=nodes[["x","y"]]
d["citykm"]=cdist(xy,xy[nodes.city]).min(axis=1)[d.idx]/1000
d["coastkm"]=cdist(xy,xy[nodes.coast]).min(axis=1)[d.idx]/1000
loss=d[["idx","loss"]].drop_duplicates().set_index("idx").loss
eligible=nodes.loss.notna()&~nodes.id.isin(["31905_AM","19103_AM","40149_AM"])
rank=pd.Series(nodes.loss[eligible].rank(pct=True),index=nodes.index[eligible])
d["rank"]=d.idx.map(rank)
d["quintile"]=d.idx.map(pd.qcut(nodes.loss[eligible].rank(method="first"),5,labels=False)+1)
print("node data",len(nodes),sum(nodes.city),sum(nodes.coast),eligible.sum(),sum(eligible&(nodes.n_mun==1)))
for night in days:
    z=d[d.date==night]
    print(night,"overall importers",z.importer.sum(),z.importer.mean(),
          "Madrid",z.loc[z.id=="28079","net"].iloc[0],
          "quintiles",z.groupby("quintile").importer.agg(["sum","size","mean"]).to_dict(),
          "loss >= 0.5 inland",z[z["rank"].notna() & ~z.coast & (z.loss>=0.5)].importer.agg(["sum","size","mean"]).to_dict(),
          "single",z[z["rank"].notna() & ~z.coast & (z.loss>=0.5)&z.n_mun.eq(1)].importer.agg(["sum","size","mean"]).to_dict())
    edges=pd.read_parquet(audit/f"edges_{night}.parquet")
    velez=edges[(edges.j==nodes.idx[nodes.id=="29094"].iloc[0])&(edges.w>0)].merge(nodes,left_on="i",right_on="idx")
    print("Velez:",velez.w.sum(),len(velez),velez.loc[velez.prov=="29","w"].sum()/velez.w.sum(),
          list(velez.groupby("prov").w.sum().nlargest(8).index),list(velez.nlargest(8,"w").name))
    g=ig.Graph(len(nodes),edges.loc[edges.w>=10,["i","j"]].to_numpy(),directed=True)
    mask=d.date==night
    d.loc[mask,"core"]=np.array(g.as_undirected().coreness())[d.loc[mask,"idx"]]
    d.loc[mask,"incore"]=np.array(g.coreness(mode="in"))[d.loc[mask,"idx"]]

for scope,select in [("all", d["rank"].notna()),("single",d["rank"].notna()&d.n_mun.eq(1))]:
    sdf=d[select].copy()
    sdf["log_vis"]=np.log1p(sdf.visitors)
    sdf["log_in"]=np.log1p(sdf.in_degree)
    sdf["away_log"]=np.log(np.maximum(sdf.travelers,1)/sdf["pop"])
    sdf["vis_log_rate"]=np.log(np.maximum(sdf.visitors,1)/sdf["pop"])
    ref=sdf[sdf.date=="2024-02-14"].sort_values("idx")
    x=pd.concat([ref[["rank","coast"]],np.log(ref[["pop","n_mun","area"]]),
                 np.log1p(ref[["citykm","coastkm"]]),pd.get_dummies(ref.prov,drop_first=True)],axis=1).to_numpy(float)
    x=np.column_stack((np.ones(len(x)),x))
    def estimate(y):
        valid=np.isfinite(y)
        return np.linalg.lstsq(x[valid],y[valid],rcond=None)[0][1]
    for night in days:
        current=sdf[sdf.date==night].sort_values("idx")
        outcomes={}
        for var in ("importer","net","net_rate","away_log","vis_log_rate","log_vis","log_in","core","incore","in_degree"):
            level=current[var].to_numpy(dtype=float)
            delta=level-ref[var].to_numpy(dtype=float)
            outcomes[var]=[round(estimate(level),5),round(estimate(delta),5)]
        print("R1",scope,night,outcomes)

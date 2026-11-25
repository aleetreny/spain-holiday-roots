from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

audit=Path(__file__).resolve().parent
node=pd.read_parquet(audit/"nodes.parquet")
data=pd.read_parquet(audit/"node_daily.parquet")
states=node.loss.notna()&~node.id.isin(["31905_AM","19103_AM","40149_AM"])
xy=node[["x","y"]]
node["citydist"]=cdist(xy,xy[node.city]).min(axis=1)/1000
node["coastdist"]=cdist(xy,xy[node.coast]).min(axis=1)/1000
node["rank"]=node.loss[states].rank(pct=True)
data.date=pd.to_datetime(data.date).dt.strftime("%Y-%m-%d")
variables={
  "importer":data.visitors.gt(data.travelers).astype(float),
  "net":(data.visitors-data.travelers)/(data.visitors+data.travelers),
  "net_rate":(data.visitors-data.travelers)/node["pop"].to_numpy()[data.idx],
  "log_vis":np.log1p(data.visitors),
  "log_degree_in":np.log1p(data.in_degree),
  "log_away":np.log(data.travelers.clip(lower=1)/node["pop"].to_numpy()[data.idx]),
}
for name,series in variables.items():
    data[name]=series
rng=np.random.default_rng(44129)
W=rng.choice(np.r_[-np.sqrt([.5,1,1.5]),np.sqrt([.5,1,1.5])],size=(9999,47))

def estimate(f,y):
    X=np.column_stack([f["rank"],np.log(f[["pop","area"]]),f.coast,np.log1p(f[["citydist","coastdist"]])])
    if f.n_mun.nunique()>1:
        X=np.column_stack([X,np.log(f.n_mun)])
    keep=np.isfinite(y)
    v=pd.DataFrame(np.column_stack([X[keep],np.asarray(y)[keep]]))
    codes=f.prov[keep].to_numpy()
    v=v-v.groupby(codes).transform("mean")
    X=v.iloc[:,:-1].to_numpy()
    y=v.iloc[:,-1].to_numpy()
    inv=np.linalg.inv(X.T@X)
    b=inv@(X.T@y)
    group=[np.flatnonzero(codes==p) for p in np.unique(codes)]
    XX=np.stack([X[i].T@X[i] for i in group])
    adj=len(group)/(len(group)-1)*(len(y)-1)/(len(y)-X.shape[1])
    def boot(residual):
        xu=np.stack([X[i].T@residual[i] for i in group])
        draws=W[:,:len(group)]
        change=(draws@xu)@inv
        contributions=draws*(xu@inv[:,0])-change@(XX@inv[:,0]).T
        se=np.sqrt(adj*(contributions*contributions).sum(axis=1))
        sample_se=np.sqrt(adj*np.square(xu@inv[:,0]).sum())
        return change[:,0]/se,sample_se
    draws,se=boot(y-X@b)
    q=np.quantile(draws,[.025,.975])
    null=y-X[:,1:]@np.linalg.lstsq(X[:,1:],y,rcond=None)[0]
    null_draws,_=boot(null)
    p=(np.count_nonzero(np.abs(null_draws)>=abs(b[0]/se))+1)/(len(draws)+1)
    return b[0],b[0]-q[1]*se,b[0]-q[0]*se,p

main="2024-02-14"
samples={"all":states,"single":states&node.n_mun.eq(1)}
outputs=[]
for scope,mask in samples.items():
    f=node[mask].copy()
    for d in [main,"2024-08-15","2024-11-01","2024-12-25","2024-03-29"]:
        yday=data[data.date==d].set_index("idx").loc[f.idx]
        ybase=data[data.date==main].set_index("idx").loc[f.idx]
        for term in variables:
            for what in ("level","change") if d!=main else ("level",):
                y=yday[term].to_numpy()
                if what=="change":
                    y=y-ybase[term].to_numpy()
                b,lo,hi,p=estimate(f,y)
                outputs.append([scope,d,term,what,b,lo,hi,p])
pd.DataFrame(outputs,columns=["scope","date","term","type","estimate","lo","hi","p"]).to_csv(audit/"R1_wcr.csv",index=False)

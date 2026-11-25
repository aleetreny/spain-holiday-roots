from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import nbinom
import igraph as ig
import random

audit=Path(__file__).resolve().parent
nodes=pd.read_parquet(audit/"nodes.parquet")
N=len(nodes)
xy=nodes[["x","y"]]
dist=cdist(xy,xy)/1000
bins=np.searchsorted([0,5,10,15,20,25,30,40,50,60,80,100,125,150,200,250,300,400,500,600,800,1200],dist,side="right")-1
np.fill_diagonal(bins,21)
stack=bins.ravel()
province=pd.get_dummies(nodes.prov).to_numpy(float).T
eligible=nodes.loss.notna()&~nodes.id.isin(["31905_AM","19103_AM","40149_AM"])
nodes["rank"]=nodes.loss[eligible].rank(pct=True)
nodes["ldcity"]=np.log1p(dist[:,nodes.city].min(axis=1))
nodes["ldcoast"]=np.log1p(dist[:,nodes.coast].min(axis=1))
nodes["lpop"]=np.log(nodes["pop"])
nodes["ln_mun"]=np.log(nodes.n_mun)
nodes["larea"]=np.log(nodes.area)
vars=["rank","lpop","ln_mun","larea","coast","ldcity","ldcoast"]

def coef(X,y):
    keep=np.isfinite(y)&np.isfinite(X).all(axis=1)
    return np.linalg.lstsq(X[keep],y[keep],rcond=None)[0][0]

def ipf(W):
    a=W.sum(axis=1)
    b=W.sum(axis=0)
    target=np.bincount(stack,W.ravel(),minlength=22)
    estimate=np.ones_like(W)
    np.fill_diagonal(estimate,0)
    def ratio(a,b):
        return np.divide(a,b,where=b>0,out=np.zeros_like(a))
    for k in range(1,201):
        estimate*=ratio(a,estimate.sum(axis=1))[:,None]
        estimate*=ratio(b,estimate.sum(axis=0))[None,:]
        estimate*=ratio(target,np.bincount(stack,estimate.ravel(),minlength=22))[bins]
        gap=np.max(np.r_[np.abs(ratio(estimate.sum(axis=1)-a,a)),
                             np.abs(ratio(estimate.sum(axis=0)-b,b))])
        if gap<1e-8:
            break
    alpha=(((W-estimate)**2).sum()-estimate.sum())/np.square(estimate).sum()
    prob=nbinom.sf(9,1/alpha,1/(1+alpha*estimate))
    return estimate,prob,k,gap,alpha

def graph(edge):
    return ig.Graph(N,np.column_stack(np.nonzero(edge)),directed=True)

def diversity(weight):
    v=province@weight
    return np.divide(v.sum(axis=0)**2,(v*v).sum(axis=0),
                     out=np.full(N,np.nan),where=v.sum(axis=0)>0)

rng=np.random.default_rng(172151)
ig.set_random_number_generator(random.Random(882514))
metrics={}
stats=[]
for date in ("2024-02-14","2024-08-15","2024-11-01"):
    frame=pd.read_parquet(audit/f"edges_{date}.parquet")
    W=np.zeros((N,N))
    W[frame.i,frame.j]=frame.w
    g=graph(W>=10)
    u=g.as_undirected()
    T,s,it,gap,alpha=ipf(W)
    print(date,it,gap,alpha,s.sum(),flush=True)
    gm=[]
    wm=[]
    for trial in range(16):
        gm.append(graph(rng.random(s.shape)<s).coreness(mode="in"))
        h=u.copy()
        h.rewire(n=10*u.ecount(),allowed_edge_types="simple")
        wm.append(h.coreness())
    metrics[date]=dict(expected_in=np.log1p(s.sum(axis=0)),observed_in=np.log1p(g.indegree()),
                      expected_incore=np.mean(gm,axis=0),observed_incore=g.coreness(mode="in"),
                      expected_div=np.log(diversity(T)),observed_div=np.log(diversity(W)),
                      expected_core=np.mean(wm,axis=0),observed_core=u.coreness())
    stats.append(dict(date=date,iter=it,margin=gap,alpha=alpha,exp_edges=s.sum(),obs_edges=g.ecount()))
pd.DataFrame(stats).to_csv(audit/"gravity_fit.csv",index=False)
out=[]
for date in ("2024-08-15","2024-11-01"):
    for name in ("in","incore","core","div"):
        for scope,keep in [("all",eligible),("single",eligible&nodes.n_mun.eq(1))]:
            f=nodes[keep]
            X=pd.concat([f[vars],pd.get_dummies(f.prov)],axis=1).to_numpy(float)
            slopes={}
            for typ in ("expected","observed"):
                key=f"{typ}_{name}"
                y=np.asarray(metrics[date][key])-np.asarray(metrics["2024-02-14"][key])
                slopes[typ]=coef(X,y[keep])
            out.append(dict(date=date,outcome=name,scope=scope,
                            excess=slopes["observed"]-slopes["expected"],**slopes))
pd.DataFrame(out).to_csv(audit/"gravity_slopes.csv",index=False)

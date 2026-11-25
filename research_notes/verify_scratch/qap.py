from pathlib import Path
import duckdb
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import runpy

audit=Path(__file__).resolve().parent
base=runpy.run_path(str(audit/"roots.py"))
nodes=base["nodes"]
dates=base["days"]
p=nodes[["prov","pop","x","y"]].copy()
p[["x","y"]]=p[["x","y"]].mul(p["pop"],axis=0)
p=p.groupby("prov").sum()
center=p[["x","y"]].div(p["pop"],axis=0)
province_order=center.index
city=pd.read_csv("results/r2_city_set.csv",dtype={"ine_code":str,"prov":str}).sort_values("ine_code")
city_nodes=city.merge(nodes[["id","idx","x","y"]],left_on="ine_code",right_on="id")
def matrix_city(data,val,index="ine_code",dest="prov_birth"):
    return data.pivot_table(index=index,columns=dest,values=val,aggfunc="sum").reindex(index=city_nodes.ine_code,columns=province_order).fillna(0).to_numpy()
def matrix_prov(data,val,index="prov_residence",dest="prov_birth"):
    return data.pivot_table(index=index,columns=dest,values=val,aggfunc="sum").reindex(index=province_order,columns=province_order).fillna(0).to_numpy()
roots=pd.read_csv("results/r2_roots_city_prov.csv",dtype={"ine_code":str,"prov_birth":str})
proot=pd.read_csv("results/r2_roots_prov_prov.csv",dtype=str)
proot.persons=pd.to_numeric(proot.persons)
daily=base["flows"]
mat_city=[]; mat_prov=[]
orig_prov=nodes.set_index("idx").prov
for d in dates:
    evening=daily[daily.date==d]
    c=evening.merge(city_nodes[["idx","ine_code"]],left_on="i",right_on="idx")
    mat_city.append(matrix_city(c,"w",dest="prov"))
    p=evening.assign(origin=evening.i.map(orig_prov))
    mat_prov.append(matrix_prov(p,"w",index="origin",dest="prov"))

def ipf(W,B,keep):
    Y=W*keep
    O,D=Y.sum(axis=1),Y.sum(axis=0)
    distance=np.bincount(B.ravel(),Y.ravel())
    T=keep.astype(float)
    div=lambda a,b:np.divide(a,b,out=np.zeros_like(a),where=b>0)
    for iteration in range(600):
        T*=div(O,T.sum(axis=1))[:,None]
        T*=div(D,T.sum(axis=0))[None,:]
        T*=div(distance,np.bincount(B.ravel(),T.ravel(),minlength=len(distance)))[B]
        if max(np.max(np.abs(div(T.sum(axis=1)-O,O))),np.max(np.abs(div(T.sum(axis=0)-D,D))))<1e-9:
            return T
    raise RuntimeError("IPF")

def standardize(x,keep):
    return (x-x[keep].mean())/x[keep].std()

def residual(V,yhat,variant):
    if variant=="pearson":
        return (V-yhat)/np.maximum(np.sqrt(yhat),.01)
    if variant=="log":
        return np.log1p(V)-np.log1p(yhat)
    if variant=="asinh":
        return np.arcsinh(V)-np.arcsinh(yhat)
    return V-yhat

def run():
    for label,R,flows,origin,xy in [
        ("city",matrix_city(roots,"born"),mat_city,city_nodes.prov.to_numpy(),city_nodes[["x","y"]]),
        ("province",matrix_prov(proot,"persons"),mat_prov,province_order.to_numpy(),center)]:
        keep=origin[:,None]!=province_order.to_numpy()[None,:]
        dist=cdist(xy,center)/1000
        b=np.digitize(dist,np.quantile(dist[keep],np.linspace(0,1,21 if label=='city' else 16)[1:-1]))
        expected=[ipf(mat,b,keep) for mat in [R,*flows]]
        def decompose(X):
            X=X.copy()
            if label=="city":
                for j in np.unique(origin):
                    m=origin==j
                    X[m]-=X[m].mean(axis=0)
            return np.where(keep,X-X[keep].mean(),0.)
        for variant in ["log"]:
            X,Z=[decompose(residual(m,t,variant)) for m,t in zip([R,*flows],expected)][0],[
                decompose(residual(m,t,variant)) for m,t in zip(flows,expected[1:])]
            r=lambda a,b: np.corrcoef(a[keep],b[keep])[0,1]
            rng=np.random.default_rng(33913)
            permutations=np.repeat(np.arange(len(X))[None,:],2000,axis=0)
            group_ids=np.unique(origin)
            for row in permutations:
                if label=="province":
                    rng.shuffle(row)
                else:
                    for code in group_ids:
                        section=np.flatnonzero(origin==code)
                        row[section]=rng.permutation(section)
            shuffled=X[permutations] if label=="city" else X[permutations[:,:,None],permutations[:,None,:]]
            def p_diff(a):
                obs=(a*X).sum()
                draws=np.einsum("bij,ij->b",shuffled,a)
                return (1+(draws>=obs).sum())/(len(draws)+1)
            def partial(y):
                control=Z[0]
                norm=np.square(control).sum()
                return r(y-control*(y*control).sum()/norm,X-control*(X*control).sum()/norm)
            print(label,variant,[round(r(z,X),5) for z in Z],
                  "diff",[round(r(z-Z[0],X),5) for z in Z[1:]],
                  "partial",[round(partial(z),5) for z in Z[1:]],
                  "p_diff",[p_diff(z-Z[0]) for z in Z[1:]],flush=True)

if __name__=="__main__":
    run()

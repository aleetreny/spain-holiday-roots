from pathlib import Path
import igraph as ig
import numpy as np
import pandas as pd
import random
from scipy.stats import spearmanr

audit=Path(__file__).resolve().parent
c=pd.read_parquet(audit/"calendar.parquet").pivot(index="date",columns="kind",values="dist100")
for basis,t in [("file",c.index),("evening",c.index-pd.Timedelta(days=1))]:
    print(basis,len(c), sum(t.month==8),sum((t.month==8)&(c.coast>c.inland)),sum((t.month!=8)&(c.inland>c.coast)))
print(c.loc[["2024-08-01","2024-08-02","2024-08-15","2024-08-31"]])
bins=np.unique(np.round(np.geomspace(1,2436,30)).astype(int))
ig.set_random_number_generator(random.Random(770022))
for label in ("Feb","Aug","Nov"):
    g=ig.Graph(2435,np.loadtxt(audit/f"und_{label}.txt",dtype=int))
    degree=np.array(g.degree())
    cuts=[(lo,hi) for lo,hi in zip(bins[:-1],bins[1:]) if ((degree>=lo)&(degree<hi)).sum()>=5]
    def summarize(graph):
        knn=np.array(graph.knn()[0])
        return np.array([knn[(degree>=lo)&(degree<hi)].mean() for lo,hi in cuts])
    real=summarize(g)
    rows=[]
    for _ in range(10):
        h=g.copy()
        h.rewire(n=10*g.ecount(),allowed_edge_types="simple")
        rows.append(summarize(h))
    rewired=np.mean(rows,axis=0)
    print("knn",label,len(cuts),"Spearman",spearmanr(np.arange(len(cuts)),real).statistic,
          "strict increases",sum(np.diff(real)>0),"median fractional difference",
          np.median(abs(real-rewired)/rewired))

v=pd.read_csv(audit/"r1_ancova.csv")
print(v[v.date.isin(["2024-08-15","2024-11-01"])].to_string(index=False))
v=pd.read_csv(audit/"r1_level.csv")
for name,d in v.groupby("scope"):
    t=pd.to_datetime(d.date)
    m=t.dt.month
    dd=t.dt.day
    z=(m==2)|((m==3)&(dd<=22))|((m==10)&(~dd.isin([11,12,13])))|((m==11)&(dd>=4))|((m==12)&(dd<=20)&~dd.isin(range(5,10)))
    z=z&t.dt.year.eq(2024)
    for k in ["log_share","net_rate"]:
        X=np.column_stack((np.ones(len(d)),np.log(d.volume)))
        B=np.linalg.lstsq(X[z],d.loc[z,k],rcond=None)[0]
        resid=d[k]-X@B
        out=resid/resid[z].std(ddof=2)
        a=d.loc[m.eq(8)&t.dt.year.eq(2024)].volume>d.loc[z].volume.max()
        important=list(zip(d.date,out))
        keep={"2024-08-15","2024-11-01","2024-12-25"}
        print("volume",name,k,z.sum(),1-resid[z].var()/d.loc[z,k].var(),[(date,r) for date,r in important if date in keep],
               "extrapolated Aug",a.sum(),len(a),
               "August 10–18 residuals",np.percentile(out[m.eq(8)&dd.between(10,18)&t.dt.year.eq(2024)],[0,50,100]))

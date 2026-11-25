# Impact of the per-pair 'sum/count(DISTINCT date)' window-mean bug on p72 PPML and p93 neighbouring-city QAP.
import duckdb, numpy as np, pandas as pd, pyfixest as pf, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0,'explore'); from common import nodes
con=duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"
B=pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv',dtype={'ine_code':str,'prov_birth':str})
B=B[(B.year==2024)&B.prov_birth.str.match(r'^\d\d$')]
N=nodes().set_index('id'); cities=sorted(set(B.ine_code)&set(N.index)); cl=",".join(f"'{c}'" for c in cities)
M=N[~N.foreign].copy(); M['w']=M.population.fillna(0)+1
pc=M.groupby('prov').apply(lambda g: pd.Series({'x':np.average(g.x,weights=g.w),'y':np.average(g.y,weights=g.w)}))
grid=pd.MultiIndex.from_product([cities,sorted(pc.index)],names=['c','p']).to_frame(index=False)
grid=grid[grid.p.isin([f'{i:02d}' for i in range(1,51)])]; grid['cprov']=grid.c.str[:2]; grid=grid[grid.cprov!=grid.p]
grid['km']=np.hypot(N.loc[grid.c,'x'].values-pc.loc[grid.p,'x'].values,N.loc[grid.c,'y'].values-pc.loc[grid.p,'y'].values)/1000
grid['dbin']=np.digitize(grid.km,np.r_[0:100:10,100:400:20,400:1200:50])
grid['lborn']=np.log1p(B.set_index(['ine_code','prov_birth']).persons.reindex(list(zip(grid.c,grid.p))).fillna(0).values)
grid['pairFE']=grid.cprov+'_'+grid.p
sl={'Feb weekday':("date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",12),
    'ordinary Sat':("date IN ('2024-02-10','2024-02-17','2024-02-24','2024-10-19','2024-11-16','2024-11-23')",6),
    'mid-Aug':("date BETWEEN '2024-08-10' AND '2024-08-18'",9)}
print('city count',len(cities))
for lab,(cond,nd) in sl.items():
    for mode,den in [('buggy','count(DISTINCT date)'),('correct',str(nd))]:
        f=con.execute(f"SELECT r c, substr(s,1,2) p, sum(p)/{den} w FROM {q} WHERE {cond} AND r IN ({cl}) AND r<>s GROUP BY 1,2").df().set_index(['c','p']).w
        g=grid.copy(); g['w']=f.reindex(list(zip(g.c,g.p))).fillna(0).values
        t1=pf.fepois("w ~ lborn | c + p + dbin",data=g,vcov={'CRV1':'c'}).tidy(); t2=pf.fepois("w ~ lborn | c + pairFE",data=g,vcov={'CRV1':'c'}).tidy()
        print(f"{lab:12s} {mode:8s} total {g.w.sum()/1e6:.3f}M distFE {t1.loc['lborn','Estimate']:.3f} pairFE {t2.loc['lborn','Estimate']:.3f} ({t2.loc['lborn','Std. Error']:.3f})",flush=True)

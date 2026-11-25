# Zone-level (city x destination zone) PPML: impact of window-mean bug on the 'Feb weekday 0.29 vs holiday 0.5' contrast (p84 Part B)
import duckdb, numpy as np, pandas as pd, pyfixest as pf, warnings; warnings.filterwarnings('ignore')
tt=pd.read_parquet('explore/out_p84_city_zone_tt.parquet'); print(tt.columns.tolist(), len(tt), tt.c.nunique(), tt.z.nunique())
B=pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv',dtype={'ine_code':str,'prov_birth':str})
B=B[(B.year==2024)&B.prov_birth.str.match(r'^\d\d$')]; Bm=B.set_index(['ine_code','prov_birth']).persons
cz=tt[['c','z','eucl_km','tt_adj']].copy(); cz['p']=cz.z.str[:2]; cz['cprov']=cz.c.str[:2]; cz=cz[cz.cprov!=cz.p]
cz['lborn']=np.log1p(Bm.reindex(list(zip(cz.c,cz.p))).fillna(0).values); cz['pairFE']=cz.cprov+'_'+cz.p
cz['ttbin']=np.digitize(cz.tt_adj,np.r_[0:120:10,120:360:15,360:1000:30])
con=duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"; cl=",".join(f"'{c}'" for c in cz.c.unique())
for lab,cond,nd in [('Feb weekday',"date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",12),('Feb 14 single night',"date='2024-02-14'",1),('mid-Aug',"date BETWEEN '2024-08-10' AND '2024-08-18'",9),('Aug 15 single',"date='2024-08-15'",1)]:
    for mode,den in ([('buggy','count(DISTINCT date)'),('correct',str(nd))] if nd>1 else [('single','1')]):
        f=con.execute(f"SELECT r c, s z, sum(p)/{den} w FROM {q} WHERE {cond} AND r IN ({cl}) AND r<>s GROUP BY 1,2").df().set_index(['c','z']).w
        g=cz.copy(); g['w']=f.reindex(list(zip(g.c,g.z))).fillna(0).values
        t=pf.fepois("w ~ lborn | c + z + ttbin",data=g,vcov={'CRV1':'c'}).tidy(); t2=pf.fepois("w ~ lborn | c + z + pairFE + ttbin",data=g,vcov={'CRV1':'c'}).tidy()
        print(f"{lab:20s} {mode:8s} nz={int((g.w>0).sum())} c+z+ttbin {t.loc['lborn','Estimate']:.3f} ({t.loc['lborn','Std. Error']:.3f})  c+z+pairFE+ttbin {t2.loc['lborn','Estimate']:.3f} ({t2.loc['lborn','Std. Error']:.3f})",flush=True)

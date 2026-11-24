# Summer territories: bipartite city -> destination overnight network, communities under (a) Barber/BiCM-like null and (b) spatial gravity null.
# Also: what explains colonies beyond distance? same province / same region / birthplace roots (PPML).
import pandas as pd, numpy as np, duckdb, sys, pyfixest as pf
sys.path.insert(0,'explore'); from common import nodes, PROV2CCAA
from sklearn.metrics import normalized_mutual_info_score as nmi
N = nodes().set_index('id'); M = N[(~N.foreign)&(~N.island)&(~N.prov.isin(['51','52']))]
nm = pd.read_csv('raw_data_cache/v2/zonificacion/zonificacion_municipios/nombres_municipios.csv', sep='|', dtype=str).set_index('ID')['name']
name = lambda z: nm.get(z.replace('_AM',''), z)
con = duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"
NO = 100
origins = M.sort_values('population', ascending=False).head(NO).index.tolist()
print('origins cover', round(100*M.loc[origins,'population'].sum()/M.population.sum(),1), '% of mainland population')
dests = M.index.tolist()
R = pd.read_csv('_local/external/ine_birth_prov/birth_x_residence_2021_long.csv', dtype={'birth':str,'res':str}).set_index(['birth','res']).n
full = pd.MultiIndex.from_product([origins, dests], names=['o','d']).to_frame(index=False); full = full[full.o!=full.d].reset_index(drop=True)
dx = M.loc[full.o,['x','y']].values - M.loc[full.d,['x','y']].values
full['km'] = np.hypot(dx[:,0],dx[:,1])/1000
full['dbin'] = np.digitize(full.km, np.r_[0:30:2, 30:100:5, 100:300:10, 300:1200:25])
po, pd_ = full.o.str[:2], full.d.str[:2]
full['same_prov'] = (po.values==pd_.values).astype(int)
full['same_ccaa'] = (po.map(PROV2CCAA).values==pd_.map(PROV2CCAA).values).astype(int)
full['lroots'] = np.log1p(R.reindex(list(zip(pd_, po))).fillna(0).values)   # born in destination province, living in origin province
full.loc[full.same_prov==1,'lroots'] = 0
def ipf(df, it=150):
    E = np.ones(len(df)); w = df.w.values
    keys = [df.o.values, df.d.values, df.dbin.values]; tgt = [pd.Series(w).groupby(k).sum() for k in keys]
    for _ in range(it):
        for k,t in zip(keys,tgt):
            cur = pd.Series(E).groupby(k).sum(); E = E*(t/cur).fillna(0).reindex(k).values
    return E
def louvain(Bm, seed=0, passes=30):
    n=Bm.shape[0]; r=np.random.default_rng(seed); Bc=Bm.copy(); memb=np.arange(n)
    while True:
        nn=Bc.shape[0]; comm=np.arange(nn); improved=True; it=0
        while improved and it<passes:
            improved=False; it+=1
            for i in r.permutation(nn):
                ci=comm[i]; s=np.bincount(comm,weights=Bc[i],minlength=nn); s[ci]-=Bc[i,i]
                b=np.argmax(s)
                if s[b]>s[ci]+1e-12 and b!=ci: comm[i]=b; improved=True
        u,comm=np.unique(comm,return_inverse=True)
        if len(u)==nn: break
        Mm=np.zeros((nn,len(u))); Mm[np.arange(nn),comm]=1; Bc=Mm.T@Bc@Mm; memb=comm[memb]
    return memb
slices = {'mid-Aug': "date BETWEEN '2024-08-10' AND '2024-08-18'", 'All Saints': "date='2024-11-01'",
          'Feb weekday': "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)"}
ol = ",".join(f"'{x}'" for x in origins)
out = {}
for lab, cond in slices.items():
    f = con.execute(f"SELECT r o, s d, sum(p)/count(DISTINCT date) w FROM {q} WHERE {cond} AND r IN ({ol}) AND r<>s GROUP BY 1,2").df().set_index(['o','d']).w
    df = full.copy(); df['w'] = f.reindex(list(zip(df.o, df.d))).fillna(0).values
    df['E'] = ipf(df)
    # PPML: what explains excess beyond distance?
    m = pf.fepois("w ~ same_prov + same_ccaa + lroots | o + d + dbin", data=df, vcov={'CRV1':'o'})
    t = m.tidy()
    print(f"\n=== {lab}: PPML (o,d,dist-bin FE): same_prov x{np.exp(t.loc['same_prov','Estimate']):.2f}  same_ccaa x{np.exp(t.loc['same_ccaa','Estimate']):.2f}"
          f"  roots elasticity {t.loc['lroots','Estimate']:.3f} (se {t.loc['lroots','Std. Error']:.3f})")
    # bipartite modularity matrices
    oi = {o:i for i,o in enumerate(origins)}; di = {d:i for i,d in enumerate(dests)}
    no, nd_ = len(origins), len(dests); W = df.w.sum()
    O = np.zeros((no,nd_)); E = np.zeros((no,nd_))
    O[df.o.map(oi).values, df.d.map(di).values] = df.w.values; E[df.o.map(oi).values, df.d.map(di).values] = df.E.values
    keep = O.sum(0)>0; O, E = O[:,keep], E[:,keep]; dl = np.array(dests)[keep]
    ko, kd = O.sum(1), O.sum(0); Ebar = np.outer(ko,kd)/W
    res = {}
    for nul, Ex in [('Barber (no space)', Ebar), ('spatial gravity', E)]:
        Bb = O - Ex; n = no+len(dl); B = np.zeros((n,n)); B[:no,no:] = Bb; B[no:,:no] = Bb.T
        best=None
        for s in range(4):
            mb = louvain(B, seed=s); Q = sum(Bb[mb[:no][:,None]==mb[no:][None,:]])/W
            if best is None or Q>best[0]: best=(Q,mb)
        Q, mb = best
        co, cd = mb[:no], mb[no:]
        ncom = len(set(co))
        # how much of each origin's flow stays in its own territory
        inside = sum(O[i, cd==co[i]].sum() for i in range(no))/W
        dprov = pd.Series(dl).str[:2]; nmi_ccaa = nmi(dprov.map(PROV2CCAA), cd); nmi_prov = nmi(dprov, cd)
        print(f"  [{nul}] Q_B={Q:.3f}  territories(with cities)={ncom}  flow inside own territory {100*inside:.0f}%  NMI(dest comm, CCAA)={nmi_ccaa:.2f} prov={nmi_prov:.2f}")
        groups = pd.Series(origins).groupby(co).apply(lambda s: [name(x) for x in s])
        for c, cities in groups.items():
            if len(cities)>=1:
                top = pd.Series(O[:, cd==c].sum(0), index=dl[cd==c]).sort_values(ascending=False).head(4)
                print(f"     T{c}: {', '.join(cities[:8])}{' +%d'%(len(cities)-8) if len(cities)>8 else ''}  || dest n={int((cd==c).sum())}: {', '.join(name(z) for z in top.index)}")
        res[nul] = dict(Q=Q, co=co, cd=cd, dl=dl)
    out[lab] = res
pd.to_pickle(out, 'explore/out_summer_territories.pkl')

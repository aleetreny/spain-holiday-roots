# "Does each big city have its own private holiday territory?" Bipartite city -> destination overnight network vs gravity null.
import pandas as pd, numpy as np, duckdb, sys, pyfixest as pf
sys.path.insert(0,'explore'); from common import nodes
N = nodes().set_index('id'); M = N[(~N.foreign)&(~N.island)&(~N.prov.isin(['51','52']))]
nm = pd.read_csv('raw_data_cache/v2/zonificacion/zonificacion_municipios/nombres_municipios.csv', sep='|', dtype=str).set_index('ID')['name']
con = duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"
origins = M.sort_values('population', ascending=False).head(60).index.tolist()
def flows(cond):
    o = ",".join(f"'{x}'" for x in origins)
    df = con.execute(f"SELECT r o, s d, sum(p)/count(DISTINCT date) w FROM {q} WHERE {cond} AND r IN ({o}) AND r<>s GROUP BY 1,2").df()
    return df[df.d.isin(M.index)]
slices = {'mid-Aug': "date BETWEEN '2024-08-10' AND '2024-08-18'", 'All Saints': "date='2024-11-01'",
          'Feb weekday': "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)"}
dests = M.index.tolist()
full = pd.MultiIndex.from_product([origins, dests], names=['o','d']).to_frame(index=False)
full = full[full.o!=full.d]
dx = M.loc[full.o,['x','y']].values - M.loc[full.d,['x','y']].values
full['km'] = np.hypot(dx[:,0],dx[:,1])/1000
edges = np.r_[0:30:2, 30:100:5, 100:300:10, 300:1200:25]
full['dbin'] = np.digitize(full.km, edges)
res = {}
def ipf(df, it=200):
    # doubly-constrained gravity with binned deterrence = PPML with o, d, dbin fixed effects (ML), fitted by IPF
    E = np.ones(len(df)); w = df.w.values
    keys = [df.o.values, df.d.values, df.dbin.values]
    tgt = [pd.Series(w).groupby(k).sum() for k in keys]
    for _ in range(it):
        for k,t in zip(keys,tgt):
            cur = pd.Series(E).groupby(k).sum()
            fac = (t/cur).fillna(0)
            E = E*fac.reindex(k).values
    return E
for lab, cond in slices.items():
    f = flows(cond).set_index(['o','d']).w
    df = full.copy(); df['w'] = f.reindex(list(zip(df.o, df.d))).fillna(0).values
    df['E'] = ipf(df)
    # destination-level concentration: share of top origin, observed vs expected
    g = df.groupby('d')
    tot = g.w.transform('sum'); totE = g.E.transform('sum')
    df['sh'] = df.w/tot; df['shE'] = df.E/totE
    hh = df.groupby('d').apply(lambda x: pd.Series({'W':x.w.sum(),'HHI':(x.sh**2).sum(),'HHI_E':(x.shE**2).sum(),
                                                    'top':x.o.values[x.w.values.argmax()], 'topsh':x.sh.max(),
                                                    'topsh_E': x.shE.values[x.w.values.argmax()]}))
    hh = hh[hh.W>=200]
    ratio = (hh.HHI/hh.HHI_E)
    print(f"\n=== {lab}: destinations with >=200 visitors from the 60 big cities: {len(hh)}  total {df.w.sum()/1e6:.2f}M")
    print(f"  visitor-weighted HHI obs {np.average(hh.HHI,weights=hh.W):.3f} vs gravity-null {np.average(hh.HHI_E,weights=hh.W):.3f}"
          f"  (median ratio {ratio.median():.2f}; share of destinations with ratio>1.5: {100*(ratio>1.5).mean():.0f}%)")
    print(f"  destinations where ONE city supplies >50% of big-city visitors: obs {100*(hh.topsh>0.5).mean():.0f}%  expected-under-null for same city {100*(hh.topsh_E>0.5).mean():.0f}%")
    # excess concentration of flow: share of flow in cells with O/E>2
    print(f"  share of flow in (city,dest) cells with obs/exp>2: {100*df.w[df.w/df.E>2].sum()/df.w.sum():.1f}%  ; >5: {100*df.w[df.w/df.E>5].sum()/df.w.sum():.1f}%")
    res[lab] = (df, hh)
    # strongest colonies beyond gravity (min 500 sleepers)
    c = df[(df.w>=500)].assign(ratio=lambda x: x.w/x.E).sort_values('ratio', ascending=False).head(12)
    print('  top colonies (city -> dest, obs, obs/exp, km):')
    for r in c.itertuples():
        print(f"    {nm.get(r.o.replace('_AM',''), r.o):22s} -> {nm.get(r.d.replace('_AM',''), r.d)+' ('+r.d[:2]+')':32s} {r.w:8.0f}  x{r.ratio:5.1f}  {r.km:5.0f} km")
pd.to_pickle(res, 'explore/out_colonies.pkl')

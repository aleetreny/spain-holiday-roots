# Holiday-excess overnight network vs ordinary away-nights: distance, size mixing, communities (contiguity)
import duckdb, pandas as pd, numpy as np, sys, igraph as ig
sys.path.insert(0,'explore'); from common import nodes
con = duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"
N = nodes().set_index('id'); N = N[~N.foreign]
# reference: ordinary Fridays (no holiday/puente) Feb, Mar(early), Oct, Nov
ref_fri = ['2024-02-02','2024-02-09','2024-02-16','2024-02-23','2024-03-08','2024-03-15','2024-10-18','2024-11-15','2024-11-22']
ref_wd  = "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)"
def mat(cond):
    return con.execute(f"SELECT r, s, sum(p)/count(DISTINCT date) p FROM {q} WHERE {cond} AND r<>s GROUP BY 1,2").df()
reff = mat("date IN (" + ",".join(f"'{d}'" for d in ref_fri) + ")").set_index(['r','s']).p
def excess(cond):
    h = mat(cond).set_index(['r','s']).p
    e = (h - reff.reindex(h.index).fillna(0)).clip(lower=0)
    return e[e>0].rename('w').reset_index()
X = N[['x','y','population','prov']]
def describe(df, lab):
    df = df[df.r.isin(X.index)&df.s.isin(X.index)].copy()
    dx = X.loc[df.r,['x','y']].values - X.loc[df.s,['x','y']].values
    df['km'] = np.hypot(dx[:,0],dx[:,1])/1000
    w = df.w
    pr, ps = np.log10(X.loc[df.r,'population'].values), np.log10(X.loc[df.s,'population'].values)
    # weighted Pearson of log sizes (size assortativity of flows)
    def wcorr(a,b,w):
        ma, mb = np.average(a,weights=w), np.average(b,weights=w)
        return np.average((a-ma)*(b-mb),weights=w)/np.sqrt(np.average((a-ma)**2,weights=w)*np.average((b-mb)**2,weights=w))
    print(f"{lab:14s} total {w.sum()/1e6:5.2f}M  edges {len(df):7d}  mean km {np.average(df.km,weights=w):6.1f}  share>100km {100*w[df.km>100].sum()/w.sum():5.1f}%  "
          f"share other prov {100*w[X.loc[df.r,'prov'].values!=X.loc[df.s,'prov'].values].sum()/w.sum():5.1f}%  size-corr(res,dest) {wcorr(pr,ps,w):+.2f}  "
          f"share into <5k zones {100*w[X.loc[df.s,'population'].values<5000].sum()/w.sum():5.1f}%")
    return df
base = mat(ref_wd).rename(columns={'p':'w'})
describe(base, 'ordinary wd')
nets = {}
for lab, cond in [('Good Fri',"date='2024-03-29'"),('Aug 15',"date='2024-08-15'"),('All Saints',"date='2024-11-01'"),('Christmas Eve',"date='2024-12-24'")]:
    nets[lab] = describe(excess(cond), lab)
# communities of the All Saints excess network (undirected, weighted), contiguity check
import geopandas as gpd
g = gpd.read_file('clean_data/v2/zones/municipios_mitma.gpkg').set_index('id')
for lab in ['All Saints','Aug 15']:
    df = nets[lab]; df = df[df.w>=5]
    u = df.assign(a=np.minimum(df.r,df.s), b=np.maximum(df.r,df.s)).groupby(['a','b']).w.sum().reset_index()
    G = ig.Graph.TupleList(u.itertuples(index=False), directed=False, weights=True)
    best=None
    for seed in range(5):
        cm = G.community_leiden(objective_function='modularity', weights='weight', n_iterations=-1, resolution=1)
        Q = G.modularity(cm, weights='weight')
        if best is None or Q>best[1]: best=(cm,Q)
    cm,Q = best
    memb = pd.Series(cm.membership, index=G.vs['name'])
    big = memb.value_counts(); big = big[big>=10].index
    # contiguity: number of connected spatial pieces per community (touching polygons)
    pieces=[]
    for c in big:
        z = memb[memb==c].index
        sub = g.loc[g.index.intersection(z)].geometry
        pieces.append(sub.buffer(200).union_all().geoms.__len__() if sub.buffer(200).union_all().geom_type=='MultiPolygon' else 1)
    from sklearn.metrics import normalized_mutual_info_score as nmi
    print(f"\n{lab}: N={G.vcount()} E={G.ecount()} Q={Q:.3f}, comms>=10 nodes: {len(big)}, NMI(prov)={nmi(memb.index.str[:2], memb.values):.2f}")
    print('  spatial pieces per community (1 = contiguous):', sorted(pieces, reverse=True)[:20])
    memb.to_csv(f'explore/out_comm_{lab.replace(" ","_")}.csv')

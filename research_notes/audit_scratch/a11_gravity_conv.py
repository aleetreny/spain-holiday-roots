# Convergence of p93's fit_dc_gravity_rect (copied verbatim) on the matrices it is used for.
import duckdb, numpy as np, pandas as pd
def fit_dc_gravity_rect(W, bin_mat, mask, max_iter=300):
    nr, nc = W.shape; a = np.ones(nr); b = np.ones(nc); n_bins = int(bin_mat.max()) + 1
    f = np.ones(n_bins); f[0] = 0.0; W_m = np.where(mask, W, 0.0); r_tot = W_m.sum(1); c_tot = W_m.sum(0)
    bm = np.where(mask, bin_mat, 0); W_bin = np.bincount(bm.ravel(), weights=W_m.ravel(), minlength=n_bins)
    for it in range(max_iter):
        F = np.where(mask, f[bm], 0.0)
        for _ in range(5):
            a = np.where(r_tot > 0, r_tot / np.maximum(F @ b, 1e-300), 0.0)
            b = np.where(c_tot > 0, c_tot / np.maximum(F.T @ a, 1e-300), 0.0)
        E = np.outer(a, b) * F
        E_bin = np.bincount(bm.ravel(), weights=E.ravel(), minlength=n_bins)
        fn = np.where(E_bin > 0, f * W_bin / np.maximum(E_bin, 1e-300), 0.0); fn[0] = 0.0
        if np.max(np.abs(fn - f) / np.maximum(f, 1e-12)) < 1e-6: break
        f = fn
    E = np.outer(a, b) * np.where(mask, f[bm], 0.0)
    rel = lambda x,y: np.abs(x-y).sum()/W_m.sum()
    return E, it+1, rel(E.sum(1),r_tot), rel(E.sum(0),c_tot), rel(np.bincount(bm.ravel(),weights=E.ravel(),minlength=n_bins),W_bin)
nd=pd.read_parquet('data/derived/nodes_municipios.parquet'); nd=nd[~nd.foreign&~nd.island&~nd.prov.isin(['51','52'])].reset_index(drop=True)
idx={k:i for i,k in enumerate(nd.id)}; provs=sorted(nd.prov.unique()); p2i={p:i for i,p in enumerate(provs)}; n=len(nd)
bm=pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv',dtype={'ine_code':str,'prov_birth':str})
cities=sorted(set(bm[bm.year==2024].ine_code)&set(idx)); c2i={c:i for i,c in enumerate(cities)}; cp=np.array([c[:2] for c in cities])
print('mainland cities matched', len(cities))
nd['w']=nd.population.fillna(0)+1
pc=nd.groupby('prov').apply(lambda g: pd.Series({'x':np.average(g.x,weights=g.w),'y':np.average(g.y,weights=g.w)}),include_groups=False).loc[provs]
cx=nd.set_index('id').loc[cities,'x'].values; cy=nd.set_index('id').loc[cities,'y'].values
D=np.hypot(cx[:,None]-pc.x.values[None,:],cy[:,None]-pc.y.values[None,:])/1000; mask=cp[:,None]!=np.array(provs)[None,:]
qe=np.quantile(D[mask],np.linspace(0,1,21)); qe[0]=-1; qe[-1]+=1; Bc=np.digitize(D,qe); Bc[~mask]=0
con=duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"; cl=",".join(f"'{c}'" for c in cities); ml=",".join(f"'{z}'" for z in nd.id)
df=con.execute(f"SELECT r c, substr(s,1,2) ps, sum(p)/9 w FROM {q} WHERE date BETWEEN '2024-08-10' AND '2024-08-18' AND r IN ({cl}) AND s IN ({ml}) AND substr(r,1,2)<>substr(s,1,2) GROUP BY 1,2").df()
M=np.zeros((len(cities),len(provs))); M[df.c.map(c2i),df.ps.map(p2i)]=df.w
E,it,er,ec,eb=fit_dc_gravity_rect(M,Bc,mask); print('city x prov Mid-Aug: iters',it,'rel margin err rows %.2e cols %.2e bins %.2e'%(er,ec,eb))
# city x zone (as in p93 Part 2, max_iter=200)
Dz=np.hypot(cx[:,None]-nd.x.values[None,:],cy[:,None]-nd.y.values[None,:])/1000
mz=np.ones((len(cities),n),bool)
for i,c in enumerate(cities): mz[i,idx[c]]=False
Bz=np.digitize(Dz,np.r_[0:30:2,30:100:5,100:300:10,300:1200:25]); Bz[~mz]=0
df=con.execute(f"SELECT r c, s d, sum(p)/9 w FROM {q} WHERE date BETWEEN '2024-08-10' AND '2024-08-18' AND r IN ({cl}) AND s IN ({ml}) AND r<>s GROUP BY 1,2").df()
O=np.zeros((len(cities),n)); O[df.c.map(c2i),df.d.map(idx)]=df.w
E,it,er,ec,eb=fit_dc_gravity_rect(O,Bz,mz,max_iter=200); print('city x zone Mid-Aug: iters',it,'rel margin err rows %.2e cols %.2e bins %.2e'%(er,ec,eb))
# overdispersion: Pearson chi2/df of O vs E on active dests
keep=O.sum(0)>=10; Ok=O[:,keep]; Ek=E[:,keep]; mm=Ek>0
print('Pearson dispersion (O-E)^2/E mean over cells: %.1f (Poisson=1)'%(((Ok[mm]-Ek[mm])**2/Ek[mm]).mean()))

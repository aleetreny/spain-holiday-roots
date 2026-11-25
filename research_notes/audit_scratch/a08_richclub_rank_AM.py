import duckdb, numpy as np, pandas as pd, igraph as ig
from scipy import sparse
con=duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"
nd=pd.read_parquet('data/derived/nodes_municipios.parquet')
nd=nd[~nd.foreign&~nd.island&~nd.prov.isin(['51','52'])&(nd.id.str.len()<=8)].reset_index(drop=True)
idx={k:i for i,k in enumerate(nd.id)}; ml=",".join(f"'{z}'" for z in nd.id); n=len(nd)
def build(cond):
    df=con.execute(f"SELECT r,s,sum(p) w FROM {q} WHERE {cond} AND r<>s AND r IN ({ml}) AND s IN ({ml}) GROUP BY 1,2").df()
    df=df[df.w>=10]; ri=df.r.map(idx).values; si=df.s.map(idx).values
    u=pd.DataFrame({'a':np.minimum(ri,si),'b':np.maximum(ri,si)}).drop_duplicates()
    return ig.Graph(n=n,edges=list(zip(u.a,u.b))), np.bincount(si,minlength=n)
def A_(G):
    e=np.array(G.get_edgelist()); return sparse.csr_matrix((np.ones(2*len(e)),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(n,n))
def phi_top(A,deg,N):
    o=np.argsort(-deg,kind='stable')[:N]; return A[o][:,o].nnz/(N*(N-1))
rng=np.random.default_rng(2)
res={}
for lab,cond in [('Feb14',"date='2024-02-14'"),('Aug15',"date='2024-08-15'")]:
    G,kin=build(cond); deg=np.array(G.degree()); A=A_(G)
    tops=[10,25,50,100,200]
    obs={N:phi_top(A,deg,N) for N in tops}; rew={N:[] for N in tops}
    for _ in range(10):
        Gr=G.copy(); Gr.rewire(n=10*G.ecount()); Ar=A_(Gr)
        for N in tops: rew[N].append(phi_top(Ar,deg,N))
    print(lab,'k at rank 50/100/200:',np.sort(deg)[::-1][[49,99,199]].tolist(),' phi_obs/phi_rew/rho by top-N:',{N:(round(obs[N],3),round(np.mean(rew[N]),3),round(obs[N]/np.mean(rew[N]),3)) for N in tops})
    res[lab]=(G,kin)
# _AM confound for k-core promotion
H=pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv',dtype={'ine_code':str})
rel=pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv',sep='|',dtype=str)[['municipio_ine','municipio_mitma']].drop_duplicates()
cols=['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960','pop_2021']
Z=H.merge(rel,left_on='ine_code',right_on='municipio_ine').groupby('municipio_mitma')[cols].sum(min_count=1)
loss=1-Z.pop_2021/Z[cols[:7]].max(axis=1); nmun=rel.groupby('municipio_mitma').size()
cf=pd.read_csv('explore/out_p86_coastal_flags.csv',dtype={'id':str}).set_index('id')
T=nd[['id','population']].copy(); T['cf']=res['Feb14'][0].coreness(); T['ca']=res['Aug15'][0].coreness(); T['ka']=res['Aug15'][1]; T['kf']=res['Feb14'][1]
T['loss']=T.id.map(loss); T['nmun']=T.id.map(nmun); T['AM']=T.id.str.endswith('_AM'); T['coast']=T.id.map(cf.coastal_true)
I=T[(T.coast==0)&(T.population>=200)&T.loss.notna()].copy()
S=I[I.population<6000].copy(); S['t']=pd.qcut(S.loss,3,labels=['low','mid','extreme'])
print(S.groupby('t',observed=True).agg(n=('id','size'),share_AM=('AM','mean'),mean_nmun=('nmun','mean'),med_nmun=('nmun','median'),kc_feb=('cf','mean'),kc_aug=('ca','mean'),kin_aug=('ka','mean')).round(2).to_string())
import statsmodels.formula.api as smf
I['lrank']=I.loss.rank(pct=True); I['lpop']=np.log(I.population); I['lnmun']=np.log(I.nmun); I['prov']=I.id.str[:2]
for f in ['ca ~ lrank + lpop + cf','ca ~ lrank + lpop + cf + lnmun','ca ~ lrank + lpop + cf + lnmun + C(prov)','np.log1p(ka) ~ lrank + lpop + np.log1p(kf) + lnmun + C(prov)']:
    m=smf.ols(f,data=I).fit(cov_type='HC1'); print(f,'-> loss_rank b=%.3f se=%.3f'%(m.params['lrank'],m.bse['lrank']), 'n',int(m.nobs))
# within single-municipality zones only
J=I[I.nmun==1]
m=smf.ols('ca ~ lrank + lpop + cf + C(prov)',data=J).fit(cov_type='HC1'); print('single-municipality zones only n=%d loss_rank b=%.3f se=%.3f'%(m.nobs,m.params['lrank'],m.bse['lrank']))

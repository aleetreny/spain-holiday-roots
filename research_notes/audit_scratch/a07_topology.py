# Re-run p92 core metrics: node set (2484 vs 2435), coastal flag (p73 buffer vs p86 coastal_true),
# window means (buggy per-pair day count vs correct), rewiring acceptance.
import duckdb, numpy as np, pandas as pd, igraph as ig, json
from scipy import sparse
con=duckdb.connect(); q="'data/derived/pernoct_mun_2024_6m.parquet'"
nd=pd.read_parquet('data/derived/nodes_municipios.parquet')
nd=nd[~nd.foreign&~nd.island&~nd.prov.isin(['51','52'])].reset_index(drop=True)
cf=pd.read_csv('explore/out_p86_coastal_flags.csv',dtype={'id':str}).set_index('id')
H=pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv',dtype={'ine_code':str})
rel=pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv',sep='|',dtype=str)[['municipio_ine','municipio_mitma']].drop_duplicates()
cols=['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960','pop_2021']
Z=H.merge(rel,left_on='ine_code',right_on='municipio_ine').groupby('municipio_mitma')[cols].sum(min_count=1)
loss=1-Z.pop_2021/Z[cols[:7]].max(axis=1)
rng=np.random.default_rng(1)
def build(nodes, cond, ndays=None):
    idx={k:i for i,k in enumerate(nodes.id)}; ml=",".join(f"'{z}'" for z in nodes.id)
    den = f"{ndays}" if ndays else "count(DISTINCT date)"
    df=con.execute(f"SELECT r,s,sum(p)/{den} w FROM {q} WHERE {cond} AND r<>s AND r IN ({ml}) AND s IN ({ml}) GROUP BY 1,2").df()
    df=df[df.w>=10]; ri=df.r.map(idx).values; si=df.s.map(idx).values
    a,b=np.minimum(ri,si),np.maximum(ri,si); u=pd.DataFrame({'a':a,'b':b}).drop_duplicates()
    G=ig.Graph(n=len(nodes),edges=list(zip(u.a,u.b)))
    kin=np.bincount(si,minlength=len(nodes))
    return G,kin,len(df)
def csr(G):
    e=np.array(G.get_edgelist()); n=G.vcount()
    return sparse.csr_matrix((np.ones(2*len(e)),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(n,n))
def phi(A,k):
    d=np.asarray(A.sum(1)).ravel(); m=d>k; nk=m.sum(); return A[m][:,m].nnz/(nk*(nk-1)) if nk>1 else np.nan
out={}
full=nd; real=nd[nd.id.str.len()<=8].reset_index(drop=True)
for lab,cond,nd_ in [('Feb14',"date='2024-02-14'",None),('Aug15',"date='2024-08-15'",None)]:
    for nset_name,nodes in [('2484',full),('2435',real)]:
        G,kin,Ed=build(nodes,cond); n=G.vcount(); E=G.ecount()
        res=dict(N=n,E_und=E,E_dir=Ed,density=2*E/(n*(n-1)),k_mean=float(np.mean(G.degree())),
                 isolated=int(sum(1 for d in G.degree() if d==0)),
                 C_avg_zero=G.transitivity_avglocal_undirected(mode='zero'),r=G.assortativity_degree(directed=False))
        if nset_name=='2435':
            orig=set(map(tuple,np.sort(np.array(G.get_edgelist()),1).tolist()))
            rr=[];ov=[];ph={k:[] for k in (50,100,200)}
            A=csr(G); ph0={k:phi(A,k) for k in ph}
            for rep in range(10):
                Gr=G.copy(); Gr.rewire(n=10*E)
                assert Gr.degree()==G.degree()
                new=set(map(tuple,np.sort(np.array(Gr.get_edgelist()),1).tolist()))
                ov.append(len(orig&new)/E); rr.append(Gr.assortativity_degree(directed=False))
                Ar=csr(Gr)
                for k in ph: ph[k].append(phi(Ar,k))
            res.update(r_rew_mean=float(np.mean(rr)),r_rew_sd=float(np.std(rr)),edge_overlap_after_rewire=float(np.mean(ov)),
                       rho={k:(ph0[k]/np.mean(v), float(np.std(np.array(v)/np.mean(v))*ph0[k]/np.mean(v))) for k,v in ph.items()},
                       kmax=int(max(G.degree())))
            # expected overlap of random graph with same degrees ~ sum k_i k_j/(2E) over edges -> report
            out[lab+'_graph']=(G,kin,nodes)
        out[f'{lab}_{nset_name}']=res
        print(lab,nset_name,{k:(round(v,4) if isinstance(v,float) else v) for k,v in res.items()},flush=True)
# window mean bug impact on network properties
for lab,cond,nd_days in [('MidAug_mean',"date BETWEEN '2024-08-10' AND '2024-08-18'",9),('Feb_wd_mean',"date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",12)]:
    for mode,nds in [('buggy',None),('correct',nd_days)]:
        G,kin,Ed=build(real,cond,nds); n=G.vcount(); E=G.ecount()
        print(lab,mode,'E_und',E,'density %.4f'%(2*E/(n*(n-1))),'r %.3f'%G.assortativity_degree(directed=False),'kcore max',max(G.coreness()),flush=True)
# k-core quintile table: buggy coastal vs true coastal
Gf,kf,nodes=out['Feb14_graph']; Ga,ka,_=out['Aug15_graph']
T=pd.DataFrame({'id':nodes.id,'pop':nodes.population.values,'cf':np.array(Gf.coreness()),'ca':np.array(Ga.coreness()),'kf':kf,'ka':ka})
T['loss']=T.id.map(loss); T=T.join(cf,on='id')
for flag in ['coastal_p73','coastal_true']:
    I=T[(T[flag]==0)&(T['pop']>=200)&T.loss.notna()].copy(); I['q']=pd.qcut(I.loss,5,labels=False)
    t=I.groupby('q').agg(n=('pop','size'),med_pop=('pop','median'),med_loss=('loss','median'),kc_feb=('cf','mean'),kc_aug=('ca','mean'),kin_feb=('kf','mean'),kin_aug=('ka','mean'))
    t['kc_x']=t.kc_aug/t.kc_feb; t['kin_x']=t.kin_aug/t.kin_feb; t['kc_diff']=t.kc_aug-t.kc_feb
    print('\n',flag); print(t.round(2).to_string())
    # size-matched: pop<6000 terciles
    S=I[I['pop']<6000].copy(); S['t']=pd.qcut(S.loss,3,labels=False)
    print(S.groupby('t').agg(n=('pop','size'),med_pop=('pop','median'),mean_pop=('pop','mean'),med_loss=('loss','median'),kc_feb=('cf','mean'),kc_aug=('ca','mean'),kin_feb=('kf','mean'),kin_aug=('ka','mean')).round(2).to_string())
    # regression: coreness_aug ~ loss_rank + log pop + coreness_feb (OLS)
    import statsmodels.formula.api as smf
    I['lpop']=np.log(I['pop']); I['lrank']=I.loss.rank(pct=True)
    f=smf.ols('ca ~ lrank + lpop + cf + C(id.str[:2])',data=I.assign(**{'id.str[:2]':I.id.str[:2]})) if False else smf.ols('ca ~ lrank + lpop + cf',data=I).fit(cov_type='HC1')
    print('OLS kcore_aug ~ loss_rank + lpop + kcore_feb: b=%.2f se=%.2f'%(f.params['lrank'],f.bse['lrank']))

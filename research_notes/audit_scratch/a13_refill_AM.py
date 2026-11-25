import pandas as pd, numpy as np, pyfixest as pf
H=pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv',dtype={'ine_code':str})
rel=pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv',sep='|',dtype=str)[['municipio_ine','municipio_mitma']].drop_duplicates()
cols=['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960','pop_2021']
Z=H.merge(rel,left_on='ine_code',right_on='municipio_ine').groupby('municipio_mitma')[cols].sum(min_count=1)
Z['loss']=1-Z.pop_2021/Z[cols[:7]].max(axis=1); Z['nmun']=rel.groupby('municipio_mitma').size()
cf=pd.read_csv('explore/out_p86_coastal_flags.csv',dtype={'id':str}).set_index('id')
R=pd.read_csv('explore/out_refill.csv',dtype={'z':str}).set_index('z')
D=R.join(Z[['loss','nmun']],how='inner').join(cf,how='inner'); D['prov']=D.index.str[:2]
D=D[(D.resid>=200)&~D.prov.isin(['07','35','38','51','52'])].copy()
D['lmult_aug']=np.log(D.mult_aug); D['lres']=np.log(D.resid); D['loss_rank']=D.loss.rank(pct=True); D['lnmun']=np.log(D.nmun)
I=D[D.coastal_true==0]
for f in ["lmult_aug ~ loss_rank + lres | prov","lmult_aug ~ loss_rank + lres + lnmun | prov"]:
    t=pf.feols(f,data=I,vcov={'CRV1':'prov'}).tidy(); print(f,'n',len(I),'b=%.3f se=%.3f'%(t.loc['loss_rank','Estimate'],t.loc['loss_rank','Std. Error']))
J=I[I.nmun==1]; t=pf.feols("lmult_aug ~ loss_rank + lres | prov",data=J,vcov={'CRV1':'prov'}).tidy(); print('single-municipality zones n',len(J),'b=%.3f se=%.3f'%(t.loc['loss_rank','Estimate'],t.loc['loss_rank','Std. Error']))
J=J.copy(); J['q']=pd.qcut(J.loss,5,labels=False); print(J.groupby('q').agg(n=('resid','size'),med_loss=('loss','median'),mult_aug=('mult_aug','median'),mult_nov1=('mult_nov1','median'),mult_feb=('mult_feb','median')).round(2).to_string())

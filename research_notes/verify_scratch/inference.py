import runpy
from pathlib import Path
from itertools import product
import sys

audit=Path(__file__).resolve().parent
r=runpy.run_path(str(audit/"roots.py"))
globals().update({k:v for k,v in r.items() if not k.startswith("__")})

def prep(f, p, terms, date):
    y=p[date].to_numpy()
    active=pd.Series(y,index=f.index).groupby(f.origin+f.prov).transform("sum").gt(0).to_numpy()
    ff=f[active]
    return ff,design(ff,terms),y[active]

def score_test(frame0, frame1, k, beta, webb):
    f0,x0,y0=frame0
    f1,x1,y1=frame1
    n0,n1=len(y0),len(y1)
    all_y=np.r_[y0,y1]
    tgt=np.r_[np.zeros(n0),x1[:,k].toarray().ravel()]
    base=sparse.block_diag([x0[:,np.arange(x0.shape[1])!=k],x1[:,np.arange(x1.shape[1])!=k]],format="csr")
    shared=sparse.vstack([x0[:,k],x1[:,k]],format="csr")
    X=sparse.hstack([base,shared],format="csr")
    _,mu,_=solve(X,all_y,offset=beta*tgt)
    weighted=X.T@X.multiply(mu[:,None])
    partial=tgt-X@spsolve(weighted.tocsc(),X.T@(tgt*mu))
    cluster=pd.factorize(pd.concat([f0.origin,f1.origin]),sort=True)[0]
    score=np.bincount(cluster,partial*(all_y-mu))
    weights=webb[:,:len(score)]
    stat=np.abs(score.sum())/np.linalg.norm(score)
    null=np.abs(weights@score)/np.sqrt((weights**2)@(score**2))
    return (1+np.count_nonzero(null>=stat))/(len(null)+1)

def run_wcr(name, sister_only=False, model_only=None):
    r2=name=="R2"
    f=controls(r2_data() if r2 else census_data(),"pop")
    scopes=["all","inland","coast"]
    panels=panel(f,scopes)
    fits=pd.read_csv(audit/f"{name}_pop.csv",dtype={"drop":str}).query("drop.isna()",engine="python")
    main={"R2":["ln_roots"]} if r2 else {"cohort":["pre1981","recent"],"birth":["village","urban"],
           "age":["pre65","post65","pre_u65","recent_u65"],"horserace":["pre1981","recent","born65"],
           "youngvillage":["village_u65","village65","urban_u65","urban65"]}
    webb=np.random.default_rng(573162).choice(np.sqrt([.5,1.,1.5]),size=(9999,47))
    webb*=np.random.default_rng(73911).choice([-1,1],size=webb.shape)
    out=[]
    file=audit/f"{name}_wcr{'_sisters' if sister_only else ''}{'_'+model_only if model_only else ''}.csv"
    for scope,p in panels.items():
        for label,terms in main.items():
            if model_only and label!=model_only:
                continue
            augmented=label not in ("R2","cohort","birth")
            if augmented and scope=="coast":
                continue
            if sister_only:
                if augmented or scope!="all":
                    continue
                pairs=list(product(dates_aug,dates_feb))
            else:
                pairs=[(night,days[0]) for night in (days[1:2] if augmented else days[1:])]
                if r2:
                    pairs += [(night,"oct_mean") for night in days[1:]]
            prep_fits={date:prep(f,p,terms,date) for date in set(sum(([h,o] for h,o in pairs),[]))}
            fitted=fits[(fits.model==label)&(fits.scope==scope)].set_index("date")
            for h,o in pairs:
                for k,term in enumerate(terms):
                    estimate=fitted.loc[h,term]-fitted.loc[o,term]
                    test=lambda b:score_test(prep_fits[o],prep_fits[h],k,b,webb)
                    pval=test(0)
                    lo,hi=np.nan,np.nan
                    if not sister_only and h==days[1] and o==days[0]:
                        test95=lambda b:test(b)-.05
                        lo=brentq(test95,estimate-.45,estimate,xtol=.00015)
                        hi=brentq(test95,estimate,estimate+.45,xtol=.00015)
                    out.append(dict(scope=scope,model=label,term=term,holiday=h,ordinary=o,
                                    delta=estimate,p=pval,lo=lo,hi=hi))
                pd.DataFrame(out).to_csv(file,index=False)
            print(name,scope,label,len(out),flush=True)

if __name__=="__main__":
    run_wcr(sys.argv[1], "sisters" in sys.argv, sys.argv[3] if len(sys.argv)>3 else None)

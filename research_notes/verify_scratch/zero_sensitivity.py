from pathlib import Path
import runpy
import numpy as np
import pandas as pd

audit=Path(__file__).resolve().parent
r=runpy.run_path(str(audit/"roots.py"))
transforms={
    "n0.1":lambda n,N:np.log(n+.1),
    "n0.5":lambda n,N:np.log(n+.5),
    "n1.0":lambda n,N:np.log(n+1),
    "share1e-4":lambda n,N:np.log(n/N+1e-4),
    "share1e-5":lambda n,N:np.log(n/N+1e-5),
    "asinh1000share":lambda n,N:np.arcsinh(1000*n/N),
}
out=[]
for spec,fn in transforms.items():
    f=r["controls"](r["census_data"](fn),"pop")
    panel=r["panel"](f,["all"])["all"]
    for terms in (["pre1981","recent"],["village","urban"]):
        fits=[]
        for night in r["days"][:2]:
            y=panel[night].to_numpy()
            keep=pd.Series(y,index=f.index).groupby(f.origin+f.prov).transform("sum").gt(0)
            beta,_,_=r["solve"](r["design"](f[keep],terms),y[keep])
            fits.append(beta[:len(terms)])
        for term,slope in zip(terms,fits[1]-fits[0]):
            out.append(dict(spec=spec,term=term,delta=slope))
    pd.DataFrame(out).to_csv(audit/"v1_zeros.csv",index=False)

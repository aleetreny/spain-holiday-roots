from pathlib import Path
import runpy

scratch=Path(__file__).resolve().parent
methods=runpy.run_path(str(scratch/"inference.py"))
np=methods["np"]
pd=methods["pd"]
f=methods["controls"](methods["census_data"](),"pop")
panels=methods["panel"](f,["inland","coast"])
generator=np.random.default_rng(19816040)
webb=generator.choice(np.r_[-np.sqrt([.5,1,1.5]), np.sqrt([.5,1,1.5])],size=(9999,44))
o,h="2024-02-14","2024-08-15"
out=[]
for scope,p in panels.items():
    for terms,primary in [(["pre1981","recent"],[0,1]),(["village","urban"],[0])]:
        ff=[methods["prep"](f,p,terms,t) for t in (o,h)]
        beta=[methods["solve"](X,y)[0] for _,X,y in ff]
        for k in primary:
            d=beta[1][k]-beta[0][k]
            test=lambda null:methods["score_test"](*ff,k,null,webb)
            test95=lambda null:test(null)-0.05
            bounds=[methods["brentq"](test95,min(a,d),max(a,d),xtol=.0002) for a in (d-.35,d+.35)]
            out.append([scope,terms[k],d,test(0),*bounds])
            df=pd.DataFrame(out,columns=["scope","term","delta","p","lo","hi"])
            df.to_csv(scratch/"figure3_inference.csv",index=False)

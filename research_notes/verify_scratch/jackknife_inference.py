from pathlib import Path
import runpy
import sys

scratch=Path(__file__).resolve().parent
m=runpy.run_path(str(scratch/"inference.py"))
np=m["np"]
pd=m["pd"]
f=m["controls"](m["census_data"](),"pop")
full=m["panel"](f,["all"])["all"]
fit=pd.read_csv(scratch/"V1_pop.csv",dtype={"drop":str})
label=sys.argv[1]
terms,key=(["pre1981","recent"],"pre1981") if label=="cohort" else (["village","urban"],"village")
provinces=sorted(f.origin.unique())
generator=np.random.default_rng(931641)
webb=generator.choice(np.r_[-np.sqrt([.5,1,1.5]),np.sqrt([.5,1,1.5])],size=(9999,len(provinces)))
out=[]
for province in provinces:
    keep=f.origin.ne(province)
    data=[m["prep"](f[keep].reset_index(drop=True),full.loc[keep.to_numpy()],terms,t)
          for t in ["2024-02-14","2024-08-15"]]
    block=fit[(fit["drop"]==province)&(fit.scope=="all")&(fit.model==label)].set_index("date")
    delta=block.loc["2024-08-15",key]-block.loc["2024-02-14",key]
    p=m["score_test"](*data,0,0,webb)
    out.append([province,key,delta,p])
    pd.DataFrame(out,columns=["drop","term","delta","p"]).to_csv(scratch/f"v1_jack_{label}.csv",index=False)

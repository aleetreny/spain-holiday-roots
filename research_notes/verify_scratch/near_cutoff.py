import runpy
from pathlib import Path
import numpy as np

audit=Path(__file__).resolve().parent
i=runpy.run_path(str(audit/"inference.py"))
f=i["controls"](i["census_data"](),"pop")
p=i["panel"](f,["all"])["all"]
terms=["pre1981","recent"]
ordinary=i["prep"](f,p,terms,"2024-02-28")
rng=np.random.default_rng(468551)
webb=rng.choice(np.r_[-np.sqrt([.5,1,1.5]),np.sqrt([.5,1,1.5])],size=(149999,44))
for day in ["2024-08-12","2024-08-15"]:
    holiday=i["prep"](f,p,terms,day)
    print(day,i["score_test"](ordinary,holiday,1,0,webb),flush=True)

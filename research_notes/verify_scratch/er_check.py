from pathlib import Path
import random
import numpy as np
import pandas as pd
import igraph as ig

scratch=Path(__file__).resolve().parent
ig.set_random_number_generator(random.Random(6245471))
samples=[]
for d in range(100):
    g=ig.Graph.Erdos_Renyi(n=2435,m=101996)
    samples.append(g.average_path_length())
    if (d+1)%10==0:
        pd.Series(samples,name="l").to_csv(scratch/"aug_er_check.csv",index=False)
        print(d+1,np.mean(samples),np.std(samples,ddof=1),flush=True)

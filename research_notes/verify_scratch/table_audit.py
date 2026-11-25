from pathlib import Path
import re
import numpy as np
import pandas as pd

scratch=Path(__file__).resolve().parent
real=pd.read_csv(scratch/"network_real.csv").set_index("label")
rand=pd.read_csv(scratch/"network_null.csv")
constraints=pd.read_csv(scratch/"constraints.csv")
tails=pd.read_csv(scratch/"tails.csv").set_index(["label","kind"])
history=pd.read_csv("results/b02_basic_properties.csv")
history=history.set_index(["date","backbone","metric","null"])
struct=pd.read_csv("results/b02_structural_constraints.csv").set_index(["date","backbone"])
labels=["Feb","Aug","Nov","AugEq"]

def avg(sub,key):
    return sub[key].mean()

def val(label, metric, null=None):
    if null is None:
        return real.loc[label,metric]
    z=rand.loc[rand.label.eq(label)&rand.model.eq(null)]
    if label=="Aug" and metric=="path" and null=="er":
        return pd.concat([z[metric],pd.read_csv(scratch/"aug_er_check.csv").l]).mean()
    return avg(z,metric)

def r_avg(label):
    sub=constraints[constraints.label.eq(label)]
    return sub[sub.model.eq("real")].R_avg.iloc[0], avg(sub[sub.model.eq("phi")],"R_avg")

def bounds(label):
    d=constraints[constraints.label.eq(label)]
    return d.loc[d.model.eq("min"),"r"].min(),d.loc[d.model.eq("max"),"r"].max()

def tail(label):
    out=[]
    for kind in ("in","out"):
        r=tails.loc[(label,kind)]
        out.append("log-normal" if r.p<.05 and r.R<0 else "tie")
    return " / ".join(out)

def calculate(label, key):
    row=real.loc[label]
    mappings={
        "giant":f"{row.giant:,}", "N":f"{row.N:,} ({row.giant:,})",
        "L":f"{row.L:,}", "density":f"{100*row.density:.2f}",
        "kavg":f"{row.k_mean:.1f}", "kmed":f"{row.k_median:.0f}",
        "max":f"{row.k_max:,} ({100*row.k_max/row.N:.0f} %)",
        "med":f"{row.med_out:.0f} / {row.med_in:.0f}",
        "retained":f"{100*row.retained:.2f}","weight":f"{row.full/1e6:.6f}",
        "reciprocity":f"{row.recip:.3f} ({val(label,'recip','rew'):.4f})",
        "r":f"{row.r:.4f} ({val(label,'r','rew'):.5f})",
        "phi":f"{row.phi50:.4f} ({val(label,'phi50','rew'):.4f})",
        "prefix":f"{row.clique_prefix:.0f}",
        "clique":f"{row.clique_prefix:.0f} ({val(label,'clique_prefix','rew'):.2f})",
        "path":f"{row.path:.4f}",
        "paths":f"{row.path:.4f} ({val(label,'path','rew'):.4f} / {val(label,'path','er'):.5f})",
        "pathstar":f"{np.log(row.giant)/np.log(2*row.L/row.giant):.4f}",
        "clustering":f"{row.C:.4f} ({val(label,'C','rew'):.5f} / {val(label,'C','er'):.5f})",
        "core":f"{row.core:.0f} ({val(label,'core','rew'):.2f} / {val(label,'core','er'):.2f})",
        "diam":f"{row.diam:.0f} ({val(label,'diam','rew'):.2f})",
        "touch":f"{100*row.touch121:.2f}",
        "tail":tail(label),
        "range":"[{:.6f}, {:.6f}]".format(*bounds(label)),
        "Ravg":"{:.5f} ({:.5f})".format(*r_avg(label)),
    }
    return mappings[key]

rows1={5:"giant",6:"L",7:"density",8:"kavg",9:"max",10:"med",11:"reciprocity",
       12:"r",13:"phi",14:"prefix",15:"path",16:"pathstar",17:"clustering"}
rowsA={6:"weight",7:"retained",8:"N",9:"L",10:"density",12:"kavg",13:"kmed",14:"max",15:"med",
       16:"tail",18:"paths",19:"pathstar",20:"diam",21:"clustering",23:"r",24:"range",
       25:"reciprocity",26:"Ravg",28:"phi",29:"clique",30:"touch",31:"core"}
refs={
    "reciprocity":[("recip","rew","reciprocity","rew_dir")],
    "r":[("r","rew","r_und","rew_und")],
    "phi":[("phi50","rew","phi_top50","rew_und")],
    "paths":[("path","rew","L","rew_und"),("path","er","L","er_und")],
    "clustering":[("C","rew","C_avg","rew_und"),("C","er","C_avg","er_und")],
    "core":[("core","rew","kcore_max","rew_und"),("core","er","kcore_max","er_und")],
    "diam":[("diam","rew","diameter","rew_und")],
    "clique":[("clique_prefix","rew","n_clique","rew_und")],
    "Ravg":[("R_avg","phi","rcpres_R_avg","phi")],
}
notes={
    "range":"Feasible witnesses under the fixed degree sequence and connected edge-bearing subgraph. Endpoints are heuristic estimates.",
    "giant":"Largest component. Other Table 1 denominators follow the full 2,435-node projection.",
    "N":"Includes isolates. Weak giant equals largest undirected component.",
    "L":"Unique unordered pairs.",
    "density":"Denominator is all 2,435 choose 2 possible pairs.",
    "max":"The integer and percentage agree. Table 1 uses a slash in place of parentheses.",
    "med":"Both medians are from the directed backbone.",
    "tail":"Positive in/out degrees. p ≥ 0.05 is a failure to distinguish the fits.",
    "touch":"Counts an undirected edge once if either endpoint is in the top 121.",
    "prefix":"Largest fully connected prefix of the degree ranking. Maximum clique within the top 50 differs.",
    "clique":"Largest fully connected prefix, not the maximum clique.",
    "pathstar":"Computed within the giant component, confirmed against b02_zhou_table.csv.",
    "clustering":"Local coefficient is set to zero for nodes of degree 0 or 1.",
    "diam":"Source means 5.06 / 3.22 / 4.00 / 4.88. The current formatter rounds to integers before adding a decimal.",
    "Ravg":"Uses zero for empty neighbour sets. Null keeps every degree and the whole degree-threshold φ(k) curve.",
}
def note_diagnose(label,key):
    if key not in refs:
        return ("within noise" if key=="range" else "yes"), notes.get(key,"")
    row=real.loc[label]
    pieces=[]
    beyond=False
    for ours, kind, stored, storedkind in refs[key]:
        if kind=="phi":
            source=struct.loc[(row.date,row["mode"])]
            original,sd,n=source.rcpres_R_avg_mean,source.rcpres_R_avg_sd,source.rcpres_n
            x=constraints.loc[constraints.label.eq(label)&constraints.model.eq("phi"),ours]
        else:
            source=history.loc[(row.date,row["mode"],stored,storedkind)]
            original,sd,n=source.null_mean,source.null_sd,source.n_null
            x=rand.loc[rand.label.eq(label)&rand.model.eq(kind),ours]
        if label=="Aug" and ours=="path" and kind=="er":
            x=pd.concat([x,pd.read_csv(scratch/"aug_er_check.csv").l])
        diff=x.mean()-original
        se=np.sqrt(x.var(ddof=1)/len(x)+sd**2/n)
        beyond|=abs(diff)>2.58*se+1e-12
        pieces.append(f"{kind} Δ {diff:+.4g}, MC SE {se:.3g}")
    if key=="diam" and label!="Nov":
        return "NO",notes[key]+" "+"; ".join(pieces)+"."
    return ("NO" if beyond else "within noise"),(notes.get(key,"")+" "+"; ".join(pieces)+".").strip()

def main():
    out=["## Tables 1 and A1, cell by cell","",
         "Each row below is one cell in the manuscript. Dates are file dates. AugEq keeps 35,812 directed links. "
         "For null means, MC SE is the combined Monte Carlo standard error of the independent and reported means. "
         "I use 2.58 combined SE as a screening threshold across the many cells. "
         "Null estimates use 16 draws per graph and null family, with 20 for the φ-preserving null and 116 for the August ER path check.","",
         "| Claim | Where it appears | Reported value | Independent value | Match | Note |",
         "|---|---|---:|---:|---|---|"]
    for file, mapping in [("table1.md",rows1),("tableA1.md",rowsA)]:
        source=(Path("reports/tables")/file).read_text().splitlines()
        for line,key in mapping.items():
            cells=[v.strip() for v in source[line-1].strip("|").split("|")]
            for label,txt in zip(labels,cells[1:]):
                status,note=note_diagnose(label,key)
                out.append(f"| {cells[0]} ({label}) | [{file}:{line}](../reports/tables/{file}#L{line}) | {txt} | {calculate(label,key)} | {status} | {note} |")
    (scratch/"table_cells.md").write_text("\n".join(out)+"\n")

if __name__=="__main__":
    main()

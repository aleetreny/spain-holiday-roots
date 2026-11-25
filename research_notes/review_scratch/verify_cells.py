"""Compare every TeX data cell with its corresponding CSV display value."""
from pathlib import Path
import re
import pandas as pd

snap = Path("research_notes/review_scratch/snapshot")
t = pd.read_csv("results/t01_table1.csv")
rows = t.row.unique().tolist()

def canon(s):
    s = s.replace(r"$^\dagger$", "†").replace(r"$^{\circ}$", "°")
    s = s.replace("−", "-").replace(r"\%", "%").replace("$", "")
    return re.sub(r"\s+", "", s).lower().replace("\\", "")

def latex_cells(file, line):
    row = (snap/file).read_text().splitlines()[line-1]
    return [canon(s) for s in row.split("&")[1:]]

def source_cells(i):
    return list(t[t.row.eq(rows[i])].cell)

a1_rows = [12,13,14,15,16,19,20,21,22,23,26,27,28,29,32,33,34,35,38,39,40,41]
for i, line in enumerate(a1_rows):
    assert latex_cells("tableA1.tex", line) == list(map(canon, source_cells(i))), (i,line)
t1_rows = [2,3,4,5,7,8,16,14,18,19,10,11,13]
for line, i in enumerate(t1_rows,11):
    out = source_cells(i)[:3]
    if i == 2:
        out = [s.split("(")[1].strip(")") for s in out]
    if i == 7:
        out = [s.replace("(", "/ ").strip(")") for s in out]
    if i in {10,19}:
        out = [s.split("(")[0].strip() for s in out]
    assert latex_cells("table1.tex", line) == list(map(canon, out)), (i,line)

v = pd.read_csv("results/v1_c1_zero_sensitivity.csv")
for line, (_, group) in enumerate(v.groupby("spec", sort=False),10):
    if line == 17:
        break
    vals = group.set_index("term").loc[["pre1981", "recent_2016_23", "village_le2k"]]
    cells = [f"{r.delta_aug_feb:+.3f}" + ("°" if r.p_wcr_delta >= .05 else "") for _,r in vals.iterrows()]
    assert latex_cells("tableA2.tex", line) == list(map(canon, cells)), line
print("Matched 39 Table 1, 88 Table A1, and 21 Table A2 display cells.")

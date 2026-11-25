"""Table A2: sensitivity of the generation contrasts (Figure 3) to the treatment of empty city-province cells.
Reads results/v1_c1_zero_sensitivity.csv (pipeline/v1_verify_n1.py); writes reports/tables/tableA2.tex.
Run: .venv/bin/python -I pipeline/t02_zero_cells.py
"""
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results" / "v1_c1_zero_sensitivity.csv")
SPECS = [("N1_haldane_0.5_sample", "Add 0.5 to each count (main)"),
         ("plus_1.0_sample_eq_log1p", "Add 1 to each count"),
         ("plus_0.1_sample_eq_pop_plus1", "Add 0.1 to each count"),
         ("log_share_plus_1e4", r"Add $10^{-4}$ to each share"),
         ("log_share_plus_1e5", r"Add $10^{-5}$ to each share"),
         ("asinh_1000_share", r"$\operatorname{asinh}(1000\times\text{share})$"),
         ("R2_112_large_cities_G22_haldane0.5", "112 large cities, 22 provinces"),]
TERMS = [("pre1981", "Arrived before 1981"), ("recent_2016_23", "Arrived 2016--2023"), ("village_le2k", "Born in a village")]
def cell(r):
    star = "" if r.p_wcr_delta < 0.05 else r"$^{\circ}$"
    return f"${r.delta_aug_feb:+.3f}$".replace("+", "+").replace("-", "-") + star
rows = []
for spec, lab in SPECS:
    s = d[d.spec == spec].set_index("term")
    rows.append(lab + " & " + " & ".join(cell(s.loc[t]) for t, _ in TERMS) + r" \\")
tex = r"""\begin{table}[t]
\centering
\small
\caption{Change in the pull of the roots, 14--15 August minus 13--14 February, under different treatments of empty city--province cells}
\label{tab:zero}
\begin{tabular}{lrrr}
\toprule
Treatment & """ + " & ".join(l for _, l in TERMS) + r""" \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\par\smallskip
\begin{minipage}{\linewidth}\footnotesize
All destinations; groups entered together as in Figure~\ref{fig:roots}. Between 46\,\% and 60\,\% of city--province cells have no sampled resident of a group, so the regressor needs a rule for empty cells. $^{\circ}$\,Wild cluster bootstrap $p\geq0.05$. Source: results/v1\_c1\_zero\_sensitivity.csv.
\end{minipage}
\end{table}
"""
(ROOT / "reports" / "tables" / "tableA2.tex").write_text(tex)
print(tex)

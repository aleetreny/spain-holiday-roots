"""Table 1 of the COMP0123 report: basic properties of the overnight-stay networks.

Reads only finished pipeline outputs (no recomputation):
  results/b01_networks/manifest.csv     size, people away, share kept by the backbone
  results/b02_basic_properties.csv      real values and null means (ER, degree-preserving rewiring)
  results/b02_zhou_table.csv            l* = ln N / ln <k> (giant component), top-5 % link share
  results/b02_structural_constraints.csv  [r_min, r_max] allowed by P(k); R_avg null preserving P(k) and phi(k)
  results/b02_powerlaw.csv              Clauset-Shalizi-Newman power law vs log-normal
  results/b02_representativeness.csv    sister-night agreement (dagger marks)
Writes results/t01_table1.csv (every cell with its source), reports/tables/table1.{md,tex} (main text,
13 rows, Zhou & Mondragon 2007 style) and reports/tables/tableA1.{md,tex} (appendix, all rows + equal-links column).
Run: .venv/bin/python -I pipeline/t01_table1.py
"""
from pathlib import Path
import math
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = ROOT / "reports" / "tables"; OUT.mkdir(parents=True, exist_ok=True)

# Columns: (key, file date, backbone, header). Date convention: file d = night (d-1 -> d).
COLS = [
    ("ordinary", "2024-02-14", "w10", "Ordinary weekday", "13–14 Feb"),
    ("summer", "2024-08-15", "w10", "Summer holiday", "14–15 Aug"),
    ("allsaints", "2024-11-01", "w10", "All Saints", "31 Oct–1 Nov"),
    ("summer_eq", "2024-08-15", "equal_edges", "Summer, same directed link count", "14–15 Aug"),
]
GROUP = {"ordinary": "ordinary", "summer": "summer", "allsaints": "allsaints", "summer_eq": "summer"}

man = pd.read_csv(RES / "b01_networks" / "manifest.csv", dtype={"date": str})
bp = pd.read_csv(RES / "b02_basic_properties.csv", dtype={"date": str})
zt = pd.read_csv(RES / "b02_zhou_table.csv", dtype={"date": str})
sc = pd.read_csv(RES / "b02_structural_constraints.csv", dtype={"date": str})
pl = pd.read_csv(RES / "b02_powerlaw.csv", dtype={"date": str})
rep = pd.read_csv(RES / "b02_representativeness.csv")


def one(df, **kw):
    m = pd.Series(True, index=df.index)
    for k, v in kw.items():
        m &= df[k] == v
    sub = df[m]
    assert len(sub) == 1, (kw, len(sub))
    return sub.iloc[0]


def B(c, metric, null="none", field="real"):
    _, d, bb, *_ = c
    if field == "real":  # the real value is repeated on every null row; some metrics have no 'none' row
        sub = bp[(bp.date == d) & (bp.backbone == bb) & (bp.metric == metric)]
        assert len(sub) >= 1 and (sub.real.max() - sub.real.min()) <= 1e-9 * max(1.0, abs(sub.real.iloc[0])), (d, bb, metric)
        return float(sub.real.iloc[0]), f"b02_basic_properties.csv[{metric},real]"
    v = one(bp, date=d, backbone=bb, metric=metric, null=null)[field]
    return float(v), f"b02_basic_properties.csv[{metric},{null},{field}]"


def M(c, field):
    _, d, bb, *_ = c
    return float(one(man, date=d, backbone=bb)[field]), f"b01_networks/manifest.csv[{field}]"


def Z(c, field):
    _, d, bb, *_ = c
    return float(one(zt, date=d, backbone=bb)[field]), f"b02_zhou_table.csv[{field}]"


def S(c, field):
    _, d, bb, *_ = c
    return float(one(sc, date=d, backbone=bb)[field]), f"b02_structural_constraints.csv[{field}]"


def dagger(c, metric):
    """Holiday columns: mark the cell if the change from the ordinary night does not keep its
    direction in every sister-night pair (b02_representativeness: sister_pairs_agree < sister_pairs)."""
    key, d, bb, *_ = c
    if key == "ordinary":
        return ""
    r = rep[(rep.backbone == bb) & (rep.group == GROUP[key]) & (rep.metric == metric)]
    if r.empty:
        return ""
    r = r.iloc[0]
    return "†" if r.sister_pairs_agree < r.sister_pairs else ""


def f_int(x):
    return f"{int(round(x)):,}"


def f(x, nd):
    s = f"{x:.{nd}f}"
    return s.replace("-", "−")


def csn(c):
    """Tail of the in- and out-degree distributions: power law vs log-normal (CSN likelihood ratio)."""
    _, d, bb, *_ = c
    out = []
    for kind in ("in", "out"):
        r = one(pl, date=d, backbone=bb, kind=kind)
        if r.p_pl_vs_lognormal >= 0.05:
            out.append("tie")
        else:
            out.append("log-normal" if r.R_pl_vs_lognormal < 0 else "power law")
    return " / ".join(out), "b02_powerlaw.csv[R_pl_vs_lognormal,p_pl_vs_lognormal]"


# Each row: (section, label, function c -> (text, source, raw value or None))
def row_defs():
    R = []

    def add(section, label, fn):
        R.append((section, label, fn))

    def realnull(metric, nulls, nd, scale=1.0, null_nd=None):
        null_nd = nd if null_nd is None else null_nd

        def fn(c):
            v, s = B(c, metric)
            parts = [f(v * scale, nd) + dagger(c, metric)]
            srcs = [s]
            for nl in nulls:
                m, s2 = B(c, metric, nl, "null_mean")
                parts.append(f(m * scale, null_nd)); srcs.append(s2)
            return parts[0] + (" (" + " / ".join(parts[1:]) + ")" if nulls else ""), "; ".join(srcs), v
        return fn

    # --- size
    add("Size", "People sleeping away from home (millions)",
        lambda c: (f(M(c, "total_weight_full")[0] / 1e6, 2), M(c, "total_weight_full")[1], M(c, "total_weight_full")[0]))
    add("Size", "Share of them on kept links (%)",
        lambda c: (f(100 * M(c, "share_weight_retained")[0], 0), M(c, "share_weight_retained")[1], M(c, "share_weight_retained")[0]))
    add("Size", "Nodes N (in the giant component)",
        lambda c: (f"{f_int(B(c, 'N')[0])} ({f_int(B(c, 'giant_weak')[0])})", B(c, "giant_weak")[1], B(c, "N")[0]))
    add("Size", "Links L (undirected)", lambda c: (f_int(B(c, "E_und")[0]) + dagger(c, "E_und"), B(c, "E_und")[1], B(c, "E_und")[0]))
    add("Size", "Density (%)", realnull("density_und", [], 2, 100))
    # --- degree
    add("Degree", "Mean degree ⟨k⟩", realnull("k_mean", [], 1))
    add("Degree", "Median degree", realnull("k_median", [], 0))
    add("Degree", "Maximum degree k_max (% of N)",
        lambda c: (f"{f_int(B(c, 'k_max')[0])} ({100 * B(c, 'k_max')[0] / B(c, 'N')[0]:.0f} %)" + dagger(c, "k_max"),
                   B(c, "k_max")[1], B(c, "k_max")[0]))
    add("Degree", "Median out-degree / in-degree",
        lambda c: (f"{f_int(B(c, 'kout_median')[0])} / {f_int(B(c, 'kin_median')[0])}",
                   B(c, "kout_median")[1] + "; " + B(c, "kin_median")[1], None))
    add("Degree", "Tail fit, in / out: power law vs log-normal", lambda c: (*csn(c), None))
    # --- small world
    add("Small world", "Average path length ℓ (rewired / ER)", realnull("L", ["rew_und", "er_und"], 2))
    add("Small world", "ℓ* = ln N / ln⟨k⟩", lambda c: (f(Z(c, "l_star")[0], 2), Z(c, "l_star")[1], Z(c, "l_star")[0]))
    add("Small world", "Diameter (rewired)", realnull("diameter", ["rew_und"], 0, null_nd=1))
    add("Small world", "Clustering ⟨C⟩ (rewired / ER)", realnull("C_avg", ["rew_und", "er_und"], 3))
    # --- mixing
    add("Mixing", "Degree assortativity r (rewired)", realnull("r_und", ["rew_und"], 3))
    add("Mixing", "Range of r reached by rewiring [r_min, r_max]",
        lambda c: (f"[{f(S(c, 'r_min_maxdisassortative')[0], 3)}, {f(S(c, 'r_max_maxassortative')[0], 3)}]",
                   S(c, "r_min_maxdisassortative")[1] + "; " + S(c, "r_max_maxassortative")[1], None))
    add("Mixing", "Reciprocity (rewired, directed)", realnull("reciprocity", ["rew_dir"], 2))
    add("Mixing", "2nd-order mixing R_avg (null keeping P(k) and φ)",
        lambda c: (f"{f(S(c, 'rcpres_R_avg_real')[0], 3)}{dagger(c, 'R_avg')} ({f(S(c, 'rcpres_R_avg_mean')[0], 3)})",
                   S(c, "rcpres_R_avg_real")[1] + "; " + S(c, "rcpres_R_avg_mean")[1], S(c, "rcpres_R_avg_real")[0]))
    # --- rich club and core
    add("Rich club and core", "Rich club φ, top 50 zones (rewired)", realnull("phi_top50", ["rew_und"], 2))
    add("Rich club and core", "Rich-club clique size (rewired)", realnull("n_clique", ["rew_und"], 0, null_nd=1))
    add("Rich club and core", "Links touching the top 5 % of nodes (%)",
        lambda c: (f(100 * Z(c, "share_l_top5_any")[0], 0) + dagger(c, "l_top5_any"), Z(c, "share_l_top5_any")[1], Z(c, "share_l_top5_any")[0]))
    add("Rich club and core", "Largest k of a non-empty k-core (rewired / ER)", realnull("kcore_max", ["rew_und", "er_und"], 0))
    return R


rows, tidy = [], []
for section, label, fn in row_defs():
    cells = []
    for c in COLS:
        txt, src, raw = fn(c)
        cells.append(txt)
        tidy.append({"section": section, "row": label, "column": c[0], "date": c[1], "backbone": c[2],
                     "cell": txt, "value": raw, "source": src})
    rows.append((section, label, cells))

pd.DataFrame(tidy).to_csv(RES / "t01_table1.csv", index=False)

NOTE_FULL = (
    "Each column is the network of one night of 2024: node = one of 2,435 MITMA zones of mainland Spain; "
    "a directed link i→j joins zones when at least 10 residents of i sleep in j (last column: the summer night cut "
    "to the same number of directed links as the ordinary night, 35,812, by keeping the heaviest). Properties are "
    "computed on the undirected projection unless marked directed. In brackets: mean of 50 networks with the same "
    "degree sequence (degree-preserving rewiring, 10 swaps per link) and, where shown, of 50 Erdős–Rényi networks "
    "with the same N and L. [r_min, r_max]: the most disassortative and most assortative values reached by "
    "a heuristic rewiring search that keeps P(k) and a single connected component (Zhou & Mondragón 2007). R_avg: second-order mixing "
    "(Zhou, Cox & Hansen 2017), null keeping P(k) and φ(k). Rich-club clique: the largest n such that the n "
    "best-connected zones are all linked to one another (Zhou & Mondragón 2007). Tail fit: Clauset–Shalizi–Newman likelihood-ratio test; "
    "'tie' means p ≥ 0.05. † The change from the ordinary night does not keep its direction in every pair of "
    "sister nights (117 summer pairs, 39 All Saints pairs). 31 Oct–1 Nov is the quietest of its three sister nights. "
    ""
)
NOTE_MAIN = (
    "Undirected networks of 2,435 zones, one per night of 2024; median degrees and reciprocity use the "
    "directed links. In brackets: mean of 50 degree-preserving rewirings (for clustering, also 50 Erdős–Rényi "
    "graphs). † The direction of the change is not the same on every pair of sister nights. "
    "Definitions and the full version: Table A1."
)

# Main text: Zhou & Mondragón (2007) Table 1 style, three networks, 12 rows, no sections.
# (label in the full table, label in the main table, keep brackets?)
MAIN_ROWS = [
    ("Nodes N (in the giant component)", "Nodes N in the giant component", True),
    ("Links L (undirected)", "Links L", True),
    ("Density (%)", "Density (%)", True),
    ("Mean degree ⟨k⟩", "Average degree ⟨k⟩", True),
    ("Maximum degree k_max (% of N)", "Maximum degree k_max / share of all zones", True),
    ("Median out-degree / in-degree", "Median out- / in-degree", True),
    ("Reciprocity (rewired, directed)", "Reciprocity (rewired)", True),
    ("Degree assortativity r (rewired)", "Assortativity r (rewired)", True),
    ("Rich club φ, top 50 zones (rewired)", "Rich club φ, top 50 (rewired)", True),
    ("Rich-club clique size (rewired)", "Rich-club clique n_clique", False),
    ("Average path length ℓ (rewired / ER)", "Average path length ℓ", False),
    ("ℓ* = ln N / ln⟨k⟩", "Random-graph estimate ℓ* = ln N / ln⟨k⟩", True),
    ("Clustering ⟨C⟩ (rewired / ER)", "Clustering ⟨C⟩ (rewired / random)", True),
]


def tex(s):
    rep_ = {"%": r"\%", "⟨": r"$\langle$", "⟩": r"$\rangle$", "ℓ*": r"$\ell^*$", "ℓ": r"$\ell$", "φ": r"$\phi$",
            "†": r"$^\dagger$", "−": r"$-$", "→": r"$\to$", "≥": r"$\geq$", "&": r"\&"}
    for k, v in {"n_clique": r"$n_{\mathrm{clique}}$", "k_max": r"$k_{\max}$", "r_min": r"$r_{\min}$",
                 "r_max": r"$r_{\max}$", "R_avg": r"$R_{\mathrm{avg}}$"}.items():
        s = s.replace(k, v)
    for k, v in rep_.items():
        s = s.replace(k, v)
    return s


def write(name, title, label, cols, rows_, note, sections):
    hdr = "| | " + " | ".join(f"{c[3]}<br>{c[4]}" for c in cols) + " |"
    md = [f"# {title}", "", hdr, "|---|" + "---:|" * len(cols)]
    tx = [r"\begin{table}[" + ("!h" if name == "tableA1" else "t") + "]", r"\centering", r"\small", r"\setlength{\tabcolsep}{3pt}", rf"\caption{{{title.split('. ', 1)[1]}}}",
          rf"\label{{{label}}}", r"\resizebox{\linewidth}{!}{\begin{tabular}{l" + "r" * len(cols) + "}", r"\toprule",
          " & " + " & ".join(rf"\shortstack[r]{{{c[3]}\\{c[4]}}}" for c in cols) + r" \\", r"\midrule"]
    last = None
    for section, lab, cells in rows_:
        if sections and section != last:
            md.append(f"| **{section}** |" + " |" * len(cols))
            if last is not None:
                tx.append(r"\addlinespace")
            tx.append(rf"\multicolumn{{{len(cols) + 1}}}{{l}}{{\textit{{{section}}}}} \\")
            last = section
        md.append(f"| {lab} | " + " | ".join(cells) + " |")
        tx.append((r"\quad " if sections else "") + tex(lab) + " & " + " & ".join(tex(x) for x in cells) + r" \\")
    md += ["", note, ""]
    tx += [r"\bottomrule", r"\end{tabular}}", r"\par\smallskip", r"\begin{minipage}{\linewidth}\footnotesize",
           tex(note.replace(" Source of every cell: results/t01_table1.csv.", "")), r"\end{minipage}", r"\end{table}"]
    (OUT / f"{name}.md").write_text("\n".join(md))
    (OUT / f"{name}.tex").write_text("\n".join(tx) + "\n")
    return "\n".join(md)


write("tableA1", "Table A1. Basic properties of the overnight-stay networks, full version", "tab:basic_full",
      COLS, rows, NOTE_FULL, sections=True)

by_label = {lab: cells for _, lab, cells in rows}
main = []
for full, short, brackets in MAIN_ROWS:
    cells = by_label[full][:3]
    if full.startswith("Nodes"):  # giant component only, as in Zhou & Mondragon (2007); brackets mean nulls here
        cells = [x.partition(" (")[2].rstrip(")") for x in cells]
    elif full.startswith("Maximum degree"):
        cells = [x.replace(" (", " / ").replace(" %)", " %") for x in cells]
    elif not brackets:
        cells = [x.partition(" (")[0] for x in cells]
    main.append(("", short, cells))
print(write("table1", "Table 1. Basic properties of the overnight-stay networks", "tab:basic",
            COLS[:3], main, NOTE_MAIN, sections=False))

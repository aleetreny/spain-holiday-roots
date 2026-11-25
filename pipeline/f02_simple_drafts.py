"""Second round of figure drafts: one simple chart per message, no link maps.

Visual system: serif type; warm grey = ordinary night; sea blue = coast / beach;
terracotta = inland / roots; dark ink for text. Run with .venv/bin/python -I.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "figures_drafts"; OUT.mkdir(parents=True, exist_ok=True)
RES, NET = ROOT / "results", ROOT / "results" / "b01_networks"
BLUE, TERRA, GREY, INK, MUTED, GRID = "#1f6e8c", "#c0623b", "#a8a39b", "#1d1d1b", "#6b6862", "#e6e2db"
plt.rcParams.update({"font.family": "STIX Two Text", "mathtext.fontset": "stix", "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "savefig.dpi": 300})
nodes = pd.read_parquet(NET / "nodes.parquet").set_index("idx").sort_index()
X, Y, COAST = nodes.x.values, nodes.y.values, nodes.coastal_true.values.astype(bool)

# ---- Fig 1: who sleeps far from home, and where (coast vs inland) ----
def fig_where():
    nights = [("2024-02-14", "Ordinary Wednesday, 14 Feb"), ("2024-03-29", "Good Friday, 29 Mar"),
              ("2024-08-15", "Summer holiday, 15 Aug"), ("2024-11-01", "All Saints, 1 Nov"),
              ("2024-12-25", "Christmas Eve, 24–25 Dec")]
    rows = []
    for d, lab in nights:
        e = pd.read_parquet(NET / f"edges_{d}__full.parquet"); e = e[e.i != e.j]
        km = np.hypot(X[e.i] - X[e.j], Y[e.i] - Y[e.j]) / 1000
        f = e[km >= 100]
        rows.append((lab, f.w[COAST[f.j]].sum() / 1e6, f.w[~COAST[f.j]].sum() / 1e6))
    df = pd.DataFrame(rows, columns=["night", "coast", "inland"])
    fig, ax = plt.subplots(figsize=(6.3, 2.6))
    y = np.arange(len(df))[::-1]
    ax.barh(y, df.inland, color=TERRA, height=0.62, label="inland")
    ax.barh(y, df.coast, left=df.inland + 0.03, color=BLUE, height=0.62, label="at the coast")
    for yi, r in zip(y, df.itertuples()):
        tot = r.coast + r.inland
        ax.text(tot + 0.08, yi, f"{tot:.1f} M  ({r.inland/tot:.0%} inland)", va="center", fontsize=8.5, color=INK)
    ax.set_yticks(y, df.night); ax.tick_params(axis="y", length=0, labelcolor=INK)
    ax.set_xlabel("People sleeping at least 100 km from home (millions)")
    ax.xaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True); ax.set_xlim(0, 5.6)
    ax.spines["left"].set_visible(False)
    ax.legend(frameon=False, loc="lower right", ncol=2, fontsize=8.5, bbox_to_anchor=(1, 1.0))
    fig.savefig(OUT / "G1_where_people_sleep.png", bbox_inches="tight"); plt.close(fig)
    print(df)

# ---- Fig 2: senders and receivers (k_in vs k_out) ----
def fig_split():
    m = pd.read_parquet(RES / "p_node_metrics.parquet")
    big = nodes.population.values >= 100_000
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 3.3), sharex=True, sharey=True)
    for ax, (d, lab) in zip(axs, [("2024-02-14", "Ordinary night (14 Feb)"), ("2024-08-15", "Holiday night (15 Aug)")]):
        s = m[m.date == d].set_index("idx").reindex(nodes.index)
        kin, kout = s.kin_w10.values + 1, s.kout_w10.values + 1
        for mask, c, lbl, z in ((~COAST & ~big, TERRA, "inland zone", 2), (COAST & ~big, BLUE, "coastal zone", 3),
                                (big, INK, "city > 100,000", 4)):
            ax.scatter(kout[mask], kin[mask], s=5 if lbl != "city > 100,000" else 12, color=c, alpha=0.45 if z < 4 else 0.9,
                       lw=0, label=lbl, zorder=z)
        ax.plot([1, 3000], [1, 3000], color=MUTED, lw=0.7, ls="--", zorder=1)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1, 3000); ax.set_ylim(1, 3000)
        rho = pd.Series(kin).corr(pd.Series(kout), method="spearman")
        ax.set_title(f"{lab}\nSpearman correlation {rho:.2f}", loc="left", fontsize=9.5)
        ax.set_xlabel("zones it sends people to (out-degree + 1)")
    axs[0].set_ylabel("zones it receives people from (in-degree + 1)")
    axs[1].text(1500, 4, "senders\n(cities)", ha="right", fontsize=8.5, color=INK)
    axs[1].text(1.6, 1300, "receivers\n(coast, villages)", va="top", fontsize=8.5, color=INK)
    axs[0].legend(frameon=False, fontsize=8, loc="upper left", markerscale=2.5)
    fig.savefig(OUT / "G2_senders_receivers.png", bbox_inches="tight"); plt.close(fig)

# ---- Fig 3: the generation effect (dot + 95% CI) ----
def fig_roots():
    c = pd.read_csv(RES / "v1_c2_main_contrasts.csv")
    c = c[c.contrast == "Aug15_vs_Feb14"]
    terms = [("pre1981_exodus", "Arrived before 1981\n(the exodus generation)"),
             ("born_village_le2k", "Born in a village\n(2,000 people or fewer)"),
             ("born_urban_gt10k", "Born in a town\n(> 10,000 people)"),
             ("recent_2016_2024", "Arrived 2016–2023\n(recent migrants)")]
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)
    for ax, (smp, lab, col) in zip(axs, [("inland", "Trips to inland destinations", TERRA),
                                         ("coastal", "Trips to coastal destinations", BLUE)]):
        for k, (t, _) in enumerate(terms):
            r = c[(c["sample"] == smp) & (c.term == t)].iloc[0]
            yk = len(terms) - 1 - k
            ax.plot([r.ci_wcr_lo, r.ci_wcr_hi], [yk, yk], color=col, lw=2, solid_capstyle="round")
            ax.scatter([r.delta], [yk], s=36, color=col, edgecolor="white", lw=1.2, zorder=3)
        ax.axvline(0, color=MUTED, lw=0.8)
        ax.set_title(lab, loc="left", fontsize=9.5); ax.set_xlim(-0.25, 0.3)
        ax.xaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.set_xlabel("Extra pull on 15 Aug vs 14 Feb (elasticity)")
        ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    axs[0].set_yticks(range(len(terms))[::-1], [l for _, l in terms], color=INK)
    fig.savefig(OUT / "G3_roots_by_generation.png", bbox_inches="tight"); plt.close(fig)

# ---- Fig 4: villages switch from senders to receivers ----
def fig_villages():
    r = pd.read_csv(RES / "r1_net_balance_by_quintile.csv")
    r = r[(r["sample"] == "all") & (r.role == "main") & r.cls.str.match(r"Q\d")]
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    x = np.arange(5)
    for d, lab, col in [("2024-02-14", "Ordinary night, 14 Feb", GREY), ("2024-08-15", "Summer, 15 Aug", BLUE),
                        ("2024-12-25", "Christmas Eve", TERRA)]:
        s = r[r.date == d].sort_values("cls")
        ax.plot(x, 100 * s.share_net_importer, color=col, lw=2, marker="o", ms=6, mec="white", mew=1.2)
        yl = {"2024-08-15": 93, "2024-12-25": 104}.get(d, 100 * s.share_net_importer.iloc[-1])
        ax.text(4.15, yl, lab, va="center", fontsize=8.5, color=INK)
    ax.axhline(50, color=MUTED, lw=0.7, ls="--")
    ax.text(2.0, 51.5, "above: more arrive than leave", fontsize=7.5, color=MUTED)
    ax.text(2.0, 33, "below: more leave than arrive", fontsize=7.5, color=MUTED)
    ax.set_xticks(x, ["least\nemptied", "", "", "", "most\nemptied"])
    ax.set_xlabel("Zones grouped by population lost since its peak (fifths)")
    ax.set_ylabel("Zones receiving more people\nthan they send (%)")
    ax.set_ylim(0, 105); ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    fig.savefig(OUT / "G4_villages_switch.png", bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    for f in (fig_where, fig_split, fig_roots, fig_villages): f()

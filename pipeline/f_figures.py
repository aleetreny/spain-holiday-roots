"""Final figures for the COMP0123 report (one function per figure).

Visual system (all figures): STIX Two Text (serif, like the text); sea blue = coast;
terracotta = inland towns and villages ("roots"); light stone = inland cities;
warm grey = ordinary night; dark ink for text. Text never wears series colour.
Date convention: the MITMA file dated d holds the night ending on d (d-1 -> d).
Run: .venv/bin/python -I pipeline/f_figures.py [fig1 ...]
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "figures"; OUT.mkdir(parents=True, exist_ok=True)
RES, NET = ROOT / "results", ROOT / "results" / "b01_networks"
# Style copied from Nature / Nature Communications figures (Wesolowski 2017; Schlapfer 2021):
# Helvetica, small type, closed thin frame, outward ticks, bold lowercase panel letters,
# frameless legends inside the panel, dashed vertical lines for holidays.
BLUE, TERRA, STONE, GREY = "#2b6f9e", "#c0623b", "#c9c9c9", "#a8a39b"
INK, MUTED, GRID = "#222222", "#666666", "#ececec"
plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 8.5, "axes.titlesize": 8.5, "axes.labelsize": 8.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.linewidth": 0.6, "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": True, "axes.spines.right": True, "legend.frameon": False,
    "savefig.dpi": 600, "pdf.fonttype": 42, "lines.solid_capstyle": "butt"})

def panel_letter(ax, L, x=-0.08, y=1.02):
    ax.text(x, y, L, transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="bottom", ha="right")
TEXTWIDTH = 6.3  # inches, A4 with 2.5 cm margins

nodes = pd.read_parquet(NET / "nodes.parquet").set_index("idx").sort_index()
X, Y = nodes.x.values, nodes.y.values
COAST = nodes.coastal_true.values.astype(bool)
CITY = nodes.population.values >= 100_000

def save(fig, name):
    # PDF for the report: vector lines and text; raster layers (dense scatters, matrices) at 1200 dpi.
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", dpi=1200)
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=600)
    plt.close(fig)

# ---------------------------------------------------------------- Figure 1
FIG1_NIGHTS = [  # (file date, label: the night itself)
    ("2024-02-14", "An ordinary weekday\nnight of 13–14 February"),
    ("2024-03-29", "Easter\nnight of 28–29 March"),
    ("2024-08-15", "Summer holiday\nnight of 14–15 August"),
    ("2024-11-01", "All Saints\nnight of 31 October–1 November"),
    ("2024-12-25", "Christmas Eve\nnight of 24–25 December"),
]
FIG1_MIN_KM = 100

def fig1_data():
    rows = []
    for d, lab in FIG1_NIGHTS:
        e = pd.read_parquet(NET / f"edges_{d}__full.parquet"); e = e[e.i != e.j]
        e = e[np.hypot(X[e.i] - X[e.j], Y[e.i] - Y[e.j]) / 1000 >= FIG1_MIN_KM]
        j = e.j.values
        rows.append(dict(date=d, label=lab,
                         towns=e.w[~COAST[j] & ~CITY[j]].sum() / 1e6,
                         cities=e.w[~COAST[j] & CITY[j]].sum() / 1e6,
                         coast=e.w[COAST[j]].sum() / 1e6))
    df = pd.DataFrame(rows); df["total"] = df[["towns", "cities", "coast"]].sum(axis=1)
    df.to_csv(OUT / "fig1_data.csv", index=False)
    return df

FIG1_SHORT = ["Ordinary weekday\n13–14 Feb", "Easter\n28–29 Mar", "Summer holiday\n14–15 Aug",
              "All Saints\n31 Oct–1 Nov", "Christmas Eve\n24–25 Dec"]
SEGS = [("towns", TERRA, "Inland towns and villages", "white"),
        ("cities", STONE, "Inland cities (>100,000 people)", INK),
        ("coast", BLUE, "Coast", "white")]

def bars_abs(ax, df, labels, label_min=0.18, total_fmt="{:.1f} M"):
    y = np.arange(len(df))[::-1]
    left = np.zeros(len(df))
    for col, colour, lab, txt in SEGS:
        ax.barh(y, df[col], left=left, height=0.66, color=colour, label=lab,
                edgecolor="white", linewidth=0.8)
        for yi, l, v in zip(y, left, df[col]):
            if v >= label_min:
                ax.text(l + v / 2, yi, f"{v:.1f}", ha="center", va="center", fontsize=8, color=txt)
        left = left + df[col].values
    for yi, tot in zip(y, df.total):
        ax.text(tot + 0.06, yi, total_fmt.format(tot), va="center", ha="left", fontsize=8)
    ax.set_yticks(y, labels); ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, len(df) - 0.4)
    ax.set_xlim(0, 5.3); ax.set_xticks(range(6))
    ax.set_xlabel(f"People sleeping at least {FIG1_MIN_KM} km from home (millions)")

def fig1():
    """Option 1: one panel, how many people and where they slept, five nights."""
    df = fig1_data()
    fig, ax = plt.subplots(figsize=(TEXTWIDTH, 2.6))
    bars_abs(ax, df, FIG1_SHORT)
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.02), handlelength=1.0, handleheight=1.0,
              borderaxespad=0.4, labelspacing=0.35)
    save(fig, "fig1_option1_bars")

def fig1_calendar_data():
    f = OUT / "fig1_calendar_data.csv"
    d = pd.read_csv(f, parse_dates=["date"])
    d["night"] = d.date - pd.Timedelta(days=1)   # file date d = night (d-1 -> d)
    for c in ("towns", "cities", "coast"): d[c] = d[c] / 1e6
    return d

def fig1b():
    """Option 2: a = the 2024 calendar night by night; b = the five nights studied."""
    df, cal = fig1_data(), fig1_calendar_data()
    blocks = [("2024-01-31", "2024-03-30"), ("2024-07-31", "2024-08-30"), ("2024-09-30", "2024-12-30")]
    days = [(pd.Timestamp(b_) - pd.Timestamp(a_)).days + 1 for a_, b_ in blocks]
    fig = plt.figure(figsize=(TEXTWIDTH, 4.9))
    gs = fig.add_gridspec(2, 3, width_ratios=days, height_ratios=[1.15, 1], wspace=0.05, hspace=0.48)
    axs = [fig.add_subplot(gs[0, k]) for k in range(3)]
    marks = {"2024-02-13": ("Ordinary\nweekday", "center"), "2024-03-28": ("Easter", "right"),
             "2024-08-14": ("Summer\nholiday", "center"), "2024-10-31": ("All Saints", "center"),
             "2024-12-24": ("Christmas\nEve", "right")}
    for k, (ax, (a_, b_)) in enumerate(zip(axs, blocks)):
        s = cal[(cal.night >= a_) & (cal.night <= b_)]
        ax.fill_between(s.night, s.towns, s.coast, where=s.towns >= s.coast, interpolate=True,
                        color=TERRA, alpha=0.16, lw=0)
        ax.fill_between(s.night, s.towns, s.coast, where=s.towns < s.coast, interpolate=True,
                        color=BLUE, alpha=0.20, lw=0)
        ax.plot(s.night, s.coast, color=BLUE, lw=1.1, label="Coast")
        ax.plot(s.night, s.towns, color=TERRA, lw=1.1, label="Inland towns and villages")
        for dstr, (lab, ha) in marks.items():
            t = pd.Timestamp(dstr)
            if pd.Timestamp(a_) <= t <= pd.Timestamp(b_):
                ax.axvline(t, color=INK, lw=0.6, ls=(0, (3, 2)), zorder=0)
                tx = t + pd.Timedelta(days=0.8 if ha == "right" else 0)
                ax.text(tx, 2.63, lab, ha=ha, va="bottom", fontsize=7.2, linespacing=1.0)
        ax.set_xlim(pd.Timestamp(a_) - pd.Timedelta(days=0.5), pd.Timestamp(b_) + pd.Timedelta(days=0.5))
        ax.set_ylim(0, 2.6)
        ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
        ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator(bymonthday=15))
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("1 %b"))
        if k > 0:
            ax.tick_params(axis="y", left=False, labelleft=False)
            ax.spines["left"].set_visible(False)
        if k < 2:
            ax.spines["right"].set_visible(False)
        # break marks between blocks
        kw = dict(transform=ax.transAxes, color=INK, clip_on=False, lw=0.6)
        if k < 2:
            for yy in (0, 1): ax.plot([1 - 0.012 * 92 / days[k], 1 + 0.012 * 92 / days[k]], [yy - 0.03, yy + 0.03], **kw)
        if k > 0:
            for yy in (0, 1): ax.plot([-0.012 * 92 / days[k], 0.012 * 92 / days[k]], [yy - 0.03, yy + 0.03], **kw)
    axs[0].set_ylabel(f"People sleeping at least\n{FIG1_MIN_KM} km from home (millions)")
    axs[1].set_xlabel("Night of 2024 (dated by the evening it starts; gaps = no data)", labelpad=4)
    h, l = axs[2].get_legend_handles_labels()
    axs[2].legend(h[::-1], l[::-1], loc="upper left", bbox_to_anchor=(0.36, 0.97), handlelength=1.0, handletextpad=0.4, borderaxespad=0.2, frameon=True, facecolor="white", edgecolor="none", framealpha=1,
                  labelspacing=0.35).set_zorder(5)
    for ax in axs:
        ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5])
    axb = fig.add_subplot(gs[1, :])
    bars_abs(axb, df, [l.replace("\n", ", ") for l in FIG1_SHORT])
    axb.set_xlim(0, 5.0)
    axb.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=3, handlelength=1.0, handleheight=1.0,
               borderaxespad=0.0, columnspacing=1.6)
    x0 = -0.075
    fig.text(x0, axs[0].get_position().y1 + 0.045, "a", fontsize=10.5, fontweight="bold", va="bottom")
    fig.text(x0, axb.get_position().y1 + 0.012, "b", fontsize=10.5, fontweight="bold", va="bottom")
    save(fig, "fig1_option2_calendar")

# ---------------------------------------------------------------- Figure 2
def fig2():
    """a, b: in- vs out-degree per zone (Feb vs Aug); c: Zhou's structural-constraint range of r."""
    m = pd.read_parquet(RES / "p_node_metrics.parquet")
    sc = pd.read_csv(RES / "b02_structural_constraints.csv").query("backbone == 'w10'")
    big = nodes.population.values >= 100_000
    groups = [(~big & ~COAST, TERRA, "Inland towns and villages", 4, 0.30, 2),
              (~big & COAST, BLUE, "Coastal towns and villages", 7, 0.75, 3),
              (big, INK, "Cities over 100,000 people", 16, 0.95, 4)]
    rng = np.random.default_rng(20261124)
    jit = lambda k: k * np.exp(rng.uniform(-0.12, 0.12, len(k)))  # display-only jitter (log scale)
    fig = plt.figure(figsize=(TEXTWIDTH, 5.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.2, 0.9], hspace=0.45, wspace=0.08)
    axa, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    lab_pts = {"Madrid": "Madrid", "Barcelona": "Barcelona", "Benidorm": "Benidorm (resort)"}
    for ax, d, ttl in ((axa, "2024-02-14", "Ordinary weekday (13–14 Feb)"),
                       (axb, "2024-08-15", "Summer holiday (14–15 Aug)")):
        x = m[m.date == d].set_index("idx").reindex(nodes.index)
        kout, kin = x.kout_w10.values + 1, x.kin_w10.values + 1
        rho = pd.Series(kin).corr(pd.Series(kout), method="spearman")
        jo, ji = jit(kout), jit(kin)
        for mask, c, lab, sz, al, z in groups:
            ax.scatter(np.where(big, kout, jo)[mask], np.where(big, kin, ji)[mask], s=sz, color=c, alpha=al, lw=0.3 if z == 4 else 0,
                       edgecolor="white", label=lab, zorder=z, rasterized=True)
        ax.plot([1, 4000], [1, 4000], color=INK, lw=0.6, ls=(0, (3, 2)), zorder=1)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.8, 4000); ax.set_ylim(0.8, 4000)
        ax.set_aspect("equal")
        ax.set_title(f"{ttl}\nin–out correlation (Spearman) = {rho:.2f}", fontsize=8.5, linespacing=1.4)
        ax.set_xlabel("Out-degree + 1")
        for i in nodes.index[nodes.name.isin(lab_pts)]:
            if d == "2024-02-14" and nodes.name[i] == "Benidorm":
                continue
            nm = nodes.name[i]
            if d != "2024-02-14" and nm == "Barcelona":
                xy_t, ha, va = (3600, 38), "right", "center"
            elif d != "2024-02-14" and nm == "Benidorm":
                xy_t, ha, va = (110, 676), "right", "center"
            else:
                xy_t, ha, va = (kout[i] / 1.6, kin[i] * 1.45), "right", "baseline"
            ax.annotate(lab_pts[nm], (kout[i], kin[i]), xytext=xy_t, textcoords="data",
                        ha=ha, va=va, fontsize=7.5, zorder=6,
                        arrowprops=dict(arrowstyle="-", color=INK, lw=0.5, shrinkA=1, shrinkB=2))
    axa.set_ylabel("In-degree + 1")
    axb.tick_params(axis="y", labelleft=False)
    axb.text(0.97, 0.05, "senders: residents sleep\nin far more zones than\nthe zones sending people here", transform=axb.transAxes, ha="right",
             va="bottom", fontsize=7.5, color=MUTED)
    axb.text(0.04, 0.96, "receivers: host people\nfrom many more zones\nthan they send to", transform=axb.transAxes, ha="left",
             va="top", fontsize=7.5, color=MUTED)
    from matplotlib.lines import Line2D
    gh = [Line2D([], [], marker="o", ls="", color=c, mec="white", mew=0.5, ms=m_, label=lab)
          for (_, c, lab, _, _, _), m_ in zip(groups, (5.5, 5.5, 7))]
    fig.legend(handles=gh, loc="lower center", ncol=3, bbox_to_anchor=(0.53, 0.985), handletextpad=0.3,
               columnspacing=1.8)
    # c: structural constraints (Zhou & Mondragon 2007)
    axc = fig.add_subplot(gs[1, :])
    rows = [("ordinary_2024-02-14__w10", "Ordinary weekday"), ("allsaints_2024-11-01__w10", "All Saints"),
            ("summer_2024-08-15__w10", "Summer holiday")]
    for k, (lab_, name) in enumerate(rows):
        r = sc[sc.label == lab_].iloc[0]; yk = len(rows) - 1 - k
        axc.plot([r.r_min_maxdisassortative, r.r_max_maxassortative], [yk, yk], color=STONE, lw=7,
                 solid_capstyle="butt", zorder=1)
        axc.plot([r.maxrand_conn_r_mean] * 2, [yk - 0.22, yk + 0.22], color=MUTED, lw=1.2, zorder=2)
        axc.scatter([r.r_real], [yk], s=30, color=INK, zorder=3, edgecolor="white", lw=0.8)
    axc.set_yticks(range(len(rows))[::-1], [n_ for _, n_ in rows]); axc.tick_params(axis="y", length=0)
    axc.set_ylim(-0.7, len(rows) - 0.3); axc.set_xlim(-0.45, -0.05)
    axc.set_xlabel("Degree assortativity $r$ (below 0: hubs link mainly to lower-degree zones)")
    from matplotlib.lines import Line2D
    hh = [Line2D([], [], color=STONE, lw=6, label="range reached by rewiring, same degrees"),
          Line2D([], [], color=MUTED, lw=1.2, label="random rewiring (mean of 50)"),
          Line2D([], [], marker="o", ls="", color=INK, mec="white", ms=5, label="real network")]
    axc.legend(handles=hh, loc="lower left", ncol=3, bbox_to_anchor=(0.0, 1.0), handlelength=1.6,
               columnspacing=1.2, borderaxespad=0.2)
    x0 = -0.075
    for ax, L in ((axa, "a"), (axb, "b"), (axc, "c")):
        pos = ax.get_position()
        fig.text(pos.x0 - (0.085 if L != "b" else 0.02), pos.y1 + (0.055 if L != "c" else 0.06), L,
                 fontsize=10.5, fontweight="bold", va="bottom")
    save(fig, "fig2_senders_receivers")
    print(sc[["label", "r_real", "r_min_maxdisassortative", "r_max_maxassortative", "maxrand_conn_r_mean"]])

# ---------------------------------------------------------------- Figure 3
def fig3():
    """The roots pull inland in August, and it belongs to the generation that left (v1 results)."""
    mc = pd.read_csv(RES / "v1_c2_main_contrasts.csv").query("contrast == 'Aug15_vs_Feb14'")
    ag = pd.read_csv(RES / "v1_c2_age_confound.csv").query("spec == 'E2a_split65_4way_exclusive' and sample == 'all'")
    sp = pd.read_csv(RES / "v1_c2_sister_pairs.csv").query("grid == 'summer_x_feb_all'")
    XLAB = "Change in the pull of the province of birth, 14–15 Aug vs 13–14 Feb\n(+0.1 ≈ 7% larger expected flow for twice the group's share)"
    fig = plt.figure(figsize=(TEXTWIDTH, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1], width_ratios=[1, 1], hspace=0.72, wspace=0.62)
    gtop = gs[0, :].subgridspec(1, 2, wspace=0.06)
    axa = fig.add_subplot(gtop[0, 0]); axa2 = fig.add_subplot(gtop[0, 1], sharey=axa)
    rows = [("pre1981_exodus", "Arrived before 1981\n(the exodus generation)"),
            ("born_village_le2k", "Born in a village\n(2,000 people or fewer)"),
            ("recent_2016_2024", "Arrived 2016–2023\n(recent migrants)")]
    def zero_guides(ax, y_top):
        ax.axvline(0, color=INK, lw=0.6)
        ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    for ax, smp, col, ttl in ((axa, "inland", TERRA, "To inland destinations"),
                              (axa2, "coastal", BLUE, "To coastal destinations")):
        for k, (t, lab) in enumerate(rows):
            yk = len(rows) - 1 - k
            r = mc[(mc["sample"] == smp) & (mc.term == t)].iloc[0]
            ax.plot([r.ci_wcr_lo, r.ci_wcr_hi], [yk] * 2, color=col, lw=1.6)
            ax.scatter([r.delta], [yk], s=30, color=col, edgecolor="white", lw=0.8, zorder=3)
        ax.set_xlim(-0.16, 0.24); ax.set_xticks([-0.1, 0, 0.1, 0.2])
        ax.set_ylim(-0.6, len(rows) - 0.4); zero_guides(ax, len(rows))
        ax.set_title(ttl, fontsize=8.5)
        ax.text(0.235, -0.5, "pulls more\nin August", ha="right", va="bottom", fontsize=7.2, color=MUTED)
        ax.text(-0.155, -0.5, "pulls less", ha="left", va="bottom", fontsize=7.2, color=MUTED)
    axa.set_yticks(range(len(rows))[::-1], [l for _, l in rows]); axa.tick_params(axis="y", length=0)
    axa2.tick_params(axis="y", length=0, labelleft=False)
    fig.text(0.5 * (axa.get_position().x0 + axa2.get_position().x1), axa.get_position().y0 - 0.048, XLAB,
             ha="center", va="top", fontsize=8.5)
    # b: the age test
    axb = fig.add_subplot(gs[1, 0])
    arows = [("pre1981_ge65", "Arrived before 1981,\naged 65+"), ("pre1981_lt65_age43_64", "Arrived before 1981,\naged 43–64"),
             ("post1981_ge65", "Arrived 1981 or later,\naged 65+"), ("recent2016_lt65", "Arrived 2016–2023,\nunder 65")]
    for k, (t, lab) in enumerate(arows):
        r = ag[ag.term == t].iloc[0]; yk = len(arows) - 1 - k
        axb.plot([r.ci_wcr_lo, r.ci_wcr_hi], [yk] * 2, color=INK, lw=1.4)
        axb.scatter([r.delta_aug_feb], [yk], s=26, color=INK, edgecolor="white", lw=0.8, zorder=3)
    axb.set_yticks(range(len(arows))[::-1], [l for _, l in arows]); axb.tick_params(axis="y", length=0, labelsize=7.5)
    axb.set_ylim(-0.6, len(arows) - 0.4); axb.set_xlim(-0.13, 0.15); zero_guides(axb, 4)
    axb.set_xlabel("Extra pull on 14–15 Aug\n(all destinations)")
    axb.set_title("Age or generation?", fontsize=8.5)
    # c: all 117 summer x ordinary night pairs
    axc = fig.add_subplot(gs[1, 1])
    crow = [("pre1981_exodus", "Arrived\nbefore 1981"), ("born_village_le2k", "Born in\na village"),
            ("recent_2016_2024", "Arrived\n2016–2023")]
    rng = np.random.default_rng(20261124)
    for k, (t, lab) in enumerate(crow):
        v = sp[sp.term == t].delta.values; yk = len(crow) - 1 - k
        axc.scatter(v, yk + rng.uniform(-0.18, 0.18, len(v)), s=5, color=INK, alpha=0.45, lw=0, rasterized=True)
        axc.plot([np.median(v)] * 2, [yk - 0.3, yk + 0.3], color=INK, lw=1.4)
        n_same = int((np.sign(v) == np.sign(np.median(v))).sum())
        axc.text(0.205, yk, f"{n_same}/{len(v)}", va="center", ha="right", fontsize=7.5, color=MUTED)
    axc.set_yticks(range(len(crow))[::-1], [l for _, l in crow]); axc.tick_params(axis="y", length=0, labelsize=7.5)
    axc.set_ylim(-0.6, len(crow) - 0.4); axc.set_xlim(-0.13, 0.21); zero_guides(axc, 3)
    axc.set_xticks([-0.1, 0, 0.1])
    axc.text(0.205, 2.33, "same sign", va="bottom", ha="right", fontsize=7.2, color=MUTED)
    axc.set_xlabel("Extra pull, each of 117 summer ×\nordinary night pairs (all destinations)")
    axc.set_title("Every pair of nights", fontsize=8.5)
    for ax, L, dx in ((axa, "a", 0.20), (axb, "b", 0.20), (axc, "c", 0.13)):  # a spans both top panels
        pos = ax.get_position()
        fig.text(pos.x0 - dx, pos.y1 + 0.02, L, fontsize=10.5, fontweight="bold", va="bottom")
    save(fig, "fig3_roots_generation")

def fig3b():
    """Slope-chart version: the pull of the home province on 13-14 Feb and 14-15 Aug."""
    mc = pd.read_csv(RES / "v1_c2_main_contrasts.csv").query("contrast == 'Aug15_vs_Feb14'")
    ag = pd.read_csv(RES / "v1_c2_age_confound.csv").query("spec == 'E2a_split65_4way_exclusive' and sample == 'all'")
    fig, axs = plt.subplots(1, 3, figsize=(TEXTWIDTH, 3.0), sharey=True, gridspec_kw=dict(wspace=0.95))
    def slope(ax, b0, b1, lo, hi, colour, lw, label, mk="o"):
        sig = lo > 0 or hi < 0
        ax.plot([0, 1], [b0, b1], color=colour, lw=lw, zorder=2, solid_capstyle="round")
        ax.scatter([0, 1], [b0, b1], s=22, color=colour, edgecolor="white", lw=0.7, zorder=3, marker=mk)
        d = b1 - b0
        txt = f"{label}\n{d:+.2f}".replace("-", "−") + ("" if sig else " (n.s.)")
        return (b1, txt)
    def finish(ax, labels, ylim, title, ylabel=None):
        # spread right-hand labels so they do not collide
        labels = sorted(labels, key=lambda t: t[0])
        ys, last = [], -1e9
        gap = (ylim[1] - ylim[0]) * 0.135
        for y, _ in labels:
            y = max(y, last + gap); ys.append(y); last = y
        for (y0, txt), y in zip(labels, ys):
            ax.annotate(txt, (1, y0), xytext=(1.16, y), textcoords="data", va="center", ha="left",
                        fontsize=7.2, linespacing=1.15, annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", color=GREY, lw=0.4, shrinkA=0, shrinkB=3))
        ax.set_xlim(-0.15, 1.08); ax.set_ylim(*ylim)
        ax.set_xticks([0, 1], ["13–14 Feb\nordinary", "14–15 Aug\nholiday"])
        ax.set_title(title, fontsize=8.5, pad=8)
        ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
        if ylabel: ax.set_ylabel(ylabel)
    groups = [("pre1981_exodus", "Exodus generation\n(arrived before 1981)", 1.8),
              ("born_village_le2k", "Born in a village", 1.2),
              ("recent_2016_2024", "Recent migrants\n(arrived 2016–23)", 0.9)]
    for ax, smp, col, ttl in ((axs[0], "inland", TERRA, "To inland destinations"),
                              (axs[1], "coastal", BLUE, "To coastal destinations")):
        labs = []
        for t, lab, lw in groups:
            r = mc[(mc["sample"] == smp) & (mc.term == t)].iloc[0]
            c_ = GREY if t.startswith("recent") else col
            labs.append(slope(ax, r.beta_ordinary, r.beta_holiday, r.ci_wcr_lo, r.ci_wcr_hi, c_, lw, lab.split("\n")[0]))
        finish(ax, labs, (0, 0.3), ttl, "Pull of the home province\n(elasticity)" if smp == "inland" else None)
    labs = []
    for t, lab, c_, lw in (("pre1981_ge65", "65+, arrived before 1981", INK, 1.8),
                           ("post1981_ge65", "65+, arrived after 1981", MUTED, 1.2),
                           ("recent2016_lt65", "Recent, under 65", GREY, 0.9)):
        r = ag[ag.term == t].iloc[0]
        labs.append(slope(axs[2], r.beta_feb14, r.beta_aug15, r.ci_wcr_lo, r.ci_wcr_hi, c_, lw, lab))
    finish(axs[2], labs, (0, 0.3), "Age or generation?\n(all destinations)")
    for ax in axs[1:]:
        ax.tick_params(axis="y", labelleft=False)
    for ax, L in zip(axs, "abc"):
        ax.text(-0.42 if L == "a" else -0.12, 1.1, L, transform=ax.transAxes, fontsize=10.5, fontweight="bold")
    save(fig, "fig3_roots_generation_slopes")

# ---------------------------------------------------------------- Figure 4
SEND, RECV, MID = "#6a3d9a", "#1b7837", "#efefef"   # diverging pair for net balance (not used elsewhere)

def _zones_geo():
    import geopandas as gpd
    g = gpd.read_file(ROOT / "clean_data/v2/zones/municipios_mitma.gpkg")[["id", "geometry"]]
    g = g[g.id.isin(nodes.id)].copy()
    from shapely.geometry import Polygon
    out = g.geometry.buffer(300).union_all().buffer(-300).simplify(600)
    main = max(getattr(out, "geoms", [out]), key=lambda q: q.area)   # mainland only, no holes
    g["geometry"] = g.geometry.simplify(400)
    return g, gpd.GeoSeries([Polygon(main.exterior)], crs=g.crs).boundary

def fig4():
    """a: maps, shaded where a zone receives more people than it sends; b: share of such zones by emptying.
    Greys only: blue/terracotta/ink are reserved for coast/inland/cities in every other figure."""
    from matplotlib.colors import ListedColormap, BoundaryNorm
    import matplotlib.cm as cm
    rn = pd.read_parquet(RES / "r1_node_metrics.parquet")
    rq = pd.read_csv(RES / "r1_net_balance_by_quintile.csv")
    g, outline = _zones_geo()
    bounds = [-1, 0, 0.1, 0.25, 0.5, 1]
    cols = ["#ffffff", "#e4e1db", "#bdb8af", "#8c867c", "#57524b"]
    cmap, norm = ListedColormap(cols), BoundaryNorm(bounds, len(cols))
    # (date, label, line colour, marker, line width, line style, marker face)
    nights = [("2024-02-14", "Ordinary weekday\n13–14 Feb", "#9a948a", "o", 1.3, (0, (4, 2)), "white"),
              ("2024-08-15", "Summer holiday\n14–15 Aug", INK, "s", 1.6, "-", "white"),
              ("2024-12-25", "Christmas Eve\n24–25 Dec", INK, "D", 1.6, "-", INK)]
    fig = plt.figure(figsize=(TEXTWIDTH, 5.9))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.42, wspace=0.02)
    for k, (d, ttl, *_ ) in enumerate(nights):
        ax = fig.add_subplot(gs[0, k])
        v = rn[rn.date == d][["id", "net_bal"]]
        gg = g.merge(v, on="id", how="left")
        gg.plot(column="net_bal", cmap=cmap, norm=norm, ax=ax, linewidth=0.15, edgecolor="face")  # vector map
        outline.plot(ax=ax, color=INK, lw=0.4)
        ax.set_axis_off(); ax.set_aspect("equal")
        ax.set_xlim(-30000, 1_150_000); ax.set_ylim(3_960_000, 4_880_000)
        share = (v.net_bal > 0).mean()
        ax.set_title(ttl, fontsize=8.5, linespacing=1.3)
        ax.text(0.5, -0.02, f"{share:.0%} of zones receive\nmore than they send", transform=ax.transAxes,
                ha="center", va="top", fontsize=7.5)
        if k == 0: axmap0 = ax
    from matplotlib.patches import Rectangle
    lax = fig.add_axes([0.10, 0.465, 0.84, 0.03]); lax.set_axis_off(); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
    lax.add_patch(Rectangle((0.0, 0.35), 0.035, 0.45, facecolor="white", edgecolor=INK, lw=0.4))
    lax.text(0.045, 0.57, "sends more people\nthan it receives", va="center", fontsize=7.2, linespacing=1.1)
    x0, w = 0.30, 0.095
    for i, c in enumerate(cols[1:]):
        lax.add_patch(Rectangle((x0 + i * w, 0.35), w, 0.45, facecolor=c, edgecolor=INK, lw=0.4))
    for i, t in enumerate(["0", "0.1", "0.25", "0.5", "1"]):
        lax.text(x0 + i * w, 0.22, t, ha="center", va="top", fontsize=7.2)
    lax.text(x0 + 4 * w + 0.015, 0.57, "receives more people\nthan it sends", va="center", fontsize=7.2, linespacing=1.1)
    lax.text(x0 + 2 * w, 1.35, "Net balance = (in − out) / (in + out) people sleeping away", ha="center",
             va="bottom", fontsize=7.5)
    axd = fig.add_subplot(gs[1, :])
    q = rq[(rq["sample"] == "all") & (rq.role == "main") & rq.cls.str.match(r"Q\d")]
    x = np.arange(5)
    for d, lab, col, mk, lw, ls, mfc in nights:
        sq = q[q.date == d].sort_values("cls")
        y = 100 * sq.share_net_importer.values
        axd.plot(x, y, color=col, lw=lw, ls=ls, marker=mk, ms=5.5, mfc=mfc, mec=col, mew=1.1, zorder=3)
        dy = {"2024-08-15": -6, "2024-12-25": 5}.get(d, 0)
        axd.text(4.3, y[-1] + dy, lab.replace("\n", ", "), va="center", fontsize=8,
                 color=col if d == "2024-02-14" else INK)
    axd.axhline(50, color=MUTED, lw=0.6, zorder=1)
    axd.text(4.3, 51.5, "above: most zones in the\ngroup receive more than they send", fontsize=7.2, color=MUTED, va="bottom")
    axd.text(4.3, 48.5, "below: most send more\nthan they receive", fontsize=7.2, color=MUTED, va="top")
    axd.set_xticks(x, ["least emptied\nfifth", "", "", "", "most emptied\nfifth"])
    axd.set_xlim(-0.2, 5.5); axd.set_ylim(0, 105)
    axd.set_xlabel("Zones grouped by population lost since their historical peak (1900–1960 to 2021)")
    axd.set_ylabel("Zones receiving more\npeople than they send (%)")
    axd.yaxis.grid(True, color=GRID, lw=0.5); axd.set_axisbelow(True)
    axd.spines["right"].set_visible(False); axd.spines["top"].set_visible(False)
    pa = axmap0.get_position()
    fig.text(0.07, pa.y1 + 0.005, "a", fontsize=10.5, fontweight="bold")
    pos = axd.get_position(); fig.text(0.07, pos.y1 + 0.02, "b", fontsize=10.5, fontweight="bold")
    save(fig, "fig4_villages")

# ---------------------------------------------------------------- Basic properties (section 3.2)
# Zhou & Mondragon (2007) Fig. 2 layout: P(k), knn(k), rich club, one curve per network, log-log.
# Same night styles as Figure 4: dashed grey = ordinary; black = holidays. Undirected degree, as Table 1.
DEG_NIGHTS = [("2024-02-14", "Ordinary weekday, 13–14 Feb", "#9a948a", "o", 1.3, (0, (4, 2)), "white"),
              ("2024-11-01", "All Saints, 31 Oct–1 Nov", "#6f6a62", "^", 1.4, "-", "#6f6a62"),
              ("2024-08-15", "Summer holiday, 14–15 Aug", INK, "s", 1.6, "-", "white")]

def _marked(ax, x, y, col, lw, ls, mk, mfc, n_marks=10, label=None):
    """Line plus markers at log-spaced x only, so that overlapping curves stay readable."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ax.plot(x, y, color=col, lw=lw, ls=ls, zorder=3, label=label)
    idx = np.unique(np.searchsorted(x, np.geomspace(x.min(), x.max(), n_marks)).clip(0, len(x) - 1))
    ax.plot(x[idx], y[idx], ls="", marker=mk, ms=4.2, mfc=mfc, mec=col, mew=1.0, zorder=4)

def figbasic():
    """a: CCDF of degree; b: knn(k) real vs rewired; c: rich club phi(r/N) real vs rewired (w10 backbone)."""
    cc = pd.read_parquet(RES / "b02_degree_ccdf.parquet").query("backbone == 'w10' and kind == 'und' and k >= 1")
    kn = pd.read_csv(RES / "b03_knn.csv", dtype={"date": str}).query("n_nodes >= 5")
    rc = pd.read_csv(RES / "b02_rich_club_zhou.csv", dtype={"date": str}).query("backbone == 'w10'")
    fig, axs = plt.subplots(1, 3, figsize=(TEXTWIDTH, 2.45), gridspec_kw=dict(wspace=0.62))
    axa, axb, axc = axs
    for d, lab, col, mk, lw, ls, mfc in DEG_NIGHTS:
        s = cc[cc.date == d].sort_values("k")
        _marked(axa, s.k, s.ccdf, col, lw, ls, mk, mfc)
        real = kn[(kn.date == d) & (kn.series == "real")].sort_values("k_lo")
        rew = (kn[(kn.date == d) & (kn.series == "rewired")].groupby("k_lo")
               .agg(k=("k_center", "mean"), knn=("knn", "mean")).reset_index().sort_values("k_lo"))
        axb.plot(rew.k, rew.knn, color=col, lw=0.8, ls=(0, (1, 1.2)), zorder=2)
        _marked(axb, real.k_center, real.knn, col, lw, ls, mk, mfc, n_marks=8)
        r = rc[rc.date == d].sort_values("r_norm")
        axc.plot(r.r_norm, r.phi_rew_mean, color=col, lw=0.8, ls=(0, (1, 1.2)), zorder=2)
        _marked(axc, r.r_norm, r.phi_real, col, lw, ls, mk, mfc, n_marks=8)
    axa.axvline(len(nodes) - 1, color=MUTED, lw=0.6, ls=(0, (1, 1.5)), zorder=1)
    axa.text(len(nodes) * 0.85, 0.75, "linked to\nall other\nzones", ha="right", va="top", fontsize=7.2, color=MUTED)
    axa.set_xlim(0.8, 4000); axa.set_ylim(2.5e-4, 1.6)
    axa.set_xlabel("Degree $k$"); axa.set_ylabel("Share of zones with degree ≥ $k$")
    axa.set_title("Degree distribution", fontsize=8.5)
    axb.set_xlim(4, 3000); axb.set_ylim(20, 1000)
    axb.set_xlabel("Degree $k$"); axb.set_ylabel("Neighbours' mean degree $k_{nn}$")
    axb.set_title("Mixing", fontsize=8.5)
    axc.set_xlim(5e-4, 1.2); axc.set_ylim(5e-3, 1.4)
    axc.set_xlabel("Rank $r/N$ (best connected first)")
    axc.set_ylabel("Rich-club connectivity $\\phi$")
    axc.set_title("Rich club", fontsize=8.5)
    for ax in axs:
        ax.set_xscale("log"); ax.set_yscale("log")
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
    axb.yaxis.set_major_locator(FixedLocator([30, 100, 300, 1000]))
    axb.yaxis.set_minor_formatter(NullFormatter())
    axb.set_yticklabels(["30", "100", "300", "1,000"])
    from matplotlib.lines import Line2D
    hh = [Line2D([], [], color=col, lw=lw, ls=ls, marker=mk, ms=4.2, mfc=mfc, mec=col, label=lab)
          for d, lab, col, mk, lw, ls, mfc in DEG_NIGHTS]
    hh.append(Line2D([], [], color=INK, lw=0.8, ls=(0, (1, 1.2)), label="same degrees, randomly rewired"))
    fig.legend(handles=hh, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.99), handlelength=2.6,
               columnspacing=1.6)
    for ax, L in zip(axs, "abc"):
        panel_letter(ax, L, x=-0.22, y=1.04)
    save(fig, "fig_basic_properties")

def figbosam():
    """BOSAM (Guo, Chen & Zhou 2008), three main nights, w10 backbone, undirected projection.
    Top row: whole matrix; with 2,435 zones a pixel cannot hold one link at print size, so each grey level is
    the share of links present in a block of BLOCK x BLOCK cells. Bottom row: the 5 % best-connected zones
    (Zhou & Mondragon's top 5 %), one cell per pair, with a strip giving the kind of zone (project colours)."""
    from matplotlib.colors import ListedColormap, LogNorm
    from matplotlib.patches import Rectangle
    nights = [("ordinary", "2024-02-14", "Ordinary, 13–14 Feb"),
              ("allsaints", "2024-11-01", "All Saints, 31 Oct–1 Nov"),
              ("summer", "2024-08-15", "Summer, 14–15 Aug")]
    n = len(nodes); BLOCK = 25; top = int(np.floor(0.05 * n))
    fig = plt.figure(figsize=(TEXTWIDTH, 5.6))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.07], hspace=0.62, wspace=0.10, bottom=0.13)
    kind_cmap = ListedColormap([TERRA, BLUE, INK])
    pairs = []
    for c, (g, d, ttl) in enumerate(nights):
        o = pd.read_csv(RES / f"b02_bosam_order_{g}_{d}__w10.csv")
        e = pd.read_parquet(RES / f"b02_bosam_edges_{g}_{d}__w10.parquet")
        A = np.zeros((n, n), dtype=np.uint8)
        A[e.row.values, e.col.values] = 1; A[e.col.values, e.row.values] = 1
        m = (n // BLOCK) * BLOCK
        D = A[n - m:, n - m:].reshape(m // BLOCK, BLOCK, m // BLOCK, BLOCK).mean(axis=(1, 3))
        ax = fig.add_subplot(gs[0, c])
        im = ax.imshow(np.clip(D, 1e-3, 1), origin="lower", cmap="Greys", norm=LogNorm(1e-3, 1),
                       interpolation="nearest", extent=(n - m, n, n - m, n), rasterized=True)
        import matplotlib.patheffects as pe
        ax.add_patch(Rectangle((n - top, n - top), top, top, fill=False, ec="white", lw=0.9,
                               path_effects=[pe.Stroke(linewidth=2.0, foreground=INK), pe.Normal()]))
        ax.set_xlim(0, n); ax.set_ylim(0, n); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([0, n], ["1", f"{n:,}"]) if c == 0 else ax.set_yticks([])
        ax.tick_params(axis="y", labelsize=7.5)
        ax.set_title(ttl, fontsize=8.5)
        phi = A[n - top:, n - top:].sum() / (top * (top - 1))
        # zoom: top 5 %, one cell per pair
        az = fig.add_subplot(gs[1, c])
        az.imshow(A[n - top:, n - top:], origin="lower", cmap="Greys", vmin=0, vmax=1, interpolation="nearest",
                  extent=(n - top, n, n - top, n), rasterized=True)
        az.set_xticks([]); az.set_yticks([n - top, n], [f"{n - top + 1:,}", f"{n:,}"]) if c == 0 else az.set_yticks([])
        az.tick_params(axis="y", labelsize=7.5)
        az.set_title(f"{100 * phi:.0f} % of pairs linked", fontsize=8)
        idx = o.idx.values[n - top:]
        kind = np.where(CITY[idx], 2, np.where(COAST[idx], 1, 0))[None, :]
        axs = fig.add_subplot(gs[2, c])
        axs.imshow(kind, cmap=kind_cmap, vmin=0, vmax=2, aspect="auto", interpolation="nearest",
                   extent=(n - top, n, 0, 1), rasterized=True)
        axs.set_yticks([]); axs.set_xticks([])
        if c == 1:
            ax.set_xlabel("2,425 of 2,435 zones in BOSAM order (degree, then neighbours), least connected first", labelpad=2)
            axs.set_xlabel(f"Its {top} best-connected zones (outlined square above), same order", labelpad=3)
        panel_letter(ax, "abc"[c], x=-0.04 if c else -0.22, y=1.13)
        panel_letter(az, "def"[c], x=-0.04 if c else -0.22, y=1.09)
        if c == 2:
            top_ax = ax
        pairs.append((az, axs))
    fig.canvas.draw()  # square zooms shrink to their aspect; put each strip right under its zoom
    for az, axs in pairs:
        pz = az.get_position()
        axs.set_position([pz.x0, pz.y0 - 0.050, pz.width, 0.026])
    pt = top_ax.get_position()
    cax = fig.add_axes([pt.x1 + 0.012, pt.y0, 0.012, pt.height])
    cb = fig.colorbar(im, cax=cax, ticks=[1e-3, 1e-2, 1e-1, 1])
    cb.ax.set_yticklabels(["0.1 %", "1 %", "10 %", "100 %"]); cb.ax.tick_params(labelsize=7.2, length=2)
    cb.set_label("Pairs linked in each\n25 × 25 block", fontsize=7.5); cb.outline.set_linewidth(0.5)
    from matplotlib.patches import Patch
    hh = [Patch(color=TERRA, label="Inland towns and villages"), Patch(color=BLUE, label="Coastal towns and villages"),
          Patch(color=INK, label="Cities over 100,000 people")]
    fig.legend(handles=hh, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.975), handlelength=1.2,
               columnspacing=1.6)
    save(fig, "fig_bosam")

FIGS = {"fig1": fig1, "fig1b": fig1b, "fig2": fig2, "fig3": fig3, "fig3b": fig3b, "fig4": fig4,
        "figbasic": figbasic, "figbosam": figbosam}
if __name__ == "__main__":
    for name in (sys.argv[1:] or FIGS):
        FIGS[name]()

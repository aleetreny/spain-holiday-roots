"""Figure 1 options: forms that show quantity, in two visual themes. Drafts."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "figures_drafts" / "fig1_options"; OUT.mkdir(parents=True, exist_ok=True)
d1 = pd.read_csv(ROOT / "reports/figures/fig1_data.csv")
SHORT = ["Ordinary weekday\n13–14 Feb", "Easter\n28–29 Mar", "Summer holiday\n14–15 Aug",
         "All Saints\n31 Oct–1 Nov", "Christmas Eve\n24–25 Dec"]

THEMES = {
  "A": dict(font="STIX Two Text", size=9, towns="#c0623b", cities="#d9d2c5", coast="#1f6e8c",
            ink="#1d1d1b", muted="#6b6862", grid="#e6e2db"),
  "B": dict(font="Helvetica", size=7.5, towns="#d08c2f", cities="#cfcfcf", coast="#2b6f9e",
            ink="#222222", muted="#666666", grid="#ebebeb"),
}
def setup(t):
    T = THEMES[t]
    plt.rcParams.update({"font.family": T["font"], "font.size": T["size"], "axes.spines.top": False,
        "axes.spines.right": False, "axes.edgecolor": T["muted"], "axes.labelcolor": T["ink"],
        "xtick.color": T["muted"], "ytick.color": T["muted"], "text.color": T["ink"], "savefig.dpi": 300,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    return T
CATS = [("towns", "inland towns and villages"), ("cities", "inland cities (>100,000)"), ("coast", "coast")]

def stacked_abs(ax, T, label_min=0.35):
    y = np.arange(len(d1))[::-1]; left = np.zeros(len(d1))
    for c, lab in CATS:
        ax.barh(y, d1[c] - 0.02, left=left, height=0.62, color=T[c], label=lab)
        for yi, l, v in zip(y, left, d1[c]):
            if v >= label_min:
                ax.text(l + v / 2, yi, f"{v:.1f}", ha="center", va="center",
                        color="white" if c != "cities" else T["ink"])
        left += d1[c].values
    for yi, tot in zip(y, d1.total):
        ax.text(tot + 0.06, yi, f"{tot:.1f} M", va="center", color=T["ink"])
    ax.set_yticks(y, SHORT); ax.tick_params(axis="y", length=0, labelcolor=T["ink"])
    ax.spines["left"].set_visible(False); ax.xaxis.grid(True, color=T["grid"], lw=0.6); ax.set_axisbelow(True)
    ax.set_xlim(0, 5.2); ax.set_xlabel("People sleeping ≥100 km from home (millions)".replace("≥", "at least "))

def legend_top(ax, T, ncol=3):
    ax.legend(frameon=False, ncol=ncol, loc="lower left", bbox_to_anchor=(-0.02, 1.0), handlelength=1.2,
              columnspacing=1.2)

# Option 1: absolute stacked bars (theme A and B)
for t in "AB":
    T = setup(t)
    fig, ax = plt.subplots(figsize=(6.3, 2.7)); stacked_abs(ax, T); legend_top(ax, T)
    fig.savefig(OUT / f"opt1_absolute_bars_theme{t}.png", bbox_inches="tight"); plt.close(fig)

# Option 2: dot plot, three dots per night (theme B)
T = setup("B")
fig, ax = plt.subplots(figsize=(6.3, 2.6)); y = np.arange(len(d1))[::-1]
for yi, r in zip(y, d1.itertuples()):
    vals = [r.towns, r.cities, r.coast]
    ax.plot([min(vals), max(vals)], [yi, yi], color=T["grid"], lw=2.5, zorder=1, solid_capstyle="round")
for c, lab in CATS:
    ax.scatter(d1[c], y, s=46, color=T[c], edgecolor="white", lw=1, zorder=3, label=lab)
ax.set_yticks(y, SHORT); ax.tick_params(axis="y", length=0, labelcolor=T["ink"]); ax.spines["left"].set_visible(False)
ax.xaxis.grid(True, color=T["grid"], lw=0.6); ax.set_axisbelow(True); ax.set_xlim(0, 2.4)
ax.set_xlabel("People sleeping at least 100 km from home (millions)"); legend_top(ax, T)
fig.savefig(OUT / "opt2_dotplot_themeB.png", bbox_inches="tight"); plt.close(fig)

# Option 3: the calendar, every night with data in 2024 (theme B)
dm = pd.read_csv(ROOT / "results/p_daily_metrics.csv"); dm = dm[dm.date.str.startswith("2024")].copy()
dm["date"] = pd.to_datetime(dm.date)
dm["coast"] = dm.away_persons * dm.share_dest_coastal / 1e6
dm["inland"] = dm.away_persons * (1 - dm.share_dest_coastal) / 1e6
blocks = [("2024-02-01", "2024-03-31", "February–March"), ("2024-08-01", "2024-08-31", "August"),
          ("2024-10-01", "2024-12-31", "October–December")]
widths = [60, 31, 92]
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.5), sharey=True, gridspec_kw=dict(width_ratios=widths, wspace=0.06))
marks = {"2024-02-14": "ordinary\nweekday", "2024-03-29": "Easter", "2024-08-15": "15 Aug",
         "2024-11-01": "All Saints", "2024-12-25": "Christmas\nEve"}
for ax, (a, b, ttl) in zip(axs, blocks):
    s = dm[(dm.date >= a) & (dm.date <= b)]
    ax.plot(s.date, s.inland, color=T["towns"], lw=1.4, label="slept inland")
    ax.plot(s.date, s.coast, color=T["coast"], lw=1.4, label="slept at the coast")
    for dstr, lab in marks.items():
        dd = pd.Timestamp(dstr)
        if pd.Timestamp(a) <= dd <= pd.Timestamp(b):
            v = s.loc[s.date == dd, "inland"].iloc[0]
            ax.annotate(lab, (dd, v), xytext=(0, 14), textcoords="offset points", ha="center", fontsize=6.5,
                        color=T["ink"], arrowprops=dict(arrowstyle="-", color=T["muted"], lw=0.5))
    ax.set_title(ttl, loc="left", fontsize=7.5, color=T["muted"])
    ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.yaxis.grid(True, color=T["grid"], lw=0.6); ax.set_axisbelow(True)
for ax in axs[1:]:
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
axs[0].set_ylabel("People sleeping away\nfrom home (millions)"); axs[0].set_ylim(0, 6)
axs[0].legend(frameon=False, loc="upper left")
fig.savefig(OUT / "opt3_calendar_themeB.png", bbox_inches="tight"); plt.close(fig)

# Option 4: two panels, a = how many, b = where (shares), theme B, Nature-style letters
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.7), sharey=True, gridspec_kw=dict(width_ratios=[1, 1.15], wspace=0.08))
y = np.arange(len(d1))[::-1]
axs[0].barh(y, d1.total, height=0.62, color=T["muted"])
for yi, tot in zip(y, d1.total): axs[0].text(tot + 0.08, yi, f"{tot:.1f} M", va="center")
axs[0].set_xlim(0, 5.3); axs[0].set_xlabel("People sleeping at least 100 km\nfrom home (millions)")
sh = d1[["towns", "cities", "coast"]].div(d1.total, axis=0); left = np.zeros(len(d1))
for c, lab in CATS:
    axs[1].barh(y, sh[c] - 0.004, left=left, height=0.62, color=T[c], label=lab)
    for yi, l, v in zip(y, left, sh[c]):
        if v >= 0.07: axs[1].text(l + v / 2, yi, f"{v:.0%}", ha="center", va="center", color="white" if c != "cities" else T["ink"])
    left += sh[c].values
axs[1].set_xlim(0, 1); axs[1].set_xticks([0, .5, 1], ["0%", "50%", "100%"]); axs[1].set_xlabel("Where they slept (share)")
for ax, L in zip(axs, "ab"):
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    ax.xaxis.grid(True, color=T["grid"], lw=0.6); ax.set_axisbelow(True)
    ax.text(-0.02 if L == "b" else -0.42, 1.04, L, transform=ax.transAxes, fontweight="bold", fontsize=9)
axs[0].set_yticks(y, SHORT); axs[0].tick_params(axis="y", labelcolor=T["ink"])
axs[1].legend(frameon=False, ncol=3, loc="lower left", bbox_to_anchor=(-0.02, 1.06), handlelength=1.2, columnspacing=1.0)
fig.savefig(OUT / "opt4_two_panels_themeB.png", bbox_inches="tight"); plt.close(fig)
print("ok")

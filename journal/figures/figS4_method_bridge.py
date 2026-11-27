"""Paired estimates on the six nights recalculated under both MITMA methods."""

import pandas as pd

from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from style import A2, INK, MARK, ORDINARY, decimals, figure, panel, place, read, save


def night_label(date):
    day = (date - pd.Timedelta(days=1)).day
    return f"{day}–{date.day} {date:%b %Y}"


def main():
    data = read(A2, "method_bridge_dates").sort_values("file_date")
    assert len(data) == 6
    fig = figure(60)
    for i, (term, title) in enumerate([("pooled", "Pooled roots"), ("pre1981", "Arrived before 1981")]):
        ax = place(fig, 30 + i * 86, 10, 52, 33)
        panel(ax, "ab"[i], title, at=2 + i * 86)
        for y, row in enumerate(data.to_dict("records")):
            old, new = row[f"old_{term}"], row[f"new_{term}"]
            ax.plot([old, new], [y - .10, y + .10], color=ORDINARY, linewidth=.85, zorder=1)
            # Stagger endpoints on the same date to preserve near-identical estimates.
            ax.plot(old, y - .10, "o", color=INK, mfc="white", mew=.75, ms=MARK, zorder=3)
            ax.plot(new, y + .10, "s", color=INK, ms=MARK - .4, zorder=4)
        labels = [night_label(date) for date in pd.to_datetime(data.file_date)]
        ax.set_yticks(range(len(data)), labels)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.set_ylim(len(data) - .45, -.55)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.xaxis.set_major_formatter(decimals(2 if i == 0 else 3))
        ax.set_xlabel("Roots elasticity")
    fig.legend(handles=[
        Line2D([], [], ls="", marker="o", ms=MARK, mfc="white", mec=INK, label="Old method"),
        Line2D([], [], ls="", marker="s", ms=MARK - .4, color=INK, label="New method"),
    ], loc="upper center", bbox_to_anchor=(.5, 1 - 54 / 60), ncol=2)
    save(fig, "S4", data)


if __name__ == "__main__":
    main()

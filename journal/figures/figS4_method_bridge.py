"""Paired estimates on the six nights recalculated under both MITMA methods."""

import pandas as pd

from style import A2, MUTED, NEUTRAL, figure, panel, read, save


def night_label(date):
    day = (date - pd.Timedelta(days=1)).day
    return f"{day}–{date.day} {date:%b %Y}"


def main():
    data = read(A2, "method_bridge_dates").sort_values("file_date")
    assert len(data) == 6
    fig = figure(63)
    for i, (term, title) in enumerate([("pooled", "Pooled roots"), ("pre1981", "Arrived before 1981")]):
        ax = fig.add_axes([.185 + i * .425, .205, .365, .626])
        panel(ax, "ab"[i], title, x=-.11)
        for y, row in enumerate(data.to_dict("records")):
            old, new = row[f"old_{term}"], row[f"new_{term}"]
            ax.plot([old, new], [y - .10, y + .10], color="#999999", linewidth=.85, zorder=1)
            # Stagger endpoints on the same date to preserve near-identical estimates.
            ax.plot(old, y - .10, "o", color=NEUTRAL, mfc="white", mew=.75, ms=4, zorder=3)
            ax.plot(new, y + .10, "s", color=NEUTRAL, ms=3.6, zorder=4)
        labels = [night_label(date) for date in pd.to_datetime(data.file_date)]
        ax.set_yticks(range(len(data)), labels if i == 0 else [""] * len(data))
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.set_ylim(len(data) - .45, -.55)
        ax.set_xlabel("Roots elasticity")
    fig.text(.58, .034, "○  Old method     ■  New method", ha="center", color=MUTED, fontsize=7)
    save(fig, "S4", data)


if __name__ == "__main__":
    main()

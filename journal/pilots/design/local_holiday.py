"""Audit whether Castilla y León local holidays support a nightly event study."""
import re
import sys
import unicodedata
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd

from pilot_data import OUT, birth_ties, municipalities

HOLIDAYS = Path("/Users/alejandrotreny/Proyectos/spain-holiday-roots-journal/_local/external_journal/jcyl_fiestas/fiestas.csv")
# 2024 BOE national / Castilla y León regional holidays and 29–31 October DANA disruption.
CONFLICTS = pd.to_datetime([
    "2024-01-01", "2024-01-06", "2024-03-28", "2024-03-29", "2024-04-23", "2024-05-01",
    "2024-08-15", "2024-10-12", "2024-10-29", "2024-10-30", "2024-10-31", "2024-11-01",
    "2024-12-06", "2024-12-08", "2024-12-09", "2024-12-25",
])


def name_key(name):
    text = unicodedata.normalize("NFKD", name.upper()).encode("ascii", "ignore").decode()
    return " ".join(sorted(re.findall(r"[A-Z]+", text)))


def reconcile(calendar, villages):
    events = calendar.loc[calendar["event_date"].dt.year.eq(2024) & calendar["dest"].isin(villages.index)].copy()
    # Some official rows give the name of a submunicipal nucleus for the same INE code.
    events = events.loc[events["municipio"].map(name_key).eq(events["dest"].map(villages["LAU_NAME"].map(name_key)))]
    return events[["dest", "event_date"]].drop_duplicates()


def clean_windows(events, offsets):
    windows = pd.DataFrame({shift: events["event_date"] + pd.Timedelta(days=shift) for shift in offsets})
    local = pd.MultiIndex.from_frame(events)
    good = pd.Series(True, index=events.index)
    for night in windows:
        for side in [-1, 0, 1]:
            adjacent = windows[night] + pd.Timedelta(days=side)
            good &= ~adjacent.isin(CONFLICTS)
            overlap = pd.MultiIndex.from_arrays([events["dest"], adjacent]).isin(local)
            good &= ~(overlap & adjacent.ne(events["event_date"]))
    return windows.loc[good]


def main():
    calendar = pd.read_csv(HOLIDAYS, sep=";", dtype={"ine": str}, parse_dates=["fecha_fiesta"])
    calendar = calendar.rename(columns={"fecha_fiesta": "event_date"})
    calendar["dest"] = calendar["ine"].str.zfill(5)
    calendar["year"] = calendar["event_date"].dt.year.astype("Int64").astype("string").fillna("missing_date")
    calendar.groupby("year").agg(raw_rows=("event_date", "size"), distinct_ine_codes=("dest", "nunique")).to_csv(
        OUT / "jcyl_calendar_coverage.csv"
    )
    daily = pd.read_csv(OUT / "mitma_daily_inventory.csv", parse_dates=["file_date"])
    published = daily.loc[daily["standard_listed"], "file_date"]
    local = daily.loc[daily["derived_local"] | daily["raw_local"], "file_date"]

    rows = []
    for cutoff in [2000, 5000]:
        cities, villages = municipalities(cutoff)
        villages = villages.loc[villages["n_municipalities"].eq(1) & villages.index.to_series().eq(villages["zone"])]
        events = reconcile(calendar, villages)
        ties = birth_ties(cities.index, events["dest"].unique()).reset_index()
        ties = ties.loc[ties["orig"].str[:2].ne(ties["dest"].str[:2]) & ties["g2_households"].gt(0)]
        tied_villages = set(ties["dest"])
        for label, offsets in [
            ("main_plus_minus_7", [0, 1, -7, -6, 7, 8]),
            ("main_and_lead_placebo", [0, 1, -7, -6, 7, 8, -14, -13]),
        ]:
            windows = clean_windows(events, offsets)
            for coverage, kept in [
                ("matched_names_and_zones", events.index),
                ("no_other_holiday_or_DANA", windows.index),
                ("complete_rss_windows", windows.index[windows.isin(published.tolist()).all(axis=1)]),
                ("complete_local_windows", windows.index[windows.isin(local.tolist()).all(axis=1)]),
            ]:
                subset = events.loc[kept]
                exposed = subset["dest"].isin(tied_villages)
                needed = set(windows.loc[windows.index.intersection(kept)].stack()) - set(local)
                rows.append({
                    "population_max": cutoff, "window": label, "coverage": coverage,
                    "events": len(subset), "villages": subset["dest"].nunique(),
                    "provinces": subset["dest"].str[:2].nunique(),
                    "events_with_strict_parent_tie": int(exposed.sum()),
                    "distinct_window_days_not_local": len(needed),
                })
    pd.DataFrame(rows).to_csv(OUT / "local_holiday_feasibility.csv", index=False)


if __name__ == "__main__":
    main()

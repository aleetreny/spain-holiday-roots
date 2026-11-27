"""a06_descriptive: Descriptive calendar extras (clearly labelled descriptive, no formal tests).

Answers Task 6 & PREREGISTRATION.md §5 descriptive outputs:
  - Highest and lowest beta_t nights of 2022-2024 (pooled all, pooled inland, pre1981, village-born)
  - July-August trajectory (week by week and day by day: is 15 Aug special or all of August?)
  - Weekends (Fri/Sat/Sun evenings) vs Mon-Wed ordinary vs Thursday evenings
  - Long weekends (puentes: Constitución/Inmaculada Dec 6-8, 1 May, 12 Oct, 1 Nov, etc.)
  - Christmas & New Year window (20 Dec - 7 Jan: Christmas Eve vs New Year's Eve vs Epiphany)
  - Emptied-village receiving-end daily summary (Q5 receiver share, arrivals & departures per resident)

Run:
  /Users/alejandrotreny/Proyectos/spain-holiday-roots/.venv/bin/python -I journal/analysis/A_calendar/a06_descriptive.py
"""
from __future__ import annotations

import datetime as dt
import sys
import time
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a00_paths import B0, C, OUT_DIR, PANEL_DIR, ensure_dirs, make_logger  # noqa: E402
from a01_calendar import load_target_holiday_dates  # noqa: E402


def tag_descriptive_period(df: pd.DataFrame, hols: set[dt.date]) -> pd.Series:
    """Assign mutually exclusive descriptive calendar buckets (evening-date based)."""
    labels = []
    for row in df.itertuples(index=False):
        d = str(row.date)
        mmdd = d[5:10]
        e = dt.date.fromisoformat(str(row.evening_date))
        e_mmdd = e.isoformat()[5:10]
        e_dow = e.weekday()  # 0=Mon .. 6=Sun
        nc = row.night_class

        if nc == "AUG":
            labels.append("AUG_peak_0810_0818")
        elif nc == "XMAS":
            labels.append("XMAS_eve_1225")
        elif nc == "EASTER":
            labels.append("EASTER_good_friday")
        elif nc == "ALLSAINTS":
            labels.append("ALLSAINTS_1101")
        elif nc == "O":
            labels.append("O_ordinary_MonTueWed")
        elif mmdd == "01-01":
            labels.append("NEW_YEAR_eve_0101")
        elif mmdd == "01-06":
            labels.append("EPIPHANY_eve_0106")
        elif e_mmdd >= "12-20" or e_mmdd <= "01-07":
            labels.append("XMAS_NEWYEAR_other")
        elif "07-01" <= e_mmdd <= "07-31":
            labels.append("JULY_all")
        elif "08-01" <= e_mmdd <= "08-08":
            labels.append("AUG_early_0801_0808")
        elif "08-18" <= e_mmdd <= "08-31":
            labels.append("AUG_late_0818_0831")
        elif "12-05" <= e_mmdd <= "12-09":
            labels.append("PUENTE_Dec_Constitucion")
        elif ("10-10" <= e_mmdd <= "10-14") or ("04-29" <= e_mmdd <= "05-03") or ("11-01" <= e_mmdd <= "11-03"):
            labels.append("PUENTE_other_national")
        elif e in hols or (e + dt.timedelta(days=1)) in hols:
            labels.append("OTHER_holiday_or_eve")
        elif e_dow == 3:
            labels.append("NONHOL_Thursday_eve")
        elif e_dow == 4:
            labels.append("WEEKEND_Friday_eve")
        elif e_dow == 5:
            labels.append("WEEKEND_Saturday_eve")
        elif e_dow == 6:
            labels.append("WEEKEND_Sunday_eve")
        else:
            labels.append("OTHER_shoulder_or_holyweek")
    return pd.Series(labels, index=df.index)


def main() -> None:
    t0 = time.time()
    log = make_logger(t0)
    ensure_dirs()

    man = pd.read_csv(OUT_DIR / "night_manifest.csv")
    np_df = pd.read_csv(OUT_DIR / "nightly_pooled.csv")
    nc_df = pd.read_csv(OUT_DIR / "nightly_cohorts.csv")
    hols = load_target_holiday_dates()

    # Build wide nightly summary joining pooled (all, inland, coastal) + cohort pre1981/recent + village-born
    p_all = np_df[np_df["sample"] == "all"][["date", "evening_date", "year", "evening_weekday", "split", "night_class", "beta", "se_crv1", "ci_lo", "ci_hi", "flow_total"]].rename(
        columns={"beta": "beta_pooled_all", "se_crv1": "se_pooled_all", "ci_lo": "ci_lo_all", "ci_hi": "ci_hi_all", "flow_total": "flow_112_all"}
    )
    p_inl = np_df[np_df["sample"] == "inland"][["date", "beta"]].rename(columns={"beta": "beta_pooled_inland"})
    p_cst = np_df[np_df["sample"] == "coastal"][["date", "beta"]].rename(columns={"beta": "beta_pooled_coastal"})
    c_pre = nc_df[(nc_df["sample"] == "all") & (nc_df["model"] == "cohort_joint") & (nc_df["term"] == "pre1981")][["date", "beta"]].rename(columns={"beta": "beta_pre1981_all"})
    c_rec = nc_df[(nc_df["sample"] == "all") & (nc_df["model"] == "cohort_joint") & (nc_df["term"] == "recent_2016_2023")][["date", "beta"]].rename(columns={"beta": "beta_recent_all"})
    c_vil = nc_df[(nc_df["sample"] == "all") & (nc_df["model"] == "birthplace_joint") & (nc_df["term"] == "born_le_2000")][["date", "beta"]].rename(columns={"beta": "beta_vil2k_all"})
    c_pre_i = nc_df[(nc_df["sample"] == "inland") & (nc_df["model"] == "cohort_joint") & (nc_df["term"] == "pre1981")][["date", "beta"]].rename(columns={"beta": "beta_pre1981_inland"})

    wide = p_all.merge(p_inl, on="date").merge(p_cst, on="date").merge(c_pre, on="date").merge(c_rec, on="date").merge(c_vil, on="date").merge(c_pre_i, on="date")
    wide["desc_bucket"] = tag_descriptive_period(wide, hols)
    wide["month"] = wide["date"].str[5:7].astype(int)

    # 1. Top 25 & bottom 15 nights by pooled beta_all and by pre1981_all
    top_p = wide.sort_values("beta_pooled_all", ascending=False).head(25).assign(ranking="top25_pooled_all")
    bot_p = wide.sort_values("beta_pooled_all", ascending=True).head(15).assign(ranking="bot15_pooled_all")
    top_pre = wide.sort_values("beta_pre1981_all", ascending=False).head(25).assign(ranking="top25_pre1981_all")
    top_vil = wide.sort_values("beta_vil2k_all", ascending=False).head(25).assign(ranking="top25_vil2k_all")
    pd.concat([top_p, bot_p, top_pre, top_vil], ignore_index=True).to_csv(OUT_DIR / "descriptive_top_nights.csv", index=False)

    # 2. Summary by descriptive bucket (overall + by year)
    agg_cols = ["beta_pooled_all", "beta_pooled_inland", "beta_pooled_coastal", "beta_pre1981_all", "beta_pre1981_inland", "beta_recent_all", "beta_vil2k_all", "flow_112_all"]
    bucket_rows = []
    for bkt, g in wide.groupby("desc_bucket", sort=True):
        rec = {"desc_bucket": bkt, "n_nights": len(g)}
        for c in agg_cols:
            rec[f"{c}_mean"] = float(g[c].mean())
            rec[f"{c}_median"] = float(g[c].median())
            rec[f"{c}_p10"] = float(g[c].quantile(0.10))
            rec[f"{c}_p90"] = float(g[c].quantile(0.90))
        bucket_rows.append(rec)
    df_bkt = pd.DataFrame(bucket_rows).sort_values("beta_pooled_all_mean", ascending=False)
    df_bkt.to_csv(OUT_DIR / "descriptive_calendar_summary.csv", index=False)

    # 3. July-August weekly & day-of-summer trajectory (June 15 .. Sep 15)
    ja = wide[(wide["date"].str[5:10] >= "06-15") & (wide["date"].str[5:10] <= "09-15")].copy()
    ja["mmdd"] = ja["date"].str[5:10]
    ja_day = ja.groupby("mmdd")[["beta_pooled_all", "beta_pooled_inland", "beta_pooled_coastal", "beta_pre1981_all", "beta_vil2k_all", "flow_112_all"]].agg(["count", "mean", "min", "max"])
    ja_day.columns = [f"{a}_{b}" for a, b in ja_day.columns]
    ja_day.reset_index().to_csv(OUT_DIR / "descriptive_summer_trajectory.csv", index=False)

    # 4. Monthly mean/median across all 12 months (all valid nights vs O-only nights)
    m_all = wide.groupby("month")[["beta_pooled_all", "beta_pooled_inland", "beta_pre1981_all", "beta_recent_all", "beta_vil2k_all"]].agg(["count", "mean", "median"])
    m_all.columns = [f"{a}_{b}" for a, b in m_all.columns]
    m_all.reset_index().assign(subset="all_valid").to_csv(OUT_DIR / "descriptive_monthly_summary.csv", index=False)
    m_o = wide[wide["night_class"] == "O"].groupby("month")[["beta_pooled_all", "beta_pooled_inland", "beta_pre1981_all", "beta_recent_all", "beta_vil2k_all"]].agg(["count", "mean", "median"])
    m_o.columns = [f"{a}_{b}" for a, b in m_o.columns]
    pd.concat([
        m_all.reset_index().assign(subset="all_valid"),
        m_o.reset_index().assign(subset="O_only"),
    ], ignore_index=True).to_csv(OUT_DIR / "descriptive_monthly_summary.csv", index=False)

    # 5. Emptied-village receiving-end daily summary (PREREGISTRATION.md §5 descriptive output)
    nodes = B0.load_nodes()
    elig = nodes["loss_peak_2021"].notna() & (~nodes["id"].isin(C.ARTEFACT_ZONES))
    q_series = pd.Series(np.nan, index=nodes.index)
    q_series.loc[elig] = pd.qcut(nodes.loc[elig, "loss_peak_2021"], 5, labels=False) + 1
    nodes["loss_q"] = q_series
    q5_ids = nodes.loc[nodes["loss_q"] == 5, ["id", "population"]].copy()
    assert len(q5_ids) == 475, f"Expected 475 Q5 zones, got {len(q5_ids)}"

    con = duckdb.connect()
    con.register("q5", q5_ids)
    pz = PANEL_DIR / "panel_zones2435_daily.parquet"
    q5_daily = con.execute(f"""
        SELECT z.date,
               count(*) AS n_q5,
               avg(CASE WHEN z.visitors_in > z.residents_away THEN 1.0 ELSE 0.0 END) AS q5_share_net_receiver,
               sum(CASE WHEN z.visitors_in > z.residents_away THEN 1 ELSE 0 END) AS q5_n_net_receiver,
               avg(z.visitors_in / q5.population) AS q5_mean_arrivals_per_res,
               avg(z.residents_away / q5.population) AS q5_mean_departures_per_res,
               avg((z.visitors_in - z.residents_away) / q5.population) AS q5_mean_net_per_res,
               median(z.visitors_in / q5.population) AS q5_med_arrivals_per_res,
               median(z.residents_away / q5.population) AS q5_med_departures_per_res
        FROM read_parquet('{pz}') z
        JOIN q5 ON z.zone_id = q5.id
        GROUP BY 1
        ORDER BY 1
    """).df()
    q5_daily = q5_daily.merge(wide[["date", "evening_date", "year", "evening_weekday", "split", "night_class", "desc_bucket"]], on="date", how="left")
    q5_daily.to_csv(OUT_DIR / "descriptive_q5_villages_daily.csv", index=False)

    log("Top 10 highest pooled beta nights (descriptive):")
    print(top_p[["date", "evening_weekday", "night_class", "desc_bucket", "beta_pooled_all", "beta_pooled_inland", "beta_pre1981_all", "flow_112_all"]].head(10).to_string(index=False))
    log("\nBucket summary (descriptive, sorted by pooled mean):")
    print(df_bkt[["desc_bucket", "n_nights", "beta_pooled_all_mean", "beta_pooled_inland_mean", "beta_pre1981_all_mean", "beta_vil2k_all_mean"]].to_string(index=False))
    log(f"Descriptive analysis finished in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()

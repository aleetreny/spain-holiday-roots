"""a03_manifest_and_panel: Build the 2022-2024 night manifest with pre-registered quality rules
and BOE calendar classes, and construct the three valid-night panels:
  (a) panel_cities137_prov.parquet  (137 census cities x 46 other mainland provinces, inland/coastal/all)
  (b) panel_cities675_prov.parquet  (675 cohort cities x 46 other mainland provinces, inland/coastal/all)
  (c) panel_zones2435_daily.parquet (2,435 mainland zones: residents_home, residents_away, visitors_in)

Run:
  .venv/bin/python -I journal/analysis/A_calendar/a03_manifest_and_panel.py
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a00_paths import (  # noqa: E402
    B0,
    C,
    CW_EXT_DIR,
    OUT_DIR,
    PANEL_DIR,
    PARQUET_DIR,
    ensure_dirs,
    make_logger,
)
from a01_calendar import (  # noqa: E402
    DANA_DATES,
    MITMA_INCIDENT_DATES,
    classify_file_dates,
    load_target_holiday_dates,
)
from a02_build_dataset import discover_raw_files  # noqa: E402


def get_city_sets(nodes: pd.DataFrame) -> tuple[list[str], list[str]]:
    """Exact 137 census cities (R2) and 675 single-zone >10k municipalities (V1)."""
    # 137 mainland census cities from table 68545 (2024)
    bm = pd.read_csv(
        CW_EXT_DIR / "ine_birth_residence" / "mun_residence_x_prov_birth_censo_2021_2025.csv",
        dtype={"ine_code": str, "prov_birth": str},
    )
    bm24 = bm[bm["year"] == 2024]
    cities137 = sorted(c for c in bm24["ine_code"].unique() if c[:2] not in C.EXCLUDED_PROVINCES)
    assert len(cities137) == 137 and set(cities137) <= set(nodes["id"])

    # 675 identifying single-zone >10k municipalities from Censo 2024 10% microdata (exact V1 rule)
    con = duckdb.connect()
    pq = CW_EXT_DIR / "ine_censo_microdatos_2024" / "md_CensoAnual_muestra_2024.parquet"
    spec_muni = "('991','992','993','666','966','777','888','999')"
    m_tot = con.execute(f"""
        SELECT CPRO || CMUN AS mun, CPRO AS prov
        FROM read_parquet('{pq}')
        WHERE CMUN NOT IN {spec_muni}
          AND CPRO NOT IN ('07','35','38','51','52')
        GROUP BY 1, 2
    """).df()
    single_set = set(nodes.loc[nodes["n_municipalities"] == 1, "id"])
    m_single = m_tot[m_tot["mun"].isin(single_set)]
    pcounts = m_single.groupby("prov").size()
    id_provs = set(pcounts[pcounts >= 2].index)
    cities675 = sorted(m_single.loc[m_single["prov"].isin(id_provs), "mun"].unique())
    assert len(cities675) == 675 and len(id_provs) == 44
    return cities137, cities675


def build_manifest(nodes: pd.DataFrame, threads: int = 4) -> pd.DataFrame:
    log = make_logger()
    ensure_dirs()

    raw_df = discover_raw_files()
    split_map = dict(zip(raw_df["file_date"], raw_df["split"]))

    # Find all partitioned parquet files currently built
    part_files = sorted(PARQUET_DIR.glob("date=*/data.parquet"))
    assert part_files, f"No parquet partitions found under {PARQUET_DIR}"
    part_dates = [p.parent.name.split("=", 1)[1] for p in part_files]
    log(f"Scanning {len(part_files)} partitioned dates for manifest quality metrics...")

    con = duckdb.connect()
    con.execute(f"SET threads TO {threads}")
    con.register("ml_nodes", nodes[["id", "prov", "coastal_true"]])

    # Compute per-date quality metrics in a single DuckDB scan
    file_list_sql = "[" + ",".join(f"'{p}'" for p in part_files) + "]"
    q_df = con.execute(f"""
        WITH x AS (
            SELECT CAST(date AS VARCHAR) AS file_date, r, s, p
            FROM read_parquet({file_list_sql})
        ),
        ml AS (
            SELECT x.file_date, x.r, x.s, x.p,
                   (nr.id IS NOT NULL) AS r_ml,
                   (ns.id IS NOT NULL) AS s_ml
            FROM x
            LEFT JOIN ml_nodes nr ON x.r = nr.id
            LEFT JOIN ml_nodes ns ON x.s = ns.id
        )
        SELECT file_date,
               count(*) AS n_rows_all,
               sum(p) AS total_p_all,
               sum(p) FILTER (WHERE r_ml AND s_ml) AS total_p_mainland,
               sum(p) FILTER (WHERE r_ml AND s_ml AND r != s) AS total_p_mainland_away,
               sum(p) FILTER (WHERE r_ml) AS total_p_res_mainland,
               count(DISTINCT r) FILTER (WHERE r_ml AND p > 0) AS n_res_zones_ml,
               count(DISTINCT s) FILTER (WHERE s_ml AND p > 0) AS n_sleep_zones_ml,
               count(DISTINCT r) FILTER (WHERE r_ml AND s_ml AND p > 0) AS n_res_zones_ml_in_ml
        FROM ml
        GROUP BY 1
        ORDER BY 1
    """).df()

    # Full 2022-2024 calendar (1,096 days) so missing MITMA dates are also explicit in manifest
    all_dates = [d.strftime("%Y-%m-%d") for d in pd.date_range("2022-01-01", "2024-12-31")]
    hols = load_target_holiday_dates()
    cal_df = classify_file_dates(all_dates, hols)

    man = cal_df.merge(q_df, on="file_date", how="left")
    man["file_present"] = man["file_date"].isin(set(part_dates))
    man["split"] = man["file_date"].map(split_map).fillna("missing")

    # Quality rules (PREREGISTRATION.md §2)
    man["q_fail_missing_file"] = ~man["file_present"]
    man["q_fail_ministry_incident"] = man["file_date"].isin(MITMA_INCIDENT_DATES)
    man["q_fail_dana"] = man["file_date"].isin(DANA_DATES)

    # Rule (c): fewer than 99% of median number of mainland residence zones report residents
    med_res_zones = float(man.loc[man["file_present"], "n_res_zones_ml"].median())
    man["median_res_zones_ml"] = med_res_zones
    man["res_zone_coverage"] = man["n_res_zones_ml"] / med_res_zones
    man["q_fail_low_coverage"] = man["file_present"] & (man["res_zone_coverage"] < 0.99)

    # Rule (d): total mainland persons (home + away) deviate by >10% from 28-night centred rolling median
    # Computed on the chronological calendar series (missing days ignored by rolling median)
    s_tot = man["total_p_mainland"].astype(float)
    # Use 28-night centred rolling median over available nights (also compute 29-night check)
    avail_idx = man.index[man["file_present"]]
    roll28 = pd.Series(np.nan, index=man.index, dtype=float)
    roll28.loc[avail_idx] = s_tot.loc[avail_idx].rolling(28, center=True, min_periods=14).median()
    man["rolling_med28_total_p_mainland"] = roll28
    man["rel_dev_total_p_mainland"] = (man["total_p_mainland"] - man["rolling_med28_total_p_mainland"]).abs() / man["rolling_med28_total_p_mainland"]
    man["q_fail_total_deviation"] = man["file_present"] & (man["rel_dev_total_p_mainland"] > 0.10)

    # Diagnostic check: 29-night centred rolling median vs 28-night
    roll29 = pd.Series(np.nan, index=man.index, dtype=float)
    roll29.loc[avail_idx] = s_tot.loc[avail_idx].rolling(29, center=True, min_periods=14).median()
    diff_28_29 = (
        (man["file_present"] & (((man["total_p_mainland"] - roll29).abs() / roll29) > 0.10))
        != man["q_fail_total_deviation"]
    ).sum()
    log(f"Rolling median check: 28-night vs 29-night disagreement count = {int(diff_28_29)}")

    man["is_valid"] = (
        man["file_present"]
        & (~man["q_fail_ministry_incident"])
        & (~man["q_fail_dana"])
        & (~man["q_fail_low_coverage"])
        & (~man["q_fail_total_deviation"])
    )

    def excl_reason(row) -> str:
        reasons = []
        if row.q_fail_missing_file:
            reasons.append("missing_file")
        if row.q_fail_ministry_incident:
            reasons.append("rule_a_ministry_incident")
        if row.q_fail_dana:
            reasons.append("rule_b_dana_window")
        if row.q_fail_low_coverage:
            reasons.append("rule_c_low_res_coverage")
        if row.q_fail_total_deviation:
            reasons.append("rule_d_total_dev_gt10pct")
        return ";".join(reasons)

    man["exclusion_reasons"] = [excl_reason(r) for r in man.itertuples(index=False)]
    man["night_class"] = np.where(man["is_valid"], man["night_class_calendar"], "excluded")

    out_csv = OUT_DIR / "night_manifest.csv"
    man.to_csv(out_csv, index=False)
    log(f"Saved {out_csv} ({len(man)} calendar days; {int(man['file_present'].sum())} files present; "
        f"{int(man['is_valid'].sum())} valid; {int((man['file_present'] & ~man['is_valid']).sum())} present-but-excluded)")

    excl_pres = man[man["file_present"] & ~man["is_valid"]]
    if len(excl_pres):
        log("Present files excluded by quality rules:")
        for r in excl_pres.itertuples(index=False):
            log(f"  {r.file_date} ({r.split}, cal={r.night_class_calendar}): {r.exclusion_reasons} "
                f"[res_cov={r.res_zone_coverage:.4f}, rel_dev={r.rel_dev_total_p_mainland:.4f}]")

    log("Valid night counts by split x night_class:")
    print(man[man["is_valid"]].groupby(["split", "night_class"]).size().unstack(fill_value=0).to_string())
    log("Valid night counts by year x split x night_class:")
    print(man[man["is_valid"]].groupby(["year", "split", "night_class"]).size().unstack(fill_value=0).to_string())
    return man


def build_panels(nodes: pd.DataFrame, man: pd.DataFrame, threads: int = 4) -> None:
    log = make_logger()
    ensure_dirs()
    cities137, cities675 = get_city_sets(nodes)
    valid_dates = sorted(man.loc[man["is_valid"], "file_date"].tolist())
    valid_files = [str(PARQUET_DIR / f"date={d}" / "data.parquet") for d in valid_dates]
    log(f"Building panels over {len(valid_dates)} valid nights...")

    con = duckdb.connect()
    con.execute(f"SET threads TO {threads}")
    con.register("ml_nodes", nodes[["id", "prov", "coastal_true"]])
    con.register("c137", pd.DataFrame({"city": cities137}))
    con.register("c675", pd.DataFrame({"city": cities675}))
    con.register("art_zones", pd.DataFrame({"id": list(C.ARTEFACT_ZONES)}))
    con.register("valid_dates_tbl", pd.DataFrame({"date": valid_dates}))

    file_list_sql = "[" + ",".join(f"'{f}'" for f in valid_files) + "]"

    # (a) 137 census cities -> other mainland provinces (split inland/coastal/all)
    # Note: matches b00.build_edges / r2_roots.py exactly (all 2,435 mainland dest zones, r != s, oprov != dprov).
    # Also includes flow_inland_noart (excluding the 3 artefact destination zones) for sensitivity transparency.
    t1 = time.time()
    p137 = PANEL_DIR / "panel_cities137_prov.parquet"
    con.execute(f"""
        COPY (
            SELECT CAST(x.date AS VARCHAR) AS date,
                   x.r AS city,
                   ns.prov AS dprov,
                   sum(CASE WHEN NOT ns.coastal_true THEN x.p ELSE 0.0 END) AS flow_inland,
                   sum(CASE WHEN ns.coastal_true THEN x.p ELSE 0.0 END) AS flow_coastal,
                   sum(x.p) AS flow_all,
                   sum(CASE WHEN NOT ns.coastal_true AND az.id IS NULL THEN x.p ELSE 0.0 END) AS flow_inland_noart
            FROM read_parquet({file_list_sql}) x
            JOIN c137 ON x.r = c137.city
            JOIN ml_nodes ns ON x.s = ns.id
            LEFT JOIN art_zones az ON x.s = az.id
            WHERE x.r != x.s
              AND substr(x.r, 1, 2) != ns.prov
            GROUP BY 1, 2, 3
            ORDER BY 1, 2, 3
        ) TO '{p137}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    n137 = con.execute(f"SELECT count(*), count(DISTINCT date) FROM read_parquet('{p137}')").fetchone()
    log(f"Saved {p137.name}: {n137[0]:,} rows across {n137[1]} nights ({time.time() - t1:.1f}s)")

    # (b) 675 cohort cities -> other mainland provinces (split inland/coastal/all)
    t2 = time.time()
    p675 = PANEL_DIR / "panel_cities675_prov.parquet"
    con.execute(f"""
        COPY (
            SELECT CAST(x.date AS VARCHAR) AS date,
                   x.r AS city,
                   ns.prov AS dprov,
                   sum(CASE WHEN NOT ns.coastal_true THEN x.p ELSE 0.0 END) AS flow_inland,
                   sum(CASE WHEN ns.coastal_true THEN x.p ELSE 0.0 END) AS flow_coastal,
                   sum(x.p) AS flow_all,
                   sum(CASE WHEN NOT ns.coastal_true AND az.id IS NULL THEN x.p ELSE 0.0 END) AS flow_inland_noart
            FROM read_parquet({file_list_sql}) x
            JOIN c675 ON x.r = c675.city
            JOIN ml_nodes ns ON x.s = ns.id
            LEFT JOIN art_zones az ON x.s = az.id
            WHERE x.r != x.s
              AND substr(x.r, 1, 2) != ns.prov
            GROUP BY 1, 2, 3
            ORDER BY 1, 2, 3
        ) TO '{p675}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    n675 = con.execute(f"SELECT count(*), count(DISTINCT date) FROM read_parquet('{p675}')").fetchone()
    log(f"Saved {p675.name}: {n675[0]:,} rows across {n675[1]} nights ({time.time() - t2:.1f}s)")

    # (c) Per-zone daily totals for all 2,435 mainland zones (for Analysis D)
    # Exact R1/b00 definitions: mainland-to-mainland network (r in ML and s in ML):
    #   residents_home = p where r == z and s == z
    #   residents_away = sum(p) where r == z and s in ML and s != z  (s_out in R1)
    #   visitors_in    = sum(p) where s == z and r in ML and r != z  (s_in in R1)
    #   residents_away_all_s = sum(p) where r == z and s != z (incl. islands/abroad)
    t3 = time.time()
    pz = PANEL_DIR / "panel_zones2435_daily.parquet"
    con.execute(f"""
        COPY (
            WITH x AS (
                SELECT CAST(date AS VARCHAR) AS date, r, s, p
                FROM read_parquet({file_list_sql})
            ),
            r_agg AS (
                SELECT x.date, x.r AS zone_id,
                       sum(CASE WHEN x.r = x.s THEN x.p ELSE 0.0 END) AS residents_home,
                       sum(CASE WHEN x.r != x.s AND ns.id IS NOT NULL THEN x.p ELSE 0.0 END) AS residents_away,
                       sum(CASE WHEN x.r != x.s THEN x.p ELSE 0.0 END) AS residents_away_all_s
                FROM x
                JOIN ml_nodes nr ON x.r = nr.id
                LEFT JOIN ml_nodes ns ON x.s = ns.id
                GROUP BY 1, 2
            ),
            s_agg AS (
                SELECT x.date, x.s AS zone_id,
                       sum(x.p) AS visitors_in
                FROM x
                JOIN ml_nodes nr ON x.r = nr.id
                JOIN ml_nodes ns ON x.s = ns.id
                WHERE x.r != x.s
                GROUP BY 1, 2
            )
            SELECT d.date,
                   n.id AS zone_id,
                   coalesce(r.residents_home, 0.0) AS residents_home,
                   coalesce(r.residents_away, 0.0) AS residents_away,
                   coalesce(s.visitors_in, 0.0) AS visitors_in,
                   coalesce(r.residents_home, 0.0) + coalesce(r.residents_away, 0.0) AS residents_mainland_total,
                   coalesce(r.residents_away_all_s, 0.0) AS residents_away_all_s,
                   coalesce(s.visitors_in, 0.0) - coalesce(r.residents_away, 0.0) AS net_flow,
                   (az.id IS NOT NULL) AS is_artefact_zone
            FROM valid_dates_tbl d
            CROSS JOIN ml_nodes n
            LEFT JOIN r_agg r ON d.date = r.date AND n.id = r.zone_id
            LEFT JOIN s_agg s ON d.date = s.date AND n.id = s.zone_id
            LEFT JOIN art_zones az ON n.id = az.id
            ORDER BY 1, 2
        ) TO '{pz}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    nz = con.execute(f"SELECT count(*), count(DISTINCT date), count(DISTINCT zone_id) FROM read_parquet('{pz}')").fetchone()
    assert nz[0] == len(valid_dates) * 2435 and nz[2] == 2435
    log(f"Saved {pz.name}: {nz[0]:,} rows ({nz[1]} nights x {nz[2]} zones) ({time.time() - t3:.1f}s)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--manifest-only", action="store_true")
    args = ap.parse_args()

    t0 = time.time()
    nodes = B0.load_nodes()
    man = build_manifest(nodes, threads=args.threads)
    if not args.manifest_only:
        build_panels(nodes, man, threads=args.threads)
    print(f"Total elapsed: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()

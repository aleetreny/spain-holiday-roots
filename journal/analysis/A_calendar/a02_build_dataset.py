"""a02_build_dataset: Convert all available 2022-2024 MITMA raw daily CSV.gz files
(discovery + confirmation) into one date-partitioned Parquet dataset (columns: date, r, s, p)
using the exact processing behind the coursework parquets (explore/p60, p61, p88),
and validate every discovery date against the coursework parquets.

Idempotent: skips dates whose partition parquet already exists unless --force is passed.
Run:
  /Users/alejandrotreny/Proyectos/spain-holiday-roots/.venv/bin/python -I journal/analysis/A_calendar/a02_build_dataset.py
"""
from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a00_paths import (
    CONFIRMATION_RAW_DIR,
    CW_PERNOCT_2024,
    CW_PERNOCT_2223,
    CW_REL_CSV,
    DISCOVERY_RAW_DIR,
    OUT_DIR,
    PARQUET_DIR,
    ensure_dirs,
    make_logger,
)

FILE_RE = re.compile(r"^(\d{8})_Pernoctaciones_municipios\.csv\.gz$")


def discover_raw_files() -> pd.DataFrame:
    """List all raw 2022-2024 MITMA daily files across discovery and confirmation dirs."""
    rows = []
    seen: dict[str, str] = {}
    for split, folder in (("discovery", DISCOVERY_RAW_DIR), ("confirmation", CONFIRMATION_RAW_DIR)):
        if not folder.exists():
            continue
        for p in sorted(folder.glob("*_Pernoctaciones_municipios.csv.gz")):
            m = FILE_RE.match(p.name)
            if not m:
                continue
            ymd = m.group(1)
            if ymd[:4] not in ("2022", "2023", "2024"):
                continue
            d_iso = f"{ymd[:4]}-{ymd[4:6]}-{ymd[6:8]}"
            if d_iso in seen:
                raise RuntimeError(f"Duplicate raw file date {d_iso} in both {seen[d_iso]} and {split}")
            seen[d_iso] = split
            rows.append({
                "file_date": d_iso,
                "split": split,
                "raw_path": str(p),
                "raw_bytes": p.stat().st_size,
            })
    return pd.DataFrame(rows).sort_values("file_date").reset_index(drop=True)


def convert_all(raw_df: pd.DataFrame, force: bool = False, threads: int = 4) -> dict:
    """Convert raw CSV.gz files to date-partitioned Parquet under PARQUET_DIR.

    Layout: PARQUET_DIR / f"date={d_iso}" / "data.parquet", with columns (date, r, s, p)
    written explicitly inside each file AND matching the Hive partition folder name.
    Using exact p60/p61/p88 logic:
      - read_csv(delim='|', header=true, columns={'fecha':'VARCHAR','zona_residencia':'VARCHAR',
                                                  'zona_pernoctacion':'VARCHAR','personas':'DOUBLE'})
      - inner join DISTINCT (distrito_mitma -> municipio_mitma) from relacion_ine_zonificacionMitma.csv
      - sum(p) grouped by (date, r, s)
    """
    log = make_logger()
    ensure_dirs()
    con = duckdb.connect()
    con.execute(f"SET threads TO {threads}")
    con.execute(f"""
        CREATE TABLE rel AS
        SELECT DISTINCT distrito_mitma AS dz, municipio_mitma AS mz
        FROM read_csv('{CW_REL_CSV}', delim='|', header=true, all_varchar=true)
    """)
    n_rel, n_dz = con.execute("SELECT count(*), count(DISTINCT dz) FROM rel").fetchone()
    assert n_rel == n_dz, f"Non-unique distrito_mitma mapping: {n_rel} vs {n_dz}"
    log(f"Loaded MITMA crosswalk: {n_rel} unique distrito_mitma -> municipio_mitma pairs")

    todo = []
    skipped = 0
    for row in raw_df.itertuples(index=False):
        part_dir = PARQUET_DIR / f"date={row.file_date}"
        out_file = part_dir / "data.parquet"
        if out_file.exists() and out_file.stat().st_size > 0 and not force:
            skipped += 1
        else:
            todo.append((row.file_date, row.raw_path, part_dir, out_file))

    log(f"Raw files found: {len(raw_df)} (already converted: {skipped}, to convert: {len(todo)})")

    unmapped_records = []
    t_conv = time.time()
    # Process in small batches or per-file so each date lands atomically
    for idx, (d_iso, raw_path, part_dir, out_file) in enumerate(todo, 1):
        part_dir.mkdir(parents=True, exist_ok=True)
        tmp_file = part_dir / "data.parquet.tmp"
        # Check for unmapped residence codes & verify single file date
        stats_row = con.execute(f"""
            WITH raw AS (
                SELECT strptime(CAST(fecha AS VARCHAR), '%Y%m%d')::DATE AS date,
                       CAST(zona_residencia AS VARCHAR) AS r_dz,
                       CAST(zona_pernoctacion AS VARCHAR) AS s,
                       CAST(personas AS DOUBLE) AS p
                FROM read_csv('{raw_path}', delim='|', header=true,
                              columns={{'fecha':'VARCHAR','zona_residencia':'VARCHAR',
                                        'zona_pernoctacion':'VARCHAR','personas':'DOUBLE'}})
            )
            SELECT count(*) AS raw_rows,
                   count(DISTINCT date) AS n_dates,
                   CAST(min(date) AS VARCHAR) AS min_d,
                   count(*) FILTER (WHERE rel.mz IS NULL) AS unmapped_rows,
                   coalesce(sum(raw.p) FILTER (WHERE rel.mz IS NULL), 0.0) AS unmapped_p,
                   coalesce(sum(raw.p), 0.0) AS raw_total_p
            FROM raw LEFT JOIN rel ON raw.r_dz = rel.dz
        """).fetchone()
        raw_rows, n_dates, min_d, unmapped_rows, unmapped_p, raw_total_p = stats_row
        assert n_dates == 1 and min_d == d_iso, f"Date mismatch in {raw_path}: {n_dates}, {min_d} vs {d_iso}"
        if unmapped_rows > 0:
            unmapped_records.append((d_iso, int(unmapped_rows), float(unmapped_p), float(raw_total_p)))

        con.execute(f"""
            COPY (
                SELECT strptime(CAST(x.fecha AS VARCHAR), '%Y%m%d')::DATE AS date,
                       rel.mz AS r,
                       CAST(x.zona_pernoctacion AS VARCHAR) AS s,
                       sum(CAST(x.personas AS DOUBLE)) AS p
                FROM read_csv('{raw_path}', delim='|', header=true,
                              columns={{'fecha':'VARCHAR','zona_residencia':'VARCHAR',
                                        'zona_pernoctacion':'VARCHAR','personas':'DOUBLE'}}) x
                JOIN rel ON CAST(x.zona_residencia AS VARCHAR) = rel.dz
                GROUP BY 1, 2, 3
                ORDER BY 2, 3
            ) TO '{tmp_file}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """)
        tmp_file.replace(out_file)
        if idx % 50 == 0 or idx == len(todo):
            log(f"  converted {idx}/{len(todo)} dates ({time.time() - t_conv:.1f}s)")

    if unmapped_records:
        log(f"Note: {len(unmapped_records)} newly converted dates had unmapped residence codes "
            f"(e.g. foreign/special codes dropped by inner join, matching p61/p88). "
            f"Max unmapped_p = {max(x[2] for x in unmapped_records):.1f}")
    return {
        "n_raw": len(raw_df),
        "n_skipped": skipped,
        "n_converted": len(todo),
        "elapsed_s": time.time() - t_conv,
    }


def validate_against_coursework(threads: int = 4) -> pd.DataFrame:
    """Validate EVERY discovery date in PARQUET_DIR against coursework parquets."""
    log = make_logger()
    con = duckdb.connect()
    con.execute(f"SET threads TO {threads}")

    rows = []
    for label, cw_path in [("2022_2023", CW_PERNOCT_2223), ("2024_6m", CW_PERNOCT_2024)]:
        cw_dates = [
            str(r[0]) for r in con.execute(
                "SELECT DISTINCT CAST(date AS VARCHAR) FROM read_parquet(?) ORDER BY 1",
                [str(cw_path)],
            ).fetchall()
        ]
        log(f"Validating {len(cw_dates)} discovery dates against {cw_path.name}...")
        part_files = [str(PARQUET_DIR / f"date={d}" / "data.parquet") for d in cw_dates]
        missing_parts = [f for f in part_files if not Path(f).exists()]
        assert not missing_parts, f"Missing partitioned files for discovery dates: {missing_parts[:5]}"

        # Exact full outer join per (date, r, s)
        df_val = con.execute(f"""
            WITH cw AS (
                SELECT date, r, s, sum(p) AS p_cw
                FROM read_parquet('{cw_path}')
                GROUP BY 1, 2, 3
            ),
            nw AS (
                SELECT date, r, s, sum(p) AS p_new
                FROM read_parquet({part_files})
                GROUP BY 1, 2, 3
            ),
            j AS (
                SELECT coalesce(cw.date, nw.date) AS date,
                       coalesce(cw.r, nw.r) AS r,
                       coalesce(cw.s, nw.s) AS s,
                       cw.p_cw,
                       nw.p_new,
                       abs(coalesce(cw.p_cw, 0.0) - coalesce(nw.p_new, 0.0)) AS abs_diff
                FROM cw FULL OUTER JOIN nw
                  ON cw.date = nw.date AND cw.r = nw.r AND cw.s = nw.s
            )
            SELECT CAST(date AS VARCHAR) AS date,
                   count(*) AS n_pairs,
                   count(*) FILTER (WHERE p_cw IS NULL) AS extra_in_new,
                   count(*) FILTER (WHERE p_new IS NULL) AS missing_in_new,
                   sum(p_cw) AS sum_p_cw,
                   sum(p_new) AS sum_p_new,
                   max(abs_diff) AS max_abs_diff,
                   sum(abs_diff) AS sum_abs_diff
            FROM j
            GROUP BY 1
            ORDER BY 1
        """).df()
        df_val.insert(0, "source_parquet", label)
        rows.append(df_val)
        log(
            f"  {label}: {len(df_val)} dates checked | total pairs={df_val['n_pairs'].sum():,} | "
            f"missing_in_new={int(df_val['missing_in_new'].sum())} | extra_in_new={int(df_val['extra_in_new'].sum())} | "
            f"max_abs_diff={df_val['max_abs_diff'].max():.3e}"
        )

    out = pd.concat(rows, ignore_index=True)
    assert len(out) == 306, f"Expected 306 discovery dates, validated {len(out)}"
    assert int(out["missing_in_new"].sum()) == 0, "Some coursework (r,s) pairs missing in new dataset!"
    assert int(out["extra_in_new"].sum()) == 0, "Extra (r,s) pairs found in new dataset!"
    assert out["max_abs_diff"].max() < 1e-7, f"Max abs diff too large: {out['max_abs_diff'].max()}"

    out_csv = OUT_DIR / "dataset_validation.csv"
    out.to_csv(out_csv, index=False)
    log(f"VALIDATION PASSED on all 306 discovery dates! Max abs diff = {out['max_abs_diff'].max():.3e}. Saved {out_csv}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="Re-convert existing partition parquets")
    ap.add_argument("--skip-validate", action="store_true")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    t0 = time.time()
    raw_df = discover_raw_files()
    info = convert_all(raw_df, force=args.force, threads=args.threads)
    print(f"Conversion summary: {info}")
    if not args.skip_validate:
        validate_against_coursework(threads=args.threads)
    print(f"Total elapsed: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()

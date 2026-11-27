"""Build the 2025 and bridge inputs with Analysis A. Run with python -I -B."""
from __future__ import annotations

import datetime as dt
import hashlib
import sys
from pathlib import Path
from unittest.mock import patch

import duckdb
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
JOURNAL = Path("/Users/alejandrotreny/Proyectos/spain-holiday-roots-journal")
A_DIR = JOURNAL / "journal/analysis/A_calendar"
A_HEAVY = JOURNAL / "_local/journal_data/A_calendar"
HEAVY = JOURNAL / "_local/journal_data/A2_confirmation2025"
RAW = JOURNAL / "_local/external_journal/mitma_pernoct_2025"
STANDARD = HEAVY / "standard"
RECALC = HEAVY / "recalc"
A_SOURCE = HERE.parents[2] / "_local/A_reference"
sys.path.insert(0, str(A_SOURCE if A_SOURCE.is_dir() else A_DIR))

import a00_paths as paths
import a01_calendar as cal
import a02_build_dataset as dataset
import a03_manifest_and_panel as panel
import a04_nightly_series as nightly
from a05_solver import fit_single_block, make_block_design

AUG = {"2025-08-10", "2025-08-14", "2025-08-15", "2025-08-16", "2025-08-17"}
EASTER = "2025-04-18"
BRIDGE = ["2024-10-15", "2024-10-18", "2024-10-19", "2025-01-13", "2025-01-19", "2025-02-08"]
BLACKOUT_FILES = {"2025-04-28", "2025-04-29"}
BOE_ID = "BOE-A-2024-21316"
BOE_TXT = f"https://www.boe.es/diario_boe/txt.php?id={BOE_ID}"
BOE_XML = f"https://www.boe.es/diario_boe/xml.php?id={BOE_ID}"


def file_list() -> pd.DataFrame:
    files = pd.read_csv(RAW / "manifest.csv", dtype={"file_date": str})
    assert files.groupby("kind").size().to_dict() == {"recalc": 6, "standard": 196}
    assert files["status"].eq("OK").all()
    assert not files.duplicated(["kind", "file_date"]).any()
    files["file_date"] = pd.to_datetime(files["file_date"], format="%Y%m%d").dt.strftime("%Y-%m-%d")
    files["raw_path"] = [
        str(RAW / kind / f"{d.replace('-', '')}_Pernoctaciones_municipios.csv.gz")
        for kind, d in zip(files["kind"], files["file_date"])
    ]
    for row in files.itertuples(index=False):
        source = Path(row.raw_path)
        assert source.stat().st_size == row.bytes, source
        assert hashlib.sha256(source.read_bytes()).hexdigest() == row.sha256, source
    assert set(files.loc[files["kind"].eq("recalc"), "file_date"]) == set(BRIDGE)
    files["method"] = np.where(
        files["kind"].eq("recalc") | files["file_date"].ge("2025-07-01"), "new", "old"
    )
    files["split"] = "confirmation_ii"
    return files


def calendar_2025(_dates, holidays: set[dt.date]) -> pd.DataFrame:
    dates = pd.date_range("2025-01-01", "2025-12-31").strftime("%Y-%m-%d").tolist()
    with (
        patch.dict(cal.HOLY_WEEK_EVENING_RANGES, {2025: (dt.date(2025, 4, 11), dt.date(2025, 4, 21))}),
        patch.object(cal, "EASTER_FILE_DATES", {EASTER}),
    ):
        frame = cal.classify_file_dates(dates, holidays)
    ordinary_months = frame["file_date"].str[5:7].isin(["01", "02", "03", "04", "05", "06", "10"])
    classes = frame["night_class_calendar"]
    frame.loc[classes.eq("O") & ~ordinary_months, "night_class_calendar"] = "other"
    frame.loc[classes.eq("AUG") & ~frame["file_date"].isin(AUG), "night_class_calendar"] = "other"
    frame.loc[classes.isin(["XMAS", "ALLSAINTS"]), "night_class_calendar"] = "other"
    return frame


def mainland_without_artefacts(nodes: pd.DataFrame) -> pd.DataFrame:
    artefacts = nodes["id"].isin(paths.C.ARTEFACT_ZONES)
    assert artefacts.sum() == len(paths.C.ARTEFACT_ZONES)
    assert not nodes.loc[artefacts, "coastal_true"].any()
    return nodes.loc[~artefacts]


def build_standard_manifest(nodes: pd.DataFrame, files: pd.DataFrame) -> pd.DataFrame:
    xml_file = HEAVY / "boe_2025.xml"
    boe_source = {
        "boe_id": BOE_ID,
        "url": BOE_TXT,
        "xml_url": BOE_XML,
        "eli_url": "https://www.boe.es/eli/es/res/2024/10/15/(1)",
        "accessed": "2026-11-27",
    }
    with patch.dict(cal.BOE_SOURCES, {2025: boe_source}):
        holidays = cal.parse_boe_year(2025, xml_file.read_bytes() if xml_file.exists() else cal.fetch_boe_xml(2025))
    assert "2025-04-18" in set(holidays.loc[holidays["is_national_all_spain"], "date"])
    holidays_path = OUT / "boe_holidays_2025.csv"
    holidays.to_csv(holidays_path, index=False)
    with (
        patch.object(panel, "discover_raw_files", lambda: files),
        patch.object(panel, "load_target_holiday_dates", lambda: cal.load_target_holiday_dates(holidays_path)),
        patch.object(panel, "classify_file_dates", calendar_2025),
        patch.object(panel, "MITMA_INCIDENT_DATES", BLACKOUT_FILES),
        patch.object(panel, "PARQUET_DIR", STANDARD / "parquet"),
        patch.object(panel, "OUT_DIR", STANDARD),
        patch.object(panel, "ensure_dirs", lambda: None),
    ):
        manifest = panel.build_manifest(mainland_without_artefacts(nodes))
    manifest = manifest.merge(
        files[["file_date", "method", "url", "bytes", "sha256"]], on="file_date", how="left", validate="1:1"
    )
    assert set(manifest.loc[manifest["file_present"], "file_date"]) == set(files["file_date"])
    class_names = {"AUG": "AUG2025", "O": "O2025", "EASTER": "EASTER2025"}
    for col in ("night_class", "night_class_calendar"):
        manifest[col] = manifest[col].replace(class_names)
    manifest.to_csv(OUT / "night_manifest.csv", index=False)
    return manifest


def build_panels(nodes: pd.DataFrame, manifest: pd.DataFrame, directory: Path) -> None:
    target = directory / "panel"
    target.mkdir(parents=True, exist_ok=True)
    with (
        patch.object(panel, "PARQUET_DIR", directory / "parquet"),
        patch.object(panel, "PANEL_DIR", target),
        patch.object(panel, "ensure_dirs", lambda: None),
    ):
        panel.build_panels(nodes, manifest)


def bridge_flow_stats(nodes: pd.DataFrame, dates: list[str]) -> pd.DataFrame:
    """Compare matched mainland zone pairs, treating an unpublished pair as zero."""
    rows = []
    with duckdb.connect() as con:
        con.register("mainland", nodes[["id"]])
        for date in dates:
            prior = A_HEAVY / "pernoct_mun_2022_2024" if date.startswith("2024") else STANDARD / "parquet"
            old = prior / f"date={date}" / "data.parquet"
            new = RECALC / "parquet" / f"date={date}" / "data.parquet"
            counts = con.execute("""
                WITH matched AS (
                    SELECT coalesce(o.r, n.r) AS r, coalesce(o.s, n.s) AS s,
                           coalesce(o.p, 0) AS old_p, coalesce(n.p, 0) AS new_p
                    FROM read_parquet(?) o FULL JOIN read_parquet(?) n USING (r, s)
                ),
                ml AS (
                    SELECT *, r != s AS away, abs(old_p - new_p) > 1e-8 AS changed,
                           old_p > 0 OR new_p > 0 AS in_union
                    FROM matched JOIN mainland a ON r = a.id JOIN mainland b ON s = b.id
                )
                SELECT sum(old_p) AS old_total, sum(new_p) AS new_total,
                       sum(old_p) FILTER (WHERE away) AS old_away,
                       sum(new_p) FILTER (WHERE away) AS new_away,
                       count(*) FILTER (WHERE away AND in_union) AS n_flows,
                       count(*) FILTER (WHERE away AND in_union AND changed) AS n_flows_changed,
                       count(*) FILTER (WHERE away AND old_p = 0 AND new_p > 0) AS n_flows_added,
                       count(*) FILTER (WHERE away AND old_p > 0 AND new_p = 0) AS n_flows_removed
                FROM ml
            """, [str(old), str(new)]).df().iloc[0].to_dict()
            rows.append({"file_date": date, **counts})
    out = pd.DataFrame(rows)
    count_cols = [c for c in out if c.startswith("n_flows")]
    out[count_cols] = out[count_cols].astype("int64")
    out["share_flows_changed"] = out["n_flows_changed"] / out["n_flows"]
    out["relative_change_total"] = out["new_total"] / out["old_total"] - 1
    out["relative_change_away"] = out["new_away"] / out["old_away"] - 1
    return out


def point_estimates(nodes: pd.DataFrame, directory: Path, dates: list[str], use_noart: bool) -> pd.DataFrame:
    out = pd.DataFrame(index=pd.Index(dates, name="date"))
    for term, city_set, design, x_key in (
        ("pooled", 137, nightly.build_pooled_meta, "Xr"),
        ("pre1981", 675, nightly.build_cohort_meta, "X_COH"),
    ):
        meta = design(nodes)
        flows = nightly.load_panel_flows(
            directory / f"panel_cities{city_set}_prov.parquet", meta, dates, use_noart=use_noart,
        )
        for d, flow in flows.items():
            block = make_block_design(
                d, int(d[:4]), flow["all"], meta[x_key], meta["cov_dist"],
                meta["c_code"], meta["p_code"], meta["cl_code"], slope_raw=meta["slope_raw"],
            )
            out.loc[d, term] = fit_single_block(block, meta["G_tot"])["beta"][0]
    assert set(out.index) == set(dates) and out.notna().all().all()
    return out


def bridge(nodes: pd.DataFrame, out_dir: Path = OUT, use_noart: bool = True) -> None:
    old = pd.concat([
        point_estimates(nodes, STANDARD / "panel", [d for d in BRIDGE if d.startswith("2025")], use_noart),
        point_estimates(nodes, A_HEAVY / "panel", [d for d in BRIDGE if d.startswith("2024")], use_noart),
    ])
    new = point_estimates(nodes, RECALC / "panel", BRIDGE, use_noart)
    compare = old.add_prefix("old_").join(new.add_prefix("new_")).rename_axis("file_date")
    summary = []
    for term, tolerance in (("pooled", 0.035), ("pre1981", 0.024)):
        compare[f"difference_{term}"] = compare[f"new_{term}"] - compare[f"old_{term}"]
        compare[f"abs_difference_{term}"] = compare[f"difference_{term}"].abs()
        mean_absolute_difference = compare[f"abs_difference_{term}"].mean()
        summary.append({
            "term": term, "n_paired_dates": len(compare),
            "mean_absolute_difference": mean_absolute_difference,
            "tolerance": tolerance, "pass": mean_absolute_difference < tolerance,
        })
    flow_nodes = mainland_without_artefacts(nodes) if use_noart else nodes
    compare.reset_index().merge(
        bridge_flow_stats(flow_nodes, BRIDGE), on="file_date", validate="1:1"
    ).sort_values("file_date").to_csv(out_dir / "method_bridge_dates.csv", index=False)
    frame = pd.DataFrame(summary)
    frame["joint_pass"] = frame["pass"].all()
    frame.to_csv(out_dir / "method_bridge_summary.csv", index=False)
    print(frame.to_string(index=False), flush=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    nodes = paths.B0.load_nodes()
    with patch.object(nightly, "PANEL_DIR", A_HEAVY / "panel"), patch.object(nightly, "OUT_DIR", OUT):
        repro = nightly.run_reproduction_check(nodes)
    a_repro = pd.read_csv(A_DIR / "out/reproduction_check.csv")
    assert repro["check"].tolist() == a_repro["check"].tolist()
    np.testing.assert_allclose(repro["reproduced"], a_repro["reproduced"], rtol=0, atol=1e-10)
    files = file_list()
    for kind, directory in (("standard", STANDARD), ("recalc", RECALC)):
        directory.mkdir(parents=True, exist_ok=True)
        sample = files.loc[files["kind"].eq(kind)]
        with (
            patch.object(dataset, "PARQUET_DIR", directory / "parquet"),
            patch.object(dataset, "ensure_dirs", lambda: None),
        ):
            dataset.convert_all(sample)
        manifest = (
            build_standard_manifest(nodes, sample) if kind == "standard"
            else sample.assign(is_valid=True)
        )
        build_panels(nodes, manifest, directory)
    bridge(nodes)


if __name__ == "__main__":
    main()

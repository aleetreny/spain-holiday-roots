"""Rerun Analysis A with the §1 zone exclusion, using A as a read-only dependency."""

import argparse
import os
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

OUT = Path(__file__).resolve().parent / "out"
DEFAULT_A = Path("/Users/alejandrotreny/Proyectos/spain-holiday-roots-journal/journal/analysis/A_calendar")
ARTEFACT_ZONES = {"31905_AM", "19103_AM", "40149_AM"}
PANEL_NAMES = ("panel_cities137_prov.parquet", "panel_cities675_prov.parquet")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-a", type=Path, default=DEFAULT_A)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")

    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.dont_write_bytecode = True
    a_dir = args.analysis_a.resolve()
    sys.path.insert(0, str(a_dir))

    import duckdb
    import pandas as pd
    import a00_paths as paths
    import a03_manifest_and_panel as panel
    import a04_nightly_series as nightly
    import a05_preregistered_tests as tests

    nodes = paths.B0.load_nodes()
    excluded = nodes[nodes["id"].isin(ARTEFACT_ZONES)]
    if set(excluded["id"]) != ARTEFACT_ZONES or excluded["coastal_true"].any():
        raise ValueError("Expected all three artefact zones in the mainland frame and inland")
    if set(paths.C.ARTEFACT_ZONES) != ARTEFACT_ZONES:
        raise ValueError("Analysis A and §1 artefact zone lists differ")

    # Both saved panels already carry flow_inland_noart. This is sufficient for
    # both arms because the three excluded destinations are inland.
    with duckdb.connect() as con:
        for name, cities in zip(PANEL_NAMES, panel.get_city_sets(nodes), strict=True):
            file = str(paths.PANEL_DIR / name)
            if ARTEFACT_ZONES.intersection(cities):
                raise ValueError(f"An artefact zone belongs to the original city frame for {name}")
            saved_cities = {row[0] for row in con.execute(
                "SELECT DISTINCT city FROM read_parquet(?)", [file]
            ).fetchall()}
            if saved_cities != set(cities):
                raise ValueError(f"{name} does not match the original city frame")
            print(f"{name}: {len(cities)} origins, zero artefact origins", flush=True)

    manifest = pd.read_csv(a_dir / "out" / "night_manifest.csv")
    needed = sorted(manifest.loc[
        manifest["is_valid"] & manifest["night_class"].isin(["O", "AUG", "XMAS", "EASTER", "ALLSAINTS"]),
        "file_date",
    ])
    pooled = nightly.build_pooled_meta(nodes)
    common = nightly.build_pooled_meta(nodes, only_dest_provs=set(pooled["COASTAL_PROVS"]))
    cohort = nightly.build_cohort_meta(nodes)
    print(f"Loading artefact-excluded flows for {len(needed)} classified nights", flush=True)
    flows = [
        nightly.load_panel_flows(paths.PANEL_DIR / name, meta, dates=needed, use_noart=True)
        for name, meta in zip((PANEL_NAMES[0], *PANEL_NAMES), (pooled, common, cohort), strict=True)
    ]

    OUT.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="analysis_a_noart_", dir=OUT) as staging:
        staging = Path(staging)
        with patch.object(tests, "OUT_DIR", staging):
            _, pooled_blocks, cohort_blocks, p3_blocks = tests.run_all_hypotheses(
                manifest, pooled, common, cohort, *flows, workers=args.workers
            )
            (staging / "preregistered_tests.csv").replace(OUT / "preregistered_tests_noart.csv")
            tests.run_placebo_windows(
                manifest, pooled_blocks, cohort_blocks, p3_blocks, pooled, cohort,
                workers=min(args.workers, 6),
            )
        for name in ("placebo_summary", "placebo_windows"):
            (staging / f"{name}.csv").replace(OUT / f"{name}_noart.csv")


if __name__ == "__main__":
    main()

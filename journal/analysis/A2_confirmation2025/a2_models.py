"""Fit Confirmation II tests, diagnostics and nightly series using Analysis A."""
from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a2_data import A_DIR, AUG, EASTER, OUT, RECALC, STANDARD, nightly, paths
import a05_preregistered_tests as tests

HYPOTHESES = ("P1", "P2", "P3", "S2", "S4", "S5", "S6")
AUG_HYPOTHESES = tuple(h for h in HYPOTHESES if h != "S2")
RECALC_REFERENCE = "2024-10-15"
PANEL_POOLED = "panel_cities137_prov.parquet"
PANEL_COHORT = "panel_cities675_prov.parquet"


def analysis_samples(manifest: pd.DataFrame, bridge_pass: bool) -> dict[str, list[str]]:
    valid = manifest.loc[manifest["is_valid"]]
    aug = sorted(valid.loc[valid["night_class"].eq("AUG2025"), "file_date"])
    assert set(aug) == AUG
    assert valid.loc[valid["file_date"].eq(EASTER), "night_class"].tolist() == ["EASTER2025"]
    ordinary = valid.loc[valid["night_class"].eq("O2025")]
    old = sorted(ordinary.loc[ordinary["method"].eq("old"), "file_date"])
    new = sorted(ordinary.loc[ordinary["method"].eq("new"), "file_date"])
    all_o = sorted(old + new)
    return {
        "aug": aug,
        "aug_reference": all_o if bridge_pass else [RECALC_REFERENCE, *new],
        "easter_reference": all_o if bridge_pass else old,
    }


def load_flows(filename: str, meta: dict, dates: list[str], use_noart: bool) -> dict:
    standard_dates = [d for d in dates if d != RECALC_REFERENCE]
    flows = nightly.load_panel_flows(STANDARD / "panel" / filename, meta, standard_dates, use_noart=use_noart)
    if RECALC_REFERENCE in dates:
        flows.update(nightly.load_panel_flows(RECALC / "panel" / filename, meta, [RECALC_REFERENCE], use_noart=use_noart))
    assert set(flows) == set(dates)
    return flows


def design_blocks(root_year: int, kind: str, nodes: pd.DataFrame, dates: list[str], use_noart: bool):
    if kind in ("cohort_joint", "birthplace_joint"):
        meta = nightly.build_cohort_meta(nodes)
        flows = load_flows(PANEL_COHORT, meta, dates, use_noart)
        blocks = tests._prep_blocks_cohort(dates, flows, {2024: meta, 2025: meta}, kind)
        return {b.label: replace(b, year=2025) for b in blocks}

    meta = nightly.build_pooled_meta(nodes, roots_year=root_year)
    if kind == "pooled":
        flows = load_flows(PANEL_POOLED, meta, dates, use_noart)
        blocks = tests._prep_blocks_pooled(dates, flows, {2024: meta, 2025: meta})
        return {b.label: replace(b, year=2025) for b in blocks}

    assert kind == "p3"
    common = nightly.build_pooled_meta(nodes, roots_year=root_year, only_dest_provs=set(meta["COASTAL_PROVS"]))
    flows = load_flows(PANEL_POOLED, common, dates, use_noart)
    blocks, _ = tests._prep_blocks_p3(dates, flows, common)
    return {
        d: (replace(bc, year=2025), replace(bi, year=2025))
        for d, bc, bi in zip(dates, blocks[::2], blocks[1::2], strict=True)
    }


def sample_manifest(holiday_dates: list[str], ordinary_dates: list[str], holiday: str) -> pd.DataFrame:
    """Assign the addendum's optional 2024 recalculation to the 2025 contrast."""
    return pd.DataFrame({
        "file_date": holiday_dates + ordinary_dates,
        "night_class": [holiday] * len(holiday_dates) + ["O"] * len(ordinary_dates),
        "is_valid": True,
        "year": 2025,
    })


def eval_test(job: tuple[str, int, pd.DataFrame, tuple, np.ndarray, np.ndarray]) -> dict:
    hyp, roots_year, sample, blocks, vb22, vb44 = job
    tests._init_hyp_worker((sample, *blocks, vb22, vb44))
    out = tests._eval_one_hyp((hyp, "all_nights", None, True))
    assert out is not None and np.isfinite([out["est"], out["ci_wcr_lo"], out["ci_wcr_hi"], out["p_wcr"]]).all()
    out["split"] = "confirmation_ii"
    out["roots_year"] = roots_year
    print(f"{roots_year} {hyp}: {out['est']:+.5f}, p={out['p_wcr']:.4f} ({out['elapsed_s']:.1f}s)", flush=True)
    return out


def run_placebos(samples: dict, blocks: dict) -> None:
    aug, ordinary = samples["aug"], samples["aug_reference"]
    if len(ordinary) < 18:
        pd.DataFrame([{
            "year": 2025, "n_ordinary_reference": len(ordinary), "n_placebo_windows": 0,
            "status": "insufficient_ordinary_reference_nights_for_analysis_A_placebos",
        }]).to_csv(OUT / "placebo_summary.csv", index=False)
        return
    _, _, _, _, rows, summary = tests._placebo_one_year((
        "confirmation_ii", 2025, aug, ordinary, blocks["pooled"],
        blocks["cohort_joint"], blocks["p3"], 22, 44,
    ))
    pd.DataFrame(rows).to_csv(OUT / "placebo_windows.csv", index=False)
    pd.DataFrame(summary).to_csv(OUT / "placebo_summary.csv", index=False)


def run_model_tests(manifest: pd.DataFrame, samples: dict, bridge_pass: bool, workers: int, mode: str, use_noart: bool) -> None:
    aug, aug_o, easter_o = (samples[k] for k in ("aug", "aug_reference", "easter_reference"))
    all_dates = sorted(set(aug + aug_o + [EASTER] + easter_o))
    aug_dates = sorted(set(aug + aug_o))
    aug_sample = sample_manifest(aug, aug_o, "AUG")
    easter_sample = sample_manifest([EASTER], easter_o, "EASTER")
    nodes = paths.B0.load_nodes()

    prep = {
        (2024, "pooled"): all_dates,
        (2024, "cohort_joint"): aug_dates,
        (2024, "birthplace_joint"): aug_dates,
        (2024, "p3"): aug_dates,
    }
    if mode == "tests":
        prep[(2025, "pooled")] = all_dates
        prep[(2025, "p3")] = aug_dates
    else:
        prep.pop((2024, "birthplace_joint"))

    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending = {
            pool.submit(design_blocks, year, kind, nodes, dates, use_noart): (year, kind)
            for (year, kind), dates in prep.items()
        }
        built = {}
        for f in as_completed(pending):
            key = pending[f]
            built[key] = f.result()
            print(f"Prepared {key[0]} {key[1]} blocks", flush=True)

        if mode == "placebos":
            run_placebos(samples, {kind: val for (_, kind), val in built.items()})
            return

        vb22 = tests.webb_draws("confirmation_ii_pooled_G22", 22)
        vb44 = tests.webb_draws("confirmation_ii_cohort_G44", 44)
        assert vb22.shape == (9999, 22) and vb44.shape == (9999, 44)
        # Submit only the design used by each test to limit worker memory.
        kinds = ("pooled", "cohort_joint", "birthplace_joint", "p3")
        kind_of = {
            "P1": "pooled", "S2": "pooled", "P3": "p3", "P2": "cohort_joint",
            "S4": "birthplace_joint", "S5": "cohort_joint", "S6": "cohort_joint",
        }
        pending = []
        for roots_year, hyps in ((2024, HYPOTHESES), (2025, ("P1", "P3", "S2"))):
            for hyp in hyps:
                test_sample = easter_sample if hyp == "S2" else aug_sample
                arguments = tuple(built[roots_year, k] if k == kind_of[hyp] else {} for k in kinds)
                pending.append(pool.submit(eval_test, (hyp, roots_year, test_sample, arguments, vb22, vb44)))
        results = pd.DataFrame([f.result() for f in as_completed(pending)])

    # The cohort models use the fixed 2024 microdata for both versions of the pooled exposure.
    fixed = results.loc[results["roots_year"].eq(2024) & results["hypothesis"].isin(["P2", "S4", "S5", "S6"])]
    sensitivity = pd.concat([results.loc[results["roots_year"].eq(2025)], fixed], ignore_index=True)
    for frame, filename in (
        (results.loc[results["roots_year"].eq(2024)].copy(), "preregistered_tests.csv"),
        (sensitivity, "roots_year_sensitivity.csv"),
    ):
        frame["p_holm"] = frame.groupby("family")["p_wcr"].transform(lambda p: tests.holm_adjust(p.to_numpy()))
        frame["alternative_sign"] = np.where(frame["hypothesis"].eq("S5"), -1, 1)
        frame["passes_holm_in_direction"] = frame["p_holm"].le(0.05) & frame["est"].mul(frame["alternative_sign"]).gt(0)
        frame["aug_reference_rule"] = "all_O2025" if bridge_pass else "new_method_only_low_power"
        frame.sort_values("hypothesis").to_csv(OUT / filename, index=False)
    manifest = manifest[["file_date", "method"]].set_index("file_date")
    sample_rows = []
    for label, sample in (("AUG", aug_sample), ("EASTER", easter_sample)):
        for row in sample.itertuples(index=False):
            sample_rows.append({
                "contrast": label, "file_date": row.file_date, "role": "reference" if row.night_class == "O" else "holiday",
                "method": "new" if row.file_date == RECALC_REFERENCE else manifest.loc[row.file_date, "method"],
                "contrast_year": row.year,
            })
    pd.DataFrame(sample_rows).to_csv(OUT / "contrast_nights.csv", index=False)


def nightly_series(nodes: pd.DataFrame, manifest: pd.DataFrame, workers: int, use_noart: bool) -> None:
    load = nightly.load_panel_flows
    with (
        patch.object(nightly, "OUT_DIR", OUT),
        patch.object(nightly, "PANEL_DIR", STANDARD / "panel"),
        patch.object(nightly, "load_panel_flows", lambda *a, **kw: load(*a, **kw, use_noart=use_noart)),
    ):
        nightly.run_nightly_series(nodes, workers=workers)
    lookup = manifest.set_index("file_date")["method"]
    for suffix in ("pooled", "cohorts"):
        path = OUT / f"nightly_{suffix}.csv"
        frame = pd.read_csv(path)
        frame["method"] = frame["date"].map(lookup)
        assert frame["method"].notna().all()
        frame.to_csv(path, index=False)


def combined_decision() -> None:
    a = pd.read_csv(A_DIR / "out/preregistered_tests.csv")
    a = a.query("split == 'confirmation_O' and year_scope == 'equal_year_avg'")
    b = pd.read_csv(OUT / "preregistered_tests.csv")
    cols = ["hypothesis", "est", "ci_wcr_lo", "ci_wcr_hi", "p_holm"]
    joint = a[cols].merge(b[cols + ["alternative_sign"]], on="hypothesis", suffixes=("_A", "_A2"), validate="1:1")
    joint["same_sign"] = np.sign(joint["est_A"]) == np.sign(joint["est_A2"])
    for label in ("A", "A2"):
        joint[f"passes_{label}"] = (
            joint[f"est_{label}"].mul(joint["alternative_sign"]).gt(0)
            & joint[f"p_holm_{label}"].le(0.05)
        )
    joint["abstract_eligible"] = (
        joint["hypothesis"].isin(AUG_HYPOTHESES)
        & joint["same_sign"] & (joint["passes_A"] | joint["passes_A2"])
    )
    joint.to_csv(OUT / "combined_decision.csv", index=False)


def main() -> None:
    global OUT
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["tests", "nightly", "placebos"])
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    parser.add_argument("--include-artefact-zones", action="store_true",
                        help="Use the original Analysis A city–province panel as a sensitivity.")
    args = parser.parse_args()
    manifest = pd.read_csv(OUT / "night_manifest.csv")
    bridge = pd.read_csv(OUT / "method_bridge_summary.csv")
    assert len(bridge) == 2 and bridge["n_paired_dates"].eq(6).all()
    bridge_pass = bool(bridge["joint_pass"].all())
    use_noart = not args.include_artefact_zones
    if not use_noart and args.output_dir.resolve() == OUT:
        parser.error("--include-artefact-zones requires a separate --output-dir")
    OUT = args.output_dir
    OUT.mkdir(parents=True, exist_ok=True)
    samples = analysis_samples(manifest, bridge_pass)
    print(
        f"Bridge pass={bridge_pass}: {len(samples['aug'])} August nights, "
        f"{len(samples['aug_reference'])} August references, "
        f"{len(samples['easter_reference'])} Easter references", flush=True,
    )
    if args.stage == "nightly":
        if OUT != Path(__file__).resolve().parent / "out":
            manifest.to_csv(OUT / "night_manifest.csv", index=False)
        nightly_series(paths.B0.load_nodes(), manifest, args.workers, use_noart)
    else:
        run_model_tests(manifest, samples, bridge_pass, args.workers, args.stage, use_noart)
        if args.stage == "tests":
            combined_decision()


if __name__ == "__main__":
    main()

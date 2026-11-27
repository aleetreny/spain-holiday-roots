"""Run the untabulated urban-born complement of S4 using the same stacked PPML.

Use the second term of the joint birthplace model and the S4 night samples.
Compute 95% CRV1 intervals by origin province and write only figure inputs.
"""

import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from style import A, A2, ANALYSIS, HEAVY, OUT, read

sys.path.insert(0, str(ANALYSIS / A2))
from a2_models import analysis_samples, design_blocks, paths, tests


def fit(blocks, ordinary, summer):
    ordered = [blocks[d] for d in ordinary + summer]
    holidays = [False] * len(ordinary) + [True] * len(summer)
    target, other = tests.build_equal_year_loadings(
        ordered, holidays, K=2, target_vec=np.array([0., 1.]),
        other_vecs=[np.array([1., 0.])],
    )
    estimate = tests.FastStackedMultiNight(ordered, target, other, 44).fit_unrestricted()
    fields = ("est", "se_crv1", "G", "ci_crv_lo", "ci_crv_hi", "p_crv1", "n_obs", "n_blocks", "flow_total")
    return {field: estimate[field] for field in fields}


def main():
    nodes = paths.B0.load_nodes()
    meta = tests.build_cohort_meta(nodes)
    manifest = read(A, "night_manifest")
    summer, ordinary, _ = tests._select_matched_dates(manifest, "AUG", "confirmation_O")
    dates = ordinary + summer
    flows = tests.load_panel_flows(
        HEAVY / A / "panel/panel_cities675_prov.parquet", meta, dates,
    )
    blocks = tests._prep_blocks_cohort(
        dates, flows, {year: meta for year in (2022, 2023, 2024)}, "birthplace_joint",
    )
    rows = [{"period": "2022–24", **fit({b.label: b for b in blocks}, ordinary, summer)}]
    del blocks, flows

    bridge = read(A2, "method_bridge_summary").joint_pass.all()
    sample = analysis_samples(read(A2, "night_manifest"), bridge)
    summer, ordinary = sample["aug"], sample["aug_reference"]
    blocks = design_blocks(2024, "birthplace_joint", nodes, sorted(ordinary + summer), True)
    rows.append({"period": "2025", **fit(blocks, ordinary, summer)})
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).assign(term="born_gt_10000_x_AUG").to_csv(
        OUT / "fig3_urban_contrast.csv", index=False,
    )


if __name__ == "__main__":
    main()

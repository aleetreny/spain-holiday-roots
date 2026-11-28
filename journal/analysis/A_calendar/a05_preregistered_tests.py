"""a05_preregistered_tests: Pre-registered hypothesis tests (P1-P3, S1-S6), placebo windows
diagnostic, and §7 robustness suite for Analysis A (Roots Calendar).

Parallelised across CPU cores for fast exact execution of all 9,999-draw Webb WCR score
bootstraps and exact test-inversion 95% WCR confidence intervals.

Implements PREREGISTRATION.md & design.yml verbatim:
  - Stacked PPML over holiday and ordinary (O) nights matched within calendar year,
    every FE (city, origin_prov x dest_prov) and covariate (log km, log km^2, income x dest)
    interacted with night (and arm for P3).
  - Year-specific contrasts averaged with equal weight across years.
  - Exact Kline & Santos (2012) Wild Cluster Restricted (WCR) score bootstrap, Webb 6-point
    weights, 9,999 draws, clustered by origin province jointly across all stacked nights & years
    (22 clusters pooled, 44 cohort), two-sided p-values, exact test-inversion 95% WCR CIs.
  - Holm-Bonferroni within family P (P1-P3) and within family S (S1-S6).

Run:
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    .venv/bin/python -I \
    journal/analysis/A_calendar/a05_preregistered_tests.py --workers 10
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a00_paths import B0, OUT_DIR, PANEL_DIR, V1, ensure_dirs, make_logger  # noqa: E402
from a04_nightly_series import build_cohort_meta, build_pooled_meta, load_panel_flows  # noqa: E402
from a05_solver import (  # noqa: E402
    BlockDesign,
    FastStackedMultiNight,
    build_equal_year_loadings,
    build_p3_inland_coast_loadings,
    drop_zero_groups,
    fit_single_block,
    holm_adjust,
    make_block_design,
)

N_BOOT = 9999
UNSEEN_PPML_AUG = {
    f"{yr}-08-{d:02d}" for yr in (2022, 2023) for d in range(10, 19) if d != 15
}


def webb_draws(tag: str, G: int, n_boot: int = N_BOOT) -> np.ndarray:
    return V1.rng_for("analysis_A", tag, G).choice(V1.WEBB, size=(n_boot, G))


def _prep_blocks_pooled(
    dates: list[str],
    flows: dict[str, dict[str, np.ndarray]],
    meta_by_year: dict[int, dict],
    sample: str = "all",
    equal_night_weight: bool = False,
) -> list[BlockDesign]:
    blocks = []
    for d in dates:
        yr = int(d[:4])
        m = meta_by_year[yr]
        y_raw = flows[d][sample]
        scale = 1.0 / float(y_raw.sum()) if (equal_night_weight and y_raw.sum() > 0) else 1.0
        blk = make_block_design(
            d, yr, y_raw, m["Xr"], m["cov_dist"],
            m["c_code"], m["p_code"], m["cl_code"], slope_raw=m["slope_raw"], scale_y=scale,
        )
        blk.mu_single = fit_single_block(blk, m["G_tot"])["mu"]
        blocks.append(blk)
    return blocks


def _prep_blocks_cohort(
    dates: list[str],
    flows: dict[str, dict[str, np.ndarray]],
    meta_by_year: dict[int, dict],
    model: str,
    sample: str = "all",
    equal_night_weight: bool = False,
) -> list[BlockDesign]:
    xkey = "X_COH" if model == "cohort_joint" else "X_BS"
    blocks = []
    for d in dates:
        yr = int(d[:4])
        m = meta_by_year[yr]
        y_raw = flows[d][sample]
        scale = 1.0 / float(y_raw.sum()) if (equal_night_weight and y_raw.sum() > 0) else 1.0
        blk = make_block_design(
            d, yr, y_raw, m[xkey], m["cov_dist"],
            m["c_code"], m["p_code"], m["cl_code"], slope_raw=m["slope_raw"], scale_y=scale,
        )
        blk.mu_single = fit_single_block(blk, m["G_tot"])["mu"]
        blocks.append(blk)
    return blocks


def _prep_blocks_p3(
    dates: list[str],
    flows_common: dict[str, dict[str, np.ndarray]],
    meta_common: dict,
    common_support: bool = True,
    flows_full: dict[str, dict[str, np.ndarray]] | None = None,
    meta_full: dict | None = None,
    equal_night_weight: bool = False,
) -> tuple[list[BlockDesign], list[bool]]:
    blocks: list[BlockDesign] = []
    is_inland: list[bool] = []
    for d in dates:
        yr = int(d[:4])
        if common_support:
            m = meta_common
            y_cst = flows_common[d]["coastal"]
            y_inl = flows_common[d]["inland"]
            keep = np.ones(len(y_cst), dtype=bool)
            while True:
                k1 = drop_zero_groups(np.where(keep, y_cst, 0.0), [m["c_code"], m["p_code"]]) & keep
                k2 = drop_zero_groups(np.where(k1, y_inl, 0.0), [m["c_code"], m["p_code"]]) & k1
                if np.array_equal(k2, keep):
                    break
                keep = k2
            sc_c = 1.0 / float(y_cst[keep].sum()) if equal_night_weight else 1.0
            sc_i = 1.0 / float(y_inl[keep].sum()) if equal_night_weight else 1.0
            bc = make_block_design(
                f"{d}_coast", yr, y_cst, m["Xr"], m["cov_dist"],
                m["c_code"], m["p_code"], m["cl_code"], slope_raw=m["slope_raw"],
                extra_keep=keep, scale_y=sc_c,
            )
            bi = make_block_design(
                f"{d}_inland", yr, y_inl, m["Xr"], m["cov_dist"],
                m["c_code"], m["p_code"], m["cl_code"], slope_raw=m["slope_raw"],
                extra_keep=keep, scale_y=sc_i,
            )
            assert bc.n == bi.n
            bc.mu_single = fit_single_block(bc, m["G_tot"])["mu"]
            bi.mu_single = fit_single_block(bi, m["G_tot"])["mu"]
            blocks.extend([bc, bi])
            is_inland.extend([False, True])
        else:
            assert flows_full is not None and meta_full is not None
            m = meta_full
            y_cst = flows_full[d]["coastal"]
            y_inl = flows_full[d]["inland"]
            sc_c = 1.0 / float(y_cst.sum()) if equal_night_weight else 1.0
            sc_i = 1.0 / float(y_inl.sum()) if equal_night_weight else 1.0
            bc = make_block_design(
                f"{d}_coast", yr, y_cst, m["Xr"], m["cov_dist"],
                m["c_code"], m["p_code"], m["cl_code"], slope_raw=m["slope_raw"], scale_y=sc_c,
            )
            bi = make_block_design(
                f"{d}_inland", yr, y_inl, m["Xr"], m["cov_dist"],
                m["c_code"], m["p_code"], m["cl_code"], slope_raw=m["slope_raw"], scale_y=sc_i,
            )
            bc.mu_single = fit_single_block(bc, m["G_tot"])["mu"]
            bi.mu_single = fit_single_block(bi, m["G_tot"])["mu"]
            blocks.extend([bc, bi])
            is_inland.extend([False, True])
    return blocks, is_inland


def _run_stacked_test(
    blocks: list[BlockDesign],
    l_tgt: list[np.ndarray],
    L_oth: list[np.ndarray],
    G_tot: int,
    VB: np.ndarray,
    do_ci: bool = True,
) -> dict:
    st = FastStackedMultiNight(blocks, l_tgt, L_oth, G_tot)
    unr = st.fit_unrestricted()
    p_wcr, t0, _ = st.p_wcr(0.0, VB, mu_start=unr["mus"])
    ci_lo, ci_hi = (np.nan, np.nan)
    if do_ci:
        ci_lo, ci_hi = st.invert_wcr_ci(unr["est"], unr["se_crv1"], VB, mu_start=unr["mus"])
    return {
        "est": unr["est"],
        "se_crv1": unr["se_crv1"],
        "G": unr["G"],
        "ci_crv_lo": unr["ci_crv_lo"],
        "ci_crv_hi": unr["ci_crv_hi"],
        "p_crv1": unr["p_crv1"],
        "t_wcr": t0,
        "p_wcr": p_wcr,
        "ci_wcr_lo": ci_lo,
        "ci_wcr_hi": ci_hi,
        "n_obs": unr["n_obs"],
        "n_blocks": unr["n_blocks"],
        "flow_total": unr["flow_total"],
    }


def _select_matched_dates(
    man: pd.DataFrame,
    hol_class: str,
    split_mode: str,
    o_months: set[int] | None = None,
) -> tuple[list[str], list[str], list[int]]:
    v = man[man["is_valid"]].copy()
    if split_mode == "confirmation_O":
        h_df = v[v["night_class"] == hol_class]
        o_df = v[(v["night_class"] == "O") & (v["split"] == "confirmation")]
    elif split_mode == "strict_confirmation":
        h_df = v[(v["night_class"] == hol_class) & (v["split"] == "confirmation")]
        o_df = v[(v["night_class"] == "O") & (v["split"] == "confirmation")]
    elif split_mode == "unseen_aug_vs_conf_O":
        h_df = v[(v["night_class"] == "AUG") & (v["file_date"].isin(UNSEEN_PPML_AUG))]
        o_df = v[(v["night_class"] == "O") & (v["split"] == "confirmation")]
    elif split_mode == "discovery":
        h_df = v[(v["night_class"] == hol_class) & (v["split"] == "discovery")]
        o_df = v[(v["night_class"] == "O") & (v["split"] == "discovery")]
    elif split_mode == "all_nights":
        h_df = v[v["night_class"] == hol_class]
        o_df = v[v["night_class"] == "O"]
    else:
        raise ValueError(f"Unknown split_mode: {split_mode}")

    if o_months is not None:
        o_df = o_df[o_df["evening_date"].str[5:7].astype(int).isin(o_months)]

    years = sorted(set(h_df["year"]) & set(o_df["year"]))
    h_dates = sorted(h_df.loc[h_df["year"].isin(years), "file_date"].tolist())
    o_dates = sorted(o_df.loc[o_df["year"].isin(years), "file_date"].tolist())
    return h_dates, o_dates, years


# =============================================================================
# Parallel Block Prep & Hypothesis Worker
# =============================================================================
def _prep_cohort_chunk(args):
    dates, fl_sub, c_meta, model, sample = args
    by_yr = {2022: c_meta, 2023: c_meta, 2024: c_meta}
    blks = _prep_blocks_cohort(dates, fl_sub, by_yr, model=model, sample=sample)
    return [(b.label, b) for b in blks]


_G_HYP_STATE = None


def _init_hyp_worker(state):
    global _G_HYP_STATE
    _G_HYP_STATE = state


def _eval_one_hyp(job):
    hyp, split_mode, year_filter, do_ci = job
    man, blk_p_all, blk_coh_all, blk_bs_all, p3_pair_map, VB22, VB44 = _G_HYP_STATE
    hol_map = {
        "P1": "AUG", "P2": "AUG", "P3": "AUG",
        "S1": "XMAS", "S2": "EASTER", "S3": "ALLSAINTS",
        "S4": "AUG", "S5": "AUG", "S6": "AUG",
    }
    fam = "P" if hyp.startswith("P") else "S"
    hol_class = hol_map[hyp]
    h_dates, o_dates, years = _select_matched_dates(man, hol_class, split_mode)
    if year_filter is not None:
        h_dates = [d for d in h_dates if int(d[:4]) == year_filter]
        o_dates = [d for d in o_dates if int(d[:4]) == year_filter]
        years = [year_filter] if (h_dates and o_dates) else []
    if not h_dates or not o_dates:
        return None

    t_h = time.time()
    if hyp in ("P1", "S1", "S2", "S3"):
        blks = [blk_p_all[d] for d in o_dates] + [blk_p_all[d] for d in h_dates]
        is_h = [False] * len(o_dates) + [True] * len(h_dates)
        lt, Lo = build_equal_year_loadings(blks, is_h, K=1, target_vec=np.array([1.0]))
        res = _run_stacked_test(blks, lt, Lo, 22, VB22, do_ci=do_ci)
        model_lbl = "pooled_112"
        term_lbl = f"roots_x_{hol_class}"
    elif hyp == "P3":
        blks, is_h, is_i = [], [], []
        for d in o_dates:
            bc, bi = p3_pair_map[d]; blks.extend([bc, bi]); is_h.extend([False, False]); is_i.extend([False, True])
        for d in h_dates:
            bc, bi = p3_pair_map[d]; blks.extend([bc, bi]); is_h.extend([True, True]); is_i.extend([False, True])
        lt, Lo = build_p3_inland_coast_loadings(blks, is_h, is_i)
        res = _run_stacked_test(blks, lt, Lo, 22, VB22, do_ci=do_ci)
        model_lbl = "pooled_112_common_coast19"
        term_lbl = "roots_x_AUG_x_inland"
    elif hyp in ("P2", "S5", "S6"):
        blks = [blk_coh_all[d] for d in o_dates] + [blk_coh_all[d] for d in h_dates]
        is_h = [False] * len(o_dates) + [True] * len(h_dates)
        if hyp == "P2":
            tvec, ovec = np.array([1.0, 0.0]), [np.array([0.0, 1.0])]
            term_lbl = "pre1981_x_AUG"
        elif hyp == "S5":
            tvec, ovec = np.array([0.0, 1.0]), [np.array([1.0, 0.0])]
            term_lbl = "recent_x_AUG"
        else:
            tvec, ovec = np.array([0.5, -0.5]), [np.array([0.5, 0.5])]
            term_lbl = "pre1981_minus_recent_x_AUG"
        lt, Lo = build_equal_year_loadings(blks, is_h, K=2, target_vec=tvec, other_vecs=ovec)
        res = _run_stacked_test(blks, lt, Lo, 44, VB44, do_ci=do_ci)
        model_lbl = "cohort_joint_675"
    elif hyp == "S4":
        blks = [blk_bs_all[d] for d in o_dates] + [blk_bs_all[d] for d in h_dates]
        is_h = [False] * len(o_dates) + [True] * len(h_dates)
        tvec, ovec = np.array([1.0, 0.0]), [np.array([0.0, 1.0])]
        lt, Lo = build_equal_year_loadings(blks, is_h, K=2, target_vec=tvec, other_vecs=ovec)
        res = _run_stacked_test(blks, lt, Lo, 44, VB44, do_ci=do_ci)
        model_lbl = "birthplace_joint_675"
        term_lbl = "born_le_2000_x_AUG"
    else:
        raise ValueError(hyp)

    return {
        "hypothesis": hyp,
        "family": fam,
        "split": split_mode,
        "year_scope": str(year_filter) if year_filter is not None else "equal_year_avg",
        "years_included": ",".join(str(y) for y in years),
        "n_years": len(years),
        "holiday_class": hol_class,
        "n_holiday_nights": len(h_dates),
        "n_ordinary_nights": len(o_dates),
        "model": model_lbl,
        "term": term_lbl,
        "elapsed_s": round(time.time() - t_h, 2),
        **res,
    }


def run_all_hypotheses(
    man: pd.DataFrame,
    p_meta: dict,
    p_meta_coast19: dict,
    c_meta: dict,
    p_flows: dict,
    p_flows_coast19: dict,
    c_flows: dict,
    workers: int = 10,
):
    log = make_logger()
    VB22 = webb_draws("pooled_G22", p_meta["G_tot"])
    VB44 = webb_draws("cohort_G44", c_meta["G_tot"])
    p_meta_by_yr = {yr: p_meta for yr in (2022, 2023, 2024)}

    all_valid = sorted(man.loc[man["is_valid"] & (man["night_class"].isin(["O", "AUG", "XMAS", "EASTER", "ALLSAINTS"])), "file_date"].tolist())
    log(f"Pre-building single-night BlockDesigns for {len(all_valid)} classified valid nights ({workers} workers)...")
    t_prep = time.time()
    blk_p_all = {b.label: b for b in _prep_blocks_pooled(all_valid, p_flows, p_meta_by_yr, sample="all")}

    # Parallel build of 675-city cohort & birthplace blocks
    chunks_d = [all_valid[i::workers] for i in range(workers)]
    tasks_coh = [(ch, {d: c_flows[d] for d in ch}, c_meta, "cohort_joint", "all") for ch in chunks_d if ch]
    tasks_bs = [(ch, {d: c_flows[d] for d in ch}, c_meta, "birthplace_joint", "all") for ch in chunks_d if ch]
    blk_coh_all, blk_bs_all = {}, {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(_prep_cohort_chunk, tasks_coh):
            blk_coh_all.update(res)
        for res in ex.map(_prep_cohort_chunk, tasks_bs):
            blk_bs_all.update(res)

    aug_o_dates = sorted(man.loc[man["is_valid"] & (man["night_class"].isin(["O", "AUG"])), "file_date"].tolist())
    p3_blocks_list, _ = _prep_blocks_p3(aug_o_dates, p_flows_coast19, p_meta_coast19, common_support=True)
    p3_pair_map: dict[str, tuple[BlockDesign, BlockDesign]] = {}
    for idx in range(0, len(p3_blocks_list), 2):
        p3_pair_map[aug_o_dates[idx // 2]] = (p3_blocks_list[idx], p3_blocks_list[idx + 1])
    log(f"BlockDesigns cached in {time.time() - t_prep:.1f}s")

    # Assemble all jobs (heaviest 675-city CIs scheduled first for optimal load balancing!)
    jobs = []
    for split_mode in ("all_nights", "confirmation_O", "unseen_aug_vs_conf_O", "discovery"):
        hyps = ("P2", "S4", "S5", "S6", "P1", "P3") if split_mode == "unseen_aug_vs_conf_O" else ("P2", "S4", "S5", "S6", "P1", "P3", "S1", "S2", "S3")
        for hyp in hyps:
            jobs.append((hyp, split_mode, None, True))
    for hyp in ("S1", "S2"):
        jobs.append((hyp, "strict_confirmation", None, True))
    for split_mode in ("confirmation_O", "all_nights"):
        for hyp in ("P1", "P2", "P3", "S1", "S2", "S3", "S4", "S5", "S6"):
            for yr in (2022, 2023, 2024):
                jobs.append((hyp, split_mode, yr, False))

    log(f"Running {len(jobs)} stacked PPML + WCR bootstrap hypothesis jobs across {workers} workers...")
    t_run = time.time()
    state = (man, blk_p_all, blk_coh_all, blk_bs_all, p3_pair_map, VB22, VB44)
    rows = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_hyp_worker, initargs=(state,)) as ex:
        for rec in ex.map(_eval_one_hyp, jobs):
            if rec is not None:
                rows.append(rec)
                if rec["year_scope"] == "equal_year_avg":
                    log(f"  {rec['hypothesis']:3s} [{rec['split']:20s}] "
                        f"nH={rec['n_holiday_nights']:2d} nO={rec['n_ordinary_nights']:3d} "
                        f"est={rec['est']:+.4f} [{rec['ci_wcr_lo']:+.4f}, {rec['ci_wcr_hi']:+.4f}] "
                        f"p_wcr={rec['p_wcr']:.4f} ({rec['elapsed_s']:.1f}s)")
    log(f"All {len(rows)} hypothesis jobs completed in {time.time() - t_run:.1f}s")

    df = pd.DataFrame(rows)
    # Sort neatly
    hyp_ord = {h: i for i, h in enumerate(["P1", "P2", "P3", "S1", "S2", "S3", "S4", "S5", "S6"])}
    sp_ord = {s: i for i, s in enumerate(["confirmation_O", "strict_confirmation", "unseen_aug_vs_conf_O", "discovery", "all_nights"])}
    df["_ho"] = df["hypothesis"].map(hyp_ord)
    df["_so"] = df["split"].map(sp_ord)
    df = df.sort_values(["_so", "year_scope", "_ho"]).drop(columns=["_ho", "_so"]).reset_index(drop=True)

    # Holm-Bonferroni within each (split, family) for equal_year_avg rows
    df["p_holm"] = np.nan
    for (sp, fam), idx in df[df["year_scope"] == "equal_year_avg"].groupby(["split", "family"]).groups.items():
        df.loc[idx, "p_holm"] = holm_adjust(df.loc[idx, "p_wcr"].to_numpy())

    # Composite Family-S Holm (m=6) using strict_confirmation for S1/S2 and confirmation_O for S3-S6
    df["p_holm_hybrid_S"] = np.nan
    hyb_idx = []
    for h in ("S1", "S2"):
        hyb_idx.extend(df[(df["hypothesis"] == h) & (df["split"] == "strict_confirmation") & (df["year_scope"] == "equal_year_avg")].index.tolist())
    for h in ("S3", "S4", "S5", "S6"):
        hyb_idx.extend(df[(df["hypothesis"] == h) & (df["split"] == "confirmation_O") & (df["year_scope"] == "equal_year_avg")].index.tolist())
    if len(hyb_idx) == 6:
        df.loc[hyb_idx, "p_holm_hybrid_S"] = holm_adjust(df.loc[hyb_idx, "p_wcr"].to_numpy())

    out_csv = OUT_DIR / "preregistered_tests.csv"
    df.to_csv(out_csv, index=False)
    log(f"Saved {out_csv} ({len(df)} rows)")
    return df, blk_p_all, blk_coh_all, p3_pair_map


# =============================================================================
# 2. Placebo Windows Diagnostic (PREREGISTRATION.md §3)
# =============================================================================
def _fit_pooled_slope_subset(blocks: list[BlockDesign], K: int, G_tot: int, mu_warm: list[np.ndarray] | None = None):
    l_tgt = [np.eye(K)[0] for _ in blocks]
    L_oth = [np.eye(K)[:, 1:] for _ in blocks]
    st = FastStackedMultiNight(blocks, l_tgt, L_oth, G_tot)
    if mu_warm is not None:
        st.mu_init = [m.copy() for m in mu_warm]
    L_full = [np.eye(K) for _ in blocks]
    beta, mus, _, _ = st._solve_irls(L_full)
    return beta, mus


def _placebo_one_year(args_tuple):
    split_mode, yr, aug_y, o_y, p_blks, c_blks, p3_pairs, G_p, G_c = args_tuple
    t_y = time.time()

    o_p_blks = [p_blks[d] for d in o_y]
    o_c_blks = [c_blks[d] for d in o_y]
    o_p3_c = [p3_pairs[d][0] for d in o_y]
    o_p3_i = [p3_pairs[d][1] for d in o_y]

    b_o_p1, mu_o_p1 = _fit_pooled_slope_subset(o_p_blks, 1, G_p)
    b_o_p2, mu_o_p2 = _fit_pooled_slope_subset(o_c_blks, 2, G_c)
    b_o_3c, mu_o_3c = _fit_pooled_slope_subset(o_p3_c, 1, G_p)
    b_o_3i, mu_o_3i = _fit_pooled_slope_subset(o_p3_i, 1, G_p)

    b_h_p1, _ = _fit_pooled_slope_subset([p_blks[d] for d in aug_y], 1, G_p)
    b_h_p2, _ = _fit_pooled_slope_subset([c_blks[d] for d in aug_y], 2, G_c)
    b_h_3c, _ = _fit_pooled_slope_subset([p3_pairs[d][0] for d in aug_y], 1, G_p)
    b_h_3i, _ = _fit_pooled_slope_subset([p3_pairs[d][1] for d in aug_y], 1, G_p)

    aug_p1 = float(b_h_p1[0] - b_o_p1[0])
    aug_p2 = float(b_h_p2[0] - b_o_p2[0])
    aug_p3 = float((b_h_3i[0] - b_o_3i[0]) - (b_h_3c[0] - b_o_3c[0]))

    rows = [{
        "split": split_mode, "year": yr, "window_type": "AUG", "window_idx": -1,
        "start_date": aug_y[0], "end_date": aug_y[-1],
        "n_win": len(aug_y), "n_ref": len(o_y),
        "delta_P1": aug_p1, "delta_P2": aug_p2, "delta_P3": aug_p3,
    }]

    plac_p1, plac_p2, plac_p3 = [], [], []
    n_win = len(o_y) - 9 + 1
    for k in range(n_win):
        w_idx = list(range(k, k + 9))
        r_idx = list(range(0, k)) + list(range(k + 9, len(o_y)))

        bw1, _ = _fit_pooled_slope_subset([o_p_blks[i] for i in w_idx], 1, G_p)
        br1, _ = _fit_pooled_slope_subset([o_p_blks[i] for i in r_idx], 1, G_p, mu_warm=[mu_o_p1[i] for i in r_idx])
        dp1 = float(bw1[0] - br1[0])

        bw2, _ = _fit_pooled_slope_subset([o_c_blks[i] for i in w_idx], 2, G_c)
        br2, _ = _fit_pooled_slope_subset([o_c_blks[i] for i in r_idx], 2, G_c, mu_warm=[mu_o_p2[i] for i in r_idx])
        dp2 = float(bw2[0] - br2[0])

        bw3c, _ = _fit_pooled_slope_subset([o_p3_c[i] for i in w_idx], 1, G_p)
        br3c, _ = _fit_pooled_slope_subset([o_p3_c[i] for i in r_idx], 1, G_p, mu_warm=[mu_o_3c[i] for i in r_idx])
        bw3i, _ = _fit_pooled_slope_subset([o_p3_i[i] for i in w_idx], 1, G_p)
        br3i, _ = _fit_pooled_slope_subset([o_p3_i[i] for i in r_idx], 1, G_p, mu_warm=[mu_o_3i[i] for i in r_idx])
        dp3 = float((bw3i[0] - br3i[0]) - (bw3c[0] - br3c[0]))

        plac_p1.append(dp1); plac_p2.append(dp2); plac_p3.append(dp3)
        rows.append({
            "split": split_mode, "year": yr, "window_type": "placebo", "window_idx": k,
            "start_date": o_y[k], "end_date": o_y[k + 8],
            "n_win": 9, "n_ref": len(r_idx),
            "delta_P1": dp1, "delta_P2": dp2, "delta_P3": dp3,
        })

    sum_rows = []
    for hyp, aval, pvals in (("P1", aug_p1, plac_p1), ("P2", aug_p2, plac_p2), ("P3", aug_p3, plac_p3)):
        arr = np.array(pvals)
        n_ge = int(np.sum(arr >= aval))
        rank_desc = n_ge + 1
        sum_rows.append({
            "split": split_mode, "year": yr, "hypothesis": hyp,
            "aug_delta": aval, "n_placebo_windows": len(arr),
            "n_placebo_ge_aug": n_ge, "aug_rank_of_total": f"{rank_desc}/{len(arr) + 1}",
            "empirical_upper_share": float((n_ge + 1) / (len(arr) + 1)),
            "placebo_min": float(arr.min()), "placebo_median": float(np.median(arr)),
            "placebo_max": float(arr.max()), "placebo_p95": float(np.quantile(arr, 0.95)),
        })
    return split_mode, yr, n_win, time.time() - t_y, rows, sum_rows


def run_placebo_windows(
    man: pd.DataFrame,
    blk_p_all: dict[str, BlockDesign],
    blk_coh_all: dict[str, BlockDesign],
    p3_pair_map: dict[str, tuple[BlockDesign, BlockDesign]],
    p_meta: dict,
    c_meta: dict,
    workers: int = 6,
) -> pd.DataFrame:
    log = make_logger()
    log("\n=== PLACEBO WINDOWS DIAGNOSTIC (9 consecutive O nights vs remaining O nights) ===")
    tasks = []
    for split_mode in ("confirmation_O", "all_nights"):
        for yr in (2022, 2023, 2024):
            aug_dates, o_dates, _ = _select_matched_dates(man, "AUG", split_mode)
            aug_y = sorted(d for d in aug_dates if int(d[:4]) == yr)
            o_y = sorted(d for d in o_dates if int(d[:4]) == yr)
            if len(o_y) < 18 or not aug_y:
                continue
            need = set(aug_y) | set(o_y)
            tasks.append((
                split_mode, yr, aug_y, o_y,
                {d: blk_p_all[d] for d in need},
                {d: blk_coh_all[d] for d in need},
                {d: p3_pair_map[d] for d in need},
                p_meta["G_tot"], c_meta["G_tot"],
            ))

    rows, sum_rows = [], []
    with ProcessPoolExecutor(max_workers=min(workers, max(1, len(tasks)))) as ex:
        for split_mode, yr, n_win, dt_s, r_list, s_list in ex.map(_placebo_one_year, tasks):
            rows.extend(r_list)
            sum_rows.extend(s_list)
            ranks = {s["hypothesis"]: s["aug_rank_of_total"] for s in s_list}
            log(f"  Placebo {split_mode} {yr}: {n_win} windows in {dt_s:.1f}s | AUG ranks: {ranks}")

    df_w = pd.DataFrame(rows)
    df_s = pd.DataFrame(sum_rows)
    df_w.to_csv(OUT_DIR / "placebo_windows.csv", index=False)
    df_s.to_csv(OUT_DIR / "placebo_summary.csv", index=False)
    log(f"Saved placebo_windows.csv ({len(df_w)} rows) & placebo_summary.csv ({len(df_s)} rows)")
    return df_s


# =============================================================================
# 3. Robustness Suite (PREREGISTRATION.md §7) — parallelised across workers
# =============================================================================
_G_LODO = None


def _init_lodo(payload):
    global _G_LODO
    _G_LODO = payload


def _lodo_one_prov(dp: str) -> dict:
    p_meta, c_meta, p_flows, c_flows, o_dates, aug_dates, is_h = _G_LODO
    keep_p = (p_meta["dprov_dyad"] != dp)
    blks_l = [
        make_block_design(d, int(d[:4]), p_flows[d]["all"], p_meta["Xr"], p_meta["cov_dist"],
                          p_meta["c_code"], p_meta["p_code"], p_meta["cl_code"],
                          slope_raw=p_meta["slope_raw"], extra_keep=keep_p)
        for d in o_dates + aug_dates
    ]
    lt, Lo = build_equal_year_loadings(blks_l, is_h, K=1, target_vec=np.array([1.0]))
    u_p1 = FastStackedMultiNight(blks_l, lt, Lo, p_meta["G_tot"]).fit_unrestricted()

    keep_c = (c_meta["dprov_dyad"] != dp)
    blks_cl = [
        make_block_design(d, int(d[:4]), c_flows[d]["all"], c_meta["X_COH"], c_meta["cov_dist"],
                          c_meta["c_code"], c_meta["p_code"], c_meta["cl_code"],
                          slope_raw=c_meta["slope_raw"], extra_keep=keep_c)
        for d in o_dates + aug_dates
    ]
    lt2, Lo2 = build_equal_year_loadings(blks_cl, is_h, K=2, target_vec=np.array([1.0, 0.0]), other_vecs=[np.array([0.0, 1.0])])
    u_p2 = FastStackedMultiNight(blks_cl, lt2, Lo2, c_meta["G_tot"]).fit_unrestricted()

    return {
        "dropped_dest_prov": dp,
        "is_coastal_prov": dp in p_meta["COASTAL_PROVS"],
        "P1_est": u_p1["est"], "P1_se_crv1": u_p1["se_crv1"], "P1_p_crv1": u_p1["p_crv1"],
        "P2_est": u_p2["est"], "P2_se_crv1": u_p2["se_crv1"], "P2_p_crv1": u_p2["p_crv1"],
    }


_G_ROB = None


def _init_rob(payload):
    global _G_ROB
    _G_ROB = payload


def _run_rob_task(task_spec):
    """Run one independent robustness check spec inside a worker process."""
    (nodes, man, p_meta, p_meta_coast19, c_meta, p_flows, p_flows_coast19, c_flows,
     o_dates, aug_dates, is_h, ih_p3, VB22, VB44) = _G_ROB
    split_mode = "confirmation_O"
    all_used = sorted(set(aug_dates) | set(o_dates))
    p_meta_yr = {yr: p_meta for yr in (2022, 2023, 2024)}
    c_meta_yr = {yr: c_meta for yr in (2022, 2023, 2024)}

    kind = task_spec[0]
    out = []

    def emit(check, hyp, res, note=""):
        out.append({
            "check": check, "hypothesis": hyp, "split": split_mode,
            "n_holiday_nights": res.get("n_h", len(aug_dates)),
            "n_ordinary_nights": res.get("n_o", len(o_dates)),
            "est": res["est"], "se_crv1": res["se_crv1"], "G": res["G"],
            "ci_crv_lo": res["ci_crv_lo"], "ci_crv_hi": res["ci_crv_hi"],
            "p_crv1": res["p_crv1"], "p_wcr": res["p_wcr"],
            "ci_wcr_lo": res["ci_wcr_lo"], "ci_wcr_hi": res["ci_wcr_hi"],
            "note": note,
        })

    if kind == "year_matched_P1":
        p_meta_ym = {yr: build_pooled_meta(nodes, roots_year=yr, income_year=yr) for yr in (2022, 2023, 2024)}
        blks = _prep_blocks_pooled(o_dates, p_flows, p_meta_ym, "all") + _prep_blocks_pooled(aug_dates, p_flows, p_meta_ym, "all")
        lt, Lo = build_equal_year_loadings(blks, is_h, K=1, target_vec=np.array([1.0]))
        emit("1_year_matched_roots_income", "P1", _run_stacked_test(blks, lt, Lo, 22, VB22, do_ci=True), "Census 68545 & ADRH matched to year")

    elif kind == "year_matched_P2":
        c_meta_ym = {yr: build_cohort_meta(nodes, income_year=yr) for yr in (2022, 2023, 2024)}
        blks = _prep_blocks_cohort(o_dates, c_flows, c_meta_ym, "cohort_joint", "all") + _prep_blocks_cohort(aug_dates, c_flows, c_meta_ym, "cohort_joint", "all")
        lt, Lo = build_equal_year_loadings(blks, is_h, K=2, target_vec=np.array([1.0, 0.0]), other_vecs=[np.array([0.0, 1.0])])
        emit("1_year_matched_income_cohort", "P2", _run_stacked_test(blks, lt, Lo, 44, VB44, do_ci=True), "ADRH matched to year (2024 10% microdata fixed)")

    elif kind == "road_time":
        pm_r = build_pooled_meta(nodes, use_road_time=True)
        blks = _prep_blocks_pooled(o_dates, p_flows, {yr: pm_r for yr in (2022, 2023, 2024)}, "all") + _prep_blocks_pooled(aug_dates, p_flows, {yr: pm_r for yr in (2022, 2023, 2024)}, "all")
        lt, Lo = build_equal_year_loadings(blks, is_h, K=1, target_vec=np.array([1.0]))
        emit("2_road_travel_time", "P1", _run_stacked_test(blks, lt, Lo, 22, VB22, do_ci=True), "log(road time) + log(road time)^2 from p84")

        pm_rc19 = build_pooled_meta(nodes, use_road_time=True, only_dest_provs=set(p_meta["COASTAL_PROVS"]))
        b3, i3 = _prep_blocks_p3(o_dates + aug_dates, p_flows_coast19, pm_rc19, common_support=True)
        lt3, Lo3 = build_p3_inland_coast_loadings(b3, ih_p3, i3)
        emit("2_road_travel_time", "P3", _run_stacked_test(b3, lt3, Lo3, 22, VB22, do_ci=True), "road time on common 19 coastal provinces")

    elif kind == "drop_orig":
        _, lbl, drop_set = task_spec
        pm_d = build_pooled_meta(nodes, drop_orig_provs=drop_set)
        pf_d = load_panel_flows(PANEL_DIR / "panel_cities137_prov.parquet", pm_d, dates=all_used)
        blks_d = _prep_blocks_pooled(o_dates, pf_d, {yr: pm_d for yr in (2022, 2023, 2024)}, "all") + _prep_blocks_pooled(aug_dates, pf_d, {yr: pm_d for yr in (2022, 2023, 2024)}, "all")
        lt, Lo = build_equal_year_loadings(blks_d, is_h, K=1, target_vec=np.array([1.0]))
        VB_d = webb_draws(f"rob_{lbl}_p", pm_d["G_tot"])
        emit(lbl, "P1", _run_stacked_test(blks_d, lt, Lo, pm_d["G_tot"], VB_d, do_ci=True), f"dropped orig provs {sorted(drop_set)}")

        cm_d = build_cohort_meta(nodes, drop_orig_provs=drop_set)
        cf_d = load_panel_flows(PANEL_DIR / "panel_cities675_prov.parquet", cm_d, dates=all_used)
        blks_cd = _prep_blocks_cohort(o_dates, cf_d, {yr: cm_d for yr in (2022, 2023, 2024)}, "cohort_joint", "all") + _prep_blocks_cohort(aug_dates, cf_d, {yr: cm_d for yr in (2022, 2023, 2024)}, "cohort_joint", "all")
        lt2, Lo2 = build_equal_year_loadings(blks_cd, is_h, K=2, target_vec=np.array([1.0, 0.0]), other_vecs=[np.array([0.0, 1.0])])
        VB_cd = webb_draws(f"rob_{lbl}_c", cm_d["G_tot"])
        emit(lbl, "P2", _run_stacked_test(blks_cd, lt2, Lo2, cm_d["G_tot"], VB_cd, do_ci=True), f"dropped orig provs {sorted(drop_set)}")

        if lbl == "3_drop_MAD_and_BCN":
            pm_dc19 = build_pooled_meta(nodes, drop_orig_provs=drop_set, only_dest_provs=set(p_meta["COASTAL_PROVS"]))
            pf_dc19 = load_panel_flows(PANEL_DIR / "panel_cities137_prov.parquet", pm_dc19, dates=all_used)
            b3, i3 = _prep_blocks_p3(o_dates + aug_dates, pf_dc19, pm_dc19, common_support=True)
            lt3, Lo3 = build_p3_inland_coast_loadings(b3, ih_p3, i3)
            emit(lbl, "P3", _run_stacked_test(b3, lt3, Lo3, pm_dc19["G_tot"], VB_d, do_ci=True), "dropped MAD & BCN")

    elif kind == "dest_sample":
        _, samp_lbl, samp_key = task_spec
        blks_s = _prep_blocks_pooled(o_dates, p_flows, p_meta_yr, samp_key) + _prep_blocks_pooled(aug_dates, p_flows, p_meta_yr, samp_key)
        lt, Lo = build_equal_year_loadings(blks_s, is_h, K=1, target_vec=np.array([1.0]))
        emit(samp_lbl, "P1", _run_stacked_test(blks_s, lt, Lo, 22, VB22, do_ci=True), f"sample={samp_key}")

        blks_cs = _prep_blocks_cohort(o_dates, c_flows, c_meta_yr, "cohort_joint", samp_key) + _prep_blocks_cohort(aug_dates, c_flows, c_meta_yr, "cohort_joint", samp_key)
        lt2, Lo2 = build_equal_year_loadings(blks_cs, is_h, K=2, target_vec=np.array([1.0, 0.0]), other_vecs=[np.array([0.0, 1.0])])
        emit(samp_lbl, "P2", _run_stacked_test(blks_cs, lt2, Lo2, 44, VB44, do_ci=True), f"sample={samp_key}")

    elif kind == "p3_full":
        b3_f, i3_f = _prep_blocks_p3(o_dates + aug_dates, p_flows_coast19, p_meta_coast19, common_support=False, flows_full=p_flows, meta_full=p_meta)
        lt3f, Lo3f = build_p3_inland_coast_loadings(b3_f, ih_p3, i3_f)
        emit("5_P3_full_support_47_vs_19", "P3", _run_stacked_test(b3_f, lt3f, Lo3f, 22, VB22, do_ci=True), "all 47 inland zones vs 19 coastal zones")

    elif kind == "zero_rule":
        _, zrule, r112, lbl = task_spec
        cm_z = build_cohort_meta(nodes, zero_rule=zrule, restrict_to_112=r112)
        cf_z = load_panel_flows(PANEL_DIR / "panel_cities675_prov.parquet", cm_z, dates=all_used) if r112 else c_flows
        cm_z_yr = {yr: cm_z for yr in (2022, 2023, 2024)}
        VB_z = VB22 if r112 else VB44
        blks_z = _prep_blocks_cohort(o_dates, cf_z, cm_z_yr, "cohort_joint", "all") + _prep_blocks_cohort(aug_dates, cf_z, cm_z_yr, "cohort_joint", "all")
        for hyp_z, tv, ov in [
            ("P2", np.array([1.0, 0.0]), [np.array([0.0, 1.0])]),
            ("S5", np.array([0.0, 1.0]), [np.array([1.0, 0.0])]),
            ("S6", np.array([0.5, -0.5]), [np.array([0.5, 0.5])]),
        ]:
            lt_z, Lo_z = build_equal_year_loadings(blks_z, is_h, K=2, target_vec=tv, other_vecs=ov)
            emit(lbl, hyp_z, _run_stacked_test(blks_z, lt_z, Lo_z, cm_z["G_tot"], VB_z, do_ci=(hyp_z == "P2" or r112)), f"zero_rule={zrule}, r112={r112}")
        blks_zb = _prep_blocks_cohort(o_dates, cf_z, cm_z_yr, "birthplace_joint", "all") + _prep_blocks_cohort(aug_dates, cf_z, cm_z_yr, "birthplace_joint", "all")
        lt_zb, Lo_zb = build_equal_year_loadings(blks_zb, is_h, K=2, target_vec=np.array([1.0, 0.0]), other_vecs=[np.array([0.0, 1.0])])
        emit(lbl, "S4", _run_stacked_test(blks_zb, lt_zb, Lo_zb, cm_z["G_tot"], VB_z, do_ci=r112), f"zero_rule={zrule}, r112={r112}")

    elif kind == "seasonal_js":
        h_js, o_js, _ = _select_matched_dates(man, "AUG", split_mode, o_months={6, 9})
        is_h_js = [False] * len(o_js) + [True] * len(h_js)
        bp = _prep_blocks_pooled(o_js, p_flows, p_meta_yr, "all") + _prep_blocks_pooled(h_js, p_flows, p_meta_yr, "all")
        lt, Lo = build_equal_year_loadings(bp, is_h_js, K=1, target_vec=np.array([1.0]))
        r = _run_stacked_test(bp, lt, Lo, 22, VB22, do_ci=True); r["n_h"], r["n_o"] = len(h_js), len(o_js)
        emit("7_seasonal_ref_Jun_Sep_only", "P1", r, "O restricted to June & September evenings")

        bc = _prep_blocks_cohort(o_js, c_flows, c_meta_yr, "cohort_joint", "all") + _prep_blocks_cohort(h_js, c_flows, c_meta_yr, "cohort_joint", "all")
        lt2, Lo2 = build_equal_year_loadings(bc, is_h_js, K=2, target_vec=np.array([1.0, 0.0]), other_vecs=[np.array([0.0, 1.0])])
        r2 = _run_stacked_test(bc, lt2, Lo2, 44, VB44, do_ci=True); r2["n_h"], r2["n_o"] = len(h_js), len(o_js)
        emit("7_seasonal_ref_Jun_Sep_only", "P2", r2, "O restricted to June & September evenings")

        b3, i3 = _prep_blocks_p3(o_js + h_js, p_flows_coast19, p_meta_coast19, common_support=True)
        ih3 = [False] * (2 * len(o_js)) + [True] * (2 * len(h_js))
        lt3, Lo3 = build_p3_inland_coast_loadings(b3, ih3, i3)
        r3 = _run_stacked_test(b3, lt3, Lo3, 22, VB22, do_ci=True); r3["n_h"], r3["n_o"] = len(h_js), len(o_js)
        emit("7_seasonal_ref_Jun_Sep_only", "P3", r3, "O restricted to June & September evenings")

    elif kind == "equal_night":
        bp = _prep_blocks_pooled(o_dates, p_flows, p_meta_yr, "all", equal_night_weight=True) + _prep_blocks_pooled(aug_dates, p_flows, p_meta_yr, "all", equal_night_weight=True)
        lt, Lo = build_equal_year_loadings(bp, is_h, K=1, target_vec=np.array([1.0]))
        emit("8_equal_weight_by_night", "P1", _run_stacked_test(bp, lt, Lo, 22, VB22, do_ci=True), "each night rescaled to sum(y)=1")

        bc = _prep_blocks_cohort(o_dates, c_flows, c_meta_yr, "cohort_joint", "all", equal_night_weight=True) + _prep_blocks_cohort(aug_dates, c_flows, c_meta_yr, "cohort_joint", "all", equal_night_weight=True)
        lt2, Lo2 = build_equal_year_loadings(bc, is_h, K=2, target_vec=np.array([1.0, 0.0]), other_vecs=[np.array([0.0, 1.0])])
        emit("8_equal_weight_by_night", "P2", _run_stacked_test(bc, lt2, Lo2, 44, VB44, do_ci=True), "each night rescaled to sum(y)=1")

        b3, i3 = _prep_blocks_p3(o_dates + aug_dates, p_flows_coast19, p_meta_coast19, common_support=True, equal_night_weight=True)
        lt3, Lo3 = build_p3_inland_coast_loadings(b3, ih_p3, i3)
        emit("8_equal_weight_by_night", "P3", _run_stacked_test(b3, lt3, Lo3, 22, VB22, do_ci=True), "each night x arm rescaled to sum(y)=1")

    elif kind == "single_pooled_yr":
        bp = _prep_blocks_pooled(o_dates, p_flows, p_meta_yr, "all") + _prep_blocks_pooled(aug_dates, p_flows, p_meta_yr, "all")
        lt, Lo = build_equal_year_loadings(bp, is_h, K=1, target_vec=np.array([1.0]), pooled_across_years=True)
        emit("8_single_pooled_delta_across_yrs", "P1", _run_stacked_test(bp, lt, Lo, 22, VB22, do_ci=True), "single common beta_0 and delta across 2022-24")

        bc = _prep_blocks_cohort(o_dates, c_flows, c_meta_yr, "cohort_joint", "all") + _prep_blocks_cohort(aug_dates, c_flows, c_meta_yr, "cohort_joint", "all")
        lt2, Lo2 = build_equal_year_loadings(bc, is_h, K=2, target_vec=np.array([1.0, 0.0]), other_vecs=[np.array([0.0, 1.0])], pooled_across_years=True)
        emit("8_single_pooled_delta_across_yrs", "P2", _run_stacked_test(bc, lt2, Lo2, 44, VB44, do_ci=True), "single common beta_0 and delta across 2022-24")

        b3, i3 = _prep_blocks_p3(o_dates + aug_dates, p_flows_coast19, p_meta_coast19, common_support=True)
        lt3, Lo3 = build_p3_inland_coast_loadings(b3, ih_p3, i3, pooled_across_years=True)
        emit("8_single_pooled_delta_across_yrs", "P3", _run_stacked_test(b3, lt3, Lo3, 22, VB22, do_ci=True), "single common delta across 2022-24")

    elif kind == "drop_3art":
        pf_na = load_panel_flows(PANEL_DIR / "panel_cities137_prov.parquet", p_meta, dates=all_used, use_noart=True)
        blks_na = _prep_blocks_pooled(o_dates, pf_na, p_meta_yr, "all") + _prep_blocks_pooled(aug_dates, pf_na, p_meta_yr, "all")
        lt, Lo = build_equal_year_loadings(blks_na, is_h, K=1, target_vec=np.array([1.0]))
        emit("9_drop_3_artefact_dest_zones", "P1", _run_stacked_test(blks_na, lt, Lo, 22, VB22, do_ci=True), "excluding 31905_AM, 19103_AM, 40149_AM as destinations")

    return out


def run_robustness(
    nodes: pd.DataFrame,
    man: pd.DataFrame,
    p_meta: dict,
    p_meta_coast19: dict,
    c_meta: dict,
    p_flows: dict,
    p_flows_coast19: dict,
    c_flows: dict,
    workers: int = 10,
) -> pd.DataFrame:
    log = make_logger()
    log(f"\n=== ROBUSTNESS SUITE (PREREGISTRATION.md §7, {workers} workers) ===")
    split_mode = "confirmation_O"
    aug_dates, o_dates, _ = _select_matched_dates(man, "AUG", split_mode)
    is_h = [False] * len(o_dates) + [True] * len(aug_dates)
    ih_p3 = [False] * (2 * len(o_dates)) + [True] * (2 * len(aug_dates))
    VB22 = webb_draws("rob_pooled_G22", 22)
    VB44 = webb_draws("rob_cohort_G44", 44)

    # 1. Parallel §7 checks
    tasks = [
        ("year_matched_P2",),
        ("drop_orig", "3_drop_MAD_and_BCN", {"28", "08"}),
        ("drop_orig", "3_drop_BCN_only", {"08"}),
        ("drop_orig", "3_drop_MAD_only", {"28"}),
        ("dest_sample", "5_inland_provinces_only_28", "inlandprov"),
        ("dest_sample", "5_inland_zones_all47", "inland"),
        ("dest_sample", "5_coastal_zones_19", "coastal"),
        ("zero_rule", "plus_1.0", False, "6_zero_plus_1.0"),
        ("zero_rule", "plus_0.1", False, "6_zero_plus_0.1"),
        ("zero_rule", "share_plus_1e4", False, "6_zero_share_plus_1e4"),
        ("zero_rule", "share_plus_1e5", False, "6_zero_share_plus_1e5"),
        ("zero_rule", "asinh_1000_share", False, "6_zero_asinh_1000_share"),
        ("zero_rule", "haldane_0.5", True, "6_zero_112_large_cities"),
        ("equal_night",),
        ("single_pooled_yr",),
        ("seasonal_js",),
        ("year_matched_P1",),
        ("road_time",),
        ("p3_full",),
        ("drop_3art",),
    ]
    t_rob = time.time()
    payload = (nodes, man, p_meta, p_meta_coast19, c_meta, p_flows, p_flows_coast19, c_flows,
               o_dates, aug_dates, is_h, ih_p3, VB22, VB44)
    rob_rows = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_rob, initargs=(payload,)) as ex:
        for res_list in ex.map(_run_rob_task, tasks):
            for r in res_list:
                rob_rows.append(r)
                ci_str = f"[{r['ci_wcr_lo']:+.4f}, {r['ci_wcr_hi']:+.4f}]" if np.isfinite(r["ci_wcr_lo"]) else f"[{r['ci_crv_lo']:+.4f}, {r['ci_crv_hi']:+.4f}]*"
                log(f"  [ROB] {r['check']:34s} {r['hypothesis']:3s}: est={r['est']:+.4f} {ci_str} p_wcr={r['p_wcr']:.4f}")
    df_rob = pd.DataFrame(rob_rows).sort_values(["check", "hypothesis"]).reset_index(drop=True)
    df_rob.to_csv(OUT_DIR / "robustness.csv", index=False)
    log(f"Saved robustness.csv ({len(df_rob)} rows, {time.time() - t_rob:.1f}s)")

    # 2. Parallel LODO (47 destination provinces)
    log(f"  Running leave-one-destination-province-out (47 provinces, {workers} workers) for P1 & P2...")
    t_lodo = time.time()
    lodo_payload = (p_meta, c_meta, p_flows, c_flows, o_dates, aug_dates, is_h)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_lodo, initargs=(lodo_payload,)) as ex:
        lodo_rows = list(ex.map(_lodo_one_prov, list(p_meta["PROVS"])))
    df_lodo = pd.DataFrame(lodo_rows)
    df_lodo.to_csv(OUT_DIR / "robustness_lodo_dest.csv", index=False)
    log(f"  LODO done ({time.time() - t_lodo:.1f}s) | "
        f"P1 range: [{df_lodo['P1_est'].min():+.4f} (drop {df_lodo.loc[df_lodo['P1_est'].idxmin(), 'dropped_dest_prov']}), "
        f"{df_lodo['P1_est'].max():+.4f} (drop {df_lodo.loc[df_lodo['P1_est'].idxmax(), 'dropped_dest_prov']})] | "
        f"P2 range: [{df_lodo['P2_est'].min():+.4f} (drop {df_lodo.loc[df_lodo['P2_est'].idxmin(), 'dropped_dest_prov']}), "
        f"{df_lodo['P2_est'].max():+.4f} (drop {df_lodo.loc[df_lodo['P2_est'].idxmax(), 'dropped_dest_prov']})]")
    return df_rob


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-main", action="store_true")
    ap.add_argument("--skip-placebo", action="store_true")
    ap.add_argument("--skip-robustness", action="store_true")
    ap.add_argument("--workers", type=int, default=10)
    args = ap.parse_args()

    t0 = time.time()
    ensure_dirs()
    log = make_logger(t0)
    nodes = B0.load_nodes()
    man = pd.read_csv(OUT_DIR / "night_manifest.csv")

    p_meta = build_pooled_meta(nodes)
    p_meta_coast19 = build_pooled_meta(nodes, only_dest_provs=set(p_meta["COASTAL_PROVS"]))
    c_meta = build_cohort_meta(nodes)

    needed = sorted(man.loc[man["is_valid"] & man["night_class"].isin(["O", "AUG", "XMAS", "EASTER", "ALLSAINTS"]), "file_date"].tolist())
    log(f"Loading panel flows for {len(needed)} classified valid nights...")
    p_flows = load_panel_flows(PANEL_DIR / "panel_cities137_prov.parquet", p_meta, dates=needed)
    p_flows_coast19 = load_panel_flows(PANEL_DIR / "panel_cities137_prov.parquet", p_meta_coast19, dates=needed)
    c_flows = load_panel_flows(PANEL_DIR / "panel_cities675_prov.parquet", c_meta, dates=needed)

    blk_p_all = blk_coh_all = p3_pair_map = None
    if not args.skip_main:
        _, blk_p_all, blk_coh_all, p3_pair_map = run_all_hypotheses(
            man, p_meta, p_meta_coast19, c_meta, p_flows, p_flows_coast19, c_flows, workers=args.workers
        )
    if not args.skip_placebo:
        assert blk_p_all is not None
        run_placebo_windows(man, blk_p_all, blk_coh_all, p3_pair_map, p_meta, c_meta, workers=min(args.workers, 6))
    if not args.skip_robustness:
        run_robustness(nodes, man, p_meta, p_meta_coast19, c_meta, p_flows, p_flows_coast19, c_flows, workers=args.workers)
    log(f"ALL PRE-REGISTERED TESTS & DIAGNOSTICS COMPLETE in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()

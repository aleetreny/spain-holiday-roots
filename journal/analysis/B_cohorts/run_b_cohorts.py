"""Analysis B: Hardening the arrival-cohort and birthplace-size results for EPJ Data Science.

Fast multi-core staged pipeline (uses ProcessPoolExecutor with 12-14 workers on pure-NumPy
tasks after DuckDB data loading closes cleanly; every WCR p-value and CI uses 9,999 Webb draws).

Stages:
  --stage b1_b2         B.1 audit comparison vs V1 + B.2 main 2024/2022/2023 contrasts (both baselines)
  --stage b2_multinight B.2 PREREGISTRATION §3 multi-night stacked PPML (9 Aug vs 13 Feb & 12 Oct, 2022-24)
  --stage b2_sisters    B.2 pairwise sister-night grids (2024 Aug x Feb13, Aug x Oct12, Nov x Feb13, 2022/23)
  --stage b3            B.3 direct 4-block stacked inland x coast test (19 common provs + 47 vs 19)
  --stage b4            B.4 age-matched cohort & birthplace-size models within 43-64 and 65+ bands
  --stage b5            B.5 Table A2 zero-cell rules on a common predicted-flow contrast scale
  --stage b6            B.6 44-origin-prov jackknife, 47-dest-prov jackknife, drop Madrid/Barcelona/both
  --stage assemble      Assemble tidy plotting table out/fig_who_pulls_tidy.csv
  --stage all           Run all stages sequentially (each stage parallelized across CPU cores!), then assemble.

Run:
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    /Users/alejandrotreny/Proyectos/spain-holiday-roots/.venv/bin/python -I \
    journal/analysis/B_cohorts/run_b_cohorts.py --stage all
"""
from __future__ import annotations

import argparse
import gc
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))

import b_common as BC  # noqa: E402
from b_common import V1, C, MAIN_ROOT, JOURNAL_DATA_BD, N_BOOT_WCR, WEBB  # noqa: E402

T0 = time.time()
N_WORKERS = min(14, max(1, (os.cpu_count() or 4) - 2))

# Shared read-only globals populated BEFORE forking worker pool (copy-on-write on macOS/POSIX)
_G: dict = {}


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def run_pool(fn, items: list, n_workers: int = N_WORKERS) -> list:
    """Execute fn over items using forked workers (zero serialization overhead on read-only _G)."""
    gc.collect()
    nw = min(n_workers, max(1, len(items)))
    if nw <= 1:
        return [fn(it) for it in items]
    ctx_mp = mp.get_context("fork")
    with ProcessPoolExecutor(max_workers=nw, mp_context=ctx_mp) as ex:
        return list(ex.map(fn, items))


def verify_multiblock_vs_pyfixest(ctx: dict, FLOWS: dict, X_coh: np.ndarray) -> None:
    """Assert FastStackedMultiBlock matches pyfixest.fepois on both 2-night diff and 4-block inland x coast."""
    import pyfixest as pf

    no, _, _ = BC.make_design(ctx, FLOWS["2024-02-14"]["all"], X_coh)
    nh, _, _ = BC.make_design(ctx, FLOWS["2024-08-15"]["all"], X_coh)
    fo, fh = V1.fit_single(no, ctx["G_TOT"]), V1.fit_single(nh, ctx["G_TOT"])
    st_diff = BC.build_two_night_diff_stack(nh, no, 0, 1, fh["mu"], fo["mu"])
    res_diff = st_diff.fit_unrestricted(ctx["G_TOT"])
    expected_diff = float((fh["beta"][0] - fo["beta"][0]) - (fh["beta"][1] - fo["beta"][1]))
    dIF = (fh["IF"][:, 0] - fo["IF"][:, 0]) - (fh["IF"][:, 1] - fo["IF"][:, 1])
    Geff = ctx["G_TOT"]
    expected_se = float(np.sqrt(Geff / (Geff - 1) * np.sum(dIF ** 2)))
    assert abs(res_diff["beta"][0] - expected_diff) < 1e-7
    assert abs(res_diff["se"][0] - expected_se) < 1e-7

    m19 = ctx["cell_coastal_prov"]
    no_c, _, G19 = BC.make_design(ctx, FLOWS["2024-02-14"]["coastal"], X_coh, mask=m19)
    nh_c, _, _ = BC.make_design(ctx, FLOWS["2024-08-15"]["coastal"], X_coh, mask=m19)
    no_i, _, _ = BC.make_design(ctx, FLOWS["2024-02-14"]["inland"], X_coh, mask=m19)
    nh_i, _, _ = BC.make_design(ctx, FLOWS["2024-08-15"]["inland"], X_coh, mask=m19)
    foc, fhc = V1.fit_single(no_c, G19), V1.fit_single(nh_c, G19)
    foi, fhi = V1.fit_single(no_i, G19), V1.fit_single(nh_i, G19)
    st4 = BC.build_four_block_inland_coast_stack(no_c, nh_c, no_i, nh_i, foc["mu"], fhc["mu"], foi["mu"], fhi["mu"])
    r4 = st4.fit_unrestricted(G19)
    gam_exp = (fhi["beta"] - foi["beta"]) - (fhc["beta"] - foc["beta"])
    assert np.allclose(r4["beta"][:2], gam_exp, atol=1e-7)

    provs19 = ctx["COASTAL_PROVS"]
    P2I = ctx["P2I"]
    ri19, pj19 = ctx["ri_all"][m19], ctx["pj_all"][m19]
    linc19 = ctx["LINC"][ri19]
    lkm19, lkm2_19 = ctx["COV_DIST"][m19, 0], ctx["COV_DIST"][m19, 1]
    blocks = [
        ("0_oc", FLOWS["2024-02-14"]["coastal"][m19], 0, 0),
        ("1_hc", FLOWS["2024-08-15"]["coastal"][m19], 1, 0),
        ("2_oi", FLOWS["2024-02-14"]["inland"][m19], 0, 1),
        ("3_hi", FLOWS["2024-08-15"]["inland"][m19], 1, 1),
    ]
    rows_pf = []
    for btag, y_b, hol_f, inl_f in blocks:
        d_b = pd.DataFrame({
            "y": y_b,
            "x0": X_coh[m19, 0], "x1": X_coh[m19, 1],
            "x0_h": X_coh[m19, 0] * hol_f, "x1_h": X_coh[m19, 1] * hol_f,
            "x0_i": X_coh[m19, 0] * inl_f, "x1_i": X_coh[m19, 1] * inl_f,
            "x0_hi": X_coh[m19, 0] * hol_f * inl_f, "x1_hi": X_coh[m19, 1] * hol_f * inl_f,
            "c_fe": btag + "_" + ctx["MMUN"][ri19],
            "p_fe": btag + "_" + ctx["MPROV"][ri19] + "_" + ctx["PROVS"][pj19],
            "cl": ctx["MPROV"][ri19],
            "b_lkm": lkm19, "b_lkm2": lkm2_19, "btag": btag,
        })
        for p_idx in range(1, len(provs19)):
            pcode = provs19[p_idx]
            d_b[f"inc_{pcode}"] = linc19 * (pj19 == P2I[pcode])
        rows_pf.append(d_b)
    df_pf = pd.concat(rows_pf, ignore_index=True)
    cov_terms = []
    for btag, _, _, _ in blocks:
        is_b = (df_pf["btag"] == btag).astype(float)
        df_pf[f"d1_{btag}"] = df_pf["b_lkm"] * is_b
        df_pf[f"d2_{btag}"] = df_pf["b_lkm2"] * is_b
        cov_terms.extend([f"d1_{btag}", f"d2_{btag}"])
        for p_idx in range(1, len(provs19)):
            pcode = provs19[p_idx]
            col = f"inc_{pcode}_{btag}"
            df_pf[col] = df_pf[f"inc_{pcode}"] * is_b
            cov_terms.append(col)
    rhs = "x0_hi + x1_hi + x0_h + x1_h + x0_i + x1_i + x0 + x1 + " + " + ".join(cov_terms)
    m_pf = pf.fepois(f"y ~ {rhs} | c_fe + p_fe", data=df_pf, vcov={"CRV1": "cl"}, iwls_tol=1e-10, iwls_maxiter=100)
    t_pf = m_pf.tidy()
    for k_i, nm in enumerate(("x0_hi", "x1_hi")):
        assert abs(float(t_pf.loc[nm, "Estimate"]) - r4["beta"][k_i]) < 1e-5
    log("Verified FastStackedMultiBlock (2-night diff & 4-block inland x coast) against pyfixest (|diff| < 1e-5)")


def run_two_night_contrast(
    nh: V1.SingleNightDesign, fh: dict,
    no: V1.SingleNightDesign, fo: dict,
    tnames: list[str],
    VB: np.ndarray,
    do_ci: bool = True,
    diff_pair: tuple[int, int, str] | None = None,
    do_diff_ci: bool = True,
) -> list[dict]:
    st = BC.build_two_night_stack(nh, no, fh["mu"], fo["mu"])
    dIF = fh["IF"] - fo["IF"]
    Geff = int((np.abs(dIF[:, 0]) > 0).sum())
    se_d = np.sqrt(Geff / (Geff - 1) * np.sum(dIF ** 2, axis=0))
    delta = fh["beta"] - fo["beta"]
    tc = stats.t.ppf(0.975, Geff - 1)

    rows = []
    for k_i, tnm in enumerate(tnames):
        pw, _ = st.p_param(k_i, 0.0, VB)
        if do_ci:
            lo, hi = BC.invert_ci_safe(
                lambda v, m0: st.p_param(k_i, v, VB, m0, return_stats=True),
                float(delta[k_i]), float(se_d[k_i]),
            )
        else:
            lo, hi = np.nan, np.nan
        rows.append({
            "term": tnm, "n": nh.n + no.n, "G": Geff,
            "beta_holiday": float(fh["beta"][k_i]), "beta_ordinary": float(fo["beta"][k_i]),
            "delta": float(delta[k_i]), "se_crv1": float(se_d[k_i]),
            "p_crv1": float(2.0 * stats.t.sf(abs(delta[k_i] / se_d[k_i]), Geff - 1)),
            "ci_crv_lo": float(delta[k_i] - tc * se_d[k_i]), "ci_crv_hi": float(delta[k_i] + tc * se_d[k_i]),
            "p_wcr": pw, "ci_wcr_lo": lo, "ci_wcr_hi": hi,
        })

    if diff_pair is not None:
        ip, ineg, dname = diff_pair
        st_d = BC.build_two_night_diff_stack(nh, no, ip, ineg, fh["mu"], fo["mu"])
        d_est = float(delta[ip] - delta[ineg])
        d_if = dIF[:, ip] - dIF[:, ineg]
        d_se = float(np.sqrt(Geff / (Geff - 1) * np.sum(d_if ** 2)))
        pw, _ = st_d.p_param(0, 0.0, VB)
        if do_ci and do_diff_ci:
            lo, hi = BC.invert_ci_safe(
                lambda v, m0: st_d.p_param(0, v, VB, m0, return_stats=True), d_est, d_se
            )
        else:
            lo, hi = np.nan, np.nan
        rows.append({
            "term": dname, "n": nh.n + no.n, "G": Geff,
            "beta_holiday": float(fh["beta"][ip] - fh["beta"][ineg]),
            "beta_ordinary": float(fo["beta"][ip] - fo["beta"][ineg]),
            "delta": d_est, "se_crv1": d_se,
            "p_crv1": float(2.0 * stats.t.sf(abs(d_est / d_se), Geff - 1)),
            "ci_crv_lo": float(d_est - tc * d_se), "ci_crv_hi": float(d_est + tc * d_se),
            "p_wcr": pw, "ci_wcr_lo": lo, "ci_wcr_hi": hi,
        })
    return rows


def init_shared_data(load_v1dest: bool = False):
    """Load all census & cached MITMA data cleanly into _G BEFORE any fork()."""
    if _G.get("initialized"):
        return
    ctx = BC.load_census_and_design()
    ns = BC.cached_night_sets()
    G_TOT = ctx["G_TOT"]
    VB_9999 = BC.rng_for("b_webb_9999").choice(WEBB, size=(N_BOOT_WCR, G_TOT))
    MATS, NTOT = ctx["MATS"], ctx["NTOT"]
    ri_all, pj_all = ctx["ri_all"], ctx["pj_all"]
    lsh05 = lambda col: np.log((MATS[col] + 0.5) / NTOT)[ri_all, pj_all]  # noqa: E731
    X_COH = np.column_stack([lsh05("n_pre1981"), lsh05("n_rec2016")])
    X_BS = np.column_stack([lsh05("n_vil2k_dated"), lsh05("n_urb10k_dated")])
    X_BS_RAW = np.column_stack([lsh05("n_vil2k_raw"), lsh05("n_urb10k_raw")])
    MODELS = {
        "cohort_joint": (X_COH, ["pre1981", "recent_2016_2023"], (0, 1, "diff_pre1981_minus_recent")),
        "birthsize_joint": (X_BS, ["born_village_le2k", "born_urban_gt10k"], (0, 1, "diff_village_minus_urban")),
    }
    FLOWS = BC.load_mitma_flows(ctx, ns["all_needed"], exclude_artefact_dest=True)
    for k, src in [
        ("2024-02-mean13", ns["feb24_sis"]),
        ("2024-08-mean9", ns["aug24_sis"]),
        ("2024-10-mean12", ns["oct24_clean12"]),
        ("2024-10-mean15", ns["oct24_all15"]),
        ("2024-10-mean10_noVal", ns["oct24_no_val10"]),
        ("2022-02-mean12", ns["feb22_sis"]),
        ("2022-08-mean9", ns["aug22_sis"]),
        ("2023-02-mean12", ns["feb23_sis"]),
        ("2023-08-mean9", ns["aug23_sis"]),
    ]:
        BC.add_window_mean(FLOWS, k, src)

    FLOWS_V1DEST = None
    if load_v1dest:
        FLOWS_V1DEST = BC.load_mitma_flows(
            ctx, sorted(set(["2024-02-14", "2024-08-15"] + ns["oct24_all15"])), exclude_artefact_dest=False
        )
        BC.add_window_mean(FLOWS_V1DEST, "2024-10-mean15", ns["oct24_all15"])
        BC.add_window_mean(FLOWS_V1DEST, "2024-10-mean12", ns["oct24_clean12"])

    _G.update({
        "initialized": True, "ctx": ctx, "ns": ns, "G_TOT": G_TOT, "VB_9999": VB_9999,
        "X_COH": X_COH, "X_BS": X_BS, "X_BS_RAW": X_BS_RAW, "MODELS": MODELS,
        "FLOWS": FLOWS, "FLOWS_V1DEST": FLOWS_V1DEST,
    })
    gc.collect()


# =============================================================================
# Worker functions (top-level for pickling / fork compatibility)
# =============================================================================
def _w_fit_single(task):
    mname, sample, d = task
    ctx = _G["ctx"]
    Xmat = _G["MODELS"][mname][0]
    nt, _, _ = BC.make_design(ctx, _G["FLOWS"][d][sample], Xmat)
    f = V1.fit_single(nt, _G["G_TOT"])
    return (mname, sample, d), (nt, f)


def _w_b1_task(task):
    vname, use_v1dest, use_raw_bs, o_date, v1_cname, sample, mname = task
    ctx, G_TOT, VB = _G["ctx"], _G["G_TOT"], _G["VB_9999"]
    fdict = _G["FLOWS_V1DEST"] if use_v1dest else _G["FLOWS"]
    if mname == "cohort_joint":
        Xmat, tnames = _G["X_COH"], ["pre1981", "recent_2016_2023"]
    else:
        Xmat = _G["X_BS_RAW"] if use_raw_bs else _G["X_BS"]
        tnames = ["born_village_le2k", "born_urban_gt10k"]
    v1_csv = _G["v1_csv"]
    v1_map = {
        ("cohort_joint", "pre1981"): "pre1981_exodus",
        ("cohort_joint", "recent_2016_2023"): "recent_2016_2024",
        ("birthsize_joint", "born_village_le2k"): "born_village_le2k",
        ("birthsize_joint", "born_urban_gt10k"): "born_urban_gt10k",
    }
    no, _, _ = BC.make_design(ctx, fdict[o_date][sample], Xmat)
    nh, _, _ = BC.make_design(ctx, fdict["2024-08-15"][sample], Xmat)
    fo, fh = V1.fit_single(no, G_TOT), V1.fit_single(nh, G_TOT)
    st = BC.build_two_night_stack(nh, no, fh["mu"], fo["mu"])
    dIF = fh["IF"] - fo["IF"]
    se_d = np.sqrt(G_TOT / (G_TOT - 1) * np.sum(dIF ** 2, axis=0))
    out = []
    for k_i, tnm in enumerate(tnames):
        delta = float(fh["beta"][k_i] - fo["beta"][k_i])
        pw, _ = st.p_param(k_i, 0.0, VB)
        v1_row = v1_csv[
            (v1_csv["contrast"] == v1_cname) & (v1_csv["sample"] == sample)
            & (v1_csv["model"] == mname) & (v1_csv["term"] == v1_map[(mname, tnm)])
        ].iloc[0]
        v1_delta = float(v1_row["delta"])
        if vname in ("v1_exact_repro", "v1_exact_repro_oct15"):
            assert abs(delta - v1_delta) < 1e-8
        out.append({
            "variant": vname, "ordinary": o_date, "sample": sample, "model": mname, "term": tnm,
            "beta_ordinary": float(fo["beta"][k_i]), "beta_aug15": float(fh["beta"][k_i]),
            "delta": delta, "se_crv1": float(se_d[k_i]), "p_wcr": pw,
            "v1_contrast_ref": v1_cname, "v1_delta": v1_delta, "diff_vs_v1": delta - v1_delta,
        })
    return out


def _w_b2_contrast(task):
    cname, h, o, sample, mname, do_ci = task
    _, tnames, dpair = _G["MODELS"][mname]
    nh, fh = _G["FITS"][(mname, sample, h)]
    no, fo = _G["FITS"][(mname, sample, o)]
    res_l = run_two_night_contrast(nh, fh, no, fo, tnames, _G["VB_9999"], do_ci=do_ci, diff_pair=dpair)
    for r in res_l:
        r.update({"contrast": cname, "holiday": h, "ordinary": o, "sample": sample, "model": mname})
    return res_l


def _w_b2_multinight(task):
    sname, h_dates, o_dates, sample, mname, do_ci = task
    G_TOT, VB = _G["G_TOT"], _G["VB_9999"]
    _, tnames, (ip, ineg, dname) = _G["MODELS"][mname]
    d_ord = [_G["FITS"][(mname, sample, d)][0] for d in o_dates]
    d_hol = [_G["FITS"][(mname, sample, d)][0] for d in h_dates]
    mu0_l = [_G["FITS"][(mname, sample, d)][1]["mu"] for d in o_dates] + [_G["FITS"][(mname, sample, d)][1]["mu"] for d in h_dates]
    stm = BC.build_multinight_stack(d_ord, d_hol, mu0_l)
    fu = stm.fit_unrestricted(G_TOT)
    k = len(tnames)
    tc = stats.t.ppf(0.975, fu["G"] - 1)
    out = []
    for k_i, tnm in enumerate(tnames):
        est, se = float(fu["beta"][k_i]), float(fu["se"][k_i])
        pw, _ = stm.p_param(k_i, 0.0, VB, mu0=fu["mu"])
        if do_ci:
            lo, hi = BC.invert_ci_safe(lambda v, m0: stm.p_param(k_i, v, VB, m0, return_stats=True), est, se)
        else:
            lo, hi = np.nan, np.nan
        out.append({
            "spec": sname, "n_hol_nights": len(h_dates), "n_ord_nights": len(o_dates),
            "sample": sample, "model": mname, "term": tnm,
            "n_stacked": fu["n"], "G": fu["G"],
            "beta_ordinary_pooled": float(fu["beta"][k + k_i]),
            "beta_holiday_pooled": float(fu["beta"][k + k_i] + est),
            "delta": est, "se_crv1": se,
            "p_crv1": float(2.0 * stats.t.sf(abs(est / se), fu["G"] - 1)),
            "ci_crv_lo": float(est - tc * se), "ci_crv_hi": float(est + tc * se),
            "p_wcr": pw, "ci_wcr_lo": lo, "ci_wcr_hi": hi,
        })
    V_hol, V_com = stm.V[:, :k], stm.V[:, k:]
    V_diff = np.column_stack([V_hol[:, ip], V_hol[:, ip] + V_hol[:, ineg], V_com])
    stm_d = BC.FastStackedMultiBlock(stm.designs, V_diff, mu0_warm=fu["mu"])
    d_est = float(fu["beta"][ip] - fu["beta"][ineg])
    d_if = fu["IF"][:, ip] - fu["IF"][:, ineg]
    d_se = float(np.sqrt(fu["G"] / (fu["G"] - 1) * np.sum(d_if ** 2)))
    pw_d, _ = stm_d.p_param(0, 0.0, VB, mu0=fu["mu"])
    if do_ci:
        lo_d, hi_d = BC.invert_ci_safe(lambda v, m0: stm_d.p_param(0, v, VB, m0, return_stats=True), d_est, d_se)
    else:
        lo_d, hi_d = np.nan, np.nan
    out.append({
        "spec": sname, "n_hol_nights": len(h_dates), "n_ord_nights": len(o_dates),
        "sample": sample, "model": mname, "term": dname,
        "n_stacked": fu["n"], "G": fu["G"],
        "beta_ordinary_pooled": float(fu["beta"][k + ip] - fu["beta"][k + ineg]),
        "beta_holiday_pooled": float((fu["beta"][k + ip] + fu["beta"][ip]) - (fu["beta"][k + ineg] + fu["beta"][ineg])),
        "delta": d_est, "se_crv1": d_se,
        "p_crv1": float(2.0 * stats.t.sf(abs(d_est / d_se), fu["G"] - 1)),
        "ci_crv_lo": float(d_est - tc * d_se), "ci_crv_hi": float(d_est + tc * d_se),
        "p_wcr": pw_d, "ci_wcr_lo": lo_d, "ci_wcr_hi": hi_d,
    })
    return out


def _w_b2_sister_slice(task):
    gname, sample, mname, h, o_list = task
    _, tnames, _ = _G["MODELS"][mname]
    VB = _G["VB_9999"]
    nh, fh = _G["FITS"][(mname, sample, h)]
    rows = []
    for o in o_list:
        no, fo = _G["FITS"][(mname, sample, o)]
        st = BC.build_two_night_stack(nh, no, fh["mu"], fo["mu"])
        dIF = fh["IF"] - fo["IF"]
        Geff = int((np.abs(dIF[:, 0]) > 0).sum())
        se_d = np.sqrt(Geff / (Geff - 1) * np.sum(dIF ** 2, axis=0))
        delta = fh["beta"] - fo["beta"]
        for k_i, tnm in enumerate(tnames):
            pw, _ = st.p_param(k_i, 0.0, VB)
            rows.append({
                "grid": gname, "sample": sample, "model": mname, "term": tnm,
                "holiday": h, "ordinary": o,
                "delta": float(delta[k_i]), "se_crv1": float(se_d[k_i]),
                "p_crv1": float(2.0 * stats.t.sf(abs(delta[k_i] / se_d[k_i]), Geff - 1)),
                "p_wcr": pw,
            })
    return rows


def _w_b3_task(task):
    support_name, use_mask19_inl, use_stratum_cov, strict_cell_intersect, cname, h, o, mname, do_ci = task
    ctx, FLOWS, G_TOT, VB = _G["ctx"], _G["FLOWS"], _G["G_TOT"], _G["VB_9999"]
    Xmat, tnames, (ip, ineg, dname) = _G["MODELS"][mname]
    m19 = ctx["cell_coastal_prov"]
    use_mask_inl = m19 if use_mask19_inl else None
    cov_i = ctx["COV_DIST_INL"] if use_stratum_cov else ctx["COV_DIST"]
    cov_c = ctx["COV_DIST_CST"] if use_stratum_cov else ctx["COV_DIST"]

    no_c, _, G_c = BC.make_design(ctx, FLOWS[o]["coastal"], Xmat, mask=m19, cov_dist=cov_c)
    nh_c, _, _ = BC.make_design(ctx, FLOWS[h]["coastal"], Xmat, mask=m19, cov_dist=cov_c)
    no_i, _, G_i = BC.make_design(ctx, FLOWS[o]["inland"], Xmat, mask=use_mask_inl, cov_dist=cov_i)
    nh_i, _, _ = BC.make_design(ctx, FLOWS[h]["inland"], Xmat, mask=use_mask_inl, cov_dist=cov_i)
    assert G_c == G_TOT and G_i == G_TOT

    if strict_cell_intersect:
        keep4 = no_c.keep & nh_c.keep & no_i.keep & nh_i.keep
        m_exact = np.zeros(len(m19), dtype=bool)
        m_exact[np.where(m19)[0][keep4]] = True
        no_c, _, _ = BC.make_design(ctx, FLOWS[o]["coastal"], Xmat, mask=m_exact, cov_dist=cov_c)
        nh_c, _, _ = BC.make_design(ctx, FLOWS[h]["coastal"], Xmat, mask=m_exact, cov_dist=cov_c)
        no_i, _, _ = BC.make_design(ctx, FLOWS[o]["inland"], Xmat, mask=m_exact, cov_dist=cov_i)
        nh_i, _, _ = BC.make_design(ctx, FLOWS[h]["inland"], Xmat, mask=m_exact, cov_dist=cov_i)
        assert no_c.n == nh_c.n == no_i.n == nh_i.n

    foc, fhc = V1.fit_single(no_c, G_TOT), V1.fit_single(nh_c, G_TOT)
    foi, fhi = V1.fit_single(no_i, G_TOT), V1.fit_single(nh_i, G_TOT)
    st4 = BC.build_four_block_inland_coast_stack(no_c, nh_c, no_i, nh_i, foc["mu"], fhc["mu"], foi["mu"], fhi["mu"])

    d_inl = fhi["beta"] - foi["beta"]
    d_cst = fhc["beta"] - foc["beta"]
    gamma = d_inl - d_cst
    IF_gam = (fhi["IF"] - foi["IF"]) - (fhc["IF"] - foc["IF"])
    Geff = int((np.abs(IF_gam[:, 0]) > 0).sum())
    se_gam = np.sqrt(Geff / (Geff - 1) * np.sum(IF_gam ** 2, axis=0))
    tc = stats.t.ppf(0.975, Geff - 1)

    rows = []
    for k_i, tnm in enumerate(tnames):
        pw, _ = st4.p_param(k_i, 0.0, VB)
        lo, hi = BC.invert_ci_safe(
            lambda v, m0: st4.p_param(k_i, v, VB, m0, return_stats=True),
            float(gamma[k_i]), float(se_gam[k_i]),
        ) if do_ci else (np.nan, np.nan)
        rows.append({
            "support": support_name, "contrast": cname, "holiday": h, "ordinary": o,
            "model": mname, "term": tnm,
            "n_4block": st4.N, "n_inland_ord": no_i.n, "n_inland_hol": nh_i.n,
            "n_coast_ord": no_c.n, "n_coast_hol": nh_c.n, "G": Geff,
            "delta_inland": float(d_inl[k_i]), "delta_coast": float(d_cst[k_i]),
            "gamma_inland_minus_coast": float(gamma[k_i]), "se_crv1": float(se_gam[k_i]),
            "p_crv1": float(2.0 * stats.t.sf(abs(gamma[k_i] / se_gam[k_i]), Geff - 1)),
            "ci_crv_lo": float(gamma[k_i] - tc * se_gam[k_i]),
            "ci_crv_hi": float(gamma[k_i] + tc * se_gam[k_i]),
            "p_wcr": pw, "ci_wcr_lo": lo, "ci_wcr_hi": hi,
        })

    st4_d = BC.build_four_block_diff_stack(no_c, nh_c, no_i, nh_i, ip, ineg, foc["mu"], fhc["mu"], foi["mu"], fhi["mu"])
    g_diff = float(gamma[ip] - gamma[ineg])
    if_d = IF_gam[:, ip] - IF_gam[:, ineg]
    se_d = float(np.sqrt(Geff / (Geff - 1) * np.sum(if_d ** 2)))
    pw_d, _ = st4_d.p_param(0, 0.0, VB)
    lo_d, hi_d = BC.invert_ci_safe(
        lambda v, m0: st4_d.p_param(0, v, VB, m0, return_stats=True), g_diff, se_d
    ) if do_ci else (np.nan, np.nan)
    rows.append({
        "support": support_name, "contrast": cname, "holiday": h, "ordinary": o,
        "model": mname, "term": dname,
        "n_4block": st4.N, "n_inland_ord": no_i.n, "n_inland_hol": nh_i.n,
        "n_coast_ord": no_c.n, "n_coast_hol": nh_c.n, "G": Geff,
        "delta_inland": float(d_inl[ip] - d_inl[ineg]),
        "delta_coast": float(d_cst[ip] - d_cst[ineg]),
        "gamma_inland_minus_coast": g_diff, "se_crv1": se_d,
        "p_crv1": float(2.0 * stats.t.sf(abs(g_diff / se_d), Geff - 1)),
        "ci_crv_lo": float(g_diff - tc * se_d), "ci_crv_hi": float(g_diff + tc * se_d),
        "p_wcr": pw_d, "ci_wcr_lo": lo_d, "ci_wcr_hi": hi_d,
    })
    return rows


def _w_b4_task(task):
    spec_idx, cname, o_date, sample = task
    aspec, X_a, tnames, use_agecity, diff_pairs = _G["AGE_MODELS"][spec_idx]
    ctx, FLOWS, G_TOT, VB = _G["ctx"], _G["FLOWS"], _G["G_TOT"], _G["VB_9999"]
    sf_a = ctx["SLOPE_AGECITY"] if use_agecity else ctx["SLOPE_MAIN"]

    no_a, _, _ = BC.make_design(ctx, FLOWS[o_date][sample], X_a, slope_raw=sf_a)
    nh_a, _, _ = BC.make_design(ctx, FLOWS["2024-08-15"][sample], X_a, slope_raw=sf_a)
    fo_a, fh_a = V1.fit_single(no_a, G_TOT), V1.fit_single(nh_a, G_TOT)
    st_a = BC.build_two_night_stack(nh_a, no_a, fh_a["mu"], fo_a["mu"])
    dIF = fh_a["IF"] - fo_a["IF"]
    Geff = int((np.abs(dIF[:, 0]) > 0).sum())
    se_d = np.sqrt(Geff / (Geff - 1) * np.sum(dIF ** 2, axis=0))
    delta = fh_a["beta"] - fo_a["beta"]
    tc = stats.t.ppf(0.975, Geff - 1)

    rows = []
    for k_i, tnm in enumerate(tnames):
        pw, _ = st_a.p_param(k_i, 0.0, VB)
        lo, hi = BC.invert_ci_safe(
            lambda v, m0: st_a.p_param(k_i, v, VB, m0, return_stats=True),
            float(delta[k_i]), float(se_d[k_i]),
        )
        rows.append({
            "spec": aspec, "contrast": cname, "sample": sample, "term": tnm,
            "is_within_band_diff": False, "n": nh_a.n + no_a.n, "G": Geff,
            "beta_ordinary": float(fo_a["beta"][k_i]), "beta_aug15": float(fh_a["beta"][k_i]),
            "delta": float(delta[k_i]), "se_crv1": float(se_d[k_i]),
            "p_crv1": float(2.0 * stats.t.sf(abs(delta[k_i] / se_d[k_i]), Geff - 1)),
            "ci_crv_lo": float(delta[k_i] - tc * se_d[k_i]),
            "ci_crv_hi": float(delta[k_i] + tc * se_d[k_i]),
            "p_wcr": pw, "ci_wcr_lo": lo, "ci_wcr_hi": hi,
        })
    for ip, ineg, dname in diff_pairs:
        st_d = BC.build_two_night_diff_stack(nh_a, no_a, ip, ineg, fh_a["mu"], fo_a["mu"])
        d_est = float(delta[ip] - delta[ineg])
        d_if = dIF[:, ip] - dIF[:, ineg]
        d_se = float(np.sqrt(Geff / (Geff - 1) * np.sum(d_if ** 2)))
        pw_d, _ = st_d.p_param(0, 0.0, VB)
        lo_d, hi_d = BC.invert_ci_safe(
            lambda v, m0: st_d.p_param(0, v, VB, m0, return_stats=True), d_est, d_se
        )
        rows.append({
            "spec": aspec, "contrast": cname, "sample": sample, "term": dname,
            "is_within_band_diff": True, "n": nh_a.n + no_a.n, "G": Geff,
            "beta_ordinary": float(fo_a["beta"][ip] - fo_a["beta"][ineg]),
            "beta_aug15": float(fh_a["beta"][ip] - fh_a["beta"][ineg]),
            "delta": d_est, "se_crv1": d_se,
            "p_crv1": float(2.0 * stats.t.sf(abs(d_est / d_se), Geff - 1)),
            "ci_crv_lo": float(d_est - tc * d_se), "ci_crv_hi": float(d_est + tc * d_se),
            "p_wcr": pw_d, "ci_wcr_lo": lo_d, "ci_wcr_hi": hi_d,
        })
    return rows


def _w_b5_task(task):
    cname, o_date, sample, r_idx, mtype, cols, tnames = task
    rname, gfn, smask = _G["ZERO_RULES"][r_idx]
    ctx, FLOWS, RAW_SHARE_STATS = _G["ctx"], _G["FLOWS"], _G["RAW_SHARE_STATS"]
    MATS, NTOT, ri_all, pj_all = ctx["MATS"], ctx["NTOT"], ctx["ri_all"], ctx["pj_all"]
    N_cell = NTOT[ri_all, 0]

    m_idx = np.asarray(smask, dtype=bool)
    _, cl_u = pd.factorize(ctx["MPROV"][ri_all][m_idx], sort=True)
    G_sub = len(cl_u)
    VB_sub = BC.rng_for("b5_webb", G_sub).choice(WEBB, size=(N_BOOT_WCR, G_sub))
    hd_sub = V1.FastHDFE2(pd.factorize(ctx["c_code"][m_idx])[0], pd.factorize(ctx["p_code"][m_idx])[0])
    ones_sub = np.ones(int(m_idx.sum()))
    N_sub = N_cell[m_idx]

    y_ord_s = FLOWS[o_date][sample][m_idx]
    y_aug_s = FLOWS["2024-08-15"][sample][m_idx]
    w_flow = 0.5 * (y_ord_s / y_ord_s.sum() + y_aug_s / y_aug_s.sum())

    s_cols = [(MATS[c] / NTOT)[ri_all, pj_all][m_idx] for c in cols]
    g_cols = [gfn(sc, N_sub) for sc in s_cols]
    X_full = np.column_stack([gfn((MATS[c] / NTOT)[ri_all, pj_all], N_cell) for c in cols])
    no_s, _, _ = BC.make_design(ctx, FLOWS[o_date][sample], X_full, mask=m_idx)
    nh_s, _, _ = BC.make_design(ctx, FLOWS["2024-08-15"][sample], X_full, mask=m_idx)
    fo_s, fh_s = V1.fit_single(no_s, G_sub), V1.fit_single(nh_s, G_sub)
    st_s = BC.build_two_night_stack(nh_s, no_s, fh_s["mu"], fo_s["mu"])
    dIF = fh_s["IF"] - fo_s["IF"]
    se_d = np.sqrt(G_sub / (G_sub - 1) * np.sum(dIF ** 2, axis=0))
    delta = fh_s["beta"] - fo_s["beta"]

    rows = []
    for k_i, (col, tnm) in enumerate(zip(cols, tnames)):
        pw, _ = st_s.p_param(k_i, 0.0, VB_sub)
        lo, hi = BC.invert_ci_safe(
            lambda v, m0: st_s.p_param(k_i, v, VB_sub, m0, return_stats=True),
            float(delta[k_i]), float(se_d[k_i]),
        )
        ds_w = RAW_SHARE_STATS[col]["sd_within_fe"]
        ds_raw = RAW_SHARE_STATS[col]["sd_raw"]
        dg_w_vec = gfn(s_cols[k_i] + ds_w, N_sub) - g_cols[k_i]
        dg_raw_vec = gfn(s_cols[k_i] + ds_raw, N_sub) - g_cols[k_i]
        dg_w_fw = float(np.sum(w_flow * dg_w_vec))
        dg_w_unw = float(np.mean(dg_w_vec))
        dg_raw_fw = float(np.sum(w_flow * dg_raw_vec))
        g_dm = hd_sub.demean(g_cols[k_i], ones_sub, tol=1e-11)
        sd_g_w = float(np.std(g_dm))

        dg_aug, dg_ord = dg_w_vec[nh_s.keep], dg_w_vec[no_s.keep]
        mult_aug = float(np.sum(fh_s["mu"] * np.exp(fh_s["beta"][k_i] * dg_aug)) / np.sum(fh_s["mu"]))
        mult_ord = float(np.sum(fo_s["mu"] * np.exp(fo_s["beta"][k_i] * dg_ord)) / np.sum(fo_s["mu"]))
        pct_agg_exact = 100.0 * (mult_aug / mult_ord - 1.0)

        pct_from = lambda d_val, scale: float(100.0 * (np.exp(d_val * scale) - 1.0))  # noqa: E731
        rows.append({
            "contrast": cname, "sample": sample, "rule": rname, "model": mtype, "term": tnm,
            "G": G_sub, "n_aug": nh_s.n, "n_ord": no_s.n,
            "zero_cell_share": float(np.mean(MATS[col][ri_all, pj_all][m_idx] == 0)),
            "beta_ordinary": float(fo_s["beta"][k_i]), "beta_aug15": float(fh_s["beta"][k_i]),
            "delta_raw": float(delta[k_i]), "se_crv1": float(se_d[k_i]),
            "p_wcr": pw, "ci_wcr_lo": lo, "ci_wcr_hi": hi,
            "ds_within_raw_share": ds_w,
            "dg_flow_weighted_1sd_within_share": dg_w_fw,
            "dg_unweighted_1sd_within_share": dg_w_unw,
            "sd_within_fe_g": sd_g_w,
            "pct_flow_1sd_rawshare_fw": pct_from(delta[k_i], dg_w_fw),
            "pct_flow_1sd_rawshare_fw_ci_lo": pct_from(lo, dg_w_fw),
            "pct_flow_1sd_rawshare_fw_ci_hi": pct_from(hi, dg_w_fw),
            "pct_flow_1sd_rawshare_agg_exact": pct_agg_exact,
            "pct_flow_1sd_rawshare_unw": pct_from(delta[k_i], dg_w_unw),
            "pct_flow_1sd_rawshare_unw_ci_lo": pct_from(lo, dg_w_unw),
            "pct_flow_1sd_rawshare_unw_ci_hi": pct_from(hi, dg_w_unw),
            "pct_flow_1sd_overall_share_fw": pct_from(delta[k_i], dg_raw_fw),
            "pct_flow_1sd_overall_share_fw_ci_lo": pct_from(lo, dg_raw_fw),
            "pct_flow_1sd_overall_share_fw_ci_hi": pct_from(hi, dg_raw_fw),
            "pct_flow_1sd_g_within": pct_from(delta[k_i], sd_g_w),
            "pct_flow_1sd_g_within_ci_lo": pct_from(lo, sd_g_w),
            "pct_flow_1sd_g_within_ci_hi": pct_from(hi, sd_g_w),
        })
    return rows


def _w_b6_origin_jk(task):
    cl_drop_idx, p_drop = task
    ctx, FLOWS, G_TOT, VB = _G["ctx"], _G["FLOWS"], _G["G_TOT"], _G["VB_9999"]
    n_muns_p = int((ctx["MPROV"] == p_drop).sum())
    rows = []
    for cname, o_date in (("Aug15_vs_Feb14", "2024-02-14"), ("Aug15_vs_OctMean12", "2024-10-mean12")):
        for sample in ("all", "inland"):
            for mname, (Xmat, tnames, _) in _G["MODELS"].items():
                no_j, _, _ = BC.make_design(ctx, FLOWS[o_date][sample], Xmat, drop_cl=cl_drop_idx)
                nh_j, _, _ = BC.make_design(ctx, FLOWS["2024-08-15"][sample], Xmat, drop_cl=cl_drop_idx)
                fo_j, fh_j = V1.fit_single(no_j, G_TOT), V1.fit_single(nh_j, G_TOT)
                st_j = BC.build_two_night_stack(nh_j, no_j, fh_j["mu"], fo_j["mu"])
                dIF = fh_j["IF"] - fo_j["IF"]
                Geff = int((np.abs(dIF[:, 0]) > 0).sum())
                assert Geff == 43
                se_d = np.sqrt(Geff / (Geff - 1) * np.sum(dIF ** 2, axis=0))
                delta = fh_j["beta"] - fo_j["beta"]
                for k_i, tnm in enumerate(tnames):
                    pw, _ = st_j.p_param(k_i, 0.0, VB)
                    rows.append({
                        "contrast": cname, "sample": sample, "model": mname, "term": tnm,
                        "dropped_prov": p_drop, "dropped_n_munis": n_muns_p, "G_rem": Geff,
                        "beta_ordinary": float(fo_j["beta"][k_i]), "beta_aug15": float(fh_j["beta"][k_i]),
                        "delta": float(delta[k_i]), "se_crv1": float(se_d[k_i]), "p_wcr": pw,
                    })
    return rows


def _w_b6_drop_mb(task):
    dlabel, drop_provs, cname, o_date, sample, mname = task
    ctx, FLOWS = _G["ctx"], _G["FLOWS"]
    Xmat, tnames, dpair = _G["MODELS"][mname]
    oprov_cell = ctx["MPROV"][ctx["ri_all"]]
    dmask = ~np.isin(oprov_cell, drop_provs)

    if sample != "common_19_provs_inland_x_coast":
        no_m, _, G_m = BC.make_design(ctx, FLOWS[o_date][sample], Xmat, mask=dmask)
        nh_m, _, _ = BC.make_design(ctx, FLOWS["2024-08-15"][sample], Xmat, mask=dmask)
        fo_m, fh_m = V1.fit_single(no_m, G_m), V1.fit_single(nh_m, G_m)
        VB_m = BC.rng_for("b6_drop_webb", dlabel, G_m).choice(WEBB, size=(N_BOOT_WCR, G_m))
        res_l = run_two_night_contrast(nh_m, fh_m, no_m, fo_m, tnames, VB_m, do_ci=True, diff_pair=dpair)
        for r in res_l:
            r.update({"drop_spec": dlabel, "contrast": cname, "sample": sample, "model": mname, "estimand": "delta_aug_ord"})
        return res_l

    m_19_drop = dmask & ctx["cell_coastal_prov"]
    no_c, _, G_m = BC.make_design(ctx, FLOWS[o_date]["coastal"], Xmat, mask=m_19_drop)
    nh_c, _, _ = BC.make_design(ctx, FLOWS["2024-08-15"]["coastal"], Xmat, mask=m_19_drop)
    no_i, _, _ = BC.make_design(ctx, FLOWS[o_date]["inland"], Xmat, mask=m_19_drop)
    nh_i, _, _ = BC.make_design(ctx, FLOWS["2024-08-15"]["inland"], Xmat, mask=m_19_drop)
    foc, fhc = V1.fit_single(no_c, G_m), V1.fit_single(nh_c, G_m)
    foi, fhi = V1.fit_single(no_i, G_m), V1.fit_single(nh_i, G_m)
    st4 = BC.build_four_block_inland_coast_stack(no_c, nh_c, no_i, nh_i, foc["mu"], fhc["mu"], foi["mu"], fhi["mu"])
    gamma = (fhi["beta"] - foi["beta"]) - (fhc["beta"] - foc["beta"])
    IF_g = (fhi["IF"] - foi["IF"]) - (fhc["IF"] - foc["IF"])
    se_g = np.sqrt(G_m / (G_m - 1) * np.sum(IF_g ** 2, axis=0))
    tc = stats.t.ppf(0.975, G_m - 1)
    VB_m = BC.rng_for("b6_drop_webb", dlabel, G_m).choice(WEBB, size=(N_BOOT_WCR, G_m))
    rows = []
    for k_i, tnm in enumerate(tnames):
        pw, _ = st4.p_param(k_i, 0.0, VB_m)
        lo, hi = BC.invert_ci_safe(
            lambda v, m0: st4.p_param(k_i, v, VB_m, m0, return_stats=True),
            float(gamma[k_i]), float(se_g[k_i]),
        )
        rows.append({
            "drop_spec": dlabel, "contrast": cname, "sample": "common_19_provs_inland_x_coast",
            "model": mname, "estimand": "gamma_inland_minus_coast", "term": tnm,
            "n": st4.N, "G": G_m,
            "beta_holiday": float(fhi["beta"][k_i] - foi["beta"][k_i]),
            "beta_ordinary": float(fhc["beta"][k_i] - foc["beta"][k_i]),
            "delta": float(gamma[k_i]), "se_crv1": float(se_g[k_i]),
            "p_crv1": float(2.0 * stats.t.sf(abs(gamma[k_i] / se_g[k_i]), G_m - 1)),
            "ci_crv_lo": float(gamma[k_i] - tc * se_g[k_i]),
            "ci_crv_hi": float(gamma[k_i] + tc * se_g[k_i]),
            "p_wcr": pw, "ci_wcr_lo": lo, "ci_wcr_hi": hi,
        })
    return rows


def _w_b6_dest_jk(task):
    p_idx, p_code_str = task
    ctx, FLOWS, VB = _G["ctx"], _G["FLOWS"], _G["VB_9999"]
    dmask = ctx["pj_all"] != p_idx
    rows = []
    for cname, o_date in (("Aug15_vs_Feb14", "2024-02-14"), ("Aug15_vs_OctMean12", "2024-10-mean12")):
        for sample in ("all", "inland"):
            for mname, (Xmat, tnames, _) in _G["MODELS"].items():
                no_d, _, G_d = BC.make_design(ctx, FLOWS[o_date][sample], Xmat, mask=dmask)
                nh_d, _, _ = BC.make_design(ctx, FLOWS["2024-08-15"][sample], Xmat, mask=dmask)
                fo_d, fh_d = V1.fit_single(no_d, G_d), V1.fit_single(nh_d, G_d)
                st_d = BC.build_two_night_stack(nh_d, no_d, fh_d["mu"], fo_d["mu"])
                dIF = fh_d["IF"] - fo_d["IF"]
                se_d = np.sqrt(G_d / (G_d - 1) * np.sum(dIF ** 2, axis=0))
                delta = fh_d["beta"] - fo_d["beta"]
                for k_i, tnm in enumerate(tnames):
                    pw, _ = st_d.p_param(k_i, 0.0, VB)
                    rows.append({
                        "contrast": cname, "sample": sample, "model": mname, "term": tnm,
                        "dropped_dest_prov": p_code_str, "G": G_d,
                        "delta": float(delta[k_i]), "se_crv1": float(se_d[k_i]), "p_wcr": pw,
                    })
    return rows


# =============================================================================
# Stage drivers
# =============================================================================
def stage_b1_b2():
    log(f"=== STAGE b1_b2 (workers={N_WORKERS}) ===")
    init_shared_data(load_v1dest=True)
    verify_multiblock_vs_pyfixest(_G["ctx"], _G["FLOWS"], _G["X_COH"])

    # B.1 Audit diff vs V1
    _G["v1_csv"] = pd.read_csv(MAIN_ROOT / "results" / "v1_c2_main_contrasts.csv")
    b1_tasks = []
    for vname, use_v1d, use_raw_bs, o_date, v1_cname in [
        ("v1_exact_repro", True, True, "2024-02-14", "Aug15_vs_Feb14"),
        ("drop_7777_only", True, False, "2024-02-14", "Aug15_vs_Feb14"),
        ("drop_artefact_dest_only", False, True, "2024-02-14", "Aug15_vs_Feb14"),
        ("prereg_binding_both", False, False, "2024-02-14", "Aug15_vs_Feb14"),
        ("v1_exact_repro_oct15", True, True, "2024-10-mean15", "Aug15_vs_OctMean15"),
        ("prereg_binding_oct15", False, False, "2024-10-mean15", "Aug15_vs_OctMean15"),
        ("prereg_binding_oct12_noDANA", False, False, "2024-10-mean12", "Aug15_vs_OctMean15"),
    ]:
        for sample in ("all", "inland", "coastal"):
            for mname in ("cohort_joint", "birthsize_joint"):
                b1_tasks.append((vname, use_v1d, use_raw_bs, o_date, v1_cname, sample, mname))
    b1_rows = [r for sub in run_pool(_w_b1_task, b1_tasks) for r in sub]
    pd.DataFrame(b1_rows).to_csv(OUT / "b1_v1_audit_comparison.csv", index=False)
    log(f"B.1 complete ({len(b1_rows)} rows) -> out/b1_v1_audit_comparison.csv")

    # Fit key dates in parallel
    KEY_DATES = sorted({
        "2024-02-14", "2024-08-15", "2024-11-01", "2024-12-25", "2024-03-29",
        "2024-02-mean13", "2024-08-mean9", "2024-10-mean12", "2024-10-mean15", "2024-10-mean10_noVal",
        "2022-02-16", "2022-08-15", "2022-11-01", "2022-02-mean12", "2022-08-mean9",
        "2023-02-15", "2023-08-15", "2023-02-mean12", "2023-08-mean9",
    })
    fit_tasks = [(m, s, d) for m in _G["MODELS"] for s in ("all", "inland", "coastal", "inlandprov") for d in KEY_DATES]
    FITS = dict(run_pool(_w_fit_single, fit_tasks))
    _G["FITS"] = FITS
    sn_rows = []
    for (mname, sample, d), (nt, f) in FITS.items():
        for k_i, tnm in enumerate(_G["MODELS"][mname][1]):
            sn_rows.append({
                "model": mname, "sample": sample, "date": d, "term": tnm,
                "n": nt.n, "G": f["G"], "beta": float(f["beta"][k_i]), "se_crv1": float(f["se"][k_i]),
            })
    pd.DataFrame(sn_rows).to_csv(OUT / "b2_key_single_nights.csv", index=False)

    MAIN_CONTRASTS = [
        ("Aug15_vs_Feb14", "2024-08-15", "2024-02-14", ("all", "inland", "coastal", "inlandprov"), True),
        ("Aug15_vs_OctMean12", "2024-08-15", "2024-10-mean12", ("all", "inland", "coastal", "inlandprov"), True),
        ("Aug15_vs_FebMean13", "2024-08-15", "2024-02-mean13", ("all", "inland", "coastal"), True),
        ("Aug15_vs_OctMean15", "2024-08-15", "2024-10-mean15", ("all", "inland", "coastal"), False),
        ("Aug15_vs_OctMean10_noVal", "2024-08-15", "2024-10-mean10_noVal", ("all", "inland", "coastal"), True),
        ("AugMean9_vs_FebMean13", "2024-08-mean9", "2024-02-mean13", ("all", "inland", "coastal"), True),
        ("AugMean9_vs_OctMean12", "2024-08-mean9", "2024-10-mean12", ("all", "inland", "coastal"), True),
        ("Nov01_vs_Feb14", "2024-11-01", "2024-02-14", ("all", "inland", "coastal"), True),
        ("Nov01_vs_OctMean12", "2024-11-01", "2024-10-mean12", ("all", "inland"), True),
        ("Dec25_vs_Feb14", "2024-12-25", "2024-02-14", ("all", "inland", "coastal"), True),
        ("Dec25_vs_OctMean12", "2024-12-25", "2024-10-mean12", ("all", "inland"), True),
        ("Mar29_vs_Feb14", "2024-03-29", "2024-02-14", ("all", "inland", "coastal"), True),
        ("Mar29_vs_OctMean12", "2024-03-29", "2024-10-mean12", ("all", "inland"), True),
        ("Y2022_Aug15_vs_Feb16", "2022-08-15", "2022-02-16", ("all", "inland", "coastal"), True),
        ("Y2022_AugMean9_vs_FebMean12", "2022-08-mean9", "2022-02-mean12", ("all", "inland"), True),
        ("Y2022_Nov01_vs_Feb16", "2022-11-01", "2022-02-16", ("all", "inland"), True),
        ("Y2023_Aug15_vs_Feb15", "2023-08-15", "2023-02-15", ("all", "inland", "coastal"), True),
        ("Y2023_AugMean9_vs_FebMean12", "2023-08-mean9", "2023-02-mean12", ("all", "inland"), True),
        ("Placebo_OctMean12_vs_FebMean13", "2024-10-mean12", "2024-02-mean13", ("all", "inland", "coastal"), True),
    ]
    c_tasks = [(cn, h, o, s, m, ci) for cn, h, o, ss, ci in MAIN_CONTRASTS for s in ss for m in _G["MODELS"]]
    b2_rows = [r for sub in run_pool(_w_b2_contrast, c_tasks) for r in sub]
    pd.DataFrame(b2_rows)[[
        "contrast", "holiday", "ordinary", "sample", "model", "term", "n", "G",
        "beta_holiday", "beta_ordinary", "delta", "se_crv1", "p_crv1",
        "ci_crv_lo", "ci_crv_hi", "p_wcr", "ci_wcr_lo", "ci_wcr_hi",
    ]].to_csv(OUT / "b2_main_contrasts.csv", index=False)
    log(f"STAGE b1_b2 COMPLETE ({len(b2_rows)} contrast rows).")


def stage_b2_multinight_and_sisters():
    log(f"=== STAGE b2_multinight_and_sisters (workers={N_WORKERS}) ===")
    init_shared_data(load_v1dest=False)
    ns = _G["ns"]
    MULTINIGHT_SPECS = [
        ("2024_Aug9_vs_Feb13_stacked", ns["aug24_sis"], ns["feb24_sis"], ("all", "inland", "coastal"), True),
        ("2024_Aug9_vs_Oct12_stacked", ns["aug24_sis"], ns["oct24_clean12"], ("all", "inland", "coastal"), True),
        ("2022_Aug9_vs_Feb12_stacked", ns["aug22_sis"], ns["feb22_sis"], ("all", "inland"), True),
        ("2023_Aug9_vs_Feb12_stacked", ns["aug23_sis"], ns["feb23_sis"], ("all", "inland"), True),
    ]
    SISTER_GRIDS = [
        ("2024_aug9_x_feb13", ns["aug24_sis"], ns["feb24_sis"], ("all", "inland", "coastal")),
        ("2024_aug9_x_oct12", ns["aug24_sis"], ns["oct24_clean12"], ("all", "inland", "coastal")),
        ("2024_nov3_x_feb13", ns["nov24_sis"], ns["feb24_sis"], ("all", "inland")),
        ("2022_aug9_x_feb12", ns["aug22_sis"], ns["feb22_sis"], ("all", "inland")),
        ("2023_aug9_x_feb12", ns["aug23_sis"], ns["feb23_sis"], ("all", "inland")),
    ]
    needed_d = sorted({d for _, hl, ol, _ in SISTER_GRIDS for d in hl + ol})
    fit_tasks = [(m, s, d) for m in _G["MODELS"] for s in ("all", "inland", "coastal") for d in needed_d]
    FITS = dict(run_pool(_w_fit_single, fit_tasks))
    _G["FITS"] = FITS
    log(f"  Built {len(FITS)} single-night fits in parallel")

    # Stacked multi-night PPML (20 tasks across 14 workers!)
    mn_tasks = [(sn, hl, ol, s, m, ci) for sn, hl, ol, ss, ci in MULTINIGHT_SPECS for s in ss for m in _G["MODELS"]]
    mn_rows = [r for sub in run_pool(_w_b2_multinight, mn_tasks) for r in sub]
    pd.DataFrame(mn_rows).to_csv(OUT / "b2_stacked_multinight.csv", index=False)
    log(f"  b2_stacked_multinight done ({len(mn_rows)} rows)")

    # Pairwise sister grids (parallelized by (grid, sample, model, holiday_night)!)
    sis_tasks = [
        (gname, s, m, h, o_list)
        for gname, h_list, o_list, ss in SISTER_GRIDS
        for s in ss for m in _G["MODELS"] for h in h_list
    ]
    sis_rows = [r for sub in run_pool(_w_b2_sister_slice, sis_tasks) for r in sub]
    df_sis = pd.DataFrame(sis_rows)
    JOURNAL_DATA_BD.mkdir(parents=True, exist_ok=True)
    df_sis.to_parquet(JOURNAL_DATA_BD / "b2_sister_pairs.parquet", index=False)

    sis_sum = []
    for (gname, sample, mname, tnm), g in df_sis.groupby(["grid", "sample", "model", "term"], sort=False):
        vals, pws = g["delta"].to_numpy(), g["p_wcr"].to_numpy()
        med = float(np.median(vals))
        sgn = np.sign(med) if med != 0 else 1.0
        sis_sum.append({
            "grid": gname, "sample": sample, "model": mname, "term": tnm,
            "n_pairs": len(vals), "median_delta": med,
            "q25_delta": float(np.quantile(vals, 0.25)), "q75_delta": float(np.quantile(vals, 0.75)),
            "min_delta": float(vals.min()), "max_delta": float(vals.max()),
            "same_sign_count": int((np.sign(vals) == sgn).sum()),
            "same_sign_share": float((np.sign(vals) == sgn).mean()),
            "sig_wcr_05_count": int(((pws < 0.05) & (np.sign(vals) == sgn)).sum()),
            "sig_wcr_05_share": float(((pws < 0.05) & (np.sign(vals) == sgn)).mean()),
            "median_p_wcr": float(np.median(pws)), "max_p_wcr": float(pws.max()),
        })
    pd.DataFrame(sis_sum).to_csv(OUT / "b2_sister_summary.csv", index=False)
    log(f"STAGE b2_multinight_and_sisters COMPLETE ({len(df_sis)} sister pairs).")


def stage_b3():
    log(f"=== STAGE b3 (workers={N_WORKERS}) ===")
    init_shared_data(load_v1dest=False)
    B3_CONTRASTS = [
        ("Aug15_vs_Feb14", "2024-08-15", "2024-02-14", True),
        ("Aug15_vs_OctMean12", "2024-08-15", "2024-10-mean12", True),
        ("Aug15_vs_FebMean13", "2024-08-15", "2024-02-mean13", True),
        ("AugMean9_vs_FebMean13", "2024-08-mean9", "2024-02-mean13", True),
        ("AugMean9_vs_OctMean12", "2024-08-mean9", "2024-10-mean12", True),
        ("Y2022_Aug15_vs_Feb16", "2022-08-15", "2022-02-16", True),
        ("Y2023_Aug15_vs_Feb15", "2023-08-15", "2023-02-15", True),
        ("Nov01_vs_Feb14", "2024-11-01", "2024-02-14", True),
        ("Dec25_vs_Feb14", "2024-12-25", "2024-02-14", True),
        ("Mar29_vs_Feb14", "2024-03-29", "2024-02-14", True),
    ]
    tasks = []
    for support_name, use_m19, use_scov, strict_int in [
        ("common_19_provs", True, False, False),
        ("common_19_provs_exact_cell_intersect", True, False, True),
        ("common_19_provs_stratum_centroid", True, True, False),
        ("full_47_inland_vs_19_coast", False, False, False),
    ]:
        cl = B3_CONTRASTS if support_name in ("common_19_provs", "full_47_inland_vs_19_coast") else B3_CONTRASTS[:2]
        for cn, h, o, ci in cl:
            for m in _G["MODELS"]:
                tasks.append((support_name, use_m19, use_scov, strict_int, cn, h, o, m, ci))
    b3_rows = [r for sub in run_pool(_w_b3_task, tasks) for r in sub]
    pd.DataFrame(b3_rows).to_csv(OUT / "b3_inland_coast_interaction.csv", index=False)
    log(f"STAGE b3 COMPLETE ({len(b3_rows)} rows).")


def stage_b4():
    log(f"=== STAGE b4 (workers={N_WORKERS}) ===")
    init_shared_data(load_v1dest=False)
    ctx = _G["ctx"]
    MATS, NTOT, ri_all, pj_all = ctx["MATS"], ctx["NTOT"], ctx["ri_all"], ctx["pj_all"]
    lsh05 = lambda col: np.log((MATS[col] + 0.5) / NTOT)[ri_all, pj_all]  # noqa: E731
    _G["AGE_MODELS"] = [
        ("A0_v1_E2a_raw_with_7777",
         np.column_stack([lsh05("n_pre1981_43_64"), lsh05("n_pre1981_ge65"), lsh05("n_post1981_ge65_raw_v1"), lsh05("n_rec2016_lt65")]),
         ["pre1981_43_64", "pre1981_ge65", "post1981_ge65_raw_v1", "recent2016_lt65"], False,
         [(1, 2, "within_65plus__pre1981_minus_post1981_raw")]),
        ("A1_v1_E2a_audit_dated",
         np.column_stack([lsh05("n_pre1981_43_64"), lsh05("n_pre1981_ge65"), lsh05("n_post1981_ge65_dated"), lsh05("n_rec2016_lt65")]),
         ["pre1981_43_64", "pre1981_ge65", "post1981_ge65_dated", "recent2016_lt65"], False,
         [(1, 2, "within_65plus__pre1981_minus_post1981")]),
        ("A2_age_matched_binary_43_64_and_65plus",
         np.column_stack([lsh05("n_pre1981_43_64"), lsh05("n_post1981_43_64_dated"), lsh05("n_pre1981_ge65"), lsh05("n_post1981_ge65_dated"), lsh05("n_post1981_lt43_dated")]),
         ["pre1981_43_64", "post1981_43_64", "pre1981_ge65", "post1981_ge65", "post1981_lt43"], False,
         [(0, 1, "within_43_64__pre1981_minus_post1981"), (2, 3, "within_65plus__pre1981_minus_post1981")]),
        ("A3_age_matched_extremes_43_64_and_65plus",
         np.column_stack([lsh05("n_pre1981_43_64"), lsh05("n_rec2016_43_64"), lsh05("n_pre1981_ge65"), lsh05("n_rec2016_ge65"), lsh05("n_rec2016_lt43")]),
         ["pre1981_43_64", "recent_43_64", "pre1981_ge65", "recent_ge65", "recent_lt43"], False,
         [(0, 1, "within_43_64__pre1981_minus_recent"), (2, 3, "within_65plus__pre1981_minus_recent")]),
        ("A4_age_matched_binary_43_59_and_60plus",
         np.column_stack([lsh05("n_pre1981_43_59"), lsh05("n_post1981_43_59_dated"), lsh05("n_pre1981_ge60"), lsh05("n_post1981_ge60_dated"), lsh05("n_post1981_lt43_dated")]),
         ["pre1981_43_59", "post1981_43_59", "pre1981_ge60", "post1981_ge60", "post1981_lt43"], False,
         [(0, 1, "within_43_59__pre1981_minus_post1981"), (2, 3, "within_60plus__pre1981_minus_post1981")]),
        ("A5_birthsize_age_matched_3bands",
         np.column_stack([lsh05("n_vil2k_lt43_dated"), lsh05("n_urb10k_lt43_dated"), lsh05("n_vil2k_43_64_dated"), lsh05("n_urb10k_43_64_dated"), lsh05("n_vil2k_ge65_dated"), lsh05("n_urb10k_ge65_dated")]),
         ["village_lt43", "urban_lt43", "village_43_64", "urban_43_64", "village_ge65", "urban_ge65"], False,
         [(0, 1, "within_lt43__village_minus_urban"), (2, 3, "within_43_64__village_minus_urban"), (4, 5, "within_65plus__village_minus_urban")]),
        ("A6_cohort_ctrl_city_share65_slopes", _G["X_COH"], ["pre1981", "recent_2016_2023"], True, [(0, 1, "diff_pre1981_minus_recent")]),
        ("A7_birthsize_ctrl_city_share65_slopes", _G["X_BS"], ["born_village_le2k", "born_urban_gt10k"], True, [(0, 1, "diff_village_minus_urban")]),
    ]
    tasks = []
    for idx, (aspec, _, _, _, _) in enumerate(_G["AGE_MODELS"]):
        for cname, o_date in (("Aug15_vs_Feb14", "2024-02-14"), ("Aug15_vs_OctMean12", "2024-10-mean12")):
            if aspec in ("A0_v1_E2a_raw_with_7777", "A4_age_matched_binary_43_59_and_60plus") and cname != "Aug15_vs_Feb14":
                continue
            for sample in ("all", "inland"):
                tasks.append((idx, cname, o_date, sample))
    b4_rows = [r for sub in run_pool(_w_b4_task, tasks) for r in sub]
    pd.DataFrame(b4_rows).to_csv(OUT / "b4_age_matched.csv", index=False)
    log(f"STAGE b4 COMPLETE ({len(b4_rows)} rows).")


def stage_b5():
    log(f"=== STAGE b5 (workers={N_WORKERS}) ===")
    init_shared_data(load_v1dest=False)
    ctx = _G["ctx"]
    MATS, NTOT, ri_all, pj_all = ctx["MATS"], ctx["NTOT"], ctx["ri_all"], ctx["pj_all"]
    hd_675 = V1.FastHDFE2(ctx["c_code"], ctx["p_code"])
    ones_675 = np.ones(len(ri_all))
    _G["RAW_SHARE_STATS"] = {}
    for col in ("n_pre1981", "n_rec2016", "n_mid81_15", "n_vil2k_dated", "n_urb10k_dated", "n_mid2_10k_dated"):
        s_vec = (MATS[col] / NTOT)[ri_all, pj_all]
        s_dm = hd_675.demean(s_vec, ones_675, tol=1e-11)
        _G["RAW_SHARE_STATS"][col] = {
            "sd_raw": float(np.std(s_vec)), "sd_within_fe": float(np.std(s_dm)),
            "mean_raw": float(np.mean(s_vec)), "zero_share": float(np.mean(MATS[col][ri_all, pj_all] == 0)),
        }
    _G["ZERO_RULES"] = [
        ("haldane_0.5_main", lambda s, N: np.log(s + 0.5 / N), np.ones(len(ri_all), bool)),
        ("plus_1.0_log1p", lambda s, N: np.log(s + 1.0 / N), np.ones(len(ri_all), bool)),
        ("plus_0.1_pop_plus1", lambda s, N: np.log(s + 0.1 / N), np.ones(len(ri_all), bool)),
        ("log_share_plus_1e4", lambda s, N: np.log(s + 1e-4), np.ones(len(ri_all), bool)),
        ("log_share_plus_1e5", lambda s, N: np.log(s + 1e-5), np.ones(len(ri_all), bool)),
        ("asinh_1000_share", lambda s, N: np.arcsinh(1000.0 * s), np.ones(len(ri_all), bool)),
        ("R2_112_large_cities_haldane0.5", lambda s, N: np.log(s + 0.5 / N), ctx["cell_112"]),
    ]
    tasks = []
    for cname, o_date in (("Aug15_vs_Feb14", "2024-02-14"), ("Aug15_vs_OctMean12", "2024-10-mean12")):
        for sample in ("all", "inland"):
            for r_idx, (rname, _, _) in enumerate(_G["ZERO_RULES"]):
                mspecs = [
                    ("cohort_joint", ["n_pre1981", "n_rec2016"], ["pre1981", "recent_2016_2023"]),
                    ("birthsize_joint", ["n_vil2k_dated", "n_urb10k_dated"], ["born_village_le2k", "born_urban_gt10k"]),
                ]
                if rname == "haldane_0.5_main":
                    mspecs += [
                        ("cohort_complete_3term", ["n_pre1981", "n_mid81_15", "n_rec2016"], ["pre1981", "mid_1981_2015", "recent_2016_2023"]),
                        ("birthsize_complete_3term", ["n_vil2k_dated", "n_mid2_10k_dated", "n_urb10k_dated"], ["born_village_le2k", "born_mid_2k_10k", "born_urban_gt10k"]),
                    ]
                for mtype, cols, tnames in mspecs:
                    tasks.append((cname, o_date, sample, r_idx, mtype, cols, tnames))
    b5_rows = [r for sub in run_pool(_w_b5_task, tasks) for r in sub]
    pd.DataFrame(b5_rows).to_csv(OUT / "b5_zero_cells_common_scale.csv", index=False)
    log(f"STAGE b5 COMPLETE ({len(b5_rows)} rows).")


def stage_b6():
    log(f"=== STAGE b6 (workers={N_WORKERS}) ===")
    init_shared_data(load_v1dest=False)
    ctx = _G["ctx"]

    # 44-origin-province jackknife
    jk_tasks = list(enumerate(ctx["CL_LEVELS"]))
    jk_rows = [r for sub in run_pool(_w_b6_origin_jk, jk_tasks) for r in sub]
    df_jk = pd.DataFrame(jk_rows)
    df_jk.to_csv(OUT / "b6_jackknife_44provs.csv", index=False)
    jk_sum = []
    for (cname, sample, mname, tnm), g in df_jk.groupby(["contrast", "sample", "model", "term"], sort=False):
        i_min, i_max, i_maxp = g["delta"].idxmin(), g["delta"].idxmax(), g["p_wcr"].idxmax()
        jk_sum.append({
            "contrast": cname, "sample": sample, "model": mname, "term": tnm,
            "min_delta": float(g.loc[i_min, "delta"]), "min_prov": str(g.loc[i_min, "dropped_prov"]),
            "max_delta": float(g.loc[i_max, "delta"]), "max_prov": str(g.loc[i_max, "dropped_prov"]),
            "median_delta": float(g["delta"].median()),
            "max_p_wcr": float(g.loc[i_maxp, "p_wcr"]), "max_p_prov": str(g.loc[i_maxp, "dropped_prov"]),
            "n_sig_05": int((g["p_wcr"] < 0.05).sum()), "n_provs": len(g),
        })
    pd.DataFrame(jk_sum).to_csv(OUT / "b6_jackknife_summary.csv", index=False)
    log(f"  b6 44-origin-prov jackknife done ({len(df_jk)} rows)")

    # Drop Madrid, Barcelona, both (with WCR CIs)
    mb_tasks = [
        (dlabel, dprovs, cname, o_date, sample, mname)
        for dlabel, dprovs in [
            ("drop_Madrid", ["28"]),
            ("drop_Barcelona", ["08"]),
            ("drop_Madrid_and_Barcelona", ["08", "28"]),
        ]
        for cname, o_date in (("Aug15_vs_Feb14", "2024-02-14"), ("Aug15_vs_OctMean12", "2024-10-mean12"))
        for sample in ("all", "inland", "coastal", "common_19_provs_inland_x_coast")
        for mname in _G["MODELS"]
    ]
    mb_rows = [r for sub in run_pool(_w_b6_drop_mb, mb_tasks) for r in sub]
    pd.DataFrame(mb_rows).to_csv(OUT / "b6_drop_mad_bcn.csv", index=False)
    log(f"  b6 drop Madrid/Barcelona done ({len(mb_rows)} rows)")

    # 47-destination-province jackknife
    djk_tasks = list(enumerate(ctx["PROVS"]))
    djk_rows = [r for sub in run_pool(_w_b6_dest_jk, djk_tasks) for r in sub]
    pd.DataFrame(djk_rows).to_csv(OUT / "b6_dest_jackknife_47provs.csv", index=False)
    log(f"STAGE b6 COMPLETE ({len(djk_rows)} dest-jackknife rows).")


def stage_assemble():
    log("=== STAGE assemble -> out/fig_who_pulls_tidy.csv ===")
    df_b2 = pd.read_csv(OUT / "b2_main_contrasts.csv")
    df_mn = pd.read_csv(OUT / "b2_stacked_multinight.csv")
    df_b3 = pd.read_csv(OUT / "b3_inland_coast_interaction.csv")
    df_b4 = pd.read_csv(OUT / "b4_age_matched.csv")
    df_b5 = pd.read_csv(OUT / "b5_zero_cells_common_scale.csv")
    df_mb = pd.read_csv(OUT / "b6_drop_mad_bcn.csv")

    fig_rows = []
    for _, r in df_b2[df_b2["contrast"].isin([
        "Aug15_vs_Feb14", "Aug15_vs_OctMean12", "Aug15_vs_FebMean13",
        "Y2022_Aug15_vs_Feb16", "Y2023_Aug15_vs_Feb15",
        "Nov01_vs_Feb14", "Dec25_vs_Feb14", "Mar29_vs_Feb14",
        "Placebo_OctMean12_vs_FebMean13",
    ])].iterrows():
        fig_rows.append({
            "panel": "A_main_contrasts", "group": r["term"], "model": r["model"],
            "baseline": r["contrast"], "dest_type": r["sample"], "estimand": "delta_elasticity",
            "unit": "log_flow_per_log_share",
            "point": r["delta"], "ci_lo": r["ci_wcr_lo"], "ci_hi": r["ci_wcr_hi"],
            "p_wcr": r["p_wcr"], "n": r["n"], "G": r["G"], "note": f"{r['holiday']}_vs_{r['ordinary']}",
        })
    for _, r in df_mn.iterrows():
        fig_rows.append({
            "panel": "A_stacked_multinight", "group": r["term"], "model": r["model"],
            "baseline": r["spec"], "dest_type": r["sample"], "estimand": "delta_elasticity_stacked",
            "unit": "log_flow_per_log_share",
            "point": r["delta"], "ci_lo": r["ci_wcr_lo"], "ci_hi": r["ci_wcr_hi"],
            "p_wcr": r["p_wcr"], "n": r["n_stacked"], "G": r["G"],
            "note": f"{r['n_hol_nights']}hol_vs_{r['n_ord_nights']}ord",
        })
    for _, r in df_b3[df_b3["support"].isin(["common_19_provs", "full_47_inland_vs_19_coast"])].iterrows():
        fig_rows.append({
            "panel": "B_inland_x_coast", "group": r["term"], "model": r["model"],
            "baseline": r["contrast"], "dest_type": r["support"], "estimand": "gamma_inland_minus_coast",
            "unit": "log_flow_per_log_share",
            "point": r["gamma_inland_minus_coast"], "ci_lo": r["ci_wcr_lo"], "ci_hi": r["ci_wcr_hi"],
            "p_wcr": r["p_wcr"], "n": r["n_4block"], "G": r["G"],
            "note": f"inland={r['delta_inland']:+.4f},coast={r['delta_coast']:+.4f}",
        })
    for _, r in df_b4[df_b4["spec"].isin([
        "A1_v1_E2a_audit_dated", "A2_age_matched_binary_43_64_and_65plus",
        "A3_age_matched_extremes_43_64_and_65plus", "A5_birthsize_age_matched_3bands",
    ])].iterrows():
        fig_rows.append({
            "panel": "C_age_matched", "group": r["term"], "model": r["spec"],
            "baseline": r["contrast"], "dest_type": r["sample"],
            "estimand": "within_band_diff" if r["is_within_band_diff"] else "delta_elasticity",
            "unit": "log_flow_per_log_share",
            "point": r["delta"], "ci_lo": r["ci_wcr_lo"], "ci_hi": r["ci_wcr_hi"],
            "p_wcr": r["p_wcr"], "n": r["n"], "G": r["G"], "note": "",
        })
    for _, r in df_b5[df_b5["model"].isin(["cohort_joint", "birthsize_joint"])].iterrows():
        fig_rows.append({
            "panel": "D_zero_cells_common_scale", "group": r["term"], "model": r["rule"],
            "baseline": r["contrast"], "dest_type": r["sample"], "estimand": "pct_flow_1sd_within_raw_share_fw",
            "unit": "percent_flow_change_aug_vs_ord",
            "point": r["pct_flow_1sd_rawshare_fw"],
            "ci_lo": r["pct_flow_1sd_rawshare_fw_ci_lo"],
            "ci_hi": r["pct_flow_1sd_rawshare_fw_ci_hi"],
            "p_wcr": r["p_wcr"], "n": r["n_aug"] + r["n_ord"], "G": r["G"],
            "note": f"raw_delta={r['delta_raw']:+.4f},pct_1sd_g={r['pct_flow_1sd_g_within']:+.2f}%",
        })
    for _, r in df_mb.iterrows():
        fig_rows.append({
            "panel": "E_drop_madrid_barcelona", "group": r["term"], "model": r["model"],
            "baseline": r["contrast"], "dest_type": r["sample"], "estimand": r["estimand"],
            "unit": "log_flow_per_log_share",
            "point": r["delta"], "ci_lo": r["ci_wcr_lo"], "ci_hi": r["ci_wcr_hi"],
            "p_wcr": r["p_wcr"], "n": r["n"], "G": r["G"], "note": r["drop_spec"],
        })
    pd.DataFrame(fig_rows).to_csv(OUT / "fig_who_pulls_tidy.csv", index=False)
    log(f"Saved out/fig_who_pulls_tidy.csv ({len(fig_rows)} rows). ALL ANALYSIS B STAGES COMPLETE.")


STAGES = {
    "b1_b2": stage_b1_b2,
    "b2_multinight_and_sisters": stage_b2_multinight_and_sisters,
    "b3": stage_b3,
    "b4": stage_b4,
    "b5": stage_b5,
    "b6": stage_b6,
    "assemble": stage_assemble,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=["all"] + list(STAGES.keys()))
    args = ap.parse_args()
    if args.stage == "all":
        for nm in ("b1_b2", "b2_multinight_and_sisters", "b3", "b4", "b5", "b6", "assemble"):
            STAGES[nm]()
    else:
        STAGES[args.stage]()


if __name__ == "__main__":
    main()

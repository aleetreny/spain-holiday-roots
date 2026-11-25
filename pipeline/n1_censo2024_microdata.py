"""n1_censo2024_microdata: Test Censo Anual 2024 10% microdata (4.86M residents, 780 municipalities >10k).

Fixes three limitations of R2's table-68545 neighbouring-city design simultaneously:
  1. Expands from 112 identifying cities in 22 provinces -> ~650 mainland municipalities (>10k) in ~45-47 provinces,
     resolving Prof. Zhou's M2 (few clusters, Barcelona dependence).
  2. Separates PRE-1981 rural-exodus arrivals (VARANOP == '0001') from RECENT migrants (VARANOP in 2016..2023)
     AT THE NEIGHBOURING-CITY LEVEL (with city FE + origin_prov x dest_prov FE), resolving G5.
  3. Separates VILLAGE-BORN (CMUN_NAC == '991', <=2,000 pop) from URBAN-BORN (>10k city) natives of province p.
  4. Tests true >10k municipality x >10k municipality roots on INE municipal OD.

Run: .venv/bin/python -I pipeline/n1_censo2024_microdata.py 2>&1 | tee results/n1_censo2024.log
Outputs:
  results/n1_censo2024_expanded_city_prov.csv   Neighbouring-city PPML on ~650 municipalities (~46 clusters) + WCR + Barcelona jackknife
  results/n1_censo2024_cohorts.csv              Pre-1981 rural exodus vs recent (2016-23) & village-born vs urban-born
  results/n1_censo2024_mun10k_od.csv            True >10k municipality x >10k municipality roots on INE OD (with cp_fe + dest FE)
"""
from __future__ import annotations

import gzip
import sys
import time
import zlib
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pyfixest as pf
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import b00_common as B0  # noqa: E402
from n1_new_data import HDFE, ppml_hdfe, ppml_crv1_and_wcr, parse_es_num, seed_for, PROV_NAME_TO_CODE  # noqa: E402

T0 = time.time()
EXT = C.ROOT / "_local" / "external"
OUT = C.RESULTS
NET = C.RESULTS / "b01_networks"
EXCL = set(C.EXCLUDED_PROVINCES)

CEN24_PQ = EXT / "ine_censo_microdatos_2024" / "md_CensoAnual_muestra_2024.parquet"
ADRH_CSV = EXT / "ine_adrh" / "adrh_municipios.csv"
INE_OD_DIR = EXT / "ine_turismo_moviles_od"
INE_TM_RAW = EXT / "ine_turismo_moviles" / "raw"


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only-cdiff", action="store_true")
    args = ap.parse_args()
    log("Loading Censo Anual 2024 10% microdata (4.86M rows)...")
    cols = ["CPRO", "CMUN", "CPRO_NAC", "CMUN_NAC", "VARANOP"]
    df = pd.read_parquet(CEN24_PQ, columns=cols)
    for c in cols:
        df[c] = df[c].astype(str).str.strip()

    mitma_nodes = B0.load_nodes()
    NID = mitma_nodes["id"].values
    NID_SET = set(NID)
    PROVS = np.array(sorted(mitma_nodes["prov"].unique()))
    P2I = {p: k for k, p in enumerate(PROVS)}
    NPI = np.array([P2I[p] for p in mitma_nodes["prov"].values])
    NCOAST = mitma_nodes["coastal_true"].values.astype(bool)
    SINGLE_ZONES = set(mitma_nodes.loc[mitma_nodes["n_municipalities"] == 1, "id"])

    # Exact residence municipality (>10k): not in 991,992,993,999,666
    df["m_res"] = df["CPRO"].str.zfill(2) + df["CMUN"].str.zfill(3)
    df["m_nac"] = df["CPRO_NAC"].str.zfill(2) + df["CMUN_NAC"].str.zfill(3)
    bad_cmun = {"991", "992", "993", "999", "666", "966", "777", "888", ""}
    ok_res = (~df["CPRO"].isin(EXCL)) & (~df["CMUN"].isin(bad_cmun))
    sub = df[ok_res].copy()

    # Keep mainland >10k municipalities that are ALSO exact single-municipality MITMA zones AND have ADRH 2024 income
    A = pd.read_csv(ADRH_CSV, dtype={"ine_code": str})
    inc24 = A[(A["year"] == 2024) & (A["indicator"] == "Renta neta media por persona")].set_index("ine_code")["value"].dropna()

    cand_muns = sorted(set(sub["m_res"].unique()) & SINGLE_ZONES & set(inc24.index))
    pcount = Counter(m[:2] for m in cand_muns)
    id_muns = np.array(sorted(m for m in cand_muns if pcount[m[:2]] >= 2))
    id_provs = sorted({m[:2] for m in id_muns})
    log(
        f"Expanded >10k single-zone municipalities with ADRH income: {len(cand_muns)} total, "
        f"{len(id_muns)} identifying in {len(id_provs)} provinces (vs 112 in 22 in R2)!"
    )

    # Filter microdata to identifying municipalities
    sub = sub[sub["m_res"].isin(set(id_muns))].copy()
    M2R = {m: i for i, m in enumerate(id_muns)}
    MPROV = np.array([m[:2] for m in id_muns])
    M_NODE = np.array([np.where(NID == m)[0][0] for m in id_muns])
    NODE2ROW = -np.ones(len(NID), dtype=int)
    NODE2ROW[M_NODE] = np.arange(len(id_muns))

    # Build municipality x birth-province matrices (multiply 10% sample counts by 10 for approximate persons;
    # scale cancels in shares/logs anyway)
    tot_m = sub.groupby("m_res").size().reindex(id_muns).values.astype(float)
    sp = sub[sub["CPRO_NAC"].isin(PROVS)]

    def make_mat(df_s):
        g = df_s.groupby(["m_res", "CPRO_NAC"]).size()
        M = np.zeros((len(id_muns), 47), dtype=float)
        for (m, p), cnt in g.items():
            M[M2R[m], P2I[p]] = float(cnt)
        return M

    B_ALL = make_mat(sp)
    B_PRE81 = make_mat(sp[sp["VARANOP"] == "0001"])                                    # Arrived in province BEFORE 1981 (rural exodus!)
    B_MID = make_mat(sp[sp["VARANOP"].isin(["0002", "0003", "0004", "0005"])])         # 1981-2015
    B_REC16 = make_mat(sp[sp["VARANOP"].isin([str(y) for y in range(2016, 2025)])])    # Recent 2016-2024
    B_VIL = make_mat(sp[sp["CMUN_NAC"] == "991"])                                      # Born in village <= 2,000 pop!
    B_SM = make_mat(sp[sp["CMUN_NAC"].isin(["991", "992"])])                           # Born in small mun <= 5,000 pop
    B_URB = make_mat(sp[~sp["CMUN_NAC"].isin(bad_cmun)])                               # Born in >10k city!

    OWN = MPROV[:, None] == PROVS[None, :]
    MASK = ~OWN
    ri_all, pj_all = np.where(MASK)

    # Distance & income controls (exact R2 structure)
    pw = mitma_nodes["population"].values.astype(float)
    PCX = np.array([np.average(mitma_nodes["x"].values[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    PCY = np.array([np.average(mitma_nodes["y"].values[NPI == k], weights=pw[NPI == k]) for k in range(47)])
    MX, MY = mitma_nodes["x"].values[M_NODE], mitma_nodes["y"].values[M_NODE]
    KM = np.hypot(MX[:, None] - PCX[None, :], MY[:, None] - PCY[None, :]) / 1000.0
    lkm = np.log(np.maximum(KM[ri_all, pj_all], 1.0))
    lkm_c = lkm - lkm.mean()

    LINC = np.log(inc24.reindex(id_muns).values.astype(float))
    LINC = LINC - LINC.mean()
    inc_cols = [np.where(pj_all == k, LINC[ri_all], 0.0) for k in range(1, 47)]
    X_ctrl = np.column_stack([lkm_c, lkm_c ** 2] + inc_cols)

    # Log shares (with 0.5-person pseudocount in 10% sample -> 0.5/tot_m for zeros)
    def lshare(M):
        return np.log((M[ri_all, pj_all] + 0.5) / tot_m[ri_all])

    x_all = lshare(B_ALL)
    x_pre81 = lshare(B_PRE81)
    x_rec16 = lshare(B_REC16)
    x_vil = lshare(B_VIL)
    x_urb = lshare(B_URB)

    oprov_i = np.array([P2I[p] for p in MPROV[ri_all]])
    pair_fe = oprov_i * 47 + pj_all
    fe_codes = [ri_all, pair_fe]
    cl_oprov = oprov_i

    # Check zero rates & within-FE correlations
    hd0 = HDFE([np.unique(c, return_inverse=True)[1] for c in fe_codes])
    w1 = np.ones(len(x_all))
    d_pre = hd0.demean(x_pre81, w1)
    d_rec = hd0.demean(x_rec16, w1)
    d_vil = hd0.demean(x_vil, w1)
    d_urb = hd0.demean(x_urb, w1)
    log(
        f"  Within-FE corr(pre81, rec16)={np.corrcoef(d_pre, d_rec)[0,1]:.3f}; "
        f"corr(village<=2k, urban>10k)={np.corrcoef(d_vil, d_urb)[0,1]:.3f}; "
        f"zero rates: all={(B_ALL[MASK]==0).mean():.1%}, pre81={(B_PRE81[MASK]==0).mean():.1%}, rec16={(B_REC16[MASK]==0).mean():.1%}"
    )

    # Load MITMA main/appendix nights (+ Oct mean baseline) and INE 53000 monthly
    datasets = {}
    for dlabel, ddate in [
        ("MITMA_Feb14", "2024-02-14"),
        ("MITMA_Aug15", "2024-08-15"),
        ("MITMA_Nov01", "2024-11-01"),
        ("MITMA_Dec25", "2024-12-25"),
        ("MITMA_Mar29", "2024-03-29"),
    ]:
        f = NET / f"edges_{ddate}__full.parquet"
        e = pd.read_parquet(f) if f.exists() else B0.build_edges(ddate, mitma_nodes)[0]
        i, j, w = e["i"].values, e["j"].values, e["w"].values
        r = NODE2ROW[i]
        pj = NPI[j]
        for kind, keep in [("all", np.ones(len(j), bool)), ("inland", ~NCOAST[j])]:
            M = np.zeros((len(id_muns), 47))
            sel = (r >= 0) & keep
            np.add.at(M, (r[sel], pj[sel]), w[sel])
            datasets[f"{dlabel}_{kind}"] = M[ri_all, pj_all]

    # INE 53000 monthly (Feb, Mar, Aug, Nov, Dec 2024)
    recs = []
    target_per = {"2024M02", "2024M03", "2024M08", "2024M11", "2024M12"}
    with gzip.open(INE_TM_RAW / "turismo_53000_munorig_x_provdest.csv.gz", "rt", encoding="utf-8") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(";")
            if len(p) < 7 or p[5].strip() not in target_per:
                continue
            oc = p[0][:5]
            if oc not in M2R:
                continue
            dp = PROV_NAME_TO_CODE.get(p[3].strip())
            if dp not in P2I:
                continue
            v = parse_es_num(p[6])
            if v is not None:
                recs.append((p[5].strip(), M2R[oc], P2I[dp], v))
    d53 = pd.DataFrame(recs, columns=["per", "ri", "pj", "val"])
    assert len(d53) > 30000, f"53000 parse too small: {len(d53)}"
    for per in sorted(target_per):
        M = np.zeros((len(id_muns), 47))
        sub53 = d53[d53["per"] == per]
        M[sub53["ri"].values, sub53["pj"].values] = sub53["val"].values
        datasets[f"INE53000_{per}"] = M[ri_all, pj_all]

    # --- 1. Expanded neighbouring-city design (~650 muns, ~46 clusters) + Barcelona jackknife ---
    exp_rows = []
    if args.only_cdiff:
        cdiff_rows = []
        for h_key, o_key, lbl in [
            ("MITMA_Aug15_all", "MITMA_Feb14_all", "MITMA_Aug15_vs_Feb14_all"),
            ("MITMA_Aug15_inland", "MITMA_Feb14_inland", "MITMA_Aug15_vs_Feb14_inland"),
            ("MITMA_Dec25_all", "MITMA_Feb14_all", "MITMA_Dec25_vs_Feb14_all"),
            ("INE53000_2024M08", "INE53000_2024M02", "INE53000_Aug_vs_Feb"),
        ]:
            y_st = np.concatenate([datasets[o_key], datasets[h_key]])
            hol = np.concatenate([np.zeros(len(ri_all)), np.ones(len(ri_all))])
            fe_st = [np.concatenate([ri_all, ri_all + 1000]), np.concatenate([pair_fe, pair_fe + 5000])]
            cl_st = np.concatenate([cl_oprov, cl_oprov])
            Xc_st = np.vstack([np.hstack([X_ctrl, np.zeros_like(X_ctrl)]), np.hstack([np.zeros_like(X_ctrl), X_ctrl])])
            for grp_name, xa, xb, la, lb in [
                ("cohort_joint", x_pre81, x_rec16, "pre1981_exodus", "recent_2016_2024"),
                ("birthsize_joint", x_vil, x_urb, "born_village_le2k", "born_urban_gt10k"),
            ]:
                xa_st = np.concatenate([xa, xa])
                xb_st = np.concatenate([xb, xb])
                X_st = np.column_stack([xa_st * hol, xb_st * hol, xa_st, xb_st, Xc_st])
                for fidx, tlbl in [(0, la), (1, lb)]:
                    res = ppml_crv1_and_wcr(y_st, X_st, fe_st, cl_st, focus_idx=fidx, n_boot=4999, seed=seed_for("cdiff", lbl, grp_name, tlbl), do_ci=True)
                    cdiff_rows.append({
                        "contrast": lbl, "model": grp_name, "term": f"{tlbl}_x_holiday",
                        "n": res["n"], "G": res["G"],
                        "delta": res["beta"], "se_crv1": res["se_crv1"], "p_crv1": res["p_crv1"],
                        "p_wcr": res["p_wcr"], "ci_wcr_lo": res["ci_wcr_lo"], "ci_wcr_hi": res["ci_wcr_hi"],
                    })
                    log(f"  cdiff {lbl} {tlbl}: delta={res['beta']:+.4f} [{res['ci_wcr_lo']:+.4f}, {res['ci_wcr_hi']:+.4f}] p_wcr={res['p_wcr']:.4f}")
        pd.DataFrame(cdiff_rows).to_csv(OUT / "n1_censo2024_cohort_diffs.csv", index=False)
        log("Saved n1_censo2024_cohort_diffs.csv")
        return
    for dname in [
        "MITMA_Feb14_all", "MITMA_Aug15_all", "MITMA_Nov01_all", "MITMA_Dec25_all", "MITMA_Mar29_all",
        "MITMA_Feb14_inland", "MITMA_Aug15_inland", "MITMA_Dec25_inland",
        "INE53000_2024M02", "INE53000_2024M03", "INE53000_2024M08", "INE53000_2024M11", "INE53000_2024M12",
    ]:
        do_ci = dname in ("MITMA_Feb14_all", "MITMA_Aug15_all", "MITMA_Aug15_inland", "MITMA_Dec25_all", "INE53000_2024M02", "INE53000_2024M08", "INE53000_2024M12")
        res = ppml_crv1_and_wcr(
            datasets[dname], np.column_stack([x_all, X_ctrl]), fe_codes, cl_oprov,
            focus_idx=0, n_boot=4999, seed=seed_for("exp", dname), do_ci=do_ci,
        )
        exp_rows.append({
            "dataset": dname, "sample": "all_id_provs", "n": res["n"], "G": res["G"],
            "beta": res["beta"], "se_crv1": res["se_crv1"], "p_crv1": res["p_crv1"],
            "ci_crv_lo": res["ci_crv_lo"], "ci_crv_hi": res["ci_crv_hi"],
            "p_wcr": res["p_wcr"], "ci_wcr_lo": res["ci_wcr_lo"], "ci_wcr_hi": res["ci_wcr_hi"],
        })
        log(f"  exp done {dname}: beta={res['beta']:.4f} p_wcr={res['p_wcr']:.4f}")

    # Stacked differences (with & without Barcelona!)
    bcn_i = P2I["08"]
    for h_key, o_key, lbl in [
        ("MITMA_Aug15_all", "MITMA_Feb14_all", "MITMA_Aug15_vs_Feb14_all"),
        ("MITMA_Aug15_inland", "MITMA_Feb14_inland", "MITMA_Aug15_vs_Feb14_inland"),
        ("MITMA_Dec25_all", "MITMA_Feb14_all", "MITMA_Dec25_vs_Feb14_all"),
        ("MITMA_Nov01_all", "MITMA_Feb14_all", "MITMA_Nov01_vs_Feb14_all"),
        ("INE53000_2024M08", "INE53000_2024M02", "INE53000_Aug_vs_Feb"),
        ("INE53000_2024M12", "INE53000_2024M02", "INE53000_Dec_vs_Feb"),
    ]:
        for s_lbl, keep_m in [("all_id_provs", np.ones(len(ri_all), bool)), ("drop_Barcelona", oprov_i != bcn_i)]:
            if s_lbl == "drop_Barcelona" and lbl not in ("MITMA_Aug15_vs_Feb14_all", "MITMA_Aug15_vs_Feb14_inland", "INE53000_Aug_vs_Feb"):
                continue
            y0, y1 = datasets[o_key][keep_m], datasets[h_key][keep_m]
            y_st = np.concatenate([y0, y1])
            hol = np.concatenate([np.zeros(len(y0)), np.ones(len(y1))])
            ri_s, pf_s, cl_s = ri_all[keep_m], pair_fe[keep_m], cl_oprov[keep_m]
            fe_st = [np.concatenate([ri_s, ri_s + 1000]), np.concatenate([pf_s, pf_s + 5000])]
            cl_st = np.concatenate([cl_s, cl_s])
            Xc_s = X_ctrl[keep_m]
            Xc_st = np.vstack([np.hstack([Xc_s, np.zeros_like(Xc_s)]), np.hstack([np.zeros_like(Xc_s), Xc_s])])
            x_s = np.concatenate([x_all[keep_m], x_all[keep_m]])
            X_st = np.column_stack([x_s * hol, x_s, Xc_st])
            res = ppml_crv1_and_wcr(y_st, X_st, fe_st, cl_st, focus_idx=0, n_boot=4999, seed=seed_for("expdiff", lbl, s_lbl), do_ci=True)
            exp_rows.append({
                "dataset": lbl, "sample": s_lbl, "n": res["n"], "G": res["G"],
                "beta": res["beta"], "se_crv1": res["se_crv1"], "p_crv1": res["p_crv1"],
                "ci_crv_lo": res["ci_crv_lo"], "ci_crv_hi": res["ci_crv_hi"],
                "p_wcr": res["p_wcr"], "ci_wcr_lo": res["ci_wcr_lo"], "ci_wcr_hi": res["ci_wcr_hi"],
            })
            log(f"  expdiff done {lbl} ({s_lbl}): delta={res['beta']:.4f} p_wcr={res['p_wcr']:.4f}")
    exp_df = pd.DataFrame(exp_rows)
    exp_df.to_csv(OUT / "n1_censo2024_expanded_city_prov.csv", index=False)
    log("Saved n1_censo2024_expanded_city_prov.csv:\n" + exp_df.to_string(index=False))

    # --- 2. Cohort decomposition (Pre-1981 exodus vs Recent 2016-24) & Village-born vs Urban-born ---
    coh_rows = []
    for dname in [
        "MITMA_Feb14_all", "MITMA_Aug15_all", "MITMA_Aug15_inland", "MITMA_Nov01_all", "MITMA_Dec25_all",
        "INE53000_2024M02", "INE53000_2024M08", "INE53000_2024M12",
    ]:
        yv = datasets[dname]
        for sname, Xm, fidx, tlbl in [
            ("cohort_joint__pre1981", np.column_stack([x_pre81, x_rec16, X_ctrl]), 0, "pre1981_exodus"),
            ("cohort_joint__rec2016", np.column_stack([x_pre81, x_rec16, X_ctrl]), 1, "recent_2016_2024"),
            ("birthsize_joint__village_le2k", np.column_stack([x_vil, x_urb, X_ctrl]), 0, "born_village_le2k"),
            ("birthsize_joint__urban_gt10k", np.column_stack([x_vil, x_urb, X_ctrl]), 1, "born_urban_gt10k"),
        ]:
            do_ci = dname in ("MITMA_Feb14_all", "MITMA_Aug15_all", "MITMA_Aug15_inland", "MITMA_Dec25_all", "INE53000_2024M08")
            res = ppml_crv1_and_wcr(yv, Xm, fe_codes, cl_oprov, focus_idx=fidx, n_boot=4999, seed=seed_for("coh", dname, sname), do_ci=do_ci)
            coh_rows.append({
                "dataset": dname, "spec": sname, "term": tlbl, "n": res["n"], "G": res["G"],
                "beta": res["beta"], "se_crv1": res["se_crv1"], "p_crv1": res["p_crv1"],
                "ci_crv_lo": res["ci_crv_lo"], "ci_crv_hi": res["ci_crv_hi"],
                "p_wcr": res["p_wcr"], "ci_wcr_lo": res["ci_wcr_lo"], "ci_wcr_hi": res["ci_wcr_hi"],
            })
        log(f"  cohorts done {dname}")
    coh_df = pd.DataFrame(coh_rows)
    coh_df.to_csv(OUT / "n1_censo2024_cohorts.csv", index=False)
    log("Saved n1_censo2024_cohorts.csv:\n" + coh_df[["dataset", "term", "beta", "ci_wcr_lo", "ci_wcr_hi", "p_wcr"]].to_string(index=False))

    # --- 2b. Formal stacked difference tests (holiday - ordinary) for pre-1981 vs recent AND village<=2k vs urban>10k ---
    cdiff_rows = []
    for h_key, o_key, lbl in [
        ("MITMA_Aug15_all", "MITMA_Feb14_all", "MITMA_Aug15_vs_Feb14_all"),
        ("MITMA_Aug15_inland", "MITMA_Feb14_inland", "MITMA_Aug15_vs_Feb14_inland"),
        ("MITMA_Dec25_all", "MITMA_Feb14_all", "MITMA_Dec25_vs_Feb14_all"),
        ("INE53000_2024M08", "INE53000_2024M02", "INE53000_Aug_vs_Feb"),
    ]:
        y_st = np.concatenate([datasets[o_key], datasets[h_key]])
        hol = np.concatenate([np.zeros(len(ri_all)), np.ones(len(ri_all))])
        fe_st = [np.concatenate([ri_all, ri_all + 1000]), np.concatenate([pair_fe, pair_fe + 5000])]
        cl_st = np.concatenate([cl_oprov, cl_oprov])
        Xc_st = np.vstack([np.hstack([X_ctrl, np.zeros_like(X_ctrl)]), np.hstack([np.zeros_like(X_ctrl), X_ctrl])])
        for grp_name, xa, xb, la, lb in [
            ("cohort_joint", x_pre81, x_rec16, "pre1981_exodus", "recent_2016_2024"),
            ("birthsize_joint", x_vil, x_urb, "born_village_le2k", "born_urban_gt10k"),
        ]:
            xa_st = np.concatenate([xa, xa])
            xb_st = np.concatenate([xb, xb])
            X_st = np.column_stack([xa_st * hol, xb_st * hol, xa_st, xb_st, Xc_st])
            for fidx, tlbl in [(0, la), (1, lb)]:
                res = ppml_crv1_and_wcr(y_st, X_st, fe_st, cl_st, focus_idx=fidx, n_boot=4999, seed=seed_for("cdiff", lbl, grp_name, tlbl), do_ci=True)
                cdiff_rows.append({
                    "contrast": lbl, "model": grp_name, "term": f"{tlbl}_x_holiday",
                    "n": res["n"], "G": res["G"],
                    "delta": res["beta"], "se_crv1": res["se_crv1"], "p_crv1": res["p_crv1"],
                    "p_wcr": res["p_wcr"], "ci_wcr_lo": res["ci_wcr_lo"], "ci_wcr_hi": res["ci_wcr_hi"],
                })
                log(f"  cdiff {lbl} {tlbl}: delta={res['beta']:+.4f} [{res['ci_wcr_lo']:+.4f}, {res['ci_wcr_hi']:+.4f}] p_wcr={res['p_wcr']:.4f}")
    pd.DataFrame(cdiff_rows).to_csv(OUT / "n1_censo2024_cohort_diffs.csv", index=False)
    log("Saved n1_censo2024_cohort_diffs.csv")

    # --- 3. True >10k municipality x >10k municipality roots (Censo 2024 10%) on INE OD ---
    # Check exact >10k x >10k interprovincial pairs (~700 x ~700 = ~480k cells) with (orig x dest_prov) FE + dest FE!
    m10k_ml = sorted(set(sub["m_res"].unique()) & set(mitma_nodes["id"]))
    m10k_set = set(m10k_ml)
    sub_b = sub[(~sub["CPRO_NAC"].isin(EXCL)) & (~sub["CMUN_NAC"].isin(bad_cmun)) & (sub["CPRO"] != sub["CPRO_NAC"]) & sub["m_nac"].isin(m10k_set)]
    b10 = sub_b.groupby(["m_res", "m_nac"]).size().rename("born10").reset_index()
    b10.columns = ["orig", "dest", "born10"]

    od24 = pd.read_parquet(INE_OD_DIR / "ine_od_mun_2024.parquet")
    od10 = od24[od24["orig"].isin(m10k_set) & od24["dest"].isin(m10k_set)].copy()
    # Active origins/dests among m10k
    o_arr = np.array(m10k_ml)
    o_prov = np.array([x[:2] for x in o_arr])
    nd_idx = mitma_nodes.set_index("id")
    ox = nd_idx.loc[o_arr, "x"].values
    oy = nd_idx.loc[o_arr, "y"].values
    ii, jj = np.where(o_prov[:, None] != o_prov[None, :])
    dist10 = np.hypot(ox[ii] - ox[jj], oy[ii] - oy[jj]) / 1000.0
    sk = pd.DataFrame({
        "orig": o_arr[ii], "dest": o_arr[jj], "dprov": o_prov[jj],
        "cp_fe": pd.Series(o_arr[ii]) + "_" + pd.Series(o_prov[jj]),
        "lkm": np.log(np.maximum(dist10, 1.0)),
    })
    sk = sk.merge(b10, on=["orig", "dest"], how="left").fillna({"born10": 0.0})
    sk["lborn"] = np.log1p(sk["born10"] * 10.0)
    m10_rows = []
    for mm in ["2024-02", "2024-03", "2024-08", "2024-10", "2024-11", "2024-12"]:
        f = od10[od10["mes"] == mm].set_index(["orig", "dest"])["turistas"]
        sk["y"] = sk.set_index(["orig", "dest"]).index.map(f).fillna(0.0).values
        m = pf.fepois("y ~ lborn + lkm | cp_fe + dest", data=sk, vcov={"CRV1": "dprov"}, iwls_tol=1e-7, iwls_maxiter=40)
        r = m.tidy().loc["lborn"]
        m10_rows.append({
            "mes": mm, "n": int(m._N), "n_muns": len(m10k_ml), "beta": float(r["Estimate"]),
            "se": float(r["Std. Error"]), "t": float(r["t value"]), "p": float(r["Pr(>|t|)"]),
            "ci_lo": float(r["2.5%"]), "ci_hi": float(r["97.5%"]),
        })
    m10_df = pd.DataFrame(m10_rows)
    m10_df.to_csv(OUT / "n1_censo2024_mun10k_od.csv", index=False)
    log("Saved n1_censo2024_mun10k_od.csv:\n" + m10_df.to_string(index=False))
    log("DONE.")


if __name__ == "__main__":
    main()

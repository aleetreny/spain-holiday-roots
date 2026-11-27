"""Companion diagnostic for Analysis D: inland-only controlled arrivals/departures
decomposition (showing why all-zone linear arrivals/res is pulled down by coastal resorts)
and multi-year / 306-night summary statistics for REPORT_BD.md.

Run:
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    /Users/alejandrotreny/Proyectos/spain-holiday-roots/.venv/bin/python -I \
    journal/analysis/D_villages/d_inland_decomp.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))

import run_d_villages as RD  # noqa: E402


def main():
    nd = RD.load_fixed_nodes_and_housing()
    nd_el = nd[nd["eligible"]].copy().reset_index(drop=True)
    panel = pd.read_parquet(RD.JOURNAL_DATA_BD / "d_zone_daily_panel.parquet")

    # Build profiles for key nights
    el_ids = nd_el["id"].tolist()
    pop_el = nd_el["population"].to_numpy(dtype=float)
    all_dates = sorted(panel["date"].unique())
    d2row = {d: i for i, d in enumerate(all_dates)}
    ps = panel.sort_values(["date", "id"])
    N_D, N_Z = len(all_dates), len(el_ids)
    SIN = ps["s_in"].to_numpy(float).reshape(N_D, N_Z)
    SOUT = ps["s_out"].to_numpy(float).reshape(N_D, N_Z)

    oct12 = [d for d in all_dates if d.startswith("2024-10") and pd.Timestamp(d).day_of_week in (1, 2, 3) and d not in RD.DANA_DATES]
    PROF = {}
    for k, dl in [
        ("2024-02-14", ["2024-02-14"]), ("2024-10-mean12", oct12),
        ("2024-08-15", ["2024-08-15"]), ("2024-11-01", ["2024-11-01"]),
        ("2024-12-25", ["2024-12-25"]), ("2024-03-29", ["2024-03-29"]),
        ("2022-02-16", ["2022-02-16"]), ("2022-08-15", ["2022-08-15"]),
        ("2023-02-15", ["2023-02-15"]), ("2023-08-15", ["2023-08-15"]),
    ]:
        ii = [d2row[d] for d in dl]
        si, so = SIN[ii].mean(0), SOUT[ii].mean(0)
        PROF[k] = {"s_in": si, "s_out": so, "arr": si / pop_el, "dep": so / pop_el,
                   "net": (si - so) / pop_el, "imp": (si > so).astype(float)}

    inland_m = ~nd_el["coastal_true"].to_numpy(bool)
    single_m = (nd_el["n_municipalities"] == 1).to_numpy(bool)
    SAMPLES = {"inland": inland_m, "inland_single_mun": inland_m & single_m}
    CTRLS_INL = ["lpop", "lnmun", "larea", "ldcity", "ldcoast"]

    rows = []
    for sname, smask in SAMPLES.items():
        d_sub = nd_el.loc[smask].reset_index(drop=True)
        x_loss = d_sub["loss_pct"].to_numpy(float)
        Z = d_sub[CTRLS_INL].to_numpy(float)
        provs_s = d_sub["prov"].to_numpy(str)
        b_idx = RD.make_zone_bootstrap_indices(provs_s, RD.N_BOOT_ZONE, RD.rng_for("d_inl_zboot", sname))
        idx_p, picks = RD.make_prov_bootstrap_draws(provs_s, RD.N_BOOT_PROV, RD.rng_for("d_inl_pboot", sname))

        col_meta, Y_cols = [], []
        for kn in ("2024-02-14", "2024-10-mean12", "2024-08-15", "2024-12-25"):
            p = PROF[kn]
            for oname, yv in [("arrivals_per_res", p["arr"][smask]), ("departures_per_res", p["dep"][smask]),
                              ("net_per_res", p["net"][smask]), ("net_importer", p["imp"][smask])]:
                col_meta.append(("level", kn, oname))
                Y_cols.append(yv)
        for pname, h, o in [
            ("Aug15_minus_Feb14", "2024-08-15", "2024-02-14"),
            ("Aug15_minus_OctMean12", "2024-08-15", "2024-10-mean12"),
            ("Dec25_minus_Feb14", "2024-12-25", "2024-02-14"),
            ("Dec25_minus_OctMean12", "2024-12-25", "2024-10-mean12"),
            ("Y2022_Aug15_minus_Feb16", "2022-08-15", "2022-02-16"),
            ("Y2023_Aug15_minus_Feb15", "2023-08-15", "2023-02-15"),
        ]:
            ph, po = PROF[h], PROF[o]
            for oname, yv in [
                ("arrivals_per_res_diff", (ph["arr"] - po["arr"])[smask]),
                ("departures_per_res_diff", (ph["dep"] - po["dep"])[smask]),
                ("net_per_res_diff", (ph["net"] - po["net"])[smask]),
                ("log_s_in_ratio", np.log((ph["s_in"][smask] + 1.0) / (po["s_in"][smask] + 1.0))),
                ("net_importer_diff", (ph["imp"] - po["imp"])[smask]),
            ]:
                col_meta.append(("change", pname, oname))
                Y_cols.append(yv)

        Y_mat = np.column_stack(Y_cols)
        b_full, b_zone = RD.fast_ols_fe_zone_boot(Y_mat, x_loss, Z, provs_s, b_idx)
        b_prov = RD.fast_ols_fe_prov_boot(Y_mat, x_loss, Z, idx_p, picks)
        for j, (ml, pl, oname) in enumerate(col_meta):
            sz = RD.summarize_boot(b_full[j], b_zone[:, j])
            sp = RD.summarize_boot(b_full[j], b_prov[:, j])
            r_wcb = RD.R1.fe_ols(Y_mat[:, j], x_loss, Z, provs_s, provs_s, n_boot=RD.N_WCB_OLS,
                                 seed=int(RD.rng_for("d_inl_wcb", sname, ml, pl, oname).integers(1, 2**31)))
            rows.append({
                "sample": sname, "mode": ml, "period": pl, "outcome": oname,
                "n": len(d_sub), "G": r_wcb["G"], "coef": float(b_full[j]),
                "zone_boot_ci_lo": sz["ci_lo"], "zone_boot_ci_hi": sz["ci_hi"], "zone_boot_p": sz["p_boot"],
                "prov_boot_ci_lo": sp["ci_lo"], "prov_boot_ci_hi": sp["ci_hi"],
                "p_wcr_prov": float(r_wcb["p_wcr"]),
                "ci_wcu_prov_lo": float(r_wcb["ci_wcu_lo"]), "ci_wcu_prov_hi": float(r_wcb["ci_wcu_hi"]),
            })
    df_inl = pd.DataFrame(rows)
    df_inl.to_csv(OUT / "d_decomp_controlled_inland.csv", index=False)
    print(df_inl[df_inl["period"].isin(["2024-08-15", "2024-12-25", "Aug15_minus_Feb14", "Dec25_minus_Feb14"])][
        ["sample", "mode", "period", "outcome", "coef", "zone_boot_ci_lo", "zone_boot_ci_hi", "ci_wcu_prov_lo", "ci_wcu_prov_hi"]
    ].to_string(index=False))

    # Also summarise all 306 cached nights by calendar period & year for Q5 (all & single_mun)
    qa = pd.read_csv(OUT / "d_quintiles_all_nights.csv")
    q5 = qa[qa["cls"] == "Q5"].copy()
    print("\n=== Q5 across all 306 cached nights by year-month (all sample) ===")
    q5_all = q5[q5["sample"] == "all"].copy()
    q5_all["ym"] = q5_all["date"].str[:7]
    print(q5_all.groupby("ym").agg(
        n_nights=("date", "count"),
        imp_min=("share_net_importer", "min"),
        imp_med=("share_net_importer", "median"),
        imp_max=("share_net_importer", "max"),
        arr_med=("mean_arrivals_per_res", "median"),
        dep_med=("mean_departures_per_res", "median"),
        net_med=("mean_net_per_res", "median"),
    ).to_string())

    # Check 2022 & 2023 main nights & sister nights from d_quintiles_main_boot.csv
    qb = pd.read_csv(OUT / "d_quintiles_main_boot.csv")
    print("\n=== Q5 levels on 2022, 2023, 2024 key nights (all & single_mun) ===")
    sub = qb[(qb["mode"] == "level") & (qb["cls"] == "Q5") & (qb["period"].isin([
        "2024-02-14", "2024-02-mean13", "2024-10-mean12", "2024-08-15", "2024-08-mean9",
        "2024-11-01", "2024-12-25", "2024-03-29",
        "2022-02-16", "2022-08-15", "2022-11-01", "2023-02-15", "2023-08-15",
    ])) & (qb["metric"].isin(["share_net_importer", "arrivals_per_res", "departures_per_res", "net_per_res"]))]
    piv = sub.pivot_table(index=["sample", "period"], columns="metric", values="mean")
    print(piv[["share_net_importer", "arrivals_per_res", "departures_per_res", "net_per_res"]].to_string())


if __name__ == "__main__":
    main()

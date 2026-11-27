"""Read-only checks of eight proposed headline results and the INE workbook import.

Run with the coursework environment's .venv/bin/python -I. Prints aggregate
results only. Pass a different absolute coursework checkout as the first argument.
"""
from pathlib import Path
import sys
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import duckdb
import numpy as np
import pandas as pd
import pyfixest as pf


ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else
            "/Users/alejandrotreny/Proyectos/spain-holiday-roots")
RES = ROOT / "results"
EXT = ROOT / "_local" / "external"
NET = RES / "b01_networks"
NODES = pd.read_parquet(NET / "nodes.parquet").sort_values("idx")
FEB, AUG, XMAS = "2024-02-14", "2024-08-15", "2024-12-25"
CON = duckdb.connect()
CON.execute("SET threads TO 2")


def csv(name, **kwargs):
    return pd.read_csv(RES / (name + ".csv"), **kwargs)


def edges(date, backbone="full"):
    name = f"edges_{date}__{backbone}.parquet"
    path = NET / name
    return pd.read_parquet(path if path.exists() else RES / "r1_cache" / name)


def graph_stats(date, backbone):
    edge = edges(date, backbone)
    directed = pd.MultiIndex.from_frame(edge[["i", "j"]])
    reverse = pd.MultiIndex.from_arrays([edge.j, edge.i])
    undirected = np.unique(np.sort(edge[["i", "j"]].values, axis=1), axis=0)
    n = len(NODES)
    degree = np.bincount(undirected.ravel(), minlength=n)
    ends = degree[undirected]
    assortativity = ((ends[:, 0] * ends[:, 1]).mean() - ends.mean() ** 2) / ends.var()
    kin = pd.Series(np.bincount(edge.j, minlength=n))
    kout = pd.Series(np.bincount(edge.i, minlength=n))
    return 2 * len(undirected) / (n * (n - 1)), kin.corr(kout, method="spearman"), directed.isin(reverse).mean(), assortativity


def make_dyads(cities):
    provinces = sorted(NODES.prov.unique())
    d = cities[["city", "income"]].merge(pd.DataFrame({"dest": provinces}), how="cross")
    d["orig"] = d.city.str[:2]
    d = d[d.orig != d.dest].copy()
    d["pair"] = d.orig + "_" + d.dest
    nodes = NODES.set_index("id")
    centroids = {
        province: np.average(s[["x", "y"]], weights=s.population, axis=0)
        for province, s in NODES.groupby("prov")
    }
    xy = nodes.loc[d.city, ["x", "y"]].to_numpy()
    dest_xy = np.array([centroids[p] for p in d.dest])
    d["ld"] = np.log(np.linalg.norm(xy - dest_xy, axis=1) / 1000)
    d["ld2"] = d.ld ** 2
    income = np.log(d.income)
    income -= income.mean()
    slopes = [f"inc_{p}" for p in provinces[1:]]
    d = d.assign(**{col: income * (d.dest == province)
                    for col, province in zip(slopes, provinces[1:])})
    return d, ["ld", "ld2", *slopes]


def fit_nights(d, terms, covariates, dates, october_mean=False):
    formula = "flow ~ " + " + ".join([*terms, *covariates]) + " | city + pair"
    pairs = pd.MultiIndex.from_frame(d[["city", "dest"]])
    estimates = []
    outcomes = []
    for date in dates:
        edge = edges(date)
        city = NODES.id.to_numpy()[edge.i]
        prov = NODES.prov.to_numpy()[edge.j]
        flow = edge.groupby([city, prov]).w.sum()
        outcomes.append(flow.reindex(pairs, fill_value=0).to_numpy())
    if october_mean:
        raw = ROOT / "data/derived/pernoct_mun_2024_6m.parquet"
        october = CON.execute("""
            SELECT r AS city, left(s, 2) AS dest, sum(p) AS persons,
                   count(DISTINCT date) AS populated_nights
            FROM read_parquet(?)
            WHERE date >= DATE '2024-10-01' AND date < DATE '2024-11-01'
              AND dayofweek(date) IN (2, 3, 4) AND r IN ? AND s IN ? AND r != s
            GROUP BY 1, 2
        """, [str(raw), list(d.city.unique()), NODES.id.tolist()]).df()
        assert october.populated_nights.max() == 15
        outcomes.append(october.set_index(["city", "dest"]).persons.reindex(pairs, fill_value=0).to_numpy() / 15)
    for outcome in outcomes:
        per_night = d.assign(flow=outcome)
        estimate = pf.fepois(formula, per_night, iwls_tol=1e-10, iwls_maxiter=100).coef()
        estimates.append(estimate.loc[terms].to_numpy())
    return np.array(estimates)


def check_roots():
    r2 = csv("r2_city_set", dtype={"ine_code": str})
    cities = r2[r2.identifying].rename(columns={"ine_code": "city", "income_2024": "income"})
    dyads, cov = make_dyads(cities)
    roots = csv("r2_roots_city_prov", dtype={"ine_code": str, "prov_birth": str})
    dyads = dyads.merge(roots.rename(columns={"ine_code": "city", "prov_birth": "dest"})[
        ["city", "dest", "share"]], validate="one_to_one")
    dyads["roots"] = np.log(dyads.share)
    nights = csv("r2_ppml_nights").query("sample == 'all' and spec == 'main'").set_index("date")
    beta = fit_nights(dyads, ["roots"], cov, [FEB, AUG], october_mean=True).ravel()
    expected = nights.loc[[FEB, AUG, "2024-10-mean15"], "beta"].to_numpy()
    np.testing.assert_allclose(beta, expected, atol=1e-6)
    print("2. Pooled PPML (Feb, Aug):", beta[:2], "Aug-Feb:", beta[1] - beta[0])
    diff = csv("r2_diff_tests").query("sample == 'all' and spec == 'main' and holiday == @AUG")
    for baseline, ordinary in (("Feb", FEB), ("Oct", "2024-10-mean15")):
        row = diff[diff.ordinary == ordinary].iloc[0]
        assert np.isclose(row.delta, row.beta_holiday - row.beta_ordinary)
        print("   ", baseline, row.delta, row.p_wcr_delta0, "(published WCR p-value)")


def check_cohorts():
    census = EXT / "ine_censo_microdatos_2024/md_CensoAnual_muestra_2024.parquet"
    counts = CON.execute("""
        SELECT CPRO || CMUN AS city, CPRO_NAC AS dest, count(*) AS people,
            count(*) FILTER (WHERE VARANOP = '0001') AS early,
            count(*) FILTER (WHERE VARANOP BETWEEN '2016' AND '2023') AS recent,
            count(*) FILTER (WHERE CMUN_NAC = '991') AS village,
            count(*) FILTER (WHERE CMUN_NAC NOT IN
                ('991','992','993','666','966','777','888','999')) AS urban,
            count(*) FILTER (WHERE VAREDAD < 65 AND VARANOP = '0001') AS early_young,
            count(*) FILTER (WHERE VAREDAD >= 65 AND VARANOP = '0001') AS early_older,
            count(*) FILTER (WHERE VAREDAD >= 65 AND VARANOP != '0001') AS late_older,
            count(*) FILTER (WHERE VAREDAD >= 65 AND (VARANOP BETWEEN '0002' AND '0005'
                OR VARANOP BETWEEN '2016' AND '2023')) AS late_older_dated,
            count(*) FILTER (WHERE VAREDAD < 65 AND VARANOP BETWEEN '2016' AND '2023') AS recent_young
        FROM read_parquet(?)
        WHERE CMUN NOT IN ('991','992','993','666','966','777','888','999')
        GROUP BY 1, 2
    """, [str(census)]).df()
    eligible = NODES.loc[NODES.n_municipalities == 1, ["id", "prov"]]
    eligible = eligible[eligible.id.isin(counts.city)]
    cities = eligible[eligible.groupby("prov").id.transform("size") >= 2].rename(columns={"id": "city"})
    assert len(cities) == 675
    income = pd.read_csv(EXT / "ine_adrh/adrh_municipios.csv", dtype={"ine_code": str})
    income = income.query("year == 2024 and indicator == 'Renta neta media por persona'")
    cities["income"] = cities.city.map(income.set_index("ine_code").value)
    dyads, cov = make_dyads(cities)
    dyads = dyads.merge(counts.drop(columns="people"), how="left", validate="one_to_one")
    terms = ["early", "recent", "village", "urban"]
    age_terms = ["early_young", "early_older", "late_older", "recent_young"]
    late_all = int(dyads.late_older.sum())
    late_dated = int(dyads.late_older_dated.sum())
    city_population = dyads.city.map(counts.groupby("city").people.sum())
    for term in [*terms, *age_terms, "late_older_dated"]:
        dyads[term] = np.log((dyads[term].fillna(0) + 0.5) / city_population)
    result = []
    for group in (terms[:2], terms[2:]):
        feb, aug = fit_nights(dyads, group, cov, [FEB, AUG])
        result.extend(aug - feb)
    table = csv("v1_c2_main_contrasts").query("contrast == 'Aug15_vs_Feb14' and sample == 'all'")
    np.testing.assert_allclose(result, table.delta, atol=1e-6)
    print("3. Cohort and birthplace-size PPML (early, recent, village, urban):", np.array(result))
    feb, aug = fit_nights(dyads, age_terms, cov, [FEB, AUG])
    original = aug - feb
    age_table = csv("v1_c2_age_confound").query("spec == 'E2a_split65_4way_exclusive' and sample == 'all'")
    np.testing.assert_allclose(original, age_table.delta_aug_feb, atol=1e-6)
    dated_terms = ["early_young", "early_older", "late_older_dated", "recent_young"]
    feb, aug = fit_nights(dyads, dated_terms, cov, [FEB, AUG])
    print("   Age split (published coding / dated arrivals only):", original, aug - feb)
    print("   Older comparison, special codes:", late_all - late_dated, "/", late_all)


def check_villages():
    eligible = NODES.loss_peak_2021.notna() & ~NODES.id.isin(["31905_AM", "19103_AM", "40149_AM"])
    q5 = pd.qcut(NODES.loc[eligible, "loss_peak_2021"], 5, labels=False) == 4
    q5 = q5[q5].index
    table = csv("r1_net_balance_by_quintile").query("sample == 'all' and role == 'main' and cls == 'Q5'").set_index("date")
    for d in [FEB, AUG, XMAS]:
        e = edges(d)
        incoming = np.bincount(e.j, weights=e.w, minlength=len(NODES))
        outgoing = np.bincount(e.i, weights=e.w, minlength=len(NODES))
        count = (incoming[q5] > outgoing[q5]).sum()
        assert np.isclose(count / len(q5), table.loc[d, "share_net_importer"])
        print("4. Q5 net receivers:", d, count, "/", len(q5), "=", count / len(q5))


def check_workbook():
    folder = EXT / "ine_turismo_moviles_od"
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    rel = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
    with ZipFile(folder / "raw/exp_tmov_interno_mun_2024.xlsx") as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = {entry.attrib["Id"]: entry.attrib["Target"]
                         for entry in ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))}
        sheets = {entry.attrib["name"].replace("_", "-"): relationships[entry.attrib[rel + "id"]]
                  for entry in workbook.iter(ns + "sheet")}
        strings = ["".join(entry.itertext())
                   for entry in ET.fromstring(archive.read("xl/sharedStrings.xml"))]
        for month in ["2024-02", "2024-08"]:
            count, total = 0, 0
            with archive.open("xl/" + sheets[month]) as raw:
                column = None
                for _, element in ET.iterparse(raw, events=("end",)):
                    if element.tag != ns + "row":
                        continue
                    row = element.attrib["r"]
                    if row == "1":
                        column = next(c.attrib["r"][:-1] for c in element
                                      if c.attrib.get("t") == "s"
                                      and strings[int(c.findtext(ns + "v"))] == "turistas")
                    else:
                        value = element.findtext(f"{ns}c[@r='{column}{row}']/{ns}v")
                        assert value is not None
                        count += 1
                        total += int(value)
                    element.clear()
            cached = CON.execute("SELECT count(*), sum(turistas) FROM read_parquet(?) WHERE mes = ?",
                                 [str(folder / "ine_od_mun_2024.parquet"), month]).fetchone()
            print("INE workbook", month, "rows/travellers", (count, total), "parquet", cached)
            assert (count, total) == cached


def check_municipal_support():
    municipal = RES / "n1_mun_nodes.parquet"
    trips = EXT / "ine_turismo_moviles_od/ine_od_mun_2024.parquet"
    roots = RES / "n1_roots_city_village_ecepov2021.parquet"
    support = CON.execute("""
        WITH nodes AS (SELECT * FROM read_parquet(?)),
        cities AS (SELECT ine_code AS orig, prov FROM nodes WHERE is_census_city),
        trip AS (
            SELECT t.orig, t.dest, t.mes FROM read_parquet(?) t
            JOIN cities c USING (orig)
            JOIN nodes n ON t.dest = n.ine_code
            WHERE mes IN ('2024-02', '2024-03', '2024-08', '2024-10', '2024-11', '2024-12')
        ),
        departures AS (
            SELECT orig, dest,
                bool_or(mes = '2024-02') AS feb,
                bool_or(mes = '2024-08') AS aug
            FROM trip GROUP BY 1, 2
        )
        SELECT coalesce(r.born_g2, 0) > 0 AS sampled_parental_tie,
            count(*) AS eligible_dyads,
            count(*) FILTER (WHERE d.feb) AS feb_published,
            count(*) FILTER (WHERE d.aug) AS aug_published,
            count(*) FILTER (WHERE d.feb AND d.aug) AS both_published
        FROM cities c
        CROSS JOIN (SELECT DISTINCT dest FROM trip) a
        JOIN nodes n ON a.dest = n.ine_code
        LEFT JOIN read_parquet(?) r ON c.orig = r.orig AND a.dest = r.dest
        LEFT JOIN departures d ON c.orig = d.orig AND a.dest = d.dest
        WHERE c.prov != n.prov AND NOT n.coastal_true
        GROUP BY 1 ORDER BY 1
    """, [str(municipal), str(trips), str(roots)]).df()
    print("   Municipal inland disclosure diagnostic on N1's pre-separation frame:")
    print(support.to_string(index=False))


def check_households():
    folder = EXT / "ine_ecepov_2021"
    hh = pd.read_parquet(folder / "ecepov_vivienda_2021.parquet",
                         columns=["IDEN", "SEGUNRESI", "LUGSEGUNRESI", "OTRACCAA", "FACTOR"])
    hh = hh[hh.SEGUNRESI.eq("1") & hh.LUGSEGUNRESI.eq("4")].copy()
    cols = ["IDEN", "NACIM", "PROVN", "NACIMP1", "PROVNP1", "NACIMP2", "PROVNP2"]
    persons = pd.read_parquet(folder / "ecepov_hogar_2021.parquet", columns=cols)
    merged = persons.merge(hh[["IDEN", "OTRACCAA"]], validate="many_to_one")
    provinces = [
        "04 11 14 18 21 23 29 41", "22 44 50", "33", "07", "35 38", "39",
        "05 09 24 34 37 40 42 47 49", "02 13 16 19 45", "08 17 25 43", "03 12 46",
        "06 10", "15 27 32 36", "28", "30", "31", "01 20 48", "26", "51", "52",
    ]
    community = {prov: f"{i:02d}" for i, group in enumerate(provinces, 1) for prov in group.split()}
    hit = np.logical_or.reduce([
        merged[nativity].eq("1") & merged[birth].map(community).eq(merged.OTRACCAA)
        for nativity, birth in [("NACIM", "PROVN"), ("NACIMP1", "PROVNP1"), ("NACIMP2", "PROVNP2")]
    ])
    hh["match"] = hh.IDEN.isin(merged.loc[hit, "IDEN"])
    table = csv("v1_c4_ecepov").set_index("subgroup")
    for sample, name in [(hh, "all_cross_ccaa"),
                         (hh[hh.OTRACCAA.isin(["02", "07", "08", "11", "13", "15", "17"])],
                          "inland_dest_ccaa_7_incl_Madrid")]:
        share = np.average(sample.match, weights=sample.FACTOR)
        assert len(sample) == table.loc[name, "n_hh"]
        assert np.isclose(share, table.loc[name, "match_any_weighted"])
        print("6. ECEPOV weighted household match:", name, len(sample), share)


def main():
    table = csv("b02_basic_properties")
    for date, bb in [(FEB, "w10"), (AUG, "w10"), (AUG, "equal_edges")]:
        density, spearman, recip, r = graph_stats(date, bb)
        check = table.query("date == @date and backbone == @bb").groupby("metric").real.first()
        np.testing.assert_allclose([recip, r], check.loc[["reciprocity", "r_und"]], atol=1e-9)
        print("1. Graph:", date, bb, density, spearman, recip, r)
    sc = csv("b02_structural_constraints").query("backbone == 'w10'").set_index("date")
    shifts = sc.loc[[FEB, AUG], ["r_real", "maxrand_conn_r_mean"]].diff()
    print("   Assortativity shift share:", shifts.loc[AUG, "maxrand_conn_r_mean"] / shifts.loc[AUG, "r_real"])
    check_roots()
    check_cohorts()
    check_villages()
    housing = csv("v1_c3_mediation").set_index("spec")
    for label in ["MITMA_Aug15_vs_Feb14_log_sin", "MITMA_Dec25_vs_Feb14_log_sin"]:
        h = housing.loc[label]
        share = 1 - h.beta_after_med / h.beta_base
        assert np.isclose(share, h.decomp_share)
        print("5. Housing slope attenuation:", label, share)
    check_households()
    municipal = csv("n1_city_village_roots").query("term == 'lg2_h' and spec.str.endswith('inland_g2')")
    print("7. Municipal second-generation contrasts (CSV trace):",
          municipal.set_index("mes").beta.to_dict())
    check_workbook()
    check_municipal_support()
    calendar = csv("p_daily_metrics")
    available = set()
    for name in ["pernoct_mun_2024_6m.parquet", "pernoct_mun_2022_2023.parquet"]:
        raw = ROOT / "data/derived" / name
        available.update(row[0] for row in CON.execute(
            "SELECT DISTINCT CAST(date AS VARCHAR) FROM read_parquet(?)", [str(raw)]
        ).fetchall())
    assert set(calendar.date) == available and calendar.date.is_unique
    print("8. Calendar counts by year:", calendar.date.str[:4].value_counts().to_dict())


if __name__ == "__main__":
    main()

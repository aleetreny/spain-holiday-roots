"""p_patterns: exploratory pattern mining for question B ("Roots or beaches?").

Wave-2 agent P. Descriptive, systematic search for patterns that could change,
strengthen or threaten the story. Every number in results/p_patterns.md comes
from a table written here.

Stages (run all, or name some on the command line):
  daily     one network per night for EVERY 2024 night in the data (180) plus the 2022-23
            nights; network + flow metrics per night -> results/p_daily_metrics.csv
  zones     per-zone per-night totals (residents, stays, in/out away strength)
            -> results/p_zone_daily.parquet (cache) and artefact tables
  hubs      degree distribution, in/out hubs, PageRank, betweenness, net balance on the
            main/sister/year nights -> results/p_hubs_*.csv, p_centrality_*.csv
  cities    rooted vs beach ranking of the 137 mainland INE cities -> results/p_cities_*.csv
  holidays  Easter / Christmas / All Saints / weekends vs August (excess similarity)
            -> results/p_holiday_*.csv
  distance  distance distribution of stays and of the holiday excess -> results/p_distance_*.csv
  plots     quick diagnostic plots -> results/p_diag/

Conventions
-----------
* Loaders: b00_common (node table N = 2,435, build_edges, backbones, to_igraph).
* Weights are persons (MITMA expanded). Edge i->j: residents of i sleeping in j, i != j.
* "Ordinary baseline" for a pair = MEAN over the 13 Feb-2024 Tue/Wed/Thu sister nights with a
  missing pair counted as 0 (b00.window_mean_from_frame), never sum/count(distinct date) (audit C3).
  Node-level comparisons with February use the per-node MEDIAN over the same 13 nights.
  "Excess" of a night = w(night) - w_base, per pair (can be negative).
* Coastal flag: coastal_true (config). "Emptied" = loss_peak_2021 >= 0.5 (lost at least half
  of its 1900-1960 peak). Descriptive label only; R1 owns the node-level test.
* Distances: Euclidean between zone representative points (EPSG:3042), km.

Run: .venv/bin/python -I pipeline/p_patterns.py [stage ...]
"""
from __future__ import annotations

import sys
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path

import duckdb
import igraph as ig
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import b00_common as B  # noqa: E402

RES = C.RESULTS
DIAG = RES / "p_diag"
NETS = RES / "b01_networks"
EMPTIED_LOSS = 0.5
W_MIN = C.BACKBONE_MAIN["w_min"]
E_EQUAL = 35812  # equal-edges E (decisions.md, W1-core)

# ---- 2024 calendar (national holidays + the main regional ones, and the eves) ------------
# Source: BOE resolution of 2023-10-05 (fiestas laborales 2024). Labels are descriptive.
HOLIDAYS_2024 = {
    "2024-02-13": "Carnival Tue (local)",
    "2024-02-28": "Andalucia Day (AND)",
    "2024-03-19": "San Jose (VAL, MUR, partly)",
    "2024-03-28": "Holy Thursday (all but CAT, VAL)",
    "2024-03-29": "Good Friday (national)",
    "2024-08-15": "Assumption (national)",
    "2024-10-12": "National Day (Sat)",
    "2024-11-01": "All Saints (national, Fri)",
    "2024-12-06": "Constitution Day (Fri)",
    "2024-12-08": "Immaculate Conception (Sun; Mon 9 in some regions)",
    "2024-12-25": "Christmas Day",
    "2024-12-24": "Christmas Eve (not a holiday)",
    "2024-12-31": "New Year's Eve (not a holiday)",
}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def nodes() -> pd.DataFrame:
    return B.load_nodes()


def gini(x) -> float:
    x = np.sort(np.asarray(x, float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return np.nan
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def wquantile(v, w, q) -> float:
    o = np.argsort(v)
    v, w = np.asarray(v)[o], np.asarray(w)[o]
    c = np.cumsum(w) / w.sum()
    return float(np.interp(q, c, v))


def dist_km(i, j, nd) -> np.ndarray:
    x, y = nd["x"].values, nd["y"].values
    return np.hypot(x[i] - x[j], y[i] - y[j]) / 1000.0


def spearman(a, b) -> float:
    a = pd.Series(a).rank().values
    b = pd.Series(b).rank().values
    return float(np.corrcoef(a, b)[0, 1])


def load_edges(date: str, kind: str = "full") -> pd.DataFrame:
    """Saved b01 edges when they exist, else rebuilt with b00.build_edges."""
    f = NETS / f"edges_{date}__{kind}.parquet"
    if f.exists():
        return pd.read_parquet(f)
    e, _ = B.build_edges(date, nodes())
    if kind == "full":
        return e
    if kind == "w10":
        return B.backbone_threshold(e, W_MIN)
    if kind == "equal_edges":
        return B.backbone_equal_edges(e, E_EQUAL)
    raise ValueError(kind)


def all_dates() -> list[str]:
    d24 = B.available_dates(C.PERNOCT_2024)
    d23 = B.available_dates(C.PERNOCT_2022_23)
    return sorted(set(d24) | set(d23))


def ordinary_sisters() -> list[str]:
    return sorted(pd.read_csv(NETS / "networks.csv").query("group=='ordinary' and role=='sister'").date.tolist())


NATIONAL_2024 = {"2024-03-28", "2024-03-29", "2024-08-15", "2024-10-12", "2024-11-01",
                 "2024-12-06", "2024-12-08", "2024-12-25"}   # Holy Thursday: all regions but CAT, VAL


def day_type(date: str) -> str:
    """Label of the calendar date as stored in the file (see the date-convention check)."""
    if date in NATIONAL_2024:
        return "holiday"
    return ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"][pd.Timestamp(date).day_of_week]


def period(date: str) -> str:
    d = pd.Timestamp(date)
    if d.year != 2024:
        return f"{d.year}-{d.month:02d}"
    if d.month == 3 and d.day >= 23:
        return "Easter week"
    if d.month == 12 and d.day >= 21:
        return "Christmas"
    if d.month == 8:
        return "August"
    if (d.month == 10 and d.day == 31) or (d.month == 11 and d.day <= 3):
        return "All Saints"
    return {2: "Feb", 3: "Mar", 10: "Oct", 11: "Nov", 12: "Dec"}[d.month]


# ---------------------------------------------------------------------------
# roots matrix (city x birth province), used for a descriptive daily index and the city ranking
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def roots_matrix(year: int = 2024) -> pd.DataFrame:
    """Long table: city zone id, city idx, prov_birth, persons (Census 68545), mainland cities only,
    Spanish provinces only (ABROAD / FORMER_ES dropped), own province kept (flag own)."""
    nd = nodes()
    m = pd.read_csv(C.ROOT / "_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv",
                    dtype={"ine_code": str, "prov_birth": str})
    m = m[(m.year == year) & m.prov_birth.str.fullmatch(r"\d\d")]
    m = m[~m.ine_code.str[:2].isin(C.EXCLUDED_PROVINCES) & ~m.prov_birth.isin(C.EXCLUDED_PROVINCES)]
    pos = dict(zip(nd.id, nd.idx))
    m = m[m.ine_code.isin(pos)].copy()
    m["c"] = m.ine_code.map(pos)
    m["own"] = m.ine_code.str[:2] == m.prov_birth
    assert m.ine_code.nunique() == 137, m.ine_code.nunique()
    return m[["ine_code", "municipality_name", "c", "prov_birth", "persons", "own"]]


@lru_cache(maxsize=1)
def prov_centroids() -> pd.DataFrame:
    nd = nodes()
    w = nd.population.fillna(nd.population.median())
    g = pd.DataFrame({"prov": nd.prov, "x": nd.x * w, "y": nd.y * w, "w": w}).groupby("prov").sum()
    return pd.DataFrame({"x": g.x / g.w, "y": g.y / g.w})


def city_prov_flows(e: pd.DataFrame) -> pd.DataFrame:
    """Flows from the 137 cities to each destination province, own province excluded."""
    nd = nodes()
    R = roots_matrix()
    cities = R[["c", "ine_code"]].drop_duplicates()
    d = e[e.i.isin(cities.c)].copy()
    d["p"] = nd.prov.values[d.j.values]
    f = d.groupby(["i", "p"], as_index=False).w.sum().rename(columns={"i": "c", "p": "prov_birth", "w": "flow"})
    full = R[~R.own][["c", "prov_birth", "persons"]].merge(f, on=["c", "prov_birth"], how="left").fillna({"flow": 0.0})
    pc = prov_centroids()
    cx, cy = nd.x.values[full.c.values], nd.y.values[full.c.values]
    full["dist_km"] = np.hypot(cx - pc.x.reindex(full.prov_birth).values, cy - pc.y.reindex(full.prov_birth).values) / 1000
    return full


def roots_index(e: pd.DataFrame) -> dict:
    """Descriptive roots alignment of one night: correlation over (city, other province) pairs of
    log(1+flow) with log(1+born), raw and after removing city means and log distance (OLS
    residuals). Not a test (R2 owns the inference)."""
    f = city_prov_flows(e)
    y = np.log1p(f.flow.values)
    x = np.log1p(f.persons.values)
    raw = float(np.corrcoef(x, y)[0, 1])
    # within-city, distance-partialled
    X = np.column_stack([np.log(f.dist_km.values + 1)])
    def resid(v):
        v = v - pd.Series(v).groupby(f.c.values).transform("mean").values
        Xd = X - pd.DataFrame(X).groupby(f.c.values).transform("mean").values
        b = np.linalg.lstsq(Xd, v, rcond=None)[0]
        return v - Xd @ b
    part = float(np.corrcoef(resid(x), resid(y))[0, 1])
    share = float((f.flow * f.persons / f.groupby("c").persons.transform("sum")).sum() / f.flow.sum())
    return {"roots_corr_raw": raw, "roots_corr_partial": part, "roots_weighted_share": share}


# ---------------------------------------------------------------------------
# stage: daily
# ---------------------------------------------------------------------------
def _undirected_from(eb: pd.DataFrame, n: int) -> ig.Graph:
    return B.to_igraph(eb, n, directed=False)


def night_metrics(date: str) -> dict:
    nd = nodes()
    n = len(nd)
    e, info = B.build_edges(date, nd)
    w = e.w.values
    i, j = e.i.values, e.j.values
    dkm = dist_km(i, j, nd)
    coast = nd.coastal_true.values
    loss = nd.loss_peak_2021.values
    pop = nd.population.values
    prov = nd.prov.values
    ccaa = nd.ccaa.values
    emptied_inland = (loss >= EMPTIED_LOSS) & ~coast
    small_inland = (pop < C.SMALL_VILLAGE_POP) & ~coast
    away = w.sum()
    out = {"date": date, **{k: v for k, v in info.items() if k != "date"}}
    out.update({
        "away_persons": away,
        "n_pairs": len(e),
        "dist_mean_km": float((w * dkm).sum() / away),
        "dist_median_km": wquantile(dkm, w, 0.5),
        "dist_p90_km": wquantile(dkm, w, 0.9),
        "share_lt50km": float(w[dkm < 50].sum() / away),
        "share_50_150km": float(w[(dkm >= 50) & (dkm < 150)].sum() / away),
        "share_150_300km": float(w[(dkm >= 150) & (dkm < 300)].sum() / away),
        "share_ge300km": float(w[dkm >= 300].sum() / away),
        "share_dest_coastal": float(w[coast[j]].sum() / away),
        "share_dest_emptied_inland": float(w[emptied_inland[j]].sum() / away),
        "share_dest_small_inland": float(w[small_inland[j]].sum() / away),
        "share_same_prov": float(w[prov[i] == prov[j]].sum() / away),
        "share_same_ccaa": float(w[ccaa[i] == ccaa[j]].sum() / away),
        # "downhill": host zone smaller than residence zone (population hierarchy)
        "share_downhill": float(w[np.nan_to_num(pop[j], nan=0) < np.nan_to_num(pop[i], nan=0)].sum() / away),
    })
    # weighted reciprocity on the full network: sum min(w_ij, w_ji) / sum w
    key = i.astype(np.int64) * n + j
    rkey = j.astype(np.int64) * n + i
    s = pd.Series(w, index=key)
    wr = s.reindex(rkey).fillna(0).values
    out["recip_weighted_full"] = float(np.minimum(w, wr).sum() / away)
    # strengths
    s_in = np.bincount(j, w, n)
    s_out = np.bincount(i, w, n)
    out["gini_s_in"] = gini(s_in)
    out["gini_s_out"] = gini(s_out)
    out["top10_out_share"] = float(np.sort(s_out)[-10:].sum() / away)
    out["top10_in_share"] = float(np.sort(s_in)[-10:].sum() / away)
    # backbones
    for bb, eb in (("w10", B.backbone_threshold(e, W_MIN)), ("eq", B.backbone_equal_edges(e, E_EQUAL))):
        Gd = B.to_igraph(eb, n, directed=True)
        Gu = _undirected_from(eb, n)
        kin = np.array(Gd.indegree())
        kout = np.array(Gd.outdegree())
        ku = np.array(Gu.degree())
        out.update({
            f"{bb}_E": Gd.ecount(),
            f"{bb}_E_und": Gu.ecount(),
            f"{bb}_density_und": 2 * Gu.ecount() / (n * (n - 1)),
            f"{bb}_share_weight": float(eb.w.sum() / away),
            f"{bb}_kin_max": int(kin.max()),
            f"{bb}_kout_max": int(kout.max()),
            f"{bb}_kin_hub": nd.name.values[kin.argmax()],
            f"{bb}_kout_hub": nd.name.values[kout.argmax()],
            f"{bb}_gini_kin": gini(kin),
            f"{bb}_gini_kout": gini(kout),
            f"{bb}_spearman_kin_kout": spearman(kin, kout),
            f"{bb}_r_und": Gu.assortativity_degree(directed=False),
            f"{bb}_recip": Gd.reciprocity(),
            f"{bb}_C_avg": Gu.transitivity_avglocal_undirected(mode="zero"),
            f"{bb}_kmax_core": int(max(Gu.coreness())),
            f"{bb}_n_isolated": int((ku == 0).sum()),
        })
    out.update(roots_index(e))
    return out


def stage_daily() -> pd.DataFrame:
    dates = all_dates()
    with Pool(8) as pool:
        rows = pool.map(night_metrics, dates, chunksize=2)
    df = pd.DataFrame(rows).sort_values("date")
    df["weekday"] = [pd.Timestamp(d).day_name()[:3] for d in df.date]
    df["day_type"] = [day_type(d) for d in df.date]
    df["period"] = [period(d) for d in df.date]
    df["holiday_label"] = df.date.map(HOLIDAYS_2024).fillna("")
    df.to_csv(RES / "p_daily_metrics.csv", index=False)
    return df


# ---------------------------------------------------------------------------
# stage: zones (per-zone daily totals) + artefacts
# ---------------------------------------------------------------------------
def stage_zones() -> None:
    nd = nodes()
    ids = list(nd.id)
    con = duckdb.connect()
    con.register("ids", pd.DataFrame({"id": ids}))
    frames = []
    for path in (C.PERNOCT_2024, C.PERNOCT_2022_23):
        q = f"""
        WITH d AS (SELECT date, r, s, sum(p) p FROM read_parquet('{path}') GROUP BY 1,2,3)
        SELECT date, r AS id,
               sum(p) AS residents_all,                                 -- residents of r anywhere (incl. home)
               sum(CASE WHEN r = s THEN p ELSE 0 END) AS at_home,
               sum(CASE WHEN r <> s AND s IN (SELECT id FROM ids) THEN p ELSE 0 END) AS out_away_ml
        FROM d WHERE r IN (SELECT id FROM ids) GROUP BY 1, 2"""
        a = con.execute(q).df()
        q2 = f"""
        WITH d AS (SELECT date, r, s, sum(p) p FROM read_parquet('{path}') GROUP BY 1,2,3)
        SELECT date, s AS id, sum(p) AS sleepers_all,
               sum(CASE WHEN r <> s AND r IN (SELECT id FROM ids) THEN p ELSE 0 END) AS in_away_ml
        FROM d WHERE s IN (SELECT id FROM ids) GROUP BY 1, 2"""
        b = con.execute(q2).df()
        frames.append(a.merge(b, on=["date", "id"], how="outer"))
    z = pd.concat(frames).fillna(0)
    z["date"] = z.date.astype(str).str[:10]
    z.to_parquet(RES / "p_zone_daily.parquet", index=False)

    # --- artefact table 1: residents_all vs padron population, and day-to-day stability
    z24 = z[z.date.str.startswith("2024")]
    g = z24.groupby("id").residents_all
    st = pd.DataFrame({"res_median": g.median(), "res_min": g.min(), "res_max": g.max(),
                       "res_cv": g.std() / g.mean(), "n_days": g.size()})
    st = st.join(nd.set_index("id")[["name", "population", "prov", "is_am", "n_municipalities", "coastal_true"]])
    st["res_over_pop"] = st.res_median / st.population
    # month-level jump (calibration): median Feb-Mar vs median Aug-Dec
    fm = z24[z24.date < "2024-04"].groupby("id").residents_all.median()
    ad = z24[z24.date > "2024-07"].groupby("id").residents_all.median()
    st["calib_jump_pct"] = 100 * (ad / fm - 1)
    # max deviation of a single day from the zone's own median (robust)
    med = z24.groupby("id").residents_all.transform("median")
    z24 = z24.assign(dev=(z24.residents_all / med - 1))
    worst = z24.loc[z24.groupby("id").dev.apply(lambda s: s.abs().idxmax())][["id", "date", "dev"]].set_index("id")
    st = st.join(worst.rename(columns={"date": "worst_day", "dev": "worst_day_dev"}))
    st.reset_index().to_csv(RES / "p_artefact_zone_residents.csv", index=False)

    # --- artefact table 2: daily national totals and days that break the weekday pattern
    d = pd.read_csv(RES / "p_daily_metrics.csv")
    d24 = d[d.date.str.startswith("2024")].copy()
    d24["month"] = d24.date.str[5:7]
    for col in ("n_pairs", "away_persons", "total_p_all"):
        ref = d24[d24.day_type != "holiday"].groupby(["month", "weekday"])[col].median()
        mad = d24[d24.day_type != "holiday"].groupby(["month", "weekday"])[col].apply(lambda s: (s - s.median()).abs().median())
        k = list(zip(d24.month, d24.weekday))
        d24[f"{col}_rel_to_month_weekday"] = d24[col].values / ref.reindex(k).values - 1
        d24[f"{col}_robust_z"] = (d24[col].values - ref.reindex(k).values) / (1.4826 * mad.reindex(k).values)
    d24.to_csv(RES / "p_artefact_days.csv", index=False)

    # --- artefact table 3: smallest reported flows (is there a reporting floor?)
    con = duckdb.connect()
    fl = con.execute(f"""SELECT date, min(p) pmin, quantile_cont(p, 0.01) p01, quantile_cont(p, 0.5) p50,
                          avg(CASE WHEN p < 5 THEN 1 ELSE 0 END) share_lt5, count(*) n
                          FROM read_parquet('{C.PERNOCT_2024}') WHERE r <> s GROUP BY 1 ORDER BY 1""").df()
    fl["date"] = fl.date.astype(str).str[:10]
    fl.to_csv(RES / "p_artefact_flow_floor.csv", index=False)


# ---------------------------------------------------------------------------
# stage: hubs / degree distribution / centrality / net balance
# ---------------------------------------------------------------------------
def hub_dates() -> pd.DataFrame:
    nets = pd.read_csv(NETS / "networks.csv")
    return nets.drop_duplicates("date")[["date", "group", "role"]].sort_values("date")


def node_metrics_night(date: str) -> pd.DataFrame:
    nd = nodes()
    n = len(nd)
    e = load_edges(date, "full")
    out = pd.DataFrame({"idx": nd.idx, "date": date})
    out["s_in"] = np.bincount(e.j, e.w, n)
    out["s_out"] = np.bincount(e.i, e.w, n)
    out["balance"] = (out.s_in - out.s_out) / (out.s_in + out.s_out).replace(0, np.nan)
    for bb in ("w10", "equal_edges"):
        eb = load_edges(date, bb)
        Gd = B.to_igraph(eb, n, directed=True)
        Gu = B.to_igraph(eb, n, directed=False)
        out[f"kin_{bb}"] = Gd.indegree()
        out[f"kout_{bb}"] = Gd.outdegree()
        out[f"k_{bb}"] = Gu.degree()
        if bb == "w10":
            out["btw_w10"] = Gu.betweenness(directed=False)
            # reciprocated share of each node's out-links
            key = set(zip(eb.i, eb.j))
            rec = np.array([(b, a) in key for a, b in zip(eb.i, eb.j)])
            out["recip_out_w10"] = pd.Series(rec.astype(float)).groupby(eb.i.values).mean().reindex(range(n)).values
    Gw = ig.Graph(n=n, edges=list(zip(e.i.tolist(), e.j.tolist())), directed=True)
    out["pagerank_w"] = Gw.pagerank(weights=e.w.tolist(), damping=0.85)
    return out


def stage_hubs() -> None:
    nd = nodes()
    hd = hub_dates()
    with Pool(6) as pool:
        frames = pool.map(node_metrics_night, hd.date.tolist())
    nm = pd.concat(frames).merge(hd, on="date")
    nm.to_parquet(RES / "p_node_metrics.parquet", index=False)

    attrs = nd.set_index("idx")[["id", "name", "population", "prov", "coastal_true", "loss_peak_2021", "is_am", "n_municipalities"]]
    mains = C.MAIN_NIGHTS
    # ---- top-20 hubs per main night (in, out, pagerank, betweenness)
    rows = []
    for g, d in mains.items():
        x = nm[nm.date == d].set_index("idx").join(attrs)
        for col in ("kin_w10", "kout_w10", "kin_equal_edges", "kout_equal_edges", "pagerank_w", "btw_w10", "s_in"):
            t = x.sort_values(col, ascending=False).head(20)
            for r, (idx, rr) in enumerate(t.iterrows(), 1):
                rows.append({"group": g, "date": d, "metric": col, "rank": r, "id": rr.id, "name": rr["name"],
                             "value": rr[col], "population": rr.population, "coastal": rr.coastal_true,
                             "loss": rr.loss_peak_2021})
    pd.DataFrame(rows).to_csv(RES / "p_hubs_top20.csv", index=False)

    # ---- hub composition and stability across sister / year nights
    rows = []
    for d, sub in nm.groupby("date"):
        x = sub.set_index("idx").join(attrs)
        g = sub.group.iloc[0]
        for col in ("kin_w10", "kout_w10", "kin_equal_edges", "pagerank_w", "btw_w10"):
            top = x.sort_values(col, ascending=False).head(50)
            rows.append({"date": d, "group": g, "role": sub.role.iloc[0], "metric": col,
                         "top50_share_coastal": top.coastal_true.mean(),
                         "top50_share_pop_gt100k": (top.population > 1e5).mean(),
                         "top50_share_emptied_inland": ((top.loss_peak_2021 >= EMPTIED_LOSS) & ~top.coastal_true).mean(),
                         "top50_median_pop": top.population.median(),
                         "top50_ids": ";".join(top.id)})
    comp = pd.DataFrame(rows)
    # Jaccard of top-50 sets vs the main night of the same group
    jac = []
    for _, r in comp.iterrows():
        main = comp[(comp.date == mains[r.group]) & (comp.metric == r.metric)].top50_ids.iloc[0].split(";")
        a = set(r.top50_ids.split(";"))
        jac.append(len(a & set(main)) / len(a | set(main)))
    comp["jaccard_top50_vs_main"] = jac
    comp.drop(columns="top50_ids").to_csv(RES / "p_hubs_composition.csv", index=False)

    # ---- degree distribution summaries per night (in vs out heterogeneity, decoupling)
    rows = []
    for d, sub in nm.groupby("date"):
        for bb in ("w10", "equal_edges"):
            kin, kout = sub[f"kin_{bb}"].values, sub[f"kout_{bb}"].values
            rows.append({"date": d, "group": sub.group.iloc[0], "role": sub.role.iloc[0], "backbone": bb,
                         "kin_mean": kin.mean(), "kin_max": kin.max(), "kout_max": kout.max(),
                         "kin_p99": np.percentile(kin, 99), "kout_p99": np.percentile(kout, 99),
                         "kin_cv": kin.std() / kin.mean(), "kout_cv": kout.std() / kout.mean(),
                         "gini_kin": gini(kin), "gini_kout": gini(kout),
                         "spearman_kin_kout": spearman(kin, kout),
                         "n_kout_gt_half_N": int((kout > len(kin) / 2).sum()),
                         "n_kin_gt_half_N": int((kin > len(kin) / 2).sum())})
    pd.DataFrame(rows).to_csv(RES / "p_degree_summary.csv", index=False)

    # ---- rank shifts holiday vs ordinary, and what explains the gain (zone attributes)
    feb = nm[nm.group.eq("ordinary") & nm.date.str.startswith("2024")].groupby("idx").median(numeric_only=True)
    rows, reg = [], []
    for d, sub in nm.groupby("date"):
        x = sub.set_index("idx")
        for col in ("kin_w10", "kout_w10", "pagerank_w", "btw_w10", "s_in", "balance"):
            rows.append({"date": d, "group": x.group.iloc[0], "role": x.role.iloc[0], "metric": col,
                         "spearman_vs_feb_median": spearman(x[col], feb[col])})
        # OLS: log gain of in-degree / pagerank on zone attributes (incl. aggregation controls, audit C1/M4)
        A = attrs.copy()
        ok = A.population.notna() & A.loss_peak_2021.notna()
        for col in ("kin_w10", "pagerank_w", "s_in"):
            tr = np.log if col == "pagerank_w" else np.log1p   # PageRank > 0 and << 1: log, not log1p
            y = tr(x[col]) - tr(feb[col])
            X = pd.DataFrame({"const": 1.0, "coastal": A.coastal_true.astype(float),
                              "emptied": ((A.loss_peak_2021 >= EMPTIED_LOSS) & ~A.coastal_true).astype(float),
                              "log_pop": np.log(A.population), "log_nmun": np.log(A.n_municipalities),
                              "log_feb": tr(feb[col])})[ok]
            yy = y[ok].values
            b, *_ = np.linalg.lstsq(X.values, yy, rcond=None)
            res = yy - X.values @ b
            XtX = np.linalg.inv(X.values.T @ X.values)
            # HC1 robust se
            meat = (X.values * res[:, None] ** 2).T @ X.values
            se = np.sqrt(np.diag(XtX @ meat @ XtX) * len(yy) / (len(yy) - X.shape[1]))
            for k, name in enumerate(X.columns):
                reg.append({"date": d, "group": x.group.iloc[0], "role": x.role.iloc[0], "outcome": f"dlog_{col}",
                            "term": name, "coef": b[k], "se_hc1": se[k], "n": len(yy),
                            "r2": 1 - res.var() / yy.var()})
    pd.DataFrame(rows).to_csv(RES / "p_centrality_rank_stability.csv", index=False)
    pd.DataFrame(reg).to_csv(RES / "p_centrality_gain_ols.csv", index=False)

    # ---- net balance (one-way flows) vs zone attributes on every night
    rows = []
    for (d, subset), x in ((k, v) for k0, sub in nm.groupby("date")
                           for k, v in (((k0, "all"), sub.set_index("idx").join(attrs)),
                                        ((k0, "single_municipality"),
                                         sub.set_index("idx").join(attrs).query("n_municipalities == 1")))):
        for lab, m in (("coastal", x.coastal_true),
                       ("emptied_inland", (x.loss_peak_2021 >= EMPTIED_LOSS) & ~x.coastal_true),
                       ("other_inland_small", (x.loss_peak_2021 < EMPTIED_LOSS) & ~x.coastal_true & (x.population < C.SMALL_VILLAGE_POP)),
                       ("city_gt100k", x.population > 1e5)):
            rows.append({"date": d, "group": x.group.iloc[0], "role": x.role.iloc[0], "subset": subset,
                         "class": lab, "n": int(m.sum()),
                         "median_balance": x.balance[m].median(),
                         "share_net_importer": float((x.balance[m] > 0).mean()),
                         "net_persons": float((x.s_in - x.s_out)[m].sum())})
    pd.DataFrame(rows).to_csv(RES / "p_net_balance_by_class.csv", index=False)


# ---------------------------------------------------------------------------
# stage: cities (rooted vs beach)
# ---------------------------------------------------------------------------
def _baseline(dates: list[str]) -> pd.DataFrame:
    """Mean over the ordinary nights of w per (i, j), a pair missing on a night counting as 0
    (b00.window_mean_from_frame, audit C3). The mean is used, not the median, so that the
    baseline total equals the average ordinary-night volume and excess totals add up (the
    per-pair median of sparse pairs is 0 and would understate the baseline by ~18%)."""
    fr = [load_edges(d, "full").assign(date=d).rename(columns={"i": "r", "j": "s", "w": "p"}) for d in dates]
    m = B.window_mean_from_frame(pd.concat(fr), n_nights=len(dates))
    return m.rename(columns={"r": "i", "s": "j", "p": "w_base"})


@lru_cache(maxsize=4)
def baseline_for_year(year: int) -> pd.DataFrame:
    if year == 2024:
        ds = ordinary_sisters()
    else:
        ds = [d for d in B.available_dates(C.PERNOCT_2022_23)
              if d.startswith(f"{year}-02") and pd.Timestamp(d).day_of_week in (1, 2, 3)]
    return _baseline(ds)


def city_scores(date: str) -> pd.DataFrame:
    nd = nodes()
    R = roots_matrix()
    cities = R[["c", "ine_code", "municipality_name"]].drop_duplicates()
    base = baseline_for_year(int(date[:4]))
    e = load_edges(date, "full")
    e = e[e.i.isin(cities.c)]
    m = e.merge(base[base.i.isin(cities.c)], on=["i", "j"], how="outer").fillna(0.0)
    m["x"] = m.w - m.w_base                      # holiday excess (can be negative)
    m["prov_j"] = nd.prov.values[m.j]
    m["prov_i"] = nd.prov.values[m.i]
    m["coast_j"] = nd.coastal_true.values[m.j]
    m["emptied_j"] = (nd.loss_peak_2021.values[m.j] >= EMPTIED_LOSS) & ~m.coast_j
    m["dkm"] = dist_km(m.i.values, m.j.values, nd)
    rows = []
    born = R.set_index(["c", "prov_birth"]).persons
    for c, g in m.groupby("i"):
        X = g.x.sum()
        pos = g[g.x > 0]
        own = g.prov_i.iloc[0]
        gp = g[g.prov_j != own].groupby("prov_j").x.sum()
        bp = born.loc[c].drop(own, errors="ignore").reindex(gp.index).fillna(0)
        # distance of each province from the city (for a partial correlation)
        pc = prov_centroids().reindex(gp.index)
        dd = np.hypot(nd.x.values[c] - pc.x.values, nd.y.values[c] - pc.y.values) / 1000
        lx = np.log1p(gp.clip(lower=0).values)
        lb = np.log1p(bp.values)
        ld = np.log(dd + 1)
        def res(v):
            Xd = np.column_stack([np.ones_like(ld), ld])
            return v - Xd @ np.linalg.lstsq(Xd, v, rcond=None)[0]
        rows.append({
            "date": date, "c": c, "excess_total": X, "excess_pos_total": pos.x.sum(),
            "base_total": g.w_base.sum(), "night_total": g.w.sum(),
            "ratio_night_base": g.w.sum() / g.w_base.sum(),
            "beach_share": pos.x[pos.coast_j].sum() / pos.x.sum(),
            "emptied_share": pos.x[pos.emptied_j].sum() / pos.x.sum(),
            "own_prov_share": pos.x[pos.prov_j == own].sum() / pos.x.sum(),
            "excess_mean_km": (pos.x * pos.dkm).sum() / pos.x.sum(),
            "roots_spearman": spearman(gp.values, bp.values),
            "roots_partial_r": float(np.corrcoef(res(lx), res(lb))[0, 1]),
            "top_dest_prov": gp.idxmax(),
            "top_birth_prov": bp.idxmax(),
        })
    s = pd.DataFrame(rows).merge(cities, on="c")
    s["population"] = nd.population.values[s.c]
    s["city_coastal"] = nd.coastal_true.values[s.c]
    return s


def stage_cities() -> None:
    hd = hub_dates()
    dates = hd[hd.group.isin(["summer", "allsaints"])].date.tolist() + ["2024-03-29", "2024-12-25", "2024-12-24"]
    S = pd.concat([city_scores(d) for d in dates])
    S.to_csv(RES / "p_cities_scores_by_night.csv", index=False)
    # ranking stability: Spearman of each score between the main summer night and every other night
    main = S[S.date == C.MAIN_NIGHTS["summer"]].set_index("c")
    rows = []
    for d, g in S.groupby("date"):
        g = g.set_index("c").reindex(main.index)
        for col in ("beach_share", "emptied_share", "own_prov_share", "roots_spearman", "roots_partial_r", "excess_mean_km",
                    "ratio_night_base"):
            rows.append({"date": d, "score": col, "spearman_vs_aug15": spearman(g[col], main[col]),
                         "median": g[col].median()})
    pd.DataFrame(rows).to_csv(RES / "p_cities_rank_stability.csv", index=False)
    # summary per city: median over the 9 August sisters + years
    sis = S[S.date.isin(C.SISTER_NIGHTS["summer"])]
    agg = sis.groupby(["c", "municipality_name"]).agg(
        population=("population", "first"), city_coastal=("city_coastal", "first"),
        beach_share_med=("beach_share", "median"), beach_share_min=("beach_share", "min"), beach_share_max=("beach_share", "max"),
        emptied_share_med=("emptied_share", "median"), own_prov_share_med=("own_prov_share", "median"),
        roots_spearman_med=("roots_spearman", "median"), roots_partial_med=("roots_partial_r", "median"),
        excess_km_med=("excess_mean_km", "median"), ratio_med=("ratio_night_base", "median"),
        top_dest_prov=("top_dest_prov", lambda s: s.mode().iloc[0]), top_birth_prov=("top_birth_prov", "first")).reset_index()
    for yd in C.YEAR_REPLICATION["summer"]:
        y = S[S.date == yd].set_index("c")
        agg[f"beach_share_{yd[:4]}"] = agg.c.map(y.beach_share)
        agg[f"roots_partial_{yd[:4]}"] = agg.c.map(y.roots_partial_r)
    nv = S[S.date == C.MAIN_NIGHTS["allsaints"]].set_index("c")
    agg["beach_share_nov1"] = agg.c.map(nv.beach_share)
    agg["emptied_share_nov1"] = agg.c.map(nv.emptied_share)
    agg["roots_partial_nov1"] = agg.c.map(nv.roots_partial_r)
    agg.sort_values("beach_share_med").to_csv(RES / "p_cities_rooted_vs_beach.csv", index=False)
    # what predicts a "beach" city? simple correlations (n = 137, descriptive)
    nd = nodes()
    agg["log_pop"] = np.log(agg.population)
    agg["city_is_coastal"] = agg.city_coastal.astype(float)
    born = roots_matrix()
    own = born[born.own].set_index("c").persons
    tot = born.groupby("c").persons.sum()
    agg["share_born_other_prov"] = agg.c.map(1 - own / tot)
    rows = []
    for col in ("log_pop", "city_is_coastal", "share_born_other_prov"):
        for sc in ("beach_share_med", "emptied_share_med", "roots_partial_med"):
            rows.append({"x": col, "score": sc, "spearman": spearman(agg[col], agg[sc]), "n": len(agg)})
    pd.DataFrame(rows).to_csv(RES / "p_cities_correlates.csv", index=False)


# ---------------------------------------------------------------------------
# stage: holidays (is All Saints a smaller August? Easter, Christmas, weekends)
# ---------------------------------------------------------------------------
HOLI_NIGHTS = {
    "Aug15": "2024-08-15", "Aug14": "2024-08-14", "Aug10": "2024-08-10", "Aug18": "2024-08-18",
    "Aug05": "2024-08-05", "Aug25": "2024-08-25", "Aug29": "2024-08-29",
    "HolyThu28Mar": "2024-03-28", "GoodFri29Mar": "2024-03-29", "Sat30Mar": "2024-03-30", "Sun31Mar": "2024-03-31",
    "Oct31": "2024-10-31", "Nov01": "2024-11-01", "Nov02": "2024-11-02", "Nov03": "2024-11-03",
    "Dec07": "2024-12-07", "Dec08": "2024-12-08",
    "Dec24": "2024-12-24", "Dec25": "2024-12-25", "Dec31": "2024-12-31",
    "Oct12Sat": "2024-10-12", "Sat17Feb": "2024-02-17", "Sun18Feb": "2024-02-18", "Sun20Oct": "2024-10-20",
    "Sun17Nov": "2024-11-17", "Feb13Carnival": "2024-02-13", "Feb28AndDay": "2024-02-28",
    "Aug15_2023": "2023-08-15", "Aug15_2022": "2022-08-15", "Nov01_2022": "2022-11-01",
}


def stage_holidays() -> None:
    nd = nodes()
    n = len(nd)
    base = baseline_for_year(2024)
    sb_in = np.bincount(base.j, base.w_base, n)
    coast = nd.coastal_true.values
    emptied = (nd.loss_peak_2021.values >= EMPTIED_LOSS) & ~coast
    small = (nd.population.values < C.SMALL_VILLAGE_POP) & ~coast
    prov = nd.prov.values
    vec, pairx, feats = {}, {}, []
    for lab, d in HOLI_NIGHTS.items():
        e = load_edges(d, "full")
        b = base if d.startswith("2024") else baseline_for_year(int(d[:4]))
        m = e.merge(b, on=["i", "j"], how="outer").fillna(0.0)
        m["x"] = m.w - m.w_base
        s_in = np.bincount(m.j, m.w, n)
        sbi = np.bincount(m.j, m.w_base, n)
        vec[lab] = s_in - sbi
        pairx[lab] = m.set_index(m.i.astype(np.int64) * n + m.j).x
        pos = m[m.x > 0]
        dk = dist_km(pos.i.values, pos.j.values, nd)
        X = pos.x.values
        feats.append({"label": lab, "date": d, "weekday": pd.Timestamp(d).day_name()[:3],
                      "away": m.w.sum(), "base_away": m.w_base.sum(), "excess_net": m.x.sum(),
                      "excess_pos": X.sum(),
                      "exc_share_coastal": X[coast[pos.j]].sum() / X.sum(),
                      "exc_share_emptied_inland": X[emptied[pos.j]].sum() / X.sum(),
                      "exc_share_small_inland": X[small[pos.j]].sum() / X.sum(),
                      "exc_share_same_prov": X[prov[pos.i] == prov[pos.j]].sum() / X.sum(),
                      "exc_mean_km": (X * dk).sum() / X.sum(),
                      "exc_median_km": wquantile(dk, X, 0.5),
                      "exc_share_ge300km": X[dk >= 300].sum() / X.sum(),
                      **{f"roots_{k}": v for k, v in roots_index(e).items()}})
    F = pd.DataFrame(feats)
    F.to_csv(RES / "p_holiday_features.csv", index=False)
    labs = list(HOLI_NIGHTS)
    V = pd.DataFrame(vec)
    # cosine similarity of destination excess vectors, and pair-level excess correlation
    rows = []
    for a in labs:
        for b in labs:
            va, vb = V[a].values, V[b].values
            cos = va.dot(vb) / np.linalg.norm(va) / np.linalg.norm(vb)
            # scale + fit: vb ~ alpha * va
            alpha = va.dot(vb) / va.dot(va)
            r2 = 1 - ((vb - alpha * va) ** 2).sum() / ((vb - vb.mean()) ** 2).sum()
            pa, pb = pairx[a], pairx[b]
            u = pa.index.union(pb.index)
            pa_, pb_ = pa.reindex(u).fillna(0).values, pb.reindex(u).fillna(0).values
            rows.append({"a": a, "b": b, "cos_dest_excess": cos, "alpha_b_on_a": alpha, "r2_b_on_a": r2,
                         "spearman_dest_excess": spearman(va, vb),
                         "pearson_pair_excess": float(np.corrcoef(pa_, pb_)[0, 1])})
    pd.DataFrame(rows).to_csv(RES / "p_holiday_similarity.csv", index=False)
    # where is All Saints different from August? residual of Nov01 on Aug15 by zone class
    a, b = V["Aug15"].values, V["Nov01"].values
    alpha = a.dot(b) / a.dot(a)
    res = b - alpha * a
    rows = []
    for lab, msk in (("coastal", coast), ("emptied_inland", emptied), ("small_inland", small),
                     ("city_gt100k", nd.population.values > 1e5), ("all", np.ones(n, bool))):
        rows.append({"class": lab, "n": int(msk.sum()), "alpha": alpha,
                     "aug15_excess": a[msk].sum(), "nov01_excess": b[msk].sum(),
                     "nov01_over_aug15": b[msk].sum() / a[msk].sum(),
                     "residual_sum": res[msk].sum()})
    pd.DataFrame(rows).to_csv(RES / "p_holiday_allsaints_vs_aug.csv", index=False)
    V.insert(0, "id", nd.id.values)
    V.to_parquet(RES / "p_holiday_dest_excess.parquet", index=False)


# ---------------------------------------------------------------------------
# stage: distance distributions
# ---------------------------------------------------------------------------
DIST_BINS = np.array([0, 10, 25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 1200])


def stage_distance() -> None:
    nd = nodes()
    base = baseline_for_year(2024)
    rows = []
    nights = {"Feb14": "2024-02-14", "Sat17Feb": "2024-02-17", "Aug15": "2024-08-15", "Aug05": "2024-08-05",
              "Nov01": "2024-11-01", "Nov02": "2024-11-02", "GoodFri": "2024-03-29", "Dec25": "2024-12-25",
              "Aug15_2023": "2023-08-15", "Aug15_2022": "2022-08-15"}
    coast = nd.coastal_true.values
    for lab, d in nights.items():
        e = load_edges(d, "full")
        b = base if d.startswith("2024") else baseline_for_year(int(d[:4]))
        m = e.merge(b, on=["i", "j"], how="outer").fillna(0.0)
        m["x"] = m.w - m.w_base
        m["dk"] = dist_km(m.i.values, m.j.values, nd)
        m["bin"] = pd.cut(m.dk, DIST_BINS, right=False)
        m["coast_j"] = coast[m.j]
        eb = B.backbone_threshold(e, W_MIN)
        eb = eb.assign(dk=dist_km(eb.i.values, eb.j.values, nd))
        eb["bin"] = pd.cut(eb.dk, DIST_BINS, right=False)
        g = m.groupby("bin", observed=False)
        t = pd.DataFrame({"persons": g.w.sum(), "base_persons": g.w_base.sum(), "excess": g.x.sum(),
                          "excess_to_coast": m[m.coast_j].groupby("bin", observed=False).x.sum(),
                          "excess_to_inland": m[~m.coast_j].groupby("bin", observed=False).x.sum(),
                          "w10_links": eb.groupby("bin", observed=False).size()})
        t["share_persons"] = t.persons / t.persons.sum()
        t["share_excess"] = t.excess / t.excess.sum()
        t["ratio_vs_base"] = t.persons / t.base_persons
        t = t.reset_index().assign(label=lab, date=d)
        t["bin"] = t["bin"].astype(str)
        rows.append(t)
    pd.concat(rows).to_csv(RES / "p_distance_bins.csv", index=False)


# ---------------------------------------------------------------------------
# stage: plots
# ---------------------------------------------------------------------------
def stage_plots() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    DIAG.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(RES / "p_daily_metrics.csv")
    d = d[d.date.str.startswith("2024")].copy()
    d["t"] = pd.to_datetime(d.date)
    cols = ["away_persons", "w10_density_und", "w10_recip", "w10_r_und", "eq_r_und", "eq_recip",
            "share_dest_coastal", "share_dest_emptied_inland", "share_same_prov", "dist_mean_km",
            "roots_corr_partial", "w10_spearman_kin_kout", "w10_kin_max", "w10_kout_max", "eq_C_avg", "share_downhill"]
    months = [2, 3, 8, 10, 11, 12]
    fig, axes = plt.subplots(len(cols), len(months), figsize=(22, 2.0 * len(cols)), sharey="row")
    colors = {"Sat": "tab:orange", "Sun": "tab:red", "holiday": "black"}
    for r, col in enumerate(cols):
        for c, mth in enumerate(months):
            ax = axes[r, c]
            s = d[d.t.dt.month == mth]
            ax.plot(s.t.dt.day, s[col], "-", color="0.6", lw=1)
            for typ, colr in colors.items():
                q = s[s.day_type == typ]
                ax.plot(q.t.dt.day, q[col], "o", color=colr, ms=3)
            if r == 0:
                ax.set_title({2: "Feb", 3: "Mar", 8: "Aug", 10: "Oct", 11: "Nov", 12: "Dec"}[mth])
            if c == 0:
                ax.set_ylabel(col, fontsize=7)
            ax.tick_params(labelsize=6)
    fig.suptitle("Daily metrics 2024 (orange = date is Sat, red = Sun, black = national holiday)")
    fig.tight_layout()
    fig.savefig(DIAG / "p_daily_timeseries.png", dpi=110)
    plt.close(fig)

    # metric vs volume (is a holiday just a bigger weekend?)
    fig, axes = plt.subplots(2, 4, figsize=(18, 8))
    pcol = {"Feb": "0.6", "Mar": "0.4", "Easter week": "tab:green", "August": "tab:blue", "Oct": "0.5",
            "All Saints": "tab:purple", "Nov": "0.7", "Dec": "0.3", "Christmas": "tab:red"}
    for ax, col in zip(axes.flat, ["w10_r_und", "eq_r_und", "w10_recip", "eq_recip", "share_dest_coastal",
                                   "share_dest_emptied_inland", "roots_corr_partial", "dist_mean_km"]):
        for p, colr in pcol.items():
            s = d[d.period == p]
            ax.scatter(s.away_persons / 1e6, s[col], s=12, color=colr, label=p)
        ax.set_xlabel("away persons (M)")
        ax.set_title(col, fontsize=9)
    axes[0, 0].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(DIAG / "p_metric_vs_volume.png", dpi=110)
    plt.close(fig)

    # distance distributions
    t = pd.read_csv(RES / "p_distance_bins.csv")
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
    for lab, g in t.groupby("label", sort=False):
        axes[0].plot(range(len(g)), g.share_persons, "-o", ms=3, label=lab)
        if lab != "Feb14":
            axes[1].plot(range(len(g)), g.ratio_vs_base, "-o", ms=3, label=lab)
    for ax in axes:
        ax.set_xticks(range(len(DIST_BINS) - 1))
        ax.set_xticklabels([f"{a}-{b}" for a, b in zip(DIST_BINS[:-1], DIST_BINS[1:])], rotation=60, fontsize=7)
    axes[0].set_title("share of away persons by distance (km)")
    axes[1].set_title("persons / Feb baseline (mean of 13 nights)")
    axes[1].set_yscale("log")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(DIAG / "p_distance.png", dpi=110)
    plt.close(fig)

    # degree CCDFs in vs out, Feb vs Aug vs Nov
    nm = pd.read_parquet(RES / "p_node_metrics.parquet")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, bb in zip(axes, ("w10", "equal_edges")):
        for g, dte in C.MAIN_NIGHTS.items():
            s = nm[nm.date == dte]
            for col, ls in ((f"kin_{bb}", "-"), (f"kout_{bb}", "--")):
                k = np.sort(s[col].values)
                k = k[k > 0]
                ax.loglog(k, 1 - np.arange(len(k)) / len(k), ls, label=f"{g} {col.split('_')[0]}")
        ax.set_title(f"CCDF k_in (solid) / k_out (dashed), {bb}")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(DIAG / "p_degree_ccdf_in_out.png", dpi=110)
    plt.close(fig)

    # holiday similarity heatmap
    S = pd.read_csv(RES / "p_holiday_similarity.csv")
    M = S.pivot(index="a", columns="b", values="cos_dest_excess").loc[list(HOLI_NIGHTS), list(HOLI_NIGHTS)]
    fig, ax = plt.subplots(figsize=(11, 10))
    im = ax.imshow(M.values, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(M)))
    ax.set_xticklabels(M.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(len(M)))
    ax.set_yticklabels(M.index, fontsize=7)
    for a in range(len(M)):
        for b in range(len(M)):
            ax.text(b, a, f"{M.values[a, b]:.2f}", ha="center", va="center", fontsize=5)
    fig.colorbar(im)
    ax.set_title("cosine similarity of destination excess (night - Feb mean)")
    fig.tight_layout()
    fig.savefig(DIAG / "p_holiday_similarity.png", dpi=110)
    plt.close(fig)


STRUCT_METRICS = ["eq_r_und", "w10_r_und", "eq_recip", "w10_recip", "eq_C_avg", "w10_spearman_kin_kout",
                  "eq_gini_kout", "eq_gini_kin", "recip_weighted_full", "share_dest_coastal", "share_dest_emptied_inland",
                  "share_same_prov", "share_downhill", "roots_corr_partial", "dist_mean_km"]


def _ols_r2(X: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray]:
    X = np.column_stack([np.ones(len(y)), X])
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    r = y - X @ b
    return 1 - r.var() / y.var(), b


def stage_panel() -> None:
    d = pd.read_csv(RES / "p_daily_metrics.csv")
    d24 = d[d.date.str.startswith("2024")].copy()

    # ---- 1. date convention: weekday profile of ordinary weeks and holiday eves vs holiday nights
    ordm = d24[d24.period.isin(["Feb", "Oct", "Nov"]) & (d24.day_type != "holiday")]
    wk = ordm.groupby("weekday").agg(n=("date", "size"), away_median=("away_persons", "median"),
                                     dist_mean_median=("dist_mean_km", "median"),
                                     share_same_prov_median=("share_same_prov", "median"))
    wk = wk.reindex(["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]).reset_index()
    wk.insert(0, "table", "weekday_profile_ordinary_Feb_Oct_Nov_2024")
    probe = ["2022-10-31", "2022-11-01", "2022-11-02", "2024-10-31", "2024-11-01", "2024-11-02", "2024-11-03",
             "2024-11-04", "2024-12-23", "2024-12-24", "2024-12-25", "2024-12-26", "2024-03-27", "2024-03-28",
             "2024-03-29", "2024-03-30", "2024-03-31", "2024-08-14", "2024-08-15", "2024-08-16",
             "2024-12-05", "2024-12-06", "2024-12-07", "2024-12-08", "2024-12-09", "2024-12-10"]
    pr = d[d.date.isin(probe)][["date", "weekday", "away_persons", "dist_mean_km", "share_same_prov"]].copy()
    pr.insert(0, "table", "holiday_probe")
    pr["calendar_label"] = pr.date.map(HOLIDAYS_2024).fillna("")
    pd.concat([wk, pr]).to_csv(RES / "p_date_convention.csv", index=False)

    # ---- 2. is the network shape a function of volume, or of trip length? (180 nights of 2024)
    rows = []
    la = np.log(d24.away_persons.values)
    dm = d24.dist_mean_km.values
    coast = d24.share_dest_coastal.values
    for m in STRUCT_METRICS:
        y = d24[m].values
        r2v, bv = _ols_r2(la[:, None], y)
        r2d, _ = _ols_r2(dm[:, None], y)
        r2vd, bvd = _ols_r2(np.column_stack([la, dm]), y)
        r2vdc, _ = _ols_r2(np.column_stack([la, dm, coast]), y)
        rows.append({"metric": m, "n_nights": len(y), "r2_log_volume": r2v, "r2_dist_mean": r2d,
                     "r2_volume_plus_dist": r2vd, "r2_volume_dist_coastshare": r2vdc,
                     "pearson_with_log_volume": float(np.corrcoef(la, y)[0, 1]),
                     "pearson_with_dist_mean": float(np.corrcoef(dm, y)[0, 1])})
    pd.DataFrame(rows).to_csv(RES / "p_volume_vs_structure.csv", index=False)

    # ---- 3. matched-volume contrast: weekend nights vs holiday nights with similar volume
    sun = d24[(d24.weekday == "Sun") & d24.period.isin(["Feb", "Mar", "Oct", "Nov", "Dec"])]
    hol = d24[d24.period.isin(["Easter week", "All Saints", "Christmas"]) | d24.date.isin(["2024-12-07", "2024-12-08"])]
    rows = []
    for _, h in hol.iterrows():
        k = (sun.away_persons - h.away_persons).abs().idxmin()
        s = sun.loc[k]
        r = {"holiday_date": h.date, "holiday_period": h.period, "holiday_away": h.away_persons,
             "matched_sunday": s.date, "sunday_away": s.away_persons}
        for m in STRUCT_METRICS:
            r[f"{m}__hol"] = h[m]
            r[f"{m}__sun"] = s[m]
        rows.append(r)
    pd.DataFrame(rows).to_csv(RES / "p_weekend_vs_holiday_matched.csv", index=False)

    # ---- 4. switch: how fast does each period move away from the ordinary range?
    feb = d24[d24.date.isin(ordinary_sisters())]
    windows = {"Easter": ("2024-03-11", "2024-03-31"), "August": ("2024-08-01", "2024-08-31"),
               "AllSaints": ("2024-10-21", "2024-11-08"), "Christmas": ("2024-12-10", "2024-12-31"),
               "DecPuente": ("2024-12-02", "2024-12-12")}
    rows = []
    for w, (a, b) in windows.items():
        s = d24[(d24.date >= a) & (d24.date <= b)].sort_values("date")
        for m in STRUCT_METRICS + ["away_persons"]:
            lo, hi, med = feb[m].min(), feb[m].max(), feb[m].median()
            sd = feb[m].std()
            y = s[m].values
            peak_i = int(np.argmax(np.abs(y - med)))
            peak = y[peak_i]
            frac = (y - med) / (peak - med) if peak != med else np.full(len(y), np.nan)
            outside = (y < lo) | (y > hi)
            jumps = np.abs(np.diff(y))
            for k, (dt, yy, fr, out) in enumerate(zip(s.date, y, frac, outside)):
                rows.append({"window": w, "metric": m, "date": dt, "weekday": s.weekday.iloc[k], "value": yy,
                             "ord_median": med, "ord_min": lo, "ord_max": hi, "z_vs_ord": (yy - med) / sd,
                             "frac_of_peak_shift": fr, "outside_ord_range": bool(out),
                             "peak_date": s.date.iloc[peak_i],
                             "max_daily_jump_share": float(jumps.max() / abs(peak - med)) if len(jumps) else np.nan})
    pd.DataFrame(rows).to_csv(RES / "p_switch.csv", index=False)

    # ---- 5. one curve? fit metric ~ log(volume) on NON-holiday nights of 2024 (Feb, Mar before Easter
    # week, Oct, Nov and Dec outside the holiday windows; weekdays AND weekends), then residuals of
    # every night (2022-23 included) from that curve, summarised by period. Extrapolation above
    # ~5.1 M away persons (the largest ordinary Sunday) is flagged.
    d["period"] = [period(x) for x in d.date]
    hol_dates = {"2024-12-06", "2024-12-07", "2024-12-08", "2024-12-09"}
    fit = d[d.date.str.startswith("2024") & d.period.isin(["Feb", "Mar", "Oct", "Nov", "Dec"]) & ~d.date.isin(hol_dates)]
    vmax = fit.away_persons.max()
    rows = []
    for m in ["share_dest_emptied_inland", "share_dest_small_inland", "share_dest_coastal", "eq_r_und", "eq_recip",
              "share_same_prov", "roots_corr_partial", "dist_mean_km"]:
        r2, b = _ols_r2(np.log(fit.away_persons.values)[:, None], fit[m].values)
        pred = b[0] + b[1] * np.log(d.away_persons.values)
        res = d[m].values - pred
        sd = (fit[m].values - (b[0] + b[1] * np.log(fit.away_persons.values))).std()
        for p, g in d.assign(res=res, extrap=d.away_persons > vmax).groupby("period"):
            rows.append({"metric": m, "fit_r2": r2, "fit_slope_per_log_volume": b[1], "fit_resid_sd": sd,
                         "period": p, "n_nights": len(g), "median_resid": g.res.median(),
                         "median_resid_in_sd": g.res.median() / sd, "share_extrapolated": g.extrap.mean()})
    pd.DataFrame(rows).to_csv(RES / "p_one_curve_residuals.csv", index=False)


STAGES = {"daily": stage_daily, "zones": stage_zones, "hubs": stage_hubs, "cities": stage_cities,
          "holidays": stage_holidays, "distance": stage_distance, "panel": stage_panel, "plots": stage_plots}

if __name__ == "__main__":
    todo = sys.argv[1:] or list(STAGES)
    for s in todo:
        print(f"[p_patterns] stage {s}", flush=True)
        STAGES[s]()

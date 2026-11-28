"""Inputs and survey household contributions for Analysis C."""
from pathlib import Path
import os
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import sparse
from shapely import STRtree

HERE = Path(__file__).resolve().parent
# One repository root holds code and data. SHR_MAIN_ROOT may point the coursework inputs
# (data/, raw_data_cache/, clean_data/, _local/external/, results/) at another checkout.
JOURNAL = HERE.parents[2]
MAIN = Path(os.environ.get("SHR_MAIN_ROOT", JOURNAL)).resolve()
GISCO = JOURNAL / "_local/external_journal/gisco_lau_2021"
PRIVATE = JOURNAL / "_local/journal_data/C_municipal"
OUT = HERE / "out"
MONTHS = pd.period_range("2019-07", "2026-09", freq="M").astype(str).tolist()
sys.path.insert(0, str(MAIN / "pipeline"))
import n1_new_data as n1

PARENTS = ["MUNINP1", "MUNINP2"]
THRESHOLD = 30


def municipal_geometry():
    lau = gpd.read_file(GISCO / "LAU_RG_01M_2021_3035.geojson", where="CNTR_CODE = 'ES'")
    lau = lau.assign(code=lau["LAU_ID"].str.zfill(5)).set_index("code", verify_integrity=True)
    shore = gpd.read_file(GISCO / "NUTS_BN_01M_2021_3035_LEVL_0.geojson", where="COAS_FLAG = 'T'")
    assert lau.crs.to_epsg() == shore.crs.to_epsg() == 3035
    coast_pairs = STRtree(shore.geometry.array).query(lau.geometry.array, predicate="intersects")
    centres = lau.centroid
    geo = pd.DataFrame(dict(x=centres.x, y=centres.y))
    geo["coastal"] = False
    geo.iloc[np.unique(coast_pairs[0]), geo.columns.get_loc("coastal")] = True
    pop = pd.read_csv(n1.LAU_CSV, dtype={"LAU_ID": str}).set_index("LAU_ID")["POP_2021"]
    pop.index = pop.index.str.zfill(5)
    assert set(pop.index) == set(geo.index)
    np.testing.assert_equal(pop.reindex(lau.index).to_numpy(), lau["POP_2021"].to_numpy())
    return geo.join(pop.rename("pop_2021"))


def fixed_frame():
    geo = municipal_geometry()
    cities = pd.read_csv(MAIN / "results/r2_city_set.csv", dtype={"ine_code": str})["ine_code"]
    assert len(cities) == cities.nunique() == 137
    mainland = geo.loc[~geo.index.str[:2].isin(n1.EXCL)]
    villages = mainland.loc[mainland["pop_2021"].between(0, 2000)]
    frame = pd.MultiIndex.from_product([sorted(cities), villages.index.sort_values()],
                                       names=["orig", "dest"]).to_frame(index=False)
    frame["oprov"], frame["dprov"] = frame.orig.str[:2], frame.dest.str[:2]
    frame = frame.loc[frame.oprov.ne(frame.dprov)].reset_index(drop=True)
    frame["corridor"] = frame.orig + frame.dprov
    ox, oy = frame.orig.map(geo.x), frame.orig.map(geo.y)
    dx, dy = frame.dest.map(geo.x), frame.dest.map(geo.y)
    frame["km"] = np.hypot(ox - dx, oy - dy) / 1000
    frame["log_km"] = np.log(frame.km.clip(lower=1))
    frame["coastal"] = frame.dest.map(geo.coastal)
    frame["pop"] = frame.dest.map(geo.pop_2021)

    house = pd.read_csv(n1.VIV_RAW / "censo2021_59525.csv", sep=";", dtype=str)
    house["code"] = house["Municipios"].str.extract(r"^(\d{5})\s")
    house["count"] = house["Total"].map(n1.parse_es_num)
    housing = house.dropna(subset=["code"]).pivot(
        index="code", columns="Tipo de vivienda (principal o no)", values="count")
    share = housing["Vivienda no principal"] / housing["Total"]
    frame["housing_share"] = frame.dest.map(share)
    assert frame[["km", "coastal", "pop"]].notna().all().all()
    assert frame.index.is_unique
    geo.to_parquet(PRIVATE / "municipal_geometry.parquet")
    frame.to_parquet(PRIVATE / "fixed_frame.parquet", index=False)
    return frame, geo.index


class Survey:
    """Rebuild survey-weighted exposures by retaining whole households."""

    def __init__(self, frame, geographic_codes):
        columns = ["IDEN", "IDQ_MUN", "MUNIN", *PARENTS, "FACTOR"]
        people = pd.read_parquet(n1.ECEPOV_DIR / "ecepov_hogar_2021.parquet", columns=columns)
        people = people.rename(columns={"IDQ_MUN": "orig"})
        cities = pd.Index(frame.orig.unique())
        people = people.loc[people.orig.isin(cities)].copy()
        people["household"] = people.groupby(["orig", "IDEN"], sort=True).ngroup()
        people["city_id"] = cities.get_indexer(people.orig)
        people["city_born"] = people.orig.eq(people.MUNIN)
        assert people.FACTOR.gt(0).all()
        self.city = frame.orig.map(dict(zip(cities, range(len(cities))))).to_numpy()
        hh = people.groupby("household", sort=True)
        self.house_city = hh.city_id.first().to_numpy()
        self.total = hh.FACTOR.sum().to_numpy()
        self.local = (people.FACTOR * people.city_born).groupby(people.household).sum().to_numpy()
        self.hh_by_city = [np.flatnonzero(self.house_city == c) for c in range(len(cities))]

        local = people.loc[people.city_born].copy()
        self.completeness(people, geographic_codes, "all_residents")
        self.completeness(local, geographic_codes, "city_born")
        keys = pd.MultiIndex.from_frame(frame[["orig", "dest"]])
        g1 = people[["orig", "household", "FACTOR", "MUNIN"]].rename(columns={"MUNIN": "dest"})
        g2 = local.melt(id_vars=["orig", "household", "FACTOR"],
                        value_vars=PARENTS, value_name="dest")
        g2["FACTOR"] /= 2
        self.g1, self.g1_support, self.g1_rows = self.contributions(g1, keys)
        self.g2, self.g2_support, self.g2_rows = self.contributions(g2, keys)
        self.n = len(frame)

    def completeness(self, people, geographic_codes, scope):
        people = people.copy()
        people["slots_blank"] = people[PARENTS].apply(lambda x: x.str.strip().eq("")).sum(axis=1)
        people["slots_foreign"] = people[PARENTS].eq("99999").sum(axis=1)
        people["slots_matched_2021"] = people[PARENTS].isin(geographic_codes).sum(axis=1)
        people["slots_other_code"] = (
            2 - people[["slots_blank", "slots_foreign", "slots_matched_2021"]].sum(axis=1))
        people["both_reported"] = people.slots_blank.eq(0)
        people["one_reported"] = people.slots_blank.eq(1)
        people["none_reported"] = people.slots_blank.eq(2)
        columns = ["slots_blank", "slots_foreign", "slots_matched_2021", "slots_other_code",
                   "both_reported", "one_reported", "none_reported"]
        raw = people.groupby("orig").agg(persons=("FACTOR", "size"),
                                          households=("household", "nunique"),
                                          weighted_persons=("FACTOR", "sum"))
        unweighted = people.groupby("orig")[columns].sum().add_suffix("_unweighted")
        weighted = people[columns].multiply(people.FACTOR, axis=0).groupby(people.orig).sum()
        weighted["blank_slot_rate"] = weighted.slots_blank / (2 * raw.weighted_persons)
        weighted["both_reported_rate"] = weighted.both_reported / raw.weighted_persons
        raw.join(unweighted).join(weighted).to_csv(OUT / f"parent_completeness_{scope}.csv")

    def contributions(self, ties, keys):
        rows = keys.get_indexer(pd.MultiIndex.from_frame(ties[["orig", "dest"]]))
        keep = rows >= 0
        rows = rows[keep]
        columns = ties.household.to_numpy()[keep]
        values = ties.FACTOR.to_numpy()[keep]
        active, compact = np.unique(rows, return_inverse=True)
        weighted = sparse.csr_matrix((values, (compact, columns)), shape=(len(active), len(self.total)))
        support = weighted.copy()
        support.data[:] = 1
        return weighted, support, active

    def draw(self, rng):
        counts = np.zeros(len(self.total), dtype=int)
        for hh in self.hh_by_city:
            counts[hh] = np.bincount(rng.integers(len(hh), size=len(hh)), minlength=len(hh))
        return counts

    def exposures(self, multiplicity):
        den1 = np.bincount(self.house_city, weights=self.total * multiplicity)
        den2 = np.bincount(self.house_city, weights=self.local * multiplicity)
        g1 = np.zeros(self.n)
        g2 = np.zeros(self.n)
        h1 = np.zeros(self.n, dtype=int)
        h2 = np.zeros(self.n, dtype=int)
        g1[self.g1_rows] = (self.g1 @ multiplicity) / den1[self.city[self.g1_rows]]
        g2[self.g2_rows] = (self.g2 @ multiplicity) / den2[self.city[self.g2_rows]]
        h1[self.g1_rows] = self.g1_support @ (multiplicity > 0).astype(int)
        h2[self.g2_rows] = self.g2_support @ (multiplicity > 0).astype(int)
        return g1, g2, h1, h2


def od_year(year):
    # Same curated workbook inputs and schema as n1.stage_roots_mun.
    flow = pd.read_parquet(n1.INE_OD_DIR / f"ine_od_mun_{year}.parquet",
                           columns=["orig", "dest", "mes", "turistas"])
    assert not flow.duplicated(["orig", "dest", "mes"]).any()
    assert flow.orig.str[:2].ne(flow.dest.str[:2]).all()
    assert flow.turistas.ge(THRESHOLD).all()
    assert flow.orig.str.fullmatch(r"\d{5}").all() and flow.dest.str.fullmatch(r"\d{5}").all()
    return flow


def panel(frame):
    keys = pd.MultiIndex.from_frame(frame[["orig", "dest"]])
    months = {}
    for year in range(2019, 2027):
        for month, group in od_year(year).groupby("mes", sort=True):
            series = group.set_index(["orig", "dest"])["turistas"]
            months[month] = series.reindex(keys, fill_value=-1).to_numpy(dtype=np.int32)
    assert list(months) == MONTHS
    return months

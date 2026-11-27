"""Shared inputs for the journal design pilots."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

MAIN = Path("/Users/alejandrotreny/Proyectos/spain-holiday-roots")
EXTERNAL = MAIN / "_local/external"
OUT = Path(__file__).resolve().parent / "out"
sys.path.insert(0, str(MAIN / "pipeline"))
import b00_common as nodes_and_nights


def municipalities(max_population=2000):
    """Use 2021 municipal population and the primary MITMA zone for geometry."""
    b0 = nodes_and_nights
    nodes = b0.load_nodes().set_index("id")
    crosswalk = pd.read_csv(b0.INE_MITMA_REL, sep="|", dtype=str)
    primary = (
        crosswalk.groupby(["municipio_ine", "municipio_mitma"]).size().rename("sections")
        .reset_index().sort_values(["sections", "municipio_mitma"], ascending=[False, True])
        .drop_duplicates("municipio_ine").set_index("municipio_ine")["municipio_mitma"]
    )
    lau = pd.read_csv(EXTERNAL / "eurostat_lau_2021/es_lau_2021.csv", dtype={"LAU_ID": str})
    lau["LAU_ID"] = lau["LAU_ID"].str.zfill(5)
    assert lau["LAU_ID"].is_unique
    villages = (
        lau.loc[lau["POP_2021"].between(1, max_population)]
        .set_index("LAU_ID").join(primary.rename("zone")).join(nodes, on="zone")
        .dropna(subset=["x", "y"]).copy()
    )
    villages["dprov"] = villages.index.str[:2]
    coasts = nodes.loc[nodes["coastal_true"], ["x", "y"]]
    villages["coast_km"] = cKDTree(coasts).query(villages[["x", "y"]])[0] / 1000
    cities = pd.read_csv(MAIN / "results/r2_city_set.csv", dtype={"ine_code": str})
    cities = cities.set_index("ine_code").join(nodes[["x", "y"]], how="inner")
    return cities, villages


def birth_ties(cities, villages):
    """Count weighted ties and household support for city-born residents' parents."""
    columns = ["IDQ_MUN", "MUNIN", "MUNINP1", "MUNINP2", "FACTOR", "IDEN"]
    people = pd.read_parquet(EXTERNAL / "ine_ecepov_2021/ecepov_hogar_2021.parquet", columns=columns)
    people = people.loc[people["IDQ_MUN"].isin(cities)].rename(columns={"IDQ_MUN": "orig"})
    own_birth = people.loc[people["MUNIN"].isin(villages)]
    first = own_birth.groupby(["orig", "MUNIN"])["FACTOR"].sum()
    first.index.names = ["orig", "dest"]
    city_born = people.loc[people["orig"].eq(people["MUNIN"])]
    parents = city_born.melt(
        id_vars=["orig", "IDEN", "FACTOR"], value_vars=["MUNINP1", "MUNINP2"],
        value_name="dest",
    )
    parents = parents.loc[parents["dest"].isin(villages)].copy()
    parents["parent_weight"] = parents["FACTOR"] / 2
    second = parents.groupby(["orig", "dest"]).agg(
        g2=("parent_weight", "sum"), g2_households=("IDEN", "nunique"),
    )
    return second.join(first.rename("g1"), how="outer")


def municipal_dyads():
    cities, villages = municipalities()
    base = pd.MultiIndex.from_product([cities.index, villages.index], names=["orig", "dest"]).to_frame(index=False)
    base["dprov"] = base["dest"].str[:2]
    base["oprov"] = base["orig"].str[:2]
    base = base.loc[base["oprov"].ne(base["dprov"])].copy()
    base["city_province"] = base["orig"] + "_" + base["dprov"]
    dx = base["orig"].map(cities["x"]) - base["dest"].map(villages["x"])
    dy = base["orig"].map(cities["y"]) - base["dest"].map(villages["y"])
    base["log_km"] = np.log(np.maximum(np.hypot(dx, dy) / 1000, 1))
    base["log_km_sq"] = base["log_km"] ** 2
    base["dest_pop"] = base["dest"].map(villages["POP_2021"])
    base["orig_pop"] = base["orig"].map(cities["residents_2024"])
    base["income"] = base["orig"].map(cities["income_2024"])
    base["coast_km"] = base["dest"].map(villages["coast_km"])
    base["historical_loss"] = base["dest"].map(villages["loss_peak_2021"])
    base = base.set_index(["orig", "dest"]).join(birth_ties(cities.index, villages.index))
    for column in ["g1", "g2", "g2_households"]:
        base[column] = base[column].fillna(0)
    base["parents"] = base["g2_households"].gt(0).astype(float)
    base["first"] = base["g1"].gt(0).astype(float)
    return base.reset_index()


def ine_counts(base):
    """Missing values denote unpublished flows. Preserve them for estimation."""
    dyads = pd.MultiIndex.from_frame(base[["orig", "dest"]])
    annual = {}
    for year in range(2019, 2027):
        flow = pd.read_parquet(
            EXTERNAL / f"ine_turismo_moviles_od/ine_od_mun_{year}.parquet",
            columns=["orig", "dest", "mes", "turistas"],
        )
        assert flow["turistas"].ge(30).all()
        for month, rows in flow.groupby("mes", sort=True):
            annual[str(month)] = rows.set_index(["orig", "dest"])["turistas"].reindex(dyads).to_numpy(dtype="float32")
    return annual

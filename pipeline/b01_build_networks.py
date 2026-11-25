"""b01_build_networks: build every overnight network named in config and both backbones.

Outputs (results/b01_networks/):
  nodes.parquet                     node table (row position == node index)
  edges_<date>__full.parquet        all mainland directed edges i->j (i != j), weight w = persons
  edges_<date>__<backbone>.parquet  binary backbones (columns i, j, w)
  networks.csv                      one row per (date, night type, role)
  manifest.csv                      one row per (date, backbone)

Backbones (config):
  w10          : BACKBONE_MAIN, directed edges with w >= w_min
  equal_edges  : BACKBONE_CHECK, the E heaviest directed edges, E = E_dir of the
                 w10 backbone of the ordinary MAIN night (2024-02-14), used for every date.
Density = E / (N (N - 1)) for the directed graph; density_und = 2 E_und / (N (N - 1))
for the undirected projection (edge if either direction exists). The MAX_DENSITY check
is applied to density_und, which is always >= the directed density.

Run: .venv/bin/python -I pipeline/b01_build_networks.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import b00_common as B  # noqa: E402

OUT = C.RESULTS / "b01_networks"
BACKBONES = ("w10", "equal_edges")


def resolve_ordinary_sisters() -> list[str]:
    """All Tue/Wed/Thu of Feb 2024 present in the 2024 pernoct file."""
    dates = [d for d in B.available_dates(C.PERNOCT_2024) if d.startswith("2024-02-")]
    keep = [d for d in dates if pd.Timestamp(d).day_of_week in (1, 2, 3)]  # Mon=0
    return sorted(keep)


def network_table() -> pd.DataFrame:
    """Rows (date, group, role) for every network in config. A date can have several roles."""
    sisters = dict(C.SISTER_NIGHTS)
    if isinstance(sisters["ordinary"], str):
        sisters["ordinary"] = resolve_ordinary_sisters()
    rows = []
    for g, d in C.MAIN_NIGHTS.items():
        rows.append((d, g, "main"))
    for g, ds in sisters.items():
        for d in ds:
            rows.append((d, g, "sister"))
    for g, ds in C.YEAR_REPLICATION.items():
        for d in ds:
            rows.append((d, g, "year"))
    for g, d in getattr(C, "APPENDIX_NIGHTS", {}).items():
        rows.append((d, g, "appendix"))
    t = pd.DataFrame(rows, columns=["date", "group", "role"])
    for g, d in C.MAIN_NIGHTS.items():
        assert d in set(t[(t.group == g) & (t.role == "sister")].date), f"main {d} not in its sister list"
    return t


def und_edge_count(e: pd.DataFrame) -> int:
    if e.empty:
        return 0
    a = np.minimum(e.i.values, e.j.values).astype(np.int64)
    b = np.maximum(e.i.values, e.j.values).astype(np.int64)
    return int(np.unique(a * 100000 + b).size)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    nodes = B.load_nodes()
    N = len(nodes)
    nodes.to_parquet(OUT / "nodes.parquet", index=False)
    assert N == C.EXPECTED_N_NODES

    nets = network_table()
    nets.to_csv(OUT / "networks.csv", index=False)
    dates = sorted(nets.date.unique())
    print(f"N={N}; {len(dates)} distinct dates; ordinary sisters: {resolve_ordinary_sisters()}")

    # Check every date exists before building anything.
    avail = {p: set(B.available_dates(p)) for p in (C.PERNOCT_2024, C.PERNOCT_2022_23)}
    missing = [d for d in dates if d not in avail[B.pernoct_file(d)]]
    assert not missing, f"dates missing from data: {missing}"

    full, info = {}, {}
    for d in dates:
        full[d], info[d] = B.build_edges(d, nodes)
        full[d].to_parquet(OUT / f"edges_{d}__full.parquet", index=False)

    ref_date = C.MAIN_NIGHTS[C.BACKBONE_CHECK["reference"]]
    E_ref = len(B.backbone_threshold(full[ref_date], C.BACKBONE_MAIN["w_min"]))
    print(f"equal_edges reference: {ref_date}, E_dir = {E_ref}")

    rows = []
    for d in dates:
        e = full[d]
        tot = float(e.w.sum())
        for bb in BACKBONES:
            b = B.backbone_threshold(e, C.BACKBONE_MAIN["w_min"]) if bb == "w10" else B.backbone_equal_edges(e, E_ref)
            if bb == "equal_edges":
                assert len(b) == min(E_ref, len(e))
            b.to_parquet(OUT / f"edges_{d}__{bb}.parquet", index=False)
            E = len(b)
            E_und = und_edge_count(b)
            g = nets[nets.date == d]
            rows.append(dict(
                label=d, date=d,
                groups=";".join(sorted(g.group.unique())),
                roles=";".join(sorted(g.role.unique())),
                backbone=bb, N=N, E=E, E_und=E_und,
                density=E / (N * (N - 1)), density_und=2 * E_und / (N * (N - 1)),
                total_weight=float(b.w.sum()),
                total_weight_full=tot,
                share_weight_retained=float(b.w.sum()) / tot,
                w_min_retained=float(b.w.min()) if E else np.nan,
                E_full=len(e),
                n_isolated_full=int(N - len(np.union1d(e.i.values, e.j.values))),
                n_isolated_backbone=int(N - len(np.union1d(b.i.values, b.j.values))),
                total_p_all_zones=info[d]["total_p_all"],
                total_p_mainland=info[d]["total_p_mainland"],
            ))
    man = pd.DataFrame(rows)
    man["passes_max_density"] = man.density_und < C.MAX_DENSITY
    man.to_csv(OUT / "manifest.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(man[["date", "groups", "roles", "backbone", "E", "E_und", "density", "density_und",
                   "share_weight_retained", "passes_max_density"]].to_string(index=False))

    fails = man[~man.passes_max_density]
    if len(fails):
        print("NETWORKS FAILING density_und < MAX_DENSITY:")
        print(fails[["date", "backbone", "density_und"]].to_string(index=False))
    # Hard assert on the main backbone of the three main nights; every other failure is reported above.
    main_bb = man[(man.backbone == "w10") & man.date.isin(C.MAIN_NIGHTS.values())]
    bad_main = main_bb[~main_bb.passes_max_density]
    assert bad_main.empty, f"main backbone fails MAX_DENSITY for {bad_main.date.tolist()}"


if __name__ == "__main__":
    main()

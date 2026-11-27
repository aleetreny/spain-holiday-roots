"""Publication audit on each higher-level table's actual matched support."""
import numpy as np
import pandas as pd

from inputs import MONTHS, OUT, n1, od_year


def published_margins():
    definitions = {
        "52997_national": ("turismo_52997_prov_x_prov.csv.gz", 6, [0, 3], []),
        "52997_province_pair": ("turismo_52997_prov_x_prov.csv.gz", 6, [2, 5],
                                [(2, "oprov", "province"), (5, "dprov", "province")]),
        "53000_city_total": ("turismo_53000_munorig_x_provdest.csv.gz", 4, [1],
                             [(0, "orig", "municipality")]),
        "53000_city_to_province": ("turismo_53000_munorig_x_provdest.csv.gz", 4, [3],
                                  [(0, "orig", "municipality"), (3, "dprov", "province")]),
        "53001_village_total": ("turismo_53001_provorig_x_mundest.csv.gz", 4, [0],
                               [(3, "dest", "municipality")]),
        "53001_province_to_village": ("turismo_53001_provorig_x_mundest.csv.gz", 4, [2],
                                     [(2, "oprov", "province"), (3, "dest", "municipality")]),
    }
    national_sublevels = {
        "52997_national": [1, 2, 4, 5],
        "53000_city_total": [2, 3],
        "53001_village_total": [1, 2],
    }
    cache = {}
    for label, (file, concept, select, fields) in definitions.items():
        if file not in cache:
            data = pd.read_csv(n1.INE_TM_RAW / file, sep=";", header=None, skiprows=1,
                               dtype=str, keep_default_na=False)
            cache[file] = data.loc[data[concept].eq("Turistas")]
        data = cache[file]
        if label in national_sublevels:
            keep = (data[select].eq("Total Nacional").all(axis=1)
                    & data[national_sublevels[label]].eq("").all(axis=1))
        else:
            keep = data[select].ne("").all(axis=1)
        data = data.loc[keep].copy()
        keys = []
        for field, name, kind in fields:
            values = data[field].str[:5] if kind == "municipality" else data[field].map(n1.PROV_NAME_TO_CODE)
            if kind == "province":
                unknown = data.loc[values.isna(), field].unique()
                if len(unknown):
                    raise ValueError(f"Unmapped province names in {label}: {unknown}")
            data[name] = values
            keys.append(name)
        data["month"] = data[concept + 1].str.replace("M", "-", regex=False)
        data["margin"] = pd.to_numeric(data[concept + 2].str.replace(".", "", regex=False), errors="coerce")
        yield label, keys, data.dropna(subset=["margin"])[["month", *keys, "margin"]]


def audit(frame, codes, panel, h2):
    od = pd.concat([od_year(year) for year in range(2019, 2027)], ignore_index=True)
    od = od.rename(columns={"mes": "month"})
    national_inventory = od.groupby("month").agg(published_cells=("turistas", "size"),
                                                  minimum_published=("turistas", "min"))
    unmapped = ~od.orig.isin(codes) | ~od.dest.isin(codes)
    national_inventory["published_cells_outside_lau_codes"] = od.loc[unmapped].groupby("month").size()
    national_inventory = national_inventory.fillna(0)
    national_inventory.to_csv(OUT / "od_inventory.csv")
    od["oprov"] = od.orig.str[:2]
    od["dprov"] = od.dest.str[:2]
    rows = []
    month_index = pd.Index(MONTHS, name="month")
    for label, keys, data in published_margins():
        if "orig" in keys:
            data = data.loc[data.orig.isin(frame.orig.unique())]
        if "dest" in keys:
            data = data.loc[data.dest.isin(frame.dest.unique())]
        if "oprov" in keys and "dprov" in keys:
            data = data.loc[data.oprov.ne(data.dprov)]
        join_keys = ["month", *keys]
        aligned = data.merge(od.groupby(join_keys, as_index=False)["turistas"].sum(),
                             how="left", on=join_keys, validate="one_to_one")
        aligned["published_municipal_sum"] = aligned.turistas.fillna(0)
        aligned["gap"] = aligned.margin - aligned.published_municipal_sum
        aligned["overruns"] = aligned.gap.lt(0)
        aligned["overrun_gt_10"] = aligned.gap.lt(-10)
        summary = aligned.groupby("month").agg(
            matched_margin_cells=("margin", "size"), margin_trips=("margin", "sum"),
            published_municipal_trips=("published_municipal_sum", "sum"),
            net_gap=("gap", "sum"), overrun_cells=("overruns", "sum"),
            overrun_gt_10_cells=("overrun_gt_10", "sum"), minimum_gap=("gap", "min"),
        ).reindex(month_index)
        summary["gap_share"] = summary.net_gap / summary.margin_trips
        summary["table"] = label
        rows.append(summary.reset_index())
    pd.concat(rows, ignore_index=True).to_csv(OUT / "monthly_margins.csv", index=False)

    visibility = []
    for month, flow in panel.items():
        reference = month[:4] + "-02"
        for tied in [False, True]:
            mask = (h2 > 0) == tied
            count = int(mask.sum())
            published = flow[mask] >= 0
            baseline = panel[reference][mask] >= 0 if reference in panel else None
            visibility.append(dict(
                month=month, sampled_parent_tie=tied, dyads=count,
                count_published=int(published.sum()),
                share_published=published.mean(),
                count_published_both=np.nan if baseline is None else (published & baseline).sum(),
                share_published_both=np.nan if baseline is None else (published & baseline).mean(),
            ))
    pd.DataFrame(visibility).to_csv(OUT / "visibility.csv", index=False)

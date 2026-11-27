"""Compare the ministry RSS inventory with standard overnight files already local."""
import argparse
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd

from pilot_data import MAIN, OUT, nodes_and_nights as b0

STANDARD = re.compile(r"/por-municipios/pernoctaciones/ficheros-diarios/\d{4}-\d{2}/(\d{8})_Pernoctaciones_municipios\.csv\.gz$")
FILE = re.compile(r"(\d{8})_Pernoctaciones_municipios\.csv\.gz$")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=str(MAIN / "metadata_cache/data_links_v2_2026-11-23.xml"),
                        help="Use - for XML arriving on stdin.")
    args = parser.parse_args()
    source = sys.stdin.buffer if args.manifest == "-" else args.manifest
    tree = ET.parse(source)
    links = {item.text for item in tree.iterfind(".//item/link") if item.text}
    assert links, "No file links found"
    standard = {match[1] for link in links if (match := STANDARD.search(link))}
    bridge = {match[1] for link in links if "cambio_metod" in link and (match := FILE.search(link))}
    raw = {match[1] for file in (MAIN / "_local/external/mitma_pernoct").glob("*.gz")
           if (match := FILE.fullmatch(file.name)) and file.stat().st_size > 0}
    derived = set(b0.available_dates(b0.pernoct_file("2024-01-01"))) | set(
        b0.available_dates(b0.pernoct_file("2022-01-01"))
    )
    dates = pd.date_range("2022-01-01", "2025-12-31")
    encoded = dates.strftime("%Y%m%d")
    days = pd.DataFrame({
        "file_date": dates, "standard_listed": encoded.isin(standard),
        "bridge_listed": encoded.isin(bridge), "raw_local": encoded.isin(raw),
        "derived_local": dates.strftime("%Y-%m-%d").isin(derived),
    })
    days["standard_listed_not_local"] = days["standard_listed"] & ~(days["raw_local"] | days["derived_local"])
    days["standard_listed_local"] = days["standard_listed"] & (days["raw_local"] | days["derived_local"])
    days.to_csv(OUT / "mitma_daily_inventory.csv", index=False)
    days.groupby(days["file_date"].dt.strftime("%Y-%m"))[
        ["standard_listed", "bridge_listed", "raw_local", "derived_local",
         "standard_listed_local", "standard_listed_not_local"]
    ].sum().rename_axis("month").to_csv(OUT / "mitma_monthly_inventory.csv")


if __name__ == "__main__":
    main()

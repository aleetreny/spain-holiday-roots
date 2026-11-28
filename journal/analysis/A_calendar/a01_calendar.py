"""a01_calendar: Official Spanish labour-holiday calendar (2022-2024) from BOE XMLs
and pre-registered night classification rules.

Sources (accessed 27 November 2026):
  - 2022: BOE-A-2021-17113 (Resolución de 14 de octubre de 2021, BOE núm. 252, 21/10/2021)
          https://www.boe.es/diario_boe/txt.php?id=BOE-A-2021-17113
          https://www.boe.es/eli/es/res/2021/10/14/(3)
  - 2023: BOE-A-2022-16755 (Resolución de 7 de octubre de 2022, BOE núm. 247, 14/10/2022)
          https://www.boe.es/diario_boe/txt.php?id=BOE-A-2022-16755
          https://www.boe.es/eli/es/res/2022/10/07/(2)
  - 2024: BOE-A-2023-22014 (Resolución de 23 de octubre de 2023, BOE núm. 257, 27/10/2023)
          https://www.boe.es/diario_boe/txt.php?id=BOE-A-2023-22014
          https://www.boe.es/eli/es/res/2023/10/23/(1)

MITMA incident list (accessed 27 November 2026):
  https://www.transportes.gob.es/ministerio/proyectos-singulares/estudios-de-movilidad-con-big-data/opendata-movilidad
  - Orange outage / unavailable days (2022-2024):
      2023: 26, 27, 30, 31 Oct; 1, 2, 3 Nov
      2024: 4, 18, 19 Apr; 10, 11 Nov
  - DANA Levante anomalies:
      2024: 29, 30, 31 Oct

Run:
  .venv/bin/python -I journal/analysis/A_calendar/a01_calendar.py
"""
from __future__ import annotations

import datetime as dt
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a00_paths import OUT_DIR, ensure_dirs, make_logger

BOE_SOURCES = {
    2022: {
        "boe_id": "BOE-A-2021-17113",
        "title": "Resolución de 14 de octubre de 2021, de la Dirección General de Trabajo, por la que se publica la relación de fiestas laborales para el año 2022",
        "pub_date": "2021-10-21",
        "boe_num": 252,
        "url": "https://www.boe.es/diario_boe/txt.php?id=BOE-A-2021-17113",
        "eli_url": "https://www.boe.es/eli/es/res/2021/10/14/(3)",
        "xml_url": "https://www.boe.es/diario_boe/xml.php?id=BOE-A-2021-17113",
        "accessed": "2026-11-27",
    },
    2023: {
        "boe_id": "BOE-A-2022-16755",
        "title": "Resolución de 7 de octubre de 2022, de la Dirección General de Trabajo, por la que se publica la relación de fiestas laborales para el año 2023",
        "pub_date": "2022-10-14",
        "boe_num": 247,
        "url": "https://www.boe.es/diario_boe/txt.php?id=BOE-A-2022-16755",
        "eli_url": "https://www.boe.es/eli/es/res/2022/10/07/(2)",
        "xml_url": "https://www.boe.es/diario_boe/xml.php?id=BOE-A-2022-16755",
        "accessed": "2026-11-27",
    },
    2024: {
        "boe_id": "BOE-A-2023-22014",
        "title": "Resolución de 23 de octubre de 2023, de la Dirección General de Trabajo, por la que se publica la relación de fiestas laborales para el año 2024",
        "pub_date": "2023-10-27",
        "boe_num": 257,
        "url": "https://www.boe.es/diario_boe/txt.php?id=BOE-A-2023-22014",
        "eli_url": "https://www.boe.es/eli/es/res/2023/10/23/(1)",
        "xml_url": "https://www.boe.es/diario_boe/xml.php?id=BOE-A-2023-22014",
        "accessed": "2026-11-27",
    },
}

MITMA_INCIDENT_SOURCE_URL = (
    "https://www.transportes.gob.es/ministerio/proyectos-singulares/"
    "estudios-de-movilidad-con-big-data/opendata-movilidad"
)

# Official MITMA incident / outage list for 2022-2024 (file dates, accessed 2026-11-27)
MITMA_INCIDENT_DATES: set[str] = {
    # 2023 Orange outage (all 7 also absent from RSS)
    "2023-10-26", "2023-10-27", "2023-10-30", "2023-10-31",
    "2023-11-01", "2023-11-02", "2023-11-03",
    # 2024 Orange incidents (2024-04-18 and 2024-11-10 absent from RSS;
    # 2024-04-04, 2024-04-19, 2024-11-11 present in RSS and excluded by rule a)
    "2024-04-04", "2024-04-18", "2024-04-19",
    "2024-11-10", "2024-11-11",
}

DANA_DATES: set[str] = {"2024-10-29", "2024-10-30", "2024-10-31"}

# Holy Week windows on EVENING date e = d - 1: Friday before Palm Sunday .. Easter Monday (inclusive)
HOLY_WEEK_EVENING_RANGES: dict[int, tuple[dt.date, dt.date]] = {
    2022: (dt.date(2022, 4, 8), dt.date(2022, 4, 18)),   # Easter Sunday 2022-04-17
    2023: (dt.date(2023, 3, 31), dt.date(2023, 4, 10)),  # Easter Sunday 2023-04-09
    2024: (dt.date(2024, 3, 22), dt.date(2024, 4, 1)),   # Easter Sunday 2024-03-31
}

# Pre-registered holiday file dates
EASTER_FILE_DATES: set[str] = {"2022-04-15", "2023-04-07", "2024-03-29"}

MONTH_ES = {
    "ENERO": 1, "FEBRERO": 2, "MARZO": 3, "ABRIL": 4, "MAYO": 5, "JUNIO": 6,
    "JULIO": 7, "AGOSTO": 8, "SEPTIEMBRE": 9, "OCTUBRE": 10, "NOVIEMBRE": 11, "DICIEMBRE": 12,
}

# Target 5 origin CCAA from PREREGISTRATION.md §2 & design.yml
TARGET_CCAA_MAP = {
    "Andalucía": "Andalucia",
    "Cataluña": "Cataluna",
    "Com. Madrid": "Madrid",
    "País Vasco": "Pais_Vasco",
    "Comunitat Valenciana": "Comunidad_Valenciana",
}


def fetch_boe_xml(year: int) -> bytes:
    url = BOE_SOURCES[year]["xml_url"]
    last_err = None
    for attempt in range(1, 5):
        try:
            return subprocess.check_output(
                ["curl", "-4", "-sL", "--max-time", "30", "-A", "Mozilla/5.0", url],
                timeout=35,
            )
        except Exception as e:
            last_err = e
            time.sleep(2 * attempt)
    raise RuntimeError(f"Failed to fetch {url}: {last_err}")


def parse_boe_year(year: int, xml_bytes: bytes) -> pd.DataFrame:
    """Parse official BOE labour-calendar XML table into one row per holiday date."""
    src = BOE_SOURCES[year]
    root = ET.fromstring(xml_bytes)
    tables = root.findall(".//table")
    assert len(tables) >= 1, f"No table found in {src['boe_id']}"
    tbl = tables[0]
    trs = tbl.findall(".//tr")

    # Header row listing the 19 CCAA / autonomous cities
    ccaa_headers = None
    start_row = 0
    for idx, tr in enumerate(trs):
        cells = [" ".join("".join(c.itertext()).split()) for c in tr]
        if any(c.startswith("Andalucía") for c in cells):
            ccaa_headers = cells
            start_row = idx + 1
            break
    assert ccaa_headers is not None and len(ccaa_headers) == 19, f"Unexpected headers: {ccaa_headers}"

    # Match column index for the 5 target CCAA (ignoring footnote tags like '(2)', '(3)')
    col_map: dict[str, int] = {}
    for col_idx, raw_hdr in enumerate(ccaa_headers):
        clean_hdr = re.sub(r"\s*\(\d+\)\s*$", "", raw_hdr).strip()
        if clean_hdr in TARGET_CCAA_MAP:
            col_map[TARGET_CCAA_MAP[clean_hdr]] = col_idx
    assert len(col_map) == 5, f"Missing target CCAA in {year}: {col_map}"

    cur_month = None
    # Aggregate by date because some dates appear on two rows (e.g. 2 May in 2022: Lunes siguiente + Madrid)
    by_date: dict[dt.date, dict] = {}
    for tr in trs[start_row:]:
        cells = [" ".join("".join(c.itertext()).split()) for c in tr]
        if not cells:
            continue
        first = cells[0].strip()
        if not first:
            continue
        upper = first.upper()
        if upper in MONTH_ES:
            cur_month = MONTH_ES[upper]
            continue
        m = re.match(r"^(\d{1,2})\s+(.+)$", first)
        if not m:
            continue
        assert cur_month is not None
        day = int(m.group(1))
        name = m.group(2).rstrip(".")
        d = dt.date(year, cur_month, day)
        marks = cells[1:]
        assert len(marks) == 19, f"Row {first} has {len(marks)} marks != 19"

        # Codes in BOE table: '*' national non-substitutable, '**' national kept, '***' regional
        all_19_marked = all(mk in ("*", "**", "***") for mk in marks)
        non_subst_nat = all(mk == "*" for mk in marks)

        rec = by_date.get(d)
        if rec is None:
            rec = {
                "date": d.isoformat(),
                "year": year,
                "weekday": d.strftime("%a"),
                "holiday_name": name,
                "boe_id": src["boe_id"],
                "boe_url": src["url"],
                "boe_eli_url": src["eli_url"],
                "accessed": src["accessed"],
                "is_national_non_substitutable": non_subst_nat,
                "is_national_all_spain": all_19_marked,
                "Madrid": False,
                "Cataluna": False,
                "Andalucia": False,
                "Comunidad_Valenciana": False,
                "Pais_Vasco": False,
                "raw_marks_target": {},
            }
            by_date[d] = rec
        else:
            if name not in rec["holiday_name"]:
                rec["holiday_name"] += " / " + name
            rec["is_national_non_substitutable"] |= non_subst_nat
            rec["is_national_all_spain"] |= all_19_marked

        for cname, cidx in col_map.items():
            mk = marks[cidx]
            if mk in ("*", "**", "***"):
                rec[cname] = True
                rec["raw_marks_target"][cname] = mk

    rows = []
    for d in sorted(by_date):
        r = by_date[d]
        in_any_5 = any(r[c] for c in ("Madrid", "Cataluna", "Andalucia", "Comunidad_Valenciana", "Pais_Vasco"))
        r["is_regional_target5"] = in_any_5 and not r["is_national_all_spain"]
        r["is_national_or_target5"] = r["is_national_all_spain"] or in_any_5
        # Keep only holidays relevant to National or the 5 target CCAA (plus retain a flag)
        marks_str = "; ".join(f"{k}:{v}" for k, v in sorted(r["raw_marks_target"].items()))
        r["target_marks"] = marks_str
        del r["raw_marks_target"]
        rows.append(r)

    df = pd.DataFrame(rows)
    return df


def build_boe_calendar(save_csv: bool = True) -> pd.DataFrame:
    """Fetch and parse 2022-2024 BOE labour calendars and save to out/boe_holidays_2022_2024.csv."""
    log = make_logger()
    parts = []
    for yr in (2022, 2023, 2024):
        xml_bytes = fetch_boe_xml(yr)
        df_y = parse_boe_year(yr, xml_bytes)
        parts.append(df_y)
        n_nat = int(df_y["is_national_all_spain"].sum())
        n_t5 = int(df_y["is_national_or_target5"].sum())
        log(f"Parsed BOE {yr} ({BOE_SOURCES[yr]['boe_id']}): {len(df_y)} total dates in table, "
            f"{n_nat} all-Spain national, {n_t5} national or target-5 regional")
    cal = pd.concat(parts, ignore_index=True).sort_values("date").reset_index(drop=True)

    # Sanity assertions against manually verified BOE tables
    t5_dates = set(cal.loc[cal["is_national_or_target5"], "date"])
    for must_have in [
        "2022-02-28", "2022-03-19", "2022-04-14", "2022-04-15", "2022-04-18",
        "2022-05-02", "2022-06-06", "2022-06-24", "2022-07-25", "2022-08-15",
        "2022-09-06", "2022-10-12", "2022-11-01", "2022-12-06", "2022-12-08", "2022-12-26",
        "2023-01-02", "2023-01-06", "2023-02-28", "2023-03-20", "2023-04-06", "2023-04-07",
        "2023-04-10", "2023-05-01", "2023-05-02", "2023-06-24", "2023-07-25", "2023-08-15",
        "2023-09-11", "2023-10-09", "2023-10-12", "2023-11-01", "2023-12-06", "2023-12-08",
        "2023-12-25", "2023-12-26",
        "2024-01-01", "2024-01-06", "2024-02-28", "2024-03-19", "2024-03-28", "2024-03-29",
        "2024-04-01", "2024-05-01", "2024-05-02", "2024-06-24", "2024-07-25", "2024-08-15",
        "2024-09-11", "2024-10-09", "2024-10-12", "2024-11-01", "2024-12-06", "2024-12-09",
        "2024-12-25", "2024-12-26",
    ]:
        assert must_have in t5_dates, f"Expected target holiday {must_have} missing!"

    if save_csv:
        ensure_dirs()
        # Save full BOE table (all rows with national/target-5 indicator) AND filtered summary
        out_path = OUT_DIR / "boe_holidays_2022_2024.csv"
        cal.to_csv(out_path, index=False)
        log(f"Saved {out_path} ({len(cal)} rows, {len(t5_dates)} target-5/national dates)")
    return cal


def load_target_holiday_dates(cal_csv: Path | None = None) -> set[dt.date]:
    p = cal_csv or (OUT_DIR / "boe_holidays_2022_2024.csv")
    df = pd.read_csv(p)
    sub = df[df["is_national_or_target5"]]
    return {dt.date.fromisoformat(s) for s in sub["date"]}


def classify_file_dates(file_dates: list[str], holiday_dates: set[dt.date]) -> pd.DataFrame:
    """Classify each MITMA file date d (evening e = d - 1) per PREREGISTRATION.md §2.

    Note: Does NOT apply data-quality exclusions yet; returns calendar eligibility flags
    so the manifest can combine calendar eligibility with quality rules cleanly.
    """
    # Check invariance to Sunday national holidays not listed in BOE table (1 May 2022, 25 Dec 2022, 1 Jan 2023, 8 Dec 2024)
    sunday_nat = {dt.date(2022, 5, 1), dt.date(2022, 12, 25), dt.date(2023, 1, 1), dt.date(2024, 12, 8)}

    rows = []
    for ds in sorted(file_dates):
        d = dt.date.fromisoformat(ds)
        e = d - dt.timedelta(days=1)
        mmdd = ds[5:10]
        e_mmdd = e.isoformat()[5:10]
        e_dow = e.weekday()  # 0=Mon, 1=Tue, 2=Wed, 3=Thu, 4=Fri, 5=Sat, 6=Sun

        # Holiday classes (defined on file date d per design.yml)
        if "08-10" <= mmdd <= "08-18":
            cal_class = "AUG"
        elif mmdd == "12-25":
            cal_class = "XMAS"
        elif ds in EASTER_FILE_DATES:
            cal_class = "EASTER"
        elif mmdd == "11-01":
            cal_class = "ALLSAINTS"
        else:
            cal_class = None

        # Ordinary-night calendar eligibility (defined on evening e = d - 1)
        is_mtw = e_dow in (0, 1, 2)
        is_hol_e = e in holiday_dates
        is_adj_hol_e = ((e - dt.timedelta(days=1)) in holiday_dates) or ((e + dt.timedelta(days=1)) in holiday_dates)
        is_summer_excl = "06-24" <= e_mmdd <= "09-10"
        is_xmas_excl = (e_mmdd >= "12-20") or (e_mmdd <= "01-07")
        hw = HOLY_WEEK_EVENING_RANGES.get(e.year)
        is_holy_week = (hw is not None) and (hw[0] <= e <= hw[1])

        o_eligible = is_mtw and (not is_hol_e) and (not is_adj_hol_e) and (not is_summer_excl) and (not is_xmas_excl) and (not is_holy_week)

        # Verify invariance if Sunday national holidays were also added to holiday_dates
        hol_plus_sun = holiday_dates | sunday_nat
        o_elig_alt = (
            is_mtw
            and (e not in hol_plus_sun)
            and ((e - dt.timedelta(days=1)) not in hol_plus_sun)
            and ((e + dt.timedelta(days=1)) not in hol_plus_sun)
            and (not is_summer_excl) and (not is_xmas_excl) and (not is_holy_week)
        )
        assert o_eligible == o_elig_alt, f"Sunday national holiday changed O eligibility on {ds}"

        if cal_class is None:
            cal_class = "O" if o_eligible else "other"
        else:
            assert not o_eligible, f"Holiday class {cal_class} overlapped with O on {ds}"

        o_excl_reasons = []
        if not is_mtw:
            o_excl_reasons.append("not_MonTueWed_evening")
        if is_hol_e:
            o_excl_reasons.append("holiday_evening")
        if is_adj_hol_e:
            o_excl_reasons.append("adjacent_to_holiday")
        if is_summer_excl:
            o_excl_reasons.append("summer_window_0624_0910")
        if is_xmas_excl:
            o_excl_reasons.append("winter_window_1220_0107")
        if is_holy_week:
            o_excl_reasons.append("holy_week")

        rows.append({
            "file_date": ds,
            "evening_date": e.isoformat(),
            "year": d.year,
            "evening_year": e.year,
            "evening_weekday": e.strftime("%a"),
            "file_weekday": d.strftime("%a"),
            "night_class_calendar": cal_class,
            "is_ordinary_calendar_eligible": o_eligible,
            "evening_is_target_holiday": is_hol_e,
            "evening_adjacent_target_holiday": is_adj_hol_e,
            "in_summer_excl": is_summer_excl,
            "in_winter_excl": is_xmas_excl,
            "in_holy_week": is_holy_week,
            "ordinary_excl_reasons": ";".join(o_excl_reasons),
        })
    return pd.DataFrame(rows)


def main() -> None:
    cal = build_boe_calendar(save_csv=True)
    hols = {dt.date.fromisoformat(s) for s in cal.loc[cal["is_national_or_target5"], "date"]}
    all_dates = [d.strftime("%Y-%m-%d") for d in pd.date_range("2022-01-01", "2024-12-31")]
    cls = classify_file_dates(all_dates, hols)
    print("\nCalendar class counts across all 1,096 days (2022-2024):")
    print(cls.groupby(["year", "night_class_calendar"]).size().unstack(fill_value=0).to_string())
    print("\nTarget-5 + National holiday dates:")
    t5 = cal[cal["is_national_or_target5"]][["date", "weekday", "holiday_name", "is_national_all_spain", "target_marks"]]
    print(t5.to_string(index=False))


if __name__ == "__main__":
    main()

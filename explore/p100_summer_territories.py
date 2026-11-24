"""Build 2024 overnight profiles and matched August sensitivity checks."""
import gzip
import json
from pathlib import Path

import duckdb
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString
import shapely

out = Path('data/derived/summer_territories')
out.mkdir(parents=True, exist_ok=True)
cache = Path('_local/summer_territories_downloads')
zone_dir = Path('raw_data_cache/v2/zonificacion')
rel = pd.read_csv(zone_dir / 'relacion_ine_zonificacionMitma.csv', sep='|', dtype=str)
cross = rel[['municipio_ine','municipio_mitma']].dropna().drop_duplicates()
cross.columns = ['ine_code','zone_id']
pairs = rel[['distrito_mitma','municipio_mitma']].dropna().drop_duplicates()
assert pairs.distrito_mitma.is_unique
zones = gpd.read_file('clean_data/v2/zones/municipios_mitma.gpkg').set_index('id')
zones = zones.loc[zones.index.intersection(cross.zone_id)].sort_index()
names = pd.read_csv(zone_dir / 'zonificacion_municipios/nombres_municipios.csv', sep='|', dtype=str)
zones['name'] = names.set_index('ID')['name'].reindex(zones.index).fillna(zones['name'])
with gzip.open(cache / 'coastline.json.gz', 'rt') as file:
    coast_data = json.load(file)
lines = gpd.GeoSeries([LineString(points) for points in coast_data['lines']], crs=4326)
coast = shapely.STRtree(lines.to_crs(zones.crs).array)
zones['coastal'] = False
zones.iloc[np.unique(coast.query(zones.geometry.make_valid().array, predicate='dwithin',
                                 distance=200)[0]), zones.columns.get_loc('coastal')] = True
zones['municipality_count'] = cross.groupby('zone_id').size()
zones['ine_code'] = cross[cross.zone_id.eq(cross.ine_code)].set_index('zone_id').ine_code
split_codes = cross.loc[cross.duplicated('ine_code', keep=False), 'ine_code']
zones.loc[(zones.municipality_count != 1) | zones.index.str.endswith('_AM') |
          zones.ine_code.isin(split_codes), 'ine_code'] = None
zones[['name','population','coastal','municipality_count','ine_code']].to_csv(out / 'zones.csv')
cross.to_csv(out / 'municipality_crosswalk.csv', index=False)

xy = np.radians(zones.geometry.centroid.to_crs(4326).get_coordinates().to_numpy())
lon, lat = xy.T
km = 6371.0088 * 2 * np.arcsin(np.sqrt(np.clip(
    np.sin((lat[:,None]-lat)/2)**2 + np.cos(lat[:,None])*np.cos(lat)*
    np.sin((lon[:,None]-lon)/2)**2, 0, 1)))
bins = np.digitize(km, np.r_[2:30:2, 30:100:5, 100:300:10, 300:1200:25, 1500, 2000, 2500])
prov = zones.index.str[:2].to_numpy()
same = prov[:,None] == prov
n = len(zones)
db = duckdb.connect(config={'threads': 4})
db.register('district_map', pairs)
base = "'data/derived/pernoct_mun_2024_6m.parquet'"
raw = "'data/derived/pernoct_2024_6m.parquet'"
db.execute(f'CREATE VIEW raw AS SELECT * FROM {raw}')
db.execute(f'CREATE VIEW nights AS SELECT * FROM {base}')
for prefix in ['202407','202208','202308']:
    files = sorted(cache.glob(f'{prefix}*_Pernoctaciones_municipios.csv.gz'))
    if files:
        relation = db.read_csv([str(f) for f in files], delimiter='|', header=True, all_varchar=True)
        relation.create_view('downloads')
        query = """SELECT strptime(fecha,'%Y%m%d')::DATE date, zona_residencia r,
                          zona_pernoctacion s, personas::DOUBLE p FROM downloads"""
        db.execute(f'CREATE TABLE raw_{prefix} AS {query}')
        db.execute(f"""CREATE TABLE nights_{prefix} AS SELECT date, municipio_mitma r, s, sum(p) p
                        FROM raw_{prefix} JOIN district_map ON r=distrito_mitma GROUP BY 1,2,3""")
        db.execute('DROP VIEW downloads')

periods = {
    'august_2024': ('', "month(date)=8"),
    'easter_2024': ('', "date BETWEEN '2024-03-28' AND '2024-03-31'"),
    'all_saints_2024': ('', "date BETWEEN '2024-11-01' AND '2024-11-03'"),
    'ordinary_2024': ('', "month(date)=2 AND dayofweek(date) BETWEEN 2 AND 4"),
    'mid_august_2024': ('', "month(date)=8 AND day(date) BETWEEN 12 AND 18"),
}
for prefix, label in [('202407','july_2024'),('202208','mid_august_2022'),('202308','mid_august_2023')]:
    if db.sql(f"SELECT count(*) FROM information_schema.tables WHERE table_name='nights_{prefix}'").fetchall()[0][0]:
        periods[label] = ('_'+prefix, 'true')
if 'july_2024' in periods:
    july_days = db.sql('SELECT count(DISTINCT date) FROM nights_202407').fetchall()[0][0]
    if july_days == 31:
        db.execute("CREATE VIEW nights_summer AS SELECT * FROM nights_202407 UNION ALL SELECT * FROM nights WHERE month(date)=8")
        db.execute("CREATE VIEW raw_summer AS SELECT * FROM raw_202407 UNION ALL SELECT * FROM raw WHERE month(date)=8")
        periods['summer_2024'] = ('_summer', 'true')
    else:
        del periods['july_2024']

common_days = set(range(12,19))
for label in ('mid_august_2022','mid_august_2023','mid_august_2024'):
    if label in periods:
        suffix, condition = periods[label]
        available = db.sql(f'SELECT DISTINCT day(date) FROM nights{suffix} WHERE {condition}').fetchall()
        common_days &= {d for d, in available}
if common_days:
    day_filter = ','.join(map(str, sorted(common_days)))
    for label, (suffix, condition) in list(periods.items()):
        if label.startswith('mid_'):
            periods[label] = (suffix, f'({condition}) AND day(date) IN ({day_filter})')

origins, destinations, top, colonies, quality, stability = [], [], [], [], [], []

def fit_gravity(weight):
    origin, dest = np.nonzero((weight.sum(1)>0)[:,None] & (weight.sum(0)>0) & ~np.eye(n,dtype=bool))
    keys = (origin, dest, bins[origin,dest])
    observed = weight[origin,dest]
    targets = [np.bincount(k, weights=observed) for k in keys]
    expected = np.ones(len(origin))
    for iteration in range(1, 301):
        for key, target in zip(keys, targets):
            total = np.bincount(key, weights=expected)
            expected *= np.divide(target, total, out=np.zeros_like(target), where=total>0)[key]
        error = max(np.abs(np.bincount(key, weights=expected)-target).sum()/observed.sum()
                    for key,target in zip(keys,targets))
        if error < 1e-6:
            break
    result = np.zeros_like(weight)
    result[origin,dest] = expected
    return result, error, iteration

for period, (suffix, condition) in periods.items():
    table = 'nights'+suffix
    days = db.sql(f'SELECT DISTINCT date FROM {table} WHERE {condition} ORDER BY date').df().date
    count = len(days)
    flows = db.sql(f"""SELECT r,s,sum(p) person_nights, count(*) reported_nights,
                           count(*) FILTER (WHERE p=0) zero_nights FROM {table}
                        WHERE {condition} GROUP BY 1,2""").df()
    oi, di = zones.index.get_indexer(flows.r), zones.index.get_indexer(flows.s)
    valid = (oi>=0)&(di>=0)
    matrix = np.zeros((n,n))
    matrix[oi[valid],di[valid]] = flows.person_nights[valid]/count
    resident = matrix.sum(1)
    present = matrix.sum(0)
    home = matrix.diagonal().copy()
    np.fill_diagonal(matrix, 0)
    origin_count = (resident>0).sum()
    dest_count = (present>0).sum()
    raw_info = db.sql(f"""SELECT count(*) raw_rows, count(*) FILTER (WHERE p=0) raw_zero_rows,
        count(*) FILTER (WHERE p IS NULL) raw_null_rows,
        count(*) FILTER (WHERE municipio_mitma IS NULL) unmapped_raw_rows,
        coalesce(sum(p) FILTER(WHERE municipio_mitma IS NULL),0) unmapped_raw_person_nights
        FROM raw{suffix} LEFT JOIN district_map ON r=distrito_mitma WHERE {condition}""").df().iloc[0].to_dict()
    diag = dict(period=period, days=count, dates=days.dt.strftime('%Y-%m-%d').tolist(),
                origin_zones=int(origin_count), destination_zones=int(dest_count),
                unmapped_person_nights=float(flows.person_nights[~valid].sum()),
                domestic_person_nights=float(present.sum()*count),
                observed_od_days=int(flows.reported_nights[valid].sum()),
                absent_od_days=int(n*n*count-flows.reported_nights[valid].sum()),
                **raw_info)
    if not period.startswith('mid_'):
        expected, error, it = fit_gravity(matrix)
        diag.update(gravity_relative_margin_error=error, gravity_iterations=it)
        assert error < 1e-6, diag
    quality.append(diag)
    if period.startswith('mid_'):
        stability.append(pd.DataFrame(dict(period=period, zone_id=zones.index,
            present=present, multiplier=present/zones.population.where(zones.population>0).to_numpy(),
            away=matrix.sum(1), inbound=matrix.sum(0),
            top_destination=np.where(matrix.sum(1)>0,zones.index[matrix.argmax(1)],None))))
        continue
    daily = db.sql(f'''
        SELECT zone_id, count(DISTINCT date) FILTER (WHERE side='origin') origin_days,
               count(DISTINCT date) FILTER (WHERE side='destination') destination_days
        FROM (
            SELECT date,r zone_id,'origin' side FROM {table} WHERE ({condition}) AND p>0
            UNION ALL
            SELECT date,s zone_id,'destination' side FROM {table} WHERE ({condition}) AND p>0
        ) GROUP BY 1
    ''').df().set_index('zone_id').reindex(zones.index).fillna(0)
    for i,z in enumerate(zones.index):
        row, col = matrix[i], matrix[:,i]
        away, inbound = row.sum(), col.sum()
        base_record = dict(period=period, zone_id=z, days=count)
        shares = row[row>0]/away if away else np.array([])
        incoming = col[col>0]/inbound if inbound else np.array([])
        origins.append(dict(**base_record, reporting_days=int(daily.origin_days.iloc[i]), resident_night_mean=resident[i] if resident[i] else None,
            away_night_mean=away if resident[i] else None, away_share=away/resident[i] if resident[i] else None,
            coastal_share=row[zones.coastal].sum()/away if away else None,
            inland_share=row[~zones.coastal].sum()/away if away else None,
            own_province_share=row[same[i]].sum()/away if away else None,
            other_province_share=row[~same[i]].sum()/away if away else None,
            mean_distance_km=row@km[i]/away if away else None,
            destination_entropy=-shares@np.log(shares) if away else None,
            observed_destination_count=len(shares)))
        population = zones.population.iloc[i]
        destinations.append(dict(**base_record, reporting_days=int(daily.destination_days.iloc[i]), present_night_mean=present[i] if present[i] else None,
            home_night_mean=home[i] if present[i] else None,
            resident_population=population,
            population_multiplier=present[i]/population if present[i] and population>0 else None,
            inbound_night_mean=inbound if present[i] else None,
            dominant_origin_zone_id=zones.index[col.argmax()] if inbound else None,
            dominant_origin_share=col.max()/inbound if inbound else None,
            origin_entropy=-incoming@np.log(incoming) if inbound else None,
            observed_origin_count=len(incoming)))
        for rank,j in enumerate(np.argsort(-row,kind='stable')[:10],1):
            if row[j]>0:
                top.append(dict(**base_record, rank=rank, destination_zone_id=zones.index[j],
                                night_mean=row[j], away_share=row[j]/away, distance_km=km[i,j]))
        e = expected[i]
        ratio = np.divide(row,e,out=np.zeros(n),where=e>0)
        for j in np.flatnonzero((ratio>=2)&(row>=20)&(e>=5)):
            colonies.append(dict(**base_record, destination_zone_id=zones.index[j],
                observed_night_mean=row[j], expected_night_mean=e[j], observed_expected=ratio[j],
                distance_km=km[i,j]))
    print(period, count, flush=True)

for filename, records in [('origins',origins),('destinations',destinations),
                          ('top_destinations',top),('colonies',colonies)]:
    frame = pd.DataFrame(records)
    frame.to_parquet(out / f'{filename}.parquet', index=False)
    assert not frame.duplicated(['period','zone_id'] + (
        ['destination_zone_id'] if filename in ('top_destinations','colonies') else [])).any()
stable = pd.concat(stability)
stable.to_csv(out / 'stability.csv', index=False)
reference = stable[stable.period.eq('mid_august_2024')].set_index('zone_id')
comparisons = []
for period, frame in stable.groupby('period'):
    if period == 'mid_august_2024':
        continue
    joined = frame.set_index('zone_id').merge(reference, left_index=True, right_index=True, suffixes=('', '_2024'))
    joined = joined[(joined.present>0)&(joined.present_2024>0)]
    both = joined.dropna(subset=['top_destination','top_destination_2024'])
    comparisons.append(dict(period=period, matched_zones=len(joined),
        multiplier_spearman=joined.multiplier.corr(joined.multiplier_2024,method='spearman'),
        median_absolute_multiplier_change=(joined.multiplier-joined.multiplier_2024).abs().median(),
        zones_with_two_rankings=len(both),
        top_destination_agreement=(both.top_destination==both.top_destination_2024).mean()))
pd.DataFrame(comparisons).to_csv(out / 'stability_summary.csv', index=False)
(out / 'quality.json').write_text(json.dumps(quality, indent=2, default=lambda v:v.item())+'\n')

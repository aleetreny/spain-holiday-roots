"""Compare August rankings on directly matched municipal pairs, using workbook values only."""
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import duckdb
import pandas as pd

out = Path('data/derived/summer_territories')
file = Path('_local/summer_territories_downloads/ine2024.xlsx')
ns = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
data = []
with ZipFile(file) as archive:
    strings = [''.join(item.itertext()) for item in ET.fromstring(archive.read('xl/sharedStrings.xml'))]
    with archive.open('xl/worksheets/sheet9.xml') as sheet:
        for _, row in ET.iterparse(sheet, events=['end']):
            if row.tag != ns+'row':
                continue
            fields = {}
            for cell in row:
                field = cell.attrib['r'][0]
                value = cell.findtext(ns+'v')
                if field in 'ABDF' and value is not None:
                    fields[field] = strings[int(value)] if cell.get('t') == 's' else value
            if fields.get('A') == '2024-08' and all(k in fields for k in 'BDF'):
                data.append((fields['B'],fields['D'],fields['F']))
            row.clear()
ine = pd.DataFrame(data, columns=['r','s','tourists'])
ine.tourists = pd.to_numeric(ine.tourists, errors='coerce')
zones = pd.read_csv(out/'zones.csv', dtype={'id':str,'ine_code':str})
direct = set(zones.ine_code.dropna())
ine = ine[ine.r.isin(direct)&ine.s.isin(direct)&ine.r.ne(ine.s)&ine.tourists.gt(0)]
assert not ine.duplicated(['r','s']).any()
db = duckdb.connect(config={'threads':4})
db.register('ine',ine)
matched = db.sql("""SELECT r,s,tourists,sum(p) person_nights FROM ine JOIN
    'data/derived/pernoct_mun_2024_6m.parquet' USING (r,s)
    WHERE month(date)=8 GROUP BY 1,2,3 HAVING sum(p)>0""").df()
matched.to_parquet(out/'ine_august_matched.parquet',index=False)
rows = []
for origin, group in matched.groupby('r'):
    if len(group)>=5:
        rows.append(dict(ine_code=origin, matched_pairs=len(group),
            spearman=group.tourists.corr(group.person_nights,method='spearman'),
            top5_overlap=len(set(group.nlargest(5,'tourists').s)&set(group.nlargest(5,'person_nights').s))))
result = pd.DataFrame(rows)
result.to_csv(out/'ine_august_validation.csv',index=False)
print(len(ine),len(matched),len(result),result[['spearman','top5_overlap']].median().to_dict())

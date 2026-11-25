import duckdb, pandas as pd
con=duckdb.connect()
mun="'data/derived/pernoct_mun_2024_6m.parquet'"
d=con.execute(f"select distinct date from {mun} order by 1").df().date
allnov=pd.date_range('2024-11-01','2024-11-30'); print('missing Nov 2024', sorted(set(allnov.date)-set(pd.to_datetime(d).dt.date)))
raw="'data/derived/pernoct_2024_6m.parquet'"
print('min p raw', con.execute(f"select min(p), quantile_cont(p,[0.01,0.05,0.5]) from {raw} where date='2024-08-15'").fetchall())
print('rows p<1 share', con.execute(f"select avg((p<1)::int), avg((p<10)::int) from {raw} where date='2024-08-15'").fetchall())
# pairs presence across Aug 10-18: how many days present
x=con.execute(f"""select r,s,count(distinct date) nd, sum(p) sp from {mun} where date between '2024-08-10' and '2024-08-18' and r<>s group by 1,2""").df()
print(x.nd.value_counts().sort_index())
x['mean_bug']=x.sp/x.nd; x['mean_ok']=x.sp/9
print('pairs w>=10 bug vs ok', (x.mean_bug>=10).sum(), (x.mean_ok>=10).sum())
y=con.execute(f"""select r,s,count(distinct date) nd, sum(p) sp from {mun} where date between '2024-02-05' and '2024-02-29' and dayofweek(date) in (2,3,4) and r<>s group by 1,2""").df()
ndays=con.execute(f"select count(distinct date) from {mun} where date between '2024-02-05' and '2024-02-29' and dayofweek(date) in (2,3,4)").fetchone()[0]
print('feb ndays',ndays, y.nd.value_counts().sort_index().to_dict())
print('pairs w>=10 bug vs ok', (y.sp/y.nd>=10).sum(), (y.sp/ndays>=10).sum(), 'total w bug vs ok', (y.sp/y.nd).sum()/1e6, y.sp.sum()/ndays/1e6)
print('aug total bug vs ok', x.mean_bug.sum()/1e6, x.mean_ok.sum()/1e6)

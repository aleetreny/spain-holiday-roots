# Figure: (a) exodus -> holiday refill gradient (inland zones), (b) Barcelona-area cities: born in Andalusia vs August nights in Andalusia
import pandas as pd, numpy as np, duckdb, sys, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'serif','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
exec(open('explore/p73_exodus_refill.py').read().split("print('zones'")[0])
inl = D[D.coastal==0].copy(); inl['q'] = pd.qcut(inl.loss_peak_2021, 5, labels=False)
lab = [f"{int(round(100*v))}%" for v in inl.groupby('q').loss_peak_2021.median()]
lab = [('+'+l[1:] if l.startswith('-') else '−'+l) for l in lab]  # loss>0 shown as negative change
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
x = np.arange(5); w = 0.26
for i,(col,c,l) in enumerate([('mult_feb','#bbbbbb','Ordinary February weekday'),('mult_nov1','#7fb3d5','All Saints (1 Nov)'),('mult_aug','#1b4f72','Mid-August')]):
    med = inl.groupby('q')[col].median()
    ax[0].bar(x+(i-1)*w, med.values-1, w, bottom=1, color=c, label=l)
ax[0].axhline(1, color='k', lw=0.8)
ax[0].set_xticks(x); ax[0].set_xticklabels([f"Q{i+1}\n{lab[i]}" for i in range(5)])
ax[0].set_xlabel('Inland zones by population change since their historical peak (quintiles)')
ax[0].set_ylabel('People sleeping there / residents')
ax[0].set_title('(a) The emptier a village became, the more it refills', loc='left', fontsize=11)
ax[0].legend(frameon=False, fontsize=9, loc='upper left')
# (b)
B = pd.read_csv('_local/external/ine_birth_residence/mun_residence_x_prov_birth_censo_2021_2025.csv', dtype={'ine_code':str,'prov_birth':str})
B = B[(B.year==2024)&B.prov_birth.str.match(r'^\d\d$')]
AND = ['04','11','14','18','21','23','29','41']
cities = sorted(B[B.ine_code.str[:2]=='08'].ine_code.unique())
con = duckdb.connect(); cl = ",".join(f"'{c}'" for c in cities)
f = con.execute(f"SELECT r c, substr(s,1,2) p, sum(p)/count(DISTINCT date) w FROM 'data/derived/pernoct_mun_2024_6m.parquet' WHERE date BETWEEN '2024-08-10' AND '2024-08-18' AND r IN ({cl}) AND r<>s GROUP BY 1,2").df()
f = f[f.p!='08']
y = 100*f[f.p.isin(AND)].groupby('c').w.sum()/f.groupby('c').w.sum()
xb = 100*B[B.prov_birth.isin(AND)].groupby('ine_code').persons.sum()/B.groupby('ine_code').persons.sum()
nm = B.drop_duplicates('ine_code').set_index('ine_code').municipality_name
d = pd.DataFrame({'x':xb,'y':y}).dropna()
ax[1].scatter(d.x, d.y, s=40, color='#1b4f72')
for c in ['08019','08205','08245','08101','08073','08301','08187','08279']:
    if c in d.index: ax[1].annotate(nm[c].replace(", L'","").replace('Hospitalet de Llobregat',"L'Hospitalet"), (d.x[c], d.y[c]), xytext=(4,3), textcoords='offset points', fontsize=8.5)
b = np.polyfit(d.x, d.y, 1); xs = np.linspace(d.x.min(), d.x.max(), 10); ax[1].plot(xs, b[0]*xs+b[1], color='#7fb3d5', lw=1.5)
ax[1].set_xlabel('Residents born in Andalusia (%)'); ax[1].set_ylabel('August nights spent in Andalusia (% of nights\noutside Barcelona province)')
ax[1].set_title(f'(b) Barcelona-area cities: holidays follow birthplace (r = {d.x.corr(d.y):.2f})', loc='left', fontsize=11)
fig.tight_layout(); fig.savefig('reports/figs/fig9_exodus_holidays.png', dpi=200)
print('saved', d.x.corr(d.y))

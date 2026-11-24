# Check 5: Echo test (old 1988-1995 vs recent 2015-2024 interprovincial migration)
#          + Birthplace ties vs Meta SCI friendship ties predicting 2024 holiday overnight flows.
import sys
sys.path.insert(0, 'explore')
import duckdb, numpy as np, pandas as pd, pyfixest as pf
from common import nodes

N = nodes().set_index('id')
M = N[(~N.foreign) & (~N.island) & (~N.prov.isin(['51', '52']))].copy()
M['w'] = M.population.fillna(0) + 1
pc = M.groupby('prov').apply(
    lambda g: pd.Series({'x': np.average(g.x, weights=g.w), 'y': np.average(g.y, weights=g.w), 'pop': g.population.fillna(0).sum()})
)
provs = sorted(pc.index.tolist())  # 47 mainland provinces

# 1. Load 2024 provincial overnight flows
con = duckdb.connect()
q = "'data/derived/pernoct_mun_2024_6m.parquet'"
slices = {
    'Feb weekday': "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",
    'ordinary Sat': "date IN ('2024-02-10','2024-02-17','2024-02-24','2024-10-19','2024-11-16','2024-11-23')",
    'Good Friday': "date='2024-03-29'",
    'mid-Aug': "date BETWEEN '2024-08-10' AND '2024-08-18'",
    'All Saints': "date='2024-11-01'",
    'Christmas Eve': "date='2024-12-24'",
}

def prov_flows(cond):
    df = con.execute(
        f"SELECT substr(r,1,2) rp, substr(s,1,2) sp, sum(p)/count(DISTINCT date) w "
        f"FROM {q} WHERE {cond} GROUP BY 1,2"
    ).df()
    return df.set_index(['rp', 'sp']).w

# 2. Load interprovincial migration 1988-2024
mig = pd.read_csv('_local/external/ine_interprovincial_migration/interprovincial_migration_1988_2024.csv',
                  dtype={'prov_origin': str, 'prov_dest': str})
print("Migration sources & years:")
print(mig.groupby('source')['year'].agg(['min', 'max', 'count']))

# Note: 2021 has BOTH EM_24379 and EMCR_69742.
# For clean recent periods without double-counting 2021:
#   old_88_95: EVR_microdata 1988-1995 (8 years)
#   mid_96_07: EVR 1996-2007 (12 years)
#   rec_15_24_spliced: EM 2015-2021 rescaled to EMCR level using 2021 national ratio + EMCR 2022-2024
#   rec_15_21_em: EM_24379 2015-2021 (pure EM series, no break)
#   rec_21_24_emcr: EMCR_69742 2021-2024 (pure EMCR series, no break)
m_old = mig[(mig.year >= 1988) & (mig.year <= 1995)].groupby(['prov_origin', 'prov_dest']).migrations.sum()
m_mid = mig[(mig.year >= 1996) & (mig.year <= 2007)].groupby(['prov_origin', 'prov_dest']).migrations.sum()
m_em = mig[(mig.source == 'EM_24379') & (mig.year >= 2015) & (mig.year <= 2021)]
m_emcr = mig[(mig.source == 'EMCR_69742') & (mig.year >= 2021) & (mig.year <= 2024)]
scale_21 = m_emcr[m_emcr.year == 2021].migrations.sum() / m_em[m_em.year == 2021].migrations.sum()
print(f"2021 EMCR/EM scale factor: {scale_21:.4f}")
m_rec_em = m_em.groupby(['prov_origin', 'prov_dest']).migrations.sum()
m_rec_emcr = m_emcr.groupby(['prov_origin', 'prov_dest']).migrations.sum()
m_rec_spliced = (
    m_em[m_em.year < 2021].groupby(['prov_origin', 'prov_dest']).migrations.sum() * scale_21
).add(m_rec_emcr, fill_value=0)

# 3. Load birthplace stock (1998 and 2024)
br = pd.read_csv('_local/external/ine_birth_residence/prov_birth_x_prov_residence.csv',
                 dtype={'prov_birth': str, 'prov_residence': str})
b_2024 = br[(br.year == 2024) & (br.source == 'censo_anual_68523')].set_index(['prov_birth', 'prov_residence']).persons
b_1998 = br[(br.year == 1998) & (br.source == 'padron_continuo')].set_index(['prov_birth', 'prov_residence']).persons

# 4. Load Meta SCI (47 mainland provinces map 1-to-1 with NUTS3!)
sci = pd.read_csv('_local/external/meta_sci_spain/sci_spain_nuts3.csv',
                  dtype={'user_prov_ine': str, 'friend_prov_ine': str})
sci_m = sci[sci.user_prov_ine.isin(provs) & sci.friend_prov_ine.isin(provs)]
assert len(sci_m) == 47 * 47, f"Expected 47x47 mainland SCI, got {len(sci_m)}"
sci_s = sci_m.set_index(['user_prov_ine', 'friend_prov_ine']).scaled_sci

# Build directed mainland province pairs (rp = residence/holiday origin, sp = overnight/holiday destination)
pairs = pd.DataFrame([(a, b) for a in provs for b in provs if a != b], columns=['rp', 'sp'])
pairs['km'] = np.hypot(pc.loc[pairs.rp, 'x'].values - pc.loc[pairs.sp, 'x'].values,
                       pc.loc[pairs.rp, 'y'].values - pc.loc[pairs.sp, 'y'].values) / 1000.0
pairs['kmbin'] = pd.qcut(pairs.km, 20, labels=False)

# Directionality:
#   "in" (sp -> rp): migrants from holiday-dest sp to holiday-origin rp ("roots" direction!)
#   "out" (rp -> sp): migrants from holiday-origin rp to holiday-dest sp ("reverse" direction)
def add_dir(df, ser, name):
    df[f'l_{name}_in'] = np.log1p(ser.reindex(list(zip(df.sp, df.rp))).fillna(0).values)
    df[f'l_{name}_out'] = np.log1p(ser.reindex(list(zip(df.rp, df.sp))).fillna(0).values)

add_dir(pairs, m_old, 'mig88_95')
add_dir(pairs, m_mid, 'mig96_07')
add_dir(pairs, m_rec_spliced, 'mig15_24')
add_dir(pairs, m_rec_em, 'mig15_21')
add_dir(pairs, m_rec_emcr, 'mig21_24')
add_dir(pairs, b_2024, 'birth24')
add_dir(pairs, b_1998, 'birth98')
pairs['l_sci'] = np.log(sci_s.reindex(list(zip(pairs.rp, pairs.sp))).values)

# Check collinearity (raw and residualized on rp + sp + kmbin FE)
print("\n=== Correlation matrix of log predictors (raw) ===")
cols_c = ['l_mig88_95_in', 'l_mig15_24_in', 'l_birth98_in', 'l_birth24_in', 'l_sci']
print(pairs[cols_c].corr().round(3).to_string())

resid = pd.DataFrame(index=pairs.index)
for c in cols_c + ['l_mig88_95_out', 'l_mig15_24_out', 'l_birth24_out']:
    resid[c] = pf.feols(f"{c} ~ 1 | rp + sp + kmbin", data=pairs).resid()
print("\n=== Correlation matrix AFTER residualizing on rp + sp + kmbin FE ===")
print(resid[cols_c].corr().round(3).to_string())

# Check what explains Meta SCI itself (birthplace vs old/recent migration)
m_sci = pf.feols("l_sci ~ l_birth24_in + l_birth24_out + l_mig88_95_in + l_mig15_24_in | rp + sp + kmbin", data=pairs, vcov='hetero')
m_sci_b = pf.feols("l_sci ~ l_birth24_in + l_birth24_out | rp + sp + kmbin", data=pairs, vcov='hetero')
m_sci_0 = pf.feols("l_sci ~ 1 | rp + sp + kmbin", data=pairs, vcov='hetero')
print(f"\n=== Explaining Meta SCI (OLS with rp+sp+kmbin FE) ===")
print(f"FE-only R2={m_sci_0._r2:.3f}; +birth24 R2={m_sci_b._r2:.3f} (within R2={m_sci_b._r2_within:.3f}); +mig R2={m_sci._r2:.3f} (within R2={m_sci._r2_within:.3f})")
print(m_sci.tidy()[['Estimate', 'Std. Error', 't value', 'Pr(>|t|)']].round(4).to_string())

# Helper to run PPML and extract stats
def run_ppml(df, dep, rhs):
    m0 = pf.fepois(f"{dep} ~ 1 | rp + sp + kmbin", data=df, vcov='hetero')
    m = pf.fepois(f"{dep} ~ {rhs} | rp + sp + kmbin", data=df, vcov='hetero')
    dev0, dev1 = float(np.sum(m0.deviance)), float(np.sum(m.deviance))
    t = m.tidy()
    return t, 1.0 - dev1 / dev0

ref = prov_flows(slices['Feb weekday'])

print("\n" + "="*100)
print("PART A: ECHO TEST — OLD (1988-95) vs RECENT (2015-24) INTERPROVINCIAL MIGRATION")
print("="*100)
rows_a = []
for lab, cond in slices.items():
    f = prov_flows(cond)
    p = pairs.copy()
    p['w'] = f.reindex(list(zip(p.rp, p.sp))).fillna(0).values
    p['w_ref'] = ref.reindex(list(zip(p.rp, p.sp))).fillna(0).values
    p['excess'] = (p.w - p.w_ref).clip(lower=0)
    deps = ['w'] if lab == 'Feb weekday' else ['w', 'excess']
    for dep in deps:
        # 1. Alone: old_in only vs rec_in only
        t_o, d_o = run_ppml(p, dep, "l_mig88_95_in")
        t_r, d_r = run_ppml(p, dep, "l_mig15_24_in")
        # 2. Horse race (incoming = roots direction only)
        t_hr, d_hr = run_ppml(p, dep, "l_mig88_95_in + l_mig15_24_in")
        # 3. Full bidirectional horse race (both in and out, old and recent)
        t_bi, d_bi = run_ppml(p, dep, "l_mig88_95_in + l_mig15_24_in + l_mig88_95_out + l_mig15_24_out")
        # 4. Pure EMCR 2021-24 vs 1988-95
        t_emcr, _ = run_ppml(p, dep, "l_mig88_95_in + l_mig21_24_in + l_mig88_95_out + l_mig21_24_out")
        rows_a.append(dict(
            slice=lab, dep=dep,
            dev_old_only=d_o, dev_rec_only=d_r,
            hr_old_in=t_hr.loc['l_mig88_95_in', 'Estimate'], hr_old_se=t_hr.loc['l_mig88_95_in', 'Std. Error'],
            hr_rec_in=t_hr.loc['l_mig15_24_in', 'Estimate'], hr_rec_se=t_hr.loc['l_mig15_24_in', 'Std. Error'],
            dev_hr=d_hr,
            bi_old_in=t_bi.loc['l_mig88_95_in', 'Estimate'], bi_old_in_se=t_bi.loc['l_mig88_95_in', 'Std. Error'],
            bi_rec_in=t_bi.loc['l_mig15_24_in', 'Estimate'], bi_rec_in_se=t_bi.loc['l_mig15_24_in', 'Std. Error'],
            bi_old_out=t_bi.loc['l_mig88_95_out', 'Estimate'], bi_rec_out=t_bi.loc['l_mig15_24_out', 'Estimate'],
            dev_bi=d_bi,
            emcr_old_in=t_emcr.loc['l_mig88_95_in', 'Estimate'], emcr_rec_in=t_emcr.loc['l_mig21_24_in', 'Estimate'],
        ))
df_a = pd.DataFrame(rows_a)
df_a.to_csv('explore/out_p85_echo.csv', index=False)
pd.set_option('display.width', 240)
pd.set_option('display.max_columns', 30)
print(df_a.round(3).to_string(index=False))

print("\n" + "="*100)
print("PART A2: BONUS — 1998 BIRTHPLACE STOCK vs 2024 BIRTHPLACE STOCK (and + RECENT MIGRATION)")
print("="*100)
rows_b98 = []
for lab, cond in slices.items():
    f = prov_flows(cond)
    p = pairs.copy()
    p['w'] = f.reindex(list(zip(p.rp, p.sp))).fillna(0).values
    p['w_ref'] = ref.reindex(list(zip(p.rp, p.sp))).fillna(0).values
    p['excess'] = (p.w - p.w_ref).clip(lower=0)
    for dep in (['w'] if lab == 'Feb weekday' else ['w', 'excess']):
        t1, d1 = run_ppml(p, dep, "l_birth98_in + l_birth24_in")
        t2, d2 = run_ppml(p, dep, "l_birth98_in + l_mig15_24_in + l_mig15_24_out")
        rows_b98.append(dict(
            slice=lab, dep=dep,
            b98_vs_b24_b98=t1.loc['l_birth98_in', 'Estimate'], b98_se=t1.loc['l_birth98_in', 'Std. Error'],
            b98_vs_b24_b24=t1.loc['l_birth24_in', 'Estimate'], b24_se=t1.loc['l_birth24_in', 'Std. Error'],
            b98_plus_rec_b98=t2.loc['l_birth98_in', 'Estimate'], b98_r_se=t2.loc['l_birth98_in', 'Std. Error'],
            b98_plus_rec_in=t2.loc['l_mig15_24_in', 'Estimate'], b98_plus_rec_out=t2.loc['l_mig15_24_out', 'Estimate'],
            dev_b98_rec=d2,
        ))
print(pd.DataFrame(rows_b98).round(3).to_string(index=False))

print("\n" + "="*100)
print("PART B: BIRTHPLACE TIES vs META SCI FRIENDSHIP TIES")
print("="*100)
rows_sci = []
for lab, cond in slices.items():
    f = prov_flows(cond)
    p = pairs.copy()
    p['w'] = f.reindex(list(zip(p.rp, p.sp))).fillna(0).values
    p['w_ref'] = ref.reindex(list(zip(p.rp, p.sp))).fillna(0).values
    p['excess'] = (p.w - p.w_ref).clip(lower=0)
    for dep in (['w'] if lab == 'Feb weekday' else ['w', 'excess']):
        t_b, d_b = run_ppml(p, dep, "l_birth24_in + l_birth24_out")
        t_s, d_s = run_ppml(p, dep, "l_sci")
        t_bs, d_bs = run_ppml(p, dep, "l_birth24_in + l_birth24_out + l_sci")
        t_b1s, d_b1s = run_ppml(p, dep, "l_birth24_in + l_sci")
        rows_sci.append(dict(
            slice=lab, dep=dep,
            dev_birth_only=d_b, dev_sci_only=d_s, dev_both=d_bs,
            # birth alone
            b_in_alone=t_b.loc['l_birth24_in', 'Estimate'], b_out_alone=t_b.loc['l_birth24_out', 'Estimate'],
            # sci alone
            sci_alone=t_s.loc['l_sci', 'Estimate'], sci_alone_se=t_s.loc['l_sci', 'Std. Error'],
            # horse race: birth_in + birth_out + sci
            b_in_hr=t_bs.loc['l_birth24_in', 'Estimate'], b_in_hr_se=t_bs.loc['l_birth24_in', 'Std. Error'],
            b_out_hr=t_bs.loc['l_birth24_out', 'Estimate'], b_out_hr_se=t_bs.loc['l_birth24_out', 'Std. Error'],
            sci_hr=t_bs.loc['l_sci', 'Estimate'], sci_hr_se=t_bs.loc['l_sci', 'Std. Error'],
            # 1-way birth_in + sci
            b1_in=t_b1s.loc['l_birth24_in', 'Estimate'], sci_1=t_b1s.loc['l_sci', 'Estimate'],
        ))
df_sci = pd.DataFrame(rows_sci)
df_sci.to_csv('explore/out_p85_sci.csv', index=False)
print(df_sci.round(3).to_string(index=False))

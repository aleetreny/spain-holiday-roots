# Check 5b: Echo test restricted to INLAND destination provinces (eliminating recent beach/lifestyle retirement migration confound!)
# + also test 1960-1981 historical depopulation & birthplace stock vs 1988-95 vs 2015-24 migration
import sys
sys.path.insert(0, 'explore')
import duckdb, numpy as np, pandas as pd, pyfixest as pf
from common import nodes

COASTAL_PROVS = {'03','04','08','11','12','15','17','18','20','21','27','29','30','33','36','39','43','46','48'}

N = nodes().set_index('id')
M = N[(~N.foreign) & (~N.island) & (~N.prov.isin(['51', '52']))].copy()
M['w'] = M.population.fillna(0) + 1
pc = M.groupby('prov').apply(
    lambda g: pd.Series({'x': np.average(g.x, weights=g.w), 'y': np.average(g.y, weights=g.w)})
)
provs = sorted(pc.index.tolist())

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
    return con.execute(
        f"SELECT substr(r,1,2) rp, substr(s,1,2) sp, sum(p)/count(DISTINCT date) w "
        f"FROM {q} WHERE {cond} GROUP BY 1,2"
    ).df().set_index(['rp', 'sp']).w

mig = pd.read_csv('_local/external/ine_interprovincial_migration/interprovincial_migration_1988_2024.csv',
                  dtype={'prov_origin': str, 'prov_dest': str})
m_old = mig[(mig.year >= 1988) & (mig.year <= 1995)].groupby(['prov_origin', 'prov_dest']).migrations.sum()
m_em = mig[(mig.source == 'EM_24379') & (mig.year >= 2015) & (mig.year <= 2021)]
m_emcr = mig[(mig.source == 'EMCR_69742') & (mig.year >= 2021) & (mig.year <= 2024)]
scale_21 = m_emcr[m_emcr.year == 2021].migrations.sum() / m_em[m_em.year == 2021].migrations.sum()
m_rec = (m_em[m_em.year < 2021].groupby(['prov_origin', 'prov_dest']).migrations.sum() * scale_21).add(
    m_emcr.groupby(['prov_origin', 'prov_dest']).migrations.sum(), fill_value=0
)

br = pd.read_csv('_local/external/ine_birth_residence/prov_birth_x_prov_residence.csv',
                 dtype={'prov_birth': str, 'prov_residence': str})
b_2024 = br[(br.year == 2024) & (br.source == 'censo_anual_68523')].set_index(['prov_birth', 'prov_residence']).persons
b_1998 = br[(br.year == 1998) & (br.source == 'padron_continuo')].set_index(['prov_birth', 'prov_residence']).persons

sci = pd.read_csv('_local/external/meta_sci_spain/sci_spain_nuts3.csv',
                  dtype={'user_prov_ine': str, 'friend_prov_ine': str})
sci_s = sci[sci.user_prov_ine.isin(provs) & sci.friend_prov_ine.isin(provs)].set_index(['user_prov_ine', 'friend_prov_ine']).scaled_sci

pairs = pd.DataFrame([(a, b) for a in provs for b in provs if a != b], columns=['rp', 'sp'])
pairs['km'] = np.hypot(pc.loc[pairs.rp, 'x'].values - pc.loc[pairs.sp, 'x'].values,
                       pc.loc[pairs.rp, 'y'].values - pc.loc[pairs.sp, 'y'].values) / 1000.0
pairs['kmbin'] = pd.qcut(pairs.km, 20, labels=False)
pairs['dest_inland'] = (~pairs.sp.isin(COASTAL_PROVS)).astype(int)

for ser, name in [(m_old, 'old'), (m_rec, 'rec'), (b_2024, 'b24'), (b_1998, 'b98')]:
    pairs[f'l_{name}_in'] = np.log1p(ser.reindex(list(zip(pairs.sp, pairs.rp))).fillna(0).values)
    pairs[f'l_{name}_out'] = np.log1p(ser.reindex(list(zip(pairs.rp, pairs.sp))).fillna(0).values)
pairs['l_sci'] = np.log(sci_s.reindex(list(zip(pairs.rp, pairs.sp))).values)
# Symmetry decomposition: roots = sp->rp (born/migrated FROM holiday dest TO residence) vs reverse = rp->sp
# Also net historic exodus = log((1+b98_in)/(1+b98_out)) and net old mig vs net recent mig!
pairs['l_old_sym'] = 0.5 * (pairs.l_old_in + pairs.l_old_out)
pairs['l_rec_sym'] = 0.5 * (pairs.l_rec_in + pairs.l_rec_out)
pairs['l_b24_sym'] = 0.5 * (pairs.l_b24_in + pairs.l_b24_out)

def run_ppml(df, dep, rhs):
    m0 = pf.fepois(f"{dep} ~ 1 | rp + sp + kmbin", data=df, vcov='hetero')
    m = pf.fepois(f"{dep} ~ {rhs} | rp + sp + kmbin", data=df, vcov='hetero')
    return m.tidy(), 1.0 - float(np.sum(m.deviance)) / float(np.sum(m0.deviance))

ref = prov_flows(slices['Feb weekday'])

print("=== INLAND DESTINATION PROVINCES ONLY (28 inland provs): Echo test (1988-95 vs 2015-24) & Birthplace vs Meta SCI ===")
rows = []
for lab, cond in slices.items():
    f = prov_flows(cond)
    p = pairs[pairs.dest_inland == 1].copy()
    p['w'] = f.reindex(list(zip(p.rp, p.sp))).fillna(0).values
    p['w_ref'] = ref.reindex(list(zip(p.rp, p.sp))).fillna(0).values
    p['excess'] = (p.w - p.w_ref).clip(lower=0)
    for dep in (['w'] if lab == 'Feb weekday' else ['w', 'excess']):
        # 1. Old vs recent incoming (roots direction)
        t_or, d_or = run_ppml(p, dep, "l_old_in + l_rec_in")
        _, d_o = run_ppml(p, dep, "l_old_in")
        _, d_r = run_ppml(p, dep, "l_rec_in")
        # 2. Old vs recent symmetric
        t_sym, _ = run_ppml(p, dep, "l_old_sym + l_rec_sym")
        # 3. Birthplace (in + out) vs Meta SCI
        t_bs, d_bs = run_ppml(p, dep, "l_b24_in + l_b24_out + l_sci")
        _, d_b = run_ppml(p, dep, "l_b24_in + l_b24_out")
        _, d_s = run_ppml(p, dep, "l_sci")
        # 4. Birthplace 1998 (exodus stock) + recent migration (in + out)
        t_98r, _ = run_ppml(p, dep, "l_b98_in + l_rec_in + l_rec_out")
        rows.append(dict(
            slice=lab, dep=dep,
            dev_old=d_o, dev_rec=d_r,
            old_in=t_or.loc['l_old_in','Estimate'], old_se=t_or.loc['l_old_in','Std. Error'],
            rec_in=t_or.loc['l_rec_in','Estimate'], rec_se=t_or.loc['l_rec_in','Std. Error'],
            old_sym=t_sym.loc['l_old_sym','Estimate'], rec_sym=t_sym.loc['l_rec_sym','Estimate'],
            dev_birth=d_b, dev_sci=d_s, dev_both=d_bs,
            b24_in=t_bs.loc['l_b24_in','Estimate'], b24_out=t_bs.loc['l_b24_out','Estimate'],
            sci=t_bs.loc['l_sci','Estimate'], sci_se=t_bs.loc['l_sci','Std. Error'],
            b98_ctrl_rec=t_98r.loc['l_b98_in','Estimate'], b98_se=t_98r.loc['l_b98_in','Std. Error'],
        ))
df = pd.DataFrame(rows)
df.to_csv('explore/out_p85b_echo_inland.csv', index=False)
pd.set_option('display.width', 240)
print(df.round(3).to_string(index=False))

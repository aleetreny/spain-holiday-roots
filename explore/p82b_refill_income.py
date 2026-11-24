# Check 1D: Village-refill gradient (p73) controlling for municipal ADRH income (2024)
# + Matched-income comparisons within Barcelona and Madrid metros!
import sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, 'explore')
import numpy as np, pandas as pd, pyfixest as pf

H = pd.read_csv('_local/external/ine_historical_pop/municipal_population_1900_2025.csv',
                dtype={'ine_code': str, 'province_code': str})
rel = pd.read_csv('raw_data_cache/v2/zonificacion/relacion_ine_zonificacionMitma.csv', sep='|', dtype=str)[
    ['municipio_ine', 'municipio_mitma']
].drop_duplicates()
adrh_m = pd.read_csv('_local/external/ine_adrh/adrh_municipios.csv', dtype={'ine_code': str})
inc24 = adrh_m[(adrh_m.year == 2024) & (adrh_m.indicator == 'Renta neta media por persona')].set_index('ine_code').value

H = H.merge(rel, left_on='ine_code', right_on='municipio_ine', how='inner')
H['inc24'] = H.ine_code.map(inc24)
H['pop_w'] = H.padron_2024.fillna(H.pop_2021).fillna(1.0)
H['inc_x_pop'] = H['inc24'] * H['pop_w']

cols = ['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960','pop_1970','pop_1981','pop_2001','pop_2021']
Z = H.groupby('municipio_mitma')[cols].sum(min_count=1)
# Population-weighted average ADRH income per MITMA zone
Z_inc = H.dropna(subset=['inc24']).groupby('municipio_mitma').apply(
    lambda g: np.average(g.inc24, weights=g.pop_w), include_groups=False
).rename('inc24')
Z = Z.join(Z_inc)
Z['peak'] = Z[['pop_1900','pop_1910','pop_1920','pop_1930','pop_1940','pop_1950','pop_1960']].max(axis=1)
Z['loss_peak_2021'] = 1 - Z.pop_2021 / Z.peak

R = pd.read_csv('explore/out_refill.csv', dtype={'z': str, 'prov': str}).set_index('z')
cf = pd.read_csv('explore/out_p86_coastal_flags.csv', dtype={'id': str}).set_index('id')
D = R.join(Z, how='inner').join(cf, how='inner')
D['prov'] = D.index.str[:2]
D = D[(D.resid >= 200) & ~D.prov.isin(['07','35','38','51','52']) & D.inc24.notna()].copy()
D['lmult_aug'] = np.log(D.mult_aug)
D['lmult_nov1'] = np.log(D.mult_nov1)
D['lres'] = np.log(D.resid)
D['linc'] = np.log(D.inc24)
D['loss_rank'] = D.loss_peak_2021.rank(pct=True)

print("=== VILLAGE-REFILL GRADIENT CONTROLLING FOR DESTINATION MUNICIPAL INCOME (ADRH 2024) ===")
for cflag, clab in [('coastal_p73', 'p73 inland'), ('coastal_true', 'true inland')]:
    inl = D[D[cflag] == 0].copy()
    print(f"\n[{clab}, n={len(inl)}, corr(loss_rank, linc)={inl.loss_rank.corr(inl.linc):+.3f}]")
    for y in ['lmult_aug', 'lmult_nov1']:
        m0 = pf.feols(f"{y} ~ loss_rank + lres | prov", data=inl, vcov={'CRV1':'prov'}).tidy()
        m1 = pf.feols(f"{y} ~ loss_rank + lres + linc | prov", data=inl, vcov={'CRV1':'prov'}).tidy()
        print(f"  {y:10s}: base={m0.loc['loss_rank','Estimate']:+.3f} ({m0.loc['loss_rank','Std. Error']:.3f})  ->  "
              f"+linc: loss_rank={m1.loc['loss_rank','Estimate']:+.3f} ({m1.loc['loss_rank','Std. Error']:.3f}), "
              f"linc={m1.loc['linc','Estimate']:+.3f} ({m1.loc['linc','Std. Error']:.3f})")

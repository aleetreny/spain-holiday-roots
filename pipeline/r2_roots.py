"""r2_roots: R2, edge level. Does the holiday overnight network align with the roots
(birthplace) network more than the ordinary network does, beyond distance, size,
income and beaches?

Run:  .venv/bin/python -I pipeline/r2_roots.py   (writes results/r2_*.csv, results/r2_run.log)

Design (decisions in .planning/2026-11-24-consolidation/decisions.md, tag [R2]):
- Roots: INE Censo anual 2024 (1 Jan), residents of city c born in province p
  (table 68545, 151 large municipalities) and the 47x47 province matrix (table 68523).
- City set: the 137 mainland census cities (all are single-municipality MITMA zones).
  Under origin-province x destination-province fixed effects only cities that share
  their province with another census city identify the roots coefficient: 112 cities
  in 22 provinces. Single-city provinces are perfectly fitted by their pair effects,
  so dropping them changes nothing; they are dropped from the estimation sample.
- Outcome: one night at a time (never a window mean), persons resident in city c
  sleeping in any mainland zone of province p != prov(c), from
  results/b01_networks/edges_<date>__full.parquet.
- PPML: y_cp ~ beta*log(share born in p) + log km + (log km)^2 + income x dest-prov slopes
  | city FE + (origin prov x dest prov) FE. Own PPML-HDFE solver, checked against
  pyfixest on the main nights. Clusters: origin province (G = 22).
- Inference: CRV1 with t(G-1) and the wild cluster restricted score bootstrap
  (Kline & Santos 2012) with Webb six-point weights, N_BOOT draws; CIs by test inversion.
- Difference holiday - ordinary: stacked two-night model with night-specific fixed
  effects and covariates and a common roots slope plus roots x holiday; H0: delta = 0.
- QAP: doubly-constrained gravity null (rows, columns, distance-bin totals, IPF),
  residual log((W+1)/(E+1)); province 47x47 design with simultaneous row/column
  permutations; city design with within-pair demeaning and permutation of cities
  within their province. Partial QAP = Dekker, Krackhardt & Snijders (2007) double
  semi-partialling (DSP). Difference QAP: (Z_holiday - Z_ordinary) vs Z_roots.
"""
from __future__ import annotations

import sys
import time
import zlib
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sps
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import b00_common as B0  # noqa: E402

T0 = time.time()
NET = C.RESULTS / "b01_networks"
OUT = C.RESULTS
EXT = C.ROOT / "_local" / "external"
BIRTH_MUN_CSV = EXT / "ine_birth_residence" / "mun_residence_x_prov_birth_censo_2021_2025.csv"
BIRTH_PROV_CSV = EXT / "ine_birth_residence" / "prov_birth_x_prov_residence.csv"
ADRH_CSV = EXT / "ine_adrh" / "adrh_municipios.csv"
MIG_CSV = EXT / "ine_interprovincial_migration" / "interprovincial_migration_1988_2024.csv"
SCI_CSV = EXT / "meta_sci_spain" / "sci_spain_nuts3.csv"
TT_P84 = C.ROOT / "explore" / "out_p84_city_zone_tt.parquet"   # road time city -> zone (legacy p84 output, read only)

ROOTS_YEAR = 2024
INCOME_YEAR = 2024
INCOME_IND = "Renta neta media por persona"
N_BOOT = 9999
ALPHA = 0.05
N_DBIN_CITY = 20   # quantile distance bins, gravity null city x province (as legacy p93)
N_DBIN_PROV = 15   # quantile distance bins, gravity null province x province (as legacy p93)
WEBB = np.array([-np.sqrt(1.5), -1.0, -np.sqrt(0.5), np.sqrt(0.5), 1.0, np.sqrt(1.5)])


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def rng_for(tag: str) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([C.SEED, zlib.crc32(tag.encode())]))


# =============================================================================
# 1. Nodes, nights, city set, roots, covariates
# =============================================================================
nodes = B0.load_nodes()
saved_nodes = pd.read_parquet(NET / "nodes.parquet")
assert (saved_nodes["id"].values == nodes["id"].values).all(), "b01 nodes differ from b00.load_nodes()"
assert nodes["population"].notna().all()
NID = nodes["id"].values
NPROV = nodes["prov"].values
NCOAST = nodes["coastal_true"].values.astype(bool)
PROVS = np.array(sorted(nodes["prov"].unique()))
assert len(PROVS) == 47
P2I = {p: k for k, p in enumerate(PROVS)}
NPI = np.array([P2I[p] for p in NPROV])
COASTAL_PROVS = sorted(nodes.loc[nodes["coastal_true"], "prov"].unique())
INLAND_PROVS = [p for p in PROVS if p not in COASTAL_PROVS]
assert len(COASTAL_PROVS) == 19 and len(INLAND_PROVS) == 28

# Nights. Convention (config): file date d = the night ENDING on d (d-1 -> d); main nights are holiday eves.
NW = pd.read_csv(NET / "networks.csv", dtype=str)
_dates24 = B0.available_dates(C.PERNOCT_2024)
_feb = sorted(NW.loc[(NW.group == "ordinary") & (NW.role == "sister"), "date"])
assert _feb == [d for d in _dates24 if d[5:7] == "02" and pd.Timestamp(d).day_of_week in (1, 2, 3)]
SISTERS = {
    "ordinary": _feb,                                              # Tue/Wed/Thu files of Feb 2024 (13)
    "summer": list(C.SISTER_NIGHTS["summer"]),
    "allsaints": list(C.SISTER_NIGHTS["allsaints"]),
    "christmas": [C.APPENDIX_NIGHTS["christmas"]],
    "goodfriday": [C.APPENDIX_NIGHTS["goodfriday"]],
    # robustness baseline (pattern #12): every Tue/Wed/Thu file of Oct 2024 present in the data (15)
    "october": [d for d in _dates24 if d[5:7] == "10" and pd.Timestamp(d).day_of_week in (1, 2, 3)],
}
YEARS = {g: list(v) for g, v in C.YEAR_REPLICATION.items()}
MAIN = {**dict(C.MAIN_NIGHTS), **dict(C.APPENDIX_NIGHTS)}
HOLI = ["summer", "allsaints", "christmas", "goodfriday"]
GROUP = {}
for g_, ds in list(SISTERS.items()) + list(YEARS.items()):
    for d_ in ds:
        GROUP[d_] = g_
ALL_DATES = sorted(GROUP)
assert len(SISTERS["ordinary"]) == 13 and len(SISTERS["summer"]) == 9 and len(SISTERS["october"]) == 15
assert all(MAIN[g] in SISTERS[g] for g in ["ordinary"] + HOLI)

# ---- roots: city x birth province (2024) --------------------------------------
Bm = pd.read_csv(BIRTH_MUN_CSV, dtype={"ine_code": str, "prov_birth": str})
Bm = Bm[Bm["year"] == ROOTS_YEAR]
assert Bm["persons"].notna().all()
cities_151 = sorted(Bm["ine_code"].unique())
assert len(cities_151) == 151
CITIES = np.array([c for c in cities_151 if c[:2] not in C.EXCLUDED_PROVINCES])
assert len(CITIES) == 137 and set(CITIES) <= set(NID), "mainland census cities must be MITMA zones"
CNAME = Bm.drop_duplicates("ine_code").set_index("ine_code")["municipality_name"]
CPROV = np.array([c[:2] for c in CITIES])
CPI = np.array([P2I[p] for p in CPROV])
pcount = Counter(CPROV)
ID_PROVS = sorted(p for p, k in pcount.items() if k >= 2)
IDROW = np.array([pcount[p] >= 2 for p in CPROV])
assert len(ID_PROVS) == 22 and IDROW.sum() == 112
ID2CITYROW = {c: k for k, c in enumerate(CITIES)}
CITY_NODE = np.array([np.where(NID == c)[0][0] for c in CITIES])
NODE2ROW = -np.ones(len(NID), dtype=int)
NODE2ROW[CITY_NODE] = np.arange(len(CITIES))

tot_res = Bm.groupby("ine_code")["persons"].sum()          # all birthplaces incl. abroad
bm = Bm[Bm["prov_birth"].isin(PROVS)].set_index(["ine_code", "prov_birth"])["persons"]
BORN = np.array([[bm.get((c, p), np.nan) for p in PROVS] for c in CITIES], dtype=float)
assert np.isfinite(BORN).all()
OWN = CPROV[:, None] == PROVS[None, :]                      # 137 x 47 own-province cells
MASK_C = ~OWN
assert (BORN[MASK_C] > 0).all(), "zero birthplace cell: log share undefined"
TOT = tot_res.reindex(CITIES).values.astype(float)
SHARE = BORN / TOT[:, None]

# ---- roots: province x province (2024, censo anual) -----------------------------
Bp = pd.read_csv(BIRTH_PROV_CSV, dtype={"prov_birth": str, "prov_residence": str})


def prov_birth_matrix(year: int, source: str) -> np.ndarray:
    s = Bp[(Bp.year == year) & (Bp.source == source)].set_index(["prov_residence", "prov_birth"])["persons"]
    M = np.array([[s.get((r, b), np.nan) for b in PROVS] for r in PROVS], dtype=float)   # [residence, birth]
    assert np.isfinite(M).all()
    return M


RP = prov_birth_matrix(ROOTS_YEAR, "censo_anual_68523")     # RP[r, s] = born in s living in r
RP98 = prov_birth_matrix(1998, "padron_continuo")
MASK_P = ~np.eye(47, dtype=bool)

# ---- distances ----------------------------------------------------------------------
pw = nodes["population"].values.astype(float)
PCX = np.array([np.average(nodes["x"].values[NPI == k], weights=pw[NPI == k]) for k in range(47)])
PCY = np.array([np.average(nodes["y"].values[NPI == k], weights=pw[NPI == k]) for k in range(47)])
CX = nodes["x"].values[CITY_NODE]
CY = nodes["y"].values[CITY_NODE]
KM_C = np.hypot(CX[:, None] - PCX[None, :], CY[:, None] - PCY[None, :]) / 1000.0
KM_P = np.hypot(PCX[:, None] - PCX[None, :], PCY[:, None] - PCY[None, :]) / 1000.0
assert (KM_C[MASK_C] > 1).all() and (KM_P[MASK_P] > 1).all()

# ---- road time (pop-weighted mean of p84 tt_adj over the zones of p) -------------
tt = pd.read_parquet(TT_P84, columns=["c", "z", "tt_adj"])
tt = tt[tt["z"].isin(set(NID)) & tt["c"].isin(set(CITIES))]
zpos = pd.Series(np.arange(len(NID)), index=NID)
tt["pi"] = NPI[zpos.loc[tt["z"]].values]
tt["pop"] = pw[zpos.loc[tt["z"]].values]
tt["ri"] = tt["c"].map(ID2CITYROW)
tt["wtt"] = tt["tt_adj"] * tt["pop"]
g = tt.groupby(["ri", "pi"])[["wtt", "pop"]].sum()
TT_C = np.full((len(CITIES), 47), np.nan)
TT_C[g.index.get_level_values(0), g.index.get_level_values(1)] = (g["wtt"] / g["pop"]).values
assert np.isfinite(TT_C[MASK_C]).all()

# ---- income (ADRH) ---------------------------------------------------------------------
A = pd.read_csv(ADRH_CSV, dtype={"ine_code": str})
inc = A[(A.year == INCOME_YEAR) & (A.indicator == INCOME_IND)].set_index("ine_code")["value"]
INC = inc.reindex(CITIES).values.astype(float)
assert np.isfinite(INC).all()
LINC = np.log(INC) - np.log(INC).mean()

pd.DataFrame({
    "ine_code": CITIES, "name": CNAME.reindex(CITIES).values, "prov": CPROV, "identifying": IDROW,
    "residents_2024": TOT, "income_" + str(INCOME_YEAR): INC,
}).to_csv(OUT / "r2_city_set.csv", index=False)
rows = []
for a, c in enumerate(CITIES):
    for b, p in enumerate(PROVS):
        rows.append((c, p, OWN[a, b], BORN[a, b], TOT[a], SHARE[a, b]))
pd.DataFrame(rows, columns=["ine_code", "prov_birth", "own_province", "born", "residents", "share"]).to_csv(
    OUT / "r2_roots_city_prov.csv", index=False)
pd.DataFrame([(r, s, RP[a, b]) for a, r in enumerate(PROVS) for b, s in enumerate(PROVS)],
             columns=["prov_residence", "prov_birth", "persons"]).to_csv(OUT / "r2_roots_prov_prov.csv", index=False)
log(f"city set: {len(CITIES)} mainland cities, {IDROW.sum()} identifying in {len(ID_PROVS)} provinces")

# =============================================================================
# 2. Flows per night (from b01 full edge lists)
# =============================================================================
manifest = pd.read_csv(NET / "manifest.csv", dtype={"date": str}).drop_duplicates("date").set_index("date")


BUILT_HERE = []


def load_edges(date: str) -> pd.DataFrame:
    """b01 full edge list if b01 saved it; otherwise built with the same b00 loader (b00.build_edges)."""
    f = NET / f"edges_{date}__full.parquet"
    if f.exists():
        e = pd.read_parquet(f)
        assert abs(e["w"].sum() - manifest.loc[date, "total_weight_full"]) < 1e-6 * e["w"].sum()
        return e
    BUILT_HERE.append(date)
    e, _ = B0.build_edges(date, nodes)
    return e


for d in C.MAIN_NIGHTS.values():   # the saved b01 edge lists are exactly what the b00 loader builds
    e_saved = load_edges(d)
    assert d not in BUILT_HERE, f"{d} expected in b01"
    e_new, _ = B0.build_edges(d, nodes)
    assert len(e_saved) == len(e_new) and (e_saved[["i", "j"]].values == e_new[["i", "j"]].values).all()
    assert np.allclose(e_saved["w"].values, e_new["w"].values)
log("b01 edge lists match b00.build_edges on the three main nights")

FLOWS: dict[str, dict[str, np.ndarray]] = {}
for d in ALL_DATES:
    e = load_edges(d)
    i, j, w = e["i"].values, e["j"].values, e["w"].values
    pi, pj = NPI[i], NPI[j]
    P = np.zeros((47, 47))
    np.add.at(P, (pi, pj), w)
    np.fill_diagonal(P, 0.0)
    r = NODE2ROW[i]
    sel = r >= 0
    out = {"prov": P}
    for kind, keep in (("all", np.ones(len(j), bool)), ("inland", ~NCOAST[j]), ("coastal", NCOAST[j])):
        M = np.zeros((len(CITIES), 47))
        s2 = sel & keep
        np.add.at(M, (r[s2], pj[s2]), w[s2])
        M[OWN] = 0.0
        out[kind] = M
    M = out["all"].copy()
    M[:, [P2I[p] for p in COASTAL_PROVS]] = 0.0
    out["inlandprov"] = M
    FLOWS[d] = out
# Noise-matched ordinary references: mean of the ordinary sister nights, divided by the number of nights
# (audit C3; every pair is summed over all nights with missing = 0). Robustness references only.
POOL = "2024-02-mean13"
OCTPOOL = "2024-10-mean15"
for name, grp in ((POOL, "ordinary"), (OCTPOOL, "october")):
    ds = SISTERS[grp]
    FLOWS[name] = {k: sum(FLOWS[d][k] for d in ds) / len(ds) for k in FLOWS[ALL_DATES[0]]}
    GROUP[name] = grp
POOLS = [POOL, OCTPOOL]
log("flows built for", len(FLOWS), "nights incl. 2 means; built here with b00.build_edges:", BUILT_HERE)


# =============================================================================
# 3. PPML with high-dimensional fixed effects (own solver; checked vs pyfixest)
# =============================================================================
class HDFE:
    """Weighted alternating projections over several fixed-effect factors."""

    def __init__(self, codes_list):
        self.D = []
        for c in codes_list:
            c = np.asarray(c)
            n, L = len(c), int(c.max()) + 1
            self.D.append(sps.csr_matrix((np.ones(n), (np.arange(n), c)), shape=(n, L)))

    def demean(self, X, w, tol=1e-12, maxit=20000):
        X = np.array(X, dtype=float, copy=True)
        one = X.ndim == 1
        if one:
            X = X[:, None]
        if X.shape[1] == 0:
            return X
        scale = max(1.0, float(np.abs(X).max()))
        ws = [D.T @ w for D in self.D]
        for _ in range(maxit):
            delta = 0.0
            for D, s in zip(self.D, ws):
                m = (D.T @ (w[:, None] * X)) / s[:, None]
                X -= D @ m
                delta = max(delta, float(np.abs(m).max()))
            if delta < tol * scale:
                return X[:, 0] if one else X
        raise RuntimeError("demeaning did not converge")


def poisson_dev(y, mu):
    t = np.where(y > 0, y * np.log(np.where(y > 0, y, 1.0) / mu), 0.0)
    return float(2.0 * np.sum(t - (y - mu)))


def ppml(y, X, hd, offset=None, mu0=None, maxit=300):
    n, k = len(y), X.shape[1]
    off = np.zeros(n) if offset is None else offset
    mu = (y + y.mean()) / 2.0 if mu0 is None else mu0.copy()
    eta = np.log(mu)
    dev_old, beta = np.inf, np.zeros(k)
    for _ in range(maxit):
        z = eta - off + (y - mu) / mu
        Dm = hd.demean(np.column_stack([z, X]), mu)
        zt, Xt = Dm[:, 0], Dm[:, 1:]
        if k:
            XtW = Xt.T * mu
            b_new = np.linalg.solve(XtW @ Xt, XtW @ zt)
            e = zt - Xt @ b_new
        else:
            b_new, e = beta, zt
        eta = z - e + off
        mu = np.exp(eta)
        dev = poisson_dev(y, mu)
        db = float(np.max(np.abs(b_new - beta))) if k else 0.0
        beta = b_new
        if abs(dev - dev_old) <= 1e-12 * (abs(dev) + 0.1) and db < 1e-10:
            return beta, mu, dev
        dev_old = dev
    raise RuntimeError("PPML did not converge")


def drop_zero_groups(y, codes_list):
    """Iteratively drop observations in fixed-effect groups whose outcome is all zero (separation)."""
    keep = np.ones(len(y), bool)
    while True:
        bad = np.zeros(len(y), bool)
        for c in codes_list:
            s = np.bincount(c[keep], y[keep], minlength=c.max() + 1)
            bad |= keep & (s[c] <= 0)
        if not bad.any():
            return keep
        keep &= ~bad


# ---- design for the neighbouring-city model ---------------------------------------
ri_all, pj_all = np.where(MASK_C & IDROW[:, None])
CELL = pd.DataFrame({"ri": ri_all, "pj": pj_all})
CELL["city"] = CITIES[ri_all]
CELL["oprov"] = CPROV[ri_all]
CELL["dprov"] = PROVS[pj_all]
X_LSHARE = np.log(SHARE[ri_all, pj_all])
X_LBORN1 = np.log1p(BORN[ri_all, pj_all])
LKM = np.log(KM_C[ri_all, pj_all])
LTT = np.log(TT_C[ri_all, pj_all])
lkm_c, ltt_c = LKM - LKM.mean(), LTT - LTT.mean()
INC_REF = PROVS[0]
INC_COLS = [p for p in PROVS if p != INC_REF]
XINC = np.column_stack([LINC[ri_all] * (pj_all == P2I[p]) for p in INC_COLS])
COV = {
    "main": np.column_stack([lkm_c, lkm_c ** 2, XINC]),
    "road": np.column_stack([ltt_c, ltt_c ** 2, XINC]),
    "noinc": np.column_stack([lkm_c, lkm_c ** 2]),
    "legacy": np.zeros((len(ri_all), 0)),
}
XVAR = {"main": X_LSHARE, "road": X_LSHARE, "noinc": X_LSHARE, "legacy": X_LBORN1}
CL_LEVELS = np.array(ID_PROVS)
CL_ALL = pd.Categorical(CELL["oprov"], categories=CL_LEVELS).codes.astype(int)
CITYCL_ALL = pd.factorize(CELL["city"])[0]
G = len(CL_LEVELS)
VB = rng_for("r2_webb").choice(WEBB, size=(N_BOOT, G))       # one fixed draw for every test


class Night:
    """One night, one outcome sample, one specification: data after separation drops."""

    def __init__(self, date, sample="all", spec="main", drop_prov=None):
        self.date, self.sample, self.spec = date, sample, spec
        y = FLOWS[date][sample][ri_all, pj_all]
        c_code = pd.factorize(CELL["city"])[0]
        p_code = pd.factorize(CELL["oprov"] + "_" + CELL["dprov"])[0]
        keep = drop_zero_groups(y, [c_code, p_code])
        if drop_prov is not None:                      # leave-one-origin-province-out
            keep &= (CELL["oprov"].values != drop_prov)
        self.keep = keep
        self.y = y[keep]
        self.x = XVAR[spec][keep]
        cov = COV[spec][keep]
        self.fe = [pd.factorize(c_code[keep])[0], pd.factorize(p_code[keep])[0]]
        self.hd = HDFE(self.fe)
        # drop covariates that are zero or collinear with the fixed effects / each other
        # (e.g. income slopes of destinations removed by separation; structural rank is weight-free)
        if cov.shape[1]:
            ct = self.hd.demean(cov, np.ones(len(self.x)))
            from scipy.linalg import qr
            _, Rq, piv = qr(ct, mode="economic", pivoting=True)
            dg = np.abs(np.diag(Rq))
            ok = np.sort(piv[dg > 1e-8 * dg[0]])
            cov = cov[:, ok]
        self.cov = cov
        xt = self.hd.demean(self.x, np.ones(len(self.x)))
        if cov.shape[1]:
            ct = self.hd.demean(cov, np.ones(len(self.x)))
            res = xt - ct @ np.linalg.lstsq(ct, xt, rcond=None)[0]
            assert np.sum(res ** 2) > 1e-6 * np.sum(xt ** 2), "roots regressor collinear with controls"
        self.cl = (CITYCL_ALL if spec == "legacy" else CL_ALL)[keep]
        self.n = int(keep.sum())


def fit_unrestricted(nt: Night):
    X = np.column_stack([nt.x, nt.cov])
    beta, mu, dev = ppml(nt.y, X, nt.hd)
    Xt = nt.hd.demean(X, mu)
    H = (Xt.T * mu) @ Xt
    Hinv = np.linalg.inv(H)
    sc = Xt * (nt.y - mu)[:, None]
    assert np.abs(sc.sum(0)).max() < 1e-6 * (1 + np.abs(nt.y).sum()), "score not zero at the MLE"
    Gc = int(nt.cl.max()) + 1
    Sg = np.vstack([np.bincount(nt.cl, sc[:, q], minlength=Gc) for q in range(X.shape[1])]).T
    IF = Sg @ Hinv[:, 0]                                    # cluster influence of beta_roots
    Geff = int((np.abs(IF) > 0).sum())
    se = float(np.sqrt(Geff / (Geff - 1) * np.sum(IF ** 2)))
    return {"beta": float(beta[0]), "se": se, "IF": IF, "G": Geff, "mu": mu, "dev": dev}


def score_p(y, mu, target, Xother, hd, cl):
    """Wild cluster restricted score bootstrap (Kline & Santos 2012), Webb weights."""
    tt_ = hd.demean(target, mu)
    if Xother.shape[1]:
        Xt = hd.demean(Xother, mu)
        XtW = Xt.T * mu
        tt_ = tt_ - Xt @ np.linalg.solve(XtW @ Xt, XtW @ tt_)
    s = tt_ * (y - mu)
    Sg = np.bincount(cl, s, minlength=G)
    t0 = Sg.sum() / np.sqrt(np.sum(Sg ** 2))
    tb = (VB @ Sg) / np.sqrt((VB ** 2) @ (Sg ** 2))
    return float(np.mean(np.abs(tb) >= abs(t0))), float(t0)


def p_beta(nt: Night, b0: float, mu0=None):
    beta, mu, _ = ppml(nt.y, nt.cov, nt.hd, offset=b0 * nt.x, mu0=mu0)
    p, _ = score_p(nt.y, mu, nt.x, nt.cov, nt.hd, nt.cl)
    return p, mu


def invert_ci(pfun, est, se, alpha=ALPHA, tol_frac=2e-3):
    """Bounds of {b0 : p(b0) > alpha} by bracketing and bisection on each side of est."""
    out = []
    cache = {"mu": None}

    def pf_(b):
        p, mu = pfun(b, cache["mu"])
        cache["mu"] = mu
        return p

    for sgn in (-1, 1):
        cache["mu"] = None
        inner, step = est, se
        outer = est + sgn * step
        k = 0
        while pf_(outer) > alpha:
            inner, step = outer, step * 1.6
            outer = est + sgn * step
            k += 1
            if k > 30:
                raise RuntimeError("CI bracket failed")
        while abs(outer - inner) > tol_frac * se:
            mid = 0.5 * (inner + outer)
            if pf_(mid) > alpha:
                inner = mid
            else:
                outer = mid
        out.append(0.5 * (inner + outer))
    return out[0], out[1]


class Stacked:
    """Two nights stacked: night-specific FE and covariates, common roots slope + roots x holiday."""

    def __init__(self, nh: Night, no: Night):
        assert nh.spec == no.spec and nh.sample == no.sample
        self.y = np.r_[no.y, nh.y]
        n0, n1 = no.n, nh.n
        fe = []
        for a, b in zip(no.fe, nh.fe):
            fe.append(np.r_[a, b + a.max() + 1])
        self.hd = HDFE(fe)
        k0, k1 = no.cov.shape[1], nh.cov.shape[1]
        cov = np.zeros((n0 + n1, k0 + k1))
        cov[:n0, :k0] = no.cov
        cov[n0:, k0:] = nh.cov
        self.x = np.r_[no.x, nh.x]
        self.target = np.r_[np.zeros(n0), nh.x]               # roots x holiday
        self.Xother = np.column_stack([self.x, cov])
        self.cl = np.r_[no.cl, nh.cl]

    def p_delta(self, d0, mu0=None):
        _, mu, _ = ppml(self.y, self.Xother, self.hd, offset=d0 * self.target, mu0=mu0)
        p, _ = score_p(self.y, mu, self.target, self.Xother, self.hd, self.cl)
        return p, mu

    def fit_unrestricted(self):
        X = np.column_stack([self.target, self.Xother])
        beta, mu, _ = ppml(self.y, X, self.hd)
        return float(beta[0]), float(beta[1])


def tcrit(Geff):
    return stats.t.ppf(1 - ALPHA / 2, Geff - 1)


# ---- solver check against pyfixest -----------------------------------------------
def check_pyfixest():
    import pyfixest as pf
    rows = []
    for d in MAIN.values():
        nt = Night(d)
        fit = fit_unrestricted(nt)
        df = pd.DataFrame({"y": nt.y, "x": nt.x, "c": nt.fe[0], "pair": nt.fe[1], "cl": nt.cl})
        names = []
        for q in range(nt.cov.shape[1]):
            df[f"v{q}"] = nt.cov[:, q]
            names.append(f"v{q}")
        m = pf.fepois("y ~ x + " + " + ".join(names) + " | c + pair", data=df, vcov={"CRV1": "cl"})
        t = m.tidy()
        b_pf, se_pf = float(t.loc["x", "Estimate"]), float(t.loc["x", "Std. Error"])
        assert abs(b_pf - fit["beta"]) < 1e-6, (d, b_pf, fit["beta"])
        assert m._N == nt.n
        rows.append({"date": d, "beta_own": fit["beta"], "beta_pyfixest": b_pf, "se_own_crv1": fit["se"],
                     "se_pyfixest_crv1": se_pf, "n": nt.n})
        log(f"pyfixest check {d}: beta own {fit['beta']:.6f} pf {b_pf:.6f}; se own {fit['se']:.4f} pf {se_pf:.4f}")
    pd.DataFrame(rows).to_csv(OUT / "r2_solver_check.csv", index=False)


check_pyfixest()

# =============================================================================
# 4. Single-night PPML: every night, samples x specs
# =============================================================================
SAMPLES = ["all", "inland", "coastal", "inlandprov"]
RUNS = [("all", "main"), ("all", "road"), ("all", "noinc"), ("all", "legacy"),
        ("inland", "main"), ("coastal", "main"), ("inlandprov", "main")]
YEAR_DATES = set(sum(YEARS.values(), []))
CI_NIGHTS = set(MAIN.values()) | YEAR_DATES | set(POOLS)


def role_of(d):
    if d in POOLS:
        return "mean"
    if d in YEAR_DATES:
        return "year"
    if d in MAIN.values():
        return "main"
    return "sister"


NIGHT_FITS: dict = {}
rows = []
for sample, spec in RUNS:
    for d in ALL_DATES + POOLS:
        nt = Night(d, sample, spec)
        f = fit_unrestricted(nt)
        NIGHT_FITS[(sample, spec, d)] = (nt, f)
        rec = {"sample": sample, "spec": spec, "date": d, "group": GROUP[d], "role": role_of(d),
               "beta": f["beta"], "se_crv1": f["se"], "G": f["G"], "n_obs": nt.n,
               "n_cities": len(np.unique(nt.fe[0])), "flow_total": float(nt.y.sum())}
        tc = tcrit(f["G"])
        rec["ci_t_lo"], rec["ci_t_hi"] = f["beta"] - tc * f["se"], f["beta"] + tc * f["se"]
        if spec != "legacy":
            rec["p_wcr_beta0"], _ = p_beta(nt, 0.0)
            if spec == "main" and (d in MAIN.values() or d in POOLS or (sample == "all" and d in CI_NIGHTS)):
                lo, hi = invert_ci(lambda b, m0: p_beta(nt, b, m0), f["beta"], f["se"])
                rec["ci_wcr_lo"], rec["ci_wcr_hi"] = lo, hi
        rows.append(rec)
    log(f"single nights done: {sample}/{spec}")
NIGHTS = pd.DataFrame(rows)
NIGHTS.to_csv(OUT / "r2_ppml_nights.csv", index=False)


def summarise(df, value, keys):
    out = []
    for key, sub in df.groupby(keys):
        key = key if isinstance(key, tuple) else (key,)
        for grp in ["ordinary", "october"] + HOLI:
            s = sub[sub["date"].isin(SISTERS[grp])]
            if s.empty:
                continue
            mn = MAIN.get(grp)
            v = s[value].values
            rec = {**dict(zip(keys, key)), "group": grp, "main_night": mn or "", "n_sisters": len(v),
                   "median_sisters": float(np.median(v)), "min_sisters": float(v.min()), "max_sisters": float(v.max())}
            if mn is not None and (s["date"] == mn).any():
                main = float(s.loc[s["date"] == mn, value].iloc[0])
                rec["main"] = main
                if len(v) > 1:
                    rec["main_is_extreme"] = bool((main == v.max() and (v == v.max()).sum() == 1)
                                                  or (main == v.min() and (v == v.min()).sum() == 1))
            for k_, val in sub[sub["date"].isin(YEARS.get(grp, []))].set_index("date")[value].to_dict().items():
                rec[f"year_{k_}"] = val
            for pl in POOLS:
                if GROUP[pl] == grp and (sub["date"] == pl).any():
                    rec["mean_night_value"] = float(sub.loc[sub["date"] == pl, value].iloc[0])
            out.append(rec)
    return pd.DataFrame(out)


summarise(NIGHTS, "beta", ["sample", "spec"]).to_csv(OUT / "r2_ppml_summary.csv", index=False)

# =============================================================================
# 5. Difference holiday - ordinary (stacked model)
# =============================================================================
# baseline "feb": the 13 Feb sister nights; "oct": the 15 Oct Tue/Wed/Thu nights (pattern #12)
pairs = [(h, o, "feb") for grp in HOLI for h in SISTERS[grp] for o in SISTERS["ordinary"]]
pairs += [(MAIN[grp], o, "oct") for grp in HOLI for o in SISTERS["october"]]
pairs += [(MAIN[grp], POOL, "feb_mean") for grp in HOLI] + [(MAIN[grp], OCTPOOL, "oct_mean") for grp in HOLI]
pairs += [("2022-08-15", "2022-02-16", "year"), ("2022-11-01", "2022-02-16", "year"), ("2023-08-15", "2023-02-15", "year")]
pairs += [(OCTPOOL, POOL, "placebo")]          # ordinary vs ordinary: October mean minus February mean
MAIN_PAIRS = {(MAIN[g], b) for g in HOLI for b in (MAIN["ordinary"], POOL, OCTPOOL)}
MAIN_PAIRS |= {("2022-08-15", "2022-02-16"), ("2022-11-01", "2022-02-16"), ("2023-08-15", "2023-02-15"), (OCTPOOL, POOL)}
CI_PAIRS_ALLSAMPLES = {(MAIN[g], b) for g in HOLI for b in (MAIN["ordinary"], OCTPOOL)} | {(OCTPOOL, POOL)}
WCR_ALL_PAIRS = {("all", "main"), ("inland", "main"), ("coastal", "main")}

# one-off check: the unrestricted stacked fit reproduces beta_o and beta_h - beta_o
nh_, fh_ = NIGHT_FITS[("all", "main", MAIN["summer"])]
no_, fo_ = NIGHT_FITS[("all", "main", MAIN["ordinary"])]
dlt, bo = Stacked(nh_, no_).fit_unrestricted()
assert abs(bo - fo_["beta"]) < 1e-6 and abs(dlt - (fh_["beta"] - fo_["beta"])) < 1e-6
log(f"stacked check: delta {dlt:.6f} = {fh_['beta']:.6f} - {fo_['beta']:.6f}")

rows = []
for sample, spec in RUNS:
    if spec == "legacy":
        continue
    for h, o, base in pairs:
        nh, fh = NIGHT_FITS[(sample, spec, h)]
        no, fo = NIGHT_FITS[(sample, spec, o)]
        dIF = fh["IF"] - fo["IF"]
        Geff = int((np.abs(dIF) > 0).sum())
        se = float(np.sqrt(Geff / (Geff - 1) * np.sum(dIF ** 2)))
        delta = fh["beta"] - fo["beta"]
        tc = tcrit(Geff)
        rec = {"sample": sample, "spec": spec, "holiday": h, "ordinary": o, "baseline": base, "group": GROUP[h],
               "is_main_pair": (h, o) in MAIN_PAIRS, "beta_holiday": fh["beta"], "beta_ordinary": fo["beta"],
               "delta": delta, "se_crv1": se, "G": Geff, "ci_t_lo": delta - tc * se, "ci_t_hi": delta + tc * se,
               "p_t": float(2 * stats.t.sf(abs(delta / se), Geff - 1))}
        if (sample, spec) in WCR_ALL_PAIRS or (h, o) in MAIN_PAIRS:
            st = Stacked(nh, no)
            rec["p_wcr_delta0"], _ = st.p_delta(0.0)
            if spec == "main" and (h, o) in MAIN_PAIRS and (sample == "all" or (h, o) in CI_PAIRS_ALLSAMPLES):
                lo, hi = invert_ci(st.p_delta, delta, se)
                rec["ci_wcr_lo"], rec["ci_wcr_hi"] = lo, hi
        rows.append(rec)
    log(f"differences done: {sample}/{spec}")
DIFF = pd.DataFrame(rows)
DIFF.to_csv(OUT / "r2_diff_tests.csv", index=False)

drows = []
for (sample, spec), sub in DIFF.groupby(["sample", "spec"]):
    for grp in HOLI:
        for base in ("feb", "oct"):
            s = sub[(sub["group"] == grp) & (sub["baseline"] == base)]
            if s.empty:
                continue
            ref = MAIN["ordinary"] if base == "feb" else None
            v = s["delta"].values
            rec = {"sample": sample, "spec": spec, "group": grp, "baseline": base, "n_pairs": len(s),
                   "median": float(np.median(v)), "min": float(v.min()), "max": float(v.max()),
                   "n_positive": int((v > 0).sum()), "share_p_t_lt_0.05": float((s["p_t"] < 0.05).mean())}
            if ref is not None:
                m = s[(s["holiday"] == MAIN[grp]) & (s["ordinary"] == ref)]
                rec["delta_main_vs_feb14"] = float(m["delta"].iloc[0])
            if s["p_wcr_delta0"].notna().all() if "p_wcr_delta0" in s else False:
                rec["share_p_wcr_lt_0.05"] = float((s["p_wcr_delta0"] < 0.05).mean())
                rec["median_p_wcr"] = float(s["p_wcr_delta0"].median())
                rec["min_p_wcr"] = float(s["p_wcr_delta0"].min())
                rec["max_p_wcr"] = float(s["p_wcr_delta0"].max())
            drows.append(rec)
pd.DataFrame(drows).to_csv(OUT / "r2_diff_summary.csv", index=False)

# leave-one-origin-province-out jackknife of the main differences (vs Feb 14), all and inland destinations
jrows = []
for sample in ("all", "inland"):
    fits = {}
    for pr in [None] + ID_PROVS:
        for d in [MAIN["ordinary"]] + [MAIN[g] for g in HOLI]:
            fits[(pr, d)] = fit_unrestricted(Night(d, sample, "main", drop_prov=pr))["beta"]
    for g in HOLI:
        full = fits[(None, MAIN[g])] - fits[(None, MAIN["ordinary"])]
        for pr in ID_PROVS:
            jrows.append({"sample": sample, "group": g, "dropped_origin_prov": pr,
                          "n_cities_dropped": int((CPROV == pr).sum()), "delta_full": full,
                          "delta_without": fits[(pr, MAIN[g])] - fits[(pr, MAIN["ordinary"])]})
pd.DataFrame(jrows).to_csv(OUT / "r2_diff_jackknife.csv", index=False)
log("jackknife done")

# =============================================================================
# 6. QAP under a doubly-constrained gravity null
# =============================================================================
def qbins(D, mask, nb):
    e = np.quantile(D[mask], np.linspace(0, 1, nb + 1))
    e[0], e[-1] = -1.0, e[-1] + 1.0
    b = np.digitize(D, e)
    b[~mask] = 0
    return b


def gravity_ipf(W, bins, mask, tol=1e-9, maxit=20000):
    """E_ij = a_i b_j f_bin(ij) on mask; matches row sums, column sums and distance-bin totals."""
    Wm = np.where(mask, W, 0.0)
    r, c = Wm.sum(1), Wm.sum(0)
    nb = int(bins.max()) + 1
    Wb = np.bincount(bins[mask], Wm[mask], minlength=nb)
    a, b, f = np.ones(len(r)), np.ones(len(c)), np.where(np.arange(nb) > 0, 1.0, 0.0)
    for it in range(maxit):
        F = np.where(mask, f[bins], 0.0)
        a = np.where(r > 0, r / np.maximum(F @ b, 1e-300), 0.0)
        b = np.where(c > 0, c / np.maximum(F.T @ a, 1e-300), 0.0)
        E = np.outer(a, b) * F
        Eb = np.bincount(bins[mask], E[mask], minlength=nb)
        f = np.where(Eb > 0, f * Wb / np.maximum(Eb, 1e-300), 0.0)
        E = np.outer(a, b) * np.where(mask, f[bins], 0.0)
        er = np.max(np.abs(E.sum(1) - r) / np.maximum(r, 1e-12))
        ec = np.max(np.abs(E.sum(0) - c) / np.maximum(c, 1e-12))
        eb = np.max(np.abs(np.bincount(bins[mask], E[mask], minlength=nb) - Wb)[1:] / np.maximum(Wb[1:], 1e-12))
        if max(er, ec, eb) < tol:
            return E, it + 1, max(er, ec, eb)
    raise RuntimeError("gravity IPF did not converge")


GRAV_LOG = []


def grav_resid(W, bins, mask, label):
    E, it, err = gravity_ipf(W, bins, mask)
    GRAV_LOG.append({"matrix": label, "iterations": it, "max_rel_margin_error": err})
    Z = np.zeros_like(W, dtype=float)
    Z[mask] = np.log((W[mask] + 1.0) / (E[mask] + 1.0))
    return Z


BIN_C = qbins(KM_C, MASK_C, N_DBIN_CITY)
BIN_P = qbins(KM_P, MASK_P, N_DBIN_PROV)
ROOTS_C = np.where(MASK_C, BORN, 0.0)
ROOTS_P = np.where(MASK_P, RP, 0.0)

# city design: 112 identifying cities, residuals demeaned within (origin prov x dest prov)
IDR = np.where(IDROW)[0]
MASK_ID = MASK_C[IDR]
GROUPS_ID = [np.where(CPROV[IDR] == p)[0] for p in ID_PROVS]


def demean_pairs(Z):
    Z = Z[IDR].copy()
    for gi in GROUPS_ID:
        sub = Z[gi]
        m = sub.mean(0, keepdims=True)
        Z[gi] = np.where(MASK_ID[gi], sub - m, 0.0)
    Z[~MASK_ID] = 0.0
    return Z


def centre(Z, mask):
    Z = np.where(mask, Z - Z[mask].mean(), 0.0)
    return Z


ZR_C = demean_pairs(grav_resid(ROOTS_C, BIN_C, MASK_C, "roots_city"))
ZR_P = centre(grav_resid(ROOTS_P, BIN_P, MASK_P, "roots_prov"), MASK_P)

rq = rng_for("r2_qap")
PERM_C = np.empty((C.N_QAP_PERM, len(IDR)), dtype=int)
for b in range(C.N_QAP_PERM):
    perm = np.arange(len(IDR))
    for gi in GROUPS_ID:
        perm[gi] = rq.permutation(gi)
    PERM_C[b] = perm
PERM_P = np.vstack([rq.permutation(47) for _ in range(C.N_QAP_PERM)])
RC_PERM = ZR_C[PERM_C]                                # (B, 112, 47); mask invariant under these perms
RP_PERM = ZR_P[PERM_P[:, :, None], PERM_P[:, None, :]]
assert np.all((RC_PERM != 0) <= MASK_ID[None]) and np.all(RP_PERM[:, ~MASK_P] == 0)


def qap_stats(Y, Cc, R, Rperm, only_diff=False):
    """Return dict of (obs, null array) for: simple r(Y,R); DSP partial r(Y,R|C); r(Y-C, R)."""
    def corr_null(A):
        na = np.sqrt(np.sum(A * A))
        nr = np.sqrt(np.sum(R * R))
        obs = np.sum(A * R) / (na * nr)
        null = np.einsum("bij,ij->b", Rperm, A) / (na * nr)
        return obs, null

    if only_diff:
        return {"r_diff": corr_null(Y - Cc)}
    out = {"r": corr_null(Y), "r_diff": corr_null(Y - Cc) if Cc is not None else None}
    if Cc is not None:
        cc = np.sum(Cc * Cc)
        ry = Y - (np.sum(Y * Cc) / cc) * Cc
        eR = R - (np.sum(R * Cc) / cc) * Cc
        # DSP: permute the residual of roots on the control (Dekker, Krackhardt & Snijders 2007)
        eRp = eR[PERM_C] if R.shape[0] == len(IDR) else eR[PERM_P[:, :, None], PERM_P[:, None, :]]
        num = np.einsum("bij,ij->b", eRp, ry)
        den = np.sqrt(np.sum(ry * ry) * (np.sum(eR * eR) - np.einsum("bij,ij->b", eRp, Cc) ** 2 / cc))
        out["r_partial_dsp"] = (np.sum(ry * eR) / np.sqrt(np.sum(ry * ry) * np.sum(eR * eR)), num / den)
    return out


ZC, ZP = {}, {}
for d in ALL_DATES + POOLS:
    ZC[d] = demean_pairs(grav_resid(FLOWS[d]["all"], BIN_C, MASK_C, f"city_{d}"))
    ZP[d] = centre(grav_resid(FLOWS[d]["prov"], BIN_P, MASK_P, f"prov_{d}"), MASK_P)
pd.DataFrame(GRAV_LOG).to_csv(OUT / "r2_gravity_convergence.csv", index=False)
log("gravity residuals done; max margin error", max(g["max_rel_margin_error"] for g in GRAV_LOG))


def controls_for(d):
    """(control label, control night) pairs for a holiday night."""
    if d[:4] != "2024":
        return [("feb_same_year", {"2022": "2022-02-16", "2023": "2023-02-15"}[d[:4]])]
    return [("feb14", MAIN["ordinary"]), ("feb_mean13", POOL), ("oct_mean15", OCTPOOL)]


def qrec(design, d, ctrl_label, ctrl, stat, obs, null):
    return {"design": design, "date": d, "group": GROUP[d], "control": ctrl_label, "control_night": ctrl,
            "stat": stat, "obs": float(obs), "null_mean": float(null.mean()), "null_sd": float(null.std()),
            "z": float((obs - null.mean()) / null.std()),
            "p_upper": float((1 + np.sum(null >= obs)) / (len(null) + 1)), "n_perm": len(null)}


qrows = []
for design, Zd, R, Rperm in (("city_within_province", ZC, ZR_C, RC_PERM), ("province_47", ZP, ZR_P, RP_PERM)):
    for d in ALL_DATES + POOLS:
        obs, null = qap_stats(Zd[d], None, R, Rperm)["r"]
        qrows.append(qrec(design, d, "", "", "r", obs, null))
        if GROUP[d] in HOLI and d not in POOLS:
            for lab, o in controls_for(d):
                res = qap_stats(Zd[d], Zd[o], R, Rperm)
                for stat in ("r_diff", "r_partial_dsp"):
                    qrows.append(qrec(design, d, lab, o, stat, *res[stat]))
    for grp in HOLI:
        for h in SISTERS[grp]:
            for o in SISTERS["ordinary"]:
                qrows.append(qrec(design, h, "feb_pairs", o, "r_diff_pairs",
                                  *qap_stats(Zd[h], Zd[o], R, Rperm, only_diff=True)["r_diff"]))
        for o in SISTERS["october"]:
            qrows.append(qrec(design, MAIN[grp], "oct_pairs", o, "r_diff_pairs",
                              *qap_stats(Zd[MAIN[grp]], Zd[o], R, Rperm, only_diff=True)["r_diff"]))
QAP = pd.DataFrame(qrows)
QAP.to_csv(OUT / "r2_qap.csv", index=False)
qs = []
for (design, stat, ctrl), sub in QAP.groupby(["design", "stat", "control"]):
    for grp in ["ordinary", "october"] + HOLI:
        s = sub[sub["date"].isin(SISTERS[grp])]
        if s.empty:
            continue
        mn = MAIN.get(grp)
        m = s[s["date"] == mn] if mn else s.iloc[0:0]
        if stat == "r_diff_pairs":
            m = m[m["control_night"] == MAIN["ordinary"]] if ctrl == "feb_pairs" else m.iloc[0:0]
        qs.append({"design": design, "stat": stat, "control": ctrl, "group": grp, "n": len(s),
                   "main": float(m["obs"].iloc[0]) if len(m) else np.nan,
                   "p_main": float(m["p_upper"].iloc[0]) if len(m) else np.nan,
                   "median": float(s["obs"].median()), "min": float(s["obs"].min()), "max": float(s["obs"].max()),
                   "max_p": float(s["p_upper"].max()), "share_p_lt_0.05": float((s["p_upper"] < 0.05).mean())})
    for pl in POOLS:
        s = sub[sub["date"] == pl]
        if len(s):
            qs.append({"design": design, "stat": stat, "control": ctrl, "group": pl, "n": 1,
                       "main": float(s["obs"].iloc[0]), "p_main": float(s["p_upper"].iloc[0])})
    for d in YEAR_DATES:
        s = sub[sub["date"] == d]
        if len(s):
            qs.append({"design": design, "stat": stat, "control": ctrl, "group": f"year_{d}", "n": 1,
                       "main": float(s["obs"].iloc[0]), "p_main": float(s["p_upper"].iloc[0])})
pd.DataFrame(qs).to_csv(OUT / "r2_qap_summary.csv", index=False)
log("QAP done")

# =============================================================================
# 7. Legacy reproduction (p72 / p82: 151 cities, window means, log1p(born), c + pairFE, CRV1 by city)
# =============================================================================
def legacy_reproduction():
    import duckdb
    import geopandas as gpd
    import pyfixest as pf
    gids = set(gpd.read_file(C.ZONES_GPKG, ignore_geometry=True)["id"])
    cities = sorted(set(Bm["ine_code"]) & gids)
    assert len(cities) == 151
    con = duckdb.connect()
    cl = ",".join(f"'{c}'" for c in cities)
    win = {"Feb weekday": "date BETWEEN '2024-02-05' AND '2024-02-29' AND dayofweek(date) IN (2,3,4)",
           "mid-Aug": "date BETWEEN '2024-08-10' AND '2024-08-18'"}
    bb = Bm.set_index(["ine_code", "prov_birth"])["persons"]
    grid = pd.MultiIndex.from_product([cities, [f"{i:02d}" for i in range(1, 51)]], names=["c", "p"]).to_frame(index=False)
    grid = grid[grid.c.str[:2] != grid.p].copy()
    grid["lborn"] = np.log1p(bb.reindex(list(zip(grid.c, grid.p))).fillna(0).values)
    grid["pairFE"] = grid.c.str[:2] + "_" + grid.p
    rows = []
    for lab, cond in win.items():
        nn = con.execute(f"SELECT count(DISTINCT date) FROM read_parquet('{C.PERNOCT_2024}') WHERE {cond}").fetchone()[0]
        for how, denom in (("legacy sum/count(distinct date) [C3 bug]", "count(DISTINCT date)"),
                           ("correct sum/n_nights", str(nn))):
            f = con.execute(f"SELECT r c, substr(s,1,2) p, sum(p)/{denom} w FROM read_parquet('{C.PERNOCT_2024}') "
                            f"WHERE {cond} AND r IN ({cl}) AND r<>s GROUP BY 1,2").df().set_index(["c", "p"])["w"]
            gg = grid.copy()
            gg["w"] = f.reindex(list(zip(gg.c, gg.p))).fillna(0).values
            t = pf.fepois("w ~ lborn | c + pairFE", data=gg, vcov={"CRV1": "c"}).tidy()
            rows.append({"step": f"L0/L1 151 cities, 50 dest provinces, {how}", "slice": lab, "n_nights": nn,
                         "beta": float(t.loc["lborn", "Estimate"]), "se": float(t.loc["lborn", "Std. Error"]),
                         "se_type": "CRV1 city"})
    # L2: 137 mainland cities, 47 mainland provinces, correct window mean from the b01 nights
    for lab, dates in (("Feb weekday", [d for d in SISTERS["ordinary"] if d >= "2024-02-05"]), ("mid-Aug", SISTERS["summer"])):
        frames = []
        for d in dates:
            M = FLOWS[d]["all"]
            a, b = np.where(MASK_C)
            frames.append(pd.DataFrame({"date": d, "r": CITIES[a], "s": PROVS[b], "p": M[a, b]}))
        wm = B0.window_mean_from_frame(pd.concat(frames), n_nights=len(dates)).set_index(["r", "s"])["p"]
        gg = pd.DataFrame({"c": CITIES[ri_all], "p": PROVS[pj_all], "lborn": X_LBORN1,
                           "pairFE": CELL["oprov"].values + "_" + CELL["dprov"].values})
        gg["w"] = wm.reindex(list(zip(gg.c, gg.p))).values
        t = pf.fepois("w ~ lborn | c + pairFE", data=gg, vcov={"CRV1": "c"}).tidy()
        rows.append({"step": "L2 137 mainland cities, 47 mainland provinces, correct window mean", "slice": lab,
                     "n_nights": len(dates), "beta": float(t.loc["lborn", "Estimate"]),
                     "se": float(t.loc["lborn", "Std. Error"]), "se_type": "CRV1 city"})
    for lab, d in (("Feb 14", MAIN["ordinary"]), ("Aug 15", MAIN["summer"]), ("Nov 1", MAIN["allsaints"])):
        _, f = NIGHT_FITS[("all", "legacy", d)]
        rows.append({"step": "L3 137 cities, single night, legacy spec", "slice": lab, "n_nights": 1,
                     "beta": f["beta"], "se": f["se"], "se_type": "CRV1 city"})
        _, f = NIGHT_FITS[("all", "main", d)]
        rows.append({"step": "L4 main spec (log share, log km + sq, income x dest, CRV1 origin prov)", "slice": lab,
                     "n_nights": 1, "beta": f["beta"], "se": f["se"], "se_type": "CRV1 origin province"})
    pd.DataFrame(rows).to_csv(OUT / "r2_legacy_reproduction.csv", index=False)
    log("legacy reproduction done")


legacy_reproduction()

# =============================================================================
# 8. Appendix: old vs recent migration and Meta SCI (province 47x47, single nights)
# =============================================================================
def echo_sci():
    import pyfixest as pf
    mig = pd.read_csv(MIG_CSV, dtype={"prov_origin": str, "prov_dest": str})
    old = mig[(mig.year >= 1988) & (mig.year <= 1995) & (mig.source == "EVR_microdata")]
    rec = mig[(mig.year >= 2021) & (mig.year <= 2024) & (mig.source == "EMCR_69742")]
    assert old.year.nunique() == 8 and rec.year.nunique() == 4
    old = old.groupby(["prov_origin", "prov_dest"])["migrations"].sum()
    rec = rec.groupby(["prov_origin", "prov_dest"])["migrations"].sum()
    sci = pd.read_csv(SCI_CSV, dtype={"user_prov_ine": str, "friend_prov_ine": str})
    sci = sci[sci.user_prov_ine.isin(PROVS) & sci.friend_prov_ine.isin(PROVS)]   # islands are several NUTS3
    assert len(sci) == 47 * 47
    sci = sci.set_index(["user_prov_ine", "friend_prov_ine"])["scaled_sci"]
    a, b = np.where(MASK_P)
    df = pd.DataFrame({"rp": PROVS[a], "sp": PROVS[b]})
    df["dbin"] = BIN_P[a, b]
    # "in" = toward the origin r from the destination s (people of s living in r; moves s -> r)
    df["l_birth24_in"] = np.log(RP[a, b])
    df["l_birth24_out"] = np.log(RP[b, a])
    df["l_birth98_in"] = np.log(RP98[a, b])
    df["l_mig_old_in"] = np.log1p(old.reindex(list(zip(df.sp, df.rp))).fillna(0).values)
    df["l_mig_rec_in"] = np.log1p(rec.reindex(list(zip(df.sp, df.rp))).fillna(0).values)
    df["l_mig_rec_out"] = np.log1p(rec.reindex(list(zip(df.rp, df.sp))).fillna(0).values)
    df["l_sci"] = np.log(sci.reindex(list(zip(df.rp, df.sp))).values)
    assert df.notna().all().all()
    specs = {"birth24": "l_birth24_in", "old": "l_mig_old_in", "recent": "l_mig_rec_in",
             "old_vs_recent": "l_mig_old_in + l_mig_rec_in",
             "birth98_ctrl_recent": "l_birth98_in + l_mig_rec_in + l_mig_rec_out",
             "sci": "l_sci", "birth_vs_sci": "l_birth24_in + l_birth24_out + l_sci"}
    rows = []
    dates = sorted(set().union(*[SISTERS[g] for g in ["ordinary"] + HOLI]))
    for dest in ("all", "inland"):
        for d in dates:
            g_ = df.copy()
            g_["y"] = FLOWS[d]["prov"][a, b]
            if dest == "inland":
                g_ = g_[g_.sp.isin(INLAND_PROVS)]
            m0 = pf.fepois("y ~ 1 | rp + sp + dbin", data=g_, vcov={"CRV1": "rp"})
            dev0 = float(np.sum(m0.deviance))
            for name, rhs in specs.items():
                m = pf.fepois(f"y ~ {rhs} | rp + sp + dbin", data=g_, vcov={"CRV1": "rp"})
                t = m.tidy()
                rec_ = {"dest": dest, "date": d, "group": GROUP[d], "spec": name, "n": m._N,
                        "dev_explained_beyond_fe": 1 - float(np.sum(m.deviance)) / dev0}
                for v in t.index:
                    rec_[f"b_{v}"] = float(t.loc[v, "Estimate"])
                    rec_[f"se_{v}"] = float(t.loc[v, "Std. Error"])
                rows.append(rec_)
    E = pd.DataFrame(rows)
    E.to_csv(OUT / "r2_echo_sci_nights.csv", index=False)
    vals = [c for c in E.columns if c.startswith("b_") or c == "dev_explained_beyond_fe"]
    summ = []
    for (dest, spec), sub in E.groupby(["dest", "spec"]):
        for grp in ["ordinary"] + HOLI:
            s = sub[sub["date"].isin(SISTERS[grp])]
            r_ = {"dest": dest, "spec": spec, "group": grp, "n_sisters": len(s)}
            for v in vals:
                if s[v].notna().any():
                    r_[f"{v}_main"] = float(s.loc[s["date"] == MAIN[grp], v].iloc[0])
                    r_[f"{v}_median"] = float(s[v].median())
                    r_[f"{v}_min"] = float(s[v].min())
                    r_[f"{v}_max"] = float(s[v].max())
            summ.append(r_)
    pd.DataFrame(summ).to_csv(OUT / "r2_echo_sci_summary.csv", index=False)
    log("echo / SCI appendix done")


echo_sci()


# =============================================================================
# 9. Reconciliation with the descriptive roots index of pattern mining (results/p_patterns.md #7)
# =============================================================================
N_BOOT_INDEX = 1999


def index_moments(M, sample):
    """Per-cluster second moments for the roots index, two designs.

    'between' = p_patterns.roots_index: 137 cities x other provinces, y = log1p(flow), x = log1p(born),
    city means removed, log(km + 1) partialled out; clusters = 47 origin provinces.
    'within'  = same variables, but city AND (origin prov x dest prov) means removed (two-way, unweighted):
    only the 112 identifying cities; clusters = 22 origin provinces. Moments are additive over clusters.
    """
    cols = np.ones(47, bool) if sample != "coastal" else np.isin(PROVS, COASTAL_PROVS)
    if sample == "inlandprov":
        cols = np.isin(PROVS, INLAND_PROVS)
    out = {}
    for design in ("between", "within"):
        provs = sorted(set(CPROV)) if design == "between" else ID_PROVS
        mom = np.zeros((len(provs), 6))     # Sxx Syy Sll Sxy Sxl Syl
        for k, pr in enumerate(provs):
            rows = np.where(CPROV == pr)[0]
            msk = MASK_C[rows][:, cols]
            Y = np.log1p(M[rows][:, cols]); X = np.log1p(BORN[rows][:, cols]); L = np.log(KM_C[rows][:, cols] + 1)
            if design == "between":
                vals = []
                for Z in (X, Y, L):
                    Z = np.where(msk, Z, np.nan)
                    vals.append(np.where(msk, Z - np.nanmean(Z, 1, keepdims=True), 0.0))
            else:
                ri, cj = np.where(msk)
                hd = HDFE([pd.factorize(ri)[0], pd.factorize(cj)[0]])
                D = hd.demean(np.column_stack([X[ri, cj], Y[ri, cj], L[ri, cj]]), np.ones(len(ri)))
                vals = [D[:, 0], D[:, 1], D[:, 2]]
            x, y, l_ = (v.ravel() for v in vals)
            mom[k] = [x @ x, y @ y, l_ @ l_, x @ y, x @ l_, y @ l_]
        out[design] = mom
    return out


def index_from(m):
    Sxx, Syy, Sll, Sxy, Sxl, Syl = m
    cxy, cxx, cyy = Sxy - Sxl * Syl / Sll, Sxx - Sxl ** 2 / Sll, Syy - Syl ** 2 / Sll
    return cxy / np.sqrt(cxx * cyy), cxy / cxx, np.sqrt(cyy)


IDX_DATES = ALL_DATES + POOLS
MOM = {(smp, d): index_moments(FLOWS[d][smp], smp) for smp in ("all", "inland", "coastal") for d in IDX_DATES}
irows = []
for (smp, d), mm in MOM.items():
    rec = {"sample": smp, "date": d, "group": GROUP[d], "role": role_of(d)}
    for design in ("between", "within"):
        r_, b_, sdy = index_from(mm[design].sum(0))
        rec[f"corr_{design}"], rec[f"ols_slope_{design}"], rec[f"sd_resid_logflow_{design}"] = r_, b_, sdy
    smp_ppml = {"all": "all", "inland": "inland", "coastal": "coastal"}[smp]
    rec["ppml_beta_main"] = NIGHT_FITS[(smp_ppml, "main", d)][1]["beta"]
    irows.append(rec)
IDX = pd.DataFrame(irows)
# validation: the 'between' correlation is the pattern-mining index (p_daily_metrics.roots_corr_partial)
pdm = pd.read_csv(OUT / "p_daily_metrics.csv").set_index("date")["roots_corr_partial"]
chk = IDX[(IDX["sample"] == "all") & IDX["date"].isin(pdm.index)]
dev_ = np.abs(chk["corr_between"].values - pdm.loc[chk["date"]].values).max()
log(f"roots index reproduces p_daily_metrics on {len(chk)} nights, max abs diff {dev_:.2e}")
assert dev_ < 1e-6
IDX.to_csv(OUT / "r2_index_bridge.csv", index=False)

# difference tests for the index: cluster (origin-province) pairs bootstrap, same draws for both nights
rb = rng_for("r2_index_boot")
BOOT_IDX = {des: rb.integers(0, n_, size=(N_BOOT_INDEX, n_)) for des, n_ in (("between", 47), ("within", len(ID_PROVS)))}
assert len(set(CPROV)) == 47
brows = []
idx_pairs = [(h, o, b) for h, o, b in pairs if b != "year"]
for smp in ("all", "inland", "coastal"):
    for h, o, base in idx_pairs:
        for design in ("between", "within"):
            mh, mo = MOM[(smp, h)][design], MOM[(smp, o)][design]
            rh, _, _ = index_from(mh.sum(0)); ro, _, _ = index_from(mo.sum(0))
            B_ = BOOT_IDX[design]
            cnt = np.stack([np.bincount(b, minlength=mh.shape[0]) for b in B_])     # (B, clusters)
            bh = index_from((cnt @ mh).T)[0]
            bo = index_from((cnt @ mo).T)[0]
            dd = bh - bo
            diff = rh - ro
            p2 = float(min(1.0, 2 * min(np.mean(dd <= 0), np.mean(dd >= 0))))
            brows.append({"sample": smp, "design": design, "holiday": h, "ordinary": o, "baseline": base,
                          "group": GROUP[h], "index_holiday": rh, "index_ordinary": ro, "diff": diff,
                          "ci_lo": float(np.quantile(dd, 0.025)), "ci_hi": float(np.quantile(dd, 0.975)),
                          "p_boot_two_sided": p2})
IB = pd.DataFrame(brows)
IB.to_csv(OUT / "r2_index_diff_tests.csv", index=False)
srows = []
for (smp, design, grp, base), s_ in IB.groupby(["sample", "design", "group", "baseline"]):
    srows.append({"sample": smp, "design": design, "group": grp, "baseline": base, "n_pairs": len(s_),
                  "median_diff": s_["diff"].median(), "min_diff": s_["diff"].min(), "max_diff": s_["diff"].max(),
                  "n_positive": int((s_["diff"] > 0).sum()), "share_p_lt_0.05": float((s_["p_boot_two_sided"] < 0.05).mean())})
pd.DataFrame(srows).to_csv(OUT / "r2_index_diff_summary.csv", index=False)
log("index reconciliation done")
log("R2 finished")

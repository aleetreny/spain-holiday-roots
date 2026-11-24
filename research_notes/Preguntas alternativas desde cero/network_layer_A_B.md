# Network-Science Layer for Candidates A (1833 Borders) and B (Rural Exodus & Holidays)

*Prepared for UCL COMP0123 (Complex Networks & Online Social Networks, Prof. Shi Zhou).*  
*New scripts: `explore/p90_A_basic_properties.py` – `explore/p94_check_2nd_order_and_tails.py`.*  
*New figures: `reports/figs/net_A_topology.png`, `reports/figs/net_A_modularity_race.png`, `reports/figs/net_B_topology_shift.png`, `reports/figs/net_B_qap_and_territories.png`.*

---

## Part 1 — Question A: Do the 1833 Provincial Borders Still Shape Spain's Mobility Network?

### 1.1 Basic Properties of the Mainland Mobility Network (`N = 2,484`)

All metrics are computed on the **peninsular mainland Spain municipal network** ($N = 2,484$ MITMA zones, excluding foreign NUTS3, islands, and Ceuta/Melilla; matching `P_dcgrav_wk.npy`) filtered to persistent links ($\text{persist} \ge 0.50$, retaining $98.5\text{–}99.2\%$ of all inter-zone trips while keeping density below Prof. Zhou's $5\%$ guideline on weekdays and Sundays). Rewired nulls use **30 degree-preserving rewirings** ($10E$ edge swaps each); the spatial null uses **20 binary graph draws** calibrated to the same edge count from the doubly-constrained gravity matrix $P^{\text{grav}}_{\text{wk}}$; and the 1st-order-mixing-preserving null uses **10 joint-degree-distribution (JDD) rewirings** ($6E$ within-bin swaps, preserving both $P(k)$ and $r$).

| Property | Weekday (`net_wk`) | Sunday (`net_sun`) | August (`net_aug`) | Weekday Nulls (ER / 30× Rew / 20× Spatial) |
|---|---:|---:|---:|---|
| Nodes $N$ (giant comp. $|GCC|$) | 2,484 (2,433) | 2,484 (2,433) | 2,484 (2,434) | — |
| Edges $E_{\text{und}}$ ($E_{\text{dir}}$) | 122,374 (233,163) | 146,869 (256,685) | 168,145 (301,014) | Same $E_{\text{und}}$ |
| Density $\rho_E$ | **3.97%** | **4.76%** | 5.45% | — |
| Trip volume retained (share) | 55.40M (99.2%) | 42.14M (98.5%) | 53.21M (99.0%) | — |
| Mean / median / 99% / max $k$ | 98.5 / 78 / 415 / 2,053 | 118.3 / 91 / 499 / 2,288 | 135.4 / 107 / 543 / 2,244 | Preserved in Rew |
| Degree CV / Gini / Freeman cent. | 0.865 / 0.357 / 0.788 | 0.888 / 0.367 / 0.875 | 0.828 / 0.351 / 0.850 | — |
| Mean / median / max strength $s$ | 44.6k / 19.6k / 3.45M | 33.9k / 16.4k / 2.40M | 42.8k / 20.6k / 2.66M | — |
| Reciprocity binary / weighted | 0.950 / 0.990 | 0.856 / 0.959 | 0.883 / 0.984 | — |
| Avg local clustering $C_{\text{avg}}$ ($T_{\text{glob}}$) | **0.688** (0.443) | **0.691** (0.429) | **0.665** (0.428) | ER: 0.040 · Rew: 0.116 ± 0.0002 · **Spatial: 0.427 ± 0.001** |
| Avg path length $L$ (diameter $D$) | **2.146** (3) | **2.016** (3) | **2.009** (3) | ER: 1.980 · Rew: 1.978 ± 0.001 |
| Small-world index $\sigma_{\text{ER}}$ ($\sigma_{\text{rew}}$) | 15.99 (5.49) | 14.07 (4.59) | 11.80 (4.31) | $\sigma > 1$ strongly small-world |
| 1st-order degree assortativity $r$ | **−0.026** | **−0.053** | **−0.041** | Rew: **−0.061 ± 0.001** · Spatial: +0.006 ± 0.002 |
| Weighted strength assortativity $r_s^w$ ($r_s^u$) | +0.077 (−0.028) | +0.069 (−0.037) | +0.077 (−0.034) | — |
| **2nd-order mixing $R_{\text{avg}}$ (Zhou)** | **+0.597** | **+0.381** | **+0.447** | Rew: **−0.117** · JDD-rew: **+0.073** · **Spatial: +0.634 ± 0.008** |
| **2nd-order mixing $R_{\text{max}}$ (Zhou)** | **+0.090** | +0.041 | +0.059 | Rew: −0.016 · JDD-rew: +0.015 · Spatial: +0.271 ± 0.011 |
| Max $k$-core (size; % Catalonia) | **98** (151; **98.7% CAT**) | **130** (157; 98.7% CAT) | **129** (153; 98.7% CAT) | Rew: **60.3 ± 0.5** · **Spatial: 76.9 ± 0.7** |
| Norm. rich-club $\rho(k)$: $k=100/200/300/400$ | **1.15 / 1.23 / 1.34 / 1.41** | 1.08 / 1.14 / 1.17 / 1.32 | 1.04 / 1.20 / 1.25 / 1.37 | 30 rewirings; $z(400)=+11.5$, $\phi(400)=0.954$ |
| Opsahl weighted rich-club $\rho_w(k)$ ($k=100/300/400$) | 1.56 / 8.98 / 13.74 | — | — | Weight-reshuffled null |
| Spatial-null strength rich-club $\rho_{\text{spatial}}$ (top 5/10/50/100) | 0.99 / 1.20 / 1.03 / 1.04 | — | — | **No strength rich-club beyond gravity** |

**Key complex-network findings (directly using Prof. Zhou's methods):**
1. **Heavy-tailed but not purely scale-free CCDF:** On the weekday network, Vuong log-likelihood-ratio tests on the upper tail ($k \ge 100$, $n=829$, $\hat{\alpha}=3.21$; $s \ge 30\text{k}$, $n=840$, $\hat{\alpha}=2.16$) modestly favour a power-law tail over a truncated log-normal ($z=5.0$ at $k\ge 100$, dropping to $z=1.86, p=0.063$ at $k\ge 200$), while the bulk is unimodal ($\text{median}=78$, $\text{CV}=0.865$) with a dramatic top-city hub spike ($k_{\max}=2,053$, Madrid linked to $83\%$ of all mainland zones).
2. **Structural constraints on assortativity (Zhou & Mondragón 2007):** Observed weekday assortativity $r=-0.026$ looks "neutral", yet degree-preserving rewiring yields $r_{\text{rew}}=-0.061 \pm 0.001$ (hub degree heterogeneity structurally forces slight disassortativity). Relative to its degree-constrained baseline, Spain's mobility network is actually **slightly assortative** ($\Delta r = +0.035$), and shifts more disassortative on Sundays ($-0.053$) and August ($-0.041$) as hub reach expands ($k_{\max}=2,288$).
3. **Second-order mixing as a signature of spatial embedding (Zhou, Cox & Hansen 2017):** Weekday $R_{\text{avg}} = +0.597$ collapses under both degree-preserving rewiring ($-0.117 \pm 0.002$) and **1st-order JDD-preserving rewiring** ($r=-0.028$, $R_{\text{avg}} = +0.073 \pm 0.002$), confirming it is genuinely independent of 1st-order mixing. Crucially, the **spatial configuration null reproduces $R_{\text{avg}} = +0.634 \pm 0.008$ completely!** Thus strong positive second-order mixing in mobility networks is explained by spatial embedding (neighbours share the same local density basin), whereas clustering ($0.688$ vs spatial null $0.427$) and the Catalan $k$-core ($98$ vs spatial null $76.9$) still reflect genuine meso-scale structure beyond gravity.
4. **Rich-club dissection:** Binary $\rho(k)$ rises to $1.41$ at $k>400$ ($z=11.5$, 30 rewirings) and Opsahl weighted $\rho_w(k)$ reaches $13.7$, yet under the doubly-constrained gravity null the strength rich-club disappears ($\rho_{\text{spatial}} = 0.99\text{–}1.04$ across top 5–500 nodes, except $1.20$ at top 10 due to AVE/air hub corridors).
5. **BOSAM fingerprint (`reports/figs/net_A_topology.png`, panel d):** Sorting nodes by weekday degree produces a crisp asymmetric fingerprint showing August's mid- and low-degree densification relative to the tighter weekday core.

---

### 1.2 Network-Native Modularity Test: 1833 Provinces vs Historical Maps & 1,000 Contiguous Placebos

We evaluate weighted modularity $Q = \frac{1}{2W} \sum_{i,j} (A_{ij} - P_{ij})\,\delta(g_i, g_j)$ for six partitions under three null models: **(i) Configuration Model** ($P^{\text{CM}}_{ij} = s_i s_j / 2W$), **(ii) Distance Spatial Configuration Null** ($P^{\text{dist}}_{ij} = a_i a_j f(d_{ij})$), and **(iii) Road Travel-Time Spatial Configuration Null** ($P^{\text{tt}}_{ij} = a_i a_j f(t_{ij})$ from OSM `road_tt_150km.parquet`). Each map is benchmarked against three placebo ensembles:
- **1,000 node-size-matched contiguous agglomerations** grown on the mainland polygon adjacency graph ($N=2,484, E=7,125$, 1 connected component) matching exact unit count $K$ and size distribution;
- **1,000 area-weighted contiguous partitions** (Voronoi from $K$ seeds sampled proportionally to zone area, matching geographic rather than node-density scale);
- **60 rigidly shifted & rotated copies** of the real map polygons (15 angles $\times$ 4 shifts).

| Partition (Units $K$, Cov.) | $Q_{\text{CM}}$ ($z_{\text{node}}$ / $z_{\text{area}}$ / $z_{\text{rot}}$) | $Q_{\text{dist}}$ ($z_{\text{node}}$ [$p$] / $z_{\text{area}}$ [$p$] / $z_{\text{rot}}$) | $Q_{\text{tt}}$ ($z_{\text{node}}$ [$p$] / $z_{\text{area}}$ [$p$] / $z_{\text{rot}}$ [$p$]) |
|---|---:|---:|---:|
| **1833 Provinces** ($K=47$, 100%) | **0.8659** (+9.09 / +4.74 / +5.48) | **0.0282** (**+3.13** [.004] / **+5.86** [.001] / **+6.64**) | **0.0174** (+0.53 [.301] / **+3.41** [.002] / **+3.61** [.033]) |
| **CCAA 1979–83** ($K=15$, 100%) | 0.8385 (+3.32 / +1.75 / +3.04) | 0.0144 (+2.38 [.020] / +3.77 [.001] / +3.96) | 0.0088 (+0.61 [.263] / +2.65 [.019] / +2.05 [.049]) |
| **Voronoi (47 capitals)** ($K=47$, 100%) | 0.8362 (+7.76 / +3.19 / +3.85) | 0.0053 (−0.89 [.816] / +1.10 [.130] / +0.97) | 0.0092 (−1.75 [.960] / +0.95 [.179] / +0.96 [.213]) |
| **1789 Floridablanca** ($K=28$, 69%) | 0.5375 (+2.78 / +0.11 / +2.48) | 0.0144 (+2.12 [.024] / +4.05 [.002] / +3.33) | 0.0072 (+0.85 [.219] / +2.67 [.011] / +1.82 [.082]) |
| **Partidos judiciales (1834)** ($K=401$, 85%) | 0.3998 (+5.64 / −4.78 / −4.43) | **0.0557** (+7.30 [.001] / +9.17 [.001] / +10.96) | **0.0331** (+4.71 [.001] / +4.21 [.001] / +4.43 [.016]) |
| **Dioceses** ($K=64$, 94%) | 0.7716 (+7.02 / +1.57 / +2.40) | **0.0350** (+3.82 [.001] / +6.65 [.001] / +6.73) | **0.0312** (+3.92 [.001] / **+7.47** [.001] / **+6.38** [.016]) |

**Regime stability (1833 Provinces / CCAA under Road-TT spatial null):**
- Weekday: $Q_{\text{prov}}^{\text{CM}} = 0.8659$, $Q_{\text{prov}}^{\text{tt}} = 0.0174$, $Q_{\text{ccaa}}^{\text{tt}} = 0.0088$
- Sunday: $Q_{\text{prov}}^{\text{CM}} = 0.8334$, $Q_{\text{prov}}^{\text{tt}} = 0.0148$, $Q_{\text{ccaa}}^{\text{tt}} = 0.0081$ (provincial modularity weakens on weekends!)
- August: $Q_{\text{prov}}^{\text{CM}} = 0.8543$, $Q_{\text{prov}}^{\text{tt}} = 0.0176$, $Q_{\text{ccaa}}^{\text{tt}} = 0.0091$

**Critical methodological discovery on spatial-null placebos:** Under the straight-line spatial null ($Q_{\text{dist}}$), 1833 Provinces beat **every** placebo cleanly ($z_{\text{node-contig}} = +3.13, p=0.004$; $z_{\text{area-contig}} = +5.86, p=0.001$; $z_{\text{rot}} = +6.64$), while the Voronoi of 47 provincial capitals fails completely ($z_{\text{node-contig}} = -0.89$, $z_{\text{area-contig}} = +1.10$, $z_{\text{rot}} = +0.97$). Under the **road travel-time null** ($Q_{\text{tt}}$), 1833 Provinces still strongly beat **area-weighted contiguous placebos** ($z = +3.41, p=0.002$) and **60 rotated/shifted maps** ($z = +3.61, p=0.033$), **but drop to $z = +0.53$ ($p = 0.301$) against uniform-node-seeded contiguous blobs!** Why? Diagnostic decomposition by distance band reveals an important spatial-null artefact: in a pure travel-time null, $0\text{–}10\text{ km}$ urban pairs have systematically *positive* residuals ($B_{\text{tt}} / 2W = +0.0363$) because urban speed limits make $5\text{ km}$ take as long as $15\text{ km}$ of rural highway, whereas $10\text{–}40\text{ km}$ pairs have negative residuals ($-0.0310$). Uniform-node-seeded contiguous blobs over-sample dense urban MITMA zones and form tiny $0\text{–}10\text{ km}$ metro pockets that harvest that positive urban-congestion residual—artificially inflating both node-seeded placebos ($0.0155$) and ultra-fine partitions (*Partidos judiciales*, $K=401$). This also confirms an honest nuance: **older fine-grained territorial units (*Dioceses*, $K=64$, and *Partidos judiciales*, $K=401$) achieve higher raw spatial modularity than 1833 Provinces**, while among macro-partitions the hierarchy is strictly **1833 Provinces $>$ CCAA $>$ 1789 Floridablanca $>$ Voronoi of capitals**.

---

### 1.3 Equivalence Between PPML Gravity Border Regressions and a Degree-Corrected Spatial Stochastic Block Model (DC-SSBM)

The Poisson Pseudo-Maximum Likelihood (PPML) gravity specification used to estimate border multipliers is not an external econometric add-on—it is mathematically identical to maximum-likelihood inference of a **Degree-Corrected Spatial Stochastic Block Model (DC-SSBM)** (Karrer & Newman 2011; Newman 2016; Bernaschi et al. 2022). Specifically, modelling expected trip counts between municipalities $i$ and $j$ as $\mathbb{E}[A_{ij}] = \theta_i^{\text{out}} \theta_j^{\text{in}} \, f(t_{ij}) \, \omega_{g_i, g_j}$, where $\theta_i$ are node propensities (origin and destination fixed effects), $f(t_{ij})$ is a non-parametric spatial deterrence kernel (2-minute travel-time bin fixed effects), and $\omega_{g_i, g_j} = \exp(\beta \, \mathbb{I}[g_i = g_j])$ is the planted partition's affinity matrix, yields a Poisson log-likelihood whose first-order conditions are solved identically by PPML. Under the null hypothesis $\beta = 0$, the MLE reduces exactly to the doubly-constrained spatial configuration model $P^{\text{tt}}_{ij} = \hat{\theta}_i \hat{\theta}_j \hat{f}(t_{ij})$ (Cazabet, Borgnat & Jensen 2017); moreover, taking the score (first derivative of the log-likelihood at $\beta = 0$) gives $\left.\frac{\partial \ell}{\partial \beta}\right|_{\beta=0} = \sum_{i,j} (A_{ij} - P^{\text{tt}}_{ij})\,\mathbb{I}[g_i = g_j] = 2W \, Q_{\text{tt}}(g)$—proving that **spatial-null modularity $Q_{\text{tt}}$ is the local score test of the DC-SSBM, while $\exp(\hat{\beta})$ is its exact finite-effect block affinity ratio** (Newman 2016). Framed this way, every border test in A (1833 vs 1789 new/old/erased borders, 1822 ghosts, language vs province) is a likelihood-based model selection over competing block structures on the same spatial network.

---

## Part 2 — Question B: Does the Rural Exodus Still Draw Spain's Holiday Map?

### 2.1 Three Networks on the Same Mainland Zones (`N = 2,484`) & Holiday Structural Transformation

We construct three network families on the exact same $N = 2,484$ mainland zones (plus the two birthplace-roots networks):
1. **Ordinary weekday overnight network:** Feb Tue–Thu mean (`Feb_wd_mean`) and single ordinary night Wed 14 Feb 2024 (`Feb_wd_1d`, apples-to-apples with single holiday nights).
2. **Holiday overnight networks:** Mid-August mean (Aug 10–18), peak summer night (Aug 15), All Saints (Nov 1 — a cemetery-visit, non-beach holiday), and Christmas Eve (Dec 24 — family homecoming).
3. **Birthplace "Roots" networks:** 47×47 mainland provinces (`prov_birth_x_prov_residence.csv`, 2024: 7.40M inter-provincial lifetime migrants) and 137 mainland Census cities × 47 provinces (`mun_residence_x_prov_birth_censo_2021_2025.csv`, 2024: 4.25M migrants).

Binary topological metrics use a consistent backbone threshold of $w \ge 10$ overnight stays/night so **every network stays cleanly below Prof. Zhou's $5\%$ density limit ($0.86\%\text{–}3.74\%$)** while capturing $64\%\text{–}83\%$ of single-night away-stays; weighted metrics (reach, $r_{\text{pop}}$, QAP) use the full unthresholded matrices. Rewirings use **25 degree-preserving swaps** ($10E$ steps).

| Property ($N = 2,484$, backbone $w \ge 10$) | Feb Wed Mean | **Feb Wed (Feb 14)** | **All Saints (Nov 1)** | **Christmas (Dec 24)** | **Aug 15 (Peak)** | Mid-Aug Mean | Roots 47×47 ($\ge 2\text{k}$) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Total away-nights (M) / flow in backbone | 3.79M (40%) | **2.11M (69%)** | **4.23M (76%)** | **4.45M (75%)** | **8.46M (83%)** | 11.04M (64%) | 7.40M migrants |
| Edges $E_{\text{und}}$ ($E_{\text{dir}}$) | 32,424 (42.9k) | **26,465** (35.8k) | **48,698** (64.2k) | **52,079** (68.9k) | **101,996** (116.0k) | 115,390 (131.1k) | 387 (2,162 all) |
| Density $\rho_E$ | 1.05% | **0.86%** | **1.58%** | **1.69%** | **3.31%** | 3.74% | 35.8% |
| Mean / median / max $k$ | 26.1 / 14 / 1,582 | **21.3 / 11 / 1,499** | **39.2 / 21 / 2,008** | **41.9 / 21 / 2,183** | **82.1 / 44.5 / 2,233** | 92.9 / 52 / 2,282 | — |
| Small inland village ($\text{pop}<5\text{k}$) mean $k_{\text{in}}$ | 5.8 | **4.1** | **10.5 (2.6×)** | **10.3 (2.5×)** | **23.4 (5.7×)** | 26.9 | — |
| Binary reciprocity | 0.488 | 0.522 | 0.484 | 0.487 | **0.242** | 0.239 | — |
| Clustering $C_{\text{avg}}$ (25× Rew) | 0.428 (0.209) | 0.543 (0.242) | 0.581 (0.293) | 0.587 (0.373) | 0.535 (0.380) | 0.511 (0.375) | 0.755 |
| Avg path length $L$ (diam; Rew) | 2.388 (4; 2.35) | 2.441 (4; 2.40) | 2.163 (4; 2.11) | 2.063 (4; 2.04) | **1.997** (3; 1.98) | 1.981 (3; 1.97) | — |
| **Degree assortativity $r$ (25× Rew)** | −0.111 (−0.122) | **−0.131** (−0.140) | **−0.155** (−0.161) | **−0.180** (−0.189) | **−0.214** (−0.217) | −0.213 (−0.218) | **−0.268** |
| **Directed mixing $r(k_{\text{out}}, k_{\text{in}})$** | −0.123 | **−0.135** | **−0.147** | **−0.167** | **−0.227** | −0.226 | — |
| **Weighted $\log(\text{pop})$ mixing $r_{\text{pop}}$** | **+0.069** | −0.020 | **−0.082** | −0.043 | −0.027 | +0.036 | **−0.229** |
| 2nd-order mixing $R_{\text{avg}}$ (25× Rew) | −0.205 (−0.217) | −0.232 (−0.245) | −0.269 (−0.279) | −0.312 (−0.323) | **−0.349** (−0.308) | −0.343 (−0.303) | — |
| **Normalised rich-club $\rho(100)$ / $\rho(200)$** | **1.25 / 1.28** | **1.20 / 1.19** | **1.06 / 1.13** | **1.11 / 1.15** | **0.99 / 1.02** | **0.99 / 1.03** | — |
| Max $k$-core (25× Rew) | 37 (30.3) | **37** (29.0) | **51** (45.0) | **65** (54.2) | **84** (82.9) | 92 (88.5) | — |
| Mean weighted reach (km) / $>100\text{ km}$ | 221 / 56.4% | 157 / 39.2% | 151 / 41.0% | 187 / 49.9% | **207 / 53.6%** | 239 / 60.6% | — |
| Cross-province share / into $<5\text{k}$ zones | 67.5% / 14.3% | 50.9% / 9.0% | 55.6% / 13.4% | 62.7% / 11.4% | **69.2% / 18.1%** | 74.6% / 19.7% | — |

**Four striking, strictly network-science structural shifts on holidays (`reports/figs/net_B_topology_shift.png`):**
1. **Disassortative phase shift & Zhou–Mondragón structural constraint:** Going from an ordinary February night to All Saints, Christmas Eve, and August 15 drives degree assortativity sharply more negative ($r: -0.131 \to -0.155 \to -0.180 \to -0.214$; directed $r(k_{\text{out}}, k_{\text{in}}): -0.135 \to -0.227$), moving the mobility graph directly toward the disassortative signature of the **Roots network** ($r=-0.268$, $r_{\text{pop}}=-0.229$). Even more remarkably, **degree-preserving rewiring reproduces $r$ to within $0.003\text{–}0.009$ in every holiday slice** ($r_{\text{rew}} = -0.217$ vs $-0.214$ on Aug 15): in Prof. Zhou's structural-constraint framework, **holiday disassortativity is structurally forced by the degree distribution** as urban hubs explode to $k_{\max}=2,233$ ($90\%$ of mainland Spain) and connect outward to low-degree villages! On All Saints (no beaches), weighted population mixing itself turns most negative ($r_{\text{pop}} = -0.082$ vs $+0.069$ Feb mean).
2. **Rich-club dissolution ($\rho(100): 1.25 \to 0.99$):** Ordinary weekday overnight stays exhibit a clear city-to-city business rich club ($\rho(100)=1.25$, $\rho(200)=1.28$). On holidays the rich club **dissolves completely** in mid-August ($\rho(100)=0.99$, $\rho(50)=1.00$) and weakens sharply on All Saints ($1.06$) and Christmas Eve ($1.11$), while reciprocity halves ($0.52 \to 0.24$) as flow turns one-way from cities to villages and coasts.
3. **$k$-core promotion of emptied inland villages (even at equal current population!):** Stratifying all 1,942 inland mainland zones ($\text{pop} \ge 200$) by quintile of 1900–2021 population loss reveals a dramatic network promotion of *España vaciada*:

| Inland Loss Quintile (1900–2021) | $n$ | Med. Pop | Med. Loss | $k$-core Feb 14 $\to$ Aug 15 ($\times$) | Norm. $k$-core Feb $\to$ Aug 15 | $k_{\text{in}}$ Feb $\to$ Aug 15 ($\times$) | $k_{\text{in}}$ Nov 1 / Dec 24 | Eff. Origin Provs (Aug 15) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Q1 (grew strongly) | 389 | 9,623 | −409% | 14.2 $\to$ 49.2 (**3.46×**) | 0.384 $\to$ 0.586 | 17.9 $\to$ 24.4 (**1.36×**) | 26.9 / 30.1 | 5.02 |
| Q2 | 388 | 8,454 | −84% | 12.2 $\to$ 42.0 (3.44×) | 0.330 $\to$ 0.500 | 24.1 $\to$ 39.2 (1.63×) | 37.1 / 42.2 | 6.94 |
| Q3 | 388 | 5,706 | −9% | 8.9 $\to$ 32.7 (3.70×) | 0.239 $\to$ 0.389 | 8.8 $\to$ 24.7 (2.82×) | 16.5 / 19.8 | 6.98 |
| Q4 | 388 | 5,046 | +34% | 7.9 $\to$ 33.0 (4.18×) | 0.213 $\to$ 0.392 | 6.2 $\to$ 30.0 (4.83×) | 15.2 / 16.6 | 7.49 |
| **Q5 (emptied most)** | 389 | 4,932 | **+65%** | **9.5 $\to$ 43.1 (4.55×)** | **0.256 $\to$ 0.513 (2.0×)** | **5.9 $\to$ 49.2 (8.33×)** | **18.3 / 17.2** | **7.87** |

   - **Size-matched control (small inland zones only, $\text{pop} < 6,000$, $n = 1,053$, nearly identical median pop $4,309\text{–}4,811$):** Extreme-loss villages (median $-65\%$) jump from **$k$-core $9.2 \to 41.5$** and **$k_{\text{in}}$ $5.6 \to 47.6$ ($8.5\times$)** on Aug 15 (and $k_{\text{in}} = 17.7$ on All Saints, $16.0$ on Christmas), compared to only **$k$-core $6.2 \to 22.1$** and **$k_{\text{in}}$ $4.4 \to 13.1$ ($3.0\times$)** for equal-sized growing towns! Emptied villages literally **double** their relative $k$-core depth and draw visitors from **7.8 effective origin provinces** vs 5.0 for equal-sized non-emptied towns.
4. **Population-sorted BOSAM (`reports/figs/net_B_topology_shift.png`, panel d):** On Feb 14 only the top-left city-to-city block is populated; on Aug 15 the top-ranked cities form a dark vertical "curtain" illuminating links all the way down to node 2,484.

---

### 2.2 Spatial-Null QAP Comparison & Bipartite "Summer Territories"

#### A. Quadratic Assignment Procedure (QAP) Under Doubly-Constrained Spatial Gravity Null (2,000 Permutations)

For every directed mobility matrix and reference matrix we fit a doubly-constrained spatial gravity null $E_{ij} = a_i b_j f(d_{ij})$ and compute log-ratio residuals $Z_{ij} = \log\frac{W_{ij}+1}{E_{ij}+1}$. We run 2,000 QAP permutations (simultaneous row/col permutations at the 47×47 province level; full bipartite AND **within-origin-province neighbouring-city permutations** at the 137 mainland Census cities × 47 provinces level), plus **Freedman–Lane partial QAP** controlling for the everyday February weekday spatial residual.

| Network / Level | Ordinary Feb Wed | Mid-Aug (Full / Excess) | All Saints (Full / Excess) | Christmas Eve (Full / Excess) |
|---|---|---|---|---|
| **47×47 Provinces vs Roots 2024:** | | | | |
| — Spatial-null QAP $r$ ($z$) | 0.573 ($z=16.3^{***}$) | 0.545 ($13.3^{***}$) / **0.117** ($2.9^{**}$) | 0.560 ($14.1^{***}$) / **0.134** ($3.8^{***}$) | **0.650** ($16.7^{***}$) / **0.310** ($8.5^{***}$) |
| — **Partial QAP ctrl Feb wd** $r$ ($z$) | — | **0.325** ($7.8^{***}$) / **0.085** ($1.92^{*}$) | **0.291** ($7.4^{***}$) / **0.142** ($4.0^{***}$) | **0.405** ($12.1^{***}$) / **0.251** ($7.0^{***}$) |
| **137 Cities × 47 Provs vs City Roots 2024:** | | | | |
| — Full bipartite spatial-null QAP $r$ ($z$) | 0.388 ($21.1^{***}$) | **0.486** ($19.2^{***}$) / **0.288** ($13.2^{***}$) | **0.457** ($19.4^{***}$) / **0.267** ($13.0^{***}$) | **0.489** ($23.6^{***}$) / **0.342** ($16.5^{***}$) |
| — **Neighbouring-city QAP** (within-prov perm) | 0.330 ($8.8^{***}$) | **0.559** ($11.5^{***}$) / **0.309** ($10.9^{***}$) | **0.354** ($9.5^{***}$) / **0.185** ($8.0^{***}$) | **0.393** ($10.4^{***}$) / **0.273** ($9.9^{***}$) |
| — **Neighbouring-city Partial QAP (ctrl Feb)** | — | **0.493** ($12.5^{***}$) / **0.345** ($11.2^{***}$) | **0.282** ($10.0^{***}$) / **0.238** ($9.6^{***}$) | **0.321** ($10.1^{***}$) / **0.306** ($10.5^{***}$) |

*Honest nuance from the multi-network horse race (47×47 Partial QAP on holiday excess, ctrl Feb wd):* At the coarse province level, **recent migration (2015–2024)** ($r = 0.223 / 0.288 / 0.349$ for Aug/Nov/Xmas excess, $z=6.8\text{–}12.0$) and **Facebook SCI 2026** ($r = 0.180 / 0.253 / 0.295$) outperform **1998 birthplace stock** ($r = 0.023 / 0.076 / 0.178$) and **1988–1995 migration** ($r = 0.127 / 0.204 / 0.303$). At the **137-city neighbouring-city scale**, however, birthplace roots predict holiday flows dramatically above ordinary weekdays ($r_{\text{nb,partial}} = 0.493$ in Mid-Aug, $z=12.5, p=0.0005$; $0.345$ on pure excess, $z=11.2$). Thus holiday mobility carries **both** the multi-decade rural-exodus echo (captured sharply by city-level birthplace composition and emptied-village $k$-core jumps) **and** ongoing contemporary migration/friendship ties.

#### B. Bipartite Communities ("Summer Territories") Under Exact BiCM & Spatial-Gravity Null + Contiguity Test

We cluster the bipartite network of **137 mainland Census cities $\times$ active mainland destinations** ($2,456$ in Mid-Aug) under three nulls: **(1) Exact Bipartite Configuration Model (BiCM, Saracco et al. 2015)** solved to degree error $<10^{-6}$ on binary links ($w \ge 5$), **(2) Barber weighted configuration model** ($s_i s_j / W$), and **(3) Spatial DC-gravity bipartite null** ($E_{ij} = a_i b_j f(d_{ij})$). We evaluate spatial contiguity on the exact polygon adjacency graph and benchmark Mid-August spatial-null communities against **100 Monte Carlo networks drawn from $\text{Poisson}(E_{ij})$** clustered identically.

| Slice & Null Model | $Q_B$ | Comms ($\ge 5$ dest) | Non-Contig Comms ($\ge 2$ pieces of $\ge 3$) | Diaspora Dest % | Diaspora Flow % | Roots Enrichment (in/out) | NMI (CCAA) |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Mid-Aug — Exact BiCM** (err $9.9\times 10^{-7}$) | 0.159 | 5 | 5 / 5 (tiny enclaves) | 8.8% | — | **2.80×** | 0.55 |
| **Mid-Aug — Barber Weighted CM** | 0.418 | 7 | 3 / 7 | 2.6% | — | **3.26×** | 0.71 |
| **Mid-Aug — Spatial DC-Gravity** | **0.165** | **16** | **12 / 16 (75%)** | **35.1%** | **41.3%** | **1.36×** | **0.31** |
| *— Mid-Aug Spatial-Null MC (100 draws)* | *0.0089 ± 0.0001* | *~12* | *11.8 ± 1.2* | *80.5 ± 2.0%* | *63.8 ± 4.1%* | *1.10 ± 0.036×* | — |
| **All Saints — Exact BiCM / Barber** | 0.336 / 0.514 | 5 / 7 | 4/5 / 2/7 | 3.2% / 0.7% | — | 2.99× / 3.64× | 0.65 / 0.79 |
| **All Saints — Spatial DC-Gravity** | **0.140** | **20** | **14 / 20 (70%)** | **29.9%** | **29.5%** | **1.40×** | 0.40 |
| **Feb Wed — Exact BiCM / Barber** | 0.268 / 0.390 | 5 / 6 | 4/5 / 2/6 | 6.0% / 1.9% | — | 2.75× / 3.17× | 0.61 / 0.74 |
| **Feb Wed — Spatial DC-Gravity** | 0.101 | 14 | 11 / 14 (79%) | 35.8% | 25.2% | 1.47× | 0.30 |

**Interpretation of the contiguity & Monte Carlo test (`reports/figs/net_B_qap_and_territories.png`):**
- Under non-spatial nulls (BiCM / Barber), Spain's summer map forms 5–7 large contiguous regional blocs ($\text{NMI}_{\text{CCAA}}=0.55\text{–}0.71$, $<9\%$ disconnected destinations).
- Once distance is factored out (Spatial DC-gravity null), **bipartite modularity remains overwhelmingly above sampling noise ($Q_B = 0.165$ vs MC $0.0089 \pm 0.0001$, $z = +1,168$)**, **12 of 16 territories split into major non-contiguous archipelagos carrying 41.3% of flow** (e.g. San Sebastián + Pamplona + Irún $\leftrightarrow$ Salamanca/Palencia interior villages [60% disconnected]; Bilbao/Barakaldo $\leftrightarrow$ Burgos/Zamora/Cantabria [38%]; Madrid metro $\leftrightarrow$ Torrevieja/Benidorm/Marbella/Gandia [72%]; Zaragoza/Huesca/Teruel $\leftrightarrow$ Salou/Peñíscola [25%]), and **cities grouped with destinations in the same spatial-null territory have $1.36\times$ higher birthplace-roots overlap vs $1.10 \pm 0.036\times$ under pure spatial-null noise ($z = +7.48, p < 0.01$)**.
- *Honest caveat:* Notice that `Feb Wed` spatial-null modularity ($0.101$) is lower than Mid-August ($0.165$) and has a much smaller diaspora flow share ($25.2\%$ vs $41.3\%$), yet still exhibits $1.47\times$ roots enrichment—confirming that diaspora corridors exist year-round and **amplify sharply in volume and modularity** during holidays. Notice also that pure Poisson spatial-null noise yields *even higher* random scatter ($80.5\%$ disconnected destinations), so observed spatial-null territories are **structured archipelagos** (much more coherent than random noise, much less contiguous than BiCM).

---

## Part 3 — One-Page "Network-First" Report Designs (COMP0123 Format)

### Design A (1 Page): *Do the 1833 Provincial Borders Still Partition Spain's Mobility Network?*

- **Core Research Question (1 sentence):** *Does the 1833 provincial division still leave an independent community boundary in Spain's daily mobility network beyond geography, road travel time, and pre-existing 1789 territories?*
- **Network Definition:** Undirected & directed weighted municipal mobility graphs of peninsular Spain ($N = 2,484$ nodes, persistent links $\ge 50\%$ of days: weekday $E = 122,374$, $\rho = 3.97\%$, 55.4M trips/day; compared against Sunday $E = 146,869, \rho = 4.76\%$ and August $E = 168,145, \rho = 5.45\%$). Nodes carry administrative (1833 province, 1979–83 CCAA, 1789 Floridablanca, 1834 partido judicial, diocese, language zone) and geometric attributes.
- **Null Models (Hierarchical Ladder):** (1) Erdős–Rényi; (2) Degree-preserving rewiring (30×) & 1st-order JDD-preserving rewiring (10×); (3) Newman–Girvan configuration model; (4) Doubly-constrained spatial configuration model (both Euclidean distance and OSM road travel-time deterrence); (5) Map placebos: 1,000 contiguous agglomerations (node-matched & area-matched) + 60 shifted/rotated real maps; (6) Degree-Corrected Spatial Stochastic Block Model (DC-SSBM = PPML).
- **Section-by-Section Metrics (aligned to COMP0123 grading):**
  1. *RQ & Design (2.5 pp):* Why Spain's 1833 reform (weak provincial tier vs strong 1978–83 autonomous regions + 1789/1822 historical natural experiments + language zones) is a unique laboratory for testing institutional memory in spatial networks.
  2. *Literature Review (2 pp):* Critical synthesis of 6–7 network-science papers below (bridging Expert 2011 spatial modularity, Newman 2016 modularity–SBM equivalence, and Zhou's rewiring/2nd-order mixing tools).
  3. *Methods (2.5 pp):* Persistence filtering, spatial configuration null via IPF, formal DC-SSBM $\leftrightarrow$ spatial modularity equivalence, contiguous & rotated map placebo generation, OSM travel-time matrix.
  4. *Results (4.5 pp):*
     - **§4.1 Global Topology & Zhou Diagnostics:** Table 1 + `net_A_topology.png` (CCDF, small-world $\sigma$, structural constraints on $r$, rich-club $\rho(k)$ vs spatial rich-club, Catalan $k$-core, and the novel finding that 2nd-order mixing $R_{\text{avg}}=+0.60$ is a spatial-embedding signature + BOSAM).
     - **§4.2 Partition Modularity Horse Race:** `net_A_modularity_race.png` (1833 vs CCAA vs 1789 vs Partidos vs Dioceses vs Voronoi under CM, Euclidean, and Road-TT nulls against 1,000 contiguous & 60 rotated placebos — including the urban speed-limit diagnostic).
     - **§4.3 DC-SSBM Causal Attribution:** `fig8_which_lines_matter.png` (new 1833 borders $\times 0.59$ = old 1789+1833 borders $\times 0.60$ vs erased 1789 borders $\times 0.88$; 1822 ghosts; language zones) + `fig7_halfhour_rule.png` (province rules $<30\text{ min}$, CCAA activates $>30\text{ min}$).
  5. *Discussion & Limitations (2.5 pp):* Mechanism (41% of road km managed by *Diputaciones*, health/school catchments), honest discussion of fine-grained Dioceses/Partidos and why travel-time modularity needs area-matched placebos, MITMA `_AM` nesting bias (conservative), 2025-method invariance.
- **Core Figures (4–5):** `net_A_topology.png` (4-panel topology + 2nd-order + BOSAM), `net_A_modularity_race.png` (placebo distributions & map race), `fig8_which_lines_matter.png` (1789 vs 1833 DC-SSBM), `fig7_halfhour_rule.png` (travel-time scale separation).

---

### Design B (1 Page): *Does the Rural Exodus Still Draw Spain's Holiday Network?*

- **Core Research Question (1 sentence):** *Does the 1950s–1970s rural exodus still govern the topological rewiring and bipartite community structure of Spain's holiday overnight mobility network?*
- **Three Networks on Shared Mainland Nodes ($N = 2,484$):**
  1. *Baseline Ordinary Network:* Feb weekday overnight stays (`Feb_wd_1d` $E=26,465, \rho=0.86\%$; `Feb_wd_mean` $\rho=1.05\%$ at $w \ge 10$).
  2. *Holiday Overnight Networks:* Mid-August ($\rho=3.74\%$), Aug 15 ($\rho=3.31\%$), All Saints Nov 1 ($\rho=1.58\%$, non-beach family control), Christmas Eve Dec 24 ($\rho=1.69\%$).
  3. *Historical Roots Networks:* Birthplace $\times$ residence (47×47 provinces; 137 mainland Census cities $\times$ 47 provinces) + 1900–2021 municipal population-loss attribute + 1988–2024 migration & Meta SCI comparison layers.
- **Null Models:** (1) ER & 25× degree-preserving rewiring (structural-constraint test on holiday disassortativity & rich-club collapse); (2) Doubly-constrained spatial gravity null for QAP residuals; (3) QAP (2,000 matrix & within-province neighbouring-city permutations) + Freedman–Lane partial QAP; (4) Exact BiCM (Saracco 2015), Barber weighted CM, and Spatial DC-gravity bipartite nulls + 100-draw Monte Carlo spatial-null contiguity test.
- **Section-by-Section Metrics (aligned to COMP0123 grading):**
  1. *RQ & Design (2.5 pp):* Universal network question (*how long does a historical migration shock echo in a nation's seasonal mobility network?*), why Spain's *España vaciada* + MITMA overnight stays + Census birthplace matrices provide a clean natural experiment, and why comparing ordinary vs non-beach (All Saints) vs beach (August) holidays isolates roots from tourism.
  2. *Literature Review (2 pp):* Critical synthesis of 6–7 network-science papers below (bipartite null models, seasonal mobility networks, QAP matrix comparisons, weighted rich-club & structural constraints).
  3. *Methods (2.5 pp):* Three-network construction ($w\ge 10$ sparse backbone $<5\%$ density), structural-constraint rewiring, spatial-null partial QAP (including within-province city permutation), exact BiCM fixed-point solver, polygon-contiguity diaspora test + MC spatial null.
  4. *Results (4.5 pp):*
     - **§4.1 Holiday Topological Phase Shift:** Table 2 + `net_B_topology_shift.png` (density 3.9×, disassortative shift $r:-0.13 \to -0.21$ proven structurally forced by degree sequence via Zhou–Mondragón rewiring, rich-club collapse $\rho(100):1.25 \to 0.99$, population-sorted BOSAM "curtain").
     - **§4.2 Node-Level $k$-Core Promotion of Emptied Villages:** Size-matched inland villages ($\text{pop}<6\text{k}$) with $-65\%$ historical loss jump $k$-core $9.2 \to 41.5$ and $k_{\text{in}}$ $5.6 \to 47.6$ ($8.5\times$) vs $3.0\times$ for equal-sized growing towns, holding on All Saints (no beaches).
     - **§4.3 Edge-Level Spatial-Null QAP & Neighbouring-City Test:** `net_B_qap_and_territories.png` panels a–b (137-city neighbouring-city partial QAP $r=0.493, z=12.5$; honest comparison vs recent migration & Meta SCI).
     - **§4.4 Bipartite "Summer Territories" & Diaspora Test:** `net_B_qap_and_territories.png` panels c–d (BiCM contiguous regional blocs vs Spatial-gravity non-contiguous archipelagos: 12/16 non-contiguous, 41% diaspora flow, $1.36\times$ roots enrichment at MC $z=+7.5$).
  5. *Discussion & Limitations (2.5 pp):* Distinguishing historical exodus echo from contemporary migration/SCI (transparently acknowledging recent migration dominates coarse 47×47 flows while city-level roots & village $k$-core capture the 1950s–70s exodus), income/beach controls (from parallel statistical check), MITMA home-reassignment attenuation bias, policy relevance for rural infrastructure sizing.
- **Core Figures (4):** `net_B_topology_shift.png` (4-panel topology, mixing, $k$-core promotion, BOSAM), `net_B_qap_and_territories.png` (QAP comparison + BiCM vs Spatial-null territory maps), `fig9_exodus_holidays.png` (village refill + Barcelona–Andalusia diaspora), plus an optional rich-club collapse & $k$-core spatial map.

---

## Part 4 — Verified Network-Science Literature to Review (5–8 Papers Each)

### For Question A (1833 Borders & Spatial Network Community Structure)
1. **Expert, P., Evans, T. S., Blondel, V. D., & Lambiotte, R. (2011).** *Uncovering space-independent communities in spatial networks.* PNAS, 108(19), 7663–7668. **DOI:** `10.1073/pnas.1018962108` — *Core precedent:* introduces gravity-based spatial modularity and discovers Belgium's linguistic border; we improve their population-based null to a doubly-constrained travel-time configuration null and find the *opposite* result in Spain (institutions beat language).
2. **Sarzynska, M., Leicht, E. A., Chowell, G., & Porter, M. A. (2016).** *Null models for community detection in spatially embedded, temporal networks.* Journal of Complex Networks, 4(3), 363–406. **DOI:** `10.1093/comnet/cnv027` — *Methodological foundation:* systematically compares Newman–Girvan, gravity, and radiation nulls in spatial community detection; justifies our null-model ladder.
3. **Thiemann, C., Theis, F., Grady, D., Brune, R., & Brockmann, D. (2010).** *The Structure of Borders in a Small World.* PLOS ONE, 5(11), e15422. **DOI:** `10.1371/journal.pone.0015422` — *Direct conceptual precedent:* extracts effective mobility borders in the US and contrasts them with administrative states without controlling for travel time or historical border counterfactuals.
4. **Ratti, C., Sobolevsky, S., Calabrese, F., Andris, C., Reades, J., Martino, M., Claxton, R., & Strogatz, S. H. (2010).** *Redrawing the Map of Great Britain from a Network of Human Interactions.* PLOS ONE, 5(12), e14248. **DOI:** `10.1371/journal.pone.014248` *(and Sobolevsky et al., 2013, PLOS ONE, DOI: `10.1371/journal.pone.0081707`)* — *Key critique target:* shows standard modularity recovers administrative regions, which we prove is largely driven by distance decay unless tested under a spatial configuration null.
5. **Karrer, B., & Newman, M. E. J. (2011).** *Stochastic blockmodels and community structure in networks.* Phys. Rev. E, 83, 016107. **DOI:** `10.1103/PhysRevE.83.016107` — *together with* **Newman, M. E. J. (2016).** *Equivalence between modularity optimization and maximum likelihood methods for community detection.* Phys. Rev. E, 94, 052315. **DOI:** `10.1103/PhysRevE.94.052315` *(and Bernaschi et al., 2022, Sci. Rep., DOI: `10.1038/s41598-022-22798-6`)* — *Theoretical bridge:* underpins our proof that PPML gravity border estimation is a degree-corrected spatial SBM whose score test is spatial modularity.
6. **Zhou, S., & Mondragón, R. J. (2007).** *Structural constraints in complex networks.* New Journal of Physics, 9(6), 173. **DOI:** `10.1088/1367-2630/9/6/173` *(plus Zhou & Mondragón 2004, IEEE Comm. Lett., DOI: `10.1109/LCOMM.2004.823426`)* — *Prof. Zhou's core method:* degree-preserving rewiring to separate forced degree-sequence effects from genuine assortative/rich-club organization.
7. **Zhou, S., Cox, I. J., & Hansen, L. K. (2017).** *Second-order assortative mixing in social networks.* Complex Networks VIII (CompleNet 2017), Springer, pp. 3–13. **DOI:** `10.1007/978-3-319-54241-6_1` *(plus Guo, Chen & Zhou 2009 BOSAM, DOI: `10.1007/978-3-642-02469-6_45`)* — *Direct novel extension:* we deliver the first application of Zhou's second-order mixing to a national mobility network and show it is driven by spatial embedding.
8. **Colizza, V., Flammini, A., Serrano, M. A., & Vespignani, A. (2006).** *Detecting rich-club ordering in complex networks.* Nature Physics, 2, 110–115. **DOI:** `10.1038/nphys209` *(and Opsahl et al., 2008, PRL, DOI: `10.1103/PhysRevLett.101.168702`)* — *Rich-club normalisation standard:* motivates comparing binary, weighted, and spatial-null rich-club coefficients.

### For Question B (Rural Exodus, Holiday Mobility & Bipartite Territories)
1. **Saracco, F., Di Clemente, R., Gabrielli, A., & Squartini, T. (2015).** *Randomizing bipartite networks: the case of the World Trade Web.* Scientific Reports, 5, 10595. **DOI:** `10.1038/srep10595` — *Core bipartite null:* exact maximum-entropy Bipartite Configuration Model (BiCM), which we solve to $<10^{-6}$ error and contrast with our spatial-gravity bipartite null.
2. **Barber, M. J. (2007).** *Modularity and community detection in bipartite networks.* Phys. Rev. E, 76, 066102. **DOI:** `10.1103/PhysRevE.76.066102` — *Bipartite modularity framework:* defines $Q_B$ for city-by-destination "summer territories" before and after removing distance.
3. **Louail, T., Lenormand, M., Picornell, M., García Cantú, O., Herranz, R., Frias-Martinez, E., Ramasco, J. J., & Barthelemy, M. (2015).** *Uncovering the spatial structure of mobility networks.* Nature Communications, 6, 6007. **DOI:** `10.1038/ncomms7007` — *Spanish mobile-phone network benchmark:* classifies 31 Spanish cities via hotspot OD structure; we scale from intra-urban commuting to national inter-municipal holiday rewiring.
4. **Wesolowski, A., zu Erbach-Schoenberg, E., Tatem, A. J., Lourenço, C., Viboud, C., Charu, V., Eagle, N., Engø-Monsen, K., Qureshi, T., Buckee, C. O., & Metcalf, C. J. E. (2017).** *Multinational patterns of seasonal asymmetry in human movement influence infectious disease dynamics.* Nature Communications, 8, 2069. **DOI:** `10.1038/s41467-017-02070-z` — *Seasonal holiday mobility precedent:* documents large-scale holiday asymmetry across countries, without linking destinations to multi-decade historical rural exodus or birthplace networks.
5. **Li, J., Ye, Q., Deng, X., Liu, Y., & Liu, Y. (2016).** *Spatial-Temporal Analysis on Spring Festival Travel Rush in China Based on Multisource Big Data.* Sustainability, 8(11), 1184. **DOI:** `10.3390/su8111184` — *Closest holiday-migration QAP precedent:* relates China's *Chunyun* flows to contemporary migrant workers via QAP; Spain tests a 50–70 year intergenerational echo (children/grandchildren returning to emptied villages) under a spatial DC-gravity null.
6. **Dekker, D., Krackhardt, D., & Snijders, T. A. B. (2007).** *Sensitivity of MRQAP tests to collinearity and autocorrelation conditions.* Psychometrika, 72(4), 563–581. **DOI:** `10.1007/s11336-007-9016-1` — *Methodological standard for network comparison:* justifies our Freedman–Lane double-semi-partialing QAP controlling for the everyday mobility network.
7. **Zhou, S., & Mondragón, R. J. (2007).** *Structural constraints in complex networks.* New Journal of Physics, 9(6), 173. **DOI:** `10.1088/1367-2630/9/6/173` *(plus Guo, Chen & Zhou 2009 BOSAM, DOI: `10.1007/978-3-642-02469-6_45`)* — *Explains our headline topological result:* holiday disassortativity ($r=-0.21$) is structurally forced by the degree sequence, and BOSAM makes the city-to-village transition immediately visible.
8. **Barbosa, H., Barthelemy, M., Ghoshal, G., James, C. R., Lenormand, M., Louail, T., Menezes, R., Ramasco, J. J., Simini, F., & Tomasini, M. (2018).** *Human mobility: Models and applications.* Physics Reports, 734, 1–74. **DOI:** `10.1016/j.physrep.2018.01.001` *(and Bailey et al., 2018, JEP, DOI: `10.1257/jep.32.3.259` for Facebook SCI)* — *Authoritative spatial-mobility & social-tie review* situating gravity/radiation nulls and digital friendship layers.

---

## Part 5 — Figures & Code Inventory

| Asset | Path | Content |
|---|---|---|
| Script A1 | `explore/p90_A_basic_properties.py` | Full basic properties, 30× rewiring, rich-club (binary, Opsahl, spatial), 2nd-order mixing vs spatial null, BOSAM |
| Script A2 | `explore/p91_A_modularity_placebos.py` | Multi-map modularity race (6 maps × 3 nulls × 1,000 node-contig + 1,000 area-contig + 60 rotated placebos) |
| Script A3 | `explore/p94_check_2nd_order_and_tails.py` | 10× JDD-preserving rewiring (exact Zhou 2017 test) + Vuong power-law vs log-normal tail tests |
| Script B1 | `explore/p92_B_networks_topology.py` | 6 overnight + 2 roots networks, 25× rewiring, disassortative shift, rich-club collapse, emptied-village $k$-core promotion, BOSAM |
| Script B2 | `explore/p93_B_qap_and_bipartite.py` | 47×47 & 137×47 spatial-null & partial QAP (2,000 perms), exact BiCM, Barber & spatial bipartite territories, 100× MC contiguity test |
| Fig A1 | `reports/figs/net_A_topology.png` | 4-panel: Degree CCDF, norm. rich-club (30 rewirings), 2nd-order mixing vs spatial/rewired nulls, BOSAM weekday vs Aug |
| Fig A2 | `reports/figs/net_A_modularity_race.png` | 2-panel: 1,000 contiguous placebo distributions + 6-map modularity horse race (distance & road-TT nulls) |
| Fig B1 | `reports/figs/net_B_topology_shift.png` | 4-panel: In-degree CCDF, assortativity shift, emptied-village $k$-core promotion, population-sorted BOSAM |
| Fig B2 | `reports/figs/net_B_qap_and_territories.png` | 4-panel: Spatial-null & neighbouring-city QAP, historical-layer horse race, BiCM vs Spatial-gravity summer territories maps |

---

## Part 6 — Honest Verdict: How Network-Native Are A and B Now, and Which Is Stronger for COMP0123?

| Evaluation Axis (Prof. Shi Zhou's COMP0123) | Candidate A (1833 Borders) | Candidate B (Rural Exodus & Holidays) |
|---|---|---|
| **Network-Nativeness Score (0–10)** | **7.0 / 10** | **8.5 / 10** |
| **Why this score?** | The network layer is now genuinely rich (30× rewiring, rich-club, 2nd-order mixing = spatial embedding, BOSAM, DC-SSBM $\leftrightarrow$ modularity equivalence, 1,000 contiguous placebos). **However**, the *killer causal proof* (new vs old vs erased 1789 borders, 1822 ghosts, language zones, 30-min rule) still lives in edge-level PPML/DC-SSBM regressions, while pure network modularity under the road-TT null has a subtle vulnerability: Dioceses ($z_{\text{area}}=+7.47, z_{\text{rot}}=+6.38$) and Partidos outrank 1833 Provinces, and uniform-node contiguous blobs drop $z_{\text{tt}}$ to $+0.53$ unless area-matched ($+3.41$). | **Every single headline finding is now a genuine Complex Networks phenomenon** that maps directly onto COMP0123 lectures: comparing **three distinct networks** on the same 2,484 nodes; a dramatic **topological phase shift** (density $3.9\times$, disassortative shift structurally forced by degree sequence via Zhou–Mondragón rewiring, rich-club collapse $1.25 \to 0.99$, population-sorted BOSAM curtain); **size-matched $k$-core promotion** of emptied villages ($9.2 \to 41.5$); **spatial-null & neighbouring-city partial QAP** ($z=12.5$); and **exact BiCM vs spatial-gravity bipartite communities** ($z_Q=1,168$, roots enrichment $z=+7.5$). |
| **Biggest Risk for COMP0123** | A strict network physicist could still ask: *"Why did your modularity horse race put Dioceses ahead of 1833 Provinces, and isn't the 1789/1833 core really an econometric gravity paper dressed up as a network?"* | Two honest caveats that must be handled cleanly: (1) At the coarse 47×47 level, recent migration (2015–24) and Facebook SCI beat 1998 birthplace stock—so the report must anchor the historical exodus claim on **node-level village $k$-core promotion** + **137-city neighbouring-city QAP** (where income control from the parallel agent is decisive!). (2) Bipartite spatial-null communities are somewhat noisy/mixed (beach + roots). |
| **Simplicity & "Wow" for a 15-page Report** | High intellectual depth, slightly heavier historical scaffolding (1789, 1822, 1833, CCAA, partidos, dioceses). | Effortless one-sentence question, emotionally vivid, cleaner narrative arc, visually unforgettable (BOSAM curtain + summer territories maps). |

### Bottom-Line Recommendation
- **For COMP0123 specifically, Candidate B is now the stronger Complex Networks report (8.5/10 vs 7.0/10)** — **provided** the parallel statistical audit confirms that the 137-city neighbouring-city roots effect survives controlling for municipal income! Why? Because in B, **the network topology itself transforms before your eyes** (assortativity, structural constraints, rich-club dissolution, $k$-core promotion, BOSAM, QAP, exact BiCM vs spatial bipartite territories). Prof. Shi Zhou explicitly rewards *comparing real networks*, *structural-constraint rewiring*, *rich-club*, *$k$-core*, and *BOSAM*—and B uses **all five** to answer its core question rather than bolting them on as Section 4.1 preliminaries.
- **Keep Candidate A (7.0/10) as an ultra-solid fallback** if the income control seriously weakens B's city-level roots signal: A's causal identification (1789 vs 1833 new/old/erased borders) is rock-solid and already bulletproof, and with the DC-SSBM framing + 2nd-order mixing discovery + area-matched modularity placebos it is still a comfortable Distinction-level COMP0123 report.

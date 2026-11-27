# Table A1. Basic properties of the overnight-stay networks, full version

| | Ordinary weekday<br>13–14 Feb | Summer holiday<br>14–15 Aug | All Saints<br>31 Oct–1 Nov | Summer, same directed link count<br>14–15 Aug |
|---|---:|---:|---:|---:|
| **Size** | | | | |
| People sleeping away from home (millions) | 2.11 | 8.46 | 4.23 | 8.46 |
| Share of them on kept links (%) | 69 | 83 | 76 | 67 |
| Nodes N (in the giant component) | 2,435 (2,432) | 2,435 (2,435) | 2,435 (2,434) | 2,435 (2,422) |
| Links L (undirected) | 26,465 | 101,996 | 48,698 | 31,573 |
| Density (%) | 0.89 | 3.44 | 1.64 | 1.07 |
| **Degree** | | | | |
| Mean degree ⟨k⟩ | 21.7 | 83.8 | 40.0 | 25.9 |
| Median degree | 11 | 45 | 22 | 10 |
| Maximum degree k_max (% of N) | 1,499 (62 %) | 2,233 (92 %) | 2,008 (82 %) | 1,829 (75 %) |
| Median out-degree / in-degree | 7 / 7 | 14 / 27 | 10 / 16 | 2 / 6 |
| Tail fit, in / out: power law vs log-normal | tie / tie | log-normal / log-normal | tie / tie | log-normal / log-normal |
| **Small world** | | | | |
| Average path length ℓ (rewired / ER) | 2.44 (2.41 / 2.83) | 2.00 (1.98 / 2.02) | 2.16 (2.11 / 2.49) | 2.25 (2.24 / 2.74) |
| ℓ* = ln N / ln⟨k⟩ | 2.53 | 1.76 | 2.11 | 2.39 |
| Diameter (rewired) | 4 (5.1) | 3† (3.2) | 4† (4.0) | 4† (4.9) |
| Clustering ⟨C⟩ (rewired / ER) | 0.554 (0.247 / 0.009) | 0.546† (0.387 / 0.034) | 0.593 (0.298 / 0.016) | 0.678 (0.436 / 0.011) |
| **Mixing** | | | | |
| Degree assortativity r (rewired) | −0.131 (−0.140) | −0.214 (−0.217) | −0.155 (−0.161) | −0.187 (−0.186) |
| Range of r reached by rewiring [r_min, r_max] | [−0.221, −0.092] | [−0.424, −0.114] | [−0.274, −0.098] | [−0.271, −0.148] |
| Reciprocity (rewired, directed) | 0.52 (0.13) | 0.24 (0.09) | 0.48 (0.15) | 0.24 (0.07) |
| 2nd-order mixing R_avg (null keeping P(k) and φ) | 0.000 (−0.023) | −0.289 (−0.252) | −0.129 (−0.127) | −0.243 (−0.191) |
| **Rich club and core** | | | | |
| Rich club φ, top 50 zones (rewired) | 0.80 (0.66) | 0.96 (0.92) | 0.91 (0.81) | 0.83† (0.81) |
| Rich-club clique size (rewired) | 17 (5.8) | 17† (9.2) | 22 (7.1) | 9† (7.9) |
| Links touching the top 5 % of nodes (%) | 65 | 58 | 63 | 74 |
| Largest k of a non-empty k-core (rewired / ER) | 37 (30 / 15) | 84 (83 / 67) | 51 (45 / 30) | 48 (43 / 18) |

Each column is the network of one night of 2024: node = one of 2,435 MITMA zones of mainland Spain; a directed link i→j joins zones when at least 10 residents of i sleep in j (last column: the summer night cut to the same number of directed links as the ordinary night, 35,812, by keeping the heaviest). Properties are computed on the undirected projection unless marked directed. In brackets: mean of 50 networks with the same degree sequence (degree-preserving rewiring, 10 swaps per link) and, where shown, of 50 Erdős–Rényi networks with the same N and L. [r_min, r_max]: the most disassortative and most assortative values reached by a heuristic rewiring search that keeps P(k) and a single connected component (Zhou & Mondragón 2007). R_avg: second-order mixing (Zhou, Cox & Hansen 2017), null keeping P(k) and φ(k). Rich-club clique: the largest n such that the n best-connected zones are all linked to one another (Zhou & Mondragón 2007). Tail fit: Clauset–Shalizi–Newman likelihood-ratio test; 'tie' means p ≥ 0.05. † The change from the ordinary night does not keep its direction in every pair of sister nights (117 summer pairs, 39 All Saints pairs). 1 Nov is the quietest of its three sister nights. 

# V1: independent verification of the N1 findings (exodus cohort, village-born, casa del pueblo)

*The code is `pipeline/v1_verify_n1.py`, written by the Antigravity verifier (Gemini 4 Carbon). It stopped on a quota error after Part 1. The orchestrator re-ran the full script locally on 2026-11-24 (log `results/v1_run_orchestrator.log`). Parts C1–C4 finished; C5 (the municipal INE OD network) crashed on a sheet-name lookup (`'202402' is not in list`). C5 is for the journal paper only and is left pending. This summary was written by the orchestrator from the CSV outputs.*

## Verdicts

| Claim | Reproduced | Verdict | Notes |
|---|---|---|---|
| C1a. Pre-1981 arrivals: δ(Aug 15 − Feb 14) = +0.084 | +0.084, p_WCR < 0.001, CI [0.046, 0.117] | **CONFIRMED** | Joint model with recent arrivals; 675 municipalities, G = 44 |
| C1b. Recent arrivals (2016–23): δ = −0.072 | −0.072, p = 0.007 | **CONFIRMED** | |
| C1c. Village-born (≤ 2,000): δ = +0.086 (legacy text: +0.087) | +0.086, p < 0.001 | **CONFIRMED** | |
| C1d. Urban-born (> 10k): δ = −0.062 | −0.062, p = 0.03 | CONFIRMED, fragile | Significant in only 85/117 sister pairs; not used as a result |
| C3. "The second-home stock mediates 65 %" | 64.6 % [43.8, 98.4] (INE 53001, Aug/Feb) | CONFIRMED, with caveat | A cross-sectional decomposition (Gelbach-type), **not causal mediation**. Say "accounts for" |
| C3. "…but not at Christmas" | MITMA single nights: Aug 15 48 % [22, 82], Good Friday 26 %, All Saints 13 %, Christmas Eve **1.8 %** [−1.7, 7.5] | CONFIRMED **within the MITMA data** | INE monthly December gives 43 % [13, 239]. That is a whole month including the December long weekends, not Christmas Eve, so never mix the two sources in one contrast |
| C4. ECEPOV: 58.4 % of cross-CCAA second homes lie in a member's or a parent's birth CCAA; 77.4 % of inland ones | Exact; 2.9× and 2.0–2.6× the permutation null | **CONFIRMED** | A lower bound (13.3 % of weighted cross-CCAA second-home households miss a parent birthplace). Madrid is 48 % of the weight; without Madrid, inland is 79.0 % |
| C5. Municipal INE OD node gains | Not re-run (crash) | PENDING | The verifier's notes flag that suppression of cells < 30 inflates the magnitude |

## Robustness of C1 (new; the N1 agent had not run any of this)
- **Sister nights (9 × 13 = 117 pairs):**
  - pre-1981: δ median +0.079 [0.060, 0.101]; same sign 117/117; significant 117/117.
  - village-born: +0.080 [0.055, 0.105]; 117/117 and 117/117.
  - recent: −0.075; 117/117 and 117/117.
- **Inland vs coastal destinations:**
  - pre-1981 inland +0.154 (117/117 significant); coastal −0.017 (n.s.).
  - village-born inland +0.123; coastal +0.019 (n.s.).
  - **The holiday surge goes entirely to inland destinations.**
- **Other holidays vs February:**
  - Good Friday: pre-1981 +0.057 (13/13 significant), village-born +0.068 (13/13).
  - All Saints: pre-1981 +0.031 (19/39 significant), village-born +0.039 (38/39).
  - Christmas Eve: pre-1981 +0.026 (10/13), village-born +0.030 (12/13).
  - Recent arrivals do not fall at Christmas (+0.001).
- **October baseline:** pre-1981 +0.086 (15/15 significant), village-born +0.095 (15/15).
- **Zero cells** (46–60 %; Haldane +0.5): the sign is stable across 6 treatments and in R2's 112 large cities (G = 22, few zeros: pre-1981 +0.082, p = 0.025; village-born +0.098, p < 0.001). The magnitude depends on the treatment, so do not quote ratios such as "triples".
- **Leave-one-province-out (44 provinces):** pre-1981 δ ranges 0.072–0.089 (max p 0.001); village-born 0.081–0.092 (max p < 0.001). **No single province carries the result.**
- **Age confound (the key check), Aug 15 − Feb 14, all destinations:**
  - Horse race with "born in p and aged 65+": pre-1981 stays at +0.060 (p 0.001); the 65+ term is +0.021 (p 0.21, n.s.).
  - Four-way split:
    - pre-1981 and 65+: **+0.087** (p < 0.001);
    - **post-1981 and 65+: −0.006 (p 0.73);**
    - pre-1981 and under 65: +0.013 (n.s.; inland +0.037, p 0.006);
    - recent and under 65: −0.072.
  - So older people who migrated later do **not** do this. The effect belongs to the generation that left before 1981, not to age as such.
  - Village-born under 65: +0.037 (p 0.034; inland +0.051, p 0.021). Younger village-born people also go back.
  - City age structure × destination slopes: pre-1981 +0.080 (p 0.001).

## Most defensible sentence
> Compared with an ordinary February night, on 15 August cities send disproportionately more people to the provinces their **pre-1981 migrants** came from, and to those of their **village-born** residents. The pattern holds on all 117 summer × ordinary night pairs and in every leave-one-province-out run. It goes entirely to inland destinations and is absent towards the coast. It is not an age effect: residents aged 65+ who migrated after 1981 show no such pull. Recent migrants show the opposite pattern.

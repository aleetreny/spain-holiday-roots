# Independent audit run log

All Python entry points were run with .venv/bin/python -I from the project root. No packages were installed. The small graph swap sampler was compiled locally. Source and report directories were read as inputs.

| Step | Commands and retained evidence |
|---|---|
| District and node rebuild | rebuild.py and rebuild.log. Regenerate the temporary Parquet panels with this script if necessary. |
| Graph measures and ensembles | network.py produced network_real.csv, network_null.csv and tails.csv. constraints.py and swaps.cpp generated constraints.csv. er_check.py extended the August ER path comparison to 116 independent draws. |
| Figures and equal-edge sister graphs | figures.py, diagnostics.py and equal_sisters.py, with corresponding logs. sister_network.csv and sister_eq_network.csv contain independent cross-date real statistics. |
| R1 | sisters_nodes.py produced r1_level.csv, r1_ancova.csv and r1_sisters.csv. r1_inference.py produced R1_wcr.csv. gravity.py produced gravity_fit.csv and gravity_slopes.csv. |
| R2 | roots.py R2 pop full, r2_extra.py, qap.py and inference.py R2. R2_pop.csv, r2_extra.csv, roots_index.csv, qap_reconciled.log and R2_wcr.csv. |
| V1 | roots.py V1 pop full. inference.py V1 sisters cohort and inference.py V1 sisters birth, then near_cutoff.py to examine two borderline contrasts. jackknife_inference.py with cohort and birth. figure3_inference.py independently inverted all six Figure 3 mainland/coastal intervals. zero_sensitivity.py checked six count and share transformations. V1_wcr.csv includes the complete all-destination main comparison and some inland outputs. Its expanded secondary inference grid was stopped after the required main and sister checks finished. |
| Housing and ECEPOV | secondary.py repeated descriptive, permutation-expectation, OLS and cluster pairs checks. The 19,999-draw precision check is secondary_19k.log. |
| Audit cells | table_audit.py produced table_cells.md, included in the signed-off report. 39 main table and 88 appendix entries. |

The first district aggregation retained a few zero-weight pairs. All were subsequently verified as exactly zero and explicitly excluded when comparing positive edge lists. A first secondary-data read selected subtotals in addition to the national rows. The audit now reads only national totals with the next two hierarchy fields empty. Initial rooted QAP checks used fixed rather than quantile distance bins, and a standardized difference. The audit reconciled those choices and independently recovered the reported statistics using quantile bins and original log-residual differences.

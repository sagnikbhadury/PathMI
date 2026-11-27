# Results deposited with the manuscript

These are the outputs of the study *Measuring molecular recoverability in pathology foundation models*. The study used the TCGA colorectal cohort (COAD + READ, state MSI) and the TCGA prostate cohort (PRAD, states ERG and ETS). The models were UNI v2, Virchow2, H-optimus-1 and CONCH v1.5, with embeddings extracted by Trident at 20X.

The unit of analysis is the **patient**. A patient's slides are pooled separately, and the slide vectors are then averaged, giving 693 patients (727 slides). Every estimate is in nats.

The folder holds summary tables only. It contains no tile embeddings, pooled embedding matrices, projected coordinates or model objects that could reconstruct them. The scripts that produced every file are in [`../study`](../study).

| Folder | Content |
|---|---|
| `grid/` | `mi_summary_table.xlsx`: the 192-cell grid (state × model × pooling × projection), with the estimate, the permutation-null mean and s.d., the p value and the condition number. `best_kpca_summary.xlsx`: the 48 projections selected by max-T over a joint permutation. `incremental_mi_summary.xlsx`: covariate-only, embedding-only and combined estimates and the incremental recoverability beyond measured covariates. `kpca_selected_decomposition.xlsx`: the selected decomposition per arm. |
| `per_run/<cohort>_<model>/` | The per-run tables behind the grid. `kpca_mi_summary_*.csv` has the 16 projections with nulls, raw estimates, k* and cond(H). Alongside it are `best_kpca_mi_summary_*.csv`, `kpca_selected_decomposition_*.csv` and `incremental_mi_*.csv`. |
| `record/` | The Recoverability Record: `recoverability_record.csv` and `recoverability_record.json` hold the 12 measurement arms plus the permuted-label negative control, and `json/` holds one file per Record. Each Record gives provenance, the estimate, the fixed-arm null, p_max16, the numerical tolerance ε_tol (column `eps_rep`), half-sampling s.e. and 90 % interval, incremental and usable counterparts, tile burden, status and N_evidence; the JSON (schema 1.1) adds protocol notes and the fields `supported_interpretations` and `interpretations_not_established`. `status_recheck.*` re-derives every status, margin and N_evidence under the tolerance and checks both resampling formulas against the raw draws. |
| `usable/` | Usable recoverability for all 192 cells and the control. It comes from patient-held-out cross-fitting, with standardisation, projection and readout all fitted inside the training patients. The readouts are ridge and the PathMI readout; the fixed-projection and in-sample comparators are included. Per-repeat values and summary statistics are alongside. |
| `learning/` | Patient-number learning curves for the 12 Record arms and the control. `learning_summary.csv` gives the estimate, null, usable value, finite-population s.e., lower bound, threshold and resolved flag at each n. `learning_draws.csv` has every subsample. `n_evidence.csv` gives N_evidence and the full-cohort status. |
| `robustness/` | Half-sampling (`subsampling_*`, delete-d jackknife s.e. and 90 % interval). The bootstrap with replacement (`bootstrap_*`) is kept to show its upward bias. Also here: the retained-dimension sweep, the v_beta sweep, the chain-seed and size-prior refits, contrasts in tolerance units (`resolution_units.csv`), and the numerical tolerance (`tolerance.csv` by model group, `tolerance_identities.csv` by identity and model, `tolerance_specification.csv` for the null covariate block, `tolerance.txt`). |
| `max16/` | The joint-permutation test over all 16 configurations of each model and state (`max16.csv`: winner, observed maximum, null maxima, p_max16), and the alignment checks (`alignment.csv`, `refit.csv`: 48 stored permuted fits refitted exactly). |
| `burden/` | Tissue-sampling burden: N90/N95 per arm with bootstrap ranges, burden curves, the ordering across models, Kendall's W, and the H-optimus-1 linear-arm curve. |
| `tile_maps/<cohort>/` | The tile-level decomposition of the H-optimus-1 mean-pooling / linear-PCA arm. `fit_summary_*.csv` gives the selected model, estimate and plug-in. `slide_summary_*.csv` has one row per slide (patient, label, fitted logit and probability, tile count, share of high-surprise tiles under the thresholding rule, column `frac_informative`). The log records the exact-identity checks. |
| `simulation_context/` | The hidden-context simulation: Monte Carlo truths, replicate estimates, nulls and the summary at π_A = 1, 0.5, 0.25. |
| `simulation_estimators/` | The 72-cell × 50-replicate comparison of PathMI, LASSO-CV, OLS and the linear Donsker–Varadhan bound on simulated data. |
| `prad_batch_purity/` | Prostate refits with scanning batch and ABSOLUTE purity added to the covariates. |

Model labels are `uni_v2`, `virchow2`, `hoptimus1` and `conch_v15`.

The TCGA slides are public through the GDC. The embeddings are not distributed here; they can be regenerated with Trident ([`mahmoodlab/TRIDENT`](https://github.com/mahmoodlab/TRIDENT)) under each model's own licence.

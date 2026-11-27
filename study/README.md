# Study scripts

These are the scripts that produced the analyses and figures of the manuscript, kept as they were run.

They read Trident feature files (`.h5`) and the pooled matrices built from them. **Neither is distributed in this repository.** To rerun the analyses, extract the features with Trident for the TCGA COAD, READ and PRAD slides, then point the paths at the top of each script to them. Every analysis is also available as a documented function of the package, and those functions run on any embeddings (see the main README).

| Script | Analysis | Package function |
|---|---|---|
| `estimator/PathMI_SUBHO_logistic.R` | the estimator (v3) | `pathmi_fit()`, `subho_MI_logistic()` |
| `estimator/simulation_MI_comparison.R`, `estimator/PathMI_DOE_linear.R` | estimator comparison | `simulate_estimator_design()`, `doe_linear()` |
| `pool_pass.R`, `projection_functions.R`, `infold_functions.R` | pooling, PCA and kPCA with out-of-sample extension | `pool_cohort()`, `fit_projection()` |
| `patient_level/pl_grid.R`, `pl_promote.R`, `run_patient_chain.sh` | the patient-level 192-cell grid with shared permutations and the full 30-step chain | `pathmi_null()`, `select_representation()` |
| `patient_level/max16.R` | max-T over all 16 pooling x projection configurations of each model and state (p_max16), with the alignment check (complete permutations, seed reconstruction, exact refits of stored permuted fits) | `select_representation()`, `verify_permutation_alignment()` |
| `patient_level/pl_incremental.R`, `prad_incremental_batch_purity.R` | incremental recoverability beyond measured covariates | `incremental_recoverability()` |
| `usable_infold.R`, `usable_infold_stats.R` | patient-held-out usable recoverability | `usable_recoverability()` |
| `robustness_subsampling.R`, `robustness_prior_bootstrap.R`, `robustness_dimension.R` | half-sampling, bootstrap bias, dimension and prior sweeps | `half_sampling()` |
| `tolerance.R`, `resolution_units.R`, `resolution_table.R`, `patient_level/check_refit.R` | numerical tolerance from replicate and sign-flip refits (the null covariate block reported separately as a specification check), and every contrast in tolerance units | `numerical_tolerance()`, `replicate_refit()`, `sign_flip_refit()`, `null_block_refit()` |
| `learning_curves.R`, `learning_summary.R` | learning curves and N_evidence | `learning_curve()`, `evidence_status()` |
| `status_recheck.R` | both resampling formulas and the evidence statuses re-derived from the raw draws under the tolerance | `learning_curve()`, `evidence_status()` |
| `patient_level/run_final_edits_chain.sh` | the final rerun: tolerance, learning summary, Record, status check, tolerance units and all figures | |
| `tile_maps_hoptimus1.R`, `sampling_burden_all_models.R`, `burden_*.R` | tile decomposition, specific surprise, tissue burden | `tile_decomposition()`, `tissue_burden()` |
| `recoverability_record.R` | Recoverability Record (CSV, JSON, table) | `recoverability_record()`, `record_to_json()` |
| `hidden_context_sim.R` | hidden-context simulation | `simulate_hidden_context()` |
| `figures_*.R`, `tables_main.R`, `*_table.R` | manuscript figures and tables | `plot_*()` |

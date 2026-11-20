# Study scripts

These are the scripts that produced the analyses and figures of the manuscript, kept as they were run.

They read Trident feature files (`.h5`) and the pooled matrices built from them. **Neither is distributed in this repository.** To rerun the analyses, extract the features with Trident for the TCGA COAD, READ and PRAD slides, then point the paths at the top of each script to them. Every analysis is also available as a documented function of the package, and those functions run on any embeddings (see the main README).

| Script | Analysis | Package function |
|---|---|---|
| `estimator/PathMI_SUBHO_logistic.R` | the estimator (v3) | `pathmi_fit()`, `subho_MI_logistic()` |
| `estimator/simulation_MI_comparison.R`, `estimator/PathMI_DOE_linear.R` | estimator comparison | `simulate_estimator_design()`, `doe_linear()` |
| `pool_pass.R`, `projection_functions.R`, `infold_functions.R` | pooling, PCA and kPCA with out-of-sample extension | `pool_cohort()`, `fit_projection()` |
| `patient_level/pl_grid.R`, `pl_promote.R`, `run_patient_chain.sh` | the patient-level 192-cell grid with shared permutations and the full 30-step chain | `pathmi_null()`, `select_representation()` |
| `patient_level/pl_incremental.R`, `prad_incremental_batch_purity.R` | attributable recoverability | `incremental_recoverability()` |
| `usable_infold.R`, `usable_infold_stats.R` | patient-held-out usable recoverability | `usable_recoverability()` |
| `robustness_subsampling.R`, `robustness_prior_bootstrap.R`, `robustness_dimension.R` | half-sampling, bootstrap bias, dimension and prior sweeps | `half_sampling()` |
| `resolution_units.R`, `patient_level/check_refit.R` | reproducibility bound and resolution units | `reproducibility_bound()`, `replicate_refit()`, `sign_flip_refit()` |
| `learning_curves.R`, `learning_summary.R` | learning curves and N_evidence | `learning_curve()`, `evidence_status()` |
| `tile_maps_hoptimus1.R`, `sampling_burden_all_models.R`, `burden_*.R` | tile decomposition, specific surprise, tissue burden | `tile_decomposition()`, `tissue_burden()` |
| `recoverability_record.R` | Recoverability Record (CSV, JSON, table) | `recoverability_record()`, `record_to_json()` |
| `hidden_context_sim.R` | hidden-context simulation | `simulate_hidden_context()` |
| `figures_*.R`, `tables_main.R`, `*_table.R` | manuscript figures and tables | `plot_*()` |

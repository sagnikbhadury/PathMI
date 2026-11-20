#' PathMI: measuring molecular recoverability in pathology foundation-model embeddings
#'
#' The package follows one representation from tiles to a Recoverability
#' Record:
#'
#' 1. **Representation.** [pool_tiles()], [pool_cohort()] (one vector per
#'    patient), [fit_projection()], [apply_projection()], [project_embeddings()].
#' 2. **Estimate.** [pathmi_fit()] (the estimator; engine [subho_MI_logistic()]).
#' 3. **Inference.** [draw_permutations()], [pathmi_null()],
#'    [select_representation()] (max-T), [half_sampling()],
#'    [replicate_refit()], [sign_flip_refit()], [reproducibility_bound()].
#' 4. **Counterparts.** [incremental_recoverability()] (attributable),
#'    [usable_recoverability()] (patient-held-out).
#' 5. **Evidence.** [learning_curve()], [evidence_status()].
#' 6. **Tiles.** [tile_decomposition()], [tissue_burden()].
#' 7. **Record.** [recoverability_record()], [record_to_json()].
#' 8. **Simulation and baselines.** [simulate_cohort()],
#'    [simulate_hidden_context()], [simulate_estimator_design()],
#'    [doe_linear()], [run_lasso_cv()], [compute_ols_mi()].
#'
#' @keywords internal
#' @importFrom stats median quantile rgamma rnorm
#' @importFrom graphics hist
"_PACKAGE"

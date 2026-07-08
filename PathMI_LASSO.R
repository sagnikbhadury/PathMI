################################################################
# PathMI: LASSO-CV DOE Mutual Information Estimator
#
# Parallel implementation to the SUBHO horseshoe pipeline.
# Replaces:
#   - sagnik.horseshoe()  -->  cv.glmnet (logistic, alpha=1)
#   - lpy.x() marginal likelihood walk  -->  lambda.min from CV deviance
#   - step-up selection  -->  nonzero coefficients at lambda.min
#   - MI from log Bayes factor  -->  MI from log-likelihood ratio (null vs fitted)
#
# Pipeline:
#   1. Patient embedding aggregation
#   2. LASSO logistic regression with k-fold CV
#   3. DOE-type MI from fitted model log-likelihood ratio
#   4. Permutation test: reshuffle labels, rerun full pipeline
#   5. BH FDR correction across targets
#
# Dependencies: glmnet
################################################################

library(glmnet)

################################################################
# STEP 1: Aggregate tile embeddings -> patient matrix
################################################################

aggregate_embeddings <- function(embedding_list, method = "mean") {
  stopifnot(method %in% c("mean", "median"))
  agg_fn <- if (method == "mean") colMeans else function(x) apply(x, 2, median)
  Z <- do.call(rbind, lapply(embedding_list, agg_fn))
  rownames(Z) <- names(embedding_list)
  return(Z)
}

################################################################
# STEP 2: LASSO logistic with CV
#
# Returns:
#   - selected_features: indices of nonzero coefficients at lambda.min
#   - coef_lasso:        coefficient vector at lambda.min
#   - lambda_min:        selected lambda
#   - cv_fit:            full cv.glmnet object (for diagnostics)
################################################################

run_lasso_cv <- function(Z, y, nfolds = 10, seed = 42) {
  stopifnot(is.numeric(y), all(y %in% c(0, 1)))
  set.seed(seed)

  cv_fit <- cv.glmnet(
    x        = Z,
    y        = y,
    family   = "binomial",
    alpha    = 1,           # LASSO
    nfolds   = nfolds,
    type.measure = "deviance",
    standardize  = TRUE     # glmnet standardizes internally; coefficients returned on original scale
  )

  coef_lasso       <- coef(cv_fit, s = "lambda.min")
  selected_features <- which(coef_lasso[-1] != 0)  # exclude intercept row

  list(
    selected_features = selected_features,
    coef_lasso        = coef_lasso,
    lambda_min        = cv_fit$lambda.min,
    cv_fit            = cv_fit
  )
}

################################################################
# STEP 3: DOE-type MI estimation from LASSO fitted model
#
# MI is estimated as the log-likelihood ratio between the
# fitted LASSO model and the null (intercept-only) model,
# divided by n (per-sample average, nats if log base e).
#
# I(Z_sel; Y) ≈ [loglik(fitted) - loglik(null)] / n
#
# loglik from logistic model: sum_i [ y_i * log(p_i) + (1-y_i) * log(1-p_i) ]
# Null model: p_null = mean(y) for all i
#
# This is equivalent to KL divergence between joint and marginals
# under the DOE classifier framework — here approximated via
# the parametric log-likelihood ratio of the LASSO logistic fit.
################################################################

compute_doe_mi <- function(Z, y, lasso_result) {

  sel <- lasso_result$selected_features
  coef_vec <- as.vector(lasso_result$coef_lasso)  # includes intercept at position 1

  # null log-likelihood
  p_null   <- mean(y)
  p_null   <- pmax(pmin(p_null, 1 - 1e-10), 1e-10)  # numerical guard
  ll_null  <- sum(y * log(p_null) + (1 - y) * log(1 - p_null))

  # fitted log-likelihood
  if (length(sel) == 0) {
    # no features selected: MI is zero
    return(list(MI = 0, ll_fitted = ll_null, ll_null = ll_null, n_selected = 0))
  }

  Z_sel    <- cbind(1, Z[, sel, drop = FALSE])       # design matrix with intercept column
  eta      <- Z_sel %*% coef_vec[c(1, sel + 1)]     # linear predictor
  p_fitted <- plogis(eta)
  p_fitted <- pmax(pmin(p_fitted, 1 - 1e-10), 1e-10)
  ll_fit   <- sum(y * log(p_fitted) + (1 - y) * log(1 - p_fitted))

  n  <- length(y)
  MI <- (ll_fit - ll_null) / n

  list(
    MI         = MI,
    ll_fitted  = ll_fit,
    ll_null    = ll_null,
    n_selected = length(sel)
  )
}

################################################################
# STEP 4: Permutation test
#
# Shuffles patient labels B times, reruns full LASSO-CV + MI pipeline.
# Empirical p-value: proportion of null MIs >= observed MI.
################################################################

permutation_test_lasso <- function(Z, y, nfolds = 10, B = 500, seed = 42) {

  # observed
  lasso_obs <- run_lasso_cv(Z, y, nfolds = nfolds, seed = seed)
  mi_obs    <- compute_doe_mi(Z, y, lasso_obs)$MI

  # null distribution
  null_mi <- numeric(B)
  for (b in seq_len(B)) {
    set.seed(seed + b)
    y_perm    <- sample(y)
    lasso_b   <- run_lasso_cv(Z, y_perm, nfolds = nfolds, seed = seed + b)
    null_mi[b] <- compute_doe_mi(Z, y_perm, lasso_b)$MI
  }

  pval <- mean(null_mi >= mi_obs)

  list(
    MI_observed  = mi_obs,
    MI_null      = null_mi,
    pvalue       = pval,
    lasso_result = lasso_obs
  )
}

################################################################
# STEP 5: Run across multiple molecular targets + BH FDR
#
# targets: named list of binary vectors, one per molecular outcome
# embedding_list: named list of patient tile matrices
################################################################

run_PathMI_LASSO <- function(embedding_list,
                              targets,
                              nfolds      = 10,
                              B           = 500,
                              agg_method  = "mean",
                              fdr_alpha   = 0.05,
                              seed        = 42) {

  # aggregate once
  Z <- aggregate_embeddings(embedding_list, method = agg_method)
  Z <- scale(Z)  # center + scale before LASSO

  target_names <- names(targets)
  results      <- vector("list", length(target_names))
  names(results) <- target_names

  for (tgt in target_names) {
    cat(sprintf("Processing target: %s\n", tgt))
    y <- targets[[tgt]]

    # align patients
    common_pts <- intersect(rownames(Z), names(y))
    if (length(common_pts) < 20) {
      warning(sprintf("Target %s: <20 overlapping patients, skipping.", tgt))
      next
    }
    Z_tgt <- Z[common_pts, ]
    y_tgt <- y[common_pts]

    results[[tgt]] <- permutation_test_lasso(
      Z      = Z_tgt,
      y      = y_tgt,
      nfolds = nfolds,
      B      = B,
      seed   = seed
    )
  }

  # collect p-values and apply BH FDR
  mi_vals  <- sapply(results, function(r) if (!is.null(r)) r$MI_observed else NA)
  pvals    <- sapply(results, function(r) if (!is.null(r)) r$pvalue      else NA)
  pvals_adj <- p.adjust(pvals, method = "BH")

  summary_df <- data.frame(
    target      = target_names,
    MI          = mi_vals,
    pvalue      = pvals,
    pvalue_BH   = pvals_adj,
    significant = pvals_adj < fdr_alpha,
    stringsAsFactors = FALSE
  )

  list(
    summary    = summary_df,
    raw        = results
  )
}

################################################################
# DIAGNOSTICS
#
# Plot CV deviance curve for a single target (sanity check)
################################################################

plot_cv_lasso <- function(lasso_result, title = "LASSO CV") {
  plot(lasso_result$cv_fit, main = title)
  abline(v = log(lasso_result$lambda_min), col = "blue", lty = 2)
}

################################################################
# EXAMPLE USAGE
################################################################

# embedding_list <- readRDS("patient_embeddings.rds")  # list: patient -> tile x feature matrix
# targets <- list(
#   MSI    = setNames(c(1,0,1,...), patient_ids),
#   ERG    = setNames(c(0,1,0,...), patient_ids)
# )
# 
# out <- run_PathMI_LASSO(
#   embedding_list = embedding_list,
#   targets        = targets,
#   nfolds         = 10,
#   B              = 500,
#   fdr_alpha      = 0.05
# )
# 
# print(out$summary)
# plot_cv_lasso(out$raw$MSI$lasso_result, title = "MSI - LASSO CV")

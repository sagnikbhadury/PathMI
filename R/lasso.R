## PathMI: LASSO-CV DOE mutual information estimator.
##
## Parallel implementation to the SUBHO horseshoe pipeline (see subho.R).
## Replaces:
##   - sagnik.horseshoe()  -->  cv.glmnet (logistic, alpha=1)
##   - lpy.x() marginal likelihood walk  -->  lambda.min from CV deviance
##   - step-up selection  -->  nonzero coefficients at lambda.min
##   - MI from log Bayes factor  -->  MI from log-likelihood ratio (null vs fitted)

#' LASSO logistic regression with k-fold cross-validation
#'
#' Fits an L1-penalized (LASSO) logistic regression of \code{y} on \code{Z}
#' via \code{\link[glmnet]{cv.glmnet}}, selecting the model at
#' \code{lambda.min}.
#'
#' @param Z A patients x features numeric matrix (e.g. from
#'   \code{\link{aggregate_embeddings}}).
#' @param y Binary (0/1) response vector.
#' @param nfolds Number of cross-validation folds.
#' @param seed Random seed for fold assignment.
#'
#' @return A list with \code{selected_features} (indices of nonzero
#'   coefficients at \code{lambda.min}), \code{coef_lasso}, \code{lambda_min},
#'   and \code{cv_fit} (the full \code{cv.glmnet} object).
#' @export
run_lasso_cv <- function(Z, y, nfolds = 10, seed = 42) {
  stopifnot(is.numeric(y), all(y %in% c(0, 1)))
  set.seed(seed)

  cv_fit <- glmnet::cv.glmnet(
    x = Z, y = y, family = "binomial", alpha = 1,
    nfolds = nfolds, type.measure = "deviance",
    standardize = TRUE
  )

  coef_lasso <- stats::coef(cv_fit, s = "lambda.min")
  selected_features <- which(coef_lasso[-1] != 0)

  list(
    selected_features = selected_features,
    coef_lasso = coef_lasso,
    lambda_min = cv_fit$lambda.min,
    cv_fit = cv_fit
  )
}

#' DOE-type mutual information from a fitted LASSO model
#'
#' Estimates mutual information between the selected features and \code{y}
#' as the log-likelihood ratio between the fitted LASSO logistic model and
#' the null (intercept-only) model, per sample:
#' \deqn{I(Z_{sel}; Y) \approx [\ell(\hat\beta) - \ell(\hat\beta_0)] / n}
#'
#' @param Z A patients x features numeric matrix.
#' @param y Binary (0/1) response vector.
#' @param lasso_result Output of \code{\link{run_lasso_cv}}.
#'
#' @return A list with \code{MI}, \code{ll_fitted}, \code{ll_null}, and
#'   \code{n_selected}.
#' @export
compute_doe_mi <- function(Z, y, lasso_result) {

  sel <- lasso_result$selected_features
  coef_vec <- as.vector(lasso_result$coef_lasso)

  p_null <- mean(y)
  p_null <- pmax(pmin(p_null, 1 - 1e-10), 1e-10)
  ll_null <- sum(y * log(p_null) + (1 - y) * log(1 - p_null))

  if (length(sel) == 0) {
    return(list(MI = 0, ll_fitted = ll_null, ll_null = ll_null, n_selected = 0))
  }

  Z_sel <- cbind(1, Z[, sel, drop = FALSE])
  eta <- Z_sel %*% coef_vec[c(1, sel + 1)]
  p_fitted <- stats::plogis(eta)
  p_fitted <- pmax(pmin(p_fitted, 1 - 1e-10), 1e-10)
  ll_fit <- sum(y * log(p_fitted) + (1 - y) * log(1 - p_fitted))

  n <- length(y)
  MI <- (ll_fit - ll_null) / n

  list(MI = MI, ll_fitted = ll_fit, ll_null = ll_null, n_selected = length(sel))
}

#' Permutation test for the LASSO DOE mutual information estimate
#'
#' Shuffles \code{y} \code{B} times, refitting the full LASSO-CV + MI
#' pipeline each time, to build a null distribution for the observed MI.
#'
#' @inheritParams run_lasso_cv
#' @param B Number of permutations.
#'
#' @return A list with \code{MI_observed}, \code{MI_null} (length-\code{B}
#'   vector), \code{pvalue}, and \code{lasso_result} (the observed-data fit).
#' @export
permutation_test_lasso <- function(Z, y, nfolds = 10, B = 500, seed = 42) {

  lasso_obs <- run_lasso_cv(Z, y, nfolds = nfolds, seed = seed)
  mi_obs <- compute_doe_mi(Z, y, lasso_obs)$MI

  null_mi <- numeric(B)
  for (b in seq_len(B)) {
    set.seed(seed + b)
    y_perm <- sample(y)
    lasso_b <- run_lasso_cv(Z, y_perm, nfolds = nfolds, seed = seed + b)
    null_mi[b] <- compute_doe_mi(Z, y_perm, lasso_b)$MI
  }

  pval <- mean(null_mi >= mi_obs)

  list(MI_observed = mi_obs, MI_null = null_mi, pvalue = pval, lasso_result = lasso_obs)
}

#' Run the LASSO PathMI pipeline across multiple molecular targets
#'
#' Aggregates tile embeddings to the patient level, then runs
#' \code{\link{permutation_test_lasso}} independently for each target in
#' \code{targets}, applying Benjamini-Hochberg FDR correction across targets.
#'
#' @param embedding_list Named list of per-patient tile x feature matrices.
#' @param targets Named list of binary (0/1) target vectors, one per
#'   molecular outcome, each named by patient ID.
#' @param nfolds Number of cross-validation folds for the LASSO fit.
#' @param B Number of permutations per target.
#' @param agg_method Tile-aggregation method, \code{"mean"} or \code{"median"}.
#' @param fdr_alpha Significance threshold applied to BH-adjusted p-values.
#' @param seed Random seed.
#'
#' @return A list with \code{summary} (data frame of MI, p-value, and BH
#'   q-value per target) and \code{raw} (full per-target results).
#'
#' @examples
#' \dontrun{
#' embedding_list <- readRDS("patient_embeddings.rds")
#' targets <- list(MSI = setNames(c(1, 0, 1), patient_ids))
#' out <- run_PathMI_LASSO(embedding_list, targets, nfolds = 10, B = 500)
#' print(out$summary)
#' }
#'
#' @export
run_PathMI_LASSO <- function(embedding_list,
                              targets,
                              nfolds = 10,
                              B = 500,
                              agg_method = "mean",
                              fdr_alpha = 0.05,
                              seed = 42) {

  Z <- aggregate_embeddings(embedding_list, method = agg_method)
  Z <- scale(Z)

  target_names <- names(targets)
  results <- vector("list", length(target_names))
  names(results) <- target_names

  for (tgt in target_names) {
    cat(sprintf("Processing target: %s\n", tgt))
    y <- targets[[tgt]]

    common_pts <- intersect(rownames(Z), names(y))
    if (length(common_pts) < 20) {
      warning(sprintf("Target %s: <20 overlapping patients, skipping.", tgt))
      next
    }
    Z_tgt <- Z[common_pts, ]
    y_tgt <- y[common_pts]

    results[[tgt]] <- permutation_test_lasso(Z = Z_tgt, y = y_tgt, nfolds = nfolds, B = B, seed = seed)
  }

  mi_vals <- sapply(results, function(r) if (!is.null(r)) r$MI_observed else NA)
  pvals <- sapply(results, function(r) if (!is.null(r)) r$pvalue else NA)
  pvals_adj <- stats::p.adjust(pvals, method = "BH")

  summary_df <- data.frame(
    target = target_names, MI = mi_vals, pvalue = pvals, pvalue_BH = pvals_adj,
    significant = pvals_adj < fdr_alpha, stringsAsFactors = FALSE
  )

  list(summary = summary_df, raw = results)
}

#' Plot the cross-validation deviance curve for a LASSO fit
#'
#' @param lasso_result Output of \code{\link{run_lasso_cv}}.
#' @param title Plot title.
#' @return Invisibly \code{NULL}; called for its plotting side effect.
#' @export
plot_cv_lasso <- function(lasso_result, title = "LASSO CV") {
  plot(lasso_result$cv_fit, main = title)
  graphics::abline(v = log(lasso_result$lambda_min), col = "blue", lty = 2)
}

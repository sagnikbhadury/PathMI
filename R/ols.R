## PathMI: OLS-based DOE mutual information estimator.
##
## Classifier: Gaussian linear regression (OLS), no penalty. MI is extracted
## from the log-likelihood ratio between the fitted and null models, which
## for a Gaussian linear model collapses to the Gaussian channel-capacity
## formula MI = 0.5 * log(sigma2_null / sigma2_fitted) -- a monotone function
## of R-squared.
##
## NOTE: y is binary (0/1) but treated as continuous here (a linear
## probability model). This is a deliberate approximation with no
## penalization and no model selection: all features enter the model, so it
## should only be used when p << n. For p >= n, use the LASSO version
## (see lasso.R) instead.

#' DOE-type mutual information from an unpenalized OLS fit
#'
#' Fits \code{y} on all columns of \code{Z} via OLS and estimates mutual
#' information as \eqn{0.5 \log(\sigma^2_{null} / \sigma^2_{fit})}, the
#' Gaussian channel-capacity formula. Exact for Gaussian \code{y}; for binary
#' \code{y} (as in a linear probability model) this is an approximation,
#' interpretable as the MI implied by the linear model class.
#'
#' @param Z A patients x features numeric matrix, with \code{ncol(Z) < nrow(Z)}.
#' @param y Response vector (may be binary 0/1, treated as continuous).
#'
#' @return A list with \code{MI} (nats), \code{R2}, \code{sigma2_null},
#'   \code{sigma2_fit}, \code{n_features}, and \code{fit} (the raw
#'   \code{lm.fit} object).
#' @export
compute_ols_mi <- function(Z, y) {

  n <- length(y)
  p <- ncol(Z)

  if (p >= n) stop("p >= n: OLS is rank-deficient. Use the LASSO version instead.")

  sigma2_null <- stats::var(y) * (n - 1) / n

  X_fit <- cbind(1, Z)
  fit <- stats::lm.fit(X_fit, y)
  rss_fit <- sum(fit$residuals^2)
  sigma2_fit <- rss_fit / n
  sigma2_fit <- max(sigma2_fit, 1e-10)

  MI <- 0.5 * log(sigma2_null / sigma2_fit)
  MI <- max(MI, 0)

  ss_tot <- sum((y - mean(y))^2)
  R2 <- 1 - rss_fit / ss_tot

  list(MI = MI, R2 = R2, sigma2_null = sigma2_null, sigma2_fit = sigma2_fit,
       n_features = p, fit = fit)
}

#' Permutation test for the OLS DOE mutual information estimate
#'
#' @param Z A patients x features numeric matrix.
#' @param y Response vector.
#' @param B Number of permutations.
#' @param seed Random seed.
#'
#' @return A list with \code{MI_observed}, \code{MI_null}, and \code{pvalue}.
#' @export
permutation_test_ols <- function(Z, y, B = 500, seed = 42) {

  mi_obs <- compute_ols_mi(Z, y)$MI

  null_mi <- numeric(B)
  for (b in seq_len(B)) {
    set.seed(seed + b)
    y_perm <- sample(y)
    null_mi[b] <- compute_ols_mi(Z, y_perm)$MI
  }

  pval <- mean(null_mi >= mi_obs)

  list(MI_observed = mi_obs, MI_null = null_mi, pvalue = pval)
}

#' Run the OLS PathMI pipeline across multiple molecular targets
#'
#' Aggregates tile embeddings to the patient level, then runs
#' \code{\link{permutation_test_ols}} independently for each target in
#' \code{targets}, applying Benjamini-Hochberg FDR correction across targets.
#' Only appropriate when the number of features is much smaller than the
#' number of patients for every target (see \code{\link{compute_ols_mi}}).
#'
#' @inheritParams run_PathMI_LASSO
#' @param B Number of permutations per target.
#'
#' @return A list with \code{summary} (data frame of MI, p-value, and BH
#'   q-value per target) and \code{raw} (full per-target results).
#' @export
run_PathMI_OLS <- function(embedding_list,
                            targets,
                            agg_method = "mean",
                            B = 500,
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

    results[[tgt]] <- permutation_test_ols(Z = Z_tgt, y = y_tgt, B = B, seed = seed)
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

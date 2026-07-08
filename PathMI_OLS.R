################################################################
# PathMI: OLS-based DOE Mutual Information Estimator
#
# Classifier: Gaussian linear regression (OLS), no penalty.
# MI extracted from log-likelihood ratio (fitted vs null model).
#
# I(Z; Y) ≈ [loglik(fitted) - loglik(null)] / n
#
# Under Gaussian linear model:
#   loglik(fitted) = -(n/2)*log(2*pi*sigma2_hat) - RSS_fitted/(2*sigma2_hat)
#   loglik(null)   = -(n/2)*log(2*pi*sigma2_null) - RSS_null/(2*sigma2_null)
#
# Equivalently, MI ≈ (1/2) * log(sigma2_null / sigma2_fitted)
# which is the Gaussian channel capacity formula.
# This collapses to a monotone function of R².
#
# NOTE: y is binary (0/1) but treated as continuous here (LPM).
#       This is a deliberate approximation — no penalization,
#       no model selection. All features enter the model.
#       Use only when p << n. For p >= n, use LASSO version.
#
# Dependencies: base R only
################################################################

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
# STEP 2: OLS fit + DOE MI
#
# Uses lm() on all features. sigma2 estimated from residuals.
# MI = 0.5 * log(sigma2_null / sigma2_fitted)  [nats]
#
# This is exact for Gaussian Y. For binary Y (LPM), it is an
# approximation — interpretable as MI under the linear model class.
################################################################

compute_ols_mi <- function(Z, y) {

  n <- length(y)
  p <- ncol(Z)

  if (p >= n) stop("p >= n: OLS is rank-deficient. Use LASSO version instead.")

  # null model: intercept only
  sigma2_null <- var(y) * (n - 1) / n   # MLE variance = mean squared deviation from mean

  # fitted OLS model via lm.fit (avoids formula parser overhead)
  X_fit       <- cbind(1, Z)             # design matrix with intercept
  fit         <- lm.fit(X_fit, y)
  rss_fit     <- sum(fit$residuals^2)
  sigma2_fit  <- rss_fit / n             # MLE, not unbiased

  # guard against perfect fit
  sigma2_fit <- max(sigma2_fit, 1e-10)

  # DOE MI via Gaussian channel formula
  MI <- 0.5 * log(sigma2_null / sigma2_fit)
  MI <- max(MI, 0)   # numerical guard: MI >= 0

  # also return R² for interpretability
  ss_tot <- sum((y - mean(y))^2)
  R2     <- 1 - rss_fit / ss_tot

  list(
    MI         = MI,
    R2         = R2,
    sigma2_null = sigma2_null,
    sigma2_fit  = sigma2_fit,
    n_features  = p,
    fit         = fit
  )
}

################################################################
# STEP 3: Permutation test
################################################################

permutation_test_ols <- function(Z, y, B = 500, seed = 42) {

  mi_obs <- compute_ols_mi(Z, y)$MI

  null_mi <- numeric(B)
  for (b in seq_len(B)) {
    set.seed(seed + b)
    y_perm     <- sample(y)
    null_mi[b] <- compute_ols_mi(Z, y_perm)$MI
  }

  pval <- mean(null_mi >= mi_obs)

  list(
    MI_observed = mi_obs,
    MI_null     = null_mi,
    pvalue      = pval
  )
}

################################################################
# STEP 4: Run across targets + BH FDR
################################################################

run_PathMI_OLS <- function(embedding_list,
                            targets,
                            agg_method = "mean",
                            B          = 500,
                            fdr_alpha  = 0.05,
                            seed       = 42) {

  Z <- aggregate_embeddings(embedding_list, method = agg_method)
  Z <- scale(Z)

  target_names <- names(targets)
  results      <- vector("list", length(target_names))
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

    results[[tgt]] <- permutation_test_ols(
      Z    = Z_tgt,
      y    = y_tgt,
      B    = B,
      seed = seed
    )
  }

  mi_vals   <- sapply(results, function(r) if (!is.null(r)) r$MI_observed else NA)
  pvals     <- sapply(results, function(r) if (!is.null(r)) r$pvalue      else NA)
  pvals_adj <- p.adjust(pvals, method = "BH")

  summary_df <- data.frame(
    target      = target_names,
    MI          = mi_vals,
    pvalue      = pvals,
    pvalue_BH   = pvals_adj,
    significant = pvals_adj < fdr_alpha,
    stringsAsFactors = FALSE
  )

  list(summary = summary_df, raw = results)
}

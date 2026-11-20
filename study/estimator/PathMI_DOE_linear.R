################################################################
# PathMI: Classifier-based (Donsker-Varadhan) MI estimator with a
#         linear critic — "DOE Linear" in the manuscript and figures
#
# Version 2. Replaces the v1 implementation, which was structurally
# unable to return a positive estimate. See "Why v1 returned zero".
#
# Estimator (Belghazi et al. 2018; Poole et al. 2019):
#   I(Z;Y) = D_KL( p(z,y) || p(z) p(y) )
#          >= sup_T { E_joint[T] - log E_marginal[exp T] }      (DV bound)
#
# with equality at T* = log p(z,y)/(p(z)p(y)) + const. For binary Y
#   T*(z,y) = log p(y|z) - log p(y),
# which is NOT additive in z and y: it is a different function of z
# for y = 0 and y = 1. The critic must therefore carry a z-by-y
# interaction. We use the linear critic
#
#   T(z,y) = a'z + alpha*y + y * (w'z) + b                       (1)
#
# fitted by ridge-penalised logistic regression that discriminates
# joint pairs (z_i, y_i) [label 1] from product-of-marginals pairs
# (z_i, y_sigma(i)) [label 0], where sigma is a random permutation.
# With equal class sizes the fitted log-odds is a direct estimate of
# log p_joint / p_marginal, so it is plugged into the DV bound as T.
#
# Cross-fitting: the DV bound is evaluated on critic scores computed
# out of fold (K = 5). Both rows generated from patient i (its joint
# and its marginal row) are always in the same fold, so no row's
# score comes from a critic that has seen that z_i. Without this,
# in-sample scores overfit at p >= n and the "lower bound" inflates.
#
# Ridge (alpha = 0) rather than LASSO: the estimator is meant to be
# the no-selection classifier baseline, so every predictor enters
# the critic. The penalty is needed only to make the 2p + 2 critic
# parameters identifiable from 2n rows when p >= n. The penalty is
# chosen on the estimator's own objective: for every lambda on the
# glmnet ridge path the cross-fitted DV value is computed, and the
# lambda that maximises it is kept. This is the MINE principle (the
# critic is trained to maximise the bound) restricted to a convex
# critic class. Choosing lambda by classification deviance instead
# over-shrinks the critic, because joint-vs-marginal discrimination
# is a weak-signal task for binary y (a permutation leaves roughly
# half the labels unchanged), and the bound then collapses toward
# zero for reasons unrelated to the critic class.
#
# Ceiling of the critic class (checked with the oracle logit): the
# DV bound with T = y * eta(z) + alpha * y is far below the true MI
# and can be negative at strong signal, because the exact log-ratio
# log p(y|z) - log p(y) = y*eta(z) - softplus(eta(z)) - log p(y)
# contains a nonlinear term that no linear critic represents. The
# downward bias of this estimator is therefore a property of the
# linear critic under the DV bound, which is what it is meant to
# expose, not an implementation artefact.
#
# Why v1 returned zero (the bug this version fixes):
#   v1 fitted lm(label ~ Z + y), an ADDITIVE critic T = a'z + alpha*y + b.
#   The joint block and the permuted block contain the same rows of Z
#   and the same multiset of y values, so sum_i T(z_i, y_i) equals
#   sum_i T(z_i, y_sigma(i)) exactly. Hence E_joint[T] = E_marg[T],
#   and by Jensen log E_marg[exp T] >= E_marg[T] = E_joint[T], so the
#   DV value is <= 0 for every dataset and floored to 0. The zeros in
#   the earlier simulation were a property of the critic class, not
#   of the data or of the DV bound.
#
# Dependencies: glmnet
################################################################

suppressPackageStartupMessages(library(glmnet))

################################################################
# HELPER: critic design matrix  phi(z, y) = [ z , y , z * y ]
################################################################

critic_features <- function(Z, y) {
  Z <- as.matrix(Z)
  cbind(Z, y = y, Z * y)
}

################################################################
# HELPER: joint and marginal rows
#
# Joint:    (Z_i, y_i)         label 1
# Marginal: (Z_i, y_sigma(i))  label 0, sigma a random permutation
#
# Returns the stacked feature matrix, the 0/1 label, and the index
# of the originating observation i for every row (used for
# cross-fitting so that both rows of i share a fold).
################################################################

construct_doe_dataset <- function(Z, y, seed = 42) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  set.seed(seed)
  perm <- sample(n)

  X_joint <- critic_features(Z, y)
  X_marg  <- critic_features(Z, y[perm])

  list(
    X     = rbind(X_joint, X_marg),
    label = c(rep(1L, n), rep(0L, n)),
    obs   = c(seq_len(n), seq_len(n)),
    n     = n
  )
}

################################################################
# HELPER: log-mean-exp, numerically stable
################################################################

log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

################################################################
# STEP 1-3: fit critic, cross-fit scores, evaluate the DV bound
#
# Returns:
#   MI       floored DV estimate (nats)
#   MI_raw   pre-floor DV estimate; negative values indicate that
#            the critic found no separable dependence
#   E_joint  mean out-of-fold critic score on joint rows
#   log_E_marginal  log-mean-exp of out-of-fold scores on marginal rows
#   lambda   ridge penalty selected by CV
#   scores   out-of-fold critic scores, length 2n
################################################################

doe_linear_single <- function(Z, y, seed = 42, n_folds = 5, nlambda = 40) {

  stopifnot(all(y %in% c(0, 1)))
  d <- construct_doe_dataset(Z, y, seed = seed)
  X <- d$X; lab <- d$label; obs <- d$obs; n <- d$n

  # One ridge path, shared by every fold, so the DV value is a
  # function of the same lambda grid throughout.
  full    <- glmnet(X, lab, family = "binomial", alpha = 0,
                    nlambda = nlambda, standardize = TRUE)
  lambdas <- full$lambda

  # Cross-fitting: folds are assigned to observations i, and both
  # rows of i (joint + marginal) inherit that fold.
  set.seed(seed + 1)
  fold_of_obs <- sample(rep(seq_len(n_folds), length.out = n))
  fold        <- fold_of_obs[obs]

  S <- matrix(NA_real_, length(lab), length(lambdas))
  for (k in seq_len(n_folds)) {
    tr <- fold != k
    te <- fold == k
    fit_k <- glmnet(X[tr, , drop = FALSE], lab[tr], family = "binomial",
                    alpha = 0, lambda = lambdas, standardize = TRUE)
    S[te, ] <- predict(fit_k, X[te, , drop = FALSE], type = "link")
  }

  # DV value at every lambda, from out-of-fold scores only.
  dv_path <- vapply(seq_along(lambdas), function(j) {
    mean(S[lab == 1L, j]) - log_mean_exp(S[lab == 0L, j])
  }, numeric(1))

  j_best <- which.max(dv_path)
  lambda <- lambdas[j_best]
  scores <- S[, j_best]

  E_joint    <- mean(scores[lab == 1L])
  log_E_marg <- log_mean_exp(scores[lab == 0L])

  MI_raw <- E_joint - log_E_marg
  MI     <- max(MI_raw, 0)

  list(
    MI             = MI,
    MI_raw         = MI_raw,
    E_joint        = E_joint,
    log_E_marginal = log_E_marg,
    lambda         = lambda,
    lambda_path    = lambdas,
    dv_path        = dv_path,
    scores         = scores,
    n_features     = ncol(Z)
  )
}

################################################################
# STEP 4: Permutation test
#
# Under H0 (Z independent of Y) shuffling y before constructing the
# joint block removes any dependence; the null MI concentrates near
# zero. The full critic fit (ridge path, cross-fitting and the
# lambda choice on the DV objective) is repeated inside every
# permutation, so the selection step is part of what the null
# distribution measures.
################################################################

permutation_test_doe_linear <- function(Z, y, B = 500, seed = 42) {

  obs    <- doe_linear_single(Z, y, seed = seed)
  mi_obs <- obs$MI

  null_mi <- numeric(B)
  for (b in seq_len(B)) {
    set.seed(seed + b)
    y_perm     <- sample(y)
    res_b      <- doe_linear_single(Z, y_perm, seed = seed + b)
    null_mi[b] <- res_b$MI
  }

  pval <- (1 + sum(null_mi >= mi_obs)) / (1 + B)

  list(
    MI_observed = mi_obs,
    MI_raw      = obs$MI_raw,
    MI_null     = null_mi,
    pvalue      = pval,
    lambda      = obs$lambda
  )
}

################################################################
# STEP 5: Run across targets + BH FDR
################################################################

aggregate_embeddings <- function(embedding_list, method = "mean") {
  stopifnot(method %in% c("mean", "median"))
  agg_fn <- if (method == "mean") colMeans else function(x) apply(x, 2, median)
  Z <- do.call(rbind, lapply(embedding_list, agg_fn))
  rownames(Z) <- names(embedding_list)
  return(Z)
}

run_PathMI_DOE_linear <- function(embedding_list,
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

    results[[tgt]] <- permutation_test_doe_linear(
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

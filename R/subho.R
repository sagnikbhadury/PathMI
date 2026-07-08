## PathMI: SUBHO-based DOE mutual information estimator.
## Version 3 -- Polya-Gamma logistic horseshoe.
##
## Pipeline:
##   1. Patient-level embedding aggregation (mean pooling)
##   2. Logistic SUBHO: PG-augmented horseshoe + marginal likelihood walk
##      (variables.selected determined automatically from data, capped at n-1)
##   3. DOE-type MI estimate: log Bayes factor of best vs null model, /n
##   4. Permutation test: full pipeline reruns on shuffled patient labels
##   5. FDR correction across molecular targets (BH)
##
## Polya-Gamma augmentation (Polson, Scott & Windle 2013):
##   logistic regression y|beta ~ Bernoulli(sigma(X beta))
##   is made conditionally Gaussian via latent omega_i ~ PG(1, X_i beta),
##   yielding y_tilde_i = (y_i - 0.5) / omega_i ~ N(X_i beta, 1/omega_i).
##   The Gaussian horseshoe sampler then runs exactly on this augmented system.
##   Iteration alternates: (1) sample omega | beta, (2) sample beta | omega.
##
## This file is self-contained: sagnik.horseshoe.logistic(), lpy.x.logistic(),
## and selection.logistic() below are the logistic-model replacements for an
## earlier linear-model SUBHO pipeline (not part of this package) and do not
## depend on any external source() file.
##
## Changes from v2:
##   - lpy.x.logistic now returns a "diagnostics" attribute containing the
##     log-determinant and condition number of the posterior Hessian H at
##     the MAP, so Laplace breakdown can be detected per candidate model.
##   - selection.logistic collects these diagnostics across the walk
##     (Cond.H, Cond.H.null) and returns them alongside the log marginal
##     likelihoods.
##   - subho_MI_logistic now returns MI.raw (the pre-floor estimate) and
##     cond.H.best (diagnostic at k*), instead of silently flooring at 0.
##   - pathMI_pipeline surfaces MI_raw and Cond_H in the summary table and
##     emits a warning when cond(H) at k* exceeds 1e8.
##   - New: v_beta_sensitivity() reruns the full pipeline across a grid of
##     v_beta values so model selection sensitivity to this hyperparameter
##     can be assessed directly, plus a companion plotting helper.
##   - Intercept handling is unchanged from v2 and is intentional: the
##     horseshoe runs on the full design matrix including the intercept
##     column, and selection.logistic forces the intercept to rank first
##     (the 930316 override) so it is included in every candidate model in
##     the walk. This is a design choice, not a discrepancy to fix.

sigmoid <- function(x) 1 / (1 + exp(-x))

clip_eta <- function(eta, bound = 15) pmax(pmin(eta, bound), -bound)

## ----------------------------------------------------------------
## Polya-Gamma augmented Gaussian horseshoe sampler
##
##   Samples the logistic horseshoe posterior for binary y. The intercept
##   column of X is included in the horseshoe and is shrunk like every
##   other coefficient -- this is intentional, since the intercept is
##   always retained in every candidate model during selection (see
##   selection.logistic()).
## ----------------------------------------------------------------
sagnik.horseshoe.logistic <- function(X, y,
                                       Sigma_sq = 1, tau_sq = 1,
                                       method.sigma = "sparseIG",
                                       nmc = 2000, burn = 200, thin = 10) {

  N <- nmc
  effective.ss <- (N - burn) / thin
  n <- nrow(X)
  p <- ncol(X)

  algo <- if (p < 0.8 * n) 1L else 2L   # Makalic vs Bhattacharya

  Beta <- rep(0, p)
  lambda_sq <- rep(1, p)
  nu <- rep(1, p)
  my.kappa <- rep(0.5, p)
  zeta <- rep(0, p)
  xi <- 1
  sigma_sq <- Sigma_sq
  tau_sq <- tau_sq

  omega <- rep(0.25, n)
  kappa_pg <- y - 0.5

  betaout <- matrix(0, p, effective.ss)
  lambdaout1 <- lambdaout2 <- kappaout <- zetaout <- matrix(0, p, effective.ss)

  lambda_star <- tau_sq * lambda_sq

  if (algo == 1L) Q_star <- t(X) %*% X

  for (i in seq_len(N)) {

    eta <- clip_eta(as.vector(X %*% Beta))
    omega <- BayesLogit::rpg(n, 1, eta)
    omega <- pmax(omega, 1e-6)

    z <- kappa_pg / omega
    sigma_sq <- 1

    sqW <- sqrt(omega)
    Xw <- sqW * X
    zw <- sqW * z

    if (algo == 1L) {
      Qw <- t(Xw) %*% Xw
      myA <- Qw + diag(1 / lambda_star, p, p)
      inv_A <- tryCatch(solve(myA, tol = 1e-60),
                         error = function(e) solve(myA + diag(1e-8, p)))
      beta_mu <- inv_A %*% t(Xw) %*% zw
      beta_sigma <- inv_A
      Beta <- as.vector(MASS::mvrnorm(1, mu = beta_mu, Sigma = beta_sigma))
    } else {
      D1 <- diag(lambda_star, p, p)
      u1 <- rnorm(p, 0, sqrt(lambda_star))
      del1 <- rnorm(n, 0, 1)
      v1 <- Xw %*% u1 + del1
      w1 <- solve(Xw %*% D1 %*% t(Xw) + diag(n), zw - v1)
      Beta <- as.vector(u1 + D1 %*% t(Xw) %*% w1)
    }

    s.xi <- 1 + 1 / tau_sq
    xi <- 1 / rgamma(1, 1, scale = 1 / s.xi)

    for (j in seq_len(p)) {
      s.nu <- 1 + 1 / lambda_sq[j]
      nu[j] <- 1 / rgamma(1, 1, scale = 1 / s.nu)
      s.lsq <- 1 / nu[j] + Beta[j]^2 / (2 * tau_sq)
      lambda_sq[j] <- 1 / rgamma(1, 1, scale = 1 / s.lsq)
      my.kappa[j] <- lambda_sq[j] / (1 + lambda_sq[j])
    }

    # tau_sq full conditional: shape (p+1)/2 here is correct given p =
    # ncol(X) includes the intercept (i.e. p coefficients total). In the
    # derivation's notation, where p denotes predictors excluding the
    # intercept and there are p+1 total coefficients, the equivalent
    # shape is (p+2)/2 -- the two are the same formula under the two
    # different indexing conventions. Do not change this line.
    s.tau <- 1 / xi + sum(Beta^2 / lambda_sq) / 2
    tau_sq <- 1 / rgamma(1, (p + 1) / 2, scale = 1 / s.tau)

    if (algo == 1L) {
      zeta <- as.numeric(diag(inv_A %*% Qw))
    }

    lambda_star <- tau_sq * lambda_sq

    if (i %% 500 == 0) cat(sprintf("  [HS-logistic] Iteration %d\n", i))

    if (i > burn && (i - burn) %% thin == 0) {
      idx <- (i - burn) / thin
      betaout[, idx] <- Beta
      lambdaout1[, idx] <- lambda_sq
      lambdaout2[, idx] <- log(lambda_sq)
      kappaout[, idx] <- my.kappa
      zetaout[, idx] <- zeta
    }
  }

  pMean <- rowMeans(betaout)
  pMedian <- apply(betaout, 1, median)

  c.r.alpha <- c(0.90, 0.75, 0.50, 0.25, 0.10, 0.05)
  conf.int <- function(CRalpha) {
    cri <- t(apply(betaout, 1, quantile, probs = c(CRalpha / 2, 1 - CRalpha / 2)))
    pMean * (cri[, 1] * cri[, 2] > 0)
  }
  adjusted.beta <- sapply(c.r.alpha, conf.int)

  pseudo_p <- function(s) 2 * min(mean(s > 0), mean(s < 0))
  pseudo.output <- apply(betaout, 1, pseudo_p)

  list(
    BetaHat = pMean, medianBeta = pMedian, BetaSamples = betaout,
    Adj.Beta = adjusted.beta, LambdaMean = rowMeans(lambdaout1),
    LambdaLogMean = rowMeans(lambdaout2), KappaSamples = rowMeans(kappaout),
    LambdaMedian = apply(lambdaout1, 1, median), ZetaSamples = rowMeans(zetaout),
    Pseudo_pvalue = pseudo.output
  )
}

## ----------------------------------------------------------------
## Logistic marginal likelihood via Laplace approximation
##
##     log p(y|X) ~= log p(y|X, beta_hat) + log p(beta_hat)
##                  + (p/2) log(2*pi) - (1/2) log|H|
##
##   where beta_hat = MAP estimate and H = -d^2/dbeta^2 log p(y,beta)|_{beta_hat}
##   H = t(X) diag(pi*(1-pi)) X + (1/V_beta) I
##
##   v3 change: also computes and attaches the condition number of H at
##   the MAP as a "diagnostics" attribute on the return value, so callers
##   can flag unreliable Laplace approximations without recomputing anything.
## ----------------------------------------------------------------
lpy.x.logistic <- function(y, X, V_beta, ridge = 1e-6, maxit = 25, tol = 1e-8) {

  n <- nrow(X)
  p <- ncol(X)

  prior.prec <- (1 / V_beta) * diag(p)

  beta <- rep(0, p)
  it <- 0
  for (it in seq_len(maxit)) {
    eta <- as.vector(X %*% beta)
    pi_i <- pmax(pmin(sigmoid(eta), 1 - 1e-10), 1e-10)
    W <- pi_i * (1 - pi_i)
    grad <- t(X) %*% (y - pi_i) - prior.prec %*% beta
    H <- t(X) %*% diag(W) %*% X + prior.prec
    delta <- tryCatch(solve(H, grad),
                       error = function(e) solve(H + diag(ridge, p), grad))
    beta <- beta + delta
    if (max(abs(delta)) < tol) break
  }

  eta <- as.vector(X %*% beta)
  pi_i <- pmax(pmin(sigmoid(eta), 1 - 1e-10), 1e-10)
  W <- pi_i * (1 - pi_i)

  ll <- sum(y * log(pi_i) + (1 - y) * log(1 - pi_i))

  lp <- -0.5 * p * log(2 * pi * V_beta) - 0.5 * (1 / V_beta) * sum(beta^2)

  H <- t(X) %*% diag(W) %*% X + prior.prec
  ldet <- tryCatch(
    as.numeric(determinant(H, logarithm = TRUE)$modulus),
    error = function(e) 2 * sum(log(diag(chol(H + diag(ridge, p)))))
  )

  # Diagnostic: condition number of H at the MAP. Large values signal
  # near-rank-deficiency and an unreliable Laplace approximation.
  kappa.H <- tryCatch({
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    ev <- ev[ev > 0]
    if (length(ev) < 2) NA_real_ else max(ev) / min(ev)
  }, error = function(e) NA_real_)

  laplace.correction <- 0.5 * p * log(2 * pi) - 0.5 * ldet

  log.ml <- as.numeric(ll + lp + laplace.correction)

  attr(log.ml, "diagnostics") <- list(log.det.H = ldet, cond.H = kappa.H, n.iter = it)
  log.ml
}

## ----------------------------------------------------------------
## Logistic SUBHO selection
##
##   v3 change: also collects the Hessian condition number diagnostic
##   from lpy.x.logistic at every step in the walk (Cond.H) and for the
##   null model (Cond.H.null), returned alongside log marginal likelihoods.
##
##   Intercept handling (unchanged from v2, intentional): the intercept
##   is always forced to rank first via the 930316 override, so it is
##   included in every candidate model X.intercept / design.with.selected.X
##   throughout the walk. This is by design -- the intercept should always
##   be in the model. Both stages here treat the intercept identically to
##   slope coefficients in their respective priors, which is consistent.
## ----------------------------------------------------------------
selection.logistic <- function(X, y, iter,
                                v_beta = 1,
                                ig.m = 0.1,    # kept for API parity, unused
                                ig.v = 0.1,    # kept for API parity, unused
                                method.weight = "beta.binomial",
                                nmc = 2000, burn = 200, thin = 10) {

  p <- ncol(X)
  n <- nrow(X)

  variables.selected <- min(p - 1, n - 1)

  model1 <- sagnik.horseshoe.logistic(X, y, nmc = nmc, burn = burn, thin = thin)
  incpt.BetaHat <- model1$medianBeta

  # Intentional: force intercept to rank first so it is included in
  # every candidate model below. Not a bug -- see header note.
  incpt.BetaHat[iter] <- 930316

  feat.order <- order(abs(incpt.BetaHat), decreasing = TRUE)
  feat.order <- feat.order[feat.order != iter]
  p.feat <- length(feat.order)
  variables.selected <- min(p.feat, n - 1)

  X.intercept <- matrix(X[, iter], ncol = 1)

  req.likelihood <- rep(-Inf, variables.selected)
  log.ml.only <- rep(-Inf, variables.selected)
  cond.H.vec <- rep(NA_real_, variables.selected)
  design.with.selected.X <- X.intercept
  design <- NULL

  null.fit <- lpy.x.logistic(y, X.intercept, v_beta)
  log.ml.null.score <- as.numeric(null.fit)
  cond.H.null.score <- attr(null.fit, "diagnostics")$cond.H

  if (method.weight == "basic.binomial") {
    wt.mdl <- 1 / p.feat
    for (i in seq_len(variables.selected)) {
      position <- rep(0L, p)
      position[iter] <- 1L
      position[feat.order[seq_len(i)]] <- 1L

      k <- i

      design.with.selected.X <- cbind(X.intercept, X[, feat.order[seq_len(i)], drop = FALSE])
      design <- cbind(design, position)

      log.wt <- k * log(wt.mdl) + (p.feat - k) * log(1 - wt.mdl)
      fit <- lpy.x.logistic(y, design.with.selected.X, v_beta)
      log.ml <- as.numeric(fit)

      log.ml.only[i] <- log.ml
      cond.H.vec[i] <- attr(fit, "diagnostics")$cond.H
      req.likelihood[i] <- log.wt + log.ml
    }

  } else if (method.weight == "beta.binomial") {
    a <- 1; b <- 1
    for (i in seq_len(variables.selected)) {
      position <- rep(0L, p)
      position[iter] <- 1L
      position[feat.order[seq_len(i)]] <- 1L

      k <- i

      design.with.selected.X <- cbind(X.intercept, X[, feat.order[seq_len(i)], drop = FALSE])
      design <- cbind(design, position)

      log.wt <- lgamma(a + b) + lgamma(a + k) + lgamma(p.feat + b - k) -
        (lgamma(a) + lgamma(b) + lgamma(a + b + p.feat))
      fit <- lpy.x.logistic(y, design.with.selected.X, v_beta)
      log.ml <- as.numeric(fit)

      log.ml.only[i] <- log.ml
      cond.H.vec[i] <- attr(fit, "diagnostics")$cond.H
      req.likelihood[i] <- log.wt + log.ml
    }
  }

  best.model <- which.max(req.likelihood)

  if (best.model == 0L || all(req.likelihood == -Inf)) {
    graph.estimate <- sp.beta <- rep(0L, p)
    graph.estimate[iter] <- 1L
  } else {
    graph.estimate <- design[, best.model]
    sp.beta <- graph.estimate

    myX <- as.matrix(X[, graph.estimate == 1L, drop = FALSE])
    k.sel <- ncol(myX)
    beta0 <- rep(0, k.sel)
    for (it in seq_len(15)) {
      eta.m <- pmax(pmin(sigmoid(clip_eta(as.vector(myX %*% beta0))), 1 - 1e-10), 1e-10)
      W.m <- eta.m * (1 - eta.m)
      H.m <- t(myX) %*% diag(W.m) %*% myX + (1 / v_beta) * diag(k.sel)
      g.m <- t(myX) %*% (y - eta.m) - (1 / v_beta) * beta0
      delta <- tryCatch(solve(H.m, g.m), error = function(e) rep(0, k.sel))
      beta0 <- beta0 + delta
      if (max(abs(delta)) < 1e-8) break
    }
    sp.beta[graph.estimate == 1L] <- beta0
  }

  list(
    Likelihood = req.likelihood, LogML.only = log.ml.only, LogML.null = log.ml.null.score,
    Cond.H = cond.H.vec, Cond.H.null = cond.H.null.score, Design.Matrix = design,
    DX = design.with.selected.X, Est.Graph = graph.estimate, BetaHat = incpt.BetaHat,
    Adj.Beta = model1$Adj.Beta, Pseudo_pvalue = model1$Pseudo_pvalue, Sp.BetaHat = sp.beta
  )
}

#' Logistic SUBHO mutual information estimate
#'
#' Fits the Polya-Gamma logistic horseshoe model and runs the SUBHO marginal
#' likelihood walk to select the best-supported feature subset, then
#' estimates mutual information as the per-sample log Bayes factor of the
#' best model versus the null (intercept-only) model.
#'
#' @param Z A patients x features numeric matrix.
#' @param y Binary (0/1) response vector.
#' @param v_beta Prior variance on the regression coefficients.
#' @param method.weight Model-space prior over subset size:
#'   \code{"beta.binomial"} (default) or \code{"basic.binomial"}.
#' @param nmc,burn,thin MCMC settings for the horseshoe sampler.
#'
#' @return A list with \code{MI} (floored at 0), \code{MI.raw} (pre-floor
#'   estimate; a negative raw value indicates Laplace-approximation
#'   breakdown rather than absence of signal), \code{cond.H.best},
#'   \code{cond.H.null}, \code{n.selected}, \code{selected.features},
#'   \code{BetaHat}, \code{Sp.BetaHat}, \code{LogLikelihood} (the walk's log
#'   marginal likelihood + log prior surface), and \code{pseudo.pval}.
#' @export
subho_MI_logistic <- function(Z, y,
                               v_beta = 1,
                               method.weight = "beta.binomial",
                               nmc = 2000, burn = 200, thin = 10) {

  n <- nrow(Z)
  p <- ncol(Z)
  X <- cbind(1, Z)
  iter <- 1L

  sel <- selection.logistic(X = X, y = y, iter = iter, v_beta = v_beta,
                             method.weight = method.weight, nmc = nmc, burn = burn, thin = thin)

  best.model.idx <- which.max(sel$Likelihood)
  n.selected <- best.model.idx

  log.ML.best <- sel$LogML.only[best.model.idx]
  log.ML.null <- sel$LogML.null

  # Pre-floor (raw) MI estimate in nats. Always retained: a negative raw
  # value indicates Laplace breakdown rather than absence of signal, and
  # the floor below should not be the only thing reported.
  MI.raw <- (log.ML.best - log.ML.null) / n
  MI.estimate <- max(MI.raw, 0)

  list(
    MI = MI.estimate, MI.raw = MI.raw, cond.H.best = sel$Cond.H[best.model.idx],
    cond.H.null = sel$Cond.H.null, n.selected = n.selected,
    selected.features = which(sel$Est.Graph[-1] != 0), BetaHat = sel$BetaHat[-1],
    Sp.BetaHat = sel$Sp.BetaHat[-1], LogLikelihood = sel$Likelihood,
    pseudo.pval = sel$Pseudo_pvalue[-1]
  )
}

#' Permutation test for the logistic SUBHO mutual information estimate
#'
#' Shuffles \code{y} \code{n.perm} times, refitting the full horseshoe +
#' selection + MI pipeline (\code{\link{subho_MI_logistic}}) on each
#' permutation, to build a null distribution for the observed MI.
#'
#' @inheritParams subho_MI_logistic
#' @param n.perm Number of permutations.
#' @param seed Random seed.
#' @param verbose Print progress every 25 permutations.
#'
#' @return A numeric vector of length \code{n.perm} of null MI values.
#' @export
permutation_test_logistic <- function(Z, y,
                                       n.perm = 500,
                                       v_beta = 1,
                                       method.weight = "beta.binomial",
                                       nmc = 2000, burn = 200, thin = 10,
                                       seed = 42,
                                       verbose = TRUE) {
  set.seed(seed)
  MI.null <- numeric(n.perm)

  for (b in seq_len(n.perm)) {
    y.perm <- sample(y)
    res <- tryCatch(
      subho_MI_logistic(Z = Z, y = y.perm, v_beta = v_beta, method.weight = method.weight,
                         nmc = nmc, burn = burn, thin = thin),
      error = function(e) { warning(e); list(MI = NA) }
    )
    MI.null[b] <- res$MI
    if (verbose && b %% 25 == 0) cat(sprintf("  Permutation %d / %d\n", b, n.perm))
  }
  MI.null
}

#' Full multi-target PathMI pipeline (logistic SUBHO)
#'
#' Aggregates tile embeddings to the patient level, then for each target in
#' \code{targets} fits the logistic SUBHO model
#' (\code{\link{subho_MI_logistic}}), runs a permutation test
#' (\code{\link{permutation_test_logistic}}), and applies Benjamini-Hochberg
#' FDR correction across targets.
#'
#' @param embedding_list Named list of per-patient tile x feature matrices.
#' @param targets Named list of binary (0/1) target vectors, one per
#'   molecular outcome, each named by patient ID.
#' @param n.perm Number of permutations per target.
#' @param v_beta Prior variance on the regression coefficients.
#' @param method.weight Model-space prior over subset size:
#'   \code{"beta.binomial"} (default) or \code{"basic.binomial"}.
#' @param agg.method Tile-aggregation method, \code{"mean"} or \code{"median"}.
#' @param nmc,burn,thin MCMC settings for the horseshoe sampler.
#' @param fdr.method Method passed to \code{\link[stats]{p.adjust}}.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A list with \code{summary} (data frame of MI, raw MI, Hessian
#'   condition number, permutation p-value, and BH q-value per target),
#'   \code{details} (full per-target results), and \code{Z_scaled} (the
#'   scaled patient-level embedding matrix).
#'
#' @examples
#' \dontrun{
#' result_crc <- pathMI_pipeline(
#'   embedding_list = embedding_list_crc,
#'   targets = list(MSI = setNames(msi_df$status, msi_df$patient_id)),
#'   n.perm = 500, v_beta = 1, method.weight = "beta.binomial"
#' )
#' plot_MI_null(result_crc, "MSI")
#' plot_likelihood_surface(result_crc, "MSI")
#' }
#'
#' @export
pathMI_pipeline <- function(embedding_list,
                             targets,
                             n.perm = 500,
                             v_beta = 1,
                             method.weight = "beta.binomial",
                             agg.method = "mean",
                             nmc = 2000, burn = 200, thin = 10,
                             fdr.method = "BH",
                             seed = 42,
                             verbose = TRUE) {

  cat("=== PathMI: SUBHO Logistic MI Estimation ===\n\n")

  cat("Aggregating patient-level embeddings...\n")
  Z_raw <- aggregate_embeddings(embedding_list, method = agg.method)
  Z_scaled <- preprocess_embeddings(Z_raw)
  cat(sprintf("  %d patients x %d features\n\n", nrow(Z_scaled), ncol(Z_scaled)))

  n.targets <- length(targets)
  results <- vector("list", n.targets)

  for (t in seq_len(n.targets)) {
    tname <- names(targets)[t]
    y.raw <- as.numeric(targets[[t]])

    common.pts <- intersect(rownames(Z_scaled), names(y.raw))
    if (length(common.pts) < 20)
      stop(sprintf("Target '%s': fewer than 20 matched patients.", tname))

    Z.t <- Z_scaled[common.pts, , drop = FALSE]
    y.t <- y.raw[common.pts]

    cat(sprintf("--- Target: %s | n = %d | prevalence = %.2f ---\n", tname, length(y.t), mean(y.t)))
    cat(sprintf("  Feature walk ceiling: min(p-1, n-1) = %d\n", min(ncol(Z.t) - 1, length(y.t) - 1)))

    cat("  Fitting logistic SUBHO...\n")
    obs <- subho_MI_logistic(Z = Z.t, y = y.t, v_beta = v_beta, method.weight = method.weight,
                              nmc = nmc, burn = burn, thin = thin)
    cat(sprintf("  Observed MI = %.5f (raw = %.5f) | Features selected = %d | cond(H) = %.3g\n",
                obs$MI, obs$MI.raw, obs$n.selected, obs$cond.H.best))

    if (!is.na(obs$cond.H.best) && obs$cond.H.best > 1e8) {
      cat(sprintf("  WARNING: cond(H) = %.3g at k* -- Laplace approximation may be unreliable.\n",
                  obs$cond.H.best))
    }

    cat(sprintf("  Running %d permutations...\n", n.perm))
    MI.null <- permutation_test_logistic(
      Z = Z.t, y = y.t, n.perm = n.perm, v_beta = v_beta, method.weight = method.weight,
      nmc = nmc, burn = burn, thin = thin, seed = seed + t, verbose = verbose
    )

    p.val <- max(mean(MI.null >= obs$MI, na.rm = TRUE), 1 / n.perm)
    cat(sprintf("  Permutation p = %.4f\n\n", p.val))

    results[[t]] <- list(
      target = tname, MI = obs$MI, MI.raw = obs$MI.raw, cond.H.best = obs$cond.H.best,
      n.selected = obs$n.selected, selected.features = obs$selected.features,
      BetaHat = obs$BetaHat, Sp.BetaHat = obs$Sp.BetaHat, pseudo.pval = obs$pseudo.pval,
      MI.null = MI.null, p.value = p.val
    )
  }

  p.vals <- sapply(results, `[[`, "p.value")
  p.adj <- stats::p.adjust(p.vals, method = fdr.method)

  summary.df <- data.frame(
    Target = sapply(results, `[[`, "target"),
    MI = round(sapply(results, `[[`, "MI"), 5),
    MI_raw = round(sapply(results, `[[`, "MI.raw"), 5),
    Cond_H = signif(sapply(results, `[[`, "cond.H.best"), 3),
    N_selected = sapply(results, `[[`, "n.selected"),
    P_value = round(p.vals, 4),
    P_adj_BH = round(p.adj, 4),
    Significant = p.adj < 0.05,
    stringsAsFactors = FALSE
  )

  cat("=== Summary ===\n")
  print(summary.df)

  list(summary = summary.df, details = results, Z_scaled = Z_scaled)
}

#' Plot the permutation null distribution for one target
#'
#' @param result.obj Output of \code{\link{pathMI_pipeline}}.
#' @param target.name Name of the target (as used in \code{targets}) to plot.
#' @return Invisibly \code{NULL}; called for its plotting side effect.
#' @export
plot_MI_null <- function(result.obj, target.name) {
  idx <- which(sapply(result.obj$details, `[[`, "target") == target.name)
  if (!length(idx)) stop("Target not found in result object.")
  d <- result.obj$details[[idx]]

  xlim <- range(c(d$MI.null, d$MI), na.rm = TRUE)
  hist(d$MI.null, breaks = 40, col = "grey80", border = "white",
       main = paste("Permutation null --", target.name),
       xlab = "MI (permuted labels)", xlim = xlim)
  graphics::abline(v = d$MI, col = "firebrick", lwd = 2, lty = 2)
  graphics::legend("topright", bty = "n", text.col = "firebrick",
                    legend = sprintf("Observed = %.5f\np = %.4f", d$MI, d$p.value))
}

#' Plot the SUBHO log marginal likelihood surface for one target
#'
#' @inheritParams plot_MI_null
#' @return Invisibly \code{NULL}; called for its plotting side effect.
#' @export
plot_likelihood_surface <- function(result.obj, target.name) {
  idx <- which(sapply(result.obj$details, `[[`, "target") == target.name)
  if (!length(idx)) stop("Target not found.")
  d <- result.obj$details[[idx]]

  if (is.null(d$LogLikelihood)) {
    cat("LogLikelihood surface not stored. Rerun with verbose=TRUE.\n")
    return(invisible(NULL))
  }
  ll <- d$LogLikelihood
  k <- seq_along(ll)
  plot(k, ll, type = "b", pch = 19, cex = 0.6,
       col = ifelse(k == d$n.selected, "firebrick", "steelblue"),
       xlab = "Number of features (k)", ylab = "Log marginal likelihood + log prior",
       main = paste("SUBHO likelihood surface --", target.name))
  graphics::abline(v = d$n.selected, col = "firebrick", lty = 2)
  graphics::legend("topright", bty = "n", legend = sprintf("Selected k = %d", d$n.selected))
}

#' Sensitivity of the SUBHO MI estimate to the coefficient prior variance
#'
#' Reruns \code{\link{subho_MI_logistic}} across a grid of \code{v_beta}
#' values on the same data, since the Laplace log marginal likelihood
#' depends on \code{v_beta} through
#' \eqn{-(k+1)/2 \log(v_\beta) - \|\hat\beta\|^2 / (2 v_\beta)}: model
#' selection (k*) and the resulting MI estimate can both shift with
#' \code{v_beta}, so this lets that sensitivity be assessed and reported
#' alongside the main result.
#'
#' Each grid point refits the horseshoe sampler from scratch even though it
#' does not itself depend on \code{v_beta} -- for large grids or expensive
#' horseshoe fits, consider refactoring \code{selection.logistic} to accept
#' a precomputed model1 fit and ranking, so only the Laplace scoring step
#' (which does depend on \code{v_beta}) is rerun.
#'
#' @param Z A patients x features numeric matrix.
#' @param y Binary (0/1) response vector.
#' @param v_beta.grid Numeric vector of \code{v_beta} values to sweep.
#' @param method.weight Model-space prior over subset size.
#' @param nmc,burn,thin MCMC settings for the horseshoe sampler.
#' @param verbose Print progress.
#'
#' @return A data frame with one row per \code{v_beta.grid} value:
#'   \code{v_beta}, \code{MI}, \code{MI_raw}, \code{n.selected},
#'   \code{cond.H.best}.
#' @export
v_beta_sensitivity <- function(Z, y,
                                v_beta.grid = c(1, 5, 10, 25, 100),
                                method.weight = "beta.binomial",
                                nmc = 2000, burn = 200, thin = 10,
                                verbose = TRUE) {

  out <- vector("list", length(v_beta.grid))

  for (i in seq_along(v_beta.grid)) {
    vb <- v_beta.grid[i]
    if (verbose) cat(sprintf("  [v_beta sweep] v_beta = %g\n", vb))

    res <- subho_MI_logistic(Z = Z, y = y, v_beta = vb, method.weight = method.weight,
                              nmc = nmc, burn = burn, thin = thin)

    out[[i]] <- data.frame(
      v_beta = vb, MI = res$MI, MI_raw = res$MI.raw,
      n.selected = res$n.selected, cond.H.best = res$cond.H.best
    )
  }

  result <- do.call(rbind, out)
  if (verbose) {
    cat("\n  === v_beta sensitivity summary ===\n")
    print(result)
  }
  result
}

#' Plot the results of a \code{v_beta} sensitivity sweep
#'
#' @param sens.df Output of \code{\link{v_beta_sensitivity}}.
#' @return Invisibly \code{NULL}; called for its plotting side effect.
#' @export
plot_v_beta_sensitivity <- function(sens.df) {
  op <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(op), add = TRUE)
  plot(sens.df$v_beta, sens.df$MI, type = "b", pch = 19, log = "x",
       xlab = "v_beta", ylab = "MI (nats)", main = "MI vs v_beta")
  plot(sens.df$v_beta, sens.df$n.selected, type = "b", pch = 19, log = "x",
       xlab = "v_beta", ylab = "Features selected (k*)", main = "k* vs v_beta")
}

## Usable recoverability: the out-of-patient log-likelihood gain of a readout,
## with every data-dependent transformation fitted in the training patients;
## and attributable (incremental) recoverability given covariates.

strat_folds <- function(y, K, group = NULL) {
  if (is.null(group)) { f <- integer(length(y)); for (cls in 0:1) { i <- which(y == cls); f[i] <- sample(rep_len(seq_len(K), length(i))) }; return(f) }
  g <- unique(group); gy <- y[match(g, group)]; gf <- integer(length(g))
  for (cls in 0:1) { i <- which(gy == cls); gf[i] <- sample(rep_len(seq_len(K), length(i))) }
  gf[match(group, g)]
}
loglik_gain <- function(y, eta, pbar) { p <- pmin(pmax(stats::plogis(eta), 1e-12), 1 - 1e-12); y * log(p) + (1 - y) * log(1 - p) - (y * log(pbar) + (1 - y) * log(1 - pbar)) }
auc_score <- function(s, y) { r <- rank(s); n1 <- sum(y); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1)) }

#' Usable recoverability (patient-held-out)
#'
#' Cross-fitted out-of-patient log-likelihood gain of a readout,
#' \eqn{n^{-1}\sum_i \log\{\hat p_{-i}(y_i\mid x_i)/\hat p_{-i}(y_i)\}}. Within
#' each training fold the column standardisation, the projection (see
#' [fit_projection()]) and the readout are fitted; held-out patients are
#' projected with the training fit. Folds are stratified by label and hold whole
#' patients when `group` is given. Two readouts: ridge logistic regression with
#' its penalty tuned by inner cross-validation on the training patients, and the
#' PathMI readout (ranking, walk and maximum a posteriori model refitted). Under
#' a family that contains the constant predictors, the population usable
#' information cannot exceed the represented information (Xu et al., 2020).
#'
#' @param X pooled embeddings (rows = patients), before projection.
#' @param y binary labels.
#' @param decomp,d projection, see [fit_projection()].
#' @param K,R folds and repeats.
#' @param readout `"ridge"`, `"pathmi"` or both.
#' @param group optional patient ids (rows sharing an id stay in one fold).
#' @param seed base seed.
#' @param ... passed to the PathMI readout (`nmc`, `burn`, `thin`, `v_beta`).
#' @return a list with the mean usable recoverability and AUC per readout and
#'   the per-repeat values.
#' @export
usable_recoverability <- function(X, y, decomp = "pca", d = 50L, K = 5L, R = 10L, readout = c("ridge", "pathmi"), group = NULL, seed = 42L, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("glmnet is required")
  X <- as.matrix(X); y <- as.numeric(y); n <- length(y); readout <- match.arg(readout, several.ok = TRUE)
  reps <- lapply(seq_len(R), function(r) {
    set.seed(1000L + r); fold <- strat_folds(y, K, group); out <- list()
    eta <- matrix(NA_real_, n, 2, dimnames = list(NULL, c("ridge", "pathmi"))); gain <- eta
    for (k in seq_len(K)) {
      tr <- fold != k; te <- fold == k; pr <- fit_projection(X[tr, , drop = FALSE], decomp, d)
      Ztr <- apply_projection(pr, X[tr, , drop = FALSE]); Zte <- apply_projection(pr, X[te, , drop = FALSE]); pbar <- mean(y[tr])
      if ("ridge" %in% readout) {
        set.seed(seed + 100L * r + k); nf <- max(3L, min(10L, sum(y[tr] == 1), sum(y[tr] == 0)))
        cv <- glmnet::cv.glmnet(Ztr, y[tr], family = "binomial", alpha = 0, foldid = strat_folds(y[tr], nf, if (is.null(group)) NULL else group[tr]))
        e <- as.vector(stats::predict(cv, Zte, s = "lambda.min", type = "link")); eta[te, "ridge"] <- e; gain[te, "ridge"] <- loglik_gain(y[te], e, pbar)
      }
      if ("pathmi" %in% readout) {
        a <- list(...); set.seed(seed + 100L * r + k)
        sel <- selection.logistic(X = cbind(1, Ztr), y = y[tr], iter = 1L, v_beta = if (is.null(a$v_beta)) 1 else a$v_beta, method.weight = "beta.binomial",
                                  nmc = if (is.null(a$nmc)) 2000L else a$nmc, burn = if (is.null(a$burn)) 500L else a$burn, thin = if (is.null(a$thin)) 10L else a$thin)
        b <- which.max(sel$Likelihood); floored <- (sel$LogML.only[b] - sel$LogML.null) <= 0
        e <- if (floored) rep(stats::qlogis(pbar), sum(te)) else as.vector(cbind(1, Zte) %*% sel$Sp.BetaHat)
        eta[te, "pathmi"] <- e; gain[te, "pathmi"] <- loglik_gain(y[te], e, pbar)
      }
    }
    sapply(readout, function(ro) c(I = mean(gain[, ro]), auc = auc_score(eta[, ro], y)))
  })
  arr <- simplify2array(reps)
  if (length(dim(arr)) == 2) arr <- array(arr, c(2, 1, R), dimnames = list(c("I", "auc"), readout, NULL))
  list(usable = apply(arr["I", , , drop = FALSE], 2, mean), sd = apply(arr["I", , , drop = FALSE], 2, stats::sd), auc = apply(arr["auc", , , drop = FALSE], 2, mean), repeats = arr)
}

#' Attributable (incremental) recoverability given covariates
#'
#' The gain in the PathMI estimate when the representation is added to the
#' covariates, a model-based estimate of \eqn{I(Z;Y\mid C) = I(Z,C;Y) - I(C;Y)}.
#' @param Z projected representation.
#' @param C covariate matrix (standardised).
#' @param y binary labels.
#' @param ... passed to [pathmi_fit()].
#' @return a list with the covariate-only, embedding-only and combined estimates,
#'   the incremental recoverability and the covariate share.
#' @export
incremental_recoverability <- function(Z, C, y, ...) {
  C <- as.matrix(C); colnames(C) <- if (is.null(colnames(C))) paste0("C", seq_len(ncol(C))) else colnames(C)
  fc <- pathmi_fit(C, y, ...)$estimate; fe <- pathmi_fit(Z, y, ...)$estimate; fb <- pathmi_fit(cbind(Z, C), y, ...)$estimate
  list(covariates = fc, embedding = fe, combined = fb, incremental = fb - fc, share = if (fb > 0) fc / fb else NA_real_)
}

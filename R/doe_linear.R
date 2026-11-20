## Comparison estimator: the Donsker-Varadhan bound at the richest linear
## critic T(z, y) = a'z + alpha y + y w'z + b, fitted by ridge logistic
## regression to discriminate joint from permuted pairs with cross-fitting.
## The bound is loose by construction for a binary state (the exact log-ratio
## contains a softplus term no linear critic represents); it is included as
## the variational baseline of the paper's simulation.

#' Donsker-Varadhan bound with a linear critic
#'
#' @param Z numeric matrix (rows = patients).
#' @param y binary labels.
#' @param seed seed for the permutation and the folds.
#' @param n_folds cross-fitting folds (both rows of a patient share a fold).
#' @param nlambda length of the ridge path.
#' @return a list with `estimate` (floored), `raw` and the selected `lambda`.
#' @export
doe_linear <- function(Z, y, seed = 42L, n_folds = 5L, nlambda = 40L) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("glmnet is required")
  Z <- as.matrix(Z); y <- as.numeric(y); n <- nrow(Z); set.seed(seed); perm <- sample(n)
  feat <- function(yy) cbind(Z, y = yy, Z * yy)
  X <- rbind(feat(y), feat(y[perm])); lab <- rep(1:0, each = n); obs <- rep(seq_len(n), 2)
  lme <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
  full <- glmnet::glmnet(X, lab, family = "binomial", alpha = 0, nlambda = nlambda); lam <- full$lambda
  set.seed(seed + 1); fold <- sample(rep(seq_len(n_folds), length.out = n))[obs]; S <- matrix(NA_real_, length(lab), length(lam))
  for (k in seq_len(n_folds)) { tr <- fold != k; te <- fold == k
    S[te, ] <- stats::predict(glmnet::glmnet(X[tr, , drop = FALSE], lab[tr], family = "binomial", alpha = 0, lambda = lam), X[te, , drop = FALSE], type = "link") }
  dv <- vapply(seq_along(lam), function(j) mean(S[lab == 1, j]) - lme(S[lab == 0, j]), numeric(1)); j <- which.max(dv)
  list(estimate = max(dv[j], 0), raw = dv[j], lambda = lam[j])
}

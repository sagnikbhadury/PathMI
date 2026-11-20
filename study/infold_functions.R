################################################################
# infold_functions.R
#
# Patient-held-out versions of the pipeline's projections: every
# data-dependent transformation (column standardisation, PCA or kernel
# PCA, the RBF bandwidth, kernel centring, score standardisation) is
# fitted on training patients only and then applied to held-out
# patients. Used by usable_infold.R and learning_curves.R.
#
#   fit_projection(Xtr, decomp, d)  -> object
#   apply_projection(obj, X)        -> n x k standardised scores
#
# Fitting on all patients and applying to the same patients reproduces
# projection_functions.R (and so the saved Z_scaled) exactly: checked
# by every caller with align_diff().
#
# Kernel PCA out-of-sample extension (Scholkopf, Smola and Muller
# 1998): for a new point x with kernel row k(x, X_tr), centred with the
# training Gram means,
#   kc(x)_j = k(x, x_j) - mean_l k(x, x_l) - mean_l K_lj + mean(K),
# the score on component c is kc(x) v_c / sqrt(lambda_c), where
# (lambda_c, v_c) are the eigenpairs of the centred training Gram
# matrix; on training points this equals v_c sqrt(lambda_c), the
# pipeline's training scores.
################################################################
EIG_TOL <- 1e-8

fit_projection <- function(Xtr, decomp, d = 50L) {
  mu <- colMeans(Xtr); s <- sqrt(colSums(sweep(Xtr, 2, mu)^2) / (nrow(Xtr) - 1))
  keep <- is.finite(s) & s > 0
  std <- function(X) sweep(sweep(X[, keep, drop = FALSE], 2, mu[keep]), 2, s[keep], "/")
  Xs <- std(Xtr)
  if (decomp == "pca") {
    sv <- svd(Xs, nu = 0, nv = min(d, nrow(Xs), ncol(Xs)))
    lam <- sv$d^2 / (nrow(Xs) - 1)
    k <- min(d, sum(lam > lam[1] * EIG_TOL), ncol(sv$v))
    R <- sv$v[, seq_len(k), drop = FALSE]
    raw <- function(X) std(X) %*% R
  } else {
    kern <- switch(decomp,
      cosine = { nrm <- function(X) X / sqrt(rowSums(X^2) + 1e-12); A <- nrm(Xs); function(B) tcrossprod(nrm(B), A) },
      rbf    = { sq <- rowSums(Xs^2); D2 <- outer(sq, sq, "+") - 2 * tcrossprod(Xs); D2[D2 < 0] <- 0
                 sig <- max(sqrt(median(D2[upper.tri(D2)])), 1e-6)
                 function(B) { D <- outer(rowSums(B^2), sq, "+") - 2 * tcrossprod(B, Xs); D[D < 0] <- 0; exp(-D / (2 * sig^2)) } },
      poly   = { pdim <- ncol(Xs); function(B) (tcrossprod(B, Xs) / pdim + 1)^3L },
      stop("unknown decomposition ", decomp))
    K <- kern(Xs); K <- (K + t(K)) / 2
    cm <- colMeans(K); am <- mean(K)
    Kc <- sweep(sweep(K, 2, cm), 1, cm) + am; Kc <- (Kc + t(Kc)) / 2
    ev <- eigen(Kc, symmetric = TRUE); lam <- pmax(ev$values, 0); kp <- which(lam > max(lam) * EIG_TOL)
    k <- min(d, length(kp)); lt <- lam[kp[seq_len(k)]]; V <- ev$vectors[, kp[seq_len(k)], drop = FALSE]
    W <- sweep(V, 2, sqrt(lt), "/")
    raw <- function(X) { Kx <- kern(std(X)); Kx <- sweep(Kx, 2, cm) - rowMeans(Kx) + am; Kx %*% W }
  }
  Str <- raw(Xtr); m <- colMeans(Str); sd <- apply(Str, 2, sd)
  list(decomp = decomp, k = k, apply = function(X) { S <- sweep(sweep(raw(X), 2, m), 2, sd, "/"); rownames(S) <- rownames(X); S })
}
apply_projection <- function(obj, X) obj$apply(X)

# max abs difference after aligning column signs (eigenvectors are defined up to sign)
align_diff <- function(A, B) { A <- unclass(A); B <- unclass(B); s <- sign(colSums(A * B)); s[s == 0] <- 1; max(abs(sweep(A, 2, s, `*`) - B)) }

# label-stratified folds, the same construction as usable_recoverability.R
strat_folds <- function(y, K, group = NULL) {
  if (is.null(group)) { f <- integer(length(y)); for (cls in 0:1) { idx <- which(y == cls); f[idx] <- sample(rep_len(seq_len(K), length(idx))) }; return(f) }
  # grouped: every slide of a patient goes to the same fold; patients stratified by their (first) label
  g <- unique(group); gy <- y[match(g, group)]; gf <- integer(length(g))
  for (cls in 0:1) { idx <- which(gy == cls); gf[idx] <- sample(rep_len(seq_len(K), length(idx))) }
  gf[match(group, g)]
}

# label-stratified subsample of size m that keeps the cohort prevalence
strat_subsample <- function(y, m) {
  n1 <- max(2L, round(m * mean(y))); n0 <- m - n1
  sort(c(sample(which(y == 1), n1), sample(which(y == 0), n0)))
}

clip01 <- function(p) pmin(pmax(p, 1e-12), 1 - 1e-12)
loglik_gain <- function(y, eta, pbar) { p <- clip01(plogis(eta)); y * log(p) + (1 - y) * log(1 - p) - (y * log(pbar) + (1 - y) * log(1 - pbar)) }
auc <- function(s, y) { r <- rank(s); n1 <- sum(y); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1)) }

# ridge logistic on training scores, lambda by inner label-stratified CV inside the training patients
ridge_eta <- function(Ztr, ytr, Zte, seed, gtr = NULL) {
  set.seed(seed)
  nf <- max(3L, min(10L, sum(ytr == 1), sum(ytr == 0)))
  if (!is.null(gtr)) cvr <- glmnet::cv.glmnet(Ztr, ytr, family = "binomial", alpha = 0, foldid = strat_folds(ytr, nf, gtr))   # inner folds grouped by patient
  else if (nf == 10L) cvr <- glmnet::cv.glmnet(Ztr, ytr, family = "binomial", alpha = 0, nfolds = 10)
  else cvr <- glmnet::cv.glmnet(Ztr, ytr, family = "binomial", alpha = 0, foldid = strat_folds(ytr, nf))
  as.vector(predict(cvr, Zte, s = "lambda.min", type = "link"))
}

# one K-fold cross-fit with the projection refitted inside each training fold.
# Returns per-patient ridge log-likelihood gains (and the PathMI readout's if pathmi = TRUE).
crossfit_infold <- function(Xraw, y, decomp, K, fold_seed, model_seed, d = 50L, pathmi = FALSE,
                            v_beta = 1, nmc = 2000L, burn = 500L, thin = 10L, group = NULL) {
  n <- length(y); set.seed(fold_seed); fold <- strat_folds(y, K, group)
  l_r <- eta_r <- l_p <- eta_p <- rep(NA_real_, n); ks <- kp <- integer(K); fl <- logical(K)
  for (k in seq_len(K)) {
    tr <- fold != k; te <- fold == k
    pr <- fit_projection(Xraw[tr, , drop = FALSE], decomp, d); ks[k] <- pr$k
    Ztr <- apply_projection(pr, Xraw[tr, , drop = FALSE]); Zte <- apply_projection(pr, Xraw[te, , drop = FALSE])
    pbar <- mean(y[tr])
    er <- ridge_eta(Ztr, y[tr], Zte, model_seed + k, if (is.null(group)) NULL else group[tr]); eta_r[te] <- er; l_r[te] <- loglik_gain(y[te], er, pbar)
    if (pathmi) {
      set.seed(model_seed + k)
      sel <- selection.logistic(X = cbind(1, Ztr), y = y[tr], iter = 1L, v_beta = v_beta, method.weight = "beta.binomial", nmc = nmc, burn = burn, thin = thin)
      b <- which.max(sel$Likelihood); fl[k] <- (sel$LogML.only[b] - sel$LogML.null) <= 0; kp[k] <- b
      ep <- if (fl[k]) rep(qlogis(pbar), sum(te)) else as.vector(cbind(1, Zte) %*% sel$Sp.BetaHat)
      eta_p[te] <- ep; l_p[te] <- loglik_gain(y[te], ep, pbar)
    }
  }
  list(I_ridge = mean(l_r), auc_ridge = auc(eta_r, y), I_pathmi = if (pathmi) mean(l_p) else NA_real_,
       auc_pathmi = if (pathmi) auc(eta_p, y) else NA_real_, k_proj_min = min(ks), k_walk_mean = if (pathmi) mean(kp) else NA_real_,
       n_floored = if (pathmi) sum(fl) else NA_integer_)
}

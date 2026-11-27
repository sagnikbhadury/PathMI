## Inference around the PathMI estimate: the estimate itself, a permutation
## null with full refitting, max-T selection of a representation over a joint
## permutation, patient-level uncertainty by half-sampling (delete-d
## jackknife), and the empirical numerical tolerance.

#' PathMI estimate of recoverability
#'
#' Horseshoe ranking of the coordinates under Polya-Gamma augmentation, a walk
#' over the nested models scored by their Laplace marginal likelihood, and the
#' log Bayes factor of the selected model against the label's marginal, per
#' patient, floored at zero. Thin wrapper around [subho_MI_logistic()] with the
#' settings used in the paper.
#'
#' @param Z numeric matrix of standardised coordinates (rows = patients).
#' @param y binary 0/1 vector.
#' @param v_beta working prior variance of the walk.
#' @param nmc,burn,thin Gibbs iterations, burn-in and thinning.
#' @param seed seed set before the fit.
#' @return a list with `estimate` (nats), `raw` (before the floor), `k`
#'   (selected size), `selected` (coordinates), `cond_H` and `entropy` (H(Y)).
#' @export
pathmi_fit <- function(Z, y, v_beta = 1, nmc = 2000L, burn = 500L, thin = 10L, seed = 42L) {
  Z <- as.matrix(Z); y <- as.numeric(y); stopifnot(nrow(Z) == length(y), all(y %in% 0:1))
  if (!is.null(seed)) set.seed(seed)
  r <- subho_MI_logistic(Z = Z, y = y, v_beta = v_beta, method.weight = "beta.binomial", nmc = nmc, burn = burn, thin = thin)
  p <- mean(y)
  structure(list(estimate = r$MI, raw = r$MI.raw, k = r$n.selected, selected = r$selected.features, cond_H = r$cond.H.best,
                 entropy = -(p * log(p) + (1 - p) * log(1 - p)), walk = r$LogLikelihood, n = length(y)), class = "pathmi_fit")
}

#' @export
print.pathmi_fit <- function(x, ...) {
  cat(sprintf("PathMI estimate %.4f nats (%.1f%% of H(Y) = %.3f) | k* = %d | n = %d\n", x$estimate, 100 * x$estimate / x$entropy, x$entropy, x$k, x$n))
  invisible(x)
}

#' Draw label permutations shared across representations
#' @param n number of patients.
#' @param B number of permutations.
#' @param seed seed.
#' @return an n x B integer matrix of permutation indices.
#' @export
draw_permutations <- function(n, B = 200L, seed = 42L) { set.seed(seed); sapply(seq_len(B), function(b) sample.int(n)) }

#' Permutation null for one representation
#'
#' Refits the whole estimator (ranking, walk and selection) on permuted labels.
#' Pass the same `perms` to every candidate representation so that the max-T
#' selection test ([select_representation()]) uses a joint permutation.
#'
#' @inheritParams pathmi_fit
#' @param perms permutation matrix from [draw_permutations()].
#' @param observed optional `pathmi_fit` for the observed labels.
#' @param cores number of parallel workers (`parallel::mclapply` on Unix,
#'   a PSOCK cluster on Windows).
#' @return a list with `null` (vector of permuted estimates), `mean`, `sd`,
#'   `p_value`, the observed fit, the permutation matrix `perms` and the fit
#'   `settings` (used by [select_representation()] and
#'   [verify_permutation_alignment()]).
#' @export
pathmi_null <- function(Z, y, perms = draw_permutations(length(y)), observed = NULL, cores = 1L, v_beta = 1, nmc = 2000L, burn = 500L, thin = 10L) {
  Z <- as.matrix(Z); y <- as.numeric(y); force(perms); force(v_beta); force(nmc); force(burn); force(thin)
  if (is.null(observed)) observed <- pathmi_fit(Z, y, v_beta, nmc, burn, thin)
  one <- function(b) pathmi_fit(Z, y[perms[, b]], v_beta, nmc, burn, thin, seed = 100000L + b)$estimate
  nul <- unlist(par_apply(seq_len(ncol(perms)), one, cores))
  structure(list(observed = observed, null = nul, mean = mean(nul), sd = stats::sd(nul),
                 p_value = (1 + sum(nul >= observed$estimate)) / (1 + length(nul)), perms = perms,
                 settings = list(v_beta = v_beta, nmc = nmc, burn = burn, thin = thin)), class = "pathmi_null")
}

#' Select a representation with a max-T test
#'
#' Takes the candidate with the largest estimate and tests it against the
#' maximum over all candidates under each permutation,
#' \eqn{p = (1 + \#\{b : \max_c \hat I_{b,c} \ge \max_c \hat I_c\}) / (1 + B)}.
#' The candidates can be the four projections of one pooling rule (as in the
#' paper's projection selection) or all sixteen pooling x projection
#' configurations of a model and state (the paper's \eqn{p_{\max 16}}). The test
#' is valid only when every candidate was refitted on the same permuted labels,
#' which is checked. It adjusts detection within the set of candidates; it does
#' not adjust across separate sets and does not correct the magnitude of the
#' selected estimate.
#'
#' @param nulls named list of [pathmi_null()] results computed with the same
#'   permutation matrix.
#' @param check logical; stop if the stored permutation matrices differ.
#' @return a list with the selected name, its estimate, all estimates, the
#'   max-T p value `p_maxT`, the number of candidates and the null maxima `T_null`.
#' @export
select_representation <- function(nulls, check = TRUE) {
  if (is.null(names(nulls)) || any(names(nulls) == "")) names(nulls) <- paste0("candidate", seq_along(nulls))
  if (check) {
    pm <- lapply(nulls, `[[`, "perms")
    if (any(vapply(pm, is.null, logical(1))) || !all(vapply(pm, identical, logical(1), pm[[1]])))
      stop("the candidates were not refitted on the same permuted labels; pass one draw_permutations() matrix to every pathmi_null()")
  }
  est <- vapply(nulls, function(x) x$observed$estimate, numeric(1)); best <- names(which.max(est))
  Tn <- apply(do.call(cbind, lapply(nulls, `[[`, "null")), 1, max)
  list(selected = best, estimate = est[[best]], estimates = est, p_maxT = (1 + sum(Tn >= est[[best]])) / (1 + length(Tn)),
       n_candidates = length(nulls), T_null = Tn)
}

#' Verify that permuted fits used the stored permutations
#'
#' Refits the estimator on a few permutations of the stored matrix, with the
#' seed [pathmi_null()] used for each, and compares the result with the stored
#' permuted estimate. Exact agreement shows that the stored null of every
#' checked candidate was computed on exactly these labels, so that a max-T test
#' over the candidates is a joint permutation test.
#'
#' @param nulls named list of [pathmi_null()] results.
#' @param Zs named list of the corresponding representations.
#' @param y labels.
#' @param b permutation indices to refit.
#' @return a data frame with the candidate, `b`, the stored and refitted
#'   estimates and their absolute difference.
#' @export
verify_permutation_alignment <- function(nulls, Zs, y, b = c(1L, 2L)) {
  y <- as.numeric(y)
  do.call(rbind, lapply(names(nulls), function(k) { nl <- nulls[[k]]; st <- nl$settings
    do.call(rbind, lapply(b, function(bb) {
      r <- pathmi_fit(Zs[[k]], y[nl$perms[, bb]], st$v_beta, st$nmc, st$burn, st$thin, seed = 100000L + bb)$estimate
      data.frame(candidate = k, b = bb, stored = nl$null[bb], refit = r, abs_diff = abs(r - nl$null[bb])) })) }))
}

#' Patient-level uncertainty by half-sampling
#'
#' Refits the estimator on `B` random halves of the patients (without
#' replacement) and returns the delete-d jackknife standard error
#' (Shao and Wu, 1989), which for half-samples is the standard deviation of the
#' half-sample estimates. The bootstrap with replacement is not used: duplicated
#' patients double-count evidence in the marginal likelihood and bias the
#' estimate upward.
#'
#' @inheritParams pathmi_fit
#' @param B number of half-samples.
#' @param level coverage of the reported interval.
#' @param cores parallel workers.
#' @return a list with `se`, `interval`, `draws` and `estimate`.
#' @export
half_sampling <- function(Z, y, B = 200L, level = 0.90, cores = 1L, seed = 777L, v_beta = 1, nmc = 2000L, burn = 500L, thin = 10L) {
  Z <- as.matrix(Z); y <- as.numeric(y); force(seed); force(v_beta); force(nmc); force(burn); force(thin)
  n <- length(y); m <- floor(n / 2); fit <- pathmi_fit(Z, y, v_beta, nmc, burn, thin)
  one <- function(b) { set.seed(seed + b); idx <- sample.int(n, m); pathmi_fit(Z[idx, , drop = FALSE], y[idx], v_beta, nmc, burn, thin, seed = 30000L + b)$estimate }
  dr <- unlist(par_apply(seq_len(B), one, cores)); se <- stats::sd(dr) * sqrt(m / (n - m)); z <- stats::qnorm(1 - (1 - level) / 2)
  list(estimate = fit$estimate, se = se, interval = c(max(fit$estimate - z * se, 0), fit$estimate + z * se), draws = dr, m = m)
}

#' Empirical numerical tolerance of the protocol
#'
#' The tolerance \eqn{\epsilon_{tol}} is the largest discrepancy between an
#' original fit and refits that should return the same estimate: a replicate
#' refit with an independent random-number stream ([replicate_refit()]) and a
#' refit of an equivalent representation ([sign_flip_refit()]). It therefore
#' includes the Monte Carlo error of the sampler. Contrasts are reported in
#' units of it, and a record is resolved when its lower bound exceeds the null
#' plus the tolerance ([evidence_status()]).
#'
#' A refit after adding a covariate block that carries no information
#' ([null_block_refit()]) changes the model. It is a sensitivity to model
#' specification, is passed as `specification` and is reported separately,
#' outside the tolerance. The tolerance is empirical and protocol-specific: it
#' is not a universal error bound, and it does not cover fits that are unstable
#' across chain seeds, which should be flagged separately. Compute it per model
#' group when the groups differ (in the paper, 0.0163 nats for three
#' vision-only models and 0.0215 for CONCH v1.5).
#'
#' @param ... numeric vectors of discrepancies (nats), one per numerical identity.
#' @param specification optional numeric vector of discrepancies from a
#'   specification check (reported, not included).
#' @return a list with the `tolerance` (also as `bound`), its 95th percentile,
#'   per-identity maxima, counts of exact agreements, and the `specification`
#'   summary.
#' @export
numerical_tolerance <- function(..., specification = NULL) {
  ids <- list(...); if (!length(ids)) stop("give at least one numerical identity")
  all <- abs(unlist(ids)); tol <- max(all)
  spec <- if (is.null(specification)) NULL else list(max = max(abs(specification)), exact = sum(abs(specification) < 1e-9), n = length(specification))
  list(tolerance = tol, bound = tol, q95 = unname(stats::quantile(all, 0.95)), per_identity = vapply(ids, function(d) max(abs(d)), numeric(1)),
       exact = vapply(ids, function(d) sum(abs(d) < 1e-9), numeric(1)), n = lengths(ids), specification = spec)
}

#' Reproducibility bound (deprecated)
#'
#' Superseded by [numerical_tolerance()]. An argument named `null_block` or
#' `specification` is treated as a specification check and kept out of the
#' tolerance; every other argument is a numerical identity.
#' @param ... numeric vectors of discrepancies.
#' @return as [numerical_tolerance()].
#' @export
reproducibility_bound <- function(...) {
  ids <- list(...); nm <- names(ids); if (is.null(nm)) nm <- rep("", length(ids))
  spec <- nm %in% c("null_block", "specification")
  do.call(numerical_tolerance, c(ids[!spec], list(specification = if (any(spec)) unlist(ids[spec]) else NULL)))
}

#' Refit checks for the numerical tolerance
#'
#' `replicate_refit()` refits with an independent seed; `sign_flip_refit()`
#' flips the sign of every other coordinate (an equivalent representation) and
#' refits. These two are numerical identities and define
#' [numerical_tolerance()]. `null_block_refit()` appends `q` columns drawn
#' independently of the label, which carry no information, and refits; it
#' changes the model and is a specification check, reported separately. Each
#' returns the discrepancy from the original estimate.
#' @inheritParams pathmi_fit
#' @param fit the original [pathmi_fit()].
#' @param ... passed to [pathmi_fit()].
#' @export
replicate_refit <- function(Z, y, fit, seed = 4242L, ...) pathmi_fit(Z, y, seed = seed, ...)$estimate - fit$estimate
#' @rdname replicate_refit
#' @export
sign_flip_refit <- function(Z, y, fit, ...) { s <- rep(c(1, -1), length.out = ncol(Z)); pathmi_fit(sweep(as.matrix(Z), 2, s, `*`), y, ...)$estimate - fit$estimate }

#' @rdname replicate_refit
#' @param q number of appended null columns.
#' @param seed seed for the refit or the null columns.
#' @export
null_block_refit <- function(Z, y, fit, q = 5L, seed = 777L, ...) {
  set.seed(seed); N0 <- matrix(stats::rnorm(nrow(Z) * q), nrow(Z), q); colnames(N0) <- paste0("null", seq_len(q))
  pathmi_fit(cbind(as.matrix(Z), N0), y, ...)$estimate - fit$estimate
}

# internal: lapply over a cluster or forked workers
par_apply <- function(X, FUN, cores = 1L) {
  if (cores <= 1L) return(lapply(X, FUN))
  if (.Platform$OS.type == "unix") return(parallel::mclapply(X, FUN, mc.cores = cores))
  cl <- parallel::makeCluster(cores); on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(PathMI)))
  parallel::parLapply(cl, X, FUN)
}

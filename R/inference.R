## Inference around the PathMI estimate: the estimate itself, a permutation
## null with full refitting, max-T selection of a representation over a joint
## permutation, patient-level uncertainty by half-sampling (delete-d
## jackknife), and the reproducibility bound.

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
#'   `p_value` and the observed fit.
#' @export
pathmi_null <- function(Z, y, perms = draw_permutations(length(y)), observed = NULL, cores = 1L, v_beta = 1, nmc = 2000L, burn = 500L, thin = 10L) {
  Z <- as.matrix(Z); y <- as.numeric(y); force(perms); force(v_beta); force(nmc); force(burn); force(thin)
  if (is.null(observed)) observed <- pathmi_fit(Z, y, v_beta, nmc, burn, thin)
  one <- function(b) pathmi_fit(Z, y[perms[, b]], v_beta, nmc, burn, thin, seed = 100000L + b)$estimate
  nul <- unlist(par_apply(seq_len(ncol(perms)), one, cores))
  structure(list(observed = observed, null = nul, mean = mean(nul), sd = stats::sd(nul),
                 p_value = (1 + sum(nul >= observed$estimate)) / (1 + length(nul))), class = "pathmi_null")
}

#' Select a representation with a max-T test
#'
#' @param nulls named list of [pathmi_null()] results computed with the same
#'   permutation matrix.
#' @return a list with the selected name, its estimate and the max-T p value,
#'   which controls the familywise error over the candidates.
#' @export
select_representation <- function(nulls) {
  est <- vapply(nulls, function(x) x$observed$estimate, numeric(1)); best <- names(which.max(est))
  Tn <- apply(do.call(cbind, lapply(nulls, `[[`, "null")), 1, max)
  list(selected = best, estimate = est[[best]], estimates = est, p_maxT = (1 + sum(Tn >= est[[best]])) / (1 + length(Tn)))
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

#' Reproducibility bound from quantities required to be identical
#'
#' The bound is the largest discrepancy among quantities the design requires to
#' be identical, for example (i) a replicate refit with an independent
#' random-number stream, (ii) the estimate after adding a covariate block that
#' carries no information, and (iii) a refit after the projection is recomputed
#' by an independent routine (coordinate signs may flip). It therefore includes
#' the Monte Carlo error of the sampler. Contrasts are reported in units of this
#' bound.
#'
#' @param ... numeric vectors of discrepancies (nats), one per identity.
#' @return a list with the `bound`, the 95th percentile and per-identity maxima.
#' @export
reproducibility_bound <- function(...) {
  ids <- list(...); all <- abs(unlist(ids))
  list(bound = max(all), q95 = unname(stats::quantile(all, 0.95)), per_identity = vapply(ids, function(d) max(abs(d)), numeric(1)),
       exact = vapply(ids, function(d) sum(abs(d) < 1e-9), numeric(1)), n = lengths(ids))
}

#' Identity checks for the reproducibility bound
#'
#' `replicate_refit()` refits with an independent seed; `sign_flip_refit()`
#' flips the sign of every other coordinate (an equivalent representation) and
#' refits; `null_block_refit()` appends `q` columns drawn independently of the
#' label, which carry no information, and refits. Each returns the discrepancy
#' from the original estimate.
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

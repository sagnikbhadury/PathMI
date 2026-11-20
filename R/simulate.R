## Simulated data for examples and for the two simulation studies of the paper.

#' Simulate a cohort of slides with tile embeddings
#'
#' Each patient has one slide, or two to four slides with probability
#' `p_multi`; each slide is a square grid of about `tiles_mean` tiles. A
#' contiguous region of the slide is tumour (its share varies between slides
#' around `tumour_frac`, as purity does) and the rest is stroma. Tile embeddings
#' are Gaussian around a tissue-type centre with a patient-level shift, and the
#' tumour tiles of positive patients are shifted by `effect` along `n_signal`
#' random directions, so the state leaves a trace in part of the tissue only. A
#' patient-level covariate (for example stage) is correlated with the state with
#' strength `covariate_effect`. No real embeddings are used.
#'
#' @param n_patients number of patients.
#' @param prevalence probability of the positive state.
#' @param D embedding width.
#' @param tiles_mean mean tiles per slide.
#' @param tumour_frac mean share of tumour tiles.
#' @param effect size of the shift for positive tumour tiles.
#' @param n_signal number of directions carrying the shift.
#' @param p_multi probability that a patient has several slides.
#' @param covariate_effect mean difference of the covariate between states.
#' @param tile_sd standard deviation of the tile-level noise.
#' @param seed seed.
#' @return a list with `tiles` (list of matrices, one per slide), `coords` (tile
#'   positions per slide), `tissue` (tile types per slide), `patient` (patient
#'   id of each slide), `y` (labels named by patient) and `covariate` (named by
#'   patient).
#' @examples
#' sim <- simulate_cohort(60, D = 16, tiles_mean = 30)
#' X <- pool_cohort(sim$tiles, sim$patient)
#' dim(X)
#' @export
simulate_cohort <- function(n_patients = 300L, prevalence = 0.3, D = 64L, tiles_mean = 200L, tumour_frac = 0.4, effect = 0.35,
                            n_signal = 4L, p_multi = 0.05, covariate_effect = 0.8, tile_sd = 1, seed = 1L) {
  set.seed(seed)
  centres <- list(tumour = stats::rnorm(D), stroma = stats::rnorm(D))
  dirs <- qr.Q(qr(matrix(stats::rnorm(D * n_signal), D, n_signal)))
  y <- stats::rbinom(n_patients, 1, prevalence); ids <- sprintf("P%04d", seq_len(n_patients)); names(y) <- ids
  covariate <- stats::rnorm(n_patients, covariate_effect * y); names(covariate) <- ids
  ns <- ifelse(stats::runif(n_patients) < p_multi, sample(2:4, n_patients, replace = TRUE), 1L)
  patient <- rep(ids, ns); tiles <- vector("list", length(patient)); tissue <- tiles; coords <- tiles
  pat_shift <- matrix(stats::rnorm(n_patients * D, sd = 0.3), n_patients, D)
  for (s in seq_along(patient)) {
    i <- match(patient[s], ids); w <- max(5L, round(sqrt(stats::rpois(1, tiles_mean)))); h <- max(5L, round(w * stats::runif(1, 0.7, 1.3)))
    xy <- as.matrix(expand.grid(x = seq_len(w), y = seq_len(h))); Tn <- nrow(xy)
    frac <- min(0.9, max(0.1, tumour_frac * stats::runif(1, 0.5, 1.5))); ctr <- c(stats::runif(1, 0.3, 0.7) * w, stats::runif(1, 0.3, 0.7) * h)
    dd <- sqrt(((xy[, 1] - ctr[1]) / w)^2 + ((xy[, 2] - ctr[2]) / h)^2) + stats::rnorm(Tn, sd = 0.04)
    tt <- ifelse(rank(dd) <= round(frac * Tn), "tumour", "stroma")
    E <- t(vapply(tt, function(t) centres[[t]], numeric(D))) + matrix(stats::rnorm(Tn * D, sd = tile_sd), Tn, D) + matrix(pat_shift[i, ], Tn, D, byrow = TRUE)
    if (y[i] == 1) E[tt == "tumour", ] <- E[tt == "tumour", , drop = FALSE] + matrix(effect * rowSums(dirs), sum(tt == "tumour"), D, byrow = TRUE)
    tiles[[s]] <- unname(E); tissue[[s]] <- unname(tt); coords[[s]] <- unname(xy)
  }
  list(tiles = tiles, coords = coords, tissue = tissue, patient = patient, y = y, covariate = covariate)
}

#' Simulate the hidden-context design
#'
#' Two contexts that the estimator does not observe, with the same distribution
#' of a `p`-dimensional representation (AR(1) correlation `rho`). In context A,
#' \eqn{\mathrm{logit}\,P(Y=1\mid Z) = b_0 + Z\beta}; in context B, Y is
#' independent of Z at the prevalence of context A. Because Z is independent of
#' S, \eqn{I(Z;Y) \le I(Z;Y\mid S) = \pi_A I(Z;Y\mid S=A)}: a positive pooled
#' estimate can conceal a context without information.
#' @param n total patients.
#' @param pi_A share of context A.
#' @param p dimension.
#' @param rho AR(1) correlation.
#' @param b0,beta intercept and coefficients in context A (`beta` padded with zeros).
#' @param seed seed.
#' @param mc Monte Carlo draws for the true information (0 to skip).
#' @return a list with `Z`, `y`, `context` and, if `mc > 0`, `truth`
#'   (`I_A`, `I_pooled`, `I_conditional`).
#' @export
simulate_hidden_context <- function(n = 400L, pi_A = 0.5, p = 50L, rho = 0.3, b0 = -0.8, beta = c(1.4, -1.1, 0.8), seed = 1L, mc = 2e5) {
  L <- chol(rho^abs(outer(seq_len(p), seq_len(p), "-"))); bb <- c(beta, rep(0, p - length(beta)))
  hb <- function(q) { q <- pmin(pmax(q, 1e-300), 1 - 1e-16); -(q * log(q) + (1 - q) * log(1 - q)) }
  truth <- NULL; pA <- NULL
  if (mc > 0) { set.seed(seed + 99L); Zm <- matrix(stats::rnorm(mc * p), mc, p) %*% L; pz <- stats::plogis(b0 + as.vector(Zm %*% bb)); pA <- mean(pz)
    truth <- c(I_A = hb(pA) - mean(hb(pz)), I_pooled = hb(pA) - mean(hb(pi_A * pz + (1 - pi_A) * pA)), I_conditional = pi_A * (hb(pA) - mean(hb(pz)))) }
  if (is.null(pA)) pA <- 0.37
  set.seed(seed); nA <- round(pi_A * n); Z <- matrix(stats::rnorm(n * p), n, p) %*% L; S <- rep(c("A", "B"), c(nA, n - nA))
  y <- ifelse(S == "A", stats::rbinom(n, 1, stats::plogis(b0 + as.vector(Z %*% bb))), stats::rbinom(n, 1, pA))
  list(Z = scale(Z), y = y, context = S, truth = truth)
}

#' Simulate the estimator-comparison design
#'
#' \eqn{Z \sim N_p(0,\Sigma)} with \eqn{\Sigma_{jk} = 0.3} off the diagonal and
#' \eqn{\mathrm{logit}\,P(Y=1\mid Z) = Z\beta} with `s` nonzero coefficients of
#' size `beta` (pattern 1, -1, 0.5, -0.5, 0.8 repeated). The true information is
#' computed by Monte Carlo.
#' @param n,p,s,beta design (patients, dimension, nonzero coefficients, signal size).
#' @param seed seed.
#' @param mc Monte Carlo draws.
#' @return a list with `Z`, `y` and `truth` (nats).
#' @export
simulate_estimator_design <- function(n = 200L, p = 200L, s = 5L, beta = 0.7, seed = 1L, mc = 5e4) {
  Sig <- matrix(0.3, p, p); diag(Sig) <- 1; L <- chol(Sig); b <- rep(0, p); b[seq_len(s)] <- beta * rep_len(c(1, -1, 0.5, -0.5, 0.8), s)
  hb <- function(q) { q <- pmin(pmax(q, 1e-300), 1 - 1e-16); -(q * log(q) + (1 - q) * log(1 - q)) }
  set.seed(seed + 7L); Zm <- matrix(stats::rnorm(mc * p), mc, p) %*% L; pz <- stats::plogis(as.vector(Zm %*% b)); truth <- hb(mean(pz)) - mean(hb(pz))
  set.seed(seed); Z <- matrix(stats::rnorm(n * p), n, p) %*% L; y <- stats::rbinom(n, 1, stats::plogis(as.vector(Z %*% b)))
  list(Z = Z, y = y, truth = truth)
}

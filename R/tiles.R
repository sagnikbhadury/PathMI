## Tile-level decomposition of the linear arm (mean pooling, linear PCA) and
## the tissue-sampling burden.

#' Exact tile-level decomposition of the linear arm
#'
#' Under mean pooling with linear PCA the fitted logit of a patient is an affine
#' function of the mean of its tile embeddings, so every tile contributes a
#' scalar \eqn{g(z)} and the patient's logit is the average over its slides of
#' each slide's tile mean of \eqn{a_0 + g(z)}. For every tile the function
#' returns the contribution, the as-if logit \eqn{a_0 + g(z)} and its specific
#' surprise \eqn{\mathrm{KL}\{\mathrm{Bern}(\sigma(\eta))\,\|\,\mathrm{Bern}(\bar p)\}}
#' (DeWeese and Meister, 1999).
#'
#' @param tile_list list of tile matrices, one per slide.
#' @param patient patient id of each slide.
#' @param y labels, one per patient, named by patient id (or in order of first
#'   appearance in `patient`).
#' @param d retained coordinates.
#' @param ... passed to [pathmi_fit()].
#' @return a list with `tiles` (data frame: slide, patient, tile, g, eta,
#'   surprise), `patients` (logit, fitted probability, check of the exact
#'   identity) and the `fit`.
#' @export
tile_decomposition <- function(tile_list, patient, y, d = 50L, ...) {
  X <- pool_cohort(tile_list, patient, "mean"); units <- rownames(X)
  if (!is.null(names(y))) y <- y[units]; y <- as.numeric(y)
  pca <- stats::prcomp(X, center = TRUE, scale. = TRUE); k <- min(d, ncol(pca$x), sum(pca$sdev^2 > pca$sdev[1]^2 * 1e-8))
  PC <- pca$x[, seq_len(k), drop = FALSE]; m_k <- colMeans(PC); s_k <- apply(PC, 2, stats::sd); Z <- sweep(sweep(PC, 2, m_k), 2, s_k, "/")
  fit <- pathmi_fit(Z, y, ...); sel <- fit$selected
  Xs <- cbind(1, Z[, sel, drop = FALSE]); beta <- rep(0, ncol(Xs))
  for (it in 1:50) { p <- stats::plogis(Xs %*% beta); W <- as.vector(p * (1 - p)); H <- crossprod(Xs * W, Xs) + diag(ncol(Xs))
    g <- crossprod(Xs, y - p) - beta; dlt <- solve(H, g); beta <- beta + as.vector(dlt); if (max(abs(dlt)) < 1e-10) break }
  b0 <- beta[1]; bs <- beta[-1]; v <- pca$rotation[, sel, drop = FALSE] %*% (bs / s_k[sel]); a0 <- b0 - sum(bs * m_k[sel] / s_k[sel])
  pbar <- mean(y); kl <- function(p, q) { p <- pmin(pmax(p, 1e-12), 1 - 1e-12); p * log(p / q) + (1 - p) * log((1 - p) / (1 - q)) }
  tiles <- do.call(rbind, lapply(seq_along(tile_list), function(s) { E <- sweep(sweep(as.matrix(tile_list[[s]]), 2, pca$center), 2, pca$scale, "/")
    gg <- as.vector(E %*% v); data.frame(slide = s, patient = patient[s], tile = seq_along(gg), g = gg, eta = a0 + gg, surprise = kl(stats::plogis(a0 + gg), pbar)) }))
  eta_pat <- as.vector(Xs %*% beta); names(eta_pat) <- units
  slide_mean <- tapply(tiles$eta, tiles$slide, mean); pat_from_tiles <- tapply(slide_mean, patient, mean)[units]
  list(tiles = tiles, patients = data.frame(patient = units, y = y, eta = eta_pat, p = stats::plogis(eta_pat), eta_from_tiles = as.vector(pat_from_tiles)),
       max_identity_error = max(abs(eta_pat - pat_from_tiles)), fit = fit, prevalence = pbar)
}

#' Tissue-sampling burden
#'
#' Holding the linear-arm model of [tile_decomposition()] fixed, recomputes
#' every patient's logit from `m` tiles drawn at random per slide and returns the
#' plug-in estimate relative to its whole-slide value, and \eqn{N_q}, the fewest
#' tiles per slide at which the ratio reaches `q`.
#'
#' @param decomp a [tile_decomposition()] result.
#' @param m_grid tiles per slide.
#' @param R random draws per `m`.
#' @param q levels.
#' @param seed seed.
#' @return a list with the `curve` and `N` (named by level).
#' @export
tissue_burden <- function(decomp, m_grid = c(5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200), R = 20L, q = c(0.90, 0.95), seed = 1L) {
  tl <- decomp$tiles; pt <- decomp$patients; pbar <- decomp$prevalence; a0 <- 0
  plug <- function(eta) { p <- pmin(pmax(stats::plogis(eta), 1e-12), 1 - 1e-12); mean(pt$y * log(p) + (1 - pt$y) * log(1 - p) - (pt$y * log(pbar) + (1 - pt$y) * log(1 - pbar))) }
  full <- plug(pt$eta_from_tiles); set.seed(seed); by_slide <- split(tl$eta, tl$slide); slide_pat <- tapply(tl$patient, tl$slide, `[`, 1)
  curve <- t(sapply(m_grid, function(m) { r <- replicate(R, { se <- vapply(by_slide, function(e) mean(if (length(e) > m) sample(e, m) else e), numeric(1))
    plug(tapply(se, slide_pat, mean)[pt$patient]) }); c(m = m, mean = mean(r) / full, lo = stats::quantile(r, 0.05) / full, hi = stats::quantile(r, 0.95) / full) }))
  colnames(curve) <- c("m", "ratio", "lo", "hi")
  N <- vapply(q, function(qq) { ok <- curve[, "ratio"] >= qq; if (any(ok)) curve[which(ok)[1], "m"] else NA_real_ }, numeric(1)); names(N) <- paste0("N", round(100 * q))
  list(curve = as.data.frame(curve), N = N, full = full)
}

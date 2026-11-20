## Patient-level representation: pooling of tile embeddings, averaging of a
## patient's slides, and projection to d coordinates by linear PCA or kernel
## PCA, with an out-of-sample extension so that every data-dependent step can
## be fitted on training patients only.

#' Pool the tile embeddings of one slide
#'
#' @param tiles numeric matrix, one row per tile and one column per embedding
#'   coordinate.
#' @param method one of `"mean"`, `"norm"` (tiles weighted by their L2 norm),
#'   `"moment"` (mean and elementwise standard deviation, 2D coordinates) or
#'   `"attention"` (untrained attention with weights drawn once from the
#'   Kaiming uniform distribution, hidden size 128).
#' @param seed seed for the untrained attention weights (identical for every
#'   slide when the seed is fixed).
#' @param hidden hidden size of the attention scorer.
#' @return a numeric vector.
#' @examples
#' pool_tiles(matrix(rnorm(40), 10, 4), "moment")
#' @export
pool_tiles <- function(tiles, method = c("mean", "norm", "moment", "attention"), seed = 42L, hidden = 128L) {
  method <- match.arg(method); tiles <- as.matrix(tiles)
  switch(method,
    mean = colMeans(tiles),
    norm = { w <- sqrt(rowSums(tiles^2)); colSums(tiles * w) / sum(w) },
    moment = c(colMeans(tiles), if (nrow(tiles) > 1) apply(tiles, 2, stats::sd) else rep(0, ncol(tiles))),
    attention = {
      D <- ncol(tiles); old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
      set.seed(seed); W <- matrix(stats::runif(D * hidden, -1 / sqrt(D), 1 / sqrt(D)), D, hidden); v <- stats::runif(hidden, -1 / sqrt(hidden), 1 / sqrt(hidden))
      if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
      s <- as.vector(tanh(tiles %*% W) %*% v); a <- exp(s - max(s)); a <- a / sum(a)
      colSums(tiles * a)
    })
}

#' Pool a cohort of slides into one vector per patient
#'
#' Each slide is pooled with [pool_tiles()], and the slide-level vectors of a
#' patient with several slides are averaged, so that the unit of analysis is the
#' patient. Counting a patient's slides as separate observations double-counts
#' evidence in the marginal likelihood and leaks information across folds.
#'
#' @param tile_list list of tile matrices, one per slide.
#' @param patient character vector giving the patient of each slide.
#' @param method pooling rule, see [pool_tiles()].
#' @param ... passed to [pool_tiles()].
#' @return a numeric matrix with one row per patient (row names = patient ids,
#'   in order of first appearance).
#' @export
pool_cohort <- function(tile_list, patient, method = "mean", ...) {
  stopifnot(length(tile_list) == length(patient))
  S <- t(vapply(tile_list, pool_tiles, numeric(length(pool_tiles(tile_list[[1]], method, ...))), method = method, ...))
  units <- unique(patient)
  P <- rowsum(S, patient, reorder = FALSE)[units, , drop = FALSE] / as.vector(table(patient)[units])
  rownames(P) <- units; P
}

#' Fit a projection on training rows
#'
#' Column standardisation followed by linear PCA or kernel PCA (cosine, RBF with
#' the median pairwise distance as bandwidth, or polynomial of degree 3), keeping
#' the `d` leading components (fewer if the numerical rank is lower), with the
#' scores standardised. Held-out rows are projected with [apply_projection()]:
#' the rotation for PCA and the out-of-sample extension of kernel PCA
#' (Scholkopf, Smola and Muller, 1998), with the training Gram means used for
#' centring.
#'
#' @param X numeric matrix of pooled embeddings (rows = patients).
#' @param decomp one of `"pca"`, `"cosine"`, `"rbf"`, `"poly"`.
#' @param d number of retained coordinates.
#' @return an object of class `pathmi_projection`.
#' @export
fit_projection <- function(X, decomp = c("pca", "cosine", "rbf", "poly"), d = 50L) {
  decomp <- match.arg(decomp); X <- as.matrix(X); tol <- 1e-8
  mu <- colMeans(X); s <- sqrt(colSums(sweep(X, 2, mu)^2) / (nrow(X) - 1)); keep <- is.finite(s) & s > 0
  std <- function(A) sweep(sweep(as.matrix(A)[, keep, drop = FALSE], 2, mu[keep]), 2, s[keep], "/")
  Xs <- std(X)
  if (decomp == "pca") {
    sv <- svd(Xs, nu = 0, nv = min(d, nrow(Xs), ncol(Xs))); lam <- sv$d^2 / (nrow(Xs) - 1)
    k <- min(d, sum(lam > lam[1] * tol), ncol(sv$v)); R <- sv$v[, seq_len(k), drop = FALSE]
    raw <- function(A) std(A) %*% R; spec <- lam[seq_len(k)]
  } else {
    kern <- switch(decomp,
      cosine = { nrm <- function(A) A / sqrt(rowSums(A^2) + 1e-12); A0 <- nrm(Xs); function(B) tcrossprod(nrm(B), A0) },
      rbf = { sq <- rowSums(Xs^2); D2 <- outer(sq, sq, "+") - 2 * tcrossprod(Xs); D2[D2 < 0] <- 0
              sig <- max(sqrt(stats::median(D2[upper.tri(D2)])), 1e-6)
              function(B) { Dd <- outer(rowSums(B^2), sq, "+") - 2 * tcrossprod(B, Xs); Dd[Dd < 0] <- 0; exp(-Dd / (2 * sig^2)) } },
      poly = { pdim <- ncol(Xs); function(B) (tcrossprod(B, Xs) / pdim + 1)^3 })
    K <- kern(Xs); K <- (K + t(K)) / 2; cm <- colMeans(K); am <- mean(K)
    Kc <- sweep(sweep(K, 2, cm), 1, cm) + am; Kc <- (Kc + t(Kc)) / 2
    ev <- eigen(Kc, symmetric = TRUE); lam <- pmax(ev$values, 0); kp <- which(lam > max(lam) * tol)
    k <- min(d, length(kp)); lt <- lam[kp[seq_len(k)]]; W <- sweep(ev$vectors[, kp[seq_len(k)], drop = FALSE], 2, sqrt(lt), "/")
    raw <- function(A) { Kx <- kern(std(A)); Kx <- sweep(Kx, 2, cm) - rowMeans(Kx) + am; Kx %*% W }; spec <- lt
  }
  Str <- raw(X); m <- colMeans(Str); sdv <- apply(Str, 2, stats::sd)
  structure(list(decomp = decomp, k = k, condition = spec[1] / spec[k],
                 apply = function(A) { S <- sweep(sweep(raw(A), 2, m), 2, sdv, "/"); rownames(S) <- rownames(A); colnames(S) <- paste0("Z", seq_len(ncol(S))); S }),
            class = "pathmi_projection")
}

#' @rdname fit_projection
#' @param object a `pathmi_projection` from [fit_projection()].
#' @param newX rows to project (training or held-out).
#' @export
apply_projection <- function(object, newX) object$apply(as.matrix(newX))

#' Project a pooled matrix (fit and apply on the same rows)
#' @inheritParams fit_projection
#' @return standardised scores, one row per patient, with attribute
#'   `condition` (ratio of the largest to the smallest retained eigenvalue).
#' @export
project_embeddings <- function(X, decomp = c("pca", "cosine", "rbf", "poly"), d = 50L) {
  pr <- fit_projection(X, match.arg(decomp), d); Z <- apply_projection(pr, X); attr(Z, "condition") <- pr$condition; Z
}

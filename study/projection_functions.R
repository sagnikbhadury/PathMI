################################################################
# projection_functions.R
#
# The pipeline's projection code (kPCA_SUBHO_*_trident*.R blocks 4-6)
# as functions of the retained dimension, for scripts that rebuild
# projections from the pooled matrices saved by pool_pass.R. Kept
# identical to the pipeline so that d = 50 reproduces the saved
# Z_scaled matrices (checked by every caller).
################################################################
EIG_TOL <- 1e-8

standardize_for_kernel <- function(X) {
  Xs <- scale(X, center = TRUE, scale = TRUE)
  keep <- colSums(!is.finite(Xs)) == 0
  Xs <- Xs[, keep, drop = FALSE]; rownames(Xs) <- rownames(X); Xs
}
kernel_cosine <- function(X) { Xn <- X / sqrt(rowSums(X^2) + 1e-12); tcrossprod(Xn) }
kernel_rbf <- function(X) {
  D2 <- as.matrix(dist(X, method = "euclidean"))^2
  sigma <- max(sqrt(median(D2[upper.tri(D2)])), 1e-6)
  exp(-D2 / (2 * sigma^2))
}
kernel_poly <- function(X, degree = 3L, coef0 = 1.0) { (tcrossprod(X) / ncol(X) + coef0)^degree }
centre_gram <- function(K) { n <- nrow(K); H <- diag(n) - matrix(1 / n, n, n); Kc <- H %*% K %*% H; (Kc + t(Kc)) / 2 }

# linear PCA, as apply_pca(): prcomp(center, scale) -> top d scores -> standardise
project_pca <- function(Z_raw, d) {
  pca <- prcomp(Z_raw, center = TRUE, scale. = TRUE)
  k <- min(d, ncol(pca$x)); lam <- pca$sdev^2
  list(Z_scaled = scale(pca$x[, 1:k, drop = FALSE], center = TRUE, scale = TRUE), k_use = k,
       eig_cond = lam[1] / lam[k], var_frac = sum(lam[1:k]) / sum(lam))
}
# kernel PCA, as apply_kpca(): standardise -> kernel -> centre -> eigen -> top d scaled by sqrt(lambda) -> standardise
project_kpca <- function(Z_raw, kernel = c("cosine", "rbf", "poly"), d) {
  kernel <- match.arg(kernel)
  K_c <- centre_gram(switch(kernel, cosine = kernel_cosine, rbf = kernel_rbf, poly = kernel_poly)(standardize_for_kernel(Z_raw)))
  ev <- eigen(K_c, symmetric = TRUE); lam <- pmax(ev$values, 0); keep <- which(lam > max(lam) * EIG_TOL)
  k <- min(d, length(keep)); lt <- lam[keep[seq_len(k)]]; vt <- ev$vectors[, keep[seq_len(k)], drop = FALSE]
  s <- sweep(vt, 2, sqrt(lt), `*`); rownames(s) <- rownames(Z_raw)
  list(Z_scaled = scale(s, center = TRUE, scale = TRUE), k_use = k, eig_cond = max(lt) / min(lt), var_frac = sum(lt) / (sum(lam) + 1e-12))
}
project <- function(Z_raw, decomp, d) if (decomp == "pca") project_pca(Z_raw, d) else project_kpca(Z_raw, decomp, d)

# max abs difference after aligning column signs (eigenvectors are defined up to sign)
align_diff <- function(A, B) { s <- sign(colSums(A * B)); s[s == 0] <- 1; max(abs(sweep(A, 2, s, `*`) - B)) }

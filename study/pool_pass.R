################################################################
# pool_pass.R  <crc|prad>  <uni_v2|virchow2|hoptimus1|conch_v15>
#
# One pass over the Trident h5 files of a cohort under one model,
# for the roster the pipeline used (row names of the saved Z_scaled):
#
#   1. the four pooled patient-level matrices (mean, norm-weighted,
#      moment = mean|sd, untrained attention), computed with the same
#      numpy float32 code as the pipeline's load_h5_pooled(), so that
#      projections at other retained dimensions (robustness_dimension.R)
#      start from exactly the matrices the 50-dimensional runs used;
#   2. for mean pooling, the mean of m randomly drawn tiles per slide
#      for m on a 16-point grid and R = 20 draws (float32), which gives
#      the tissue-sampling burden of every model under its linear
#      readout (sampling_burden_all_models.R) without storing tiles:
#      under mean pooling with linear PCA the slide logit from m tiles
#      is the affine map applied to the mean of those m tiles.
#
# Checks: the mean matrix -> prcomp -> scale must reproduce the saved
# pca_mean Z_scaled; the mean matrix -> cosine kernel PCA (the
# pipeline's kernel code, copied below) must reproduce the saved
# cosine_mean Z_scaled, both up to column sign.
#
# Output (large, on the results share, not in Dropbox):
#   M:/bhadury/mi_results/pool_pass/<cohort>_<model>/pooled.rds
#   M:/bhadury/mi_results/pool_pass/<cohort>_<model>/submeans.npy  (n x 16 x 20 x D float32)
#   TRIDENT/pool_pass/<cohort>_<model>.log and manifest row
################################################################
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
COHORT <- tolower(args[1]); MODEL <- args[2]
stopifnot(COHORT %in% c("crc", "prad"), MODEL %in% c("uni_v2", "virchow2", "hoptimus1", "conch_v15"))

suppressPackageStartupMessages({ library(reticulate) })
reticulate::py_require(c("h5py", "numpy"))
py_run_string("import h5py, numpy as np")

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI  <- normalizePath(dirname(this_file), winslash = "/")
LOGDIR <- file.path(TRI, "pool_pass"); dir.create(LOGDIR, showWarnings = FALSE)
OUT <- sprintf("M:/bhadury/mi_results/pool_pass/%s_%s", COHORT, MODEL); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(LOGDIR, sprintf("%s_%s.log", COHORT, MODEL))
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }

RANDOM_SEED <- 42L; ATTN_HIDDEN <- 128L; N_PC <- 50L; EIG_TOL <- 1e-8
M_GRID <- c(5L, 10L, 15L, 20L, 25L, 30L, 40L, 50L, 75L, 100L, 150L, 200L, 300L, 500L, 1000L, 2500L); R_SUB <- 20L

H5 <- list(
  crc = list(uni_v2 = "M:/kulkarnr/slideflow_test/slideflow_crc/trident_feats/20x_256px_0px_overlap/features_uni_v2/",
             virchow2 = "M:/kulkarnr/slideflow_test/slideflow_crc/trident_feats/20x_224px_0px_overlap/features_virchow2/",
             hoptimus1 = "M:/kulkarnr/slideflow_test/slideflow_crc/trident_feats/20x_224px_0px_overlap/features_hoptimus1/",
             conch_v15 = "M:/kulkarnr/slideflow_test/slideflow_crc/trident_feats/20x_512px_0px_overlap/features_conch_v15/"),
  prad = list(uni_v2 = "M:/shared_data/TCGA_PRAD/trident-feats/20x_256px_0px_overlap/features_uni_v2/",
              virchow2 = "M:/shared_data/TCGA_PRAD/trident-feats/20x_224px_0px_overlap/features_virchow2/",
              hoptimus1 = "M:/shared_data/TCGA_PRAD/trident-feats/20x_224px_0px_overlap/features_hoptimus1/",
              conch_v15 = "M:/shared_data/TCGA_PRAD/trident-feats/20x_512px_0px_overlap/features_conch_v15/"))
h5dir <- H5[[COHORT]][[MODEL]]
rds <- file.path(TRI, if (COHORT == "crc") sprintf("crc_kpca_subho_%s/kpca_subho_results_crc_%s.rds", MODEL, MODEL) else sprintf("TCGA_PRAD_kpca_subho_%s/kpca_subho_results_prad_%s.rds", MODEL, MODEL))

t0 <- Sys.time()
logf("=== pool pass | %s | %s ===", toupper(COHORT), MODEL)
res <- readRDS(rds)
Zs_mean_saved <- res$pca_mean$Z_scaled; Zs_cos_saved <- res$cosine_mean$Z_scaled
roster <- rownames(Zs_mean_saved); n <- length(roster)
ylist <- if (COHORT == "crc") list(MSI = res$pca_mean$y[roster]) else lapply(res$pca_mean$y_list, function(v) v[roster])
logf("roster n = %d | targets: %s", n, paste(names(ylist), collapse = ", "))

files <- list.files(h5dir, pattern = "[.]h5$")
stem_map <- list()
for (f in files) {
  py_run_string(sprintf("
with h5py.File(r'%s', 'r') as _f:
    _name = _f['features'].attrs.get('name', None)
", file.path(h5dir, f)))
  nm <- py$`_name`; if (!is.null(nm)) stem_map[[nm]] <- f
}
stopifnot(all(roster %in% names(stem_map)))
logf("indexed %d h5 files; all roster slides present", length(files))

# ---- python state: pooling code identical to the pipeline, plus subsample means ----
py_run_string(sprintf("
M_GRID = np.array([%s], dtype=np.int64); R_SUB = %d; N = %d
SUB = None
_rng_attn = None
", paste(M_GRID, collapse = ","), R_SUB, n))

pool_one <- function(i, fpath) {
  py_run_string(sprintf("
with h5py.File(r'%s', 'r') as _f:
    _emb = _f['features'][:].astype(np.float32)
    _a = dict(_f['coords'].attrs)
_n = int(_emb.shape[0]); _d = int(_emb.shape[1])
_ps = int(_a.get('patch_size_level0', 0)); _lw = int(_a.get('level0_width', 0)); _lh = int(_a.get('level0_height', 0))
# mean
_mean = _emb.mean(axis=0)
# norm-weighted
_norms = np.linalg.norm(_emb, axis=1); _w = _norms / _norms.sum(); _normp = (_w @ _emb)
# moment
_sd = _emb.std(axis=0, ddof=1); _moment = np.concatenate([_mean, _sd])
# attention (Kaiming uniform, seed 42, never trained; identical to the pipeline)
_rng  = np.random.default_rng(%d)
_bW   = 1.0 / np.sqrt(_d); _bv = 1.0 / np.sqrt(%d)
_W    = _rng.uniform(-_bW, _bW, (_d, %d)); _v = _rng.uniform(-_bv, _bv, (%d, 1))
_H    = np.tanh(_emb @ _W); _at = _H @ _v; _at = np.exp(_at - _at.max()); _at = _at / _at.sum()
_attn = (_at.T @ _emb).ravel()
# subsample means for the sampling burden
if SUB is None:
    SUB = np.empty((N, len(M_GRID), R_SUB, _d), dtype=np.float32)
_rs = np.random.default_rng([2026, %d])
for _mi, _m in enumerate(M_GRID):
    if _m >= _n:
        SUB[%d, _mi, :, :] = _mean
    else:
        for _r in range(R_SUB):
            _idx = _rs.choice(_n, int(_m), replace=False)
            SUB[%d, _mi, _r, :] = _emb[_idx].mean(axis=0)
del _emb, _H, _at, _W, _v, _norms, _w
", fpath, RANDOM_SEED, ATTN_HIDDEN, ATTN_HIDDEN, ATTN_HIDDEN, i, i - 1L, i - 1L))
  list(mean = as.numeric(py$`_mean`), norm = as.numeric(py$`_normp`), moment = as.numeric(py$`_moment`), attention = as.numeric(py$`_attn`),
       n_tiles = as.integer(py$`_n`), patch_level0 = as.integer(py$`_ps`), level0_w = as.integer(py$`_lw`), level0_h = as.integer(py$`_lh`))
}

pooled <- list(mean = NULL, norm = NULL, moment = NULL, attention = NULL)
tiles <- data.frame(slide = roster, n_tiles = NA_integer_, patch_level0 = NA_integer_, level0_w = NA_integer_, level0_h = NA_integer_, stringsAsFactors = FALSE)
for (i in seq_len(n)) {
  o <- pool_one(i, file.path(h5dir, stem_map[[roster[i]]]))
  if (i == 1L) for (k in names(pooled)) pooled[[k]] <- matrix(NA_real_, n, length(o[[k]]), dimnames = list(roster, NULL))
  for (k in names(pooled)) pooled[[k]][i, ] <- o[[k]]
  tiles$n_tiles[i] <- o$n_tiles; tiles$patch_level0[i] <- o$patch_level0; tiles$level0_w[i] <- o$level0_w; tiles$level0_h[i] <- o$level0_h
  if (i %% 50 == 0) logf("  %d / %d slides (%.1f min)", i, n, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
logf("pass done: D = %d | tiles per slide median %d (min %d, max %d) | patch at level 0: %s px",
     ncol(pooled$mean), as.integer(median(tiles$n_tiles)), min(tiles$n_tiles), max(tiles$n_tiles), paste(unique(tiles$patch_level0), collapse = "/"))

# ---- checks against the saved projections --------------------------------------
align <- function(A, B) { s <- sign(colSums(A * B)); s[s == 0] <- 1; max(abs(sweep(A, 2, s, `*`) - B)) }
pca <- prcomp(pooled$mean, center = TRUE, scale. = TRUE)
Zs_rec <- scale(pca$x[, 1:N_PC], center = TRUE, scale = TRUE)
d_pca <- align(Zs_rec, Zs_mean_saved)
logf("check 1: linear PCA (mean pooling) rebuilt vs saved Z_scaled: max abs diff %.2e", d_pca)

standardize_for_kernel <- function(X) { Xs <- scale(X, center = TRUE, scale = TRUE); keep <- colSums(!is.finite(Xs)) == 0; Xs <- Xs[, keep, drop = FALSE]; rownames(Xs) <- rownames(X); Xs }
kernel_cosine <- function(X) { Xn <- X / sqrt(rowSums(X^2) + 1e-12); tcrossprod(Xn) }
centre_gram <- function(K) { n <- nrow(K); H <- diag(n) - matrix(1 / n, n, n); Kc <- H %*% K %*% H; (Kc + t(Kc)) / 2 }
kpca_scores <- function(Z_raw, kernel_fn, n_pc = N_PC) {
  K_c <- centre_gram(kernel_fn(standardize_for_kernel(Z_raw)))
  ev <- eigen(K_c, symmetric = TRUE); lam <- pmax(ev$values, 0); keep <- which(lam > max(lam) * EIG_TOL)
  k_use <- min(n_pc, length(keep)); lt <- lam[keep[seq_len(k_use)]]; vt <- ev$vectors[, keep[seq_len(k_use)], drop = FALSE]
  s <- sweep(vt, 2, sqrt(lt), `*`); rownames(s) <- rownames(Z_raw); list(Z_scaled = scale(s, center = TRUE, scale = TRUE), eig_cond = max(lt) / min(lt))
}
kc <- kpca_scores(pooled$mean, kernel_cosine)
d_cos <- align(kc$Z_scaled, Zs_cos_saved)
logf("check 2: cosine kernel PCA (mean pooling) rebuilt vs saved Z_scaled: max abs diff %.2e (kappa %.1f)", d_cos, kc$eig_cond)
if (d_pca > 1e-5 || d_cos > 1e-3) logf("[WARN] a rebuilt projection does not match the saved one")

# ---- save -----------------------------------------------------------------------
saveRDS(list(cohort = COHORT, model = MODEL, roster = roster, y = ylist, tiles = tiles, pooled = pooled,
             m_grid = M_GRID, r_sub = R_SUB, check_pca = d_pca, check_cosine = d_cos),
        file.path(OUT, "pooled.rds"))
py_run_string(sprintf("np.save(r'%s', SUB)", file.path(OUT, "submeans.npy")))
manifest <- data.frame(cohort = COHORT, model = MODEL, n = n, D = ncol(pooled$mean), tiles_median = median(tiles$n_tiles),
                       patch_level0 = paste(unique(tiles$patch_level0), collapse = "/"), check_pca = d_pca, check_cosine = d_cos,
                       minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), out = OUT)
mf <- file.path(LOGDIR, "manifest.csv")
write.table(manifest, mf, sep = ",", row.names = FALSE, col.names = !file.exists(mf), append = file.exists(mf))
logf("saved pooled.rds and submeans.npy (%s) in %.1f min", OUT, manifest$minutes)

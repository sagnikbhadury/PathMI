################################################################
# kPCA_SUBHO_CRC_trident_hoptimus_conch.R
# Kernel PCA — SUBHO Logistic MI Pipeline | CRC Trident embeddings
#
# Sibling of kPCA_SUBHO_PRAD_trident.R, collapsed to a single
# target (MSI status) since CRC has no second target like PRAD's
# ERG/ETS. Embeddings come from Trident, saved as one .h5 per
# slide, same format already confirmed working for PRAD: two
# datasets, 'features' (N_tiles x D) and 'coords' (N_tiles x 2),
# with the canonical slide name stored as features.attrs['name'].
#
# EMBED_TYPE switches between "uni_v2" and "virchow2" so the same
# script drives both runs. Set EMBED_TYPE below and source the
# whole file once per run. Output dirs/file tags are suffixed
# with EMBED_TYPE so the two runs never overwrite each other.
#
# EMBED_DIM is set per the PRAD Trident run (uni_v2 = 1536,
# virchow2 = 2560). discover_h5_keys() will warn if CRC's actual
# feature dimension differs — confirm on first run rather than
# assuming CRC's Trident extraction used identical settings.
#
# Statistical fixes ported from kPCA_SUBHO_PRAD_trident.R (same
# external review, R1-R5), with R2/R3 simplified for the
# single-target case — CRC's original script never had PRAD's
# two-target cross-leakage problem (best_key chosen by MI pooled
# across ERG/ETS), so there is no per-target loop here. R2 still
# matters: the original CRC select_best_decomp() filtered
# candidates by observed p_value < 0.05 before building the
# nested-permutation null, which is the same outcome-adaptive
# selection bug as PRAD's any_sig filter, just for one target
# instead of two.
#
#   R1. Pooled embeddings are column-standardized before all
#       three kernels (cosine, RBF, poly), not just linear PCA.
#   R2. select_best_decomp() no longer filters candidate keys by
#       observed p_value before building the null maximum. maxT
#       null is computed over the full pooling-group candidate
#       set unconditionally.
#   R3. make_figures() takes a single result and its true pool
#       name, rather than the original CRC script's "pseudo-pair"
#       pattern (passing one pooling method's result into both an
#       out_mean and out_norm slot, with panels hardcoded
#       "Mean pooling" / "Norm-weighted pooling" regardless of
#       which method actually produced the data). Every call site
#       here only ever has one pooling method's result at a time,
#       so the figure now shows exactly what was computed, with
#       no faked second panel.
#   R4. compute_incremental_mi() passes method.weight =
#       "beta.binomial" to all three subho_MI_logistic() calls,
#       matching run_subho(). Output named incremental_subho_mi,
#       documented as a model-based recoverability delta, not
#       asserted to equal population conditional MI I(Y;Z|C).
#   R5. log_tile_count is standardized in build_covariate_matrix()
#       alongside stage_ordinal_z (and site_READ, if retained), so
#       all covariates carry the same effective prior weight as
#       the unit-scaled embedding PCs.
#
# Dependencies:
#   install.packages(c("MASS","BayesLogit","matrixcalc",
#                      "reticulate","ggplot2","gridExtra",
#                      "scales","dplyr"))
#   Python env (via reticulate) needs h5py and numpy.
################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(MASS)
  library(BayesLogit)
  library(matrixcalc)
  library(reticulate)
  library(ggplot2)
  library(gridExtra)
  library(scales)
  library(dplyr)
})

reticulate::py_require(c("h5py", "numpy"))
np <- reticulate::import("numpy")

source("PathMI_SUBHO_logistic.R")

################################################################
# EMBED_TYPE switch — new backbones (2026-09-28)
#
# This file is a byte-for-byte sibling of kPCA_SUBHO_CRC_trident.R
# (the UNI v2 / Virchow2 script, left untouched) with only this
# block changed, so the two new backbones run through exactly the
# same pipeline, settings and seeds as the two existing ones.
#
# Set EMBED_TYPE below and source once per run, or set the
# environment variable PATHMI_EMBED_TYPE before launching (used by
# run_backbone.R so both backbones run as separate Rscript
# processes without editing this line).
#
# Both are Trident h5 (keys 'coords' + 'features', slide name in
# features.attrs['name']); dims confirmed on disk 2026-09-28:
#   hoptimus1  1536-d   20x_224px tiles   (H-optimus-1)
#   conch_v15   768-d   20x_512px tiles   (CONCH v1.5)
################################################################
EMBED_TYPE <- Sys.getenv("PATHMI_EMBED_TYPE", unset = "hoptimus1")   # "hoptimus1" | "conch_v15"

stopifnot(EMBED_TYPE %in% c("hoptimus1", "conch_v15"))

EMBED_DIM <- switch(EMBED_TYPE, hoptimus1 = 1536L, conch_v15 = 768L)

EMBED_DIR <- switch(EMBED_TYPE,
  hoptimus1 = "M:/kulkarnr/slideflow_test/slideflow_crc/trident_feats/20x_224px_0px_overlap/features_hoptimus1/",
  conch_v15 = "M:/kulkarnr/slideflow_test/slideflow_crc/trident_feats/20x_512px_0px_overlap/features_conch_v15/"
)

META_PATH   <- "M:/kulkarnr/slideflow_test/slideflow_crc/annotations.csv"
RESULTS_DIR <- sprintf("M:/bhadury/mi_results/crc_kpca_subho_%s/", EMBED_TYPE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("EMBED_TYPE = %s | EMBED_DIM = %d | EMBED_DIR = %s\n",
            EMBED_TYPE, EMBED_DIM, EMBED_DIR))

################################################################
# CONFIG — identical values to CRC original / PRAD Trident
################################################################
N_PC        <- 50L
N_PERM      <- 200L
RANDOM_SEED <- 42L
MIN_TILES   <- 100L
V_BETA      <- 1.0
NMC         <- 2000L
BURN        <- 500L
THIN        <- 10L
EIG_TOL     <- 1e-8
ATTN_HIDDEN <- 128L

COL_SLIDE <- "slide"
COL_MSI   <- "MSIStatus"; MSI_POS <- "MSI-H"

PALETTE   <- list(mean      = "#2166AC",
                  norm      = "#D6604D",
                  moment    = "#1B7837",
                  attention = "#762A83",
                  null      = "#BDBDBD")
AGG_LABEL <- list(mean      = "Mean pooling",
                  norm      = "Norm-weighted pooling",
                  moment    = "Moment pooling (mean+SD)",
                  attention = "Attention pooling (Kaiming)")

np <- reticulate::import("numpy")

################################################################
# BLOCK 1: Index embedding directory (.h5 files)
#
# Keyed by the h5 file's internal features['name'] attribute,
# NOT by stripping the .h5 suffix from the on-disk filename.
# Confirmed on the PRAD Trident run that on-disk filename
# punctuation can drift from the metadata 'slide' convention;
# reading the attribute directly removes that dependence
# entirely, since both the index and the metadata CSV read the
# same Trident-assigned canonical name.
#
# Cost: opens every h5 once to read an attribute. Attribute
# reads do not touch the features array, so this is cheap even
# for thousands of files (no full dataset read happens here).
################################################################

index_embed_dir <- function(embed_dir) {
  files <- list.files(embed_dir, pattern = "\\.h5$", full.names = FALSE)
  if (length(files) == 0L)
    stop(sprintf("No .h5 files found in %s", embed_dir))

  stem_map     <- list()
  name_missing <- character(0)
  name_dup     <- character(0)

  for (f in files) {
    fpath <- file.path(embed_dir, f)
    nm <- tryCatch({
      reticulate::py_run_string(sprintf("
import h5py
with h5py.File(r'%s', 'r') as _f:
    _name = _f['features'].attrs.get('name', None)
", fpath))
      val <- reticulate::py$`_name`
      reticulate::py_run_string("del _name")
      val
    }, error = function(e) {
      warning(sprintf("  [ATTR READ ERROR] %s: %s", f, e$message))
      NULL
    })

    if (is.null(nm)) {
      name_missing <- c(name_missing, f)
      next
    }
    if (!is.null(stem_map[[nm]])) {
      name_dup <- c(name_dup, nm)
    }
    stem_map[[nm]] <- f
  }

  cat(sprintf("Indexed %d h5 files -> %d unique slide names (by attrs['name']).\n",
              length(files), length(stem_map)))
  if (length(name_missing) > 0)
    cat(sprintf("  [WARN] %d file(s) had no readable 'name' attribute: %s%s\n",
                length(name_missing),
                paste(head(name_missing, 5), collapse = ", "),
                ifelse(length(name_missing) > 5, ", ...", "")))
  if (length(name_dup) > 0)
    cat(sprintf("  [WARN] %d duplicate slide name(s) across files: %s\n",
                length(name_dup), paste(unique(name_dup), collapse = ", ")))

  stem_map
}

################################################################
# BLOCK 1c: Diagnose roster-vs-directory mismatch
#
# Run once after indexing, before the full roster build. Compares
# metadata 'slide' values against the stem_map keys directly so a
# naming-convention mismatch shows up immediately as a string
# diff, rather than being inferred indirectly from a low patient
# count after the expensive roster-building loop.
################################################################

diagnose_slide_matching <- function(meta_path, stem_map) {
  meta   <- read.csv(meta_path, stringsAsFactors = FALSE)
  slides <- meta[[COL_SLIDE]]
  slides <- slides[!is.na(slides)]

  in_dir  <- slides %in% names(stem_map)
  n_match <- sum(in_dir)
  n_total <- length(slides)

  cat(sprintf("\n[diagnose_slide_matching] %d / %d metadata slides found in directory index.\n",
              n_match, n_total))

  if (n_match < n_total) {
    missing <- slides[!in_dir]
    cat(sprintf("  %d metadata slide(s) NOT found. First few:\n", length(missing)))
    for (s in head(missing, 5)) cat(sprintf("    metadata: %s\n", s))

    cat("  First few directory-indexed names for comparison:\n")
    for (s in head(names(stem_map), 5)) cat(sprintf("    directory: %s\n", s))

    case_fold_match <- sum(tolower(missing) %in% tolower(names(stem_map)))
    if (case_fold_match > 0)
      cat(sprintf("  [HINT] %d of the missing slide(s) match if case is ignored — likely a casing mismatch.\n",
                  case_fold_match))
  }

  invisible(in_dir)
}

################################################################
# BLOCK 1b: One-time key discovery / sanity check
#
# Confirmed for PRAD Trident: keys are 'coords' and 'features'.
# Re-verify here rather than assume CRC's Trident extraction used
# the same key names — Trident version drift between extraction
# runs is plausible if CRC and PRAD were processed at different
# times. 'coords' is not required for pooling, only flagged.
################################################################

discover_h5_keys <- function(embed_dir, stem_map) {
  if (length(stem_map) == 0L)
    stop(sprintf("No .h5 files found in %s", embed_dir))

  probe_file <- file.path(embed_dir, stem_map[[1]])

  reticulate::py_run_string(sprintf("
import h5py
with h5py.File(r'%s', 'r') as _f:
    _keys  = list(_f.keys())
    _shape = _f['features'].shape if 'features' in _f else None
    _dtype = str(_f['features'].dtype) if 'features' in _f else None
", probe_file))

  keys  <- reticulate::py$`_keys`
  shape <- reticulate::py$`_shape`
  dtype <- reticulate::py$`_dtype`
  reticulate::py_run_string("del _keys, _shape, _dtype")

  cat(sprintf("Probed: %s\n", probe_file))
  cat(sprintf("  Keys found: %s\n", paste(keys, collapse = ", ")))

  if (!("features" %in% keys))
    stop(sprintf(
      "Expected key 'features' not found. Actual keys: %s. Update load_h5_pooled() accordingly.",
      paste(keys, collapse = ", ")))
  if (!("coords" %in% keys))
    cat("  [WARN] 'coords' key not found — not required for pooling, only flagging.\n")

  if (!is.null(shape))
    cat(sprintf("  features shape: %s | dtype: %s\n",
                paste(shape, collapse = " x "), dtype))

  if (!is.null(shape) && length(shape) == 2 && shape[[2]] != EMBED_DIM)
    cat(sprintf(
      "  [WARN] features dim (%d) != expected EMBED_DIM (%d) for EMBED_TYPE='%s'. Check EMBED_DIM.\n",
      shape[[2]], EMBED_DIM, EMBED_TYPE))

  invisible(keys)
}

################################################################
# BLOCK 2: Load and pool a single slide (.h5)
#
# Trident h5: 'features' dataset is (N_tiles x EMBED_DIM) float.
# Read directly with h5py. Pooling logic identical to PRAD Trident
# and to the original CRC .npy loader.
################################################################

load_h5_pooled <- function(embed_dir, slide_name, stem_map,
                           agg_method = "mean") {
  fname <- stem_map[[slide_name]]
  if (is.null(fname)) return(list(vec = NULL, n_tiles = NA_integer_))

  fpath <- file.path(embed_dir, fname)

  n_tiles <- tryCatch({
    reticulate::py_run_string(sprintf("
import h5py
import numpy as np
with h5py.File(r'%s', 'r') as _f:
    _emb = _f['features'][:].astype(np.float32)
_n = int(_emb.shape[0])
", fpath))
    as.integer(reticulate::py$`_n`)
  }, error = function(e) {
    warning(sprintf("  [LOAD ERROR] %s: %s", slide_name, e$message))
    NA_integer_
  })

  if (is.na(n_tiles)) return(list(vec = NULL, n_tiles = NA_integer_))

  if (n_tiles < MIN_TILES) {
    reticulate::py_run_string("del _emb, _n")
    return(list(vec = NULL, n_tiles = n_tiles))
  }

  pool_code <- switch(agg_method,
                      mean = "
_pooled = _emb.mean(axis=0)
del _emb
",
                      norm = "
_norms  = np.linalg.norm(_emb, axis=1)
_w      = _norms / _norms.sum()
_pooled = (_w @ _emb)
del _emb, _norms, _w
",
                      moment = "
_mu     = _emb.mean(axis=0)
_sd     = _emb.std(axis=0, ddof=1)
_pooled = np.concatenate([_mu, _sd])
del _emb, _mu, _sd
",
                      attention = sprintf(
                        "_rng  = np.random.default_rng(%d)
_d    = _emb.shape[1]
_bW   = 1.0 / np.sqrt(_d)
_bv   = 1.0 / np.sqrt(%d)
_W    = _rng.uniform(-_bW, _bW, (_d, %d))
_v    = _rng.uniform(-_bv, _bv, (%d, 1))
_H    = np.tanh(_emb @ _W)
_a    = _H @ _v
_a    = np.exp(_a - _a.max())
_a    = _a / _a.sum()
_pooled = (_a.T @ _emb).ravel()
del _emb, _W, _v, _H, _a, _bW, _bv",
                        as.integer(RANDOM_SEED),
                        as.integer(ATTN_HIDDEN),
                        as.integer(ATTN_HIDDEN),
                        as.integer(ATTN_HIDDEN)),
                      stop(sprintf("Unknown agg_method: '%s'", agg_method))
  )

  vec <- tryCatch({
    reticulate::py_run_string(pool_code)
    v <- as.numeric(reticulate::py$`_pooled`)
    reticulate::py_run_string("del _pooled, _n")
    v
  }, error = function(e) {
    tryCatch(reticulate::py_run_string(
      "for _v in ['_emb','_pooled','_n','_mu','_sd',
                  '_norms','_w','_W','_H','_a','_rng']:
  globals().pop(_v, None)"),
      error = function(e2) NULL)
    warning(sprintf("  [POOL ERROR] %s [%s]: %s",
                    slide_name, agg_method, e$message))
    NULL
  })

  list(vec = vec, n_tiles = n_tiles)
}

################################################################
# BLOCK 3: Build unified patient roster
#
# Single target (MSI), so this is a direct h5 port of the
# original CRC .npy roster builder: one label_vec, not a
# per-target label_list as in PRAD.
################################################################

build_unified_roster <- function(meta_path, embed_dir, stem_map) {
  meta            <- read.csv(meta_path, stringsAsFactors = FALSE)
  meta$msi_binary <- as.integer(meta[[COL_MSI]] == MSI_POS)

  cat(sprintf("Metadata: %d rows | MSI-H: %d | MSS+MSI-L: %d | NA: %d\n",
              nrow(meta), sum(meta$msi_binary == 1L, na.rm = TRUE),
              sum(meta$msi_binary == 0L, na.rm = TRUE),
              sum(is.na(meta$msi_binary))))

  emb_mean      <- list()
  emb_norm      <- list()
  emb_moment    <- list()
  emb_attention <- list()
  id_accum      <- character(0)
  label_accum   <- integer(0)
  tile_counts   <- c()

  agg_methods <- c("mean", "norm", "moment", "attention")

  for (i in seq_len(nrow(meta))) {
    slide_id <- meta[[COL_SLIDE]][i]
    label    <- meta$msi_binary[i]

    if (is.na(slide_id)) {
      cat(sprintf("  [SKIP - NA slide] row %d\n", i)); next
    }
    if (is.na(label)) {
      cat(sprintf("  [SKIP - missing label] %s\n", slide_id)); next
    }

    results <- lapply(agg_methods, function(m)
      load_h5_pooled(embed_dir, slide_id, stem_map, m))
    names(results) <- agg_methods

    not_found <- any(sapply(results, function(r) is.na(r$n_tiles)))
    if (not_found) {
      cat(sprintf("  [SKIP - not found] %s\n", slide_id)); next
    }

    any_qc_fail <- any(sapply(results, function(r) is.null(r$vec)))
    if (any_qc_fail) {
      tile_info <- paste(sapply(agg_methods, function(m)
        sprintf("%s=%s", m,
                ifelse(is.null(results[[m]]$vec),
                       sprintf("FAIL(%d)", results[[m]]$n_tiles),
                       "OK"))),
        collapse = " | ")
      cat(sprintf("  [QC FAIL] %s - %s\n", slide_id, tile_info))
      next
    }

    emb_mean[[slide_id]]      <- results$mean$vec
    emb_norm[[slide_id]]      <- results$norm$vec
    emb_moment[[slide_id]]    <- results$moment$vec
    emb_attention[[slide_id]] <- results$attention$vec
    tile_counts[slide_id]     <- results$mean$n_tiles
    id_accum                  <- c(id_accum, slide_id)
    label_accum                <- c(label_accum, as.integer(label))
  }

  # Names attached once at the end via id_accum, same pattern as
  # the PRAD fix for label vectors silently losing names.
  label_vec <- label_accum
  names(label_vec) <- id_accum

  n_ret <- length(tile_counts)
  cat(sprintf("Unified roster: %d patients | MSI-H: %d | MSS+MSI-L: %d\n",
              n_ret, sum(label_vec), sum(label_vec == 0)))
  cat(sprintf("Tile counts — min: %d | median: %.0f | max: %d\n",
              min(tile_counts), median(tile_counts), max(tile_counts)))
  cat(sprintf("Embedding dim — mean/norm/attn: %d | moment: %d\n",
              length(emb_mean[[1]]), length(emb_moment[[1]])))

  list(emb_mean      = emb_mean,
       emb_norm      = emb_norm,
       emb_moment    = emb_moment,
       emb_attention = emb_attention,
       label_vec     = label_vec,
       tile_counts   = tile_counts)
}

################################################################
# BLOCK 4: Kernel constructors
#
# R1 fix: pooled embeddings are column-standardized before all
# three kernels, not just linear PCA. Without this, RBF and poly
# are dominated by raw embedding magnitude — a real concern when
# comparing UNI v2 (1536-dim) against Virchow2 (2560-dim), which
# sit at different native scales. Cosine is recomputed on the
# same standardized matrix for consistency: this measures angle
# w.r.t. centred, unit-variance feature axes, a different
# geometric object from raw-magnitude cosine.
################################################################

standardize_for_kernel <- function(X) {
  Xs        <- scale(X, center = TRUE, scale = TRUE)
  keep_cols <- colSums(!is.finite(Xs)) == 0
  if (any(!keep_cols))
    cat(sprintf("    [standardize_for_kernel] dropped %d zero-variance column(s)\n",
                sum(!keep_cols)))
  Xs <- Xs[, keep_cols, drop = FALSE]
  rownames(Xs) <- rownames(X)
  Xs
}

kernel_cosine <- function(X) {
  X_norm <- X / sqrt(rowSums(X^2) + 1e-12)
  tcrossprod(X_norm)
}

kernel_rbf <- function(X) {
  D2    <- as.matrix(dist(X, method = "euclidean"))^2
  sigma <- sqrt(median(D2[upper.tri(D2)]))
  sigma <- max(sigma, 1e-6)
  cat(sprintf("    RBF sigma (median heuristic): %.4f\n", sigma))
  exp(-D2 / (2 * sigma^2))
}

kernel_poly <- function(X, degree = 3L, coef0 = 1.0) {
  d <- ncol(X)
  (tcrossprod(X) / d + coef0)^degree
}

centre_gram <- function(K) {
  n   <- nrow(K)
  H   <- diag(n) - matrix(1 / n, n, n)
  K_c <- H %*% K %*% H
  (K_c + t(K_c)) / 2
}

################################################################
# BLOCK 5: Standard linear PCA
################################################################

apply_pca <- function(Z_raw, n_pc = N_PC, agg_method = "unknown") {
  cat(sprintf("\n[Standard PCA — %s pooling]\n", agg_method))
  cat(sprintf("  Input: %d x %d\n", nrow(Z_raw), ncol(Z_raw)))

  pca   <- prcomp(Z_raw, center = TRUE, scale. = TRUE)
  k_use <- min(n_pc, ncol(pca$x))

  Z_pca    <- pca$x[, 1:k_use, drop = FALSE]
  Z_scaled <- scale(Z_pca, center = TRUE, scale = TRUE)

  lam       <- pca$sdev^2
  var_proxy <- sum(lam[1:k_use]) / sum(lam)
  eig_cond  <- lam[1] / lam[k_use]

  cat(sprintf("  Components retained: %d\n", k_use))
  cat(sprintf("  Variance explained (top %d / full spectrum): %.1f%%\n",
              k_use, 100 * var_proxy))
  cat(sprintf("  Eigenvalue condition (lambda_1/lambda_k): %.2e\n", eig_cond))

  list(Z_scaled    = Z_scaled,
       eigenvalues = lam,
       k_use       = k_use,
       var_proxy   = var_proxy,
       eig_cond    = eig_cond,
       agg_method  = agg_method)
}

################################################################
# BLOCK 6: Kernel PCA — standardization moved upstream into
# apply_kpca() itself (R1); eigendecomposition/scoring unchanged.
################################################################

apply_kpca <- function(Z_raw, kernel = c("cosine", "rbf", "poly"),
                       n_pc = N_PC, agg_method = "unknown",
                       eig_tol = EIG_TOL) {
  kernel <- match.arg(kernel)
  cat(sprintf("\n[kernel PCA (%s) — %s pooling]\n", kernel, agg_method))
  cat(sprintf("  Input: %d x %d\n", nrow(Z_raw), ncol(Z_raw)))

  Z_std <- standardize_for_kernel(Z_raw)

  K_raw <- switch(kernel,
                  cosine = kernel_cosine(Z_std),
                  rbf    = kernel_rbf(Z_std),
                  poly   = kernel_poly(Z_std)
  )
  K_c <- centre_gram(K_raw)

  ev   <- eigen(K_c, symmetric = TRUE)
  lam  <- pmax(ev$values, 0)
  vecs <- ev$vectors

  tol  <- max(lam) * eig_tol
  keep <- which(lam > tol)
  if (length(keep) == 0L)
    stop(sprintf("kernel PCA (%s): no eigenvalues above tolerance %.2e",
                 kernel, tol))

  k_use    <- min(n_pc, length(keep))
  lam_top  <- lam[keep[seq_len(k_use)]]
  vecs_top <- vecs[, keep[seq_len(k_use)], drop = FALSE]

  scores <- sweep(vecs_top, 2, sqrt(lam_top), `*`)
  rownames(scores) <- rownames(Z_raw)
  Z_scaled <- scale(scores, center = TRUE, scale = TRUE)

  var_proxy <- sum(lam_top) / (sum(lam) + 1e-12)
  eig_cond  <- max(lam_top) / min(lam_top)

  cat(sprintf("  Components retained: %d / %d requested\n", k_use, n_pc))
  cat(sprintf("  Within-kernel spectral fraction (top %d): %.1f%%\n",
              k_use, 100 * var_proxy))
  cat(sprintf("  Eigenvalue condition (lambda_max/lambda_min): %.2e\n", eig_cond))
  cat(sprintf("  n/k ratio: %d / %d = %.2f\n",
              nrow(Z_scaled), k_use, nrow(Z_scaled) / k_use))

  list(Z_scaled    = Z_scaled,
       eigenvalues = lam,
       k_use       = k_use,
       var_proxy   = var_proxy,
       eig_cond    = eig_cond,
       kernel      = kernel,
       agg_method  = agg_method)
}

################################################################
# BLOCK 7: SUBHO runner — single target (MSI)
################################################################

run_subho <- function(Z_scaled, y, agg_method, label = "unknown") {
  cat(sprintf("\n--- SUBHO | %s | %s pooling ---\n", label, agg_method))
  cat(sprintf("  %d patients x %d features | MSI prevalence = %.3f\n",
              nrow(Z_scaled), ncol(Z_scaled), mean(y)))

  set.seed(RANDOM_SEED)
  obs <- subho_MI_logistic(
    Z             = Z_scaled,
    y             = y,
    v_beta        = V_BETA,
    method.weight = "beta.binomial",
    nmc = NMC, burn = BURN, thin = THIN
  )
  cat(sprintf("  Observed MI = %.5f nats | Features selected = %d\n",
              obs$MI, obs$n.selected))

  cat(sprintf("  Running %d permutations...\n", N_PERM))
  MI_null <- permutation_test_logistic(
    Z             = Z_scaled,
    y             = y,
    n.perm        = N_PERM,
    v_beta        = V_BETA,
    method.weight = "beta.binomial",
    nmc = NMC, burn = BURN, thin = THIN,
    seed          = RANDOM_SEED,
    verbose       = TRUE
  )

  n_valid <- sum(!is.na(MI_null))
  p_val   <- (1 + sum(MI_null >= obs$MI, na.rm = TRUE)) / (1 + n_valid)
  cat(sprintf("  Permutation p = %.4f  (add-one corrected)\n", p_val))

  list(
    Z_scaled   = Z_scaled,
    y          = y,
    label      = label,
    agg_method = agg_method,
    result = list(
      summary = data.frame(
        Target      = "MSI",
        MI          = round(obs$MI, 5),
        N_selected  = obs$n.selected,
        P_value     = round(p_val, 4),
        Significant = p_val < 0.05,
        stringsAsFactors = FALSE
      ),
      details = list(MSI = list(
        target     = "MSI",
        MI         = obs$MI,
        n.selected = obs$n.selected,
        MI.null    = MI_null,
        p.value    = p_val
      ))
    )
  )
}

################################################################
# BLOCK 8: Covariate adjustment
#
# Same covariates as the original CRC pipeline: log(tile_count),
# AJCC stage (ordinal), site (COAD vs READ) if retained. R5
# standardizes log_tile_count alongside stage_ordinal_z, which
# the original script left unstandardized.
################################################################

STAGE_ORDINAL <- c(
  "Stage IA"   = 1, "Stage I"    = 1,
  "Stage IIA"  = 2, "Stage IIB"  = 2, "Stage IIC" = 2, "Stage II" = 2,
  "Stage IIIA" = 3, "Stage IIIB" = 3, "Stage IIIC" = 3, "Stage III" = 3,
  "Stage IVA"  = 4, "Stage IVB"  = 4, "Stage IV"  = 4
)

build_covariate_matrix <- function(meta_path, slide_ids, tile_counts) {
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  cat(sprintf("Covariate metadata columns: %s\n",
              paste(names(meta), collapse = ", ")))

  meta_sl <- meta[!is.na(meta[[COL_SLIDE]]) & !duplicated(meta[[COL_SLIDE]]),
                  c(COL_SLIDE, "ajcc_pathologic_stage", "STUDY_ORIG")]
  rownames(meta_sl) <- meta_sl[[COL_SLIDE]]

  matched   <- slide_ids[slide_ids %in% rownames(meta_sl)]
  n_missing <- length(slide_ids) - length(matched)
  if (n_missing > 0)
    cat(sprintf("  [WARN] %d slide(s) missing from covariate metadata\n",
                n_missing))

  df               <- meta_sl[matched, ]
  df$stage_ordinal <- STAGE_ORDINAL[df$ajcc_pathologic_stage]
  n_miss_stage     <- sum(is.na(df$stage_ordinal))
  if (n_miss_stage > 0) {
    df$stage_ordinal[is.na(df$stage_ordinal)] <-
      median(df$stage_ordinal, na.rm = TRUE)
    cat(sprintf("  [WARN] %d patients with unmapped AJCC stage — median imputed\n",
                n_miss_stage))
  }

  df$site_read <- as.numeric(df$STUDY_ORIG == "TCGA-READ")
  n_read       <- sum(df$site_read)

  # R5: log_tile_count standardized alongside stage_ordinal_z (the
  # original script left log_tile_count on its raw log1p scale,
  # giving it a different effective prior weight than the
  # unit-variance embedding PCs under the common V_BETA).
  if (n_read == 0) {
    cat("  [NOTE] No TCGA-READ slides retained — site covariate excluded\n")
    C <- cbind(as.numeric(scale(log1p(tile_counts[matched]))),
               as.numeric(scale(df$stage_ordinal)))
    colnames(C) <- c("log_tile_count_z", "stage_ordinal_z")
  } else {
    C <- cbind(as.numeric(scale(log1p(tile_counts[matched]))),
               as.numeric(scale(df$stage_ordinal)),
               df$site_read)
    colnames(C) <- c("log_tile_count_z", "stage_ordinal_z", "site_READ")
  }
  rownames(C) <- matched
  cat(sprintf("Covariate matrix: %d x %d\n", nrow(C), ncol(C)))
  C
}

################################################################
# compute_incremental_mi(): single-target version. R4 adds
# method.weight = "beta.binomial" to all three SUBHO calls so the
# main MI and the covariate-adjustment MI use the same model
# prior — without this, the delta below conflates a prior
# difference with a true incremental-information difference.
################################################################

compute_incremental_mi <- function(Z_scaled, C, y) {
  common <- Reduce(intersect,
                   list(rownames(Z_scaled), rownames(C), names(y)))
  if (length(common) < 20)
    stop(sprintf("Fewer than 20 patients after covariate alignment: %d",
                 length(common)))
  cat(sprintf("  Covariate alignment: %d patients\n", length(common)))

  Z  <- Z_scaled[common, , drop = FALSE]
  Cm <- C[common,        , drop = FALSE]
  yc <- y[common]

  mi_covar    <- subho_MI_logistic(Cm,           yc, v_beta = V_BETA,
                                   method.weight = "beta.binomial",
                                   nmc = NMC, burn = BURN, thin = THIN)$MI
  mi_emb      <- subho_MI_logistic(Z,            yc, v_beta = V_BETA,
                                   method.weight = "beta.binomial",
                                   nmc = NMC, burn = BURN, thin = THIN)$MI
  mi_combined <- subho_MI_logistic(cbind(Z, Cm), yc, v_beta = V_BETA,
                                   method.weight = "beta.binomial",
                                   nmc = NMC, burn = BURN, thin = THIN)$MI

  # R4: named incremental_subho_mi and documented as a model-based
  # recoverability delta — equals population conditional MI
  # I(Y;Z|C) only if subho_MI_logistic() is a coherent estimator
  # of the same MI functional across differing model dimensions,
  # which is not guaranteed under Bayesian variable selection with
  # a horseshoe prior.
  incremental_subho_mi <- mi_combined - mi_covar

  cat(sprintf(
    "  Covar only: %.4f | Emb only: %.4f | Combined: %.4f | Incremental SUBHO MI: %+.4f nats\n",
    mi_covar, mi_emb, mi_combined, incremental_subho_mi))

  list(covar_only           = mi_covar,
       emb_only              = mi_emb,
       combined               = mi_combined,
       incremental_subho_mi   = incremental_subho_mi)
}

################################################################
# BLOCK 9: Figures A, B, C — single target, single pooling method
# per call.
#
# The original CRC script called this with a "pseudo-pair":
# the same pooling method's result passed into both an out_mean
# and out_norm slot, so e.g. the attention-pooling figure showed
# identical data in two panels hardcoded as "Mean pooling" and
# "Norm-weighted pooling" — display labels that did not describe
# what was actually plotted. Since every call site in this script
# only ever has one pooling method's result to show at a time,
# make_figures() now takes that single result and its true pool
# name directly, producing a one-panel Fig A/B/C rather than
# faking a two-panel comparison that was never real.
################################################################

make_figures <- function(out, pool_name, incremental, C,
                         decomp_label = "kPCA",
                         save_dir     = RESULTS_DIR) {

  pool_colour <- PALETTE[[pool_name]]
  pool_label  <- AGG_LABEL[[pool_name]]

  s   <- out$result$summary
  det <- out$result$details[["MSI"]]
  df_summary <- data.frame(
    agg_method   = pool_name,
    observed_mi  = s$MI[s$Target == "MSI"],
    null_mi_mean = mean(det$MI.null, na.rm = TRUE),
    null_mi_sd   = sd(det$MI.null,   na.rm = TRUE),
    p_value      = s$P_value[s$Target == "MSI"],
    stringsAsFactors = FALSE
  )

  p_A <- ggplot(df_summary,
                aes(x = factor(agg_method, levels = pool_name),
                    y = observed_mi)) +
    geom_col(width = 0.50, colour = "white", linewidth = 0.4,
             fill = pool_colour) +
    geom_errorbar(aes(ymin = null_mi_mean - null_mi_sd,
                      ymax = null_mi_mean + null_mi_sd),
                  width = 0.25, colour = PALETTE$null, linewidth = 1.2) +
    geom_text(aes(y     = observed_mi * 1.04,
                  label = ifelse(p_value < 0.001, "p < 0.001",
                                 sprintf("p = %.3f", p_value))),
              size = 3.2, fontface = "italic") +
    scale_x_discrete(labels = setNames(pool_label, pool_name)) +
    scale_y_continuous(labels = number_format(accuracy = 0.001)) +
    labs(x = NULL, y = "Mutual information (nats)",
         title = sprintf("A   Recoverability — %s (fixed-decomposition p-value)  [%s, %s, %d components]",
                         pool_label, decomp_label, EMBED_TYPE, N_PC)) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 11))

  null_df <- data.frame(null_mi = det$MI.null)

  p_B <- ggplot(null_df, aes(x = null_mi)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 45, alpha = 0.65, colour = NA,
                   fill = pool_colour) +
    geom_vline(xintercept = df_summary$observed_mi,
               colour = pool_colour, linewidth = 1.5, linetype = "dashed") +
    scale_x_continuous(labels = number_format(accuracy = 0.0001)) +
    labs(x = "MI under null (nats)", y = "Density",
         title = sprintf("B   Permutation null distribution — %s (n = %d)",
                         pool_label, N_PERM)) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

  roman <- c("i", "ii", "iii")
  df_C  <- data.frame(
    bar   = roman,
    value = c(incremental$covar_only, incremental$emb_only, incremental$combined),
    stringsAsFactors = FALSE
  )
  df_C$bar <- factor(df_C$bar, levels = roman)
  delta    <- incremental$incremental_subho_mi
  y_top    <- max(df_C$value) + 0.008

  covar_label <- paste(colnames(C), collapse = ", ")

  p_C <- ggplot(df_C, aes(x = bar, y = value, alpha = bar)) +
    geom_col(width = 0.55, colour = "white", linewidth = 0.4, fill = pool_colour) +
    geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey50") +
    annotate("text", x = 2, y = y_top, vjust = 0, size = 3.0,
             fontface = "italic", colour = pool_colour,
             label = sprintf("Incremental SUBHO MI = %+.3f", delta)) +
    scale_alpha_manual(
      values = c(i = 0.38, ii = 1.00, iii = 0.68),
      labels = c("(i) Covariates only", "(ii) Embedding only",
                 "(iii) Embedding + covariates"),
      name   = NULL) +
    scale_y_continuous(labels = number_format(accuracy = 0.001)) +
    labs(x = NULL, y = "Mutual information (nats)",
         title = sprintf(
           "C   Incremental model-based MI after covariate adjustment — %s [%s]\n     Covariates: %s",
           pool_label, EMBED_TYPE, covar_label)) +
    theme_classic(base_size = 11) +
    theme(plot.title      = element_text(face = "bold", size = 11),
          legend.position = "right")

  tag    <- gsub("[^A-Za-z0-9_]", "_", tolower(decomp_label))
  out_AB <- file.path(save_dir,
                      sprintf("figures_kpca_crc_%s_%s_AB.png", EMBED_TYPE, tag))
  png(out_AB, width = 12, height = 6, units = "in", res = 300)
  grid.arrange(p_A, p_B, ncol = 2,
    top = grid::textGrob(
      sprintf("%s — SUBHO Logistic MI | TCGA CRC (MSI) | %s", EMBED_TYPE, decomp_label),
      gp = grid::gpar(fontsize = 13, fontface = "bold")))
  dev.off()
  cat(sprintf("Figures A-B -> %s\n", out_AB))

  out_C <- file.path(save_dir,
                     sprintf("figures_kpca_crc_%s_%s_C.png", EMBED_TYPE, tag))
  ggsave(out_C, plot = p_C, width = 7, height = 6, dpi = 300)
  cat(sprintf("Figure C    -> %s\n", out_C))

  invisible(list(A = p_A, B = p_B, C = p_C))
}

################################################################
# BLOCK 10: Scree comparison across all decompositions
################################################################

method_from_key <- function(nm) {
  if (grepl("^pca",    nm)) return("PCA")
  if (grepl("^cosine", nm)) return("cosine")
  if (grepl("^rbf",    nm)) return("rbf")
  if (grepl("^poly",   nm)) return("poly")
  nm
}

pool_from_key <- function(nm) {
  if (grepl("_mean$",      nm)) return("mean")
  if (grepl("_norm$",      nm)) return("norm")
  if (grepl("_moment$",    nm)) return("moment")
  if (grepl("_attention$", nm)) return("attention")
  nm
}

make_scree_comparison <- function(decomp_list, save_dir = RESULTS_DIR) {
  build_df <- function(res, key) {
    lam <- pmax(res$eigenvalues, 0)
    tot <- sum(lam)
    k   <- min(80L, length(lam))
    data.frame(
      PC         = seq_len(k),
      ind_frac   = lam[1:k] / tot,
      cum_frac   = cumsum(lam[1:k]) / tot,
      method     = method_from_key(key),
      agg_method = pool_from_key(key)
    )
  }

  df <- do.call(rbind,
                mapply(build_df, decomp_list, names(decomp_list),
                       SIMPLIFY = FALSE))

  pal <- c("PCA"    = "#2166AC",
           "cosine" = "#D6604D",
           "rbf"    = "#4DAC26",
           "poly"   = "#8B3A8B")
  lty <- c(mean      = "solid",
           norm      = "dashed",
           moment    = "dotdash",
           attention = "longdash")

  p_ind <- ggplot(df, aes(x = PC, y = ind_frac * 100,
                          colour = method, linetype = agg_method)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = pal) +
    scale_linetype_manual(values = lty) +
    labs(x = "Component",
         y = "Within-method spectral fraction (%)",
         colour = "Decomposition", linetype = "Pooling",
         title = "Individual spectrum") +
    theme_classic(base_size = 11) +
    theme(legend.position = "right")

  p_cum <- ggplot(df, aes(x = PC, y = cum_frac * 100,
                          colour = method, linetype = agg_method)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = N_PC, linetype = "dotted", colour = "grey40") +
    geom_hline(yintercept = c(90, 95, 99), linetype = "dotted",
               colour = "grey70", linewidth = 0.4) +
    scale_colour_manual(values = pal) +
    scale_linetype_manual(values = lty) +
    labs(x = "Component",
         y = "Cumulative within-method spectral fraction (%)",
         colour = "Decomposition", linetype = "Pooling",
         title = "Cumulative spectrum") +
    theme_classic(base_size = 11) +
    theme(legend.position = "right")

  out_path <- file.path(save_dir, sprintf("kpca_scree_comparison_crc_%s.png", EMBED_TYPE))
  png(out_path, width = 16, height = 6, units = "in", res = 300)
  grid.arrange(p_ind, p_cum, ncol = 2,
    top = grid::textGrob(
      sprintf("Eigenspectrum — TCGA CRC | %s | %d components | 4 pooling methods",
              EMBED_TYPE, N_PC),
      gp = grid::gpar(fontsize = 13, fontface = "bold")))
  dev.off()
  cat(sprintf("Scree comparison -> %s\n", out_path))
  invisible(list(ind = p_ind, cum = p_cum))
}

################################################################
# BLOCK 11: MI comparison plot — single target, no target facet
################################################################

make_mi_comparison <- function(results_list, save_dir = RESULTS_DIR) {
  df <- do.call(rbind, lapply(names(results_list), function(nm) {
    r   <- results_list[[nm]]
    det <- r$result$details[["MSI"]]
    data.frame(
      key        = nm,
      label      = r$label,
      decomp     = method_from_key(nm),
      agg_method = pool_from_key(nm),
      MI_obs     = det$MI,
      p_value    = det$p.value,
      null_mean  = mean(det$MI.null, na.rm = TRUE),
      null_sd    = sd(det$MI.null,   na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  decomp_order <- c("PCA", "cosine", "rbf", "poly")
  pool_order   <- c("mean", "norm", "moment", "attention")
  df$decomp     <- factor(df$decomp,     levels = decomp_order)
  df$agg_method <- factor(df$agg_method, levels = pool_order)

  pal_decomp  <- c("PCA"    = "#2166AC", "cosine" = "#D6604D",
                   "rbf"    = "#4DAC26", "poly"   = "#8B3A8B")
  pool_labels <- c(mean = "Mean", norm = "Norm",
                   moment = "Moment", attention = "Attention")

  p_mi <- ggplot(df, aes(x = decomp, y = MI_obs, fill = decomp)) +
    geom_bar(stat = "identity", position = position_dodge(0.7),
             width = 0.65, colour = "grey20", linewidth = 0.3) +
    geom_errorbar(
      aes(ymin = null_mean - null_sd, ymax = null_mean + null_sd),
      position = position_dodge(0.7), width = 0.25,
      colour = "grey30", linewidth = 0.5
    ) +
    facet_wrap(~ agg_method, nrow = 1,
               labeller = as_labeller(pool_labels)) +
    scale_fill_manual(values = pal_decomp, name = "Decomposition") +
    labs(x = "Decomposition", y = "MI (nats)",
         title = sprintf("Observed MI with permutation-null mean \u00b1 1 SD [%s]", EMBED_TYPE)) +
    theme_classic(base_size = 10) +
    theme(axis.text.x     = element_text(angle = 30, hjust = 1),
          strip.text      = element_text(face = "bold"),
          legend.position = "bottom")

  p_floor <- 1 / (N_PERM + 1)
  p_pv <- ggplot(df, aes(x = decomp,
                         y = -log10(pmax(p_value, p_floor)),
                         fill = decomp)) +
    geom_bar(stat = "identity", position = position_dodge(0.7),
             width = 0.65, colour = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "firebrick", linewidth = 0.5) +
    facet_wrap(~ agg_method, nrow = 1,
               labeller = as_labeller(pool_labels)) +
    scale_fill_manual(values = pal_decomp, name = "Decomposition") +
    labs(x = "Decomposition",
         y = expression(-log[10](p[corrected])),
         title = sprintf("Permutation p-value (add-one corrected) [%s]", EMBED_TYPE)) +
    theme_classic(base_size = 10) +
    theme(axis.text.x     = element_text(angle = 30, hjust = 1),
          strip.text      = element_text(face = "bold"),
          legend.position = "bottom")

  out_path <- file.path(save_dir, sprintf("kpca_mi_comparison_crc_%s.png", EMBED_TYPE))
  png(out_path, width = 18, height = 8, units = "in", res = 300)
  grid.arrange(p_mi, p_pv, nrow = 2,
    top = grid::textGrob(
      sprintf("MI comparison — 4 pooling x 4 decompositions | TCGA CRC | %s", EMBED_TYPE),
      gp = grid::gpar(fontsize = 13, fontface = "bold")))
  dev.off()
  cat(sprintf("MI comparison -> %s\n", out_path))

  cat("\n=== MI Comparison Summary ===\n")
  print(df[order(df$agg_method, df$decomp),
           c("decomp", "agg_method", "MI_obs", "p_value")],
        row.names = FALSE)

  invisible(df)
}

################################################################
# BLOCK 12: Select best decomposition
#
# R2 fix. The original CRC select_best_decomp() filtered
# candidates with `eligible <- sub[sub$p_value < 0.05 & ...]`
# before computing the nested-permutation null over only those
# survivors — the same outcome-adaptive selection bug as PRAD's
# any_sig filter, just expressed with one target's p_value
# instead of pooling across two targets' significance flags.
#
# Fix: eligibility is now ONLY the pre-specified eig_cond numeric
# filter (a property of the decomposition itself, fixed before
# looking at the label at all). maxT null is built over the full
# eligible candidate set unconditionally. p_maxT replaces
# p_nested as the column name to signal the corrected procedure.
################################################################

select_best_decomp <- function(results_list, decomp_list,
                               save_dir = RESULTS_DIR) {
  df <- do.call(rbind, lapply(names(results_list), function(nm) {
    r   <- results_list[[nm]]
    det <- r$result$details[["MSI"]]
    data.frame(
      key        = nm,
      label      = r$label,
      agg_method = pool_from_key(nm),
      MI_obs     = det$MI,
      p_value    = det$p.value,
      eig_cond   = decomp_list[[nm]]$eig_cond,
      stringsAsFactors = FALSE
    )
  }))

  EIG_COND_MAX <- 1e10
  best_list    <- list()

  for (pool in c("mean", "norm", "moment", "attention")) {
    sub <- df[df$agg_method == pool, ]
    if (nrow(sub) == 0L) next

    # Eligibility is purely the eig_cond numeric filter, fixed
    # before looking at any label. No p_value < 0.05 term here,
    # unlike the original script (R2).
    elig_keys <- unique(sub$key[sub$eig_cond < EIG_COND_MAX])
    if (length(elig_keys) == 0L) {
      cat(sprintf("  [%s] No method passed eig_cond < %.0e; falling back to all (exploratory).\n",
                  pool, EIG_COND_MAX))
      elig_keys <- unique(sub$key)
    }

    sub_elig <- sub[sub$key %in% elig_keys, ]
    best_row <- sub_elig[which.max(sub_elig$MI_obs), ]
    best_key <- best_row$key
    T_obs    <- best_row$MI_obs

    # maxT null over the SAME eligible candidate set used for
    # selection, unconditional on observed significance (R2).
    null_vecs <- lapply(results_list[elig_keys], function(r)
      r$result$details[["MSI"]]$MI.null)
    B        <- min(sapply(null_vecs, length))
    null_mat <- do.call(cbind, lapply(null_vecs, `[`, 1:B))
    T_null   <- apply(null_mat, 1, max, na.rm = TRUE)
    n_valid  <- sum(!is.na(T_null))
    p_maxT   <- (1 + sum(T_null >= T_obs, na.rm = TRUE)) / (1 + n_valid)

    cat(sprintf(
      "  [%-9s] Best: %-30s | MI=%.5f | p_single=%.4f | p_maxT=%.4f (n_candidates=%d)\n",
      pool, best_row$label, T_obs, best_row$p_value, p_maxT, length(elig_keys)))

    best_list[[pool]] <- data.frame(
      key          = best_key,
      label        = best_row$label,
      agg_method   = pool,
      MI_obs       = T_obs,
      eig_cond     = best_row$eig_cond,
      p_maxT       = p_maxT,
      n_candidates = length(elig_keys),
      stringsAsFactors = FALSE
    )
  }

  best_df <- do.call(rbind, best_list)
  write.csv(best_df,
            file.path(save_dir, sprintf("kpca_selected_decomposition_crc_%s.csv", EMBED_TYPE)),
            row.names = FALSE)
  cat(sprintf("  Selection written to kpca_selected_decomposition_crc_%s.csv\n", EMBED_TYPE))
  cat("  p_maxT is the selection-adjusted p-value for the per-pooling-group\n")
  cat("  decomposition pick. Fixed-decomposition p-values reported elsewhere\n")
  cat("  (run_subho summaries, Figure A) are NOT valid for whichever\n")
  cat("  decomposition ends up selected here.\n")
  best_df
}

################################################################
# BLOCK 13: Main execution
################################################################

set.seed(RANDOM_SEED)

cat("\n", strrep("=", 70), "\n")
cat(sprintf("CRC Trident embeddings | %s | Kernel PCA comparison | SUBHO Logistic MI\n", EMBED_TYPE))
cat("Target: MSIStatus (MSI-H)\n")
cat(strrep("=", 70), "\n\n")

# 13.0  Index + key sanity check before touching the full roster
cat("=== Indexing embedding directory and checking h5 keys ===\n")
stem_map <- index_embed_dir(EMBED_DIR)
discover_h5_keys(EMBED_DIR, stem_map)
diagnose_slide_matching(META_PATH, stem_map)

# 13.1  Unified patient roster
cat("\n=== Building unified patient roster ===\n")
roster      <- build_unified_roster(META_PATH, EMBED_DIR, stem_map)
y           <- roster$label_vec
tile_counts <- roster$tile_counts

cat(sprintf("\nFinal cohort: n = %d patients\n", length(y)))

# 13.2  Raw embedding matrices
Z_raw_mean      <- do.call(rbind, roster$emb_mean)
Z_raw_norm      <- do.call(rbind, roster$emb_norm)
Z_raw_moment    <- do.call(rbind, roster$emb_moment)
Z_raw_attention <- do.call(rbind, roster$emb_attention)

stopifnot(identical(rownames(Z_raw_mean), rownames(Z_raw_norm)))
stopifnot(identical(rownames(Z_raw_mean), rownames(Z_raw_moment)))
stopifnot(identical(rownames(Z_raw_mean), rownames(Z_raw_attention)))

cat("Roster check:\n")
cat(sprintf("  mean:      %d x %d\n", nrow(Z_raw_mean),      ncol(Z_raw_mean)))
cat(sprintf("  norm:      %d x %d\n", nrow(Z_raw_norm),      ncol(Z_raw_norm)))
cat(sprintf("  moment:    %d x %d  (2d — mean+SD concatenated)\n",
            nrow(Z_raw_moment), ncol(Z_raw_moment)))
cat(sprintf("  attention: %d x %d\n", nrow(Z_raw_attention), ncol(Z_raw_attention)))

y <- y[rownames(Z_raw_mean)]

# 13.3  All decompositions: 4 pooling x 4 kernels = 16 total
cat("\n=== Computing decompositions (16 total) ===\n")

raw_map <- list(mean = Z_raw_mean, norm = Z_raw_norm,
                moment = Z_raw_moment, attention = Z_raw_attention)

pooling_methods <- c("mean", "norm", "moment", "attention")
kernel_methods  <- c("pca", "cosine", "rbf", "poly")

decomp    <- list()
label_map <- c()

for (pool in pooling_methods) {
  for (kern in kernel_methods) {
    key <- sprintf("%s_%s", kern, pool)
    cat(sprintf("\n>> %s\n", key))

    if (kern == "pca") {
      decomp[[key]] <- apply_pca(raw_map[[pool]], agg_method = pool)
    } else {
      decomp[[key]] <- apply_kpca(raw_map[[pool]], kernel = kern,
                                  agg_method = pool)
    }

    kern_label <- switch(kern, pca = "Standard PCA", cosine = "Cosine kPCA",
                         rbf = "RBF kPCA",     poly   = "Poly kPCA")
    pool_label <- switch(pool, mean = "Mean", norm = "Norm",
                         moment = "Moment", attention = "Attention")
    label_map[key] <- sprintf("%s / %s", kern_label, pool_label)
  }
}

rm(Z_raw_mean, Z_raw_norm, Z_raw_moment, Z_raw_attention, raw_map); gc()

# 13.4  Eigenspectrum comparison
make_scree_comparison(decomp)

# 13.5  SUBHO on all 16 decompositions
cat(sprintf("\n=== Running SUBHO on all 16 decompositions [%s] ===\n", EMBED_TYPE))
results <- list()
for (nm in names(decomp)) {
  results[[nm]] <- run_subho(
    Z_scaled   = decomp[[nm]]$Z_scaled,
    y          = y,
    agg_method = pool_from_key(nm),
    label      = label_map[nm]
  )
}

saveRDS(results, file.path(RESULTS_DIR, sprintf("kpca_subho_results_crc_%s.rds", EMBED_TYPE)))
cat(sprintf("Results saved -> %s\n",
            file.path(RESULTS_DIR, sprintf("kpca_subho_results_crc_%s.rds", EMBED_TYPE))))

# 13.6  MI comparison
mi_df <- make_mi_comparison(results)
write.csv(mi_df, file.path(RESULTS_DIR, sprintf("kpca_mi_summary_crc_%s.csv", EMBED_TYPE)),
          row.names = FALSE)

# 13.7  Select best decomposition per pooling method
cat("\n=== Selecting best decomposition per pooling method ===\n")
best_df <- select_best_decomp(results, decomp)
write.csv(best_df, file.path(RESULTS_DIR, "best_kpca_mi_summary_crc.csv"),
          row.names = FALSE)

# 13.8  Downstream analysis (incremental MI + Figs A-C) per pooling method
cat("\n=== Downstream analysis on selected decompositions ===\n")

any_key <- best_df$key[1]
C <- build_covariate_matrix(META_PATH,
                            rownames(results[[any_key]]$Z_scaled),
                            tile_counts)

all_incremental <- list()
for (pool in best_df$agg_method) {
  best_key <- best_df$key[best_df$agg_method == pool]
  out_best <- results[[best_key]]

  cat(sprintf("\n  [%s] Incremental MI for best: %s\n",
              pool, best_df$label[best_df$agg_method == pool]))
  all_incremental[[pool]] <- compute_incremental_mi(out_best$Z_scaled, C, y)
}

write.csv(
  do.call(rbind, lapply(names(all_incremental), function(m)
    cbind(agg_method = m, as.data.frame(all_incremental[[m]])))),
  file.path(RESULTS_DIR, sprintf("incremental_mi_crc_%s.csv", EMBED_TYPE)),
  row.names = FALSE
)

# Figures A-C: one set per pooling method using its own best
# decomposition, each figure showing only that pool's true data
# (see BLOCK 9 note — no faked mean/norm pair).
for (pool in best_df$agg_method) {
  best_key <- best_df$key[best_df$agg_method == pool]
  out_best <- results[[best_key]]
  incr     <- all_incremental[[pool]]
  dlabel   <- best_df$label[best_df$agg_method == pool]

  make_figures(out_best, pool, incr, C,
               decomp_label = sprintf("%s [%s pooling]", dlabel, pool))
}

# Final summary
cat(sprintf("\n%s\n", strrep("=", 70)))
cat(sprintf("CRC Trident (%s) kernel PCA pipeline complete.\n", EMBED_TYPE))
cat("p_maxT is the selection-adjusted p-value for each pool's pick.\n")
cat("Fixed-decomposition p-values from run_subho()/Figure A are valid\n")
cat("only before this selection step.\n\n")
for (pool in best_df$agg_method)
  cat(sprintf("Best [%-9s]: %-35s | MI=%.5f | p_maxT=%.4f\n",
              pool, best_df$label[best_df$agg_method == pool],
              best_df$MI_obs[best_df$agg_method == pool],
              best_df$p_maxT[best_df$agg_method == pool]))
cat(sprintf("%s\n", strrep("=", 70)))


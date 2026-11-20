################################################################
# kPCA_SUBHO_PRAD_trident.R
# Kernel PCA — SUBHO Logistic MI Pipeline | PRAD Trident embeddings
#
# Structurally identical to kPCA_SUBHO_PRAD.R. Differences:
#   - Embeddings come from Trident, saved as one .h5 per slide.
#     Each h5 holds two datasets: 'features' (N_tiles x D) and
#     'coords' (N_tiles x 2). Key names are discovered at
#     runtime and asserted, not hardcoded blind.
#   - EMBED_TYPE switches between "uni_v2" (D=1536) and
#     "virchow2" (D=2560) so the same script drives both runs.
#     Set EMBED_TYPE below and source the whole file per run.
#   - Output dirs/file tags are suffixed with EMBED_TYPE so the
#     two runs never overwrite each other.
#
# Statistical fixes applied in this pass (external review,
# see commit notes / chat log for full discussion):
#   R1. Pooled embeddings are column-standardized before all
#       three kernels (cosine, RBF, poly), not just linear PCA.
#       Cosine on standardized features is angle w.r.t. centred,
#       unit-variance feature axes, not raw-magnitude angle.
#   R2. select_best_decomp() no longer filters candidate keys by
#       observed significance (any_sig) before building the null
#       maximum. maxT null is computed over the full pooling-group
#       candidate set unconditionally, separately per target.
#       Best key is also chosen per-target now, not by MI_obs
#       pooled across both ERG and ETS (which let one target's
#       signal pick the decomposition reported for the other).
#   R3. make_figures() takes explicit pool_names instead of
#       hardcoded "mean"/"norm" labels, fixing the bug where the
#       moment/attention figure call had its panels labeled as
#       mean/norm pooling.
#   R4. compute_incremental_mi() now passes method.weight =
#       "beta.binomial" to all three subho_MI_logistic() calls,
#       matching run_subho(). Output is named
#       incremental_subho_mi and documented as a model-based
#       recoverability delta, not asserted to equal population
#       conditional MI I(Y;Z|C).
#   R5. log_tile_count is standardized in build_covariate_matrix()
#       alongside Gleason_sum_z, so both covariates carry the same
#       effective prior weight as the unit-scaled embedding PCs.
#
# Everything else from BLOCK 4 onward is otherwise unchanged from
# kPCA_SUBHO_PRAD.R.
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
# EMBED_TYPE switch — set this and source the file once per run
################################################################
EMBED_TYPE <- "uni_v2"  #"virchow2"   # "uni_v2"  | "virchow2"

stopifnot(EMBED_TYPE %in% c("uni_v2", "virchow2"))

EMBED_DIM <- switch(EMBED_TYPE, uni_v2 = 1536L, virchow2 = 2560L)

EMBED_DIR <- switch(EMBED_TYPE,
  uni_v2   = "M:/shared_data/TCGA_PRAD/trident-feats/20x_256px_0px_overlap/features_uni_v2/",
  virchow2 = "M:/shared_data/TCGA_PRAD/trident-feats/20x_224px_0px_overlap/features_virchow2/"
)

META_PATH <- "M:/shared_data/TCGA_PRAD/orig_annotations.csv"
RESULTS_DIR <- sprintf("M:/bhadury/mi_results/TCGA_PRAD_kpca_subho_%s/", EMBED_TYPE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("EMBED_TYPE = %s | EMBED_DIM = %d | EMBED_DIR = %s\n",
            EMBED_TYPE, EMBED_DIM, EMBED_DIR))

################################################################
# CONFIG — identical to PRAD .pt pipeline
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
COL_ERG   <- "ERG_status";  ERG_POS <- "fusion"
COL_ETS   <- "ETS";         ETS_POS <- "Positive"

TARGETS_DEF <- list(
  ERG = list(col = COL_ERG, pos_level = ERG_POS),
  ETS = list(col = COL_ETS, pos_level = ETS_POS)
)
TARGETS_V <- c("ERG", "ETS")

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
# Confirmed empirically: a test file's actual filename used an
# underscore between the slide barcode and the UUID
# (...DX1_beeda9ae-...), while both its internal attrs['name']
# and the orig_annotations.csv 'slide' column used a dot
# (...DX1.beeda9ae-...). If real M:/ filenames follow the
# underscore convention while metadata uses dots, the old
# filename-stem approach silently mismatched every slide and
# logged it as [SKIP - not found] with no further explanation.
# Reading the attribute directly removes the dependence on
# filename punctuation entirely; lookup now matches whatever
# convention orig_annotations.csv uses, since both are reading
# the same canonical Trident-assigned name.
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
# Run once after build_unified_roster() reports far fewer
# patients than expected. Compares the metadata 'slide' values
# against the stem_map keys directly so a naming-convention
# mismatch shows up immediately as a string diff, rather than
# being inferred indirectly from a low patient count.
################################################################

diagnose_slide_matching <- function(meta_path, stem_map) {
  meta   <- read.csv(meta_path, stringsAsFactors = FALSE)
  slides <- meta[[COL_SLIDE]]
  slides <- slides[!is.na(slides)]

  in_dir     <- slides %in% names(stem_map)
  n_match    <- sum(in_dir)
  n_total    <- length(slides)

  cat(sprintf("\n[diagnose_slide_matching] %d / %d metadata slides found in directory index.\n",
              n_match, n_total))

  if (n_match < n_total) {
    missing <- slides[!in_dir]
    cat(sprintf("  %d metadata slide(s) NOT found. First few:\n", length(missing)))
    for (s in head(missing, 5)) cat(sprintf("    metadata: %s\n", s))

    cat("  First few directory-indexed names for comparison:\n")
    for (s in head(names(stem_map), 5)) cat(sprintf("    directory: %s\n", s))

    # Case-insensitive re-check, purely diagnostic — does NOT
    # change the matching behaviour used in build_unified_roster().
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
# Confirmed empirically: keys found are 'coords' and 'features',
# not 'coordinates'/'features' as initially guessed from the
# Trident notes. 'coords' is not required for pooling (only
# 'features' is read), so this check only flags if it's absent
# in case anything is added later that does use tile positions.
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
# Read directly with h5py — no zip/torch workaround needed, this
# is a plain HDF5 file. Pooling logic identical to the .pt loader.
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
# Identical to kPCA_SUBHO_PRAD.R, calling load_h5_pooled instead
# of load_npy_pooled.
################################################################

build_unified_roster <- function(meta_path, embed_dir, stem_map) {
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)

  for (tgt in names(TARGETS_DEF)) {
    col       <- TARGETS_DEF[[tgt]]$col
    pos_level <- TARGETS_DEF[[tgt]]$pos_level
    meta[[paste0(tgt, "_bin")]] <- as.integer(meta[[col]] == pos_level)
    cat(sprintf("%s_bin: %d pos / %d neg / %d NA\n",
                tgt,
                sum(meta[[paste0(tgt, "_bin")]] == 1L, na.rm = TRUE),
                sum(meta[[paste0(tgt, "_bin")]] == 0L, na.rm = TRUE),
                sum(is.na(meta[[paste0(tgt, "_bin")]]))))
  }

  emb_mean      <- list()
  emb_norm      <- list()
  emb_moment    <- list()
  emb_attention <- list()
  label_accum   <- lapply(TARGETS_DEF, function(x) integer(0))
  id_accum      <- character(0)
  tile_counts   <- c()

  agg_methods <- c("mean", "norm", "moment", "attention")

  for (i in seq_len(nrow(meta))) {
    slide_id <- meta[[COL_SLIDE]][i]

    if (is.na(slide_id)) {
      cat(sprintf("  [SKIP - NA slide] row %d\n", i)); next
    }

    labels_i <- sapply(names(TARGETS_DEF), function(tgt)
      meta[[paste0(tgt, "_bin")]][i])
    if (any(is.na(labels_i))) {
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
    for (tgt in names(TARGETS_DEF))
      label_accum[[tgt]] <- c(label_accum[[tgt]], as.integer(labels_i[[tgt]]))
  }

  label_list <- lapply(label_accum, function(v) {
    names(v) <- id_accum
    v
  })
  n_ret      <- length(tile_counts)

  cat(sprintf("Unified roster: %d patients\n", n_ret))
  for (tgt in names(TARGETS_DEF))
    cat(sprintf("  %s: %d pos / %d neg\n",
                tgt, sum(label_list[[tgt]]), sum(label_list[[tgt]] == 0)))
  cat(sprintf("Tile counts — min: %d | median: %.0f | max: %d\n",
              min(tile_counts), median(tile_counts), max(tile_counts)))
  cat(sprintf("Embedding dim — mean/norm/attn: %d | moment: %d\n",
              length(emb_mean[[1]]), length(emb_moment[[1]])))

  list(emb_mean      = emb_mean,
       emb_norm      = emb_norm,
       emb_moment    = emb_moment,
       emb_attention = emb_attention,
       label_list    = label_list,
       tile_counts   = tile_counts)
}

################################################################
# BLOCK 4: Kernel constructors
#
# R1 fix: pooled embeddings are column-standardized before all
# three kernels, not just linear PCA. Without this, RBF and poly
# are dominated by raw embedding magnitude (a real problem when
# comparing UNI v2, 1024-dim, against Virchow2, 2560-dim, which
# sit at different native scales). Cosine is recomputed on the
# same standardized matrix for consistency: this measures angle
# w.r.t. centred, unit-variance feature axes, which is a
# different geometric object from raw-magnitude cosine. If a
# raw-scale cosine comparison is ever wanted alongside this,
# add it as a fifth decomposition rather than reverting this one.
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
# BLOCK 5: Standard linear PCA — unchanged
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

  # R1: standardize before ANY kernel (cosine included) so all
  # three are on a comparable footing and not scale-dominated.
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
# BLOCK 7: SUBHO runner — unchanged
################################################################

run_subho <- function(Z_scaled, y_list, agg_method, label = "unknown") {
  cat(sprintf("\n--- SUBHO | %s | %s pooling ---\n", label, agg_method))
  cat(sprintf("  %d patients x %d features\n", nrow(Z_scaled), ncol(Z_scaled)))
  for (tgt in TARGETS_V)
    cat(sprintf("  %s prevalence = %.3f\n", tgt, mean(y_list[[tgt]])))

  details <- list()

  for (tgt in TARGETS_V) {
    y <- y_list[[tgt]]
    cat(sprintf("\n  [%s]\n", tgt))

    set.seed(RANDOM_SEED)
    obs <- subho_MI_logistic(
      Z             = Z_scaled,
      y             = y,
      v_beta        = V_BETA,
      method.weight = "beta.binomial",
      nmc = NMC, burn = BURN, thin = THIN
    )
    cat(sprintf("    Observed MI = %.5f nats | Features selected = %d\n",
                obs$MI, obs$n.selected))

    cat(sprintf("    Running %d permutations...\n", N_PERM))
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
    p_raw   <- (1 + sum(MI_null >= obs$MI, na.rm = TRUE)) / (1 + n_valid)
    cat(sprintf("    Permutation p = %.4f  (add-one corrected)\n", p_raw))

    details[[tgt]] <- list(
      target     = tgt,
      MI         = obs$MI,
      n.selected = obs$n.selected,
      MI.null    = MI_null,
      p.value    = p_raw
    )
  }

  p_raw_vec        <- sapply(TARGETS_V, function(tgt) details[[tgt]]$p.value)
  p_adj_vec        <- p.adjust(p_raw_vec, method = "BH")
  names(p_adj_vec) <- TARGETS_V

  summary_df <- do.call(rbind, lapply(TARGETS_V, function(tgt) {
    data.frame(
      Target      = tgt,
      MI          = round(details[[tgt]]$MI, 5),
      N_selected  = details[[tgt]]$n.selected,
      P_value     = round(details[[tgt]]$p.value, 4),
      P_adj_BH    = round(p_adj_vec[[tgt]], 4),
      Significant = p_adj_vec[[tgt]] < 0.05,
      stringsAsFactors = FALSE
    )
  }))

  cat(sprintf("\n  Summary [%s | %s pooling]:\n", label, agg_method))
  print(summary_df, row.names = FALSE)

  list(
    Z_scaled   = Z_scaled,
    y_list     = y_list,
    label      = label,
    agg_method = agg_method,
    result     = list(summary = summary_df, details = details)
  )
}

################################################################
# BLOCK 8: Covariate adjustment — unchanged
################################################################

build_covariate_matrix <- function(meta_path, slide_ids, tile_counts) {
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  cat(sprintf("Covariate metadata columns: %s\n",
              paste(names(meta), collapse = ", ")))

  meta_sl <- meta[!is.na(meta[[COL_SLIDE]]) & !duplicated(meta[[COL_SLIDE]]),
                  c(COL_SLIDE, "Reviewed_Gleason_sum")]
  rownames(meta_sl) <- meta_sl[[COL_SLIDE]]

  matched   <- slide_ids[slide_ids %in% rownames(meta_sl)]
  n_missing <- length(slide_ids) - length(matched)
  if (n_missing > 0)
    cat(sprintf("  [WARN] %d slide(s) missing from covariate metadata\n",
                n_missing))

  df           <- meta_sl[matched, ]
  n_miss_gleas <- sum(is.na(df$Reviewed_Gleason_sum))
  if (n_miss_gleas > 0) {
    med <- median(df$Reviewed_Gleason_sum, na.rm = TRUE)
    df$Reviewed_Gleason_sum[is.na(df$Reviewed_Gleason_sum)] <- med
    cat(sprintf("  [WARN] %d patients missing Gleason sum — median imputed (%.1f)\n",
                n_miss_gleas, med))
  }

  # R5: log_tile_count standardized to match Gleason_sum_z and
  # the unit-variance embedding PCs in Z_scaled, so all columns
  # entering subho_MI_logistic() carry the same effective prior
  # weight under the common V_BETA.
  C <- cbind(as.numeric(scale(log1p(tile_counts[matched]))),
             as.numeric(scale(df$Reviewed_Gleason_sum)))
  colnames(C) <- c("log_tile_count_z", "Gleason_sum_z")
  rownames(C) <- matched
  cat(sprintf("Covariate matrix: %d x %d\n", nrow(C), ncol(C)))
  C
}

################################################################
# compute_incremental_mi_one_target(): per-target worker.
#
# Needed because select_best_decomp() (R2) now selects a
# decomposition separately per target — mean pooling, say, can
# have a different best Z_scaled for ERG than for ETS. Each
# target's incremental MI must therefore be computed against its
# own selected Z_scaled, not a single shared one.
################################################################

compute_incremental_mi_one_target <- function(Z_scaled, C, y_target, tgt_label) {
  common <- Reduce(intersect, list(rownames(Z_scaled), rownames(C), names(y_target)))
  if (length(common) < 20)
    stop(sprintf("Fewer than 20 patients after covariate alignment: %d",
                 length(common)))
  cat(sprintf("  Covariate alignment: %d patients\n", length(common)))

  Z  <- Z_scaled[common, , drop = FALSE]
  Cm <- C[common,        , drop = FALSE]
  yc <- y_target[common]

  cat(sprintf("  [%s] ", tgt_label))

  # R4: method.weight added to match run_subho() exactly — the
  # main MI and the covariate-adjustment MI must use the same
  # model prior or the delta below conflates a prior difference
  # with a true incremental-information difference.
  mi_covar    <- subho_MI_logistic(Cm,           yc, v_beta = V_BETA,
                                   method.weight = "beta.binomial",
                                   nmc = NMC, burn = BURN, thin = THIN)$MI
  mi_emb      <- subho_MI_logistic(Z,            yc, v_beta = V_BETA,
                                   method.weight = "beta.binomial",
                                   nmc = NMC, burn = BURN, thin = THIN)$MI
  mi_combined <- subho_MI_logistic(cbind(Z, Cm), yc, v_beta = V_BETA,
                                   method.weight = "beta.binomial",
                                   nmc = NMC, burn = BURN, thin = THIN)$MI

  # R4: this is mi_combined - mi_covar, the increase in
  # SUBHO-estimated MI from adding the embedding to a
  # covariate-only model. It equals the population conditional
  # MI I(Y;Z|C) only if subho_MI_logistic() is a coherent
  # estimator of the same MI functional across differing model
  # dimensions and selection paths, which is not guaranteed
  # under Bayesian variable selection with a horseshoe prior.
  # Reported and named as an incremental model-based
  # recoverability score, not claimed as conditional MI.
  incremental_subho_mi <- mi_combined - mi_covar

  cat(sprintf(
    "Covar only: %.4f | Emb only: %.4f | Combined: %.4f | Incremental SUBHO MI: %+.4f nats\n",
    mi_covar, mi_emb, mi_combined, incremental_subho_mi))

  list(covar_only           = mi_covar,
       emb_only              = mi_emb,
       combined               = mi_combined,
       incremental_subho_mi   = incremental_subho_mi)
}

################################################################
# compute_incremental_mi(): multi-target wrapper, retained for
# any fixed-decomposition use case where one Z_scaled is shared
# across both targets (e.g. a non-selection exploratory run).
# Internally just calls the per-target worker for each target.
################################################################

compute_incremental_mi <- function(Z_scaled, C, y_list) {
  result <- list()
  for (tgt in TARGETS_V)
    result[[tgt]] <- compute_incremental_mi_one_target(Z_scaled, C, y_list[[tgt]], tgt)
  result
}

################################################################
# BLOCK 9: Figures A, B, C
#
# R3 fix: the function used to hardcode "mean"/"norm" as the
# panel labels regardless of which two pooling results were
# actually passed in. The second call site in the main script
# passes moment + attention results, which were silently
# relabeled as mean/norm in the output figure. pool_names now
# carries the true labels through to every facet/legend.
################################################################

make_figures <- function(out1, out2, incremental, C,
                         pool_names   = c("mean", "norm"),
                         decomp_label = "kPCA",
                         save_dir     = RESULTS_DIR) {

  stopifnot(length(pool_names) == 2L, !anyDuplicated(pool_names))
  pn1 <- pool_names[1]; pn2 <- pool_names[2]
  pool_labels_used <- setNames(
    sapply(pool_names, function(p) AGG_LABEL[[p]]),
    pool_names
  )

  TARGET_COLOURS <- c(ERG = "#2166AC", ETS = "#D6604D")

  extract_rows <- function(out, m) {
    s <- out$result$summary
    do.call(rbind, lapply(TARGETS_V, function(tgt) {
      det <- out$result$details[[tgt]]
      data.frame(
        agg_method   = m,
        target       = tgt,
        observed_mi  = det$MI,
        null_mi_mean = mean(det$MI.null, na.rm = TRUE),
        null_mi_sd   = sd(det$MI.null,   na.rm = TRUE),
        p_value      = s$P_value[s$Target == tgt],
        stringsAsFactors = FALSE
      )
    }))
  }
  df_summary <- rbind(extract_rows(out1, pn1),
                      extract_rows(out2, pn2))
  df_summary$target     <- factor(df_summary$target,     levels = TARGETS_V)
  df_summary$agg_method <- factor(df_summary$agg_method, levels = pool_names)

  p_A <- ggplot(df_summary,
                aes(x = target, y = observed_mi, fill = target)) +
    geom_col(width = 0.50, colour = "white", linewidth = 0.4) +
    geom_errorbar(aes(ymin = null_mi_mean - null_mi_sd,
                      ymax = null_mi_mean + null_mi_sd),
                  width = 0.25, colour = PALETTE$null, linewidth = 1.2) +
    geom_text(aes(y     = observed_mi * 1.04,
                  label = ifelse(p_value < 0.001, "p < 0.001",
                                 sprintf("p = %.3f", p_value))),
              size = 3.2, fontface = "italic") +
    facet_wrap(~ agg_method, labeller = as_labeller(pool_labels_used)) +
    scale_fill_manual(values = TARGET_COLOURS) +
    scale_y_continuous(labels = number_format(accuracy = 0.001)) +
    labs(x = NULL, y = "Mutual information (nats)",
         title = sprintf("A   Recoverability by target (fixed-decomposition p-value)  [%s, %s, %d components]",
                         decomp_label, EMBED_TYPE, N_PC)) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 11))

  null_df <- do.call(rbind, lapply(pool_names, function(m) {
    out <- if (m == pn1) out1 else out2
    do.call(rbind, lapply(TARGETS_V, function(tgt) {
      data.frame(agg_method = m,
                 target     = tgt,
                 null_mi    = out$result$details[[tgt]]$MI.null,
                 stringsAsFactors = FALSE)
    }))
  }))
  null_df$target     <- factor(null_df$target,     levels = TARGETS_V)
  null_df$agg_method <- factor(null_df$agg_method, levels = pool_names)

  p_B <- ggplot(null_df, aes(x = null_mi, fill = target)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 45, alpha = 0.50, colour = NA,
                   position = "identity") +
    geom_vline(data = df_summary,
               aes(xintercept = observed_mi, colour = target),
               linewidth = 1.5, linetype = "dashed") +
    facet_wrap(~ agg_method, labeller = as_labeller(pool_labels_used)) +
    scale_fill_manual(values   = TARGET_COLOURS, name = "Null") +
    scale_colour_manual(values = TARGET_COLOURS, name = "Observed") +
    scale_x_continuous(labels = number_format(accuracy = 0.0001)) +
    labs(x = "MI under null (nats)", y = "Density",
         title = sprintf("B   Permutation null distribution (n = %d)", N_PERM)) +
    theme_classic(base_size = 11) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", size = 11))

  # Fig C draws from incremental[[pool]][[target]]. Under
  # per-target selection (R2b), the caller populates only the
  # selection target's entry for each pool (the other target is
  # deliberately absent rather than computed on a decomposition
  # it didn't choose) — TARGETS_V is filtered to whatever keys
  # are actually present so the panel just shows fewer bars rather
  # than erroring.
  roman <- c("i", "ii", "iii")
  df_C  <- do.call(rbind, lapply(pool_names, function(m) {
    present_targets <- intersect(TARGETS_V, names(incremental[[m]]))
    do.call(rbind, lapply(present_targets, function(tgt) {
      v <- incremental[[m]][[tgt]]
      data.frame(agg_method = m,
                 target     = tgt,
                 bar        = roman,
                 value      = c(v$covar_only, v$emb_only, v$combined),
                 delta      = v$incremental_subho_mi,
                 stringsAsFactors = FALSE)
    }))
  }))
  df_C$bar        <- factor(df_C$bar,        levels = roman)
  df_C$target     <- factor(df_C$target,     levels = TARGETS_V)
  df_C$agg_method <- factor(df_C$agg_method, levels = pool_names)

  delta_df <- df_C %>%
    group_by(agg_method, target) %>%
    summarise(delta = first(delta),
              y_top = max(value) + 0.008, .groups = "drop")

  covar_label <- paste(colnames(C), collapse = ", ")

  p_C <- ggplot(df_C, aes(x = bar, y = value, fill = target, alpha = bar)) +
    geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
    geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey50") +
    geom_text(data = delta_df,
              aes(x = 2, y = y_top, colour = target,
                  label = sprintf("Incremental SUBHO MI = %+.3f", delta)),
              vjust = 0, size = 3.0, fontface = "italic",
              inherit.aes = FALSE) +
    facet_grid(agg_method ~ target,
               labeller = labeller(agg_method = pool_labels_used)) +
    scale_fill_manual(values   = TARGET_COLOURS) +
    scale_colour_manual(values = TARGET_COLOURS) +
    scale_alpha_manual(
      values = c(i = 0.38, ii = 1.00, iii = 0.68),
      labels = c("(i) Covariates only", "(ii) Embedding only",
                 "(iii) Embedding + covariates"),
      name   = NULL) +
    scale_y_continuous(labels = number_format(accuracy = 0.001)) +
    labs(x = NULL, y = "Mutual information (nats)",
         title = sprintf(
           "C   Incremental model-based MI after covariate adjustment [%s]\n     Covariates: %s",
           EMBED_TYPE, covar_label)) +
    theme_classic(base_size = 11) +
    theme(plot.title      = element_text(face = "bold", size = 11),
          strip.text      = element_text(face = "bold"),
          legend.position = "right")

  tag    <- gsub("[^A-Za-z0-9_]", "_", tolower(decomp_label))
  out_AB <- file.path(save_dir,
                      sprintf("figures_kpca_prad_%s_%s_AB.png", EMBED_TYPE, tag))
  png(out_AB, width = 14, height = 6, units = "in", res = 300)
  grid.arrange(p_A, p_B, ncol = 2,
               top = grid::textGrob(
                 sprintf("%s — SUBHO Logistic MI | TCGA PRAD (ERG / ETS) | %s",
                         EMBED_TYPE, decomp_label),
                 gp = grid::gpar(fontsize = 13, fontface = "bold")))
  dev.off()
  cat(sprintf("Figures A-B -> %s\n", out_AB))

  out_C <- file.path(save_dir,
                     sprintf("figures_kpca_prad_%s_%s_C.png", EMBED_TYPE, tag))
  ggsave(out_C, plot = p_C, width = 12, height = 8, dpi = 300)
  cat(sprintf("Figure C    -> %s\n", out_C))

  invisible(list(A = p_A, B = p_B, C = p_C))
}

################################################################
# BLOCK 10: Scree comparison — unchanged
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

  out_path <- file.path(save_dir, sprintf("kpca_scree_comparison_prad_%s.png", EMBED_TYPE))
  png(out_path, width = 16, height = 6, units = "in", res = 300)
  grid.arrange(p_ind, p_cum, ncol = 2,
               top = grid::textGrob(
                 sprintf("Eigenspectrum — TCGA PRAD | %s | %d components | 4 pooling methods",
                         EMBED_TYPE, N_PC),
                 gp = grid::gpar(fontsize = 13, fontface = "bold")))
  dev.off()
  cat(sprintf("Scree comparison -> %s\n", out_path))
  invisible(list(ind = p_ind, cum = p_cum))
}

################################################################
# BLOCK 11: MI comparison plot — unchanged
################################################################

make_mi_comparison <- function(results_list, save_dir = RESULTS_DIR) {
  df <- do.call(rbind, lapply(names(results_list), function(nm) {
    r <- results_list[[nm]]
    do.call(rbind, lapply(TARGETS_V, function(tgt) {
      det <- r$result$details[[tgt]]
      s   <- r$result$summary
      data.frame(
        key        = nm,
        label      = r$label,
        decomp     = method_from_key(nm),
        agg_method = pool_from_key(nm),
        target     = tgt,
        MI_obs     = det$MI,
        p_value    = det$p.value,
        p_adj_BH   = s$P_adj_BH[s$Target == tgt],
        null_mean  = mean(det$MI.null, na.rm = TRUE),
        null_sd    = sd(det$MI.null,   na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))

  decomp_order <- c("PCA", "cosine", "rbf", "poly")
  pool_order   <- c("mean", "norm", "moment", "attention")
  df$decomp     <- factor(df$decomp,     levels = decomp_order)
  df$agg_method <- factor(df$agg_method, levels = pool_order)
  df$target     <- factor(df$target,     levels = TARGETS_V)

  pal_decomp  <- c("PCA"    = "#2166AC", "cosine" = "#D6604D",
                   "rbf"    = "#4DAC26", "poly"   = "#8B3A8B")
  pool_labels <- c(mean = "Mean", norm = "Norm",
                   moment = "Moment", attention = "Attention")

  p_mi <- ggplot(df, aes(x = decomp, y = MI_obs, fill = decomp,
                         alpha = target)) +
    geom_bar(stat = "identity", position = position_dodge(0.75),
             width = 0.65, colour = "grey20", linewidth = 0.3) +
    geom_errorbar(
      aes(ymin = null_mean - null_sd, ymax = null_mean + null_sd),
      position = position_dodge(0.75), width = 0.25,
      colour = "grey30", linewidth = 0.5
    ) +
    facet_wrap(~ agg_method, nrow = 1,
               labeller = as_labeller(pool_labels)) +
    scale_fill_manual(values = pal_decomp, name = "Decomposition") +
    scale_alpha_manual(values = c(ERG = 1.0, ETS = 0.5), name = "Target") +
    labs(x = "Decomposition", y = "MI (nats)",
         title = sprintf("Observed MI with permutation-null mean \u00b1 1 SD [%s]  [ERG solid | ETS transparent]",
                         EMBED_TYPE)) +
    theme_classic(base_size = 10) +
    theme(axis.text.x     = element_text(angle = 30, hjust = 1),
          strip.text      = element_text(face = "bold"),
          legend.position = "bottom")

  p_floor <- 1 / (N_PERM + 1)
  p_pv <- ggplot(df, aes(x = decomp,
                         y = -log10(pmax(p_adj_BH, p_floor)),
                         fill = decomp, alpha = target)) +
    geom_bar(stat = "identity", position = position_dodge(0.75),
             width = 0.65, colour = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "firebrick", linewidth = 0.5) +
    facet_wrap(~ agg_method, nrow = 1,
               labeller = as_labeller(pool_labels)) +
    scale_fill_manual(values = pal_decomp, name = "Decomposition") +
    scale_alpha_manual(values = c(ERG = 1.0, ETS = 0.5), name = "Target") +
    labs(x = "Decomposition",
         y = expression(-log[10](p[BH])),
         title = sprintf("BH-adjusted p-value [%s]  [ERG solid | ETS transparent]", EMBED_TYPE)) +
    theme_classic(base_size = 10) +
    theme(axis.text.x     = element_text(angle = 30, hjust = 1),
          strip.text      = element_text(face = "bold"),
          legend.position = "bottom")

  out_path <- file.path(save_dir, sprintf("kpca_mi_comparison_prad_%s.png", EMBED_TYPE))
  png(out_path, width = 18, height = 10, units = "in", res = 300)
  grid.arrange(p_mi, p_pv, nrow = 2,
               top = grid::textGrob(
                 sprintf("MI comparison — 4 pooling x 4 decompositions x 2 targets | TCGA PRAD | %s",
                         EMBED_TYPE),
                 gp = grid::gpar(fontsize = 13, fontface = "bold")))
  dev.off()
  cat(sprintf("MI comparison -> %s\n", out_path))

  cat("\n=== MI Comparison Summary ===\n")
  print(df[order(df$agg_method, df$decomp, df$target),
           c("decomp", "agg_method", "target", "MI_obs", "p_value", "p_adj_BH")],
        row.names = FALSE)

  invisible(df)
}

################################################################
# BLOCK 12: Select best decomposition
#
# R2 fix. Two problems in the original version:
#
#   (a) eligible <- key_summary[key_summary$any_sig & ...] used
#       observed BH significance to decide which candidates were
#       allowed to compete, then took the max observed MI among
#       only those survivors. That filters on the same outcome
#       used to compute the null maximum afterward, so the
#       resulting p_nested was not a valid p-value: significance
#       was baked in before the null was built.
#
#   (b) best_key was chosen by max_MI = max(MI_obs) pooled across
#       BOTH ERG and ETS. A decomposition could be selected and
#       reported for ERG purely because it happened to score high
#       on ETS.
#
# Fix: eligibility is now ONLY the pre-specified eig_cond numeric
# filter (a property of the decomposition itself, fixed before
# looking at any label). Best-key selection and the maxT null are
# both computed per target separately, over the FULL eligible
# candidate set in the pooling group, with no significance
# pre-filtering. This is the per-target, per-pooling-group maxT
# correction: it answers "if I always pick whichever of these
# eligible decompositions gives the largest observed MI for this
# target, how extreme is that pick relative to picking the same
# way under permuted labels."
################################################################

select_best_decomp <- function(results_list, decomp_list,
                               save_dir = RESULTS_DIR) {
  df <- do.call(rbind, lapply(names(results_list), function(nm) {
    r <- results_list[[nm]]
    do.call(rbind, lapply(TARGETS_V, function(tgt) {
      det <- r$result$details[[tgt]]
      s   <- r$result$summary
      data.frame(
        key        = nm,
        label      = r$label,
        agg_method = pool_from_key(nm),
        target     = tgt,
        MI_obs     = det$MI,
        p_value    = det$p.value,
        p_adj_BH   = s$P_adj_BH[s$Target == tgt],
        eig_cond   = decomp_list[[nm]]$eig_cond,
        stringsAsFactors = FALSE
      )
    }))
  }))

  EIG_COND_MAX <- 1e10
  best_list    <- list()

  for (pool in c("mean", "norm", "moment", "attention")) {
    sub_pool <- df[df$agg_method == pool, ]
    if (nrow(sub_pool) == 0L) next

    # Eligibility is purely the eig_cond numeric filter, fixed
    # before looking at any label. No any_sig term here (R2a).
    elig_keys <- unique(sub_pool$key[sub_pool$eig_cond < EIG_COND_MAX])
    if (length(elig_keys) == 0L) {
      cat(sprintf("  [%s] No method passed eig_cond < %.0e; falling back to all (exploratory).\n",
                  pool, EIG_COND_MAX))
      elig_keys <- unique(sub_pool$key)
    }

    for (tgt in TARGETS_V) {
      sub_tgt <- sub_pool[sub_pool$key %in% elig_keys & sub_pool$target == tgt, ]

      best_row <- sub_tgt[which.max(sub_tgt$MI_obs), ]
      best_key <- best_row$key

      # maxT null over the SAME eligible candidate set used for
      # selection, unconditional on observed significance (R2b).
      null_vecs <- lapply(results_list[elig_keys], function(r)
        r$result$details[[tgt]]$MI.null)
      B        <- min(sapply(null_vecs, length))
      null_mat <- do.call(cbind, lapply(null_vecs, `[`, 1:B))
      T_null   <- apply(null_mat, 1, max, na.rm = TRUE)
      T_obs    <- best_row$MI_obs
      n_valid  <- sum(!is.na(T_null))
      p_maxT   <- (1 + sum(T_null >= T_obs, na.rm = TRUE)) / (1 + n_valid)

      cat(sprintf(
        "  [%-9s | %-3s] Best: %-32s | MI=%.5f | p_maxT=%.4f (n_candidates=%d)\n",
        pool, tgt, best_row$label, T_obs, p_maxT, length(elig_keys)))

      best_list[[paste(pool, tgt, sep = "__")]] <- data.frame(
        key        = best_key,
        label      = best_row$label,
        agg_method = pool,
        target     = tgt,
        MI_obs     = T_obs,
        eig_cond   = best_row$eig_cond,
        p_maxT     = p_maxT,
        n_candidates = length(elig_keys),
        stringsAsFactors = FALSE
      )
    }
  }

  best_df <- do.call(rbind, best_list)
  write.csv(best_df,
            file.path(save_dir, sprintf("kpca_selected_decomposition_prad_%s.csv", EMBED_TYPE)),
            row.names = FALSE)
  cat(sprintf("  Selection written to kpca_selected_decomposition_prad_%s.csv\n", EMBED_TYPE))
  cat("  p_maxT is the selection-adjusted p-value for the per-target,\n")
  cat("  per-pooling-group decomposition pick. Fixed-decomposition\n")
  cat("  p-values reported elsewhere (run_subho summaries, Figure A)\n")
  cat("  are NOT valid for whichever decomposition ends up selected here.\n")
  best_df
}

################################################################
# BLOCK 13: Main execution
################################################################

set.seed(RANDOM_SEED)

cat("\n", strrep("=", 70), "\n")
cat(sprintf("PRAD Trident embeddings | %s | Kernel PCA comparison | SUBHO Logistic MI\n", EMBED_TYPE))
cat("Targets: ERG_status (fusion) | ETS (Positive)\n")
cat(strrep("=", 70), "\n\n")

# 13.0  Index + key sanity check before touching the full roster
cat("=== Indexing embedding directory and checking h5 keys ===\n")
stem_map <- index_embed_dir(EMBED_DIR)
discover_h5_keys(EMBED_DIR, stem_map)
diagnose_slide_matching(META_PATH, stem_map)

# 13.1  Unified patient roster
cat("\n=== Building unified patient roster ===\n")
roster   <- build_unified_roster(META_PATH, EMBED_DIR, stem_map)
y_list      <- roster$label_list
tile_counts <- roster$tile_counts

cat(sprintf("\nFinal cohort: n = %d patients\n", length(tile_counts)))

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

for (tgt in TARGETS_V) {
  names(y_list[[tgt]]) <- names(roster$label_list[[tgt]])
  y_list[[tgt]]        <- y_list[[tgt]][rownames(Z_raw_mean)]
}

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

# 13.5  SUBHO on all 16 decompositions x 2 targets
cat(sprintf("\n=== Running SUBHO on all 16 decompositions (2 targets each) [%s] ===\n", EMBED_TYPE))
results <- list()
for (nm in names(decomp)) {
  pool  <- pool_from_key(nm)
  y_fit <- lapply(y_list, function(yv) yv[rownames(decomp[[nm]]$Z_scaled)])
  results[[nm]] <- run_subho(
    Z_scaled   = decomp[[nm]]$Z_scaled,
    y_list     = y_fit,
    agg_method = pool,
    label      = label_map[nm]
  )
}

saveRDS(results, file.path(RESULTS_DIR, sprintf("kpca_subho_results_prad_%s.rds", EMBED_TYPE)))
cat(sprintf("Results saved -> %s\n",
            file.path(RESULTS_DIR, sprintf("kpca_subho_results_prad_%s.rds", EMBED_TYPE))))

# 13.6  MI comparison
mi_df <- make_mi_comparison(results)
write.csv(mi_df, file.path(RESULTS_DIR, sprintf("kpca_mi_summary_prad_%s.csv", EMBED_TYPE)),
          row.names = FALSE)

# 13.7  Select best decomposition per pooling method, per target (R2)
cat("\n=== Selecting best decomposition per pooling method, per target ===\n")
best_df <- select_best_decomp(results, decomp)
write.csv(best_df, file.path(RESULTS_DIR, "best_kpca_mi_summary_prad.csv"),
          row.names = FALSE)

# best_df now has one row per (pooling method, target) combination —
# 8 rows for 4 pools x 2 targets — because selection is no longer
# pooled across ERG/ETS (R2b).

# 13.8  Incremental MI + figures, run once per target's own selection
cat("\n=== Downstream analysis on selected decompositions (per target) ===\n")

any_key <- best_df$key[1]
C <- build_covariate_matrix(META_PATH,
                            rownames(results[[any_key]]$Z_scaled),
                            tile_counts)

POOL_ORDER <- c("mean", "norm", "moment", "attention")

# all_incremental[[selection_target]][[pool]] holds the incremental
# MI computed on the decomposition that selection_target picked for
# that pool. Both ERG and ETS incremental MI are computed inside
# each pass (compute_incremental_mi_one_target is called for the
# selection target only — the non-selection target is NOT computed
# on a decomposition it didn't choose, to avoid implying it was
# evaluated on its own best representation).
all_incremental <- list()

for (sel_tgt in TARGETS_V) {
  all_incremental[[sel_tgt]] <- list()

  for (pool in POOL_ORDER) {
    row      <- best_df[best_df$agg_method == pool & best_df$target == sel_tgt, ]
    if (nrow(row) == 0L) next
    best_key <- row$key
    out_best <- results[[best_key]]
    y_fit    <- y_list[[sel_tgt]][rownames(out_best$Z_scaled)]

    cat(sprintf("\n  [selection target = %s | pool = %s] Best decomp: %s\n",
                sel_tgt, pool, row$label))
    inc <- compute_incremental_mi_one_target(out_best$Z_scaled, C, y_fit, sel_tgt)
    all_incremental[[sel_tgt]][[pool]] <- setNames(list(inc), sel_tgt)

    write.csv(
      cbind(selection_target = sel_tgt, agg_method = pool, as.data.frame(inc)),
      file.path(RESULTS_DIR,
                sprintf("incremental_mi_prad_%s_%s_%s.csv", EMBED_TYPE, pool, sel_tgt)),
      row.names = FALSE)
  }
}

# Figures: two passes (ERG-selected decomps, ETS-selected decomps),
# each pass producing the same mean+norm / moment+attention figure
# pairs as before. Within a pass, both targets still appear as
# bars/fill inside each panel (make_figures shows result$details for
# both ERG and ETS regardless of which one drove selection) — the
# figure title flags which target's selection produced the decomp,
# so the non-selection target's numbers in that pass are read as
# "evaluated on a decomposition chosen for the other target,"
# not as that target's own best result.
for (sel_tgt in TARGETS_V) {

  key_for <- function(pool) best_df$key[best_df$agg_method == pool & best_df$target == sel_tgt]

  make_figures(results[[key_for("mean")]],
               results[[key_for("norm")]],
               list(mean = all_incremental[[sel_tgt]][["mean"]],
                    norm = all_incremental[[sel_tgt]][["norm"]]),
               C,
               pool_names   = c("mean", "norm"),
               decomp_label = sprintf("Best kPCA selected for %s [mean & norm pooling]", sel_tgt))

  make_figures(results[[key_for("moment")]],
               results[[key_for("attention")]],
               list(moment    = all_incremental[[sel_tgt]][["moment"]],
                    attention = all_incremental[[sel_tgt]][["attention"]]),
               C,
               pool_names   = c("moment", "attention"),
               decomp_label = sprintf("Best kPCA selected for %s [moment & attention pooling]", sel_tgt))
}

# Final summary
cat(sprintf("\n%s\n", strrep("=", 70)))
cat(sprintf("PRAD Trident (%s) kernel PCA pipeline complete.\n", EMBED_TYPE))
cat("Selection is per-target; p_maxT is the selection-adjusted p-value\n")
cat("for each (pool, target) pick. Fixed-decomposition p-values from\n")
cat("run_subho()/Figure A are valid only before this selection step.\n\n")
for (sel_tgt in TARGETS_V) {
  cat(sprintf("--- Selected for %s ---\n", sel_tgt))
  for (pool in POOL_ORDER) {
    row <- best_df[best_df$agg_method == pool & best_df$target == sel_tgt, ]
    if (nrow(row) == 0L) next
    cat(sprintf("  [%-9s]: %-35s | MI=%.5f | p_maxT=%.4f\n",
                pool, row$label, row$MI_obs, row$p_maxT))
  }
}
cat(sprintf("%s\n", strrep("=", 70)))


################################################################
# prad_incremental_batch_purity.R
#
# Incremental recoverability for PRAD with scanning batch and tumour
# purity added to the covariate set, for all four foundation models.
#
# The original pipeline (kPCA_SUBHO_PRAD_trident*.R, block 8) used
# C0 = {standardised log tile count, standardised reviewed Gleason
# sum}; the covariate-only estimate was exactly zero in all 32 arms.
# Here the same estimator, settings and seed are run with
#   C_p   = C0 + ABSOLUTE tumour purity (median-imputed; 36 of 300
#           annotation rows lack it),
#   C_b   = C0 + scanning batch (indicator per batch with >= 8
#           patients, smaller batches pooled as "other", reference =
#           the largest batch, columns standardised like every other
#           covariate),
#   C_all = C0 + purity + batch,
# as covariate-only models (Model 1, C -> Y) and, for C_all, combined
# with the selected embedding of every (pooling, target) arm
# (Model 2, (C, Z) -> Y). Embedding-only and combined-with-C0 values
# are taken from the existing incremental CSVs so that nothing that
# was reported changes.
#
# Tile counts per slide and model are read from the h5 feature shapes
# (they were not saved by the pipeline). Output:
#   TRIDENT/prad_incremental_batch_purity/prad_incremental_batch_purity.csv
#   TRIDENT/prad_incremental_batch_purity/covariate_only_fits.csv
#   ... and a log. Copies into PRAD/ afterwards.
# Usage (from the project root, R 4.4.3):
#   Rscript TRIDENT/prad_incremental_batch_purity.R
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(reticulate) })
reticulate::py_require(c("h5py", "numpy")); py_run_string("import h5py, os")

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI  <- normalizePath(dirname(this_file), winslash = "/")
ROOT <- normalizePath(file.path(TRI, ".."), winslash = "/")
setwd(ROOT); source("PathMI_SUBHO_logistic.R")
OUT <- file.path(TRI, "prad_incremental_batch_purity"); dir.create(OUT, showWarnings = FALSE)
LOG <- file.path(OUT, "run.log")
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }

V_BETA <- 1.0; NMC <- 2000L; BURN <- 500L; THIN <- 10L; RANDOM_SEED <- 42L
META <- "M:/shared_data/TCGA_PRAD/orig_annotations.csv"
FEATS <- "M:/shared_data/TCGA_PRAD/trident-feats"
BB <- list(uni_v2    = list(dir = file.path(FEATS, "20x_256px_0px_overlap/features_uni_v2"),    label = "UNI v2"),
           virchow2  = list(dir = file.path(FEATS, "20x_224px_0px_overlap/features_virchow2"),  label = "Virchow2"),
           hoptimus1 = list(dir = file.path(FEATS, "20x_224px_0px_overlap/features_hoptimus1"), label = "H-optimus-1"),
           conch_v15 = list(dir = file.path(FEATS, "20x_512px_0px_overlap/features_conch_v15"), label = "CONCH v1.5"))
fit_mi <- function(X, y) { set.seed(RANDOM_SEED)
  r <- subho_MI_logistic(X, y, v_beta = V_BETA, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
  list(MI = r$MI, raw = r$MI.raw, k = r$n.selected, sel = paste(colnames(X)[r$selected.features], collapse = ";")) }

tile_counts <- function(d) {
  py_run_string(sprintf("
_tc = {}
for _f in os.listdir(r'%s'):
    if _f.endswith('.h5'):
        with h5py.File(os.path.join(r'%s', _f), 'r') as _h:
            _tc[_h['features'].attrs.get('name', _f)] = int(_h['features'].shape[0])
", d, d)); unlist(py$`_tc`) }

meta <- read.csv(META, stringsAsFactors = FALSE)
meta <- meta[!is.na(meta$slide) & !duplicated(meta$slide), ]; rownames(meta) <- meta$slide
t0 <- Sys.time(); logf("=== PRAD incremental recoverability with batch and purity ===")
cov_rows <- list(); arm_rows <- list()

for (bb in names(BB)) {
  res  <- readRDS(file.path(TRI, sprintf("TCGA_PRAD_kpca_subho_%s/kpca_subho_results_prad_%s.rds", bb, bb)))
  best <- read.csv(file.path(TRI, sprintf("TCGA_PRAD_kpca_subho_%s/best_kpca_mi_summary_prad.csv", bb)), stringsAsFactors = FALSE)
  roster <- rownames(res$pca_mean$Z_scaled); y_list <- lapply(res$pca_mean$y_list, function(v) v[roster])
  tc <- tile_counts(BB[[bb]]$dir); stopifnot(all(roster %in% names(tc)))
  m <- meta[roster, ]
  # covariates, all standardised (unit variance) as in the pipeline
  gle <- m$Reviewed_Gleason_sum; gle[is.na(gle)] <- median(gle, na.rm = TRUE)
  pur <- suppressWarnings(as.numeric(m$ABSOLUTE_purity)); n_pur_na <- sum(is.na(pur)); pur[is.na(pur)] <- median(pur, na.rm = TRUE)
  bt  <- as.character(m$Batch); keep <- names(which(table(bt) >= 8)); bt[!bt %in% keep] <- "other"
  ref <- names(which.max(table(bt))); lv <- setdiff(sort(unique(bt)), ref)
  Bm  <- sapply(lv, function(l) as.numeric(bt == l)); colnames(Bm) <- paste0("batch_", lv)
  C0  <- cbind(log_tile_count_z = as.numeric(scale(log1p(tc[roster]))), Gleason_sum_z = as.numeric(scale(gle)))
  Cp  <- cbind(C0, purity_z = as.numeric(scale(pur)))
  Cb  <- cbind(C0, scale(Bm)); Call <- cbind(Cp, scale(Bm))
  rownames(C0) <- rownames(Cp) <- rownames(Cb) <- rownames(Call) <- roster
  logf("[%s] n = %d | purity imputed for %d | batch levels kept %d (+other), reference %s | C_all has %d columns",
       bb, length(roster), n_pur_na, length(lv), ref, ncol(Call))
  sets <- list(C0 = C0, C_purity = Cp, C_batch = Cb, C_all = Call)
  cov_mi <- list()
  for (tg in c("ERG", "ETS")) for (s in names(sets)) {
    f <- fit_mi(sets[[s]], as.numeric(y_list[[tg]])); cov_mi[[paste(tg, s)]] <- f$MI
    cov_rows[[length(cov_rows) + 1]] <- data.frame(model = BB[[bb]]$label, target = tg, covariate_set = s, n_cov = ncol(sets[[s]]),
                                                  MI = f$MI, MI_raw = f$raw, k_star = f$k, selected = f$sel)
    logf("  [%s %s] %-9s MI = %.4f (raw %.4f) k* = %d sel = %s", bb, tg, s, f$MI, f$raw, f$k, f$sel)
  }
  for (i in seq_len(nrow(best))) {
    tg <- best$target[i]; pool <- best$agg_method[i]; key <- best$key[i]
    Z <- res[[key]]$Z_scaled[roster, , drop = FALSE]; y <- as.numeric(y_list[[tg]])
    old <- read.csv(file.path(TRI, sprintf("TCGA_PRAD_kpca_subho_%s/incremental_mi_prad_%s_%s_%s.csv", bb, bb, pool, tg)))
    f <- fit_mi(cbind(Z, Call), y)
    arm_rows[[length(arm_rows) + 1]] <- data.frame(model = BB[[bb]]$label, target = tg, pooling = pool, key = key,
      covar_C0_reported = old$covar_only, covar_C0_refit = cov_mi[[paste(tg, "C0")]], covar_purity = cov_mi[[paste(tg, "C_purity")]],
      covar_batch = cov_mi[[paste(tg, "C_batch")]], covar_all = cov_mi[[paste(tg, "C_all")]],
      emb_only = old$emb_only, combined_C0 = old$combined, combined_all = f$MI, combined_all_raw = f$raw,
      incremental_all = f$MI - cov_mi[[paste(tg, "C_all")]], k_star_combined = f$k, selected_combined = f$sel)
    logf("  [%s %s %s] combined(C_all) = %.4f (k* = %d; %s) | emb_only %.4f | combined(C0) %.4f | %.1f min",
         bb, tg, pool, f$MI, f$k, f$sel, old$emb_only, old$combined, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
}
cov_df <- do.call(rbind, cov_rows); arm_df <- do.call(rbind, arm_rows)
write.csv(cov_df, file.path(OUT, "covariate_only_fits.csv"), row.names = FALSE)
write.csv(arm_df, file.path(OUT, "prad_incremental_batch_purity.csv"), row.names = FALSE)
dst <- file.path(ROOT, "PRAD", "prad_incremental_batch_purity"); dir.create(dst, showWarnings = FALSE)
file.copy(list.files(OUT, full.names = TRUE), dst, overwrite = TRUE)
logf("done in %.1f min; %d covariate-only fits, %d combined fits; copied to PRAD/", as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(cov_df), nrow(arm_df))
print(cov_df[, c("model", "target", "covariate_set", "MI", "k_star", "selected")], row.names = FALSE)
print(arm_df[, c("model", "target", "pooling", "covar_all", "emb_only", "combined_C0", "combined_all", "incremental_all")], digits = 4, row.names = FALSE)

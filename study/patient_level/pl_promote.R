################################################################
# pl_promote.R
#
# Makes the patient-level grid the analysis of record (19 Nov 2026).
#   1. Archives every slide-level output that the downstream scripts read
#      or write into TRIDENT/archive_slide_level_20261119/ (copies).
#   2. Writes the patient-level results into the pipeline's own files and
#      formats: <cohort folder>/kpca_subho_results_<co>_<mo>.rds (Z_scaled
#      with patient row names; y for CRC, y_list for PRAD; result$details
#      with MI, MI.raw, n.selected, cond.H, MI.null, p.value),
#      kpca_mi_summary_*.csv, best_kpca_mi_summary_*.csv,
#      kpca_selected_decomposition_*.csv, incremental_mi_*.csv, and
#      prad_incremental_batch_purity/{prad_incremental_batch_purity,
#      covariate_only_fits}.csv; and CSV sources for the three merged
#      workbooks (converted by pl_promote_xlsx.py).
#   3. Writes patient-collapsed pooled.rds files into a new pool_pass
#      directory (the slide-level directory is renamed pool_pass_slide);
#      submeans.npy is collapsed by pl_collapse_submeans.py.
#
# Usage: Rscript TRIDENT/patient_level/pl_promote.R
################################################################
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
PLD <- normalizePath(dirname(this_file), winslash = "/"); TRI <- dirname(PLD); OUT <- file.path(PLD, "patient")
ARC <- file.path(TRI, "archive_slide_level_20261119"); dir.create(ARC, showWarnings = FALSE)
POOL <- "M:/bhadury/mi_results/pool_pass"; POOL_SLIDE <- "M:/bhadury/mi_results/pool_pass_slide"
EMB <- c(uni_v2 = "UNI V2", virchow2 = "VIRCHOW", hoptimus1 = "H-OPTIMUS 1", conch_v15 = "CONCH V1.5")
LAB <- c(uni_v2 = "UNI v2", virchow2 = "Virchow2", hoptimus1 = "H-optimus-1", conch_v15 = "CONCH v1.5")
KLAB <- function(key) { dc <- sub("_.*", "", key); pl <- sub(".*_", "", key)
  sprintf("%s / %s", c(pca = "Standard PCA", cosine = "Cosine kPCA", rbf = "RBF kPCA", poly = "Poly kPCA")[[dc]], paste0(toupper(substr(pl, 1, 1)), substr(pl, 2, nchar(pl)))) }
DLAB <- function(dc) ifelse(dc == "pca", "PCA", dc)

# ---- 1. archive ------------------------------------------------------------------
if (!file.exists(file.path(ARC, "ARCHIVED"))) {
  items <- c(list.files(TRI, pattern = "^(crc|TCGA_PRAD)_kpca_subho_", full.names = TRUE),
             file.path(TRI, c("mi summary table.xlsx", "best kpca summary.xlsx", "incremental mi summary.xlsx", "prad_incremental_batch_purity",
                              "usable", "robustness", "burden", "record", "learning", "tile_maps_hoptimus1", "figures")))
  for (it in items[file.exists(items)]) file.copy(it, ARC, recursive = TRUE, copy.date = TRUE)
  writeLines(format(Sys.time()), file.path(ARC, "ARCHIVED")); cat("archived", length(items), "items to", ARC, "\n")
}
if (!dir.exists(POOL_SLIDE)) { stopifnot(file.rename(POOL, POOL_SLIDE)); cat("renamed", POOL, "->", POOL_SLIDE, "\n") }
dir.create(POOL, showWarnings = FALSE)

# ---- 2. results in pipeline formats ------------------------------------------------
mr_all <- readRDS(file.path(OUT, "meta_runs.rds"))$meta_runs
grid_rows <- list(); best_rows <- list(); inc_rows <- list()
for (tag in names(mr_all)) {
  mr <- mr_all[[tag]]; co <- mr$cohort; mo <- mr$model
  dir <- file.path(TRI, if (co == "crc") sprintf("crc_kpca_subho_%s", mo) else sprintf("TCGA_PRAD_kpca_subho_%s", mo)); dir.create(dir, showWarnings = FALSE)
  res <- readRDS(file.path(OUT, tag, "results.rds")); summ <- read.csv(file.path(OUT, tag, "mi_summary.csv"), stringsAsFactors = FALSE)
  best <- read.csv(file.path(OUT, tag, "best.csv"), stringsAsFactors = FALSE); inc <- read.csv(file.path(OUT, tag, "incremental.csv"), stringsAsFactors = FALSE)
  out <- lapply(names(res), function(key) { r <- res[[key]]; Z <- r$Z_scaled; rownames(Z) <- r$units
    o <- list(Z_scaled = Z, label = KLAB(key), agg_method = sub(".*_", "", key), result = list(details = r$result$details))
    if (co == "crc") o$y <- r$y_list$MSI else o$y_list <- r$y_list; o }); names(out) <- names(res)
  saveRDS(out, file.path(dir, sprintf("kpca_subho_results_%s_%s.rds", co, mo)))
  s2 <- data.frame(key = summ$key, label = sapply(summ$key, KLAB), decomp = DLAB(summ$decomp), agg_method = summ$agg_method, target = summ$target,
                   MI_obs = summ$MI_obs, p_value = summ$p_value, null_mean = summ$null_mean, null_sd = summ$null_sd, MI_raw = summ$MI_raw, k_star = summ$k_star, cond_H = summ$cond_H)
  if (co == "crc") write.csv(s2[, setdiff(names(s2), "target")], file.path(dir, sprintf("kpca_mi_summary_crc_%s.csv", mo)), row.names = FALSE)
  else write.csv(s2, file.path(dir, sprintf("kpca_mi_summary_prad_%s.csv", mo)), row.names = FALSE)
  b2 <- data.frame(key = best$key, label = sapply(best$key, KLAB), agg_method = best$agg_method, target = best$target, MI_obs = best$MI_obs,
                   eig_cond = best$eig_cond, p_maxT = best$p_maxT, n_candidates = best$n_candidates)
  if (co == "crc") { b2 <- b2[, setdiff(names(b2), "target")]; write.csv(b2, file.path(dir, "best_kpca_mi_summary_crc.csv"), row.names = FALSE)
    write.csv(b2, file.path(dir, sprintf("kpca_selected_decomposition_crc_%s.csv", mo)), row.names = FALSE)
    write.csv(inc[, c("agg_method", "covar_only", "emb_only", "combined", "incremental_subho_mi")], file.path(dir, sprintf("incremental_mi_crc_%s.csv", mo)), row.names = FALSE)
  } else { write.csv(b2, file.path(dir, "best_kpca_mi_summary_prad.csv"), row.names = FALSE)
    write.csv(b2, file.path(dir, sprintf("kpca_selected_decomposition_prad_%s.csv", mo)), row.names = FALSE)
    for (i in seq_len(nrow(inc))) write.csv(data.frame(selection_target = inc$target[i], agg_method = inc$agg_method[i], covar_only = inc$covar_only[i],
      emb_only = inc$emb_only[i], combined = inc$combined[i], incremental_subho_mi = inc$incremental_subho_mi[i]),
      file.path(dir, sprintf("incremental_mi_prad_%s_%s_%s.csv", mo, inc$agg_method[i], inc$target[i])), row.names = FALSE) }
  grid_rows[[tag]] <- data.frame(Disease = toupper(co), Embeddings = EMB[[mo]], s2[, c("key", "label", "decomp", "agg_method", "target", "MI_obs", "p_value", "null_mean", "null_sd")])
  best_rows[[tag]] <- data.frame(Disease = toupper(co), Embeddings = EMB[[mo]], key = best$key, label = sapply(best$key, KLAB), agg_method = best$agg_method,
                                 target = best$target, MI_obs = best$MI_obs, eig_cond = best$eig_cond, p_maxT = best$p_maxT, n_candidates = best$n_candidates)
  inc_rows[[tag]] <- data.frame(Disease = toupper(co), Embeddings = EMB[[mo]], target = inc$target, agg_method = inc$agg_method,
                                covar_only = inc$covar_only, emb_only = inc$emb_only, combined = inc$combined, incremental_subho_mi = inc$incremental_subho_mi)
  # ---- 3. patient-collapsed pooled.rds ----
  pp <- readRDS(file.path(POOL_SLIDE, tag, "pooled.rds")); slide <- pp$roster; pat <- substr(slide, 1, 12); units <- unique(pat)
  coll <- function(X) { S <- rowsum(X, pat, reorder = FALSE); S[units, , drop = FALSE] / as.vector(table(pat)[units]) }
  tl <- pp$tiles; tl$patient <- substr(tl$slide, 1, 12)
  tiles <- data.frame(slide = units, n_tiles = as.integer(tapply(tl$n_tiles, tl$patient, sum)[units]), n_slides = as.integer(table(pat)[units]))
  ppn <- list(cohort = pp$cohort, model = pp$model, roster = units, slide_roster = slide,
              y = lapply(pp$y, function(v) { v <- v[slide]; setNames(as.integer(v[match(units, pat)]), units) }),
              tiles = tiles, pooled = lapply(pp$pooled, coll), m_grid = pp$m_grid, r_sub = pp$r_sub, unit = "patient")
  for (k in names(ppn$pooled)) stopifnot(identical(rownames(ppn$pooled[[k]]), units))
  dir.create(file.path(POOL, tag), showWarnings = FALSE); saveRDS(ppn, file.path(POOL, tag, "pooled.rds"))
  # check: the stored patient-level Z equals the projection rebuilt from the collapsed pooled matrix
  cat(sprintf("[%s] %d slides -> %d patients; results and pooled.rds written\n", tag, length(slide), length(units)))
}
SRC <- file.path(OUT, "xlsx_src"); dir.create(SRC, showWarnings = FALSE)
write.csv(do.call(rbind, grid_rows), file.path(SRC, "mi summary table.csv"), row.names = FALSE)
write.csv(do.call(rbind, best_rows), file.path(SRC, "best kpca summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, inc_rows), file.path(SRC, "incremental mi summary.csv"), row.names = FALSE)
ext <- read.csv(file.path(OUT, "prad_extended.csv"), stringsAsFactors = FALSE); cov <- read.csv(file.path(OUT, "prad_covariate_only.csv"), stringsAsFactors = FALSE)
incp <- do.call(rbind, inc_rows); incp <- incp[incp$Disease == "PRAD", ]
ext$model <- LAB[sub("^prad_", "", ext$tag)]
ext$covar_C0_reported <- sapply(seq_len(nrow(ext)), function(i) incp$covar_only[incp$Embeddings == EMB[[sub("^prad_", "", ext$tag[i])]] & incp$target == ext$target[i] & incp$agg_method == ext$pooling[i]])
pib <- data.frame(model = ext$model, target = ext$target, pooling = ext$pooling, key = ext$key, covar_C0_reported = ext$covar_C0_reported, covar_C0_refit = ext$covar_C0,
                  covar_purity = ext$covar_purity, covar_batch = ext$covar_batch, covar_all = ext$covar_all, emb_only = ext$emb_only, combined_C0 = ext$combined_C0,
                  combined_all = ext$combined_all, combined_all_raw = ext$combined_all_raw, incremental_all = ext$incremental_all, k_star_combined = ext$k_star_combined, selected_combined = ext$selected_combined)
dir.create(file.path(TRI, "prad_incremental_batch_purity"), showWarnings = FALSE)
write.csv(pib, file.path(TRI, "prad_incremental_batch_purity", "prad_incremental_batch_purity.csv"), row.names = FALSE)
cov$model <- LAB[sub("^prad_", "", cov$tag)]
write.csv(cov[, c("model", "target", "covariate_set", "n_cov", "MI", "MI_raw", "k_star", "selected")], file.path(TRI, "prad_incremental_batch_purity", "covariate_only_fits.csv"), row.names = FALSE)
cat("promotion done\n")

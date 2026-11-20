################################################################
# pl_incremental.R
#
# Incremental recoverability on the selected arms of the patient-level
# grid (pl_grid.R), in the formats of the original pipelines:
#   CRC  : C0 = log tile count + AJCC stage; covariates only, embedding
#          only and combined, per selected arm (one per pooling rule).
#   PRAD : the same with C0 = log tile count + Gleason sum, per target
#          and pooling rule; and the extended set C_all = C0 + purity +
#          batch (covariates only per target, combined per arm), as in
#          prad_incremental_batch_purity.R.
# Seeds: the base block (covariates, embedding, combined) runs on one
# stream started with set.seed(4242) per cohort x model, so the
# embedding-only refit is a replicate refit with an independent stream
# (the pipelines did not reset the seed here); the extended fits reset
# set.seed(42) before every fit, as prad_incremental_batch_purity.R did.
#
# Usage: Rscript TRIDENT/patient_level/pl_incremental.R [patient|slide]
# Output: <mode>/<tag>/incremental.csv, <mode>/prad_extended.csv, <mode>/prad_covariate_only.csv
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(parallel) })
args <- commandArgs(trailingOnly = TRUE); MODE <- if (length(args)) args[1] else "patient"
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
PLD <- normalizePath(dirname(this_file), winslash = "/"); TRI <- dirname(PLD); ROOT <- dirname(TRI)
setwd(ROOT); source("PathMI_SUBHO_logistic.R")
OUT <- file.path(PLD, MODE); LOG <- file.path(OUT, "pl_incremental.log")
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }
mr_all <- readRDS(file.path(OUT, "meta_runs.rds"))$meta_runs
NMC <- 2000L; BURN <- 500L; THIN <- 10L

one_tag <- function(tag) {
  mr <- mr_all[[tag]]; res <- readRDS(file.path(OUT, tag, "results.rds")); best <- read.csv(file.path(OUT, tag, "best.csv"), stringsAsFactors = FALSE)
  fit <- function(X, y) { r <- subho_MI_logistic(X, y, v_beta = 1, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
    list(MI = r$MI, raw = r$MI.raw, k = r$n.selected, sel = paste(colnames(X)[r$selected.features], collapse = ";")) }
  C0 <- mr$Csets$C0; inc <- list(); ext <- list(); cov <- list()
  set.seed(4242L)
  for (i in seq_len(nrow(best))) { b <- best[i, ]; y <- mr$ylist[[b$target]]; Z <- res[[b$key]]$Z_scaled; colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    fc <- fit(C0, y); fe <- fit(Z, y); fb <- fit(cbind(Z, C0), y)
    inc[[i]] <- data.frame(tag = tag, target = b$target, agg_method = b$agg_method, key = b$key, covar_only = fc$MI, emb_only = fe$MI, combined = fb$MI,
                           incremental_subho_mi = fb$MI - fc$MI, emb_only_raw = fe$raw, selected_MI_obs = b$MI_obs) }
  if (!is.null(mr$Csets$C_all)) {
    Ca <- mr$Csets$C_all; sets <- list(C0 = C0, C_purity = Ca[, c(colnames(C0), "purity_z")], C_batch = Ca[, setdiff(colnames(Ca), "purity_z")], C_all = Ca)
    cm <- list()
    for (tg in names(mr$ylist)) for (s in names(sets)) { set.seed(42L); f <- fit(sets[[s]], mr$ylist[[tg]]); cm[[paste(tg, s)]] <- f$MI
      cov[[length(cov) + 1]] <- data.frame(tag = tag, target = tg, covariate_set = s, n_cov = ncol(sets[[s]]), MI = f$MI, MI_raw = f$raw, k_star = f$k, selected = f$sel) }
    for (i in seq_len(nrow(best))) { b <- best[i, ]; y <- mr$ylist[[b$target]]; Z <- res[[b$key]]$Z_scaled; colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
      set.seed(42L); f <- fit(cbind(Z, Ca), y); bi <- inc[[i]]
      ext[[i]] <- data.frame(tag = tag, target = b$target, pooling = b$agg_method, key = b$key, covar_C0 = cm[[paste(b$target, "C0")]],
        covar_purity = cm[[paste(b$target, "C_purity")]], covar_batch = cm[[paste(b$target, "C_batch")]], covar_all = cm[[paste(b$target, "C_all")]],
        emb_only = bi$emb_only, combined_C0 = bi$combined, combined_all = f$MI, combined_all_raw = f$raw, incremental_all = f$MI - cm[[paste(b$target, "C_all")]],
        k_star_combined = f$k, selected_combined = f$sel) } }
  list(inc = do.call(rbind, inc), ext = if (length(ext)) do.call(rbind, ext) else NULL, cov = if (length(cov)) do.call(rbind, cov) else NULL)
}
cl <- makeCluster(length(mr_all)); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(MASS); library(BayesLogit) }))
clusterExport(cl, c("mr_all", "OUT", "NMC", "BURN", "THIN", "sagnik.horseshoe.logistic", "lpy.x.logistic", "selection.logistic", "subho_MI_logistic", "sigmoid", "clip_eta"))
t0 <- Sys.time(); o <- parLapply(cl, names(mr_all), one_tag); names(o) <- names(mr_all)
for (tag in names(o)) { write.csv(o[[tag]]$inc, file.path(OUT, tag, "incremental.csv"), row.names = FALSE)
  for (i in seq_len(nrow(o[[tag]]$inc))) with(o[[tag]]$inc[i, ], logf("[%s %s %s] covar %.4f | emb %.4f (selected %.4f, diff %+.4f) | combined %.4f | incremental %.4f", tag, target, agg_method, covar_only, emb_only, selected_MI_obs, emb_only - selected_MI_obs, combined, incremental_subho_mi)) }
ext <- do.call(rbind, lapply(o, `[[`, "ext")); cov <- do.call(rbind, lapply(o, `[[`, "cov"))
if (!is.null(ext)) write.csv(ext, file.path(OUT, "prad_extended.csv"), row.names = FALSE)
if (!is.null(cov)) write.csv(cov, file.path(OUT, "prad_covariate_only.csv"), row.names = FALSE)
logf("=== incremental done in %.1f min ===", as.numeric(difftime(Sys.time(), t0, units = "mins")))

################################################################
# usable_recoverability.R
#
# Cross-fitted usable recoverability (out-of-patient), all 192 cells.
#
# PathMI measures the information present in a representation. This
# script measures, for the same cells, how much of it a readout
# fitted on other patients extracts in a new patient:
#
#   I_usable^(F) = (1/n) sum_i log p_{-i}(y_i | z_i) - log p_{-i}(y_i)
#
# under K-fold patient-grouped cross-fitting repeated R times, where
# p_{-i} is fitted on the folds that do not contain patient i and
# p_{-i}(y) is the training-fold prevalence. Two readout families F:
#
#   pathmi  the PathMI readout itself: horseshoe ranking, marginal-
#           likelihood walk and the MAP logistic model on the selected
#           coordinates, all refitted inside the training fold; when
#           the walk's best model has lower evidence than the null
#           (the floored case) the readout returns the prevalence;
#   ridge   ridge logistic regression on all 50 coordinates, lambda by
#           10-fold inner cross-validation (glmnet).
#
# The projection (PCA / kernel PCA) was fitted without labels on all
# patients and is kept fixed; only the readout is cross-fitted.
#
# Also computed per cell: the in-sample plug-in log-likelihood ratio of
# the full-data PathMI readout, which bounds I_usable from above under
# the same readout.
#
# Usage (project root, R 4.4.3):  Rscript TRIDENT/usable_recoverability.R
# Output: TRIDENT/usable/usable_recoverability.csv (one row per cell),
#         TRIDENT/usable/usable_by_repeat.csv (one row per cell x repeat)
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(glmnet); library(parallel) })

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI  <- normalizePath(dirname(this_file), winslash = "/")
ROOT <- normalizePath(file.path(TRI, ".."), winslash = "/")
setwd(ROOT)
source("PathMI_SUBHO_logistic.R")

OUT <- file.path(TRI, "usable"); dir.create(OUT, showWarnings = FALSE)
LOG <- file.path(OUT, "usable_recoverability.log")
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }

N_WORKERS <- 10L; K_FOLDS <- 5L; R_REP <- 10L
V_BETA <- 1.0; NMC <- 2000L; BURN <- 500L; THIN <- 10L; SEED <- 42L

MODELS <- c(uni_v2 = "UNI v2", virchow2 = "Virchow2", hoptimus1 = "H-optimus-1", conch_v15 = "CONCH v1.5")
runs <- expand.grid(cohort = c("crc", "prad"), model = names(MODELS), stringsAsFactors = FALSE)

# ---- cells -------------------------------------------------------------------
cells <- list()
for (i in seq_len(nrow(runs))) {
  co <- runs$cohort[i]; mo <- runs$model[i]
  f <- if (co == "crc") sprintf("crc_kpca_subho_%s/kpca_subho_results_crc_%s.rds", mo, mo)
       else sprintf("TCGA_PRAD_kpca_subho_%s/kpca_subho_results_prad_%s.rds", mo, mo)
  res <- readRDS(file.path(TRI, f))
  for (key in names(res)) {
    a <- res[[key]]; roster <- rownames(a$Z_scaled)
    ylist <- if (co == "crc") list(MSI = a$y[roster]) else lapply(a$y_list, function(v) v[roster])
    for (tg in names(ylist)) {
      det <- a$result$details[[tg]]
      cells[[length(cells) + 1]] <- list(
        id = sprintf("%s|%s|%s|%s", toupper(co), mo, tg, key),
        cohort = toupper(co), model = MODELS[[mo]], target = tg, key = key,
        decomp = sub("_.*", "", key), pooling = sub(".*_", "", key),
        Z = unclass(a$Z_scaled)[, , drop = FALSE], y = as.numeric(ylist[[tg]]),
        MI_obs = det$MI, null_mean = mean(det$MI.null, na.rm = TRUE), null_sd = sd(det$MI.null, na.rm = TRUE),
        p_value = det$p.value, k_star_full = det$n.selected)
    }
  }
}
stopifnot(length(cells) == 192)
logf("=== usable recoverability | %d cells | K = %d folds x R = %d repeats | %d workers ===", length(cells), K_FOLDS, R_REP, N_WORKERS)

# ---- one job: one cell, one repeat (r = 0 is the full-data in-sample fit) ----
auc <- function(s, y) { r <- rank(s); n1 <- sum(y); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1)) }
clip <- function(p) pmin(pmax(p, 1e-12), 1 - 1e-12)

pathmi_readout <- function(Ztr, ytr, Zte, seed) {
  set.seed(seed)
  sel <- selection.logistic(X = cbind(1, Ztr), y = ytr, iter = 1L, v_beta = V_BETA,
                            method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
  best <- which.max(sel$Likelihood)
  pbar <- mean(ytr)
  floored <- (sel$LogML.only[best] - sel$LogML.null) <= 0
  eta <- if (floored) rep(qlogis(pbar), nrow(Zte)) else as.vector(cbind(1, Zte) %*% sel$Sp.BetaHat)
  list(eta = eta, k = best, floored = floored, pbar = pbar,
       MI = max((sel$LogML.only[best] - sel$LogML.null) / length(ytr), 0))
}

job <- function(ci, r) {
  cell <- cells[[ci]]; Z <- cell$Z; y <- cell$y; n <- length(y)
  base_of <- function(yy, pb) yy * log(pb) + (1 - yy) * log(1 - pb)
  if (r == 0L) {                                  # in-sample, full data
    ro <- pathmi_readout(Z, y, Z, SEED)
    p  <- clip(plogis(ro$eta))
    l  <- y * log(p) + (1 - y) * log(1 - p) - base_of(y, ro$pbar)
    return(data.frame(id = cell$id, rep = 0L, I_pathmi = mean(l), auc_pathmi = auc(ro$eta, y),
                      I_ridge = NA_real_, auc_ridge = NA_real_, k_mean = ro$k, k_min = ro$k, k_max = ro$k,
                      n_floored = as.integer(ro$floored), MI_refit = ro$MI, stringsAsFactors = FALSE))
  }
  set.seed(1000L + r)
  fold <- integer(n)
  for (cls in 0:1) { idx <- which(y == cls); fold[idx] <- sample(rep_len(seq_len(K_FOLDS), length(idx))) }
  l_p <- l_r <- eta_p <- eta_r <- numeric(n); ks <- integer(K_FOLDS); fl <- logical(K_FOLDS)
  for (k in seq_len(K_FOLDS)) {
    tr <- fold != k; te <- fold == k
    ro <- pathmi_readout(Z[tr, , drop = FALSE], y[tr], Z[te, , drop = FALSE], SEED + 100L * r + k)
    eta_p[te] <- ro$eta; ks[k] <- ro$k; fl[k] <- ro$floored
    p <- clip(plogis(ro$eta))
    l_p[te] <- y[te] * log(p) + (1 - y[te]) * log(1 - p) - base_of(y[te], ro$pbar)
    set.seed(SEED + 100L * r + k)
    cvr <- cv.glmnet(Z[tr, , drop = FALSE], y[tr], family = "binomial", alpha = 0, nfolds = 10)
    er  <- as.vector(predict(cvr, Z[te, , drop = FALSE], s = "lambda.min", type = "link"))
    eta_r[te] <- er; pr <- clip(plogis(er))
    l_r[te] <- y[te] * log(pr) + (1 - y[te]) * log(1 - pr) - base_of(y[te], ro$pbar)
  }
  data.frame(id = cell$id, rep = r, I_pathmi = mean(l_p), auc_pathmi = auc(eta_p, y),
             I_ridge = mean(l_r), auc_ridge = auc(eta_r, y), k_mean = mean(ks), k_min = min(ks), k_max = max(ks),
             n_floored = sum(fl), MI_refit = NA_real_, stringsAsFactors = FALSE)
}

jobs <- expand.grid(ci = seq_along(cells), r = 0:R_REP)
jobs <- jobs[sample.int(nrow(jobs)), ]            # mix cells so progress is even
logf("%d jobs", nrow(jobs))

cl <- makeCluster(N_WORKERS); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(glmnet) }))
clusterExport(cl, c("cells", "jobs", "job", "pathmi_readout", "auc", "clip", "K_FOLDS", "V_BETA", "NMC", "BURN", "THIN", "SEED",
                    "sagnik.horseshoe.logistic", "lpy.x.logistic", "selection.logistic", "sigmoid", "clip_eta"))
t0 <- Sys.time(); out <- vector("list", nrow(jobs)); done <- 0L
chunks <- split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / 120))
for (ch in chunks) {
  o <- parLapply(cl, ch, function(j) tryCatch(job(jobs$ci[j], jobs$r[j]),
                   error = function(e) data.frame(id = cells[[jobs$ci[j]]]$id, rep = jobs$r[j], I_pathmi = NA, auc_pathmi = NA,
                     I_ridge = NA, auc_ridge = NA, k_mean = NA, k_min = NA, k_max = NA, n_floored = NA, MI_refit = NA, err = conditionMessage(e))))
  for (i in seq_along(ch)) out[[ch[i]]] <- o[[i]]
  done <- done + length(ch); el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  logf("  %4d / %d jobs | %.1f min | ETA %.1f min", done, nrow(jobs), el, el / done * (nrow(jobs) - done))
}
rep_df <- do.call(rbind, lapply(out, function(d) { if (is.null(d$err)) d$err <- NA_character_; d }))
write.csv(rep_df, file.path(OUT, "usable_by_repeat.csv"), row.names = FALSE)
if (any(!is.na(rep_df$err))) logf("[WARN] %d job(s) errored: %s", sum(!is.na(rep_df$err)), paste(head(unique(rep_df$err[!is.na(rep_df$err)])), collapse = " | "))

# ---- summary per cell --------------------------------------------------------
meta <- do.call(rbind, lapply(cells, function(c) data.frame(id = c$id, cohort = c$cohort, model = c$model, target = c$target,
  key = c$key, decomp = c$decomp, pooling = c$pooling, n = length(c$y), prevalence = mean(c$y),
  MI_obs = c$MI_obs, null_mean = c$null_mean, null_sd = c$null_sd, p_value = c$p_value, k_star_full = c$k_star_full, stringsAsFactors = FALSE)))
full <- rep_df[rep_df$rep == 0L, c("id", "I_pathmi", "auc_pathmi", "MI_refit", "k_mean")]
names(full) <- c("id", "I_insample_plugin", "auc_insample", "MI_refit", "k_star_refit")
cv <- rep_df[rep_df$rep > 0L, ]
agg <- do.call(rbind, lapply(split(cv, cv$id), function(d) data.frame(id = d$id[1],
  I_usable = mean(d$I_pathmi), I_usable_sd = sd(d$I_pathmi), I_usable_lo = min(d$I_pathmi), I_usable_hi = max(d$I_pathmi),
  auc_oof = mean(d$auc_pathmi), I_usable_ridge = mean(d$I_ridge), I_usable_ridge_sd = sd(d$I_ridge), auc_oof_ridge = mean(d$auc_ridge),
  k_oof_mean = mean(d$k_mean), floored_folds = sum(d$n_floored), n_repeats = nrow(d), stringsAsFactors = FALSE)))
res <- merge(merge(meta, full, by = "id"), agg, by = "id")
res$H_Y <- with(res, -(prevalence * log(prevalence) + (1 - prevalence) * log(1 - prevalence)))
res$usable_share <- res$I_usable / res$MI_obs
res <- res[order(res$cohort, res$target, res$model, res$key), ]
write.csv(res, file.path(OUT, "usable_recoverability.csv"), row.names = FALSE)

logf("=== done in %.1f min ===", as.numeric(difftime(Sys.time(), t0, units = "mins")))
logf("refit reproduces stored estimate: max |MI_refit - MI_obs| = %.2e", max(abs(res$MI_refit - res$MI_obs)))
logf("I_usable (PathMI readout): median %.4f | range %.4f to %.4f | negative in %d cells | r(I_usable, MI_obs) = %.3f (Pearson), %.3f (Spearman)",
     median(res$I_usable), min(res$I_usable), max(res$I_usable), sum(res$I_usable < 0), cor(res$I_usable, res$MI_obs), cor(res$I_usable, res$MI_obs, method = "spearman"))
logf("I_usable (ridge): median %.4f | r with MI_obs = %.3f", median(res$I_usable_ridge), cor(res$I_usable_ridge, res$MI_obs))
for (tk in unique(paste(res$cohort, res$target))) { d <- res[paste(res$cohort, res$target) == tk, ]
  logf("  %-9s usable/present median %.2f | AUC oof median %.3f | r = %.3f", tk, median(d$usable_share), median(d$auc_oof), cor(d$I_usable, d$MI_obs)) }

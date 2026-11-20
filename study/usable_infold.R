################################################################
# usable_infold.R  (patient-level version, 19 Nov 2026)
#
# Usable recoverability with every data-dependent transformation
# confined to training patients, on the patient-level rosters (one row
# per patient; pool_pass/<tag>/pooled.rds as written by
# patient_level/pl_promote.R). For each outer fold:
#   training patients -> column standardisation -> PCA / kernel PCA
#   (bandwidth, centring, eigenvectors, score scaling) -> readout
#   (ridge logistic with lambda tuned by inner CV on training patients;
#   and the PathMI readout: ranking, walk and MAP refitted) -> held-out
#   patients (projected with the training fit).
# Also, in the same folds, the ridge readout on the projection fitted
# once on all patients (label-free), to measure what refitting the
# projection changes; and, on all patients, the in-sample plug-in of
# the PathMI readout.
# Stratified 5-fold x 10 repeats (fold seeds 1000 + r). All 192 cells,
# plus the negative control: the H-optimus-1 CRC/MSI selected arm with
# its labels permuted once (seed 20261119).
# The slide-level version is kept as usable_infold_slide_version.R.
#
# Usage (project root, R 4.4.3): Rscript TRIDENT/usable_infold.R
# Output: TRIDENT/usable/usable_infold_by_repeat.csv, usable_infold.csv,
#         and usable_recoverability.csv in the format the earlier
#         consumers read (ridge = in-fold; PathMI readout = in-fold).
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(glmnet); library(parallel); library(readxl) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI  <- normalizePath(dirname(this_file), winslash = "/"); ROOT <- normalizePath(file.path(TRI, ".."), winslash = "/")
setwd(ROOT); source("PathMI_SUBHO_logistic.R"); source(file.path(TRI, "infold_functions.R"))
OUT <- file.path(TRI, "usable"); dir.create(OUT, showWarnings = FALSE); LOG <- file.path(OUT, "usable_infold.log")
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }
N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "24")); K_FOLDS <- 5L; R_REP <- 10L
V_BETA <- 1.0; NMC <- 2000L; BURN <- 500L; THIN <- 10L; SEED <- 42L; CONTROL_SEED <- 20261119L
POOL <- "M:/bhadury/mi_results/pool_pass"
MODELS <- c(uni_v2 = "UNI v2", virchow2 = "Virchow2", hoptimus1 = "H-optimus-1", conch_v15 = "CONCH v1.5")
g <- as.data.frame(read_excel(file.path(TRI, "mi summary table.xlsx")))
emb2code <- c("UNI V2" = "uni_v2", "VIRCHOW" = "virchow2", "H-OPTIMUS 1" = "hoptimus1", "CONCH V1.5" = "conch_v15")
g$code <- emb2code[g$Embeddings]; g$id <- sprintf("%s|%s|%s|%s", g$Disease, g$code, g$target, g$key)

cells <- list()
for (co in c("crc", "prad")) for (mo in names(MODELS)) {
  pp <- readRDS(file.path(POOL, sprintf("%s_%s", co, mo), "pooled.rds")); stopifnot(identical(pp$unit, "patient"))
  for (tg in names(pp$y)) for (dc in c("pca", "cosine", "rbf", "poly")) for (pl in c("mean", "norm", "moment", "attention")) {
    id <- sprintf("%s|%s|%s|%s_%s", toupper(co), mo, tg, dc, pl); gi <- g[g$id == id, ]
    cells[[id]] <- list(id = id, cohort = toupper(co), code = mo, model = MODELS[[mo]], target = tg, decomp = dc, pooling = pl, key = paste0(dc, "_", pl),
                        y = as.numeric(pp$y[[tg]][pp$roster]), group = pp$roster, MI_obs = gi$MI_obs, null_mean = gi$null_mean, null_sd = gi$null_sd, p_value = gi$p_value, control = FALSE)
  }
}
hsel <- g$Disease == "CRC" & g$code == "hoptimus1"; ctl_key <- g$key[hsel][which.max(g$MI_obs[hsel])]
ctl <- cells[[sprintf("CRC|hoptimus1|MSI|%s", ctl_key)]]; set.seed(CONTROL_SEED); ctl$y <- sample(ctl$y)
ctl$id <- sprintf("CONTROL|hoptimus1|MSI-permuted|%s", ctl_key); ctl$control <- TRUE; ctl$MI_obs <- NA; cells[[ctl$id]] <- ctl
stopifnot(length(cells) == 193)
logf("=== in-fold usable recoverability (patient level) | %d cells | %d folds x %d repeats | control arm %s | %d workers ===", length(cells), K_FOLDS, R_REP, ctl_key, N_WORKERS)

pool_cache <- new.env()
get_X <- function(cell) { key <- sprintf("%s_%s", tolower(cell$cohort), cell$code)
  if (is.null(pool_cache[[key]])) pool_cache[[key]] <- readRDS(file.path(POOL, key, "pooled.rds"))$pooled
  pool_cache[[key]][[cell$pooling]] }
job <- function(cid, r) {
  cell <- cells[[cid]]; X <- get_X(cell); y <- cell$y
  if (r == 0L) {                                             # in-sample plug-in of the PathMI readout on all patients
    Z <- unclass(apply_projection(fit_projection(X, cell$decomp, 50L), X)); set.seed(SEED)
    sel <- selection.logistic(X = cbind(1, Z), y = y, iter = 1L, v_beta = V_BETA, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
    b <- which.max(sel$Likelihood); eta <- as.vector(cbind(1, Z) %*% sel$Sp.BetaHat)
    return(data.frame(id = cid, rep = 0L, I_ridge = NA, auc_ridge = NA, I_ridge_fixed = NA, I_pathmi = mean(loglik_gain(y, eta, mean(y))), auc_pathmi = auc(eta, y),
                      MI_refit = max((sel$LogML.only[b] - sel$LogML.null) / length(y), 0), k_walk_mean = b, n_floored = NA))
  }
  cf <- crossfit_infold(X, y, cell$decomp, K_FOLDS, fold_seed = 1000L + r, model_seed = SEED + 100L * r, pathmi = TRUE,
                        v_beta = V_BETA, nmc = NMC, burn = BURN, thin = THIN, group = cell$group)
  # ridge on the projection fitted once on all patients, same folds and seeds
  Zf <- unclass(apply_projection(fit_projection(X, cell$decomp, 50L), X)); set.seed(1000L + r); fold <- strat_folds(y, K_FOLDS, cell$group); lf <- numeric(length(y))
  for (k in seq_len(K_FOLDS)) { tr <- fold != k; te <- fold == k; er <- ridge_eta(Zf[tr, , drop = FALSE], y[tr], Zf[te, , drop = FALSE], SEED + 100L * r + k, cell$group[tr])
    lf[te] <- loglik_gain(y[te], er, mean(y[tr])) }
  data.frame(id = cid, rep = r, I_ridge = cf$I_ridge, auc_ridge = cf$auc_ridge, I_ridge_fixed = mean(lf), I_pathmi = cf$I_pathmi, auc_pathmi = cf$auc_pathmi,
             MI_refit = NA, k_walk_mean = cf$k_walk_mean, n_floored = cf$n_floored)
}
jobs <- expand.grid(cid = names(cells), r = 0:R_REP, stringsAsFactors = FALSE); set.seed(7); jobs <- jobs[sample.int(nrow(jobs)), ]
cl <- makeCluster(N_WORKERS); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(glmnet) }))
clusterExport(cl, c("cells", "jobs", "job", "get_X", "POOL", "K_FOLDS", "V_BETA", "NMC", "BURN", "THIN", "SEED",
                    "fit_projection", "apply_projection", "strat_folds", "clip01", "loglik_gain", "auc", "ridge_eta", "crossfit_infold", "EIG_TOL",
                    "sagnik.horseshoe.logistic", "lpy.x.logistic", "selection.logistic", "sigmoid", "clip_eta"))
clusterEvalQ(cl, { pool_cache <- new.env(); NULL })
t0 <- Sys.time(); out <- vector("list", nrow(jobs)); done <- 0L
for (ch in split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / 192))) {
  o <- parLapply(cl, ch, function(j) tryCatch(job(jobs$cid[j], jobs$r[j]), error = function(e) data.frame(id = jobs$cid[j], rep = jobs$r[j], err = conditionMessage(e))))
  for (i in seq_along(ch)) out[[ch[i]]] <- o[[i]]
  done <- done + length(ch); el <- as.numeric(difftime(Sys.time(), t0, units = "mins")); logf("  %4d / %d jobs | %.1f min | ETA %.1f min", done, nrow(jobs), el, el / done * (nrow(jobs) - done))
}
errs <- Filter(function(d) !is.null(d$err), out); if (length(errs)) logf("[WARN] %d job(s) errored: %s", length(errs), paste(head(unique(sapply(errs, `[[`, "err"))), collapse = " | "))
rep_df <- do.call(rbind, Filter(function(d) is.null(d$err), out)); write.csv(rep_df, file.path(OUT, "usable_infold_by_repeat.csv"), row.names = FALSE)

hb <- function(p) -(p * log(p) + (1 - p) * log(1 - p))
summ <- do.call(rbind, lapply(split(rep_df, rep_df$id), function(d) { c0 <- cells[[d$id[1]]]; cv <- d[d$rep > 0, ]; f0 <- d[d$rep == 0, ]
  data.frame(id = d$id[1], cohort = c0$cohort, model = c0$model, target = c0$target, key = c0$key, decomp = c0$decomp, pooling = c0$pooling,
    n = length(c0$y), prevalence = mean(c0$y), H_Y = hb(mean(c0$y)), MI_obs = c0$MI_obs, null_mean = c0$null_mean, null_sd = c0$null_sd, p_value = c0$p_value,
    I_usable_infold = mean(cv$I_ridge), I_usable_infold_sd = sd(cv$I_ridge), auc_infold = mean(cv$auc_ridge),
    I_usable_fixed = mean(cv$I_ridge_fixed), diff_infold_minus_fixed = mean(cv$I_ridge - cv$I_ridge_fixed),
    I_usable_pathmi_infold = mean(cv$I_pathmi), I_usable_pathmi_sd = sd(cv$I_pathmi), auc_pathmi_infold = mean(cv$auc_pathmi), floored_folds = sum(cv$n_floored), k_walk_mean = mean(cv$k_walk_mean),
    I_insample_plugin = f0$I_pathmi, MI_refit = f0$MI_refit, control = c0$control, record = FALSE, stringsAsFactors = FALSE) }))
write.csv(summ, file.path(OUT, "usable_infold.csv"), row.names = FALSE)
comp <- with(summ[!summ$control, ], data.frame(id, cohort, model, target, key, decomp, pooling, n, prevalence, MI_obs, null_mean, null_sd, p_value,
  I_insample_plugin, MI_refit, I_usable = I_usable_pathmi_infold, I_usable_sd = I_usable_pathmi_sd, auc_oof = auc_pathmi_infold,
  I_usable_ridge = I_usable_infold, I_usable_ridge_sd = I_usable_infold_sd, auc_oof_ridge = auc_infold, k_oof_mean = k_walk_mean, floored_folds, H_Y,
  usable_share = I_usable_pathmi_infold / MI_obs))
write.csv(comp, file.path(OUT, "usable_recoverability.csv"), row.names = FALSE)
s <- summ[!summ$control, ]
logf("=== done in %.1f min ===", as.numeric(difftime(Sys.time(), t0, units = "mins")))
logf("refit reproduces stored estimate: max |MI_refit - MI_obs| = %.2e", max(abs(s$MI_refit - s$MI_obs)))
logf("192 cells ridge in-fold: median %.4f | fixed projection %.4f | median in-fold minus fixed %+.4f | above null band %d | r with PathMI %.3f | exceeds PathMI in %d",
     median(s$I_usable_infold), median(s$I_usable_fixed), median(s$diff_infold_minus_fixed), sum(s$I_usable_infold > s$null_mean + 2 * s$null_sd), cor(s$I_usable_infold, s$MI_obs), sum(s$I_usable_infold > s$MI_obs))
c1 <- summ[summ$control, ]; logf("control: ridge %.4f (AUC %.3f) | PathMI readout %.4f", c1$I_usable_infold, c1$auc_infold, c1$I_usable_pathmi_infold)

################################################################
# robustness_prior_bootstrap.R
#
# Estimator robustness on the saved 50-dimensional representations.
#
#   A. Working-prior sweep. For eight configurations (the three
#      primary ones, i.e. the best configuration of each state, and
#      five deliberately difficult ones), PathMI is refitted at
#      v_beta in {0.1, 0.3, 1, 3, 10, 30, 100} with 100 label
#      permutations at every value, so that both the estimate and its
#      null floor are followed across the prior.
#   B. Chain and size-prior sensitivity, same configurations: five
#      horseshoe chain seeds at v_beta = 1, and the basic-binomial
#      model-size prior in place of the beta-binomial.
#   C. Patient bootstrap. For the 48 selected arms and the twelve
#      linear-PCA / mean-pooling arms (the tile-level readout), B = 200
#      resamples of patients with replacement; the full estimator is
#      refitted on each. Reported: bootstrap s.e. and 5/95 percentiles.
#
# Usage (project root, R 4.4.3): Rscript TRIDENT/robustness_prior_bootstrap.R
# Output: TRIDENT/robustness/{vbeta_sweep,chain_sizeprior,bootstrap_draws,bootstrap_summary}.csv
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(parallel); library(readxl) })

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI  <- normalizePath(dirname(this_file), winslash = "/")
ROOT <- normalizePath(file.path(TRI, ".."), winslash = "/")
setwd(ROOT)
source("PathMI_SUBHO_logistic.R")

OUT <- file.path(TRI, "robustness"); dir.create(OUT, showWarnings = FALSE)
LOG <- file.path(OUT, "robustness_prior_bootstrap.log")
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }

N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "12"))
NMC <- 2000L; BURN <- 500L; THIN <- 10L; SEED <- 42L
VB_GRID <- c(0.1, 0.3, 1, 3, 10, 30, 100); N_PERM_VB <- 100L
CHAIN_SEEDS <- 1:5
B_BOOT <- 200L

MODELS <- c(uni_v2 = "UNI v2", virchow2 = "Virchow2", hoptimus1 = "H-optimus-1", conch_v15 = "CONCH v1.5")
harm <- function(e) c("UNI V2" = "uni_v2", "VIRCHOW" = "virchow2", "H-OPTIMUS 1" = "hoptimus1", "CONCH V1.5" = "conch_v15")[[e]]

# ---- load a cell -------------------------------------------------------------
res_cache <- list()
load_cell <- function(cohort, model, target, key) {
  f <- if (cohort == "CRC") sprintf("crc_kpca_subho_%s/kpca_subho_results_crc_%s.rds", model, model)
       else sprintf("TCGA_PRAD_kpca_subho_%s/kpca_subho_results_prad_%s.rds", model, model)
  if (is.null(res_cache[[f]])) res_cache[[f]] <<- readRDS(file.path(TRI, f))
  a <- res_cache[[f]][[key]]; roster <- rownames(a$Z_scaled)
  y <- if (cohort == "CRC") a$y[roster] else a$y_list[[target]][roster]
  det <- a$result$details[[target]]
  list(id = sprintf("%s|%s|%s|%s", cohort, model, target, key), cohort = cohort, model = MODELS[[model]], model_code = model,
       target = target, key = key, Z = unclass(a$Z_scaled), y = as.numeric(y),
       MI_obs = det$MI, null_mean = mean(det$MI.null, na.rm = TRUE), null_sd = sd(det$MI.null, na.rm = TRUE), k_star = det$n.selected)
}

# ---- configurations ----------------------------------------------------------
sens_cfg <- data.frame(
  role   = c("primary", "primary", "primary", "difficult", "difficult", "difficult", "difficult", "difficult"),
  note   = c("best configuration, MSI", "best configuration, ERG fusion", "best configuration, ETS fusion",
             "largest departure from additivity (CONCH v1.5, moment pooling)", "lower-tier model on the hardest state",
             "near the permutation floor (CONCH v1.5, ERG)", "tile-level readout arm (linear PCA, mean pooling)",
             "the arm behind the state-dependent winner (H-optimus-1 on ERG)"),
  cohort = c("CRC", "PRAD", "PRAD", "CRC", "PRAD", "PRAD", "CRC", "PRAD"),
  model  = c("hoptimus1", "uni_v2", "hoptimus1", "conch_v15", "virchow2", "conch_v15", "hoptimus1", "hoptimus1"),
  target = c("MSI", "ERG", "ETS", "MSI", "ERG", "ERG", "MSI", "ERG"),
  key    = c("cosine_moment", "cosine_norm", "cosine_moment", "poly_moment", "cosine_moment", "poly_moment", "pca_mean", "pca_moment"),
  stringsAsFactors = FALSE)
# 19 Nov 2026 (patient-level rerun): the primary configurations are the best of each state in the current grid
.g <- as.data.frame(readxl::read_excel(file.path(TRI, "mi summary table.xlsx")))
.e2c <- c("UNI V2" = "uni_v2", "VIRCHOW" = "virchow2", "H-OPTIMUS 1" = "hoptimus1", "CONCH V1.5" = "conch_v15")
for (.st in c("MSI", "ERG", "ETS")) { .d <- .g[.g$target == .st, ]; .b <- .d[which.max(.d$MI_obs), ]
  .i <- which(sens_cfg$role == "primary" & sens_cfg$target == .st); sens_cfg$model[.i] <- .e2c[[.b$Embeddings]]; sens_cfg$key[.i] <- .b$key }
sens_cells <- lapply(seq_len(nrow(sens_cfg)), function(i) with(sens_cfg[i, ], load_cell(cohort, model, target, key)))
names(sens_cells) <- sapply(sens_cells, `[[`, "id")

best <- as.data.frame(read_excel(file.path(TRI, "best kpca summary.xlsx")))
boot_cfg <- unique(rbind(
  data.frame(cohort = best$Disease, model = sapply(best$Embeddings, harm), target = best$target, key = best$key, arm = "selected", stringsAsFactors = FALSE),
  expand.grid(cohort = "CRC", model = names(MODELS), target = "MSI", key = "pca_mean", arm = "linear_mean", stringsAsFactors = FALSE),
  expand.grid(cohort = "PRAD", model = names(MODELS), target = c("ERG", "ETS"), key = "pca_mean", arm = "linear_mean", stringsAsFactors = FALSE)))
boot_cfg$arm[duplicated(boot_cfg[, c("cohort", "model", "target", "key")]) | duplicated(boot_cfg[, c("cohort", "model", "target", "key")], fromLast = TRUE)] <- "selected+linear_mean"
boot_cfg <- boot_cfg[!duplicated(boot_cfg[, c("cohort", "model", "target", "key")]), ]
boot_cells <- lapply(seq_len(nrow(boot_cfg)), function(i) with(boot_cfg[i, ], load_cell(cohort, model, target, key)))
names(boot_cells) <- sapply(boot_cells, `[[`, "id")
logf("=== robustness | %d sensitivity configurations | %d bootstrap arms x B = %d | %d workers ===", length(sens_cells), length(boot_cells), B_BOOT, N_WORKERS)

# ---- fit helpers ---------------------------------------------------------------
fit_one <- function(Z, y, v_beta = 1, seed = SEED, method.weight = "beta.binomial") {
  set.seed(seed)
  r <- subho_MI_logistic(Z = Z, y = y, v_beta = v_beta, method.weight = method.weight, nmc = NMC, burn = BURN, thin = THIN)
  c(MI = r$MI, MI_raw = r$MI.raw, k = r$n.selected, cond = r$cond.H.best)
}

# A: (cell, v_beta, block) ; block 0 = observed, blocks 1.. = 25 permutations each
job_vb <- function(cid, vb, block) {
  cell <- sens_cells[[cid]]; Z <- cell$Z; y <- cell$y
  if (block == 0L) { f <- fit_one(Z, y, vb, SEED)
    return(data.frame(id = cid, v_beta = vb, type = "observed", b = 0L, MI = f[["MI"]], MI_raw = f[["MI_raw"]], k = f[["k"]], cond = f[["cond"]])) }
  bs <- ((block - 1L) * 25L + 1L):(block * 25L)
  do.call(rbind, lapply(bs, function(b) { set.seed(10000L + b); yp <- sample(y); f <- fit_one(Z, yp, vb, 20000L + b)
    data.frame(id = cid, v_beta = vb, type = "perm", b = b, MI = f[["MI"]], MI_raw = f[["MI_raw"]], k = f[["k"]], cond = f[["cond"]]) }))
}
# B: chain seeds and size prior
job_chain <- function(cid, what, s) {
  cell <- sens_cells[[cid]]
  f <- if (what == "seed") fit_one(cell$Z, cell$y, 1, s) else fit_one(cell$Z, cell$y, 1, SEED, "basic.binomial")
  data.frame(id = cid, what = what, seed = s, MI = f[["MI"]], MI_raw = f[["MI_raw"]], k = f[["k"]], cond = f[["cond"]])
}
# C: bootstrap, blocks of 25 draws
job_boot <- function(cid, block) {
  cell <- boot_cells[[cid]]; Z <- cell$Z; y <- cell$y; n <- length(y)
  bs <- ((block - 1L) * 25L + 1L):(block * 25L)
  do.call(rbind, lapply(bs, function(b) { set.seed(777L + b); idx <- sample.int(n, n, replace = TRUE)
    f <- fit_one(Z[idx, , drop = FALSE], y[idx], 1, 30000L + b)
    data.frame(id = cid, b = b, MI = f[["MI"]], MI_raw = f[["MI_raw"]], k = f[["k"]], prevalence = mean(y[idx])) }))
}

jobs <- rbind(
  do.call(rbind, lapply(names(sens_cells), function(c) expand.grid(kind = "vb", cid = c, a = VB_GRID, b = 0:(N_PERM_VB / 25L), stringsAsFactors = FALSE))),
  do.call(rbind, lapply(names(sens_cells), function(c) rbind(data.frame(kind = "chain", cid = c, a = CHAIN_SEEDS, b = 0L), data.frame(kind = "sizeprior", cid = c, a = 0, b = 0L)))),
  do.call(rbind, lapply(names(boot_cells), function(c) data.frame(kind = "boot", cid = c, a = 0, b = seq_len(B_BOOT / 25L)))))
jobs <- jobs[sample.int(nrow(jobs)), ]
logf("%d jobs (%d prior-sweep blocks, %d chain/size-prior fits, %d bootstrap blocks)", nrow(jobs), sum(jobs$kind == "vb"), sum(jobs$kind %in% c("chain", "sizeprior")), sum(jobs$kind == "boot"))

run_job <- function(j) {
  r <- jobs[j, ]
  tryCatch(switch(r$kind,
    vb        = job_vb(r$cid, r$a, as.integer(r$b)),
    chain     = job_chain(r$cid, "seed", as.integer(r$a)),
    sizeprior = job_chain(r$cid, "sizeprior", SEED),
    boot      = job_boot(r$cid, as.integer(r$b))),
    error = function(e) data.frame(id = r$cid, err = conditionMessage(e), kind = r$kind))
}

cl <- makeCluster(N_WORKERS); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(MASS); library(BayesLogit) }))
clusterExport(cl, c("sens_cells", "boot_cells", "jobs", "fit_one", "job_vb", "job_chain", "job_boot", "NMC", "BURN", "THIN", "SEED",
                    "sagnik.horseshoe.logistic", "lpy.x.logistic", "selection.logistic", "subho_MI_logistic", "sigmoid", "clip_eta"))
t0 <- Sys.time(); out <- vector("list", nrow(jobs)); done <- 0L
for (ch in split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / 96))) {
  o <- parLapply(cl, ch, run_job); for (i in seq_along(ch)) out[[ch[i]]] <- o[[i]]
  done <- done + length(ch); el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  logf("  %4d / %d jobs | %.1f min | ETA %.1f min", done, nrow(jobs), el, el / done * (nrow(jobs) - done))
}
errs <- Filter(function(d) !is.null(d$err), out)
if (length(errs)) logf("[WARN] %d job(s) errored: %s", length(errs), paste(head(unique(sapply(errs, `[[`, "err"))), collapse = " | "))
ok <- Filter(function(d) is.null(d$err), out)
vb_df   <- do.call(rbind, ok[sapply(ok, function(d) "v_beta" %in% names(d))])
ch_df   <- do.call(rbind, ok[sapply(ok, function(d) "what" %in% names(d))])
boot_df <- do.call(rbind, ok[sapply(ok, function(d) "prevalence" %in% names(d))])

# ---- A summary ------------------------------------------------------------------
meta_s <- do.call(rbind, lapply(sens_cells, function(c) data.frame(id = c$id, role = sens_cfg$role[match(c$id, names(sens_cells))], note = sens_cfg$note[match(c$id, names(sens_cells))],
  cohort = c$cohort, model = c$model, target = c$target, key = c$key, n = length(c$y), MI_table = c$MI_obs, null_mean_table = c$null_mean, null_sd_table = c$null_sd, k_table = c$k_star, stringsAsFactors = FALSE)))
vb_sum <- do.call(rbind, lapply(split(vb_df, list(vb_df$id, vb_df$v_beta), drop = TRUE), function(d) {
  o <- d[d$type == "observed", ]; p <- d[d$type == "perm", ]
  data.frame(id = d$id[1], v_beta = d$v_beta[1], MI = o$MI, MI_raw = o$MI_raw, k = o$k, cond = o$cond,
             null_mean = mean(p$MI), null_sd = sd(p$MI), p_perm = (1 + sum(p$MI >= o$MI)) / (1 + nrow(p)), z = (o$MI - mean(p$MI)) / sd(p$MI), n_perm = nrow(p)) }))
vb_sum <- merge(meta_s, vb_sum, by = "id"); vb_sum <- vb_sum[order(vb_sum$role, vb_sum$id, vb_sum$v_beta), ]
write.csv(vb_sum, file.path(OUT, "vbeta_sweep.csv"), row.names = FALSE)
write.csv(vb_df, file.path(OUT, "vbeta_sweep_draws.csv"), row.names = FALSE)
ch_sum <- merge(meta_s[, c("id", "role", "cohort", "model", "target", "key", "MI_table", "k_table")], ch_df, by = "id")
write.csv(ch_sum[order(ch_sum$id, ch_sum$what, ch_sum$seed), ], file.path(OUT, "chain_sizeprior.csv"), row.names = FALSE)

# ---- C summary ------------------------------------------------------------------
meta_b <- do.call(rbind, lapply(seq_along(boot_cells), function(i) { c <- boot_cells[[i]]; data.frame(id = c$id, arm = boot_cfg$arm[i], cohort = c$cohort, model = c$model, target = c$target, key = c$key,
  n = length(c$y), MI_obs = c$MI_obs, null_mean = c$null_mean, null_sd = c$null_sd, k_star = c$k_star, stringsAsFactors = FALSE) }))
boot_sum <- do.call(rbind, lapply(split(boot_df, boot_df$id), function(d) data.frame(id = d$id[1], B = nrow(d), boot_mean = mean(d$MI), boot_se = sd(d$MI),
  boot_lo05 = quantile(d$MI, .05), boot_hi95 = quantile(d$MI, .95), boot_lo025 = quantile(d$MI, .025), boot_hi975 = quantile(d$MI, .975),
  k_median = median(d$k), k_lo = quantile(d$k, .05), k_hi = quantile(d$k, .95), frac_zero = mean(d$MI == 0))))
boot_sum <- merge(meta_b, boot_sum, by = "id"); boot_sum$boot_bias <- boot_sum$boot_mean - boot_sum$MI_obs
boot_sum <- boot_sum[order(boot_sum$cohort, boot_sum$target, boot_sum$model, boot_sum$key), ]
write.csv(boot_sum, file.path(OUT, "bootstrap_summary.csv"), row.names = FALSE)
write.csv(boot_df, file.path(OUT, "bootstrap_draws.csv"), row.names = FALSE)

logf("=== done in %.1f min ===", as.numeric(difftime(Sys.time(), t0, units = "mins")))
logf("--- prior sweep: estimate (k*) and null mean by v_beta ---")
for (id in unique(vb_sum$id)) { d <- vb_sum[vb_sum$id == id, ]
  logf("%-45s %s", id, paste(sprintf("vb=%g: %.3f (k=%d, null %.4f, p=%.3f)", d$v_beta, d$MI, d$k, d$null_mean, d$p_perm), collapse = " | ")) }
logf("--- chain seeds (v_beta = 1): range of estimate over 5 seeds; size prior ---")
for (id in unique(ch_sum$id)) { d <- ch_sum[ch_sum$id == id, ]; s <- d[d$what == "seed", ]; b <- d[d$what == "sizeprior", ]
  logf("%-45s table %.4f | seeds %.4f to %.4f (k %d to %d) | basic-binomial %.4f (k=%d)", id, d$MI_table[1], min(s$MI), max(s$MI), min(s$k), max(s$k), b$MI, b$k) }
logf("--- bootstrap ---")
for (i in seq_len(nrow(boot_sum))) with(boot_sum[i, ], logf("%-45s MI %.4f | boot mean %.4f se %.4f | 90%% [%.4f, %.4f] | k* %d (90%% %g-%g)", id, MI_obs, boot_mean, boot_se, boot_lo05, boot_hi95, k_star, k_lo, k_hi))
logf("bootstrap s.e. summary: median %.4f, range %.4f to %.4f; selected arms only: median %.4f", median(boot_sum$boot_se), min(boot_sum$boot_se), max(boot_sum$boot_se), median(boot_sum$boot_se[grepl("selected", boot_sum$arm)]))

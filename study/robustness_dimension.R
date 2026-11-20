################################################################
# robustness_dimension.R
#
# Sweep of the retained dimension d of the projection, for the same
# eight configurations as robustness_prior_bootstrap.R (three primary,
# five difficult). The projection is rebuilt from the pooled matrices
# saved by pool_pass.R with the pipeline's own code at
# d in {10, 20, 30, 50, 75, 100}; d = 50 must reproduce the saved
# Z_scaled (checked). At every d PathMI is fitted once and refitted on
# 100 label permutations, so the estimate is read against the null
# floor at that dimension; the condition number of the retained
# spectrum and k* are recorded.
#
# Usage (project root, R 4.4.3): Rscript TRIDENT/robustness_dimension.R
# Output: TRIDENT/robustness/dimension_sweep.csv, dimension_sweep_draws.csv
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(parallel) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI  <- normalizePath(dirname(this_file), winslash = "/")
ROOT <- normalizePath(file.path(TRI, ".."), winslash = "/")
setwd(ROOT); source("PathMI_SUBHO_logistic.R"); source(file.path(TRI, "projection_functions.R"))
OUT <- file.path(TRI, "robustness"); dir.create(OUT, showWarnings = FALSE)
LOG <- file.path(OUT, "robustness_dimension.log")
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }

N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "12")); NMC <- 2000L; BURN <- 500L; THIN <- 10L; SEED <- 42L
D_GRID <- c(10L, 20L, 30L, 50L, 75L, 100L); N_PERM <- 100L
POOL <- "M:/bhadury/mi_results/pool_pass"
MODELS <- c(uni_v2 = "UNI v2", virchow2 = "Virchow2", hoptimus1 = "H-optimus-1", conch_v15 = "CONCH v1.5")

cfg <- data.frame(
  role   = c("primary", "primary", "primary", "difficult", "difficult", "difficult", "difficult", "difficult"),
  cohort = c("crc", "prad", "prad", "crc", "prad", "prad", "crc", "prad"),
  model  = c("hoptimus1", "uni_v2", "hoptimus1", "conch_v15", "virchow2", "conch_v15", "hoptimus1", "hoptimus1"),
  target = c("MSI", "ERG", "ETS", "MSI", "ERG", "ERG", "MSI", "ERG"),
  key    = c("cosine_moment", "cosine_norm", "cosine_moment", "poly_moment", "cosine_moment", "poly_moment", "pca_mean", "pca_moment"),
  stringsAsFactors = FALSE)
# 19 Nov 2026 (patient-level rerun): the primary configurations are the best of each state in the current grid
.g <- as.data.frame(readxl::read_excel(file.path(TRI, "mi summary table.xlsx")))
.e2c <- c("UNI V2" = "uni_v2", "VIRCHOW" = "virchow2", "H-OPTIMUS 1" = "hoptimus1", "CONCH V1.5" = "conch_v15")
for (.st in c("MSI", "ERG", "ETS")) { .d <- .g[.g$target == .st, ]; .b <- .d[which.max(.d$MI_obs), ]
  .i <- which(cfg$role == "primary" & cfg$target == .st); cfg$model[.i] <- .e2c[[.b$Embeddings]]; cfg$key[.i] <- .b$key }
cfg$id <- with(cfg, sprintf("%s|%s|%s|%s", toupper(cohort), model, target, key))

# ---- build the representations at every d -----------------------------------------
cells <- list()
for (i in seq_len(nrow(cfg))) {
  r <- cfg[i, ]; decomp <- sub("_.*", "", r$key); pooling <- sub(".*_", "", r$key)
  pp <- readRDS(file.path(POOL, sprintf("%s_%s", r$cohort, r$model), "pooled.rds"))
  rds <- file.path(TRI, if (r$cohort == "crc") sprintf("crc_kpca_subho_%s/kpca_subho_results_crc_%s.rds", r$model, r$model) else sprintf("TCGA_PRAD_kpca_subho_%s/kpca_subho_results_prad_%s.rds", r$model, r$model))
  saved <- readRDS(rds)[[r$key]]
  y <- as.numeric(pp$y[[r$target]][pp$roster]); Zraw <- pp$pooled[[pooling]]
  stopifnot(identical(rownames(Zraw), rownames(saved$Z_scaled)))
  for (d in D_GRID) {
    pr <- project(Zraw, decomp, d)
    chk <- if (d == 50L) align_diff(pr$Z_scaled, saved$Z_scaled) else NA_real_
    cells[[sprintf("%s|d=%d", r$id, d)]] <- list(id = r$id, role = r$role, cohort = toupper(r$cohort), model = MODELS[[r$model]], target = r$target, key = r$key,
      d = d, k_use = pr$k_use, eig_cond = pr$eig_cond, var_frac = pr$var_frac, Z = unclass(pr$Z_scaled), y = y,
      MI_table = saved$result$details[[r$target]]$MI, check50 = chk)
    if (d == 50L) logf("%s | d = 50 rebuilt vs saved: max abs diff %.2e (kappa %.1f)", r$id, chk, pr$eig_cond)
  }
}
logf("=== dimension sweep | %d configurations x %d dimensions | %d permutations | %d workers ===", nrow(cfg), length(D_GRID), N_PERM, N_WORKERS)

fit_one <- function(Z, y, seed) { set.seed(seed)
  r <- subho_MI_logistic(Z = Z, y = y, v_beta = 1, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
  c(MI = r$MI, MI_raw = r$MI.raw, k = r$n.selected, cond = r$cond.H.best) }
job <- function(cid, block) {
  cell <- cells[[cid]]
  if (block == 0L) { f <- fit_one(cell$Z, cell$y, SEED); return(data.frame(cid = cid, type = "observed", b = 0L, MI = f[["MI"]], MI_raw = f[["MI_raw"]], k = f[["k"]], cond = f[["cond"]])) }
  bs <- ((block - 1L) * 25L + 1L):(block * 25L)
  do.call(rbind, lapply(bs, function(b) { set.seed(10000L + b); yp <- sample(cell$y); f <- fit_one(cell$Z, yp, 20000L + b)
    data.frame(cid = cid, type = "perm", b = b, MI = f[["MI"]], MI_raw = f[["MI_raw"]], k = f[["k"]], cond = f[["cond"]]) }))
}
jobs <- expand.grid(cid = names(cells), block = 0:(N_PERM / 25L), stringsAsFactors = FALSE); jobs <- jobs[sample.int(nrow(jobs)), ]
cl <- makeCluster(N_WORKERS); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(MASS); library(BayesLogit) }))
clusterExport(cl, c("cells", "jobs", "fit_one", "job", "NMC", "BURN", "THIN", "SEED", "sagnik.horseshoe.logistic", "lpy.x.logistic", "selection.logistic", "subho_MI_logistic", "sigmoid", "clip_eta"))
t0 <- Sys.time(); out <- vector("list", nrow(jobs)); done <- 0L
for (ch in split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / 48))) {
  o <- parLapply(cl, ch, function(j) tryCatch(job(jobs$cid[j], jobs$block[j]), error = function(e) data.frame(cid = jobs$cid[j], err = conditionMessage(e))))
  for (i in seq_along(ch)) out[[ch[i]]] <- o[[i]]
  done <- done + length(ch); el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  logf("  %4d / %d jobs | %.1f min | ETA %.1f min", done, nrow(jobs), el, el / done * (nrow(jobs) - done))
}
errs <- Filter(function(d) !is.null(d$err), out); if (length(errs)) logf("[WARN] %d errored: %s", length(errs), paste(head(unique(sapply(errs, `[[`, "err"))), collapse = " | "))
draws <- do.call(rbind, Filter(function(d) is.null(d$err), out))
write.csv(draws, file.path(OUT, "dimension_sweep_draws.csv"), row.names = FALSE)
summ <- do.call(rbind, lapply(split(draws, draws$cid), function(d) { c <- cells[[d$cid[1]]]; o <- d[d$type == "observed", ]; p <- d[d$type == "perm", ]
  data.frame(id = c$id, role = c$role, cohort = c$cohort, model = c$model, target = c$target, key = c$key, d = c$d, k_use = c$k_use, eig_cond = c$eig_cond, var_frac = c$var_frac,
             MI = o$MI, MI_raw = o$MI_raw, k = o$k, cond_H = o$cond, null_mean = mean(p$MI), null_sd = sd(p$MI), p_perm = (1 + sum(p$MI >= o$MI)) / (1 + nrow(p)),
             z = (o$MI - mean(p$MI)) / sd(p$MI), MI_table_d50 = c$MI_table, check50 = c$check50, stringsAsFactors = FALSE) }))
summ <- summ[order(summ$role, summ$id, summ$d), ]
write.csv(summ, file.path(OUT, "dimension_sweep.csv"), row.names = FALSE)
logf("=== done in %.1f min ===", as.numeric(difftime(Sys.time(), t0, units = "mins")))
for (id in unique(summ$id)) { d <- summ[summ$id == id, ]
  logf("%-45s table(d=50) %.4f | %s", id, d$MI_table_d50[1], paste(sprintf("d=%d: %.4f (k=%d, null %.4f, p=%.3f, kappa %.0f)", d$d, d$MI, d$k, d$null_mean, d$p_perm, d$eig_cond), collapse = " | ")) }

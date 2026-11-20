################################################################
# learning_curves.R
#
# Patient-number learning curves for the 12 Recoverability Record arms
# and the negative control (H-optimus-1 CRC/MSI arm, labels permuted
# once; see usable_infold.R).
#
# For n in {50, 100, 150, 200, floor(N/2)} and B = 50 label-stratified
# subsamples of n patients drawn without replacement from the cohort of
# N, each subsample is analysed as if it were the whole study: the
# projection is refitted on the n patients (retained dimension
# min(50, numerical rank)), PathMI is fitted, PathMI is refitted once on
# permuted labels (the null at that n), and usable recoverability is
# cross-fitted (5 folds, projection and ridge tuning inside the
# training patients).
#
# Uncertainty at n. A subsample of n from N shares patients with every
# other subsample, so the spread of the B estimates understates the
# sampling variability of a study of n independent patients by the
# finite-population factor (1 - n/N). The s.e. at n is therefore
#     se(n) = sd_B(n) / sqrt(1 - n/N),
# which at n = N/2 is the delete-d jackknife with d = N/2 and gives the
# s.e. at N as se(N) = sd_B(N/2) * sqrt(m / (N - m)) (Shao and Wu 1989).
#
# Evidence sufficiency:
#     LCB90(n) = mean_B Ihat(n) - 1.645 se(n)
#     N_evidence = min { n : LCB90(n) > Inull(n) + eps_rep },
# with Inull(n) the mean of the permuted-label refits at n; at n = N the
# Record values are used (Ihat, 200-permutation null mean, half-sampling
# s.e. under the jackknife scaling). If no n qualifies, the arm is "not
# resolved at the available n".
#
# Usage (project root, R 4.4.3): Rscript TRIDENT/learning_curves.R
# Output: TRIDENT/learning/learning_draws.csv, learning_summary.csv,
#         n_evidence.csv
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(glmnet); library(parallel) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI  <- normalizePath(dirname(this_file), winslash = "/")
ROOT <- normalizePath(file.path(TRI, ".."), winslash = "/")
setwd(ROOT); source("PathMI_SUBHO_logistic.R"); source(file.path(TRI, "infold_functions.R"))
OUT <- file.path(TRI, "learning"); dir.create(OUT, showWarnings = FALSE); LOG <- file.path(OUT, "learning_curves.log")
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }

N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "16")); B <- 50L; N_GRID <- c(50L, 100L, 150L, 200L)
V_BETA <- 1.0; NMC <- 2000L; BURN <- 500L; THIN <- 10L; CONTROL_SEED <- 20261119L; Z90 <- qnorm(0.95)
POOL <- "M:/bhadury/mi_results/pool_pass"
MODELS <- c(uni_v2 = "UNI v2", virchow2 = "Virchow2", hoptimus1 = "H-optimus-1", conch_v15 = "CONCH v1.5"); code <- setNames(names(MODELS), MODELS)

# Record arms = each model's largest estimate per state in the grid (as in recoverability_record.R)
g <- as.data.frame(readxl::read_excel(file.path(TRI, "mi summary table.xlsx")))
emb2code <- c("UNI V2" = "uni_v2", "VIRCHOW" = "virchow2", "H-OPTIMUS 1" = "hoptimus1", "CONCH V1.5" = "conch_v15")
EPS <- c(uni_v2 = 0.0163, virchow2 = 0.0163, hoptimus1 = 0.0163, conch_v15 = 0.0215)
g$code <- emb2code[g$Embeddings]; top <- do.call(rbind, lapply(split(g, list(g$Disease, g$target, g$code), drop = TRUE), function(d) d[which.max(d$MI_obs), ]))
arms <- list()
for (i in seq_len(nrow(top))) {
  r <- top[i, ]; co <- tolower(r$Disease); mo <- r$code; dc <- sub("_.*", "", r$key); pl <- sub(".*_", "", r$key)
  pp <- readRDS(file.path(POOL, sprintf("%s_%s", co, mo), "pooled.rds"))
  id <- sprintf("%s|%s|%s|%s", r$Disease, mo, r$target, r$key)
  arms[[id]] <- list(id = id, cohort = r$Disease, model = MODELS[[mo]], state = r$target, decomp = dc, pooling = pl,
                     X = pp$pooled[[pl]], y = as.numeric(pp$y[[r$target]][pp$roster]), g = substr(pp$roster, 1, 12), eps = EPS[[mo]], control = FALSE)
}
ctl_id <- grep("^CRC\\|hoptimus1\\|MSI\\|", names(arms), value = TRUE); stopifnot(length(ctl_id) == 1)
ctl <- arms[[ctl_id]]; set.seed(CONTROL_SEED); ctl$y <- sample(ctl$y)   # identical to usable_infold.R
ctl$id <- sub("^CRC\\|hoptimus1\\|MSI\\|", "CONTROL|hoptimus1|MSI-permuted|", ctl_id); ctl$state <- "MSI (permuted)"; ctl$control <- TRUE
arms[[ctl$id]] <- ctl

job <- function(aid, n, b) {
  # 19 Nov 2026 (QC): the unit is the patient. n patients are drawn (label-stratified at the patient
  # level) with all their slides, and the usable cross-fit groups slides by patient.
  a <- arms[[aid]]; pat <- unique(a$g); py <- a$y[match(pat, a$g)]; N <- length(pat)
  set.seed(5000L + b); sp <- pat[strat_subsample(py, n)]; idx <- which(a$g %in% sp)   # same patients across arms with the same roster
  X <- a$X[idx, , drop = FALSE]; y <- a$y[idx]
  pr <- fit_projection(X, a$decomp, 50L); Z <- unclass(apply_projection(pr, X))
  set.seed(60000L + b); f <- subho_MI_logistic(Z = Z, y = y, v_beta = V_BETA, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
  set.seed(70000L + b); yp <- sample(y)
  fp <- subho_MI_logistic(Z = Z, y = yp, v_beta = V_BETA, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
  cf <- crossfit_infold(X, y, a$decomp, 5L, fold_seed = 80000L + b, model_seed = 90000L + 10L * b, pathmi = FALSE, group = a$g[idx])
  data.frame(id = aid, n = n, N = N, b = b, k_proj = pr$k, MI = f$MI, MI_raw = f$MI.raw, k_star = f$n.selected,
             MI_perm = fp$MI, I_usable = cf$I_ridge, auc = cf$auc_ridge, n_pos = sum(y), stringsAsFactors = FALSE)
}
jobs <- do.call(rbind, lapply(names(arms), function(aid) { N <- length(unique(arms[[aid]]$g)); ng <- sort(unique(c(N_GRID[N_GRID < N], floor(N / 2))))
  expand.grid(aid = aid, n = ng, b = seq_len(B), stringsAsFactors = FALSE) }))
jobs <- jobs[order(-jobs$n), ]                                   # the largest subsamples first
logf("=== learning curves | %d arms | %d jobs | B = %d | %d workers ===", length(arms), nrow(jobs), B, N_WORKERS)

cl <- makeCluster(N_WORKERS); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(glmnet) }))
clusterExport(cl, c("arms", "jobs", "job", "V_BETA", "NMC", "BURN", "THIN", "fit_projection", "apply_projection", "strat_folds", "strat_subsample",
                    "clip01", "loglik_gain", "auc", "ridge_eta", "crossfit_infold", "EIG_TOL",
                    "sagnik.horseshoe.logistic", "lpy.x.logistic", "selection.logistic", "subho_MI_logistic", "sigmoid", "clip_eta"))
t0 <- Sys.time(); out <- vector("list", nrow(jobs)); done <- 0L
for (ch in split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / 160))) {
  o <- parLapply(cl, ch, function(j) tryCatch(job(jobs$aid[j], jobs$n[j], jobs$b[j]), error = function(e) data.frame(id = jobs$aid[j], n = jobs$n[j], b = jobs$b[j], err = conditionMessage(e))))
  for (i in seq_along(ch)) out[[ch[i]]] <- o[[i]]
  done <- done + length(ch); el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  logf("  %4d / %d jobs | %.1f min | ETA %.1f min", done, nrow(jobs), el, el / done * (nrow(jobs) - done))
}
errs <- Filter(function(d) !is.null(d$err), out)
if (length(errs)) logf("[WARN] %d job(s) errored: %s", length(errs), paste(head(unique(sapply(errs, `[[`, "err"))), collapse = " | "))
draws <- do.call(rbind, Filter(function(d) is.null(d$err), out))
write.csv(draws, file.path(OUT, "learning_draws.csv"), row.names = FALSE)

# ---- full-cohort control fit (the Record arms use their stored values) ----------------
a <- arms[[ctl$id]]; set.seed(42L)
f_ctl <- subho_MI_logistic(Z = unclass(apply_projection(fit_projection(a$X, a$decomp, 50L), a$X)), y = a$y, v_beta = V_BETA, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
logf("negative control, full cohort: Ihat = %.4f (raw %.4f, k* = %d)", f_ctl$MI, f_ctl$MI.raw, f_ctl$n.selected)
saveRDS(list(MI = f_ctl$MI, MI_raw = f_ctl$MI.raw, k = f_ctl$n.selected), file.path(OUT, "control_full_fit.rds"))
logf("=== done in %.1f min ===", as.numeric(difftime(Sys.time(), t0, units = "mins")))

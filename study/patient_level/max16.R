################################################################
# max16.R  (27 Nov 2026, final-edits check 1)
#
# Joint-permutation max-T over all 16 pooling x projection
# configurations of each model-state pair, from the stored
# patient-level permutation fits (pl_grid.R, all_draws.rds).
#
#   T_b   = max over the 16 configurations of the estimate under
#           permutation b (b = 1..200);
#   T_obs = max over the 16 observed estimates;
#   p_max16 = (1 + #{b : T_b >= T_obs}) / 201.
#
# Alignment verification (the max-T null is valid only if permutation
# b relabels the same patients in every configuration):
#   A1 every configuration holds b = 0..200 exactly once;
#   A2 the shared permutation matrix is reconstructed from pl_grid.R's
#      seed (set.seed(42) per cohort x model, drawn before any other
#      random call) and is a set of permutations of the patient roster;
#   A3 exact refits: for two configurations per model-state pair (the
#      selected one and a second pooling rule) and two permutations
#      (b = 1, 137), the estimate is refitted from the pooled matrices
#      with the reconstructed labels and the stored seed 100000 + b and
#      compared with the stored value (agreement to < 1e-10 confirms
#      that the stored null used exactly these labels);
#   A4 cross-configuration correlation of the null estimates at matched
#      b versus mismatched b (shared labels induce positive dependence).
#
# Usage (project root, R 4.4.3):
#   Rscript TRIDENT/patient_level/max16.R            (A1, A2, A4 + p values)
#   Rscript TRIDENT/patient_level/max16.R refit      (adds A3, ~10 min)
# Output: TRIDENT/patient_level/max16/max16.csv, alignment.csv, refit.csv, max16.log
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(parallel) })
args <- commandArgs(trailingOnly = TRUE); DO_REFIT <- length(args) >= 1 && args[1] == "refit"
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
PLD <- normalizePath(dirname(this_file), winslash = "/"); TRI <- dirname(PLD); ROOT <- dirname(TRI)
setwd(ROOT)
OUT <- file.path(PLD, "max16"); dir.create(OUT, showWarnings = FALSE)
LOG <- file.path(OUT, "max16.log"); cat("", file = LOG)
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }

draws <- readRDS(file.path(PLD, "patient", "all_draws.rds"))
meta <- readRDS(file.path(PLD, "patient", "meta_runs.rds"))$meta_runs
N_PERM <- 200L; SEED <- 42L
draws$tag <- sub("\\|.*", "", draws$cid); draws$key <- sub(".*\\|", "", draws$cid)

rows <- list(); al <- list()
for (tag in names(meta)) {
  mr <- meta[[tag]]; n <- length(mr$units)
  summ <- read.csv(file.path(PLD, "patient", tag, "mi_summary.csv"), stringsAsFactors = FALSE)
  # A2: reconstruct the shared permutation matrix exactly as pl_grid.R drew it
  set.seed(SEED); perms <- sapply(seq_len(N_PERM), function(b) sample.int(n))
  a2 <- all(apply(perms, 2, function(p) identical(sort(p), seq_len(n))))
  for (tg in names(mr$ylist)) {
    d <- draws[draws$tag == tag & draws$target == tg, ]
    keys <- sort(unique(d$key))
    # A1: completeness and uniqueness of b in every configuration
    a1 <- all(sapply(keys, function(k) { b <- d$b[d$key == k]; length(b) == N_PERM + 1L && setequal(b, 0:N_PERM) }))
    M <- sapply(keys, function(k) { x <- d[d$key == k & d$b > 0, ]; x$MI[order(x$b)] })    # 200 x 16, row b
    obs <- sapply(keys, function(k) d$MI[d$key == k & d$b == 0])
    Tb <- apply(M, 1, max); Tobs <- max(obs); win <- keys[which.max(obs)]
    p16 <- (1 + sum(Tb >= Tobs)) / (1 + N_PERM)
    # the same with the conditioning filter used for the four-candidate max-T (eig_cond < 1e10)
    ok <- summ$key[summ$target == tg & summ$eig_cond < 1e10]; ok <- intersect(keys, ok)
    Tb_f <- apply(M[, ok, drop = FALSE], 1, max); Tobs_f <- max(obs[ok])
    p16_f <- (1 + sum(Tb_f >= Tobs_f)) / (1 + N_PERM)
    # A4: correlation of null estimates across configurations at matched vs mismatched b
    C <- suppressWarnings(cor(M)); r_match <- median(C[upper.tri(C)], na.rm = TRUE)
    set.seed(7); sh <- sapply(seq_len(ncol(M)), function(j) M[sample.int(N_PERM), j]); Cs <- suppressWarnings(cor(sh))
    r_shuf <- median(Cs[upper.tri(Cs)], na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(tag = tag, cohort = toupper(mr$cohort), model = mr$model, target = tg, n = n,
      n_configs = length(keys), winner = win, T_obs = Tobs, T_null_mean = mean(Tb), T_null_q95 = unname(quantile(Tb, .95)),
      T_null_max = max(Tb), n_exceed = sum(Tb >= Tobs), p_max16 = p16, n_configs_filtered = length(ok), p_max16_filtered = p16_f)
    al[[length(al) + 1]] <- data.frame(tag = tag, target = tg, A1_complete = a1, A2_perms_valid = a2,
      A4_median_cor_matched = r_match, A4_median_cor_shuffled = r_shuf)
    logf("[%s %s] n %d | %d configs | winner %s T_obs %.4f | null max-16 mean %.4f q95 %.4f max %.4f | p_max16 %.4f (filtered %d: %.4f) | A1 %s A2 %s | cor matched %.2f shuffled %.2f",
         tag, tg, n, length(keys), win, Tobs, mean(Tb), quantile(Tb, .95), max(Tb), p16, length(ok), p16_f, a1, a2, r_match, r_shuf)
  }
}
res <- do.call(rbind, rows); write.csv(res, file.path(OUT, "max16.csv"), row.names = FALSE)
write.csv(do.call(rbind, al), file.path(OUT, "alignment.csv"), row.names = FALSE)

if (DO_REFIT) {
  # A3: exact refits of stored permutation fits (pl_grid.R build, patient mode)
  source("PathMI_SUBHO_logistic.R"); source(file.path(TRI, "projection_functions.R"))
  POOL <- "M:/bhadury/mi_results/pool_pass"; D <- 50L; V_BETA <- 1.0; NMC <- 2000L; BURN <- 500L; THIN <- 10L
  rec <- read.csv(file.path(TRI, "record", "recoverability_record.csv"), stringsAsFactors = FALSE)
  rec <- rec[rec$record_type == "measurement", ]
  jobs <- list()
  for (i in seq_len(nrow(rec))) {
    r <- rec[i, ]; tag <- paste(tolower(r$cohort), r$trident_key, sep = "_")
    sel <- paste0(r$projection, "_", r$pooling); alt_pool <- setdiff(c("mean", "attention"), r$pooling)[1]
    for (key in c(sel, paste0(r$projection, "_", alt_pool))) for (b in c(1L, 137L))
      jobs[[length(jobs) + 1]] <- list(tag = tag, target = r$state, key = key, b = b)
  }
  build <- function(tag) {
    pp <- readRDS(file.path(POOL, tag, "pooled.rds")); slide <- pp$roster; unit <- substr(slide, 1, 12); units <- unique(unit)
    collapse <- function(X) { S <- rowsum(X, unit, reorder = FALSE); S / as.vector(table(unit)[rownames(S)]) }
    raw <- lapply(pp$pooled, collapse)
    ylist <- lapply(pp$y, function(v) { v <- v[slide]; setNames(as.numeric(v[match(units, unit)]), units) })
    list(raw = raw, ylist = ylist, units = units)
  }
  one <- function(j) {
    jb <- jobs[[j]]; B <- build(jb$tag); n <- length(B$units)
    stopifnot(identical(B$units, meta[[jb$tag]]$units))
    set.seed(42L); perms <- sapply(seq_len(200L), function(b) sample.int(n))
    dc <- sub("_.*", "", jb$key); pl <- sub(".*_", "", jb$key)
    Z <- unclass(project(B$raw[[pl]], dc, D)$Z_scaled)
    set.seed(100000L + jb$b)
    f <- subho_MI_logistic(Z = Z, y = B$ylist[[jb$target]][perms[, jb$b]], v_beta = V_BETA, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
    data.frame(tag = jb$tag, target = jb$target, key = jb$key, b = jb$b, MI_refit = f$MI)
  }
  nw <- min(length(jobs), as.integer(Sys.getenv("N_WORKERS", "16")))
  cl <- makeCluster(nw); on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, suppressPackageStartupMessages({ library(MASS); library(BayesLogit) }))
  clusterExport(cl, c("jobs", "build", "meta", "POOL", "D", "V_BETA", "NMC", "BURN", "THIN", "ROOT", "TRI"), envir = environment())
  clusterEvalQ(cl, { setwd(ROOT); source("PathMI_SUBHO_logistic.R"); source(file.path(TRI, "projection_functions.R")) })
  rf <- do.call(rbind, parLapply(cl, seq_along(jobs), one))
  rf$MI_stored <- mapply(function(tag, tg, key, b) draws$MI[draws$cid == paste(tag, key, sep = "|") & draws$target == tg & draws$b == b],
                         rf$tag, rf$target, rf$key, rf$b)
  rf$abs_diff <- abs(rf$MI_refit - rf$MI_stored)
  write.csv(rf, file.path(OUT, "refit.csv"), row.names = FALSE)
  logf("A3 exact refits: %d fits | max |refit - stored| = %.3g | all < 1e-10: %s", nrow(rf), max(rf$abs_diff), all(rf$abs_diff < 1e-10))
}

################################################################
# robustness_subsampling.R
#
# Patient-level sampling uncertainty by subsampling without
# replacement (half-sampling), for the same 57 arms as the bootstrap
# in robustness_prior_bootstrap.R.
#
# Why not the bootstrap with replacement: a resample duplicates about
# 37 percent of the patients, and a duplicated patient enters the
# marginal likelihood twice while the Occam penalty of the walk grows
# only with log n, so the walk admits more coordinates (k* reached the
# full 50 in the upper tail of the resamples) and the resampled
# estimates sit 0.02 to 0.08 nats above the point estimate
# (robustness/bootstrap_summary.csv). Subsampling without replacement
# has no duplicates. With m = floor(n/2) patients per subsample and
# B = 200 subsamples, the standard error at n is the delete-d
# jackknife with d = n - m (Shao and Wu, Ann. Statist. 1989):
#     s.e._n = s.d.(estimate_m) * sqrt(m / (n - m)),
# which is the s.d. of the half-sample estimates when m = n/2. (Until
# 19 Nov 2026 this used sqrt(m / n), the rate for m/n -> 0, which at
# m = n/2 omits the finite-population factor and is sqrt(2) too small;
# the old summary is kept as subsampling_summary_sqrt_m_over_n.csv.)
# The 90 percent interval is the point estimate +/- 1.645 s.e._n.
# Set SUMMARY_ONLY=1 to re-summarise the stored draws without refitting. The subsample mean
# is also reported: it lies below the point estimate because the
# estimator is more conservative at smaller n (simulation), which is
# why the subsample distribution is used for its spread only.
#
# The same seed draws the same patients for arms with the same roster,
# so contrasts between such arms are paired through the draws.
#
# Usage (project root, R 4.4.3): Rscript TRIDENT/robustness_subsampling.R
# Output: TRIDENT/robustness/{subsampling_summary,subsampling_draws}.csv
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(parallel); library(readxl) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI  <- normalizePath(dirname(this_file), winslash = "/")
ROOT <- normalizePath(file.path(TRI, ".."), winslash = "/")
setwd(ROOT); source("PathMI_SUBHO_logistic.R")
OUT <- file.path(TRI, "robustness"); LOG <- file.path(OUT, "robustness_subsampling.log")
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }
N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "24")); NMC <- 2000L; BURN <- 500L; THIN <- 10L; B_SUB <- 200L; FRAC <- 0.5
MODELS <- c(uni_v2 = "UNI v2", virchow2 = "Virchow2", hoptimus1 = "H-optimus-1", conch_v15 = "CONCH v1.5")
harm <- function(e) c("UNI V2" = "uni_v2", "VIRCHOW" = "virchow2", "H-OPTIMUS 1" = "hoptimus1", "CONCH V1.5" = "conch_v15")[[e]]
res_cache <- list()
load_cell <- function(cohort, model, target, key) {
  f <- if (cohort == "CRC") sprintf("crc_kpca_subho_%s/kpca_subho_results_crc_%s.rds", model, model) else sprintf("TCGA_PRAD_kpca_subho_%s/kpca_subho_results_prad_%s.rds", model, model)
  if (is.null(res_cache[[f]])) res_cache[[f]] <<- readRDS(file.path(TRI, f))
  a <- res_cache[[f]][[key]]; roster <- rownames(a$Z_scaled)
  y <- if (cohort == "CRC") a$y[roster] else a$y_list[[target]][roster]; det <- a$result$details[[target]]
  list(id = sprintf("%s|%s|%s|%s", cohort, model, target, key), cohort = cohort, model = MODELS[[model]], target = target, key = key, Z = unclass(a$Z_scaled), y = as.numeric(y),
       MI_obs = det$MI, null_mean = mean(det$MI.null, na.rm = TRUE), null_sd = sd(det$MI.null, na.rm = TRUE), k_star = det$n.selected)
}
best <- as.data.frame(read_excel(file.path(TRI, "best kpca summary.xlsx")))
cfg <- unique(rbind(
  data.frame(cohort = best$Disease, model = sapply(best$Embeddings, harm), target = best$target, key = best$key, arm = "selected", stringsAsFactors = FALSE),
  expand.grid(cohort = "CRC", model = names(MODELS), target = "MSI", key = "pca_mean", arm = "linear_mean", stringsAsFactors = FALSE),
  expand.grid(cohort = "PRAD", model = names(MODELS), target = c("ERG", "ETS"), key = "pca_mean", arm = "linear_mean", stringsAsFactors = FALSE)))
dup <- duplicated(cfg[, 1:4]) | duplicated(cfg[, 1:4], fromLast = TRUE); cfg$arm[dup] <- "selected+linear_mean"; cfg <- cfg[!duplicated(cfg[, 1:4]), ]
cells <- lapply(seq_len(nrow(cfg)), function(i) with(cfg[i, ], load_cell(cohort, model, target, key))); names(cells) <- sapply(cells, `[[`, "id")
logf("=== subsampling | %d arms | m = floor(%.2f n) | B = %d | %d workers ===", length(cells), FRAC, B_SUB, N_WORKERS)

SUMMARY_ONLY <- Sys.getenv("SUMMARY_ONLY", "0") == "1"
t0 <- Sys.time()
if (!SUMMARY_ONLY) {
job <- function(cid, block) {
  cell <- cells[[cid]]; Z <- cell$Z; y <- cell$y; n <- length(y); m <- floor(FRAC * n)
  bs <- ((block - 1L) * 25L + 1L):(block * 25L)
  do.call(rbind, lapply(bs, function(b) { set.seed(777L + b); idx <- sample.int(n, m); set.seed(30000L + b)
    r <- subho_MI_logistic(Z = Z[idx, , drop = FALSE], y = y[idx], v_beta = 1, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
    data.frame(id = cid, b = b, m = m, MI = r$MI, MI_raw = r$MI.raw, k = r$n.selected, prevalence = mean(y[idx])) }))
}
jobs <- expand.grid(cid = names(cells), block = seq_len(B_SUB / 25L), stringsAsFactors = FALSE); jobs <- jobs[sample.int(nrow(jobs)), ]
cl <- makeCluster(N_WORKERS); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(MASS); library(BayesLogit) }))
clusterExport(cl, c("cells", "jobs", "job", "FRAC", "NMC", "BURN", "THIN", "sagnik.horseshoe.logistic", "lpy.x.logistic", "selection.logistic", "subho_MI_logistic", "sigmoid", "clip_eta"))
t0 <- Sys.time(); out <- vector("list", nrow(jobs)); done <- 0L
for (ch in split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / 96))) {
  o <- parLapply(cl, ch, function(j) tryCatch(job(jobs$cid[j], jobs$block[j]), error = function(e) data.frame(id = jobs$cid[j], err = conditionMessage(e))))
  for (i in seq_along(ch)) out[[ch[i]]] <- o[[i]]
  done <- done + length(ch); el <- as.numeric(difftime(Sys.time(), t0, units = "mins")); logf("  %4d / %d jobs | %.1f min | ETA %.1f min", done, nrow(jobs), el, el / done * (nrow(jobs) - done))
}
errs <- Filter(function(d) !is.null(d$err), out); if (length(errs)) logf("[WARN] %d errored: %s", length(errs), paste(head(unique(sapply(errs, `[[`, "err"))), collapse = " | "))
draws <- do.call(rbind, Filter(function(d) is.null(d$err), out)); write.csv(draws, file.path(OUT, "subsampling_draws.csv"), row.names = FALSE)
}
draws <- read.csv(file.path(OUT, "subsampling_draws.csv"), stringsAsFactors = FALSE)
meta <- do.call(rbind, lapply(seq_along(cells), function(i) { c <- cells[[i]]; data.frame(id = c$id, arm = cfg$arm[i], cohort = c$cohort, model = c$model, target = c$target, key = c$key, n = length(c$y),
  MI_obs = c$MI_obs, null_mean = c$null_mean, null_sd = c$null_sd, k_star = c$k_star, stringsAsFactors = FALSE) }))
summ <- do.call(rbind, lapply(split(draws, draws$id), function(d) { n <- meta$n[meta$id == d$id[1]]; m <- d$m[1]; se <- sd(d$MI) * sqrt(m / (n - m))
  data.frame(id = d$id[1], B = nrow(d), m = m, sub_mean = mean(d$MI), sub_sd = sd(d$MI), se_n = se, sub_lo05 = quantile(d$MI, .05), sub_hi95 = quantile(d$MI, .95),
             k_median = median(d$k), k_lo = quantile(d$k, .05), k_hi = quantile(d$k, .95), frac_zero = mean(d$MI == 0)) }))
summ <- merge(meta, summ, by = "id")
summ$ci90_lo <- pmax(summ$MI_obs - 1.645 * summ$se_n, 0); summ$ci90_hi <- summ$MI_obs + 1.645 * summ$se_n; summ$sub_bias <- summ$sub_mean - summ$MI_obs
summ <- summ[order(summ$cohort, summ$target, summ$model, summ$key), ]
write.csv(summ, file.path(OUT, "subsampling_summary.csv"), row.names = FALSE)
logf("=== done in %.1f min ===", as.numeric(difftime(Sys.time(), t0, units = "mins")))
for (i in seq_len(nrow(summ))) with(summ[i, ], logf("%-45s MI %.4f | half-sample mean %.4f (bias %+.4f) | s.e._n %.4f | 90%% [%.4f, %.4f] | k* %d (half-sample median %g, 90%% %g-%g)", id, MI_obs, sub_mean, sub_bias, se_n, ci90_lo, ci90_hi, k_star, k_median, k_lo, k_hi))
logf("s.e. summary: median %.4f, range %.4f to %.4f | selected arms median %.4f | median s.e. in resolution units (0.0095): %.1f | median half-sample bias %+.4f",
     median(summ$se_n), min(summ$se_n), max(summ$se_n), median(summ$se_n[grepl("selected", summ$arm)]), median(summ$se_n) / 0.0095, median(summ$sub_bias))

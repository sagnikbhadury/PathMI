################################################################
# resolution_units.R
#
# Every major comparison of the Results expressed in measurement
# units: the contrast Delta I divided by (i) the reproducibility
# bound epsilon_rep (0.0095 nats over the three vision-only models,
# 0.0202 with CONCH v1.5), the resolution below which two estimates
# cannot be told apart, and (ii) the half-sampling standard error
# of the contrast (sampling uncertainty). The two denominators answer
# different questions: epsilon_rep is what the estimator itself can
# resolve, the half-sampling s.e. is how far the cohort could move it.
#
# Contrasts between selected arms of two models on the same cohort
# are paired through the half-sample draws when the two rosters are
# identical (the same seed draws the same patients); otherwise the
# s.e. of the difference is the root sum of squares (conservative).
#
# Inputs: the merged tables, robustness/subsampling_draws.csv and
# subsampling_summary.csv, usable/usable_infold.csv (in-fold ridge readout).
# Half-sampling s.e.: delete-d jackknife, s.d. x sqrt(m / (n - m)) (Shao and Wu 1989).
# Usage: Rscript TRIDENT/resolution_units.R
# Output: TRIDENT/robustness/resolution_units.csv (+ .tex rows)
################################################################
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/")
# 27 Nov 2026: numerical tolerance from tolerance.R (replicate refit + sign-flip refit; the null covariate block is a specification check, not in it)
TOL <- read.csv(file.path(TRI, "robustness", "tolerance.csv"), stringsAsFactors = FALSE)
EPS_V <- TOL$eps_tol_rounded[TOL$group == "vision-only"]; EPS_C <- TOL$eps_tol_rounded[TOL$group == "CONCH v1.5"]
BB <- c("UNI v2", "Virchow2", "H-optimus-1", "CONCH v1.5")
harm <- function(e) case_when(grepl("^UNI", e) ~ BB[1], grepl("^VIRCHOW", e) ~ BB[2], grepl("^H-OPT", e) ~ BB[3], grepl("^CONCH", e) ~ BB[4])
code <- c("UNI v2" = "uni_v2", "Virchow2" = "virchow2", "H-optimus-1" = "hoptimus1", "CONCH v1.5" = "conch_v15")
g <- read_excel(file.path(TRI, "mi summary table.xlsx")) %>% mutate(bb = harm(Embeddings), task = paste0(Disease, "/", target), decomp = tolower(decomp))
best <- read_excel(file.path(TRI, "best kpca summary.xlsx")) %>% mutate(bb = harm(Embeddings), task = paste0(Disease, "/", target))
inc <- read_excel(file.path(TRI, "incremental mi summary.xlsx")) %>% mutate(bb = harm(Embeddings), task = paste0(Disease, "/", target))
bs <- read.csv(file.path(TRI, "robustness", "subsampling_summary.csv"), stringsAsFactors = FALSE)
bd <- read.csv(file.path(TRI, "robustness", "subsampling_draws.csv"), stringsAsFactors = FALSE)
us <- read.csv(file.path(TRI, "usable", "usable_infold.csv"), stringsAsFactors = FALSE)

rows <- list()
add <- function(group, comparison, delta, eps, se = NA_real_, note = "") rows[[length(rows) + 1]] <<- data.frame(group = group, comparison = comparison, delta = delta, eps = eps,
  units_rep = delta / eps, se_half = se, units_se = ifelse(is.na(se), NA, delta / se), note = note, stringsAsFactors = FALSE)

# half-sample draws of a cell and the paired/unpaired s.e. of a difference
draws_of <- function(id) bd$MI[bd$id == id][order(bd$b[bd$id == id])]
scale_of <- function(id) { m <- bd$m[bd$id == id][1]; n <- bs$n[bs$id == id]; sqrt(m / (n - m)) }   # half-sampling, delete-d jackknife: s.e._n = s.d._m * sqrt(m/(n-m))
se_diff <- function(id1, id2, paired) { a <- draws_of(id1); b <- draws_of(id2); if (paired) sd(a - b) * scale_of(id1) else sqrt((sd(a) * scale_of(id1))^2 + (sd(b) * scale_of(id2))^2) }
roster_same <- function(cohort, m1, m2) {
  f <- function(m) rownames(readRDS(file.path(TRI, if (cohort == "CRC") sprintf("crc_kpca_subho_%s/kpca_subho_results_crc_%s.rds", m, m) else sprintf("TCGA_PRAD_kpca_subho_%s/kpca_subho_results_prad_%s.rds", m, m)))$pca_mean$Z_scaled)
  identical(f(m1), f(m2)) }
top_arm <- function(tk, b) best %>% filter(task == tk, bb == b) %>% slice_max(MI_obs, n = 1, with_ties = FALSE)
cell_id <- function(tk, b, key) sprintf("%s|%s|%s|%s", sub("/.*", "", tk), code[[b]], sub(".*/", "", tk), key)

# 1. tier contrasts: median paired difference over the 48 matched blocks
w <- g %>% select(task, decomp, agg_method, bb, MI_obs) %>% pivot_wider(names_from = bb, values_from = MI_obs)
for (pr in list(c(BB[3], BB[1]), c(BB[3], BB[2]), c(BB[1], BB[2]), c(BB[3], BB[4]), c(BB[1], BB[4]), c(BB[2], BB[4]))) {
  d <- median(w[[pr[1]]] - w[[pr[2]]]); eps <- if (BB[4] %in% pr) EPS_C else EPS_V
  add("Model tiers (48 blocks)", sprintf("%s minus %s", pr[1], pr[2]), d, eps, note = "median paired difference over 48 blocks") }

# 2. within-state H-optimus-1 vs UNI v2 (median over 16 cells) and the selected-arm contrast with half-sampling s.e.
for (tk in c("CRC/MSI", "PRAD/ERG", "PRAD/ETS")) {
  d16 <- w %>% filter(task == tk); d <- median(d16[[BB[3]]] - d16[[BB[1]]])
  add("Within-state winner (16 cells)", sprintf("%s: H-optimus-1 minus UNI v2", tk), d, EPS_V, note = "median over 16 matched cells")
  a1 <- top_arm(tk, BB[3]); a2 <- top_arm(tk, BB[1]); paired <- roster_same(sub("/.*", "", tk), "hoptimus1", "uni_v2")
  se <- se_diff(cell_id(tk, BB[3], a1$key), cell_id(tk, BB[1], a2$key), paired)
  add("Within-state winner (best arms)", sprintf("%s: best H-optimus-1 arm minus best UNI v2 arm", tk), a1$MI_obs - a2$MI_obs, EPS_V, se, note = if (paired) "paired half-samples, identical rosters" else "unpaired half-sampling s.e.")
}
# 3. projections and pooling: median paired difference over blocks (other three factors crossed)
wd <- g %>% select(task, bb, agg_method, decomp, MI_obs) %>% pivot_wider(names_from = decomp, values_from = MI_obs)
for (pr in list(c("cosine", "pca"), c("cosine", "rbf"), c("cosine", "poly"), c("pca", "poly"))) add("Projection (48 blocks)", sprintf("%s minus %s kernel PCA", pr[1], pr[2]), median(wd[[pr[1]]] - wd[[pr[2]]]), EPS_V, note = "median paired difference over 48 blocks")
wp <- g %>% select(task, bb, decomp, agg_method, MI_obs) %>% pivot_wider(names_from = agg_method, values_from = MI_obs)
for (pr in list(c("moment", "mean"), c("norm", "mean"), c("attention", "mean"))) add("Pooling (48 blocks)", sprintf("%s minus mean pooling", pr[1]), median(wp[[pr[1]]] - wp[[pr[2]]]), EPS_V, note = "median paired difference over 48 blocks")
# 4. default configuration cost (largest shortfall across states)
cfg <- g %>% group_by(task) %>% mutate(short = max(MI_obs) - MI_obs) %>% group_by(bb, decomp, agg_method) %>% summarise(cost = max(short), .groups = "drop")
for (k in list(c(BB[3], "pca", "moment"), c(BB[3], "cosine", "moment"))) { cc <- cfg %>% filter(bb == k[1], decomp == k[2], agg_method == k[3]) %>% pull(cost)
  add("Default configuration", sprintf("%s, %s, %s pooling: worst-case cost", k[1], k[2], k[3]), cc, EPS_V, note = "largest shortfall from the state's best") }
# 5. covariates
ci <- inc %>% filter(Disease == "CRC", bb != BB[4]); add("Covariates", "CRC/MSI: covariate-only (size and stage)", mean(ci$covar_only), EPS_V, note = "mean over 12 vision-only arms")
pc <- read.csv(file.path(TRI, "prad_incremental_batch_purity", "prad_incremental_batch_purity.csv"), stringsAsFactors = FALSE)
if ("covar_batch" %in% names(pc)) for (tg in c("ERG", "ETS")) add("Covariates", sprintf("PRAD/%s: batch-only covariate block", tg), mean(pc$covar_batch[pc$target == tg], na.rm = TRUE), EPS_V, note = "mean over arms")
# 6. selection inflation and peak-over-null
add("Selection", "max-T selection inflation (about one null s.d.)", mean(g$null_sd), EPS_V, note = "mean null s.d. over 192 cells")
for (tk in c("CRC/MSI", "PRAD/ERG", "PRAD/ETS")) { d <- g %>% filter(task == tk) %>% slice_max(MI_obs, n = 1, with_ties = FALSE)
  add("Signal over null", sprintf("%s: largest estimate minus its null mean", tk), d$MI_obs - d$null_mean, EPS_V, note = sprintf("%s, %s, %s", d$bb, d$decomp, d$agg_method)) }
# 7. present minus usable, best arm per state (the extraction gap)
for (tk in c("CRC/MSI", "PRAD/ERG", "PRAD/ETS")) { d <- g %>% filter(task == tk) %>% slice_max(MI_obs, n = 1, with_ties = FALSE)
  u <- us %>% filter(id == cell_id(tk, d$bb, d$key))
  if (nrow(u)) add("Present minus usable", sprintf("%s: best arm, PathMI minus cross-fitted usable (ridge readout, in-fold projection)", tk), d$MI_obs - u$I_usable_infold, EPS_V, u$I_usable_infold_sd, note = "s.e. = s.d. over 10 cross-fitting repeats") }
# 8. half-sampling s.e. of the best arm of each state, in resolution units
for (tk in c("CRC/MSI", "PRAD/ERG", "PRAD/ETS")) { d <- g %>% filter(task == tk) %>% slice_max(MI_obs, n = 1, with_ties = FALSE)
  b <- bs %>% filter(id == cell_id(tk, d$bb, d$key)); if (nrow(b)) add("Sampling uncertainty", sprintf("%s: half-sampling s.e. of the best arm", tk), b$se_n, EPS_V, note = sprintf("90%% interval %.3f to %.3f", b$ci90_lo, b$ci90_hi)) }

res <- do.call(rbind, rows) %>% mutate(resolved = ifelse(abs(units_rep) >= 1, "resolved", "below resolution"))
write.csv(res, file.path(TRI, "robustness", "resolution_units.csv"), row.names = FALSE)
print(as.data.frame(res %>% select(comparison, delta, units_rep, se_half, units_se, resolved, note)), digits = 3, right = FALSE)
cat(sprintf("\nhalf-sampling s.e. over the 48 selected arms: median %.4f (range %.4f to %.4f) = %.2f resolution units (median)\n",
            median(bs$se_n[grepl("selected", bs$arm)]), min(bs$se_n[grepl("selected", bs$arm)]), max(bs$se_n[grepl("selected", bs$arm)]), median(bs$se_n[grepl("selected", bs$arm)]) / EPS_V))

################################################################
# usable_infold_stats.R
#
# The 192-cell statements of the usable-recoverability subsection,
# computed for the stored fixed-projection ridge readout (to check they
# reproduce the earlier text) and for the in-fold readout
# (usable_infold.R), side by side.
# Output: TRIDENT/usable/usable_infold_stats.txt
################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/")
u <- read.csv(file.path(TRI, "usable", "usable_recoverability.csv"), stringsAsFactors = FALSE)
ui <- read.csv(file.path(TRI, "usable", "usable_infold.csv"), stringsAsFactors = FALSE)
d <- u %>% left_join(ui %>% select(id, I_infold = I_usable_infold, auc_infold), by = "id") %>% mutate(task = paste0(cohort, "/", target))
sink(file.path(TRI, "usable", "usable_infold_stats.txt"), split = TRUE)
for (v in c("I_usable_ridge", "I_infold")) {
  x <- d[[v]]
  cat(sprintf("\n== %s ==\nabove null band (null mean + 2 s.d.) in %d of %d | r with PathMI %.2f (Spearman %.2f) | exceeds PathMI in %d | median usable/present %.2f\n",
              v, sum(x > d$null_mean + 2 * d$null_sd), nrow(d), cor(x, d$MI_obs), cor(x, d$MI_obs, method = "spearman"), sum(x > d$MI_obs), median(x / d$MI_obs)))
  w <- d %>% mutate(val = x, block = paste(task, decomp, pooling)) %>% select(block, model, val) %>% pivot_wider(names_from = model, values_from = val)
  ft <- friedman.test(as.matrix(w[, -1])); cat(sprintf("Friedman over %d blocks: chi2 = %.1f, p = %.1e | mean ranks: %s\n", nrow(w), ft$statistic, ft$p.value,
    paste(sprintf("%s %.2f", colnames(w)[-1], colMeans(t(apply(-as.matrix(w[, -1]), 1, rank)))), collapse = ", ")))
  for (tk in unique(d$task)) { s <- d %>% mutate(val = x) %>% filter(task == tk)
    b <- s %>% group_by(model) %>% slice_max(MI_obs, n = 1, with_ties = FALSE) %>% ungroup()
    cat(sprintf("  %s: best-arm usable by model: %s | best overall arm: %s %.3f\n", tk, paste(sprintf("%s %.3f", b$model, b$val), collapse = ", "),
                s$model[which.max(s$MI_obs)], s$val[which.max(s$MI_obs)])) }
}
cat(sprintf("\nin-fold minus fixed: median %+.4f, IQR %+.4f to %+.4f, range %+.4f to %+.4f; |diff| < eps_rep (0.0095) in %d of 192\n",
            median(d$I_infold - d$I_usable_ridge), quantile(d$I_infold - d$I_usable_ridge, .25), quantile(d$I_infold - d$I_usable_ridge, .75),
            min(d$I_infold - d$I_usable_ridge), max(d$I_infold - d$I_usable_ridge), sum(abs(d$I_infold - d$I_usable_ridge) < 0.0095)))
cat(sprintf("best arms in-fold AUC: %s\n", paste(sprintf("%s %.2f", d$task, d$auc_infold)[d$id %in% c("CRC|hoptimus1|MSI|cosine_moment", "PRAD|uni_v2|ERG|cosine_norm", "PRAD|hoptimus1|ETS|cosine_moment")], collapse = "; ")))
sink()

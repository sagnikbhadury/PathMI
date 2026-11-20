T <- "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/entropy/Research - DOE Horseshoe/TRIDENT"
u <- read.csv(file.path(T, "usable", "usable_infold.csv")); s <- u[!u$control, ]
cat(sprintf("in-fold minus fixed: median %+.4f IQR %+.4f..%+.4f range %+.4f..%+.4f | |diff| < 0.0163 in %d of 192\n", median(s$diff_infold_minus_fixed),
  quantile(s$diff_infold_minus_fixed, .25), quantile(s$diff_infold_minus_fixed, .75), min(s$diff_infold_minus_fixed), max(s$diff_infold_minus_fixed), sum(abs(s$diff_infold_minus_fixed) < 0.0163)))
r <- read.csv(file.path(T, "record", "recoverability_record.csv")); r <- r[r$record_type == "measurement", ]
r$pm <- r$I_usable_pathmi / r$I_hat; print(tapply(r$pm, r$state, median)); print(r[, c("state", "model", "I_hat", "I_usable", "auc_oof", "I_usable_pathmi", "I_adjusted")], digits = 3)
cat(sprintf("usable/present over 12: median %.2f range %.2f..%.2f | adjusted/present median %.2f range %.2f..%.2f\n", median(r$I_usable / r$I_hat), min(r$I_usable / r$I_hat), max(r$I_usable / r$I_hat),
  median(r$I_adjusted / r$I_hat), min(r$I_adjusted / r$I_hat), max(r$I_adjusted / r$I_hat)))
# PathMI readout share over 192 cells by state (in-fold)
print(tapply(s$I_usable_pathmi_infold / s$MI_obs, s$target, median)); cat("negative PathMI-readout cells:", sum(s$I_usable_pathmi_infold < 0), "\n")
# batch composition at patient level
pp <- readRDS("M:/bhadury/mi_results/pool_pass/prad_hoptimus1/pooled.rds"); m <- read.csv("M:/shared_data/TCGA_PRAD/orig_annotations.csv"); m <- m[!duplicated(m$slide), ]
b <- m$Batch[match(pp$slide_roster[match(pp$roster, substr(pp$slide_roster, 1, 12))], m$slide)]; y <- pp$y$ERG[pp$roster]
tb <- table(b); cat("largest batches:\n"); print(head(sort(tb, decreasing = TRUE), 4)); cat("batch 389: n", sum(b == 389), "ERG prev", mean(y[b == 389]), "| others", mean(y[b != 389]), "| overall", mean(y), "\n")
lc <- read.csv(file.path(T, "learning", "learning_summary.csv")); cat("resolved before full cohort:", sum(lc$resolved & lc$source == "subsample", na.rm = TRUE), "\n")
cat("usable at n=50 range:", range(lc$U_mean[lc$n == 50 & !grepl("CONTROL", lc$id)]), " estimate at n=50:", range(lc$I_mean[lc$n == 50 & !grepl("CONTROL", lc$id)]), " null at 50:", range(lc$null_mean[lc$n == 50 & !grepl("CONTROL", lc$id)]), "\n")
mx <- tapply(lc$U_mean[lc$source != "x"], lc$id[lc$source != "x"], max); print(round(mx, 3))

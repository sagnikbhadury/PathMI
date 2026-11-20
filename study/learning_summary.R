################################################################
# learning_summary.R
#
# Summarises learning_curves.R into the learning curves and
# N_evidence (definitions in learning_curves.R). Also checks the
# s.e. scaling: if se(n) = sd_B(n) / sqrt(1 - n/N) is the s.e. of a
# study of n independent patients, sqrt(n) * se(n) should be roughly
# flat in n, whereas the m/N -> 0 rate sqrt(n) * sd_B(n) should fall
# as n approaches N.
#
# Inputs: learning/learning_draws.csv, learning/control_full_fit.rds,
#         record/recoverability_record.csv, robustness/subsampling_draws.csv
#         (the stored half-samples of the Record arms), usable/usable_infold.csv
# Output: learning/learning_summary.csv, learning/n_evidence.csv
################################################################
suppressPackageStartupMessages({ library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); OUT <- file.path(TRI, "learning")
Z90 <- qnorm(0.95)
code <- c("UNI v2" = "uni_v2", "Virchow2" = "virchow2", "H-optimus-1" = "hoptimus1", "CONCH v1.5" = "conch_v15")
dr <- read.csv(file.path(OUT, "learning_draws.csv"), stringsAsFactors = FALSE)
rec <- read.csv(file.path(TRI, "record", "recoverability_record.csv"), stringsAsFactors = FALSE)
if ("record_type" %in% names(rec)) rec <- rec[rec$record_type == "measurement", ]
rec <- rec %>% mutate(id = sprintf("%s|%s|%s|%s_%s", cohort, code[model], state, projection, pooling))
hs <- read.csv(file.path(TRI, "robustness", "subsampling_draws.csv"), stringsAsFactors = FALSE)
ui <- read.csv(file.path(TRI, "usable", "usable_infold.csv"), stringsAsFactors = FALSE)
ctl <- readRDS(file.path(OUT, "control_full_fit.rds"))
CTL <- unique(grep("^CONTROL", dr$id, value = TRUE)); stopifnot(length(CTL) == 1); CTL_SRC <- sub("^CONTROL\\|hoptimus1\\|MSI-permuted\\|", "CRC|hoptimus1|MSI|", CTL)

# subsample levels
lev <- dr %>% group_by(id, n, N) %>% summarise(B = n(), I_mean = mean(MI), I_sd = sd(MI), I_q05 = quantile(MI, .05), I_q95 = quantile(MI, .95),
  null_mean = mean(MI_perm), null_sd = sd(MI_perm), U_mean = mean(I_usable), U_sd = sd(I_usable), auc = mean(auc), k_star = median(k_star), frac_zero = mean(MI == 0), .groups = "drop") %>%
  mutate(fpc = sqrt(1 - n / N), se = I_sd / fpc, U_se = U_sd / fpc, source = "subsample")

# full cohort: the Record values; s.e. from the stored half-samples under the jackknife scaling sqrt(m / (N - m))
eps_of <- setNames(rec$eps_rep, rec$id); eps_of[CTL] <- eps_of[[CTL_SRC]]
full <- rec %>% transmute(id, n, N = n, I_mean = I_hat, null_mean, null_sd) %>%
  left_join(hs %>% group_by(id) %>% summarise(m = first(m), sd_half = sd(MI), .groups = "drop"), by = "id") %>%
  mutate(se = sd_half * sqrt(m / (N - m))) %>% select(-m, -sd_half) %>%
  # 19 Nov 2026 (QC): the learning-curve axis counts patients; place the full-cohort point at the patient count
  left_join(dr %>% group_by(id) %>% summarise(Np = first(N), .groups = "drop"), by = "id") %>% mutate(n = Np, N = Np) %>% select(-Np)
ctl_half <- lev %>% filter(id == CTL, n == floor(N / 2))
full <- bind_rows(full, data.frame(id = CTL, n = ctl_half$N, N = ctl_half$N, I_mean = ctl$MI,
  null_mean = rec$null_mean[rec$id == CTL_SRC], null_sd = rec$null_sd[rec$id == CTL_SRC], se = ctl_half$I_sd * sqrt(floor(ctl_half$N / 2) / (ctl_half$N - floor(ctl_half$N / 2))))) %>%
  left_join(ui %>% transmute(id, U_mean = I_usable_infold, U_sd = I_usable_infold_sd, auc = auc_infold), by = "id") %>% mutate(source = "full cohort")
curve <- bind_rows(lev, full) %>% mutate(eps = eps_of[id], LCB = I_mean - Z90 * se, threshold = null_mean + eps, resolved = LCB > threshold) %>% arrange(id, n)
write.csv(curve, file.path(OUT, "learning_summary.csv"), row.names = FALSE)

ne <- curve %>% group_by(id) %>% summarise(N = max(N), N_evidence = if (any(resolved)) min(n[resolved]) else NA_integer_,
  stable = if (any(resolved)) all(resolved[n >= min(n[resolved])]) else NA, status_full = ifelse(resolved[n == max(N)], "resolved", "unresolved"),
  LCB_full = LCB[n == max(N)], threshold_full = threshold[n == max(N)], se_full = se[n == max(N)], .groups = "drop") %>%
  mutate(label = ifelse(is.na(N_evidence), "not resolved at available n", as.character(N_evidence)))
write.csv(ne, file.path(OUT, "n_evidence.csv"), row.names = FALSE)

print(as.data.frame(ne), digits = 3)
cat("\nscaling check: median over arms of sqrt(n) * se(n) [finite-population] and sqrt(n) * sd_B(n) [m/N -> 0 rate], by n (Record arms with I > 0.05)\n")
big <- rec$id[rec$I_hat > 0.05]
print(as.data.frame(lev %>% filter(id %in% big) %>% group_by(n_level = ifelse(n == floor(N / 2), "N/2", as.character(n))) %>%
  summarise(fpc_scaled = median(sqrt(n) * se), prw_scaled = median(sqrt(n) * I_sd), .groups = "drop")), digits = 3)
cat("\nfull-cohort s.e. (jackknife scaling) vs the previous sqrt(m/n) scaling:\n")
print(as.data.frame(full %>% transmute(id, se_jackknife = se, se_previous = se / sqrt(2))), digits = 3)
cat("\nhalf-cohort level from the learning curves (projection refitted, stratified) vs the stored half-samples (fixed projection): s.d. of the estimate\n")
print(as.data.frame(lev %>% filter(n == floor(N / 2)) %>% select(id, I_sd) %>% left_join(hs %>% group_by(id) %>% summarise(sd_stored = sd(MI)), by = "id")), digits = 3)

################################################################
# status_recheck.R  (27 Nov 2026, final-edits check 3)
#
# (i)  Resampling formulas, recomputed from the raw draws:
#      full-cohort s.e. = SD_m * sqrt(m / (N - m))   (delete-d jackknife,
#        half-samples m = floor(N/2): robustness/subsampling_draws.csv);
#      learning-curve s.e. for a study of size n = SD_n / sqrt(1 - n/N)
#        (finite-population correction: learning/learning_draws.csv).
#      At n = m = N/2 the two factors are 1 and sqrt(2): they are the
#      s.e. of studies of N and of N/2 patients, and do not coincide.
# (ii) Eq. (8) with the revised tolerance (robustness/tolerance.csv) at
#      every point of every learning curve: LCB90 = I - 1.645 s.e. >
#      I_null(n) + eps_tol, with I_null the fixed-arm null; resolution
#      margins, N_evidence and the counts, compared with the stored
#      learning/learning_summary.csv; plus a sensitivity with per-model
#      tolerances (max over the two identities for each model).
#
# Run after tolerance.R and learning_summary.R.
# Output: record/status_recheck.csv, record/status_recheck.txt
################################################################
suppressPackageStartupMessages({ library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/")
TXT <- file.path(TRI, "record", "status_recheck.txt"); cat("", file = TXT)
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = TXT, append = TRUE) }
Z90 <- qnorm(0.95)
code <- c("UNI v2" = "uni_v2", "Virchow2" = "virchow2", "H-optimus-1" = "hoptimus1", "CONCH v1.5" = "conch_v15")

# ---- (i) resampling formulas -------------------------------------------------------
hd <- read.csv(file.path(TRI, "robustness", "subsampling_draws.csv"), stringsAsFactors = FALSE)
hs <- read.csv(file.path(TRI, "robustness", "subsampling_summary.csv"), stringsAsFactors = FALSE)
chk1 <- hd %>% group_by(id) %>% summarise(m = first(m), sd_m = sd(MI), .groups = "drop") %>% inner_join(hs %>% select(id, n, se_n), by = "id") %>%
  mutate(se_recomputed = sd_m * sqrt(m / (n - m)), factor = sqrt(m / (n - m)))
say("(i) full-cohort s.e. = SD_m sqrt(m/(N-m)): %d arms | max |recomputed - stored| = %.2g | factor range %.3f-%.3f (m = floor(N/2))",
    nrow(chk1), max(abs(chk1$se_recomputed - chk1$se_n)), min(chk1$factor), max(chk1$factor))
ld <- read.csv(file.path(TRI, "learning", "learning_draws.csv"), stringsAsFactors = FALSE)
ls <- read.csv(file.path(TRI, "learning", "learning_summary.csv"), stringsAsFactors = FALSE)
chk2 <- ld %>% group_by(id, n, N) %>% summarise(sd_n = sd(MI), .groups = "drop") %>% mutate(se_recomputed = sd_n / sqrt(1 - n / N)) %>%
  inner_join(ls %>% filter(source == "subsample") %>% select(id, n, se), by = c("id", "n"))
say("    learning s.e. = SD_n / sqrt(1 - n/N): %d points | max |recomputed - stored| = %.2g", nrow(chk2), max(abs(chk2$se_recomputed - chk2$se)))
grid <- ls %>% group_by(id) %>% summarise(N = max(N), grid = paste(sort(unique(n)), collapse = ", "), .groups = "drop")
say("    evaluated grid per arm: %s", paste(unique(sprintf("N = %d: {%s}", grid$N, grid$grid)), collapse = " | "))
say("    at n = N/2 the learning factor is %.3f and the jackknife factor 1.000: different quantities (studies of N/2 and of N patients)", 1 / sqrt(0.5))

# ---- (ii) Eq. (8) with the revised tolerance -----------------------------------------
TOL <- read.csv(file.path(TRI, "robustness", "tolerance.csv"), stringsAsFactors = FALSE)
eps_grp <- c(uni_v2 = "vision-only", virchow2 = "vision-only", hoptimus1 = "vision-only", conch_v15 = "CONCH v1.5")
# per-model tolerance (sensitivity only): max over the two identities for each model, from tolerance_identities.csv
ti <- read.csv(file.path(TRI, "robustness", "tolerance_identities.csv"), stringsAsFactors = FALSE)
pm <- ti %>% group_by(model) %>% summarise(m = max(max_abs), .groups = "drop")
eps_model <- setNames(round(pm$m[match(names(code), pm$model)], 4), code)
mcode <- function(id) { x <- sub("^CONTROL\\|", "", id); x <- sub("^[A-Z]+\\|", "", x); sub("\\|.*", "", x) }
rc <- ls %>% mutate(mc = mcode(id), eps_new = TOL$eps_tol_rounded[match(eps_grp[mc], TOL$group)], eps_pm = eps_model[mc],
  LCB_new = I_mean - Z90 * se, thr_new = null_mean + eps_new, resolved_new = LCB_new > thr_new, margin = LCB_new - thr_new,
  resolved_pm = LCB_new > null_mean + eps_pm)
say("(ii) Eq. (8) with eps_tol = %s: %d curve points | status changes vs stored: %d | stored eps differs at %d points",
    paste(sprintf("%s %.4f", TOL$group, TOL$eps_tol_rounded), collapse = ", "), nrow(rc), sum(rc$resolved_new != rc$resolved), sum(abs(rc$eps_new - rc$eps) > 1e-12))
ne <- rc %>% group_by(id) %>% arrange(n) %>% summarise(N = max(N), status = ifelse(resolved_new[n == max(N)], "resolved", "unresolved"),
  margin_full = margin[n == max(N)], LCB_full = LCB_new[n == max(N)], thr_full = thr_new[n == max(N)],
  N_evidence = if (any(resolved_new)) min(n[resolved_new]) else NA_integer_, status_pm = ifelse(resolved_pm[n == max(N)], "resolved", "unresolved"),
  best_sub_margin = max(margin[n < max(N)]), n_best_sub = n[n < max(N)][which.max(margin[n < max(N)])], .groups = "drop")
write.csv(ne, file.path(TRI, "record", "status_recheck.csv"), row.names = FALSE)
meas <- ne %>% filter(!grepl("^CONTROL", id))
say("    records resolved: %d of %d measurements | control: %s | per-model tolerance sensitivity: %d of %d resolved (changes: %d)",
    sum(meas$status == "resolved"), nrow(meas), ne$status[grepl("^CONTROL", ne$id)], sum(meas$status_pm == "resolved"), nrow(meas), sum(meas$status != meas$status_pm))
say("    resolved before the full cohort: %d | N_evidence values: %s", sum(!is.na(meas$N_evidence) & meas$N_evidence < meas$N), paste(sort(unique(na.omit(meas$N_evidence))), collapse = ", "))
for (i in seq_len(nrow(ne))) say("    %-40s %-10s margin %+.4f (LCB %.4f, thr %.4f) | N_evidence %s | best sub-cohort margin %+.4f at n = %d | per-model eps: %s",
  ne$id[i], ne$status[i], ne$margin_full[i], ne$LCB_full[i], ne$thr_full[i], ifelse(is.na(ne$N_evidence[i]), "none", ne$N_evidence[i]), ne$best_sub_margin[i], ne$n_best_sub[i], ne$status_pm[i])

################################################################
# recoverability_record.R
#
# The Recoverability Record: one row per (foundation model, state),
# the measurement reported as a metrology object rather than a
# leaderboard entry, plus one negative-control record:
#   provenance  model (Trident feature key), tile geometry at 20x,
#               embedding width, projection and pooling of the
#               recorded configuration, covariate set, protocol;
#   measurement n, H(Y), I_hat (represented: the model's best
#               configuration), I/H(Y), I_null (permutation mean and
#               s.d.), epsilon_rep, half-sampling 90% interval
#               (delete-d jackknife s.e.), I_adjusted (attributable:
#               incremental given covariates, same arm), I_usable
#               (usable: cross-fitted ridge readout, every
#               transformation refitted in the training patients) with
#               out-of-patient AUC, N90 and N95 (tiles, under the
#               model's linear readout);
#   evidence    status (resolved if the 90% lower bound exceeds
#               I_null + epsilon_rep) and N_evidence (patients; see
#               learning_curves.R).
# The negative control is the H-optimus-1 CRC/MSI arm with its labels
# permuted once (usable_infold.R, learning_curves.R).
#
# Inputs: merged tables, robustness/subsampling_summary.csv,
#         usable/usable_infold.csv, learning/{learning_summary,n_evidence}.csv,
#         learning/control_full_fit.rds, burden/burden_Nq_bootstrap.csv,
#         prad_incremental_batch_purity/prad_incremental_batch_purity.csv
# Output: TRIDENT/record/recoverability_record.csv, record/json/*.json,
#         record/recoverability_record.json and manuscript/table_record.tex
################################################################
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(jsonlite) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); MS <- file.path(dirname(TRI), "manuscript")
OUT <- file.path(TRI, "record"); dir.create(file.path(OUT, "json"), recursive = TRUE, showWarnings = FALSE)
BB <- c("UNI v2", "Virchow2", "H-optimus-1", "CONCH v1.5")
harm <- function(e) case_when(grepl("^UNI", e) ~ BB[1], grepl("^VIRCHOW", e) ~ BB[2], grepl("^H-OPT", e) ~ BB[3], grepl("^CONCH", e) ~ BB[4])
code <- c("UNI v2" = "uni_v2", "Virchow2" = "virchow2", "H-optimus-1" = "hoptimus1", "CONCH v1.5" = "conch_v15")
tile_px <- c("UNI v2" = 256, "Virchow2" = 224, "H-optimus-1" = 224, "CONCH v1.5" = 512)
width <- c("UNI v2" = 1536, "Virchow2" = 2560, "H-optimus-1" = 1536, "CONCH v1.5" = 768)
EPS <- c("UNI v2" = 0.0163, "Virchow2" = 0.0163, "H-optimus-1" = 0.0163, "CONCH v1.5" = 0.0215)   # 19 Nov 2026 QC: vision-only bound 0.0095 (both cohorts)
CTL <- NULL; CTL_SRC <- NULL   # set after the usable table is read (the control is the H-optimus-1 CRC/MSI selected arm, permuted)

g <- read_excel(file.path(TRI, "mi summary table.xlsx")) %>% mutate(bb = harm(Embeddings), task = paste0(Disease, "/", target), decomp = tolower(decomp))
inc <- read_excel(file.path(TRI, "incremental mi summary.xlsx")) %>% mutate(bb = harm(Embeddings))
pc <- read.csv(file.path(TRI, "prad_incremental_batch_purity", "prad_incremental_batch_purity.csv"), stringsAsFactors = FALSE)
bs <- read.csv(file.path(TRI, "robustness", "subsampling_summary.csv"), stringsAsFactors = FALSE)
ui <- read.csv(file.path(TRI, "usable", "usable_infold.csv"), stringsAsFactors = FALSE)
CTL <- grep("^CONTROL", ui$id, value = TRUE); stopifnot(length(CTL) == 1); CTL_SRC <- sub("^CONTROL\\|hoptimus1\\|MSI-permuted\\|", "CRC|hoptimus1|MSI|", CTL)
u0 <- read.csv(file.path(TRI, "usable", "usable_recoverability.csv"), stringsAsFactors = FALSE)
# first pass of the patient-level chain runs before the learning curves exist; the second pass fills status and N_evidence
HAVE_LC <- file.exists(file.path(TRI, "learning", "n_evidence.csv")) && file.exists(file.path(TRI, "learning", "control_full_fit.rds")) &&
  any(grepl("^CONTROL", read.csv(file.path(TRI, "learning", "learning_summary.csv"))$id))
ne <- if (HAVE_LC) read.csv(file.path(TRI, "learning", "n_evidence.csv"), stringsAsFactors = FALSE) else data.frame(id = character(), N_evidence = integer(), status_full = character(), LCB_full = numeric(), threshold_full = numeric())
lc <- if (HAVE_LC) read.csv(file.path(TRI, "learning", "learning_summary.csv"), stringsAsFactors = FALSE) else data.frame(id = CTL, source = "full cohort", I_mean = NA_real_, se = NA_real_, n = NA, U_mean = NA_real_, LCB = NA_real_, threshold = NA_real_, resolved = NA, null_mean = NA_real_)
nq <- read.csv(file.path(TRI, "burden", "burden_Nq_bootstrap.csv"), stringsAsFactors = FALSE) %>%
  transmute(cohort = ifelse(state == "MSI", "CRC", "PRAD"), target = state, model, N90, N95, N90_min = N90_lo, N90_max = N90_hi, N95_min = N95_lo, N95_max = N95_hi)   # 90% intervals over resamples of the tile draws

rec <- g %>% group_by(Disease, target, bb) %>% slice_max(MI_obs, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(cohort = Disease, state = target, model = bb, key, projection = decomp, pooling = agg_method, I_hat = MI_obs, null_mean, null_sd, p_value)
rec$id <- with(rec, sprintf("%s|%s|%s|%s", cohort, code[model], state, key))
rec <- rec %>% left_join(u0 %>% transmute(id, n, prevalence, H_Y, I_insample_plugin), by = "id") %>%
  left_join(ui %>% transmute(id, I_usable = I_usable_infold, I_usable_sd = I_usable_infold_sd, auc_oof = auc_infold, I_usable_fixed_projection = I_usable_fixed,
                             I_usable_pathmi = I_usable_pathmi_infold, auc_oof_pathmi = auc_pathmi_infold), by = "id") %>%
  left_join(bs %>% transmute(id, se_n, ci90_lo, ci90_hi), by = "id") %>%
  left_join(nq %>% transmute(cohort, state = target, model, N90, N95, N90_min, N90_max, N95_min, N95_max), by = c("cohort", "state", "model"))
rec <- rec %>% left_join(inc %>% filter(Disease == "CRC") %>% transmute(cohort = Disease, state = target, model = bb, pooling = agg_method, I_adj_crc = incremental_subho_mi, cov_crc = covar_only), by = c("cohort", "state", "model", "pooling")) %>%
  left_join(pc %>% transmute(cohort = "PRAD", state = target, model, pooling, I_adj_prad = incremental_all, cov_prad = covar_all), by = c("cohort", "state", "model", "pooling")) %>%
  mutate(I_adjusted = coalesce(I_adj_crc, I_adj_prad), I_covariates = coalesce(cov_crc, cov_prad),
         covariates = ifelse(cohort == "CRC", "size, stage", "size, Gleason, batch, purity"), eps_rep = EPS[model], record_type = "measurement")

# negative control: the H-optimus-1 CRC/MSI arm with permuted labels
src <- rec %>% filter(id == CTL_SRC); cfull <- lc %>% filter(id == CTL, source == "full cohort"); cu <- ui %>% filter(id == CTL)
ctl <- src %>% mutate(id = CTL, state = "MSI (permuted labels)", I_hat = cfull$I_mean, p_value = NA_real_, se_n = cfull$se, ci90_lo = max(cfull$I_mean - 1.645 * cfull$se, 0), ci90_hi = cfull$I_mean + 1.645 * cfull$se,
  I_usable = cu$I_usable_infold, I_usable_sd = cu$I_usable_infold_sd, auc_oof = cu$auc_infold, I_usable_fixed_projection = NA_real_, I_usable_pathmi = cu$I_usable_pathmi_infold, auc_oof_pathmi = cu$auc_pathmi_infold,
  I_insample_plugin = NA_real_, I_adjusted = NA_real_, I_covariates = NA_real_, covariates = "none", N90 = NA, N95 = NA, N90_min = NA, N90_max = NA, N95_min = NA, N95_max = NA, record_type = "negative control")
rec <- bind_rows(rec, ctl) %>%
  left_join(ne %>% transmute(id, N_evidence, status = status_full, LCB90 = LCB_full, threshold = threshold_full), by = "id") %>%
  mutate(share = I_hat / H_Y, tile_px = tile_px[model], width = width[model], trident_key = code[model], magnification = "20x, ~0.5 um/px", usable_share = I_usable / I_hat,
         N_evidence_label = ifelse(is.na(N_evidence), "not resolved at available n", as.character(N_evidence))) %>%
  select(record_type, cohort, state, model, trident_key, width, tile_px, magnification, projection, pooling, covariates, n, prevalence, H_Y, I_hat, share, null_mean, null_sd, p_value,
         eps_rep, se_n, ci90_lo, ci90_hi, LCB90, threshold, status, N_evidence, N_evidence_label, I_covariates, I_adjusted, I_usable, I_usable_sd, usable_share, auc_oof, I_usable_fixed_projection,
         I_usable_pathmi, auc_oof_pathmi, I_insample_plugin, N90, N90_min, N90_max, N95, N95_min, N95_max, id) %>%
  arrange(record_type != "measurement", factor(paste(cohort, state), c("CRC MSI", "PRAD ERG", "PRAD ETS")), factor(model, BB))
write.csv(rec, file.path(OUT, "recoverability_record.csv"), row.names = FALSE)

# ---- JSON records -------------------------------------------------------------------
num <- function(x, d = 4) if (is.na(x)) NULL else round(x, d)
to_json <- function(r) {
  cu <- lc %>% filter(id == r$id) %>% arrange(n)
  list(schema = "pathmi.recoverability_record/1.0", record_type = r$record_type,
    provenance = list(model = r$model, trident_feature_key = r$trident_key, embedding_width = r$width, tile_px_at_20x = r$tile_px, magnification = r$magnification,
                      cohort = paste0("TCGA-", r$cohort), state = r$state, projection = r$projection, retained_dimension = 50, pooling = r$pooling,
                      covariates = r$covariates, roster_n = r$n, prevalence = num(r$prevalence)),
    protocol = list(estimator = "PathMI (horseshoe ranking, Laplace marginal-likelihood walk, log Bayes factor per patient), v_beta = 1, 2000 Gibbs iterations, 500 burn-in, thin 10",
                    null = "200 label permutations, whole estimator refitted", uncertainty = if (r$record_type == "measurement") "half-sampling, B = 200, delete-d jackknife s.e." else "half-sampling, B = 50 label-stratified half-cohort subsamples with the projection refitted, delete-d jackknife s.e.",
                    usable = "ridge logistic, stratified 5-fold x 10 repeats, standardisation, projection and penalty fitted in training patients only",
                    evidence_rule = "resolved if Ihat - 1.645 s.e. > I_null + epsilon_rep", units = "nats"),
    measurement = list(label_entropy = num(r$H_Y), represented = num(r$I_hat), share_of_entropy = num(r$share, 3),
                       null_mean = num(r$null_mean), null_sd = num(r$null_sd), permutation_p = num(r$p_value), epsilon_rep = num(r$eps_rep),
                       se = num(r$se_n), ci90 = c(num(r$ci90_lo), num(r$ci90_hi)), attributable = num(r$I_adjusted), covariate_only = num(r$I_covariates),
                       usable = num(r$I_usable), usable_sd = num(r$I_usable_sd), usable_auc = num(r$auc_oof, 3),
                       tile_burden = if (is.na(r$N90)) NULL else list(N90 = r$N90, N90_interval = c(r$N90_min, r$N90_max), N95 = r$N95, N95_interval = c(r$N95_min, r$N95_max))),
    evidence = list(status = r$status, lower_bound_90 = num(r$LCB90), threshold = num(r$threshold), N_evidence = if (is.na(r$N_evidence)) "not resolved at available n" else r$N_evidence,
                    learning_curve = lapply(seq_len(nrow(cu)), function(j) list(n = cu$n[j], represented_mean = num(cu$I_mean[j]), se = num(cu$se[j]), null_mean = num(cu$null_mean[j]),
                                                                                    usable_mean = num(cu$U_mean[j]), resolved = cu$resolved[j]))),
    interpretation = if (isTRUE(r$status == "resolved")) "The representation carries information about the state above the estimator's null and resolution at this patient budget."
                     else "An unresolved record indicates insufficient evidence under the evaluated representation, patient budget and protocol; it is not evidence that molecular information is universally absent.")
}
recs <- lapply(seq_len(nrow(rec)), function(i) to_json(rec[i, ]))
for (i in seq_len(nrow(rec))) write_json(recs[[i]], file.path(OUT, "json", paste0(gsub("[^A-Za-z0-9]+", "_", rec$id[i]), ".json")), auto_unbox = TRUE, pretty = TRUE, digits = NA)
write_json(recs, file.path(OUT, "recoverability_record.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

# ---- LaTeX body -------------------------------------------------------------------
f3 <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x)); f2 <- function(x) ifelse(is.na(x), "--", sprintf("%.2f", x)); fs <- function(x) ifelse(x < 0, sprintf("$-$%.3f", -x), sprintf("%.3f", x))
proj_lab <- c(pca = "linear", cosine = "cosine", rbf = "RBF", poly = "poly")
st_lab <- c(MSI = "CRC/MSI", ERG = "PRAD/ERG", ETS = "PRAD/ETS")
row_tex <- function(r, first, lab) sprintf("%s & %s & %d & %s & %s & %.0f & %s$\\pm$%s & %s & %s to %s & %s & %s & %s & %s & %s & %s & %s & %s/%s \\\\",
  if (first) lab else "", r$model, r$n, f3(r$H_Y), f3(r$I_hat), 100 * r$share, f3(r$null_mean), f3(r$null_sd), sprintf("%.4f", r$eps_rep), f3(r$ci90_lo), f3(r$ci90_hi),
  f3(r$I_adjusted), fs(r$I_usable), f2(r$auc_oof), ifelse(is.na(r$N90), "--", r$N90), ifelse(is.na(r$N95), "--", r$N95),
  r$status, ifelse(is.na(r$N_evidence), "$>n$", as.character(r$N_evidence)), proj_lab[[r$projection]], r$pooling)
lines <- c()
for (st in c("MSI", "ERG", "ETS")) { d <- rec %>% filter(record_type == "measurement", state == st)
  for (i in seq_len(nrow(d))) lines <- c(lines, row_tex(d[i, ], i == 1, st_lab[[st]])); lines <- c(lines, "\\addlinespace[2pt]") }
lines <- c(lines, "\\midrule", row_tex(rec %>% filter(record_type == "negative control"), TRUE, "Control"))
tex <- c("\\begin{tabular}{llccccccccccccccl}", "\\toprule",
  "State & Model & $n$ & $H(Y)$ & $\\widehat I$ & \\%$H$ & $I_{\\mathrm{null}}$ & $\\epsilon_{\\mathrm{rep}}$ & 90\\% CI & $I_{\\mathrm{adj}}$ & $I_{\\mathrm{usable}}$ & AUC & $N_{90}$ & $N_{95}$ & Status & $N_{\\mathrm{evidence}}$ & Projection/pooling \\\\",
  "\\midrule", lines, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(MS, "table_record.tex"))
cat("=== Recoverability Record ===\n"); print(as.data.frame(rec %>% select(state, model, n, I_hat, null_mean, eps_rep, se_n, ci90_lo, ci90_hi, status, N_evidence, I_adjusted, I_usable, auc_oof, N90)), digits = 3)
m <- rec %>% filter(record_type == "measurement")
cat(sprintf("\nusable / present: median %.2f (range %.2f to %.2f) | adjusted / present: median %.2f (range %.2f to %.2f)\n", median(m$usable_share), min(m$usable_share), max(m$usable_share),
            median(m$I_adjusted / m$I_hat), min(m$I_adjusted / m$I_hat), max(m$I_adjusted / m$I_hat)))

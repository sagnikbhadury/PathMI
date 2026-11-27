################################################################
# tolerance.R  (27 Nov 2026, final-edits check 2)
#
# The empirical numerical tolerance epsilon_tol of the protocol, from
# saved outputs only. Two numerical identities enter it:
#   replicate refit   the embedding-only fit of the incremental block
#                     repeats the selected fit with an independent
#                     random-number stream (48 selected arms);
#   sign-flip refit   the full-data refit of every cell after the
#                     projection is recomputed by an independent routine
#                     (column signs may flip; 192 cells).
# The null-covariate block (combined minus embedding-only fit where the
# covariate-only estimate is zero, 32 prostate arms) adds a covariate
# block to the model. It is a sensitivity to model specification, not
# a numerical identity, and is reported separately (not in the tolerance).
# The chain-seed range (five horseshoe seeds, eight configurations) is
# also reported separately; the CONCH v1.5 poly/moment MSI arm is
# flagged as unstable, and the tolerance does not cover it.
#
# Groups: vision-only models (UNI v2, Virchow2, H-optimus-1) and
# CONCH v1.5. epsilon_tol of a group = the largest absolute discrepancy
# over the two identities among that group's arms.
#
# Inputs: best kpca summary.xlsx, incremental mi summary.xlsx,
#         usable/usable_infold.csv, robustness/chain_sizeprior.csv
# Output: robustness/tolerance.csv (group values read by the downstream
#         scripts), robustness/tolerance_identities.csv (identity x model),
#         robustness/tolerance_specification.csv, robustness/tolerance.txt
################################################################
suppressPackageStartupMessages({ library(readxl); library(dplyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); OUT <- file.path(TRI, "robustness")
TXT <- file.path(OUT, "tolerance.txt"); cat("", file = TXT)
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = TXT, append = TRUE) }
emb2code <- c("UNI V2" = "uni_v2", "VIRCHOW" = "virchow2", "H-OPTIMUS 1" = "hoptimus1", "CONCH V1.5" = "conch_v15")
lab <- c(uni_v2 = "UNI v2", virchow2 = "Virchow2", hoptimus1 = "H-optimus-1", conch_v15 = "CONCH v1.5")
grp <- function(code) ifelse(code == "conch_v15", "CONCH v1.5", "vision-only")

best <- as.data.frame(read_excel(file.path(TRI, "best kpca summary.xlsx"))); best$code <- emb2code[best$Embeddings]
inc <- as.data.frame(read_excel(file.path(TRI, "incremental mi summary.xlsx"))); inc$code <- emb2code[inc$Embeddings]
u <- read.csv(file.path(TRI, "usable", "usable_infold.csv"), stringsAsFactors = FALSE); u <- u[!u$control, ]
stopifnot(!any(is.na(best$code)), !any(is.na(inc$code)))

j <- merge(inc, best[, c("Disease", "target", "code", "agg_method", "MI_obs")], by = c("Disease", "target", "code", "agg_method"))
id_rep <- data.frame(identity = "replicate refit", code = j$code, cell = paste(j$Disease, j$code, j$target, j$agg_method, sep = "|"), d = j$emb_only - j$MI_obs)
id_flip <- data.frame(identity = "sign-flip refit", code = sub("\\|.*", "", sub("^[A-Z]+\\|", "", u$id)), cell = u$id, d = u$MI_refit - u$MI_obs)
stopifnot(nrow(id_rep) == 48, nrow(id_flip) == 192)
ids <- rbind(id_rep, id_flip); ids$group <- grp(ids$code)

per_id <- ids %>% group_by(identity, model = lab[code], group) %>%
  summarise(n = n(), exact = sum(abs(d) < 1e-9), q95_abs = unname(quantile(abs(d), .95)), max_abs = max(abs(d)), cell_at_max = cell[which.max(abs(d))], .groups = "drop") %>%
  arrange(identity, factor(model, lab))
write.csv(per_id, file.path(OUT, "tolerance_identities.csv"), row.names = FALSE)

tol <- ids %>% group_by(group) %>% summarise(eps_tol = max(abs(d)), identity_at_max = identity[which.max(abs(d))], cell_at_max = cell[which.max(abs(d))],
  n_checks = n(), exact = sum(abs(d) < 1e-9), .groups = "drop")
per_model <- ids %>% group_by(model = lab[code]) %>% summarise(max_abs = max(abs(d)), .groups = "drop")
tol$models <- ifelse(tol$group == "vision-only", "UNI v2, Virchow2, H-optimus-1", "CONCH v1.5")
tol$eps_tol_rounded <- round(tol$eps_tol, 4)
write.csv(tol, file.path(OUT, "tolerance.csv"), row.names = FALSE)

# specification sensitivity and chain-seed instability (reported separately)
spec <- inc %>% filter(covar_only == 0) %>% transmute(check = "null covariate block (specification)", model = lab[code], cell = paste(Disease, code, target, agg_method, sep = "|"), d = combined - emb_only)
cs <- read.csv(file.path(OUT, "chain_sizeprior.csv"), stringsAsFactors = FALSE) %>% filter(what == "seed") %>% group_by(id) %>%
  summarise(lo = min(MI), hi = max(MI), range = hi - lo, .groups = "drop")
write.csv(spec, file.path(OUT, "tolerance_specification.csv"), row.names = FALSE)

say("=== numerical tolerance (replicate refit + sign-flip refit; null covariate block excluded) ===")
for (i in seq_len(nrow(tol))) say("%-11s eps_tol = %.4f nats  (max at %s, %s) | %d checks, %d exact", tol$group[i], tol$eps_tol[i], tol$identity_at_max[i], tol$cell_at_max[i], tol$n_checks[i], tol$exact[i])
say("per model (max over the two identities): %s", paste(sprintf("%s %.4f", per_model$model, per_model$max_abs), collapse = "; "))
for (i in seq_len(nrow(per_id))) say("  %-16s %-12s n %3d exact %3d q95 %.4f max %.4f (%s)", per_id$identity[i], per_id$model[i], per_id$n[i], per_id$exact[i], per_id$q95_abs[i], per_id$max_abs[i], per_id$cell_at_max[i])
say("specification sensitivity (null covariate block, 32 PRAD arms; not in the tolerance): %d exact, max |d| %.4f (%s)",
    sum(abs(spec$d) < 1e-9), max(abs(spec$d)), spec$cell[which.max(abs(spec$d))])
say("chain seeds (five seeds, eight configurations): %s", paste(sprintf("%s %.4f-%.4f", cs$id, cs$lo, cs$hi), collapse = "; "))
uns <- cs %>% filter(range > max(tol$eps_tol))
say("unstable beyond the tolerance (flagged, not covered): %s", if (nrow(uns)) paste(sprintf("%s range %.3f (%.3f to %.3f)", uns$id, uns$range, uns$lo, uns$hi), collapse = "; ") else "none")

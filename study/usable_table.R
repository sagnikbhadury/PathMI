################################################################
# usable_table.R — Supplementary Table: usable recoverability by state
# and model under the two readout families (means over the 16 cells of
# each model and state; out-of-patient AUC), with the PathMI estimate
# and the in-sample plug-in for comparison. Output: manuscript/table_usable.tex
################################################################
suppressPackageStartupMessages(library(dplyr))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); MS <- file.path(dirname(TRI), "manuscript")
BB <- c("UNI v2", "Virchow2", "H-optimus-1", "CONCH v1.5")
u <- read.csv(file.path(TRI, "usable", "usable_recoverability.csv"), stringsAsFactors = FALSE) %>% mutate(task = paste0(cohort, "/", target), model = factor(model, BB)) %>%
  left_join(read.csv(file.path(TRI, "usable", "usable_infold.csv"), stringsAsFactors = FALSE) %>% select(id, I_usable_infold, I_usable_infold_sd, auc_infold), by = "id") %>%
  mutate(I_usable_ridge = I_usable_infold, I_usable_ridge_sd = I_usable_infold_sd, auc_oof_ridge = auc_infold)   # ridge readout with the projection refitted in each training fold
s <- u %>% group_by(task, model) %>% summarise(present = mean(MI_obs), insample = mean(I_insample_plugin), ridge = mean(I_usable_ridge), ridge_sd = mean(I_usable_ridge_sd), auc_r = mean(auc_oof_ridge),
  pathmi = mean(I_usable), pathmi_sd = mean(I_usable_sd), auc_p = mean(auc_oof), neg = sum(I_usable < 0), .groups = "drop") %>% arrange(factor(task, c("CRC/MSI", "PRAD/ERG", "PRAD/ETS")), model)
lines <- c(); last <- ""
for (i in seq_len(nrow(s))) with(s[i, ], { lines <<- c(lines, sprintf("%s & %s & %.3f & %.3f & %.3f & %.2f & %s%.3f & %.2f & %d \\\\",
  if (task != last) task else "", model, present, insample, ridge, auc_r, if (pathmi < 0) "$-$" else "", abs(pathmi), auc_p, neg)); last <<- task })
tex <- c("\\begin{tabular}{llccccccc}", "\\toprule",
  " & & & & \\multicolumn{2}{c}{Ridge readout} & \\multicolumn{3}{c}{PathMI readout} \\\\", "\\cmidrule(lr){5-6}\\cmidrule(lr){7-9}",
  "State & Model & $\\widehat I$ & In-sample & $I_{\\mathrm{usable}}$ & AUC & $I_{\\mathrm{usable}}$ & AUC & Negative (of 16) \\\\",
  "\\midrule", lines, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(MS, "table_usable.tex")); cat(tex, sep = "\n")

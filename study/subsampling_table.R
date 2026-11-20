################################################################
# subsampling_table.R — Supplementary Table: patient-level sampling
# uncertainty of the twelve recorded arms by half-sampling (recorded)
# and by the bootstrap with replacement (shown for the reason it is
# not used: duplicated patients inflate the walk). Output:
# manuscript/table_subsampling.tex
################################################################
suppressPackageStartupMessages(library(dplyr))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); MS <- file.path(dirname(TRI), "manuscript")
BB <- c("UNI v2", "Virchow2", "H-optimus-1", "CONCH v1.5")
rec <- read.csv(file.path(TRI, "record", "recoverability_record.csv"), stringsAsFactors = FALSE); rec <- rec[rec$record_type == "measurement", ]
code <- c("UNI v2" = "uni_v2", "Virchow2" = "virchow2", "H-optimus-1" = "hoptimus1", "CONCH v1.5" = "conch_v15")
rec$id <- with(rec, sprintf("%s|%s|%s|%s", cohort, code[model], state, paste0(projection, "_", pooling)))
ss <- read.csv(file.path(TRI, "robustness", "subsampling_summary.csv"), stringsAsFactors = FALSE)
bs <- read.csv(file.path(TRI, "robustness", "bootstrap_summary.csv"), stringsAsFactors = FALSE)
d <- rec %>% select(cohort, state, model, id, I_hat) %>% left_join(ss %>% select(id, k_star, sub_mean, se_n, ci90_lo, ci90_hi, k_median_sub = k_median), by = "id") %>%
  left_join(bs %>% select(id, boot_mean_wr = boot_mean, boot_se_wr = boot_se, k_hi_wr = k_hi), by = "id") %>% arrange(factor(state, c("MSI", "ERG", "ETS")), factor(model, BB))
st_lab <- c(MSI = "CRC/MSI", ERG = "PRAD/ERG", ETS = "PRAD/ETS"); lines <- c(); last <- ""
for (i in seq_len(nrow(d))) with(d[i, ], { lines <<- c(lines, sprintf("%s & %s & %.3f & %d & %.3f & %.4f & %.3f to %.3f & %g & %.3f & %.4f & %g \\\\",
  if (state != last) st_lab[[state]] else "", model, I_hat, k_star, sub_mean, se_n, ci90_lo, ci90_hi, k_median_sub, boot_mean_wr, boot_se_wr, k_hi_wr)); last <<- state })
tex <- c("\\begin{tabular}{llcccccccccc}", "\\toprule", " & & & & \\multicolumn{4}{c}{Half-sampling (recorded)} & \\multicolumn{3}{c}{Bootstrap with replacement} \\\\", "\\cmidrule(lr){5-8}\\cmidrule(lr){9-11}",
  "State & Model & $\\widehat I$ & $k^*$ & Mean & s.e. & 90\\% interval & $k^*$ median & Mean & s.e. & $k^*$ 95th pct. \\\\", "\\midrule", lines, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(MS, "table_subsampling.tex")); cat(tex, sep = "\n")

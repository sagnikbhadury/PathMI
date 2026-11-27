################################################################
# resolution_table.R — Supplementary Table: every major contrast in
# measurement units, from robustness/resolution_units.csv.
# Output: manuscript/table_resolution.tex (\input into the paper)
################################################################
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); MS <- file.path(dirname(TRI), "manuscript")
r <- read.csv(file.path(TRI, "robustness", "resolution_units.csv"), stringsAsFactors = FALSE)
sg <- function(x, d) ifelse(is.na(x), "", sprintf(paste0("%+.", d, "f"), x))
cmp <- gsub("pca", "linear", r$comparison); cmp <- gsub(" kernel PCA", "", cmp); cmp <- gsub("PathMI minus cross-fitted usable \\(ridge readout, in-fold projection\\)", "PathMI minus usable", cmp)
grp <- ifelse(c(TRUE, r$group[-1] != r$group[-nrow(r)]), r$group, "")
lines <- sprintf("%s & %s & %s & %s & %s & %s \\\\", grp, cmp, sg(r$delta, 4), sg(r$units_rep, 1), ifelse(is.na(r$se_half), "", sprintf("%.3f", r$se_half)), sg(r$units_se, 1))
tex <- c("\\begin{tabular}{llcccc}", "\\toprule", "Group & Comparison & $\\Delta I$ & $\\Delta I/\\epsilon_{\\mathrm{tol}}$ & s.e. & $\\Delta I/\\mathrm{s.e.}$ \\\\", "\\midrule", lines, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(MS, "table_resolution.tex")); cat(tex, sep = "\n")

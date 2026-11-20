################################################################
# robustness_tables.R — Supplementary Tables for the estimator
# robustness sweeps on the eight configurations: retained dimension
# (estimate, k*, null mean, kappa at each d) and working prior variance
# (estimate, k*, null mean at each v_beta), plus chain seeds and the
# size prior. Output: manuscript/table_dimension.tex, table_vbeta.tex
################################################################
suppressPackageStartupMessages(library(dplyr))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); MS <- file.path(dirname(TRI), "manuscript")
dim <- read.csv(file.path(TRI, "robustness", "dimension_sweep.csv"), stringsAsFactors = FALSE)
vb  <- read.csv(file.path(TRI, "robustness", "vbeta_sweep.csv"), stringsAsFactors = FALSE)
ch  <- read.csv(file.path(TRI, "robustness", "chain_sizeprior.csv"), stringsAsFactors = FALSE)
lab <- function(d) sprintf("%s/%s, %s, %s", d$cohort, d$target, d$model, gsub("_", " ", d$key))
# primaries (best of each state, MSI, ERG, ETS) then the five difficult configurations, as defined in the sweep scripts
pr <- unique(dim[dim$role == "primary", c("id", "target")]); pr <- pr$id[order(match(pr$target, c("MSI", "ERG", "ETS")))]
ord <- c(pr, c("CRC|conch_v15|MSI|poly_moment", "PRAD|virchow2|ERG|cosine_moment", "PRAD|conch_v15|ERG|poly_moment", "CRC|hoptimus1|MSI|pca_mean", "PRAD|hoptimus1|ERG|pca_moment"))
stopifnot(all(ord %in% dim$id), all(ord %in% vb$id), length(unique(ord)) == 8)
# dimension
lines <- c()
for (id in ord) { d <- dim[dim$id == id, ]; d <- d[order(d$d), ]
  lines <- c(lines, sprintf("%s & %s \\\\", lab(d[1, ]), paste(sprintf("%.3f (%d; %.4f; %.0f)", d$MI, d$k, d$null_mean, d$eig_cond), collapse = " & "))) }
tex <- c("\\begin{tabular}{lcccccc}", "\\toprule", "Configuration & $d = 10$ & $d = 20$ & $d = 30$ & $d = 50$ & $d = 75$ & $d = 100$ \\\\", "\\midrule", lines[1:3], "\\addlinespace[2pt]", lines[4:8], "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(MS, "table_dimension.tex"))
# v_beta
lines <- c()
for (id in ord) { d <- vb[vb$id == id, ]; d <- d[order(d$v_beta), ]
  lines <- c(lines, sprintf("%s & %s \\\\", lab(d[1, ]), paste(sprintf("%.3f (%d; %.4f)", d$MI, d$k, d$null_mean), collapse = " & "))) }
seeds <- ch %>% filter(what == "seed") %>% group_by(id) %>% summarise(lo = min(MI), hi = max(MI), klo = min(k), khi = max(k)); sp <- ch %>% filter(what == "sizeprior")
lines2 <- sapply(ord, function(id) { s <- seeds[seeds$id == id, ]; b <- sp[sp$id == id, ]; sprintf("%s & %.3f to %.3f (%d to %d) & %.3f (%d) \\\\", lab(vb[vb$id == id, ][1, ]), s$lo, s$hi, s$klo, s$khi, b$MI, b$k) })
tex <- c("\\begin{tabular}{lccccccc}", "\\toprule", "Configuration & $v_\\beta = 0.1$ & $0.3$ & $1$ & $3$ & $10$ & $30$ & $100$ \\\\", "\\midrule", lines[1:3], "\\addlinespace[2pt]", lines[4:8], "\\midrule",
         "\\multicolumn{8}{l}{Chain seeds (range over five seeds, $v_\\beta = 1$) and basic-binomial size prior:} \\\\", sapply(lines2, function(l) sub("\\\\\\\\$", "& & & & \\\\\\\\", l)), "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(MS, "table_vbeta.tex"))
cat(readLines(file.path(MS, "table_dimension.tex")), sep = "\n"); cat("\n"); cat(readLines(file.path(MS, "table_vbeta.tex")), sep = "\n")

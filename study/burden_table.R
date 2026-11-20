################################################################
# burden_table.R — Supplementary Table: N90 and N95 per model and state
# (tiles, with draw-bootstrap 90% intervals, and tissue area), from
# burden/burden_Nq.csv and burden/burden_Nq_bootstrap.csv.
# Output: manuscript/table_burden.tex (tabular body)
################################################################
suppressPackageStartupMessages(library(dplyr))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); MS <- file.path(dirname(TRI), "manuscript")
BB <- c("UNI v2", "Virchow2", "H-optimus-1", "CONCH v1.5")
nq <- read.csv(file.path(TRI, "burden", "burden_Nq.csv"), stringsAsFactors = FALSE)
bs <- read.csv(file.path(TRI, "burden", "burden_Nq_bootstrap.csv"), stringsAsFactors = FALSE)
d <- nq %>% select(cohort, target, model, tile_px, I_full, rel_at_5, rel_at_25, N90, N95) %>% left_join(bs %>% select(model, state, N90_lo, N90_hi, N95_lo, N95_hi), by = c("model", "target" = "state")) %>%
  mutate(a = (tile_px * 0.5 / 1000)^2, A90 = N90 * a, A95 = N95 * a) %>% arrange(factor(target, c("MSI", "ERG", "ETS")), factor(model, BB))
st_lab <- c(MSI = "CRC/MSI", ERG = "PRAD/ERG", ETS = "PRAD/ETS")
lines <- c(); last <- ""
for (i in seq_len(nrow(d))) with(d[i, ], { lines <<- c(lines, sprintf("%s & %s & %d & %.3f & %.2f & %.2f & %d [%d, %d] & %d [%d, %d] & %.2f & %.2f \\\\",
  if (target != last) st_lab[[target]] else "", model, tile_px, I_full, rel_at_5, rel_at_25, N90, N90_lo, N90_hi, N95, N95_lo, N95_hi, A90, A95)); last <<- target })
tex <- c("\\begin{tabular}{llcccccccc}", "\\toprule",
  "State & Model & Tile (px) & $I_{\\mathrm{full}}$ & $I_5/I_{\\mathrm{full}}$ & $I_{25}/I_{\\mathrm{full}}$ & $N_{90}$ & $N_{95}$ & $A_{90}$ (mm$^2$) & $A_{95}$ (mm$^2$) \\\\",
  "\\midrule", lines, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(MS, "table_burden.tex")); cat(tex, sep = "\n")

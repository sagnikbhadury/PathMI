# Validation of pl_grid.R in slide mode: observed estimates must equal the stored pipeline estimates.
suppressPackageStartupMessages(library(readxl))
R0 <- "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/entropy/Research - DOE Horseshoe"
new <- read.csv(file.path(R0, "TRIDENT/patient_level/slide/crc_hoptimus1/mi_summary.csv"))
old <- readRDS(file.path(R0, "TRIDENT/crc_kpca_subho_hoptimus1/kpca_subho_results_crc_hoptimus1.rds"))
new$MI_stored <- sapply(new$key, function(k) old[[k]]$result$details$MSI$MI)
new$diff <- new$MI_obs - new$MI_stored
print(new[, c("key", "MI_obs", "MI_stored", "diff", "k_star")], digits = 6)
cat(sprintf("max |diff| = %.2e over %d cells\n", max(abs(new$diff)), nrow(new)))

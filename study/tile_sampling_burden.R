################################################################
# tile_sampling_burden.R
#
# Tissue-sampling burden of recoverability, from the saved tile
# scores of the mean-pooling / linear-PCA H-optimus-1 arm.
#
# Under that arm the fitted slide logit is a0 + mean over tiles of
# the tile contribution g (tile_maps_hoptimus1.R). Holding the fitted
# model fixed, we recompute every slide's logit from a random subset
# of m of its tiles, form the plug-in difference-of-entropies estimate
# over patients, and repeat over R random subsets. The curve of the
# estimate against m, relative to its value with all tiles, is the
# evaluation-time sampling burden: how many tiles a slide must
# contribute for the model to recover the information it recovers
# from the whole slide. It is task dependent because the tile-level
# specific information is spread differently over the tissue in the
# three states (Results 3.9).
#
# Output: TRIDENT/tile_maps_hoptimus1/sampling_burden.csv (per task,
# m, mean and 5/95% quantiles of the plug-in estimate and of the
# in-sample AUC over R subsets), the m at which the curve reaches 90
# and 95 percent of its full-tile value, and a supplementary figure
# figures/figS_sampling_burden.pdf (copied to manuscript/figures).
################################################################
rm(list = ls()); set.seed(2026)
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(scales) })
TRI <- "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/entropy/Research - DOE Horseshoe/TRIDENT"
TM  <- file.path(TRI, "tile_maps_hoptimus1"); FIG <- file.path(TRI, "figures"); MS_FIG <- file.path(dirname(TRI), "manuscript", "figures")
M_GRID <- c(5, 10, 25, 50, 100, 250, 500, 1000, 2500)
R_REP  <- 20L
auc <- function(p, y) { r <- rank(p); n1 <- sum(y); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1)) }

burden <- function(cohort, tg) {
  tiles  <- readRDS(file.path(TM, cohort, sprintf("tile_scores_%s_hoptimus1.rds", cohort)))
  slides <- read.csv(file.path(TM, cohort, sprintf("slide_summary_%s_hoptimus1.csv", cohort)), stringsAsFactors = FALSE)
  fit    <- readRDS(file.path(TM, cohort, sprintf("fits_%s_hoptimus1.rds", cohort)))[[tg]]
  st <- slides[slides$target == tg, ]; y <- st$y; names(y) <- st$slide; pbar <- fit$pbar; a0 <- fit$a0
  g <- split(tiles[[paste0("g_", tg)]], tiles$slide)[st$slide]
  base <- y * log(pbar) + (1 - y) * log(1 - pbar)
  full_eta <- a0 + vapply(g, mean, numeric(1))
  full_mi  <- mean(y * log(plogis(full_eta)) + (1 - y) * log(1 - plogis(full_eta)) - base)
  full_auc <- auc(full_eta, y)
  out <- list()
  for (m in M_GRID) {
    mi_r <- auc_r <- numeric(R_REP)
    for (r in seq_len(R_REP)) {
      eta <- a0 + vapply(g, function(v) mean(if (length(v) <= m) v else v[sample.int(length(v), m)]), numeric(1))
      p <- plogis(eta)
      mi_r[r] <- mean(y * log(p) + (1 - y) * log(1 - p) - base); auc_r[r] <- auc(eta, y)
    }
    out[[length(out) + 1]] <- data.frame(task = paste(toupper(cohort), tg, sep = "/"), m = m,
      mi_mean = mean(mi_r), mi_lo = quantile(mi_r, .05), mi_hi = quantile(mi_r, .95),
      auc_mean = mean(auc_r), auc_lo = quantile(auc_r, .05), auc_hi = quantile(auc_r, .95),
      mi_full = full_mi, auc_full = full_auc, frac_slides_with_m = mean(lengths(g) >= m))
  }
  do.call(rbind, out)
}
res <- rbind(burden("crc", "MSI"), burden("prad", "ERG"), burden("prad", "ETS"))
res$rel <- res$mi_mean / res$mi_full
write.csv(res, file.path(TM, "sampling_burden.csv"), row.names = FALSE)

cat("=== plug-in estimate by tiles per slide (fraction of the full-tile value) ===\n")
print(res %>% select(task, m, mi_mean, rel, auc_mean, frac_slides_with_m) %>% as.data.frame(), digits = 3)
cat("\n=== full-tile plug-in estimate and AUC ===\n"); print(unique(res[, c("task", "mi_full", "auc_full")]), digits = 4)
cat("\n=== tiles per slide at which the curve first reaches 90 / 95 percent of the full value (grid points) ===\n")
for (tk in unique(res$task)) { d <- res[res$task == tk, ]
  cat(sprintf("%-9s 90%%: %s | 95%%: %s | at m = 25: %.0f%% | at m = 100: %.0f%%\n", tk,
              d$m[which(d$rel >= 0.90)[1]], d$m[which(d$rel >= 0.95)[1]], 100 * d$rel[d$m == 25], 100 * d$rel[d$m == 100])) }

p <- ggplot(res, aes(x = m, y = rel, colour = task)) +
  geom_hline(yintercept = c(0.9, 1), linetype = "22", colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = mi_lo / mi_full, ymax = mi_hi / mi_full, fill = task), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  scale_x_log10(breaks = M_GRID) +
  scale_colour_manual(values = c("CRC/MSI" = "#E69F00", "PRAD/ERG" = "#56B4E9", "PRAD/ETS" = "#009E73"), name = NULL) +
  scale_fill_manual(values = c("CRC/MSI" = "#E69F00", "PRAD/ERG" = "#56B4E9", "PRAD/ETS" = "#009E73"), name = NULL) +
  labs(x = "Tiles sampled per slide (log scale)", y = "Plug-in estimate relative to all tiles") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
PDF_DEV <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
ggsave(file.path(FIG, "figS_sampling_burden.pdf"), p, device = PDF_DEV, width = 110, height = 75, units = "mm")
ggsave(file.path(FIG, "figS_sampling_burden.png"), p, width = 110, height = 75, units = "mm", dpi = 300)
if (dir.exists(MS_FIG)) file.copy(file.path(FIG, "figS_sampling_burden.pdf"), MS_FIG, overwrite = TRUE)
cat("\nfigure written\n")

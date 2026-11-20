################################################################
# figures_results_four.R
#
# Figure 2 of the manuscript (four-model version of the former
# Figures 1b, 3a and 3b, which were built on two models and are
# retired):
#   a  selected-kernel estimate against the condition number kappa
#      of the retained spectrum, 16 selected cells per task, colour
#      by model, shape by pooling; r is Pearson over the 16 cells,
#      and in brackets over the 12 cells without CONCH v1.5
#   b  incremental decomposition for CRC/MSI, 16 rows (four models
#      by four poolings): covariate-only and incremental parts of the
#      combined estimate, open diamond = embedding-only estimate
#   c  discrepancies between quantities the design requires to be
#      identical: replicate refit (CRC, 16) and null covariate block
#      (PRAD, 32), by model; dashed lines = reproducibility bound over
#      the three vision-only models (0.0082) and over all four (0.0202)
#
# Reads the three merged tables in TRIDENT/. Console check prints the
# numbers quoted in Results 3.3, 3.6 and 3.7.
################################################################
rm(list = ls()); set.seed(2026)
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(scales) })
TAB_DIR <- "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/entropy/Research - DOE Horseshoe/TRIDENT"
FIG_DIR <- file.path(TAB_DIR, "figures"); MS_FIG <- file.path(dirname(TAB_DIR), "manuscript", "figures")
W_MM <- 180; DPI <- 600; BASE <- 9
PDF_DEV <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
BB_LV <- c("UNI v2", "Virchow2", "H-optimus-1", "CONCH v1.5")
COL_BB <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), BB_LV)
harm <- function(e) case_when(grepl("^UNI", e) ~ BB_LV[1], grepl("^VIRCHOW", e) ~ BB_LV[2], grepl("^H-OPT", e) ~ BB_LV[3], grepl("^CONCH", e) ~ BB_LV[4], TRUE ~ NA_character_)
POOL_LV <- c(mean = "Mean", norm = "Norm", moment = "Moment", attention = "Attention")
theme_pathmi <- function(base_size = BASE) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          panel.border = element_rect(colour = "grey40", linewidth = 0.3),
          strip.background = element_rect(fill = "grey92", colour = "grey40", linewidth = 0.3),
          strip.text = element_text(size = base_size - 0.5, face = "bold", margin = margin(2, 2, 2, 2)),
          axis.text = element_text(size = base_size - 1, colour = "black"), axis.title = element_text(size = base_size),
          legend.position = "bottom", legend.title = element_text(size = base_size - 1),
          legend.text = element_text(size = base_size - 1), legend.key.size = unit(3.2, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1.5, "mm"),
          plot.tag = element_text(size = base_size + 1, face = "bold"))
}
best <- read_excel(file.path(TAB_DIR, "best kpca summary.xlsx")) %>%
  mutate(bb = factor(harm(Embeddings), BB_LV), task = factor(paste0(Disease, " / ", target), c("CRC / MSI", "PRAD / ERG", "PRAD / ETS")),
         pooling = factor(POOL_LV[tolower(agg_method)], unname(POOL_LV)))
inc <- read_excel(file.path(TAB_DIR, "incremental mi summary.xlsx")) %>%
  mutate(bb = factor(harm(Embeddings), BB_LV), task = factor(paste0(Disease, " / ", target), c("CRC / MSI", "PRAD / ERG", "PRAD / ETS")),
         pooling = factor(POOL_LV[tolower(agg_method)], unname(POOL_LV)))
stopifnot(nrow(best) == 48, nrow(inc) == 48, !anyNA(best$bb), !anyNA(inc$bb))

# ---- a: conditioning ----------------------------------------------------------
rlab <- best %>% group_by(task) %>% summarise(r4 = cor(MI_obs, eig_cond), r3 = cor(MI_obs[bb != "CONCH v1.5"], eig_cond[bb != "CONCH v1.5"]), .groups = "drop") %>%
  mutate(lab = sprintf("r = %+.2f  (%+.2f without CONCH)", r4, r3))
p2a <- ggplot(best, aes(eig_cond, MI_obs)) +
  geom_point(aes(colour = bb, shape = pooling), size = 1.9, stroke = 0.5) +
  geom_text(data = rlab, aes(x = Inf, y = -Inf, label = lab), hjust = 1.05, vjust = -0.6, size = 2.3, inherit.aes = FALSE) +
  facet_wrap(~task, nrow = 1, scales = "free_x") + scale_x_log10() +
  scale_colour_manual(values = COL_BB, name = "Foundation model") + scale_shape_manual(values = c(16, 17, 15, 18), name = "Pooling") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = expression("Condition number of the retained spectrum, " * kappa * " (log scale)"), y = "Mutual information (nats)") +
  theme_pathmi() + theme(panel.spacing.x = unit(3, "mm")) + guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))

# ---- b: incremental, CRC, 16 rows ---------------------------------------------
crc <- inc %>% filter(Disease == "CRC") %>% arrange(bb, pooling) %>%
  mutate(row = paste(bb, pooling, sep = "  |  "), row = factor(row, rev(unique(row))), pct_cov = 100 * covar_only / combined)
stopifnot(nrow(crc) == 16)
crc_long <- crc %>% select(row, bb, covar_only, incremental_subho_mi) %>%
  pivot_longer(c(covar_only, incremental_subho_mi), names_to = "part", values_to = "mi") %>%
  mutate(part = factor(part, c("covar_only", "incremental_subho_mi"), c("Covariate only", "Incremental (embedding given covariates)")))
p2b <- ggplot(crc_long, aes(x = mi, y = row, fill = part)) +
  geom_col(width = 0.66, colour = "white", linewidth = 0.25) +
  geom_point(data = crc, aes(x = emb_only, y = row), inherit.aes = FALSE, shape = 23, size = 1.8, stroke = 0.45, fill = "white", colour = "black") +
  geom_text(data = crc, aes(x = pmax(combined, emb_only), y = row, label = sprintf("%.0f%%", pct_cov)), inherit.aes = FALSE, hjust = -0.25, size = 2.1) +
  scale_fill_manual(values = c("Covariate only" = "#DDCC77", "Incremental (embedding given covariates)" = "#117733"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Mutual information (nats)", y = NULL) +
  theme_pathmi() + theme(axis.text.y = element_text(size = BASE - 2.5)) + guides(fill = guide_legend(nrow = 1))

# ---- c: reproducibility -------------------------------------------------------
# 19 Nov 2026 (QC): the replicate-refit identity is taken over both cohorts (48 arms); it was CRC only (16),
# which left out a PRAD discrepancy of 0.0095 nats (H-optimus-1, ETS, cosine/norm) above the old 0.0082 bound.
d1 <- inc %>%
  inner_join(best %>% select(Disease, bb, target, agg_method, MI_obs), by = c("Disease", "bb", "target", "agg_method")) %>%
  transmute(task, bb, pooling, d = emb_only - MI_obs, src = "Replicate refit (48)")
d2 <- inc %>% filter(covar_only == 0) %>% transmute(task, bb, pooling, d = combined - emb_only, src = "Null covariate block (PRAD, 32)")
stopifnot(nrow(d1) == 48, nrow(d2) == 32)
# 20 Nov 2026 (patient-level rerun): third identity, the full-data refit of every cell after the projection is
# recomputed by an independent routine (column signs may flip), from usable/usable_infold.csv (rep 0)
uf <- read.csv(file.path(TAB_DIR, "usable", "usable_infold.csv"), stringsAsFactors = FALSE); uf <- uf[!uf$control, ]
d3 <- data.frame(task = paste0(uf$cohort, " / ", uf$target), bb = uf$model, pooling = tools::toTitleCase(uf$pooling), d = uf$MI_refit - uf$MI_obs, src = "Sign-flip refit (192)")
disc <- bind_rows(mutate(d1, bb = as.character(bb), pooling = as.character(pooling)), mutate(d2, bb = as.character(bb), pooling = as.character(pooling)), d3) %>%
  mutate(bb = factor(bb, BB_LV), pooling = factor(pooling, c("Mean", "Norm", "Moment", "Attention")))
b3 <- max(abs(disc$d[disc$bb != "CONCH v1.5"])); b4 <- max(abs(disc$d))
p2c <- ggplot(disc, aes(x = d, y = src)) +
  annotate("rect", xmin = -b3, xmax = b3, ymin = -Inf, ymax = Inf, fill = "grey88", alpha = 0.7) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.35) +
  geom_vline(xintercept = c(-b3, b3), colour = "grey55", linewidth = 0.3, linetype = "22") +
  geom_vline(xintercept = c(-b4, b4), colour = COL_BB[["CONCH v1.5"]], linewidth = 0.3, linetype = "22") +
  geom_point(aes(colour = bb, shape = pooling), position = position_jitter(height = 0.16, width = 0, seed = 1), size = 1.9, stroke = 0.5) +
  annotate("text", x = b3, y = Inf, vjust = 1.6, hjust = 1.05, label = sprintf("%.4f", b3), size = 2.2, colour = "grey25") +
  annotate("text", x = -b4, y = Inf, vjust = 1.6, hjust = -0.1, label = sprintf("%.4f", b4), size = 2.2, colour = COL_BB[["CONCH v1.5"]]) +
  scale_colour_manual(values = COL_BB, name = "Foundation model") + scale_shape_manual(values = c(16, 17, 15, 18), name = "Pooling") +
  scale_x_continuous(labels = label_number(accuracy = 0.001)) + scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  labs(x = "Discrepancy between quantities\nrequired to be identical (nats)", y = NULL) +
  theme_pathmi() + guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))

fig2 <- p2a / (p2b | p2c) + plot_layout(heights = c(1, 1.35), guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0.5, "mm"))
ggsave(file.path(FIG_DIR, "fig2_conditioning_repro.pdf"), fig2, device = PDF_DEV, width = W_MM, height = 175, units = "mm")
ggsave(file.path(FIG_DIR, "fig2_conditioning_repro.png"), fig2, width = W_MM, height = 175, units = "mm", dpi = DPI)
if (dir.exists(MS_FIG)) file.copy(file.path(FIG_DIR, "fig2_conditioning_repro.pdf"), MS_FIG, overwrite = TRUE)

cat("=== r(MI, kappa) per task: four models / without CONCH ===\n"); print(as.data.frame(rlab[, c("task", "r4", "r3")]), digits = 3)
cat("=== kappa by task x model (selected kernels, mean) ===\n"); print(round(with(best, tapply(eig_cond, list(task, bb), mean)), 1))
cat("=== CRC covariate share (percent) by model ===\n"); print(crc %>% group_by(bb) %>% summarise(lo = min(pct_cov), hi = max(pct_cov)) %>% as.data.frame(), digits = 3)
cat("=== bounds: three vision-only models", round(b3, 4), "| all four", round(b4, 4), "===\n")
cat("nonzero discrepancies:\n"); print(disc %>% filter(abs(d) > 1e-9) %>% arrange(d) %>% as.data.frame(), digits = 3)
cat("Figure written to", FIG_DIR, "\n")

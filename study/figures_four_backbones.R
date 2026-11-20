################################################################
# figures_four_backbones.R
#
# Figure 4 of the manuscript: the four-backbone comparison on the
# 192-cell grid (UNI v2, Virchow2, H-optimus-1, CONCH v1.5) and the
# configuration-level question "is there one configuration to use
# in general?".  Reads the three merged xlsx tables in TRIDENT/ and
# recomputes every statistic quoted in Results 3.8 and the Data
# analysis section, so figure and text come from one computation
# (the same statistics are printed by config_comparison_tests.R).
#
#   a  all 192 estimates by task and backbone; grey band = null mean
#      + 2 null s.d. (task average); black bar = median of the 16
#      cells of that backbone
#   b  mean within-task rank of each backbone over the 16 matched
#      (decomposition, pooling) cells; 1 = best of four; open circle
#      = pooled over the 48 blocks
#   c  the 64 (backbone, decomposition, pooling) configurations:
#      mean rank across tasks (1 = best of 64) against the largest
#      shortfall from the task's best estimate; dashed line = the
#      0.0082 nat reproducibility bound of Section 3.7
#
# Include with \includegraphics[width=\textwidth]{figures/fig4_four_backbones.pdf}
# Dependencies: readxl, dplyr, tidyr, ggplot2, patchwork, scales
################################################################

rm(list = ls()); set.seed(2026)
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(scales) })

TAB_DIR <- "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/entropy/Research - DOE Horseshoe/TRIDENT"
FIG_DIR <- file.path(TAB_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)
MS_FIG  <- file.path(dirname(TAB_DIR), "manuscript", "figures")
W_MM <- 180; DPI <- 600; BASE <- 9
PDF_DEV <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf

BB_LV <- c("UNI v2", "Virchow2", "H-optimus-1", "CONCH v1.5")
COL_BB <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), BB_LV)
harm <- function(e) case_when(grepl("^UNI", e) ~ BB_LV[1], grepl("^VIRCHOW", e) ~ BB_LV[2],
                              grepl("^H-OPT", e) ~ BB_LV[3], grepl("^CONCH", e) ~ BB_LV[4], TRUE ~ NA_character_)
KERN_LV <- c(pca = "Linear PCA", cosine = "Cosine", rbf = "RBF", poly = "Poly")
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

grid <- read_excel(file.path(TAB_DIR, "mi summary table.xlsx")) %>%
  mutate(bb = factor(harm(Embeddings), levels = BB_LV),
         task = factor(paste0(Disease, " / ", target), levels = c("CRC / MSI", "PRAD / ERG", "PRAD / ETS")),
         kernel = factor(KERN_LV[tolower(decomp)], levels = unname(KERN_LV)),
         pooling = factor(POOL_LV[tolower(agg_method)], levels = unname(POOL_LV)),
         block = paste(task, kernel, pooling), sig = p_value < 0.05)
stopifnot(nrow(grid) == 192L, !anyNA(grid$bb), !anyNA(grid$kernel), !anyNA(grid$pooling))

# ---- statistics ---------------------------------------------------------------
ranks_task <- grid %>% group_by(task, kernel, pooling) %>% mutate(r = rank(-MI_obs)) %>% ungroup()
mr_task <- ranks_task %>% group_by(task, bb) %>% summarise(mean_rank = mean(r), .groups = "drop")
mr_all  <- ranks_task %>% group_by(bb) %>% summarise(mean_rank = mean(r), .groups = "drop") %>% mutate(task = "all tasks")
fried <- grid %>% group_by(task) %>%
  summarise(p = friedman.test(MI_obs, bb, factor(paste(kernel, pooling)))$p.value, .groups = "drop")
fried_all <- friedman.test(grid$MI_obs, grid$bb, factor(grid$block))
W_kendall <- function(R) { m <- ncol(R); n <- nrow(R); 12 * sum((rowSums(R) - m * (n + 1) / 2)^2) / (m^2 * (n^3 - n)) }
Rmat <- sapply(levels(grid$task), function(tk) rank(mr_task$mean_rank[mr_task$task == tk][match(BB_LV, mr_task$bb[mr_task$task == tk])]))
W_bb <- W_kendall(Rmat)
wide <- grid %>% select(block, bb, MI_obs) %>% pivot_wider(names_from = bb, values_from = MI_obs)
pw <- combn(BB_LV, 2, function(p) { t <- wilcox.test(wide[[p[1]]], wide[[p[2]]], paired = TRUE, exact = FALSE)
  data.frame(a = p[1], b = p[2], a_higher = sum(wide[[p[1]]] > wide[[p[2]]]), median_diff = median(wide[[p[1]]] - wide[[p[2]]]), p = t$p.value) },
  simplify = FALSE) %>% bind_rows() %>% mutate(p_holm = p.adjust(p, "holm"))
fit <- lm(MI_obs ~ (task + bb + kernel + pooling)^2, data = grid); an <- anova(fit)

cfg <- grid %>% group_by(task) %>% mutate(rank64 = rank(-MI_obs), gap = max(MI_obs) - MI_obs) %>% ungroup() %>%
  group_by(bb, kernel, pooling) %>% summarise(mean_rank = mean(rank64), worst_rank = max(rank64), max_gap = max(gap), .groups = "drop")
BOUND <- 0.0163   # 20 Nov 2026: max over the three identities (replicate refit, null covariate block, sign-flip refit), vision-only models
nullband <- grid %>% group_by(task) %>% summarise(hi = mean(null_mean) + 2 * mean(null_sd), .groups = "drop")

# ---- panel a ------------------------------------------------------------------
med <- grid %>% group_by(task, bb) %>% summarise(m = median(MI_obs), .groups = "drop")
p4a <- ggplot(grid, aes(x = bb, y = MI_obs)) +
  geom_rect(data = nullband, aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = hi), fill = "grey85", alpha = 0.65, inherit.aes = FALSE) +
  geom_jitter(aes(colour = kernel, shape = pooling), width = 0.18, height = 0, size = 1.5, stroke = 0.45) +
  geom_segment(data = med, aes(x = as.numeric(bb) - 0.32, xend = as.numeric(bb) + 0.32, y = m, yend = m), linewidth = 0.55, colour = "black") +
  facet_wrap(~task, nrow = 1) +
  scale_colour_manual(values = c("Linear PCA" = "#332288", "Cosine" = "#44AA99", "RBF" = "#CC6677", "Poly" = "#AA4499"), name = "Decomposition") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Pooling") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
  labs(x = NULL, y = "Mutual information (nats)") +
  theme_pathmi() + theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.spacing.x = unit(2.5, "mm"))

# ---- panel b ------------------------------------------------------------------
mr_plot <- bind_rows(mr_task %>% mutate(task = as.character(task)), mr_all) %>%
  mutate(task = factor(task, levels = c(levels(grid$task), "all tasks")), bb = factor(bb, levels = rev(BB_LV)))
flab <- fried %>% mutate(lab = sprintf("Friedman p = %s", ifelse(p < 1e-4, "< 0.0001", sprintf("%.4f", p))))
p4b <- ggplot(mr_plot, aes(x = mean_rank, y = bb)) +
  geom_vline(xintercept = 2.5, colour = "grey75", linewidth = 0.3, linetype = "22") +
  geom_point(data = filter(mr_plot, task != "all tasks"), aes(colour = task), size = 2, position = position_dodge(width = 0.5)) +
  geom_point(data = filter(mr_plot, task == "all tasks"), shape = 21, size = 2.6, stroke = 0.7, fill = "white", colour = "black") +
  scale_colour_manual(values = c("CRC / MSI" = "#E69F00", "PRAD / ERG" = "#56B4E9", "PRAD / ETS" = "#009E73"), name = "Task") +
  scale_x_continuous(limits = c(1, 4), breaks = 1:4) +
  labs(x = "Mean within-task rank (1 = best of four)", y = NULL,
       caption = sprintf("Open circle, all 48 matched blocks (Friedman p < 0.0001).\nKendall's W across tasks = %.2f.", W_bb)) +
  theme_pathmi() + theme(plot.caption = element_text(size = BASE - 2, hjust = 0, colour = "grey30")) +
  guides(colour = guide_legend(nrow = 1))

# ---- panel c ------------------------------------------------------------------
hl <- cfg %>% filter((bb == "H-optimus-1" & kernel == "Cosine" & pooling == "Moment") |
                     (bb == "H-optimus-1" & kernel == "Linear PCA" & pooling == "Moment"))
p4c <- ggplot(cfg, aes(x = mean_rank, y = max_gap)) +
  geom_hline(yintercept = BOUND, linetype = "22", colour = "grey45", linewidth = 0.35) +
  geom_point(aes(colour = bb, shape = kernel), size = 1.7, stroke = 0.5) +
  geom_text(data = hl, aes(label = paste(kernel, pooling, sep = " / ")), size = 2.1, hjust = -0.15,
            vjust = ifelse(hl$kernel == "Cosine", -0.6, 1.5)) +
  annotate("text", x = Inf, y = BOUND, label = sprintf("reproducibility bound %.4f", BOUND), hjust = 1.02, vjust = -0.4, size = 2.2, colour = "grey30") +
  scale_colour_manual(values = COL_BB, name = "Foundation model") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Decomposition") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Mean rank across the three tasks (1 = best of 64)", y = "Largest shortfall from the task's best (nats)") +
  theme_pathmi() + theme(legend.box = "vertical", legend.spacing.y = unit(0.5, "mm")) +
  guides(colour = guide_legend(order = 1, nrow = 2, title.position = "left"),
         shape  = guide_legend(order = 2, nrow = 2, title.position = "left"))

fig4 <- p4a / (p4b | p4c) + plot_layout(heights = c(1.15, 1.05)) + plot_annotation(tag_levels = "a")
ggsave(file.path(FIG_DIR, "fig1_four_models.pdf"), fig4, device = PDF_DEV, width = W_MM, height = 172, units = "mm")
ggsave(file.path(FIG_DIR, "fig1_four_models.png"), fig4, width = W_MM, height = 172, units = "mm", dpi = DPI)
if (dir.exists(MS_FIG)) file.copy(file.path(FIG_DIR, "fig1_four_models.pdf"), MS_FIG, overwrite = TRUE)

# ---- console check: every number quoted in the text ---------------------------
cat("\n=== cells above null by task x backbone ===\n"); print(with(grid, tapply(sig, list(task, bb), sum)))
cat("\n=== grid max MI by task x backbone ===\n"); print(round(with(grid, tapply(MI_obs, list(task, bb), max)), 4))
cat("\n=== ANOVA (two-way interactions; residual = higher order) ===\n"); print(round(an[, c("Df", "F value", "Pr(>F)")], 4))
cat("\n=== Friedman: backbone over 48 blocks ===\n"); print(fried_all)
cat("mean rank, all blocks:\n"); print(mr_all)
cat("\n=== mean rank within task ===\n"); print(as.data.frame(mr_task %>% pivot_wider(names_from = bb, values_from = mean_rank)), digits = 3)
print(as.data.frame(fried), digits = 3); cat(sprintf("Kendall W (backbones across tasks) = %.3f\n", W_bb))
cat("\n=== pairwise Wilcoxon signed-rank, 48 blocks, Holm ===\n"); print(pw, digits = 3, row.names = FALSE)
cat("\n=== decomposition and pooling Friedman over 48 blocks ===\n")
for (f in c("kernel", "pooling")) { blk <- apply(grid[, setdiff(c("task", "bb", "kernel", "pooling"), f)], 1, paste, collapse = "|")
  ft <- friedman.test(grid$MI_obs, grid[[f]], factor(blk))
  mr <- grid %>% mutate(blk = blk) %>% group_by(blk) %>% mutate(r = rank(-MI_obs)) %>% group_by(.data[[f]]) %>% summarise(mr = mean(r), .groups = "drop")
  cat(sprintf("%s: chi2 = %.2f, p = %.2e | %s\n", f, ft$statistic, ft$p.value, paste(sprintf("%s %.2f", mr[[1]], mr$mr), collapse = ", "))) }
cat("\n=== configurations ===\n")
cat(sprintf("within bound in all tasks: %d | within 2 null sd in all tasks: %d\n", sum(cfg$max_gap <= BOUND),
            sum((grid %>% group_by(task) %>% mutate(g = (max(MI_obs) - MI_obs) / mean(null_sd)) %>% group_by(bb, kernel, pooling) %>% summarise(m = max(g), .groups = "drop"))$m <= 2)))
print(as.data.frame(cfg %>% arrange(mean_rank) %>% head(5)), digits = 3)
print(as.data.frame(cfg %>% arrange(max_gap) %>% head(5)), digits = 3)
cat("\nFigure written to", FIG_DIR, "\n")

################################################################
# PathMI manuscript figures  (v4)
#
# Changes from v3:
#   - Fig 1a: x axis breaks capped at 0.10 and panel spacing widened, so
#     facet boundary labels no longer collide ("0.150.00").
#   - Fig 1a: sign test annotation moved into expanded blank space below the
#     Poly block. Top right does not work; the Linear PCA arrows reach into it.
#   - Fig 1a: "p = < 0.0001" double symbol fixed.
#   - Fig 2a: stack order reversed so Decomposition reads leftmost. Legend
#     order preserved via guide reverse.
#   - Model display labels lifted into one config block (FM_UNI / FM_VIR).
#
# Include in LaTeX with:
#   \includegraphics[width=\textwidth]{figures/fig1_fm_contrast.pdf}
# Do NOT pass scale= or a fractional width. The figures are built at final
# print size; any scaling shrinks type below 7 pt.
#
# Dependencies: readxl, dplyr, tidyr, ggplot2, patchwork, scales
################################################################

rm(list = ls())
set.seed(2026)


suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---- Config -----------------------------------------------------------------
TAB_DIR <- file.path(
  "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury",
  "CODES/SBTurboCodes/entropy/Research - DOE Horseshoe/TRIDENT"
)
FIG_DIR <- file.path(TAB_DIR, "figures")

stopifnot(dir.exists(TAB_DIR))
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(FIG_DIR))

# Foundation model display labels. Change here only; everything downstream
# derives from these two strings.
#   Note: UNI and UNI v2 are different models (1024 vs 1536 dimensions), as are
#   Virchow and Virchow2. Dropping the version identifier makes the compared
#   objects ambiguous. If shortening for legend width is the goal, prefer short
#   axis labels with full names defined in the caption and methods.
FM_UNI <- "UNI v2"
FM_VIR <- "Virchow2"
FM_LV  <- c(FM_UNI, FM_VIR)

W_MM <- 180     # double column width
DPI  <- 600
BASE <- 9       # base font size in pt

# Colourblind safe.
COL_FM <- setNames(c("#0072B2", "#D55E00"), FM_LV)
COL_UP <- "#0072B2"   # Virchow higher
COL_DN <- "#D55E00"   # UNI higher

PDF_DEV <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf

# ---- Resolve file paths -----------------------------------------------------
# On-disk names use spaces, not underscores. Patterns are anchored so that
# "^mi[ _]summary" does not also catch "incremental mi summary.xlsx", and
# "^best[ _]kpca" does not catch "kpca selected decomposition.xlsx" (which is
# byte identical to the best kpca table and must not be silently substituted).
pick <- function(pat) {
  hits <- list.files(TAB_DIR, pattern = pat, full.names = TRUE,
                     ignore.case = TRUE)
  if (length(hits) != 1L) {
    stop(sprintf("Pattern '%s' matched %d files:\n  %s",
                 pat, length(hits),
                 paste(basename(hits), collapse = "\n  ")), call. = FALSE)
  }
  hits
}

F_GRID <- pick("^mi[ _]summary")        # 96 row full factorial grid
F_BEST <- pick("^best[ _]kpca")         # 24 row selected kernel table
F_INC  <- pick("^incremental[ _]mi")    # 24 row covariate decomposition

cat("Resolved input files:\n  ",
    paste(basename(c(F_GRID, F_BEST, F_INC)), collapse = "\n   "), "\n\n")

# ---- Shared theme -----------------------------------------------------------
theme_pathmi <- function(base_size = BASE) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border       = element_rect(colour = "grey40", linewidth = 0.3),
      strip.background   = element_rect(fill = "grey92", colour = "grey40",
                                        linewidth = 0.3),
      strip.text         = element_text(size = base_size - 0.5, face = "bold",
                                        margin = margin(2, 2, 2, 2)),
      axis.text          = element_text(size = base_size - 1, colour = "black"),
      axis.title         = element_text(size = base_size),
      legend.position    = "bottom",
      legend.title       = element_text(size = base_size - 1),
      legend.text        = element_text(size = base_size - 1),
      legend.key.size    = unit(3.2, "mm"),
      legend.margin      = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1.5, "mm"),
      plot.tag           = element_text(size = base_size + 2, face = "bold")
    )
}

# ---- Load and harmonise -----------------------------------------------------
# Embeddings is coded "VIRCHOW V2" in the grid table and "VIRCHOW" in the
# other two. Harmonise before any join.
harmonise <- function(d) {
  d %>%
    mutate(
      Embeddings = case_when(
        grepl("^VIRCHOW", Embeddings) ~ FM_VIR,
        grepl("^UNI",     Embeddings) ~ FM_UNI,
        TRUE ~ NA_character_
      ),
      Embeddings = factor(Embeddings, levels = FM_LV),
      task = factor(paste0(Disease, " / ", target),
                    levels = c("CRC / MSI", "PRAD / ERG", "PRAD / ETS"))
    )
}

grid <- read_excel(F_GRID) %>% harmonise()
best <- read_excel(F_BEST) %>% harmonise()
inc  <- read_excel(F_INC)  %>% harmonise()

# The tables disagree on case: grid$decomp is "PCA" but every key prefix is
# "pca". Key the lookups on lowercase and normalise at every call site.
KERN_LV <- c(pca = "Linear PCA", cosine = "Cosine", rbf = "RBF", poly = "Poly")
POOL_LV <- c(mean = "Mean", norm = "Norm", moment = "Moment",
             attention = "Attention")

grid <- grid %>%
  mutate(kernel  = factor(KERN_LV[tolower(decomp)],     levels = unname(KERN_LV)),
         pooling = factor(POOL_LV[tolower(agg_method)], levels = unname(POOL_LV)))

best <- best %>%
  mutate(kern    = tolower(sub("_.*$", "", key)),
         kernel  = factor(KERN_LV[kern],                levels = unname(KERN_LV)),
         pooling = factor(POOL_LV[tolower(agg_method)], levels = unname(POOL_LV)))

inc <- inc %>%
  mutate(pooling = factor(POOL_LV[tolower(agg_method)], levels = unname(POOL_LV)))

# ---- Validate ---------------------------------------------------------------
# Silent NA from a level lookup would drop cells without erroring. Catch it.
stopifnot(
  nrow(grid) == 96L, nrow(best) == 24L, nrow(inc) == 24L,
  !anyNA(grid$Embeddings), !anyNA(grid$task),
  !anyNA(grid$kernel),     !anyNA(grid$pooling),
  !anyNA(best$kernel),     !anyNA(best$pooling),
  !anyNA(inc$pooling)
)
cat("Input validation passed: 96 / 24 / 24 rows, no unmapped factor levels.\n\n")


################################################################
# FIGURE 1
#   a  Paired backbone contrast across all 16 matched cells per task
#   b  Estimate against spectral conditioning, selected kernels
# Supports Results 3.2 and 3.3
################################################################

# Fixed configuration order, identical across facets, so the reader can compare
# the three tasks row by row. Built kernel-major so Poly forms a block.
CFG_LV <- rev(as.vector(t(outer(
  levels(grid$kernel), levels(grid$pooling),
  function(k, p) paste(k, p, sep = " / ")
))))
stopifnot(length(CFG_LV) == 16L, !anyDuplicated(CFG_LV))

LAB_UP <- paste(FM_VIR, "higher")
LAB_DN <- paste(FM_UNI, "higher")

pair <- grid %>%
  select(task, kernel, pooling, Embeddings, MI_obs) %>%
  pivot_wider(names_from = Embeddings, values_from = MI_obs) %>%
  mutate(
    delta = .data[[FM_VIR]] - .data[[FM_UNI]],
    dir   = ifelse(delta > 0, LAB_UP, LAB_DN),
    cfg   = factor(paste(kernel, pooling, sep = " / "), levels = CFG_LV)
  )
stopifnot(nrow(pair) == 48L, !anyNA(pair$cfg), !anyNA(pair$delta))

sgn <- pair %>%
  group_by(task) %>%
  summarise(
    k   = sum(delta > 0),
    n   = n(),
    p   = binom.test(sum(delta > 0), n())$p.value,
    med = median(delta),
    .groups = "drop"
  ) %>%
  # Fixed: "p = < 0.0001" carried a double symbol. The relational operator now
  # lives inside the branch rather than being prefixed unconditionally.
  mutate(lab = sprintf("%s higher %d/%d\nmedian %s%.3f, sign test %s",
                       FM_VIR, k, n, ifelse(med > 0, "+", ""), med,
                       ifelse(p < 1e-4, "p < 0.0001", sprintf("p = %.3f", p))))

p1a <- ggplot(pair, aes(y = cfg)) +
  geom_segment(aes(x = .data[[FM_UNI]], xend = .data[[FM_VIR]],
                   yend = cfg, colour = dir),
               linewidth = 0.45,
               arrow = arrow(length = unit(1.1, "mm"), type = "closed")) +
  geom_point(aes(x = .data[[FM_UNI]]), colour = COL_FM[[FM_UNI]], size = 0.9) +
  geom_point(aes(x = .data[[FM_VIR]]), colour = COL_FM[[FM_VIR]], size = 0.9) +
  # Annotation sits in blank space created below the Poly block by the y scale
  # expansion. Top right is not usable; the Linear PCA arrows reach into it.
  geom_text(data = sgn, aes(x = Inf, y = -Inf, label = lab),
            hjust = 1.02, vjust = -0.25, size = 2.2, lineheight = 0.95,
            inherit.aes = FALSE) +
  facet_wrap(~task, nrow = 1) +
  scale_colour_manual(values = setNames(c(COL_UP, COL_DN), c(LAB_UP, LAB_DN)),
                      name = NULL) +
  # Fixed: dropping the 0.15 break stops facet boundary labels overprinting.
  # Data max is ~0.144 and remains fully visible.
  scale_x_continuous(breaks = c(0, 0.05, 0.10),
                     expand = expansion(mult = c(0.04, 0.12))) +
  scale_y_discrete(expand = expansion(add = c(1.7, 0.6))) +
  labs(x = "Mutual information (nats)", y = NULL) +
  theme_pathmi() +
  theme(panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.2),
        panel.spacing.x    = unit(3, "mm"))

rlab <- best %>%
  group_by(task) %>%
  summarise(r = cor(MI_obs, eig_cond), .groups = "drop") %>%
  mutate(lab = sprintf("r = %s%.2f", ifelse(r > 0, "+", ""), r))

p1b <- ggplot(best, aes(eig_cond, MI_obs)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey55",
              linewidth = 0.4, linetype = "22", formula = y ~ x) +
  geom_point(aes(colour = Embeddings, shape = pooling), size = 1.9,
             stroke = 0.5) +
  geom_text(data = rlab, aes(x = Inf, y = -Inf, label = lab),
            hjust = 1.15, vjust = -0.6, size = 2.4, fontface = "bold",
            inherit.aes = FALSE) +
  facet_wrap(~task, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = COL_FM, name = "Foundation model") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Pooling") +
  labs(x = expression("Condition number of retained spectrum, " * kappa),
       y = "Mutual information (nats)") +
  theme_pathmi() +
  theme(panel.spacing.x = unit(3, "mm")) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))

fig1 <- p1a / p1b +
  plot_layout(heights = c(1.55, 1)) +
  plot_annotation(tag_levels = "a")

ggsave(file.path(FIG_DIR, "fig1_fm_contrast.pdf"), fig1, device = PDF_DEV,
       width = W_MM, height = 172, units = "mm")
ggsave(file.path(FIG_DIR, "fig1_fm_contrast.png"), fig1,
       width = W_MM, height = 172, units = "mm", dpi = DPI)


################################################################
# FIGURE 2
#   a  Variance in the estimate attributable to each design factor
#   b  Full grid, estimate by decomposition and pooling, against the null
# Supports Results 3.4 and 3.5
################################################################

part_one <- function(d, f) {
  gm <- mean(d$MI_obs)
  ss <- d %>%
    group_by(.data[[f]]) %>%
    summarise(k = n(), m = mean(MI_obs), .groups = "drop") %>%
    summarise(ss = sum(k * (m - gm)^2)) %>%
    pull(ss)
  100 * ss / sum((d$MI_obs - gm)^2)
}

vp <- grid %>%
  group_by(task) %>%
  group_modify(~ tibble(
    `Foundation model` = part_one(.x, "Embeddings"),
    `Decomposition`    = part_one(.x, "kernel"),
    `Pooling`          = part_one(.x, "pooling")
  )) %>%
  ungroup() %>%
  mutate(`Residual / interaction` =
           100 - `Foundation model` - Decomposition - Pooling)

# Fixed: position_stack places the first factor level rightmost, so the old
# level order buried Decomposition (the dominant factor) at the far right.
# Reversing the levels brings it to the left; guide reverse keeps the legend
# reading Decomposition, Foundation model, Pooling, Residual.
FILL_LV <- rev(c("Decomposition", "Foundation model", "Pooling",
                 "Residual / interaction"))

vp_long <- vp %>%
  pivot_longer(-task, names_to = "factor", values_to = "pct") %>%
  mutate(factor = factor(factor, levels = FILL_LV),
         task_r = factor(task, levels = rev(levels(task))))

p2a <- ggplot(vp_long, aes(x = pct, y = task_r, fill = factor)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(pct >= 6, sprintf("%.0f", pct), "")),
            position = position_stack(vjust = 0.5), size = 2.3,
            colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("Decomposition"          = "#117733",
                               "Foundation model"       = "#88CCEE",
                               "Pooling"                = "#DDCC77",
                               "Residual / interaction" = "grey75"),
                    breaks = rev(FILL_LV), name = NULL) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  labs(x = "Percent of variance in estimated mutual information", y = NULL) +
  theme_pathmi() +
  theme(panel.grid.major.x = element_line(colour = "grey94", linewidth = 0.2)) +
  guides(fill = guide_legend(nrow = 1))

nullband <- grid %>%
  group_by(task) %>%
  summarise(hi = mean(null_mean) + 2 * mean(null_sd), .groups = "drop")

p2b <- ggplot(grid, aes(x = kernel, y = MI_obs)) +
  geom_rect(data = nullband,
            aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = hi),
            fill = "grey85", alpha = 0.65, inherit.aes = FALSE) +
  geom_line(aes(group = pooling, colour = pooling), linewidth = 0.35,
            alpha = 0.75) +
  geom_point(aes(colour = pooling, shape = pooling), size = 1.7) +
  facet_grid(Embeddings ~ task) +
  scale_colour_manual(values = c("Mean" = "#332288", "Norm" = "#44AA99",
                                 "Moment" = "#CC6677", "Attention" = "#AA4499"),
                      name = "Pooling") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Pooling") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Decomposition method", y = "Mutual information (nats)",
       caption = paste("Grey band, permutation null mean plus two null",
                       "standard deviations.")) +
  theme_pathmi() +
  theme(axis.text.x     = element_text(angle = 30, hjust = 1),
        panel.spacing.x = unit(2.5, "mm"),
        plot.caption    = element_text(size = BASE - 2, hjust = 0,
                                       colour = "grey30"))

fig2 <- p2a / p2b +
  plot_layout(heights = c(1, 2.6)) +
  plot_annotation(tag_levels = "a")

ggsave(file.path(FIG_DIR, "fig2_design_effects.pdf"), fig2, device = PDF_DEV,
       width = W_MM, height = 155, units = "mm")
ggsave(file.path(FIG_DIR, "fig2_design_effects.png"), fig2,
       width = W_MM, height = 155, units = "mm", dpi = DPI)


################################################################
# FIGURE 3
#   a  Incremental decomposition, CRC
#   b  Reproducibility: discrepancies that are required to be zero
# Supports Results 3.6 and 3.7
################################################################

crc <- inc %>%
  filter(Disease == "CRC") %>%
  arrange(Embeddings, pooling) %>%
  mutate(row = paste(Embeddings, pooling, sep = "  |  "),
         row = factor(row, levels = rev(unique(row))),
         pct_cov = 100 * covar_only / combined)
stopifnot(nrow(crc) == 8L)

crc_long <- crc %>%
  select(row, covar_only, incremental_subho_mi) %>%
  pivot_longer(c(covar_only, incremental_subho_mi),
               names_to = "part", values_to = "mi") %>%
  mutate(part = factor(part,
                       levels = c("covar_only", "incremental_subho_mi"),
                       labels = c("Covariate only",
                                  "Incremental (embedding given covariates)")))

p3a <- ggplot(crc_long, aes(x = mi, y = row, fill = part)) +
  geom_col(width = 0.66, colour = "white", linewidth = 0.25) +
  geom_point(data = crc, aes(x = emb_only, y = row), inherit.aes = FALSE,
             shape = 23, size = 1.9, stroke = 0.45, fill = "white",
             colour = "black") +
  geom_text(data = crc, aes(x = combined, y = row,
                            label = sprintf("%.1f%%", pct_cov)),
            inherit.aes = FALSE, hjust = -0.25, size = 2.2) +
  scale_fill_manual(
    values = c("Covariate only" = "#DDCC77",
               "Incremental (embedding given covariates)" = "#117733"),
    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(x = "Mutual information (nats)", y = NULL,
       caption = paste("Open diamond, embedding only estimate.",
                       "Percentage, covariate share of the combined estimate.")) +
  theme_pathmi() +
  theme(plot.caption = element_text(size = BASE - 2, hjust = 0,
                                    colour = "grey30")) +
  guides(fill = guide_legend(nrow = 1))

# Source 1. emb_only (incremental table) and MI_obs (selection table) are the
# same quantity computed twice on identical data under identical settings.
d1 <- inc %>%
  filter(Disease == "CRC") %>%
  inner_join(best %>% select(Disease, Embeddings, target, agg_method, MI_obs),
             by = c("Disease", "Embeddings", "target", "agg_method")) %>%
  transmute(task, Embeddings, pooling,
            d = emb_only - MI_obs,
            src = "Replicate refit\n(emb_only minus selected MI_obs)")
stopifnot(nrow(d1) == 8L)

# Source 2. Where covar_only is exactly zero, combined must equal emb_only.
d2 <- inc %>%
  filter(covar_only == 0) %>%
  transmute(task, Embeddings, pooling,
            d = combined - emb_only,
            src = "Null covariate block\n(combined minus emb_only)")
stopifnot(nrow(d2) == 16L)

disc  <- bind_rows(d1, d2)
bound <- max(abs(disc$d))

p3b <- ggplot(disc, aes(x = d, y = src)) +
  annotate("rect", xmin = -bound, xmax = bound, ymin = -Inf, ymax = Inf,
           fill = "grey88", alpha = 0.7) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.35) +
  geom_vline(xintercept = c(-bound, bound), colour = "grey55",
             linewidth = 0.3, linetype = "22") +
  geom_point(aes(colour = Embeddings, shape = pooling),
             position = position_jitter(height = 0.14, width = 0, seed = 1),
             size = 1.9, stroke = 0.5) +
  annotate("text", x = bound, y = Inf, vjust = 1.6, hjust = 1,
           label = sprintf("reproducibility bound  +/- %.4f nats", bound),
           size = 2.3, colour = "grey25") +
  scale_colour_manual(values = COL_FM, name = "Foundation model") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Pooling") +
  scale_x_continuous(labels = label_number(accuracy = 0.001)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  labs(x = "Discrepancy between quantities required to be identical (nats)",
       y = NULL) +
  theme_pathmi() +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))

fig3 <- p3a / p3b +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "a")

ggsave(file.path(FIG_DIR, "fig3_incremental_repro.pdf"), fig3, device = PDF_DEV,
       width = W_MM, height = 145, units = "mm")
ggsave(file.path(FIG_DIR, "fig3_incremental_repro.png"), fig3,
       width = W_MM, height = 145, units = "mm", dpi = DPI)


################################################################
# Console check. Confirms the figures encode the numbers quoted in
# the results text rather than something the plotting silently changed.
################################################################
cat("\n=== Fig 1a, paired sign tests ===\n")
print(as.data.frame(sgn[, c("task", "k", "n", "med", "p")]), digits = 4)

cat("\n=== Fig 1b, r(MI, eig_cond) by task ===\n")
print(as.data.frame(rlab[, c("task", "r")]), digits = 3)

cat("\n=== Fig 2a, variance partition (percent) ===\n")
print(as.data.frame(vp), digits = 3)

cat("\n=== Fig 3a, covariate share of combined, CRC (percent) ===\n")
print(range(crc$pct_cov), digits = 3)

cat("\n=== Fig 3b, reproducibility bound (nats) ===\n")
print(bound, digits = 4)

cat("\n=== Fig 3b, nonzero discrepancies ===\n")
print(as.data.frame(disc %>% filter(abs(d) > 1e-12) %>% arrange(d)), digits = 3)

cat("\nFigures written to:\n  ", FIG_DIR, "\n")
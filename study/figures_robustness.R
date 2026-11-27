################################################################
# figures_robustness.R
#
# Figure: usable recoverability and estimator robustness.
#   a  cross-fitted usable recoverability against PathMI, 192 cells;
#   b  retained dimension sweep, eight configurations, with the
#      permutation floor at each dimension;
#   c  working-prior sweep, same configurations;
#   d  half-sampling 90% intervals of the recorded arm of every
#      model and state against the reproducibility bound;
#   e  tissue-sampling burden N90 and N95 per model and state with
#      the range over 20 random tile draws.
# Output: figures/fig5_usable_robustness.{pdf,png} (copied to manuscript/figures)
################################################################
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork); library(scales) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); FIG <- file.path(TRI, "figures"); MS_FIG <- file.path(dirname(TRI), "manuscript", "figures")
BB <- c("UNI v2", "Virchow2", "H-optimus-1", "CONCH v1.5")
pal_m <- c("UNI v2" = "#0072B2", "Virchow2" = "#D55E00", "H-optimus-1" = "#009E73", "CONCH v1.5" = "#CC79A7")
pal_t <- c("CRC/MSI" = "#E69F00", "PRAD/ERG" = "#56B4E9", "PRAD/ETS" = "#009E73")
# 27 Nov 2026: panel d uses each record's own tolerance (0.0215 for CONCH v1.5; it used 0.0163 for every model)
th <- theme_bw(base_size = 7.5) + theme(panel.grid.minor = element_blank(), legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6.5), plot.tag = element_text(face = "bold", size = 9))

us <- read.csv(file.path(TRI, "usable", "usable_recoverability.csv"), stringsAsFactors = FALSE) %>% mutate(task = paste0(cohort, "/", target), model = factor(model, BB)) %>%
  left_join(read.csv(file.path(TRI, "usable", "usable_infold.csv"), stringsAsFactors = FALSE) %>% select(id, I_usable_infold, I_usable_infold_sd, auc_infold), by = "id") %>%
  mutate(I_usable_ridge = I_usable_infold, I_usable_ridge_sd = I_usable_infold_sd, auc_oof_ridge = auc_infold)   # ridge readout with the projection refitted in each training fold
dim <- read.csv(file.path(TRI, "robustness", "dimension_sweep.csv"), stringsAsFactors = FALSE)
vb  <- read.csv(file.path(TRI, "robustness", "vbeta_sweep.csv"), stringsAsFactors = FALSE)
rec <- read.csv(file.path(TRI, "record", "recoverability_record.csv"), stringsAsFactors = FALSE) %>% filter(record_type == "measurement") %>% mutate(task = paste0(cohort, "/", state), model = factor(model, BB))
nq  <- read.csv(file.path(TRI, "burden", "burden_Nq_bootstrap.csv"), stringsAsFactors = FALSE) %>%
  mutate(task = paste0(ifelse(state == "MSI", "CRC", "PRAD"), "/", state), model = factor(model, BB), N90_min = N90_lo, N90_max = N90_hi, N95_min = N95_lo, N95_max = N95_hi)

lab_cfg <- function(d) { d$cfg <- sprintf("%s/%s, %s, %s", d$cohort, d$target, d$model, gsub("_", " ", d$key)); d$cfg <- factor(d$cfg, unique(d$cfg[order(d$role != "primary", d$cfg)])); d }
dim <- lab_cfg(dim); vb <- lab_cfg(vb)
pal_c <- setNames(c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#999999", "#0072B2", "#F0E442")[seq_along(levels(dim$cfg))], levels(dim$cfg))

# a  (recorded usable quantity: ridge readout)
r_a <- cor(us$I_usable_ridge, us$MI_obs); r_s <- cor(us$I_usable_ridge, us$MI_obs, method = "spearman")
pa <- ggplot(us, aes(MI_obs, I_usable_ridge, colour = task, shape = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey50", linewidth = .3) + geom_hline(yintercept = 0, colour = "grey70", linewidth = .3) +
  geom_point(size = 1.4, alpha = .85) + scale_colour_manual(values = pal_t, name = NULL) + scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
  annotate("text", x = 0.003, y = max(us$I_usable_ridge) * 1.04, hjust = 0, vjust = 1, size = 2.3, label = sprintf("r = %.2f, rho = %.2f", r_a, r_s)) +
  labs(x = "Recoverability (PathMI), nats", y = "Usable recoverability, nats", tag = "a") + th + theme(legend.position = "right") + guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))

# b
pb <- ggplot(dim, aes(d, MI, colour = cfg, linetype = role)) +
  geom_ribbon(aes(ymin = pmax(null_mean - 2 * null_sd, 0), ymax = null_mean + 2 * null_sd, fill = cfg), alpha = .12, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = .5) + geom_point(size = 1.2) + geom_vline(xintercept = 50, linetype = "12", colour = "grey60", linewidth = .3) +
  scale_x_continuous(breaks = c(10, 30, 50, 75, 100)) + scale_colour_manual(values = pal_c, name = NULL) + scale_fill_manual(values = pal_c) + scale_linetype_manual(values = c(primary = "solid", difficult = "42"), guide = "none") +
  labs(x = "Retained dimension d", y = "Recoverability, nats", tag = "b") + th + theme(legend.position = "none")
# c
pc <- ggplot(vb, aes(v_beta, MI, colour = cfg, linetype = role)) +
  geom_ribbon(aes(ymin = pmax(null_mean - 2 * null_sd, 0), ymax = null_mean + 2 * null_sd, fill = cfg), alpha = .12, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = .5) + geom_point(size = 1.2) + geom_vline(xintercept = 1, linetype = "12", colour = "grey60", linewidth = .3) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")) + scale_colour_manual(values = pal_c, name = NULL) + scale_fill_manual(values = pal_c) + scale_linetype_manual(values = c(primary = "solid", difficult = "42"), guide = "none") +
  labs(x = expression(paste("Working prior variance ", v[beta])), y = NULL, tag = "c") + th + theme(legend.position = "right", legend.text = element_text(size = 5.5))
# d
rec$x <- as.numeric(rec$model) + (as.numeric(factor(rec$task, names(pal_t))) - 2) * 0.22
pd <- ggplot(rec, aes(x = x, y = I_hat, colour = task)) +
  geom_errorbar(aes(ymin = ci90_lo, ymax = ci90_hi), width = 0, linewidth = .5) + geom_point(size = 1.6) +
  geom_errorbar(aes(ymin = I_hat - eps_rep, ymax = I_hat + eps_rep), width = .18, linewidth = .3, colour = "black") +
  scale_x_continuous(breaks = 1:4, labels = c("UNI v2", "Virchow2", "H-opt-1", "CONCH")) + scale_colour_manual(values = pal_t, name = NULL) +
  labs(x = NULL, y = "Recoverability, nats", tag = "d") + th + theme(legend.position = "none", axis.text.x = element_text(size = 6.5))
# e
nq_l <- bind_rows(nq %>% transmute(task, model, q = "N90", N = N90, lo = N90_min, hi = N90_max), nq %>% transmute(task, model, q = "N95", N = N95, lo = N95_min, hi = N95_max))
nq_l$x <- as.numeric(nq_l$model) + (as.numeric(factor(nq_l$task, names(pal_t))) - 2) * 0.22
pe <- ggplot(nq_l, aes(x = x, y = N, colour = task, shape = q)) + geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .4, position = position_dodge(width = .1)) +
  geom_point(size = 1.6, position = position_dodge(width = .1)) + scale_y_log10(breaks = c(5, 10, 25, 50, 100, 250, 500)) +
  scale_x_continuous(breaks = 1:4, labels = c("UNI v2", "Virchow2", "H-opt-1", "CONCH")) + scale_colour_manual(values = pal_t, name = NULL) + scale_shape_manual(values = c(N90 = 16, N95 = 1), name = NULL) +
  labs(x = NULL, y = "Tiles per slide", tag = "e") + th + theme(legend.position = "right", axis.text.x = element_text(size = 6.5))

p <- (pa + pb + pc + plot_layout(widths = c(1.15, 1, 1.25))) / (pd + pe + plot_layout(widths = c(1, 1.3))) + plot_layout(heights = c(1, 0.85))
PDF_DEV <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
ggsave(file.path(FIG, "fig5_usable_robustness.pdf"), p, device = PDF_DEV, width = 183, height = 118, units = "mm")
ggsave(file.path(FIG, "fig5_usable_robustness.png"), p, width = 183, height = 118, units = "mm", dpi = 300)
if (dir.exists(MS_FIG)) file.copy(file.path(FIG, "fig5_usable_robustness.pdf"), MS_FIG, overwrite = TRUE)
cat("figure written\n")

################################################################
# figures_pipeline_walkthrough.R
#
# Figure 1 built from the pipeline's own outputs (replaces the drawn
# schematic, user request 19 Nov 2026). One exemplar is followed
# through every stage: the H-optimus-1 colorectal MSI Record arm
# (panels b-e, g, h) and, for the tile map, the same model's linear
# arm (mean pooling, linear PCA; the arm with an exact tile
# decomposition). Each panel title states the inference drawn from it.
#   a  one MSI-high slide: thumbnail with its tiles at 20X
#   b  patients in the first two retained coordinates of the arm
#   c  horseshoe ranking: |posterior median| of the 50 coordinates
#   d  marginal-likelihood walk over the nested models; k*
#   e  the estimate against its 200-permutation null
#   f  tile-level specific surprise on the same slide (linear arm)
#   g  learning curve: lower bound minus threshold, and usable
#      recoverability, against the number of patients
#   h  the resulting Recoverability Record
# Output: figures/fig0_pipeline.{pdf,png} (copied to manuscript/figures)
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(ggplot2); library(patchwork); library(dplyr); library(jpeg); library(scales) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); ROOT <- dirname(TRI); setwd(ROOT); source("PathMI_SUBHO_logistic.R")
FIG <- file.path(TRI, "figures"); MS_FIG <- file.path(ROOT, "manuscript", "figures")
THUMB <- "M:/kulkarnr/slideflow_test/slideflow_crc/trident_feats/thumbnails/"
BASE <- 7; COL_LAB <- c(`0` = "#0072B2", `1` = "#D55E00"); ACC <- "#1B7F79"
th <- theme_bw(base_size = BASE) + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = BASE, face = "bold", hjust = 0),
  plot.subtitle = element_text(size = BASE - 0.8, colour = "grey25", hjust = 0, lineheight = 0.95), plot.tag = element_text(size = BASE + 2, face = "bold"),
  legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = BASE - 1), legend.title = element_text(size = BASE - 1))

rec <- read.csv(file.path(TRI, "record", "recoverability_record.csv"), stringsAsFactors = FALSE)
r <- rec[rec$record_type == "measurement" & rec$cohort == "CRC" & rec$model == "H-optimus-1", ]; stopifnot(nrow(r) == 1)
key <- paste0(r$projection, "_", r$pooling)
res <- readRDS(file.path(TRI, "crc_kpca_subho_hoptimus1", "kpca_subho_results_crc_hoptimus1.rds"))
arm <- res[[key]]; Z <- arm$Z_scaled; y <- as.numeric(arm$y[rownames(Z)]); det <- arm$result$details$MSI; n <- nrow(Z)
KL <- c(pca = "linear PCA", cosine = "cosine kernel PCA", rbf = "RBF kernel PCA", poly = "polynomial kernel PCA")[[r$projection]]

# refit the arm once to expose the ranking and the walk (the pipeline's seed reproduces the stored estimate)
set.seed(42); sel <- selection.logistic(X = cbind(1, Z), y = y, iter = 1L, v_beta = 1, method.weight = "beta.binomial", nmc = 2000L, burn = 500L, thin = 10L)
kst <- which.max(sel$Likelihood); MI_refit <- max((sel$LogML.only[kst] - sel$LogML.null) / n, 0)
stopifnot(abs(MI_refit - det$MI) < 1e-10)
med <- abs(sel$BetaHat[-1]); rk <- order(med, decreasing = TRUE)

# ---- a and f: the exemplar slide ---------------------------------------------------
tl <- readRDS(file.path(TRI, "tile_maps_hoptimus1", "crc", "tile_scores_crc_hoptimus1.rds"))
ss <- read.csv(file.path(TRI, "tile_maps_hoptimus1", "crc", "slide_summary_crc_hoptimus1.csv"), stringsAsFactors = FALSE)
pos <- ss[ss$y == 1, ]; ps <- unique(ss$patch_level0)[1]
comp <- tl %>% group_by(slide) %>% summarise(fill = n() * ps^2 / ((max(x) - min(x) + ps) * (max(y) - min(y) + ps)), aspect = (max(x) - min(x) + ps) / (max(y) - min(y) + ps), .groups = "drop")
okc <- comp$slide[comp$fill >= 0.35 & comp$aspect >= 0.5 & comp$aspect <= 2]
ex <- pos[pos$slide %in% okc, ]; ex <- ex[which.max(ex$p_hat), ]; d <- tl[tl$slide == ex$slide, ]
xr <- range(d$x) + c(0, ps); yr <- range(d$y) + c(0, ps); mx <- 0.03 * diff(xr); my <- 0.03 * diff(yr)
thumb <- file.path(THUMB, paste0(ex$slide, ".jpg"))
pa <- ggplot(d, aes(x + ps / 2, -(y + ps / 2))) + { if (file.exists(thumb)) annotation_raster(readJPEG(thumb), xmin = 0, xmax = ex$level0_w, ymin = -ex$level0_h, ymax = 0) } +
  geom_tile(data = d[seq(1, nrow(d), by = 25), ], width = ps * 3, height = ps * 3, fill = NA, colour = "#1B7F79", linewidth = 0.25) +
  coord_fixed(xlim = c(xr[1] - mx, xr[2] + mx), ylim = c(-(yr[2] + my), -(yr[1] - my)), expand = FALSE) + theme_void(base_size = BASE) +
  labs(title = "Tiles", subtitle = sprintf("%s tiles at 20X;\neach embedded in 1536\ncoordinates, then pooled", format(nrow(d), big.mark = ","))) +
  theme(plot.title = element_text(size = BASE, face = "bold"), plot.subtitle = element_text(size = BASE - 0.8, colour = "grey25", lineheight = 0.95),
        panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.3))
lim <- c(0, quantile(tl$KL_MSI, 0.99))
pf <- ggplot(d, aes(x + ps / 2, -(y + ps / 2))) + { if (file.exists(thumb)) annotation_raster(readJPEG(thumb), xmin = 0, xmax = ex$level0_w, ymin = -ex$level0_h, ymax = 0) } +
  geom_tile(aes(fill = KL_MSI), width = ps, height = ps, alpha = 0.85) +
  scale_fill_viridis_c(name = "Specific surprise (nats)", limits = lim, oob = squish, guide = guide_colourbar(barwidth = unit(22, "mm"), barheight = unit(2, "mm"), title.position = "top")) +
  coord_fixed(xlim = c(xr[1] - mx, xr[2] + mx), ylim = c(-(yr[2] + my), -(yr[1] - my)), expand = FALSE) + theme_void(base_size = BASE) +
  labs(title = "Where the evidence lies", subtitle = sprintf("linear arm, p = %.2f;\n%.0f%% of tiles move the\nbelief from the prevalence", ex$p_hat, 100 * mean(d$KL_MSI > 0.05))) +
  theme(plot.title = element_text(size = BASE, face = "bold"), plot.subtitle = element_text(size = BASE - 0.8, colour = "grey25", lineheight = 0.95),
        legend.position = "bottom", legend.title = element_text(size = BASE - 1.5), legend.text = element_text(size = BASE - 1.5),
        panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.3))

# ---- b: patients in the representation ------------------------------------------
db <- data.frame(z1 = Z[, rk[1]], z2 = Z[, rk[2]], y = factor(y))
pb <- ggplot(db, aes(z1, z2, colour = y)) + geom_point(size = 0.7, alpha = 0.7) +
  scale_colour_manual(values = COL_LAB, labels = c(`0` = "MSS / MSI-L", `1` = "MSI-H"), name = NULL) +
  labs(x = sprintf("coordinate %d (rank 1)", rk[1]), y = sprintf("coordinate %d (rank 2)", rk[2]), title = "Patients",
       subtitle = sprintf("%d patients projected to 50\ncoordinates; MSI-H shifts\nalong the top-ranked axes", n)) + th +
  theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0), legend.background = element_rect(fill = alpha("white", 0.7), colour = NA))

# ---- c: ranking ------------------------------------------------------------------------
dc <- data.frame(rank = seq_along(rk), med = med[rk], sel = seq_along(rk) <= kst)
pc <- ggplot(dc, aes(rank, med, fill = sel)) + geom_col(width = 0.8) + scale_fill_manual(values = c(`TRUE` = ACC, `FALSE` = "grey75"), guide = "none") +
  labs(x = "Coordinate rank", y = "|posterior median|", title = "Horseshoe ranking", subtitle = sprintf("noise coordinates shrink\nto zero; the first %d\ncarry the signal", kst)) + th

# ---- d: walk ----------------------------------------------------------------------------
dd <- data.frame(k = seq_along(sel$Likelihood), ev = (sel$LogML.only - sel$LogML.null) / n, sc = (sel$Likelihood - max(sel$Likelihood)))
pd <- ggplot(dd, aes(k, ev)) + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) + geom_line(colour = "grey30", linewidth = 0.4) + geom_point(size = 0.6, colour = "grey30") +
  geom_point(data = dd[kst, ], colour = ACC, size = 2) + annotate("text", x = kst, y = dd$ev[kst], label = sprintf("k* = %d", kst), hjust = -0.2, vjust = -0.8, size = 2.3, colour = ACC) +
  labs(x = "Coordinates in the model, k", y = "log Bayes factor / n (nats)", title = "Marginal-likelihood walk",
       subtitle = sprintf("a coordinate enters only if\nit raises the evidence;\nestimate %.3f nats at k*", det$MI)) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) + th

# ---- e: null --------------------------------------------------------------------------------
dn <- data.frame(I = det$MI.null)
pe <- ggplot(dn, aes(I)) + geom_histogram(bins = 30, fill = "grey70", colour = NA) + geom_vline(xintercept = det$MI, colour = ACC, linewidth = 0.7) +
  annotate("text", x = det$MI, y = Inf, label = sprintf("estimate\n%.3f", det$MI), hjust = 1.05, vjust = 1.3, size = 2.2, colour = ACC) +
  scale_x_continuous(limits = c(0, det$MI * 1.08)) +
  labs(x = "Estimate under permuted labels (nats)", y = "Permutations", title = "Permutation null",
       subtitle = sprintf("estimator refitted 200 times;\nnull mean %.4f, p = 1/201;\n%.0f%% of H(Y) recovered", mean(det$MI.null), 100 * r$share)) + th

# ---- g: learning curve --------------------------------------------------------------------
lc <- read.csv(file.path(TRI, "learning", "learning_summary.csv"), stringsAsFactors = FALSE)
lid <- sprintf("CRC|hoptimus1|MSI|%s", key); dl <- lc[lc$id == lid & !(lc$source == "subsample" & lc$n == floor(lc$N / 2) & !(lc$n %in% c(50, 100, 150, 200))), ]
dl$margin <- dl$LCB - dl$threshold
sc <- max(abs(dl$margin), na.rm = TRUE) / max(dl$U_mean, na.rm = TRUE)
pg <- ggplot(dl, aes(n)) + geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(aes(y = margin), colour = ACC, linewidth = 0.5) + geom_point(aes(y = margin, shape = resolved), colour = ACC, size = 1.5) +
  geom_line(aes(y = U_mean * sc), colour = "#D55E00", linewidth = 0.4, linetype = "22") + geom_point(aes(y = U_mean * sc), colour = "#D55E00", size = 0.9) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_y_continuous(name = "LCB minus threshold", sec.axis = sec_axis(~ . / sc, name = "Usable (nats)")) +
  labs(x = "Patients n", title = "How many patients", subtitle = sprintf("resolved from %s patients;\nusable recoverability\n(dashed) rises with n",
       ifelse(is.na(r$N_evidence), "no available", as.character(r$N_evidence)))) + th +
  theme(axis.title.y.left = element_text(colour = ACC), axis.title.y.right = element_text(colour = "#D55E00"))

# ---- h: the Record -------------------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
card <- data.frame(field = c("Model, state", "Projection", "Patients", "Represented", "Null (resolution)", "90% interval", "Attributable", "Usable (AUC)", "Tiles, N90", "N evidence", "Status"),
  value = c("H-optimus-1, MSI", sprintf("%s, %s", c(pca = "linear", cosine = "cosine", rbf = "RBF", poly = "poly")[[r$projection]], r$pooling), as.character(r$n), sprintf("%s (%.0f%% of H)", fmt(r$I_hat), 100 * r$share),
            sprintf("%s (%s)", fmt(r$null_mean, 4), fmt(r$eps_rep, 4)), sprintf("%s to %s", fmt(r$ci90_lo), fmt(r$ci90_hi)), fmt(r$I_adjusted),
            sprintf("%s (%s)", fmt(r$I_usable), fmt(r$auc_oof, 2)), as.character(r$N90), ifelse(is.na(r$N_evidence), "not resolved", as.character(r$N_evidence)), r$status))
card$row <- rev(seq_len(nrow(card)))
ph <- ggplot(card) + geom_text(aes(x = 0, y = row, label = field), hjust = 0, size = 2.0, colour = "grey30") +
  geom_text(aes(x = 0.5, y = row, label = value), hjust = 0, size = 2.0, fontface = "bold") +
  scale_x_continuous(limits = c(0, 1.05)) + scale_y_continuous(limits = c(0.5, nrow(card) + 0.5)) +
  labs(title = "Recoverability Record", subtitle = "the measurement with its null,\nresolution, uncertainty,\ncounterparts and burdens") +
  theme_void(base_size = BASE) + theme(plot.title = element_text(size = BASE, face = "bold"), plot.subtitle = element_text(size = BASE - 0.8, colour = "grey25", lineheight = 0.95),
                                       panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.3), plot.margin = margin(2, 4, 2, 4))

fig <- (pa | pb | pc | pd) / (pe | pf | pg | ph) + plot_annotation(tag_levels = "a") & theme(plot.margin = margin(2, 5, 2, 2))
ggsave(file.path(FIG, "fig0_pipeline.pdf"), fig, width = 183, height = 125, units = "mm", device = cairo_pdf)
ggsave(file.path(FIG, "fig0_pipeline.png"), fig, width = 183, height = 125, units = "mm", dpi = 250)
if (dir.exists(MS_FIG)) file.copy(file.path(FIG, "fig0_pipeline.pdf"), MS_FIG, overwrite = TRUE)
cat(sprintf("exemplar: %s | slide %s | k* = %d | MI %.4f | null mean %.4f | N_evidence %s\n", key, ex$slide, kst, det$MI, mean(det$MI.null), r$N_evidence))

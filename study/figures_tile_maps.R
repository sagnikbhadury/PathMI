################################################################
# figures_tile_maps.R
#
# Figures 3 and 4 of the manuscript: tile-level recoverability maps
# from the mean-pooling / linear-PCA SUBHO fit under H-optimus-1
# (tile_maps_hoptimus1.R output) and the slide-level quantities of
# the same fit. One block per molecular state, identical in layout:
#
#   row 1  four example slides chosen from the fitted probabilities
#          (positive with the highest p, negative with the lowest p,
#          positive with the lowest p = missed, negative with the
#          highest p = false alarm), each tile coloured by the
#          posterior mean of its specific surprise over the slide
#          thumbnail, cropped to the tissue extent; the title gives
#          p, the posterior sd of the slide logit, the tile count and
#          the informative share
#   row 2  left, every slide's posterior logit (90 percent bar) against
#          its informative share, by label; right, per-patient plug-in
#          difference-of-entropies terms by label with the posterior
#          median and 90 percent interval of the plug-in estimate
#
# Figure 3 = the block for CRC / MSI (tags a-c).
# Figure 4 = the blocks for PRAD / ERG (a-c) and PRAD / ETS (d-f).
# The console check prints every number quoted in Results 3.9.
# Dependencies: MASS, dplyr, ggplot2, patchwork, scales, jpeg, readxl
################################################################

rm(list = ls()); set.seed(2026)
suppressPackageStartupMessages({ library(MASS); library(dplyr); library(ggplot2); library(patchwork); library(scales); library(jpeg) })  # dplyr after MASS so select() is dplyr's

TRI <- "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/entropy/Research - DOE Horseshoe/TRIDENT"
FIG_DIR <- file.path(TRI, "figures"); MS_FIG <- file.path(dirname(TRI), "manuscript", "figures")
TM  <- file.path(TRI, "tile_maps_hoptimus1")
W_MM <- 180; DPI <- 600; BASE <- 8
PDF_DEV <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
THUMBS <- list(crc = "M:/kulkarnr/slideflow_test/slideflow_crc/trident_feats/thumbnails/",
               prad = "M:/shared_data/TCGA_PRAD/trident-feats/thumbnails/")
COL_LAB <- c(`0` = "#0072B2", `1` = "#D55E00")

# Null band (mean + 2 sd of the permutation null) of the pca_mean arm, from the 192-row grid table.
band <- list(); {
  g <- readxl::read_excel(file.path(TRI, "mi summary table.xlsx"))
  g <- g[grepl("^H-OPT", g$Embeddings) & g$key == "pca_mean", ]
  for (i in seq_len(nrow(g))) band[[paste(g$Disease[i], g$target[i])]] <- g$null_mean[i] + 2 * g$null_sd[i]
}

theme_pathmi <- function(base_size = BASE) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          panel.border = element_rect(colour = "grey40", linewidth = 0.3),
          axis.text = element_text(size = base_size - 1, colour = "black"), axis.title = element_text(size = base_size),
          legend.position = "bottom", legend.title = element_text(size = base_size - 1),
          legend.text = element_text(size = base_size - 1), legend.key.size = unit(2.8, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
          plot.tag = element_text(size = base_size + 1, face = "bold"))
}

load_cohort <- function(cohort) {
  list(tiles  = readRDS(file.path(TM, cohort, sprintf("tile_scores_%s_hoptimus1.rds", cohort))),
       slides = read.csv(file.path(TM, cohort, sprintf("slide_summary_%s_hoptimus1.csv", cohort)), stringsAsFactors = FALSE),
       fits   = readRDS(file.path(TM, cohort, sprintf("fits_%s_hoptimus1.rds", cohort))),
       fitsum = read.csv(file.path(TM, cohort, sprintf("fit_summary_%s_hoptimus1.csv", cohort)), stringsAsFactors = FALSE),
       Zs     = readRDS(file.path(TRI, if (cohort == "crc") "crc_kpca_subho_hoptimus1/kpca_subho_results_crc_hoptimus1.rds"
                                       else "TCGA_PRAD_kpca_subho_hoptimus1/kpca_subho_results_prad_hoptimus1.rds"))$pca_mean$Z_scaled)
}

# Posterior median and 90 percent interval of the plug-in estimate (the MAP plug-in is the
# optimistic end; draws are regenerated from the saved MAP and Sigma with a fixed seed).
plugin_posterior <- function(ft, Zs, S = 500, seed = 43) {
  X <- cbind(1, Zs[, ft$sel, drop = FALSE]); y <- ft$y; pbar <- ft$pbar
  set.seed(seed); B <- mvrnorm(S, c(ft$b0, ft$beta), ft$Sigma)
  P <- plogis(X %*% t(B))
  mi <- colMeans(y * log(P) + (1 - y) * log(1 - P)) - mean(y * log(pbar) + (1 - y) * log(1 - pbar))
  c(median = median(mi), lo = unname(quantile(mi, .05)), hi = unname(quantile(mi, .95)))
}
auc <- function(p, y) { r <- rank(p); n1 <- sum(y); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1)) }

# One example-slide map, cropped to the tissue extent (2 percent margin) so that no panel
# carries the empty part of the slide.
tile_panel <- function(d, ps, cohort, sid, fill_col, title, limits) {
  thumb <- file.path(THUMBS[[cohort]], paste0(sid, ".jpg"))
  w0 <- unique(d$level0_w); h0 <- unique(d$level0_h)
  xr <- range(d$x) + c(0, ps); yr <- range(d$y) + c(0, ps)
  mx <- 0.02 * diff(xr); my <- 0.02 * diff(yr)
  p <- ggplot(d, aes(x = x + ps / 2, y = -(y + ps / 2)))
  if (file.exists(thumb)) p <- p + annotation_raster(readJPEG(thumb), xmin = 0, xmax = w0, ymin = -h0, ymax = 0)
  p + geom_tile(aes(fill = .data[[fill_col]]), width = ps, height = ps, alpha = 0.8) +
    scale_fill_viridis_c(name = "Specific surprise (nats)", limits = limits, oob = squish,
                         guide = guide_colourbar(barwidth = unit(32, "mm"), barheight = unit(2.5, "mm"),
                                                 title.position = "left", title.vjust = 0.85)) +
    coord_fixed(xlim = c(xr[1] - mx, xr[2] + mx), ylim = c(-(yr[2] + my), -(yr[1] - my)), expand = FALSE) +
    labs(title = title) + theme_void(base_size = BASE) +
    theme(plot.title = element_text(size = BASE - 3, hjust = 0.5, lineheight = 0.95), legend.position = "bottom",
          legend.title = element_text(size = BASE - 1.5), legend.text = element_text(size = BASE - 1.5),
          panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.3), plot.margin = margin(1, 1, 1, 1))
}

# Example slides are chosen by fitted probability among slides with a compact tissue
# footprint (tiles fill at least 35 percent of their bounding box, whose aspect ratio is
# within a factor of two), so that a slide made of distant tissue pieces does not force a
# wide panel and starve the other three.
compactness <- function(tiles, ps) {
  tiles %>% group_by(slide) %>% summarise(
    fill = n() * ps^2 / ((max(x) - min(x) + ps) * (max(y) - min(y) + ps)),
    aspect = (max(x) - min(x) + ps) / (max(y) - min(y) + ps), .groups = "drop")
}
pick_examples <- function(st, comp) {
  ok <- comp$slide[comp$fill >= 0.35 & comp$aspect >= 0.5 & comp$aspect <= 2]
  s <- st[st$slide %in% ok, ]
  pos <- s[s$y == 1, ]; neg <- s[s$y == 0, ]
  list(tp = pos$slide[which.max(pos$p_hat)], tn = neg$slide[which.min(neg$p_hat)],
       fn = pos$slide[which.min(pos$p_hat)], fp = neg$slide[which.max(neg$p_hat)])
}

# The block for one state: returns list(maps = patchwork row, stats = patchwork row, st, pp).
target_block <- function(coh, cohort, tg, lab_pos, lab_neg, what) {
  st <- coh$slides[coh$slides$target == tg, ]
  bv <- band[[paste(toupper(cohort), tg)]]
  kl <- paste0("KL_", tg); kllo <- paste0("KLlo_", tg)
  inf <- coh$tiles %>% group_by(slide) %>% summarise(frac_informative = mean(.data[[kllo]] > bv), .groups = "drop")
  st <- st %>% left_join(inf, by = "slide")
  stp <- if ("patient" %in% names(st)) st[!duplicated(st$patient), ] else st   # one row per patient for patient-level quantities
  pp <- plugin_posterior(coh$fits[[tg]], coh$Zs)
  ps <- unique(st$patch_level0)
  comp <- compactness(coh$tiles, ps)
  ex <- pick_examples(st, comp)
  lim <- c(0, quantile(coh$tiles[[kl]], 0.99))
  maps <- lapply(names(ex), function(k) {
    sid <- ex[[k]]; r <- st[st$slide == sid, ]
    d <- coh$tiles[coh$tiles$slide == sid, ]; d$level0_w <- r$level0_w; d$level0_h <- r$level0_h
    ttl <- sprintf("%s\n%s\np = %.2f, sd %.2f\n%s tiles, %.0f%% high-surprise", what[[k]], substr(sid, 1, 12),
                   r$p_hat, r$eta_sd, format(r$n_tiles, big.mark = ","), 100 * r$frac_informative)
    tile_panel(d, ps, cohort, sid, kl, ttl, lim) })
  maps_row <- wrap_plots(maps, nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  p_slide <- ggplot(st, aes(x = frac_informative, y = eta, colour = factor(y))) +
    geom_hline(yintercept = qlogis(mean(st$y)), linetype = "22", colour = "grey50", linewidth = 0.3) +
    geom_errorbar(aes(ymin = eta - 1.645 * eta_sd, ymax = eta + 1.645 * eta_sd), width = 0, alpha = 0.35, linewidth = 0.3) +
    geom_point(size = 1.2, stroke = 0.4) +
    scale_colour_manual(values = COL_LAB, labels = c(`0` = lab_neg, `1` = lab_pos), name = NULL) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "High-surprise tiles (share of slide)", y = "Posterior logit (90% bar)") + theme_pathmi()
  p_doe <- ggplot(stp, aes(x = doe_term, fill = factor(y))) +
    geom_histogram(bins = 40, alpha = 0.6, position = "identity", colour = NA) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    annotate("rect", xmin = pp["lo"], xmax = pp["hi"], ymin = -Inf, ymax = Inf, fill = "grey30", alpha = 0.15) +
    geom_vline(xintercept = pp["median"], colour = "grey20", linewidth = 0.5) +
    annotate("text", x = -Inf, y = Inf, vjust = 1.4, hjust = -0.05, size = 2.2,
             label = sprintf("plug-in estimate %.3f\n[%.3f, %.3f]", pp["median"], pp["lo"], pp["hi"])) +
    scale_fill_manual(values = COL_LAB, labels = c(`0` = lab_neg, `1` = lab_pos), name = NULL) +
    labs(x = "Per-patient DoE term (nats)", y = "Patients") +
    theme_pathmi() + theme(panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.2))
  list(maps = maps_row, stats = (p_slide | p_doe), st = st, stp = stp, pp = pp, ex = ex, bv = bv, kl = kl, kllo = kllo,
       n_compact = sum(comp$slide %in% st$slide & comp$fill >= 0.35 & comp$aspect >= 0.5 & comp$aspect <= 2))
}

# ---- Figure 3: CRC / MSI ------------------------------------------------------
crc <- load_cohort("crc")
b_msi <- target_block(crc, "crc", "MSI", "MSI-H", "MSS / MSI-L",
                      c(tp = "MSI-H, highest p", tn = "MSS, lowest p", fn = "missed MSI-H, lowest p", fp = "MSS false alarm, highest p"))
fig3 <- (b_msi$maps / b_msi$stats) + plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = list(c("a", "", "", "", "b", "c")))
ggsave(file.path(FIG_DIR, "fig3_tile_maps_crc.pdf"), fig3, device = PDF_DEV, width = W_MM, height = 128, units = "mm")
ggsave(file.path(FIG_DIR, "fig3_tile_maps_crc.png"), fig3, width = W_MM, height = 128, units = "mm", dpi = 300)

# ---- Figure 4: PRAD / ERG and PRAD / ETS ---------------------------------------
prad <- load_cohort("prad")
b_erg <- target_block(prad, "prad", "ERG", "ERG positive", "ERG negative",
                      c(tp = "ERG positive, highest p", tn = "ERG negative, lowest p", fn = "missed positive, lowest p", fp = "false alarm, highest p"))
b_ets <- target_block(prad, "prad", "ETS", "ETS positive", "ETS negative",
                      c(tp = "ETS positive, highest p", tn = "ETS negative, lowest p", fn = "missed positive, lowest p", fp = "false alarm, highest p"))
fig4 <- (b_erg$maps / b_erg$stats / b_ets$maps / b_ets$stats) + plot_layout(heights = c(1.15, 1, 1.15, 1)) +
  plot_annotation(tag_levels = list(c("a", "", "", "", "b", "c", "d", "", "", "", "e", "f")))
ggsave(file.path(FIG_DIR, "fig4_tile_maps_prad.pdf"), fig4, device = PDF_DEV, width = W_MM, height = 186, units = "mm")
ggsave(file.path(FIG_DIR, "fig4_tile_maps_prad.png"), fig4, width = W_MM, height = 186, units = "mm", dpi = 300)
if (dir.exists(MS_FIG)) file.copy(file.path(FIG_DIR, c("fig3_tile_maps_crc.pdf", "fig4_tile_maps_prad.pdf")), MS_FIG, overwrite = TRUE)

# ---- console check ------------------------------------------------------------
report <- function(b, tg, tiles) {
  st <- b$st; stp <- b$stp; pp <- b$pp
  cat(sprintf("\n=== %s ===\n", tg))
  cat(sprintf("n = %d | prevalence %.3f | AUC of fitted p = %.3f | mean posterior sd of logit = %.3f (range %.2f-%.2f)\n",
              nrow(st), mean(st$y), auc(st$p_hat, st$y), mean(st$eta_sd), min(st$eta_sd), max(st$eta_sd)))
  cat(sprintf("plug-in: MAP %.4f | posterior median %.4f [%.4f, %.4f]\n", mean(stp$doe_term), pp["median"], pp["lo"], pp["hi"]))
  cat(sprintf("null band for the tile rule: %.4f | tiles %d | SI median %.4f, 99th pct %.4f | informative (lower 5%% > band) %.1f%%\n",
              b$bv, length(tiles[[b$kl]]), median(tiles[[b$kl]]), quantile(tiles[[b$kl]], .99), 100 * mean(tiles[[b$kllo]] > b$bv)))
  cat(sprintf("informative share per slide: median %.1f%% (IQR %.1f-%.1f), positives %.1f%%, negatives %.1f%%\n",
              100 * median(st$frac_informative), 100 * quantile(st$frac_informative, .25), 100 * quantile(st$frac_informative, .75),
              100 * median(st$frac_informative[st$y == 1]), 100 * median(st$frac_informative[st$y == 0])))
  cat(sprintf("Spearman(share, logit) = %.2f | Spearman(sd of g, sd of logit) = %.2f | DoE term medians pos %.3f neg %.3f | negative terms %.1f%%\n",
              cor(st$frac_informative, st$eta, method = "spearman"), cor(st$tile_g_sd, st$eta_sd, method = "spearman"),
              median(stp$doe_term[stp$y == 1]), median(stp$doe_term[stp$y == 0]), 100 * mean(stp$doe_term < 0)))
  for (k in names(b$ex)) { r <- st[st$slide == b$ex[[k]], ]
    cat(sprintf("  %-3s %s  y=%d p=%.3f sd=%.2f tiles=%d informative=%.1f%%\n", k, substr(b$ex[[k]], 1, 16), r$y, r$p_hat, r$eta_sd, r$n_tiles, 100 * r$frac_informative)) }
  cat(sprintf("  (examples chosen among %d of %d slides with a compact footprint)\n", b$n_compact, nrow(st)))
}
report(b_msi, "CRC / MSI", crc$tiles); report(b_erg, "PRAD / ERG", prad$tiles); report(b_ets, "PRAD / ETS", prad$tiles)
cat("\nFigures written to", FIG_DIR, "\n")

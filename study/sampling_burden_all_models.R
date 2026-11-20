################################################################
# sampling_burden_all_models.R
#
# Tissue-sampling burden N_q(Y, M) = min{ m : I_m / I_full >= q } for
# every model and state, under each model's linear readout (mean
# pooling, linear PCA, the arm whose slide logit is an exact average
# over tiles), with q = 0.90 and 0.95, and the robustness of the
# ordering of (model, state) pairs to the random tile draw.
#
# Inputs: pool_pass.R outputs (pooled.rds and submeans.npy, the mean
# of m random tiles per slide for m on a 16-point grid and R = 20
# draws). For each model and cohort: rebuild the linear PCA from the
# pooled means (must reproduce the saved Z_scaled), refit PathMI on
# the pca_mean arm (must reproduce the stored estimate), take the MAP
# of the selected model, and apply the affine map to the subsample
# means: the logit of a slide from m tiles is a0 + g(mean of the m
# tile embeddings). The plug-in difference-of-entropies estimate is
# formed over patients at every (m, draw); I_full uses all tiles.
#
# Robustness of the ordering: N_q is recomputed per draw; we report
# the per-draw range of N_q, Kendall's W of the ranking of the twelve
# (model, state) pairs across the 20 draws, and for each pair of
# pairs the fraction of draws in which their order agrees with the
# order of the mean curves.
#
# Usage (project root, R 4.4.3): Rscript TRIDENT/sampling_burden_all_models.R
# Output: TRIDENT/burden/{burden_curves,burden_Nq,burden_Nq_by_draw,burden_ordering}.csv,
#         figures/figS_sampling_burden_all.pdf (copied to manuscript/figures)
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(reticulate); library(ggplot2); library(dplyr); library(tidyr) })
reticulate::py_require("numpy"); np <- import("numpy")
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI  <- normalizePath(dirname(this_file), winslash = "/")
ROOT <- normalizePath(file.path(TRI, ".."), winslash = "/")
setwd(ROOT); source("PathMI_SUBHO_logistic.R"); source(file.path(TRI, "projection_functions.R"))
OUT <- file.path(TRI, "burden"); dir.create(OUT, showWarnings = FALSE)
FIG <- file.path(TRI, "figures"); MS_FIG <- file.path(ROOT, "manuscript", "figures")
LOG <- file.path(OUT, "sampling_burden_all_models.log")
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }
NMC <- 2000L; BURN <- 500L; THIN <- 10L; SEED <- 42L; N_PC <- 50L
POOL <- "M:/bhadury/mi_results/pool_pass"
MODELS <- c(uni_v2 = "UNI v2", virchow2 = "Virchow2", hoptimus1 = "H-optimus-1", conch_v15 = "CONCH v1.5")
TILE_PX <- c(uni_v2 = 256, virchow2 = 224, hoptimus1 = 224, conch_v15 = 512)   # pixels at 20x, about 0.5 um per pixel
Q <- c(0.90, 0.95)

map_logit <- function(X, y, v = 1, maxit = 50, tol = 1e-10) {
  p <- ncol(X); beta <- rep(0, p); P <- diag(1 / v, p)
  for (it in seq_len(maxit)) { eta <- as.vector(X %*% beta); pr <- plogis(eta); W <- pr * (1 - pr)
    g <- t(X) %*% (y - pr) - P %*% beta; H <- t(X) %*% (X * W) + P; d <- solve(H, g); beta <- beta + as.vector(d); if (max(abs(d)) < tol) break }
  list(beta = beta, eta = as.vector(X %*% beta))
}
auc <- function(s, y) { r <- rank(s); n1 <- sum(y); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1)) }
plugin <- function(eta, y, pbar) { p <- pmin(pmax(plogis(eta), 1e-12), 1 - 1e-12); mean(y * log(p) + (1 - y) * log(1 - p) - (y * log(pbar) + (1 - y) * log(1 - pbar))) }

curves <- list(); nq_draw <- list(); fits <- list()
for (co in c("crc", "prad")) for (mo in names(MODELS)) {
  t0 <- Sys.time()
  pp <- readRDS(file.path(POOL, sprintf("%s_%s", co, mo), "pooled.rds"))
  rds <- file.path(TRI, if (co == "crc") sprintf("crc_kpca_subho_%s/kpca_subho_results_crc_%s.rds", mo, mo) else sprintf("TCGA_PRAD_kpca_subho_%s/kpca_subho_results_prad_%s.rds", mo, mo))
  arm <- readRDS(rds)$pca_mean
  Zraw <- pp$pooled$mean; n <- nrow(Zraw); stopifnot(identical(rownames(Zraw), rownames(arm$Z_scaled)))
  pca <- prcomp(Zraw, center = TRUE, scale. = TRUE)
  Zs <- scale(pca$x[, 1:N_PC], center = TRUE, scale = TRUE)
  sgn <- sign(colSums(Zs * arm$Z_scaled)); sgn[sgn == 0] <- 1
  Zs <- sweep(Zs, 2, sgn, `*`); chk <- max(abs(Zs - arm$Z_scaled))
  R <- sweep(pca$rotation[, 1:N_PC], 2, sgn, `*`); m_k <- attr(scale(pca$x[, 1:N_PC]), "scaled:center") * sgn; s_k <- attr(scale(pca$x[, 1:N_PC]), "scaled:scale")
  c_j <- pca$center; d_j <- pca$scale
  logf("%s %s | n = %d | PCA rebuilt, max abs diff vs saved %.2e", toupper(co), mo, n, chk)
  if (chk > 1e-5) logf("  [WARN] projection mismatch")
  SUB <- np$load(file.path(POOL, sprintf("%s_%s", co, mo), "submeans.npy"), mmap_mode = "r")   # n x 16 x 20 x D
  Es_full <- sweep(sweep(Zraw, 2, c_j, `-`), 2, d_j, `/`)
  m_grid <- pp$m_grid; r_sub <- pp$r_sub
  for (tg in names(pp$y)) {
    y <- as.numeric(pp$y[[tg]][pp$roster]); stored <- arm$result$details[[tg]]$MI
    set.seed(SEED)
    obs <- subho_MI_logistic(Z = Zs, y = y, v_beta = 1, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
    logf("  [%s] refit MI %.6f | stored %.6f | diff %.1e | k* %d", tg, obs$MI, stored, obs$MI - stored, obs$n.selected)
    sel <- obs$selected.features; X <- cbind(1, Zs[, sel, drop = FALSE]); mp <- map_logit(X, y)
    b0 <- mp$beta[1]; bsel <- mp$beta[-1]
    v <- as.vector(R[, sel, drop = FALSE] %*% (bsel / s_k[sel])); a0 <- b0 - sum(bsel * m_k[sel] / s_k[sel])
    eta_full <- a0 + as.vector(Es_full %*% v); stopifnot(max(abs(eta_full - mp$eta)) < 1e-8)
    pbar <- mean(y); I_full <- plugin(eta_full, y, pbar); auc_full <- auc(eta_full, y)
    # subsample logits: ((SUB - c)/d) @ v, done in numpy on the float32 array
    v32 <- np$asarray(v / d_j, dtype = "float32"); off <- a0 - sum(c_j * v / d_j)
    ETA <- py_to_r(np$tensordot(SUB, v32, axes = 1L)) + off          # n x 16 x 20
    tn <- pp$tiles$n_tiles
    for (mi in seq_along(m_grid)) {
      I_r <- auc_r <- numeric(r_sub)
      for (r in seq_len(r_sub)) { e <- ETA[, mi, r]; I_r[r] <- plugin(e, y, pbar); auc_r[r] <- auc(e, y) }
      curves[[length(curves) + 1]] <- data.frame(cohort = toupper(co), model = MODELS[[mo]], model_code = mo, target = tg, m = m_grid[mi], tile_px = TILE_PX[[mo]],
        area_mm2 = m_grid[mi] * (TILE_PX[[mo]] * 0.5 / 1000)^2, I_mean = mean(I_r), I_lo = quantile(I_r, .05), I_hi = quantile(I_r, .95), I_sd = sd(I_r),
        auc_mean = mean(auc_r), I_full = I_full, auc_full = auc_full, rel = mean(I_r) / I_full, frac_slides_with_m = mean(tn >= m_grid[mi]))
      nq_draw[[length(nq_draw) + 1]] <- data.frame(cohort = toupper(co), model = MODELS[[mo]], target = tg, m = m_grid[mi], draw = seq_len(r_sub), rel = I_r / I_full)
    }
    fits[[length(fits) + 1]] <- data.frame(cohort = toupper(co), model = MODELS[[mo]], target = tg, n = n, prevalence = pbar, k_star = obs$n.selected,
      MI_pathmi = obs$MI, I_full_plugin = I_full, auc_full = auc_full, tiles_median = median(tn), tile_px = TILE_PX[[mo]])
  }
  rm(SUB); gc()
  logf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
curves <- do.call(rbind, curves); nq_draw <- do.call(rbind, nq_draw); fits <- do.call(rbind, fits)
curves$pair <- paste(curves$model, curves$target, sep = " / "); nq_draw$pair <- paste(nq_draw$model, nq_draw$target, sep = " / ")
write.csv(curves, file.path(OUT, "burden_curves.csv"), row.names = FALSE)
write.csv(fits, file.path(OUT, "burden_fits.csv"), row.names = FALSE)
write.csv(nq_draw, file.path(OUT, "burden_rel_by_draw.csv"), row.names = FALSE)   # relative estimate per (pair, m, draw)

# ---- N_q from the mean curve and per draw -----------------------------------------
Nq_of <- function(m, rel, q) { i <- which(rel >= q); if (length(i)) m[min(i)] else NA_integer_ }
nq <- curves %>% group_by(cohort, model, target, pair, tile_px) %>% summarise(N90 = Nq_of(m, rel, .90), N95 = Nq_of(m, rel, .95), I_full = I_full[1],
  rel_at_5 = rel[m == 5], rel_at_25 = rel[m == 25], rel_at_50 = rel[m == 50], rel_at_100 = rel[m == 100], .groups = "drop")
nqd <- nq_draw %>% group_by(cohort, model, target, pair, draw) %>% summarise(N90 = Nq_of(m, rel, .90), N95 = Nq_of(m, rel, .95), .groups = "drop")
nq_rng <- nqd %>% group_by(pair) %>% summarise(N90_med = median(N90), N90_min = min(N90), N90_max = max(N90), N95_med = median(N95), N95_min = min(N95), N95_max = max(N95), .groups = "drop")
nq <- left_join(nq, nq_rng, by = "pair") %>% mutate(A90_mm2 = N90 * (tile_px * 0.5 / 1000)^2, A95_mm2 = N95 * (tile_px * 0.5 / 1000)^2) %>% arrange(cohort, target, model)
write.csv(nq, file.path(OUT, "burden_Nq.csv"), row.names = FALSE); write.csv(nqd, file.path(OUT, "burden_Nq_by_draw.csv"), row.names = FALSE)

# ---- ordering robustness ------------------------------------------------------------
W_of <- function(d, col) { M <- d %>% select(pair, draw, all_of(col)) %>% pivot_wider(names_from = pair, values_from = all_of(col)) %>% select(-draw) %>% as.matrix()
  Rk <- t(apply(M, 1, rank)); k <- ncol(Rk); m <- nrow(Rk); S <- sum((colSums(Rk) - m * (k + 1) / 2)^2); 12 * S / (m^2 * (k^3 - k)) }
mean_order <- nq %>% arrange(N90) %>% pull(pair)
pairs <- t(combn(mean_order, 2))
agree <- apply(pairs, 1, function(pr) { a <- nqd %>% filter(pair == pr[1]) %>% arrange(draw) %>% pull(N90); b <- nqd %>% filter(pair == pr[2]) %>% arrange(draw) %>% pull(N90)
  c(frac_agree = mean(a <= b), frac_strict = mean(a < b), frac_tie = mean(a == b)) })
ordering <- data.frame(first = pairs[, 1], second = pairs[, 2], N90_first = nq$N90[match(pairs[, 1], nq$pair)], N90_second = nq$N90[match(pairs[, 2], nq$pair)], t(agree))
write.csv(ordering, file.path(OUT, "burden_ordering.csv"), row.names = FALSE)
W90 <- W_of(nqd, "N90"); W95 <- W_of(nqd, "N95")
# within-model state ordering (MSI vs fusions is across cohorts; ERG vs ETS within PRAD)
logf("=== N_q (tiles) from the mean curve, with per-draw median and range ===")
for (i in seq_len(nrow(nq))) with(nq[i, ], logf("%-28s N90 %4d [%d-%d] | N95 %4d [%d-%d] | I_full %.4f | rel at 5/25/50/100: %.2f/%.2f/%.2f/%.2f | area90 %.2f mm2",
  pair, N90, N90_min, N90_max, N95, N95_min, N95_max, I_full, rel_at_5, rel_at_25, rel_at_50, rel_at_100, A90_mm2))
logf("Kendall's W of the (model, state) ranking across 20 draws: N90 %.3f | N95 %.3f", W90, W95)
logf("pairs of (model, state) whose N90 order (mean curve) holds in every draw: %d of %d; in >= 90%% of draws: %d; ties in all draws: %d",
     sum(ordering$frac_agree == 1), nrow(ordering), sum(ordering$frac_agree >= .9), sum(ordering$frac_tie == 1))
st <- nq %>% group_by(target) %>% summarise(N90 = paste(range(N90), collapse = "-"), N95 = paste(range(N95), collapse = "-"))
for (i in seq_len(nrow(st))) logf("state %s: N90 across models %s | N95 %s", st$target[i], st$N90[i], st$N95[i])
writeLines(c(sprintf("W90,%.4f", W90), sprintf("W95,%.4f", W95)), file.path(OUT, "burden_kendallW.csv"))

# ---- figure --------------------------------------------------------------------------
pal <- c("UNI v2" = "#0072B2", "Virchow2" = "#D55E00", "H-optimus-1" = "#009E73", "CONCH v1.5" = "#CC79A7")
curves$target_lab <- factor(paste0(curves$cohort, "/", curves$target), c("CRC/MSI", "PRAD/ERG", "PRAD/ETS"))
p <- ggplot(curves, aes(x = m, y = rel, colour = model, fill = model)) + facet_wrap(~ target_lab, nrow = 1) +
  geom_hline(yintercept = c(.9, .95, 1), linetype = "22", colour = "grey60", linewidth = .3) +
  geom_ribbon(aes(ymin = I_lo / I_full, ymax = I_hi / I_full), alpha = .12, colour = NA) + geom_line(linewidth = .55) + geom_point(size = 1.1) +
  scale_x_log10(breaks = c(5, 25, 100, 500, 2500)) + scale_colour_manual(values = pal, name = NULL) + scale_fill_manual(values = pal, name = NULL) +
  coord_cartesian(ylim = c(0, 1.08)) + labs(x = "Tiles sampled per slide (log scale)", y = "Plug-in estimate relative to all tiles") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
PDF_DEV <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
ggsave(file.path(FIG, "figS_sampling_burden_all.pdf"), p, device = PDF_DEV, width = 183, height = 70, units = "mm")
ggsave(file.path(FIG, "figS_sampling_burden_all.png"), p, width = 183, height = 70, units = "mm", dpi = 300)
if (dir.exists(MS_FIG)) file.copy(file.path(FIG, "figS_sampling_burden_all.pdf"), MS_FIG, overwrite = TRUE)
logf("figure written")

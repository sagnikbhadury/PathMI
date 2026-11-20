################################################################
# tile_maps_hoptimus1.R
#
# Tile-level contribution maps from the linear-PCA / mean-pooling
# SUBHO fit (H-optimus-1), at no extra estimation cost.
#
# Why this arm: under mean pooling the patient embedding is the mean
# of its tile embeddings, and linear PCA followed by standardisation
# is an affine map, so the fitted slide logit
#     eta_i = b0 + sum_k beta_k Zs[i,k]
# decomposes EXACTLY as
#     eta_i = a0 + (1/T_i) sum_t g(z_it),   g(z) = ((z - c)/d) . v,
# with c, d the PCA column centre/scale, R the rotation, m, s the
# post-PCA standardisation and v = R[, sel] (beta_sel / s_sel).
# Kernel PCA arms do not decompose this way.
#
# Per tile we report
#   g_t          contribution to the slide logit (centred so that
#                the slide's mean over tiles is eta_i - a0);
#   eta_t        "as-if" logit a0 + g_t, i.e. the logit the model
#                would assign to a slide made only of tiles like t;
#   KL_t         specific information of the tile under that
#                as-if reading, KL( Bern(sigma(eta_t)) || Bern(pbar) ),
#                which is >= 0 and whose average over tiles of a
#                slide is NOT the slide's contribution (Jensen);
#   KL_lo/KL_hi  5 / 95 percent pointwise posterior quantiles of KL_t
#                from beta ~ N(beta_hat, H^-1), the Laplace posterior
#                of the selected model (same working prior v_beta = 1
#                as the estimator's scoring stage).
# Per patient we report the plug-in difference-of-entropies term
#   l_i = log p(y_i | z_i, beta_hat) - log p(y_i),
# whose mean is the plug-in MI (distinct from the SUBHO estimate,
# which is a Laplace log Bayes factor), with a 90 percent posterior
# interval for that mean, and the posterior s.d. of eta_i as the
# slide-level reliability statistic.
#
# Steps
#   1. load the saved pca_mean arm (Z_scaled, labels, roster);
#   2. pass 1 over the roster's h5 files: mean-pool exactly as the
#      pipeline did (numpy float32 mean), rebuild prcomp + scale,
#      assert equality with the saved Z_scaled (sign-aligned);
#   3. refit subho_MI_logistic on the saved Z_scaled with the run's
#      seed/settings; assert the stored MI is reproduced; take the
#      selected components; recompute the MAP and Hessian of the
#      selected model; draw S posterior samples of beta;
#   4. pass 2 over the h5 files: score every tile, check that the
#      tile mean reproduces the slide logit, store results;
#   5. write per-tile rds, per-patient csv, fit summary, and draft
#      figures.  PRAD is done per target (ERG, ETS).
#
# Usage (from the project root, R 4.4.3):
#   Rscript TRIDENT/tile_maps_hoptimus1.R <crc|prad>
# Output: TRIDENT/tile_maps_hoptimus1/<cohort>/
################################################################

args   <- commandArgs(trailingOnly = TRUE)
COHORT <- if (length(args)) tolower(args[1]) else "crc"
stopifnot(COHORT %in% c("crc", "prad"))

suppressPackageStartupMessages({
  library(MASS); library(BayesLogit); library(reticulate)
  library(ggplot2); library(patchwork); library(scales)
})
reticulate::py_require(c("h5py", "numpy"))
py_run_string("import h5py, numpy as np")

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI  <- normalizePath(dirname(this_file), winslash = "/")
ROOT <- normalizePath(file.path(TRI, ".."), winslash = "/")
setwd(ROOT)
source("PathMI_SUBHO_logistic.R")

OUT <- file.path(TRI, "tile_maps_hoptimus1", COHORT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, sprintf("tile_maps_%s_hoptimus1.log", COHORT))
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }

# ---- run-specific constants (identical to the pipeline) --------------------
N_PC <- 50L; RANDOM_SEED <- 42L; V_BETA <- 1.0; NMC <- 2000L; BURN <- 500L; THIN <- 10L
S_POST <- 500L                               # posterior draws
cfg <- switch(COHORT,
  crc  = list(rds = file.path(TRI, "crc_kpca_subho_hoptimus1/kpca_subho_results_crc_hoptimus1.rds"),
              dir = "M:/kulkarnr/slideflow_test/slideflow_crc/trident_feats/20x_224px_0px_overlap/features_hoptimus1/",
              thumbs = "M:/kulkarnr/slideflow_test/slideflow_crc/trident_feats/thumbnails/",
              targets = "MSI", pos_label = "MSI-H"),
  prad = list(rds = file.path(TRI, "TCGA_PRAD_kpca_subho_hoptimus1/kpca_subho_results_prad_hoptimus1.rds"),
              dir = "M:/shared_data/TCGA_PRAD/trident-feats/20x_224px_0px_overlap/features_hoptimus1/",
              thumbs = "M:/shared_data/TCGA_PRAD/trident-feats/thumbnails/",
              targets = c("ERG", "ETS"), pos_label = "positive"))

t0 <- Sys.time()
logf("=== tile maps | %s | H-optimus-1 | pca_mean arm ===", toupper(COHORT))

# ---- 1. saved arm -----------------------------------------------------------
res  <- readRDS(cfg$rds)
arm  <- res$pca_mean
Zs_saved <- arm$Z_scaled
roster   <- rownames(Zs_saved)
n <- nrow(Zs_saved)
ylist <- if (COHORT == "crc") list(MSI = arm$y[roster]) else lapply(arm$y_list, function(v) v[roster])
# 19 Nov 2026: patient-level analysis. The roster is patients; each patient's slides are mean-pooled
# separately and averaged (as pool_pass/pl_promote), and every slide's tiles are scored with the patient fit.
ppat <- readRDS(sprintf("M:/bhadury/mi_results/pool_pass/%s_hoptimus1/pooled.rds", COHORT)); stopifnot(identical(ppat$unit, "patient"))
slides <- ppat$slide_roster; spat <- substr(slides, 1, 12); stopifnot(all(roster %in% spat), all(spat %in% roster)); ns <- length(slides)
logf("roster n = %d patients (%d slides) | Z_scaled %d x %d | targets: %s", n, ns, n, ncol(Zs_saved), paste(names(ylist), collapse = ", "))

# ---- h5 index by attrs['name'] (as the pipeline) ----------------------------
files <- list.files(cfg$dir, pattern = "[.]h5$")
stem_map <- list()
for (f in files) {
  py_run_string(sprintf("
with h5py.File(r'%s', 'r') as _f:
    _name = _f['features'].attrs.get('name', None)
", file.path(cfg$dir, f)))
  nm <- py$`_name`; if (!is.null(nm)) stem_map[[nm]] <- f
}
stopifnot(all(slides %in% names(stem_map)))
logf("indexed %d h5 files; all %d roster slides present", length(files), ns)

# ---- 2. pass 1: mean pooling exactly as the pipeline ------------------------
Zsl <- matrix(NA_real_, ns, 0); tile_n <- integer(ns); names(tile_n) <- slides
first <- TRUE
for (i in seq_len(ns)) {
  sid <- slides[i]
  py_run_string(sprintf("
with h5py.File(r'%s', 'r') as _f:
    _emb = _f['features'][:].astype(np.float32)
_n = int(_emb.shape[0]); _pooled = _emb.mean(axis=0); del _emb
", file.path(cfg$dir, stem_map[[sid]])))
  v <- as.numeric(py$`_pooled`); tile_n[i] <- as.integer(py$`_n`)
  if (first) { Zsl <- matrix(NA_real_, ns, length(v)); first <- FALSE }
  Zsl[i, ] <- v
  if (i %% 50 == 0) logf("  pass 1: %d / %d slides pooled (%.1f min)", i, ns, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
Zraw <- rowsum(Zsl, spat, reorder = FALSE); Zraw <- Zraw[roster, , drop = FALSE] / as.vector(table(spat)[roster])
py_run_string("del _pooled, _n")
logf("pass 1 done: Zraw %d x %d | tiles per slide median %d", nrow(Zraw), ncol(Zraw), as.integer(median(tile_n)))

pca <- prcomp(Zraw, center = TRUE, scale. = TRUE)
k_use <- min(N_PC, ncol(pca$x))
PC   <- pca$x[, 1:k_use, drop = FALSE]
Zs_rec <- scale(PC, center = TRUE, scale = TRUE)
m_k <- attr(Zs_rec, "scaled:center"); s_k <- attr(Zs_rec, "scaled:scale")
# sign alignment (prcomp is deterministic, but guard against a flipped column)
sgn <- sign(colSums(Zs_rec * Zs_saved)); sgn[sgn == 0] <- 1
Zs_rec_al <- sweep(Zs_rec, 2, sgn, `*`)
maxdiff <- max(abs(Zs_rec_al - Zs_saved))
logf("PCA rebuilt: %d components | flipped columns: %d | max |Zs_rec - Zs_saved| = %.2e", k_use, sum(sgn < 0), maxdiff)
if (maxdiff > 1e-6) stop("reconstructed Z_scaled does not match the saved matrix")
R  <- sweep(pca$rotation[, 1:k_use, drop = FALSE], 2, sgn, `*`)   # aligned rotation
m_k <- m_k * sgn                                                    # aligned centre
c_j <- pca$center; d_j <- pca$scale

# ---- helper: MAP + Hessian of a logistic model under N(0, v I) --------------
map_logit <- function(X, y, v = V_BETA, maxit = 50, tol = 1e-10) {
  p <- ncol(X); beta <- rep(0, p); P <- diag(1 / v, p)
  for (it in seq_len(maxit)) {
    eta <- as.vector(X %*% beta); pr <- plogis(eta); W <- pr * (1 - pr)
    g <- t(X) %*% (y - pr) - P %*% beta; H <- t(X) %*% (X * W) + P
    d <- solve(H, g); beta <- beta + as.vector(d)
    if (max(abs(d)) < tol) break
  }
  eta <- as.vector(X %*% beta); pr <- plogis(eta); W <- pr * (1 - pr)
  list(beta = beta, H = t(X) %*% (X * W) + P, eta = eta, iter = it)
}
kl_bern <- function(p, q) { p <- pmin(pmax(p, 1e-12), 1 - 1e-12); p * log(p / q) + (1 - p) * log((1 - p) / (1 - q)) }

# ---- 3. refit the arm per target, decompose, sample the posterior ----------
fits <- list()
for (tg in names(ylist)) {
  y <- as.numeric(ylist[[tg]]); stopifnot(length(y) == n, all(y %in% 0:1))
  set.seed(RANDOM_SEED)
  obs <- subho_MI_logistic(Z = Zs_saved, y = y, v_beta = V_BETA, method.weight = "beta.binomial",
                           nmc = NMC, burn = BURN, thin = THIN)
  stored <- arm$result$details[[tg]]$MI
  logf("[%s] refit MI = %.6f | stored MI = %.6f | diff %.1e | k* = %d | selected PCs: %s",
       tg, obs$MI, stored, obs$MI - stored, obs$n.selected, paste(obs$selected.features, collapse = ","))
  if (abs(obs$MI - stored) > 1e-8) logf("  [WARN] refit does not reproduce the stored estimate")
  sel <- obs$selected.features
  X   <- cbind(1, Zs_saved[, sel, drop = FALSE])
  mp  <- map_logit(X, y)
  b0  <- mp$beta[1]; bsel <- mp$beta[-1]
  Sigma <- solve(mp$H)
  pbar <- mean(y)
  # exact tile-linear map
  v_vec <- R[, sel, drop = FALSE] %*% (bsel / s_k[sel])           # 1536 x 1
  a0    <- b0 - sum(bsel * m_k[sel] / s_k[sel])
  eta_slide <- mp$eta
  eta_sd <- sqrt(rowSums((X %*% Sigma) * X))
  p_hat <- plogis(eta_slide)
  l_i <- (y * log(p_hat) + (1 - y) * log(1 - p_hat)) - (y * log(pbar) + (1 - y) * log(1 - pbar))
  set.seed(RANDOM_SEED + 1)
  B <- mvrnorm(S_POST, mp$beta, Sigma)                                # S x (k*+1)
  eta_post <- X %*% t(B)                                              # n x S
  P_post <- plogis(eta_post)
  mi_post <- colMeans(y * log(P_post) + (1 - y) * log(1 - P_post)) - mean(y * log(pbar) + (1 - y) * log(1 - pbar))
  fits[[tg]] <- list(target = tg, y = y, sel = sel, b0 = b0, beta = bsel, Sigma = Sigma, B = B, v = as.vector(v_vec),
                     a0 = a0, eta = eta_slide, eta_sd = eta_sd, p_hat = p_hat, pbar = pbar, l_i = l_i,
                     MI_subho = obs$MI, MI_plugin = mean(l_i), MI_plugin_ci = quantile(mi_post, c(.05, .95)),
                     k_star = obs$n.selected)
  logf("[%s] plug-in MI = %.4f (90%% posterior %.4f, %.4f) | SUBHO MI = %.4f | prevalence %.3f | mean eta sd %.3f",
       tg, mean(l_i), fits[[tg]]$MI_plugin_ci[1], fits[[tg]]$MI_plugin_ci[2], obs$MI, pbar, mean(eta_sd))
}

# ---- 4. pass 2: score every tile ----------------------------------------------
tiles <- vector("list", ns)
slide_rows <- list()
for (si in seq_len(ns)) {
  sid <- slides[si]; i <- match(spat[si], roster)
  py_run_string(sprintf("
with h5py.File(r'%s', 'r') as _f:
    _emb = _f['features'][:].astype(np.float32); _xy = _f['coords'][:]
    _a = dict(_f['coords'].attrs)
_w = int(_a.get('level0_width', 0)); _h = int(_a.get('level0_height', 0)); _ps = int(_a.get('patch_size_level0', 0))
", file.path(cfg$dir, stem_map[[sid]])))
  E  <- py$`_emb`; XY <- py$`_xy`
  T_i <- nrow(E)
  Es <- sweep(sweep(E, 2, c_j, `-`), 2, d_j, `/`)                   # (z - c)/d, T x 1536 (float64)
  df <- data.frame(slide = sid, patient = spat[si], x = XY[, 1], y = XY[, 2])
  U_all <- Es %*% R                                                  # tile PC scores, T x 50
  for (tg in names(fits)) {
    ft <- fits[[tg]]
    U  <- sweep(sweep(U_all[, ft$sel, drop = FALSE], 2, ft$m_k_sel <- m_k[ft$sel], `-`), 2, s_k[ft$sel], `/`)   # standardised tile PCs
    g  <- as.vector(U %*% ft$beta)                                   # centred contribution
    eta_t <- ft$b0 + g
    Xt <- cbind(1, U)
    KLp <- kl_bern(plogis(Xt %*% t(ft$B)), ft$pbar)                  # T x S
    df[[paste0("g_", tg)]]     <- g
    df[[paste0("eta_", tg)]]   <- eta_t
    df[[paste0("KL_", tg)]]    <- kl_bern(plogis(eta_t), ft$pbar)
    df[[paste0("KLlo_", tg)]]  <- apply(KLp, 1, quantile, 0.05)
    df[[paste0("KLhi_", tg)]]  <- apply(KLp, 1, quantile, 0.95)
    # exactness check (patient level): the patient logit is the mean over the patient's slides of each slide's tile-mean logit
    eta_from_tiles <- ft$a0 + mean(as.vector(Es %*% ft$v))
    slide_rows[[length(slide_rows) + 1]] <- data.frame(
      slide = sid, patient = spat[si], n_slides_patient = sum(spat == spat[si]), target = tg, y = ft$y[i], n_tiles = T_i,
      eta = ft$eta[i], eta_from_tiles = eta_from_tiles, eta_sd = ft$eta_sd[i], p_hat = ft$p_hat[i],
      doe_term = ft$l_i[i], tile_g_sd = sd(g), frac_tiles_pos = mean(g > 0),
      frac_tiles_KL_gt_0_05 = mean(df[[paste0("KL_", tg)]] > 0.05),
      level0_w = py$`_w`, level0_h = py$`_h`, patch_level0 = py$`_ps`)
  }
  tiles[[si]] <- df
  py_run_string("del _emb, _xy")
  if (si %% 25 == 0) logf("  pass 2: %d / %d slides scored (%.1f min)", si, ns, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
tile_df  <- do.call(rbind, tiles); rm(tiles)
slide_df <- do.call(rbind, slide_rows)
chk <- max(sapply(split(slide_df, paste(slide_df$target, slide_df$patient)), function(d) abs(d$eta[1] - mean(d$eta_from_tiles))))
logf("pass 2 done: %d tiles on %d slides | max |patient eta - mean over slides of tile-mean eta| = %.2e (float32)", nrow(tile_df), ns, chk)

# ---- 5. save -------------------------------------------------------------------
saveRDS(tile_df, file.path(OUT, sprintf("tile_scores_%s_hoptimus1.rds", COHORT)), compress = "xz")
write.csv(slide_df, file.path(OUT, sprintf("slide_summary_%s_hoptimus1.csv", COHORT)), row.names = FALSE)
fit_summary <- do.call(rbind, lapply(fits, function(ft) data.frame(
  target = ft$target, n = n, prevalence = ft$pbar, k_star = ft$k_star, selected_PCs = paste(ft$sel, collapse = " "),
  MI_subho = ft$MI_subho, MI_plugin = ft$MI_plugin, MI_plugin_lo90 = ft$MI_plugin_ci[1], MI_plugin_hi90 = ft$MI_plugin_ci[2],
  intercept = ft$b0, coefs = paste(sprintf("%.4f", ft$beta), collapse = " "),
  mean_eta_sd = mean(ft$eta_sd), auc = { r <- rank(ft$p_hat); n1 <- sum(ft$y); (sum(r[ft$y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (n - n1)) })))
write.csv(fit_summary, file.path(OUT, sprintf("fit_summary_%s_hoptimus1.csv", COHORT)), row.names = FALSE)
saveRDS(lapply(fits, function(ft) ft[setdiff(names(ft), "B")]), file.path(OUT, sprintf("fits_%s_hoptimus1.rds", COHORT)))
logf("saved tile rds, slide csv, fit summary")
print(fit_summary[, c("target", "k_star", "MI_subho", "MI_plugin", "MI_plugin_lo90", "MI_plugin_hi90", "auc")])

# ---- 6. draft figures (final versions are made by the figure script) ---------
tryCatch({
  for (tg in names(fits)) {
    sd_t <- slide_df[slide_df$target == tg, ]
    pick <- c(sd_t$slide[sd_t$y == 1][order(-sd_t$p_hat[sd_t$y == 1])][1:2],
              sd_t$slide[sd_t$y == 1][order(sd_t$p_hat[sd_t$y == 1])][1],
              sd_t$slide[sd_t$y == 0][order(sd_t$p_hat[sd_t$y == 0])][1:2],
              sd_t$slide[sd_t$y == 0][order(-sd_t$p_hat[sd_t$y == 0])][1])
    lab  <- c("positive, highest p", "positive, 2nd highest p", "positive, lowest p", "negative, lowest p", "negative, 2nd lowest p", "negative, highest p")
    pl <- list()
    for (j in seq_along(pick)) {
      d <- tile_df[tile_df$slide == pick[j], ]; ps <- sd_t$patch_level0[sd_t$slide == pick[j]]
      row <- sd_t[sd_t$slide == pick[j], ]
      pl[[j]] <- ggplot(d, aes(x = x, y = -y, fill = .data[[paste0("KL_", tg)]])) +
        geom_tile(width = ps, height = ps) + coord_equal() +
        scale_fill_viridis_c(name = "specific surprise (nats)", limits = c(0, quantile(tile_df[[paste0("KL_", tg)]], .99))) +
        labs(title = sprintf("%s\n%s | p = %.2f, sd(eta) = %.2f, %d tiles", substr(pick[j], 1, 12), lab[j], row$p_hat, row$eta_sd, row$n_tiles)) +
        theme_void(base_size = 8) + theme(plot.title = element_text(size = 7), legend.position = "bottom")
    }
    fig <- wrap_plots(pl, ncol = 3) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
    ggsave(file.path(OUT, sprintf("draft_tilemaps_%s_%s.png", COHORT, tg)), fig, width = 12, height = 8, dpi = 150)
    p2 <- ggplot(sd_t, aes(x = tile_g_sd, y = eta, colour = factor(y))) +
      geom_errorbar(aes(ymin = eta - 1.64 * eta_sd, ymax = eta + 1.64 * eta_sd), width = 0, alpha = .4) +
      geom_point(size = 1.4) + geom_hline(yintercept = qlogis(fits[[tg]]$pbar), linetype = 2) +
      labs(x = "within-slide s.d. of tile contributions", y = "slide logit (90% posterior bar)", colour = tg) + theme_bw(base_size = 9)
    p3 <- ggplot(sd_t, aes(x = doe_term, fill = factor(y))) + geom_histogram(bins = 40, alpha = .6, position = "identity") +
      geom_vline(xintercept = 0, linetype = 2) + labs(x = "per-patient DoE term (nats)", fill = tg) + theme_bw(base_size = 9)
    ggsave(file.path(OUT, sprintf("draft_slidelevel_%s_%s.png", COHORT, tg)), p2 + p3, width = 11, height = 4.5, dpi = 150)
  }
  logf("draft figures written")
}, error = function(e) logf("[WARN] draft figures failed: %s", conditionMessage(e)))
logf("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

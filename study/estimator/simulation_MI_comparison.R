################################################################
# Simulation: Comparing DOE-type MI Estimators
#
# Estimators:
#   1. Horseshoe SUBHO logistic  — subho_MI_logistic()
#   2. LASSO logistic CV         — run_lasso_cv() + compute_doe_mi()
#   3. OLS (p < n only)          — Gaussian LPM, no selection, baseline
#   4. DOE linear                — DV lower bound, joint vs marginal, baseline
#
# Sources required (load before running):
#   source("PathMI_SUBHO_logistic.R")   # subho_MI_logistic, lpy.x.logistic
#   source("PathMI_LASSO.R")            # run_lasso_cv, compute_doe_mi
#   source("PathMI_OLS.R")              # compute_ols_mi
#   source("PathMI_DOE_linear.R")       # doe_linear_single
#
# Design (fully crossed, skip s >= p):
#   n:             100, 200, 500
#   p:             50, 200, 500
#   s (signals):   5, 20, 50
#   beta strength: 0.3 (weak), 0.7 (moderate), 1.5 (strong)
#   rho (feature corr): 0.3 (compound symmetric)
#
# Ground truth MI:
#   H(Y) - H(Y|Z) via large Monte Carlo from known logistic DGP.
#   Consistent with the logistic likelihood used by both estimators.
#
# Performance metrics:
#   bias  = MI_hat - MI_true
#   RMSE  = sqrt(mean(bias^2)) across replicates
#   ratio = MI_hat / MI_true   (1 = perfect)
#   TPR   = true positive rate of selected features (selection methods only)
#   FPR   = false positive rate of selected features (selection methods only)
#
# Dependencies: MASS, glmnet, BayesLogit, matrixcalc
################################################################


suppressPackageStartupMessages({
  library(MASS)
  library(glmnet)
})

################################################################
# DATA GENERATING MECHANISM
#
# Z_i ~ N(0, Sigma),  Sigma_jk = rho (j != k), 1 (j == k)
# logit P(Y=1|Z_i) = Z_i[1:s] %*% beta_true
# beta_true[1:s]: alternating signs and magnitudes for realism
################################################################

generate_data <- function(n, p, s, beta_strength, rho = 0.3, seed = 1) {
  set.seed(seed)

  Sigma        <- matrix(rho, p, p)
  diag(Sigma)  <- 1
  Z            <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

  # alternating sign pattern: 1, -1, 0.5, -0.5, 0.8, ... repeated
  sign_mag     <- c(1, -1, 0.5, -0.5, 0.8)
  beta_true    <- rep(0, p)
  beta_true[seq_len(s)] <- beta_strength * sign_mag[((seq_len(s) - 1) %% 5) + 1]

  eta  <- as.vector(Z %*% beta_true)
  p_i  <- plogis(eta)
  y    <- rbinom(n, 1, p_i)

  list(Z = Z, y = y, beta_true = beta_true, true_signals = seq_len(s))
}

################################################################
# TRUE MI VIA MONTE CARLO
#
# I(Z;Y) = H(Y) - H(Y|Z)
# H(Y)   = binary entropy of marginal P(Y=1)
# H(Y|Z) = -E[p*log(p) + (1-p)*log(1-p)],  p = plogis(Z beta)
# N = 50000 draws for low MC variance
################################################################

compute_true_mi <- function(beta_true, p_dim, rho = 0.3, N = 50000, seed = 999) {
  set.seed(seed)
  Sigma <- matrix(rho, p_dim, p_dim); diag(Sigma) <- 1
  Z_mc  <- mvrnorm(N, mu = rep(0, p_dim), Sigma = Sigma)
  p_hat <- plogis(as.vector(Z_mc %*% beta_true))
  p_hat <- pmax(pmin(p_hat, 1 - 1e-10), 1e-10)

  H_Y_given_Z <- -mean(p_hat * log(p_hat) + (1 - p_hat) * log(1 - p_hat))
  p_marg      <- mean(p_hat)
  p_marg      <- pmax(pmin(p_marg, 1 - 1e-10), 1e-10)
  H_Y         <- -(p_marg * log(p_marg) + (1 - p_marg) * log(1 - p_marg))

  max(H_Y - H_Y_given_Z, 0)
}

################################################################
# ESTIMATOR WRAPPERS
# Standardized interface: each returns list(MI, selected, n_selected)
# selected = integer vector of feature indices (1-based, no intercept)
################################################################

## 1. Horseshoe SUBHO logistic
# subho_MI_logistic() returns: $MI, $n.selected, $selected.features
wrap_horseshoe <- function(Z, y, nmc = 2000, burn = 200, thin = 10) {
  tryCatch({
    res <- subho_MI_logistic(Z = Z, y = y,
                              v_beta        = 1,
                              method.weight = "beta.binomial",
                              nmc = nmc, burn = burn, thin = thin)
    list(MI         = res$MI,
         selected   = res$selected.features,
         n_selected = res$n.selected)
  }, error = function(e) {
    warning(sprintf("Horseshoe failed: %s", conditionMessage(e)))
    list(MI = NA, selected = integer(0), n_selected = NA)
  })
}

## 2. LASSO logistic CV
# run_lasso_cv() returns: $selected_features, $coef_lasso, $lambda_min, $cv_fit
# compute_doe_mi() returns: $MI, $ll_fitted, $ll_null, $n_selected
wrap_lasso <- function(Z, y, nfolds = 10, seed = 42) {
  tryCatch({
    lasso_res <- run_lasso_cv(Z, y, nfolds = nfolds, seed = seed)
    mi_res    <- compute_doe_mi(Z, y, lasso_res)
    list(MI         = mi_res$MI,
         selected   = lasso_res$selected_features,
         n_selected = mi_res$n_selected)
  }, error = function(e) {
    warning(sprintf("LASSO failed: %s", conditionMessage(e)))
    list(MI = NA, selected = integer(0), n_selected = NA)
  })
}

## 3. OLS baseline — only valid when p < n
# compute_ols_mi() returns: $MI, $R2, $sigma2_null, $sigma2_fit
wrap_ols <- function(Z, y) {
  if (ncol(Z) >= nrow(Z)) {
    return(list(MI = NA, selected = seq_len(ncol(Z)), n_selected = ncol(Z)))
  }
  tryCatch({
    res <- compute_ols_mi(Z, y)
    list(MI = res$MI, selected = seq_len(ncol(Z)), n_selected = ncol(Z))
  }, error = function(e) {
    list(MI = NA, selected = integer(0), n_selected = NA)
  })
}

## 4. DOE linear baseline — DV lower bound on joint vs marginal
# doe_linear_single() returns: $MI, $E_joint, $log_E_marginal
wrap_doe_linear <- function(Z, y, seed = 42) {
  tryCatch({
    res <- doe_linear_single(Z, y, seed = seed)
    list(MI = res$MI, selected = seq_len(ncol(Z)), n_selected = ncol(Z))
  }, error = function(e) {
    list(MI = NA, selected = integer(0), n_selected = NA)
  })
}

################################################################
# SELECTION METRICS (horseshoe and LASSO only)
################################################################

selection_metrics <- function(selected, true_signals, p) {
  noise <- setdiff(seq_len(p), true_signals)
  TP    <- length(intersect(selected, true_signals))
  FP    <- length(intersect(selected, noise))
  TPR   <- TP / length(true_signals)
  FPR   <- if (length(noise) > 0) FP / length(noise) else NA
  list(TPR = TPR, FPR = FPR)
}

################################################################
# SINGLE REPLICATE
################################################################

run_one_replicate <- function(n, p, s, beta_strength, rho = 0.3, seed = 1,
                               nmc = 2000, burn = 200, thin = 10,
                               nfolds = 10) {

  dat      <- generate_data(n, p, s, beta_strength, rho = rho, seed = seed)
  Z        <- dat$Z
  y        <- dat$y
  MI_true  <- compute_true_mi(dat$beta_true, p, rho = rho)
  true_sig <- dat$true_signals

  # scale once — consistent with both PathMI pipelines
  Z_sc <- scale(Z)

  res_hs  <- wrap_horseshoe(Z_sc, y, nmc = nmc, burn = burn, thin = thin)
  res_las <- wrap_lasso(Z_sc, y, nfolds = nfolds, seed = seed)
  res_ols <- wrap_ols(Z_sc, y)
  res_doe <- wrap_doe_linear(Z_sc, y, seed = seed)

  sm_hs  <- selection_metrics(res_hs$selected,  true_sig, p)
  sm_las <- selection_metrics(res_las$selected, true_sig, p)

  data.frame(
    n = n, p = p, s = s, beta = beta_strength, seed = seed,
    MI_true = MI_true,

    MI_horseshoe   = res_hs$MI,
    MI_lasso       = res_las$MI,
    MI_ols         = res_ols$MI,
    MI_doe         = res_doe$MI,

    nsel_horseshoe = res_hs$n_selected,
    nsel_lasso     = res_las$n_selected,

    bias_horseshoe = res_hs$MI  - MI_true,
    bias_lasso     = res_las$MI - MI_true,
    bias_ols       = res_ols$MI - MI_true,
    bias_doe       = res_doe$MI - MI_true,

    ratio_horseshoe = res_hs$MI  / MI_true,
    ratio_lasso     = res_las$MI / MI_true,
    ratio_ols       = res_ols$MI / MI_true,
    ratio_doe       = res_doe$MI / MI_true,

    TPR_horseshoe  = sm_hs$TPR,
    TPR_lasso      = sm_las$TPR,
    FPR_horseshoe  = sm_hs$FPR,
    FPR_lasso      = sm_las$FPR,

    stringsAsFactors = FALSE
  )
}

################################################################
# FULL SIMULATION
################################################################

run_simulation <- function(n_reps      = 100,
                            n_vec       = c(100, 200, 500),
                            p_vec       = c(50, 200, 500),
                            s_vec       = c(5, 20, 50),
                            beta_vec    = c(0.3, 0.7, 1.5),
                            rho         = 0.3,
                            nmc         = 2000,
                            burn        = 200,
                            thin        = 10,
                            nfolds      = 10,
                            output_file = "simulation_results.rds") {

  grid <- expand.grid(n = n_vec, p = p_vec, s = s_vec, beta = beta_vec,
                      stringsAsFactors = FALSE)
  grid <- grid[grid$s < grid$p, ]
  cat(sprintf("Scenarios: %d | Replicates: %d | Total runs: %d\n\n",
              nrow(grid), n_reps, nrow(grid) * n_reps))

  all_results <- list()

  for (i in seq_len(nrow(grid))) {
    cfg <- grid[i, ]
    cat(sprintf("[%d/%d] n=%d p=%d s=%d beta=%.1f\n",
                i, nrow(grid), cfg$n, cfg$p, cfg$s, cfg$beta))

    for (r in seq_len(n_reps)) {
      seed_r <- i * 10000 + r
      res <- tryCatch(
        run_one_replicate(
          n = cfg$n, p = cfg$p, s = cfg$s, beta_strength = cfg$beta,
          rho = rho, seed = seed_r,
          nmc = nmc, burn = burn, thin = thin, nfolds = nfolds
        ),
        error = function(e) {
          warning(sprintf("  Rep %d failed: %s", r, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(res)) all_results <- c(all_results, list(res))
    }
  }

  results_df <- do.call(rbind, all_results)
  saveRDS(results_df, output_file)
  cat(sprintf("\nDone. %d completed replicates -> %s\n", nrow(results_df), output_file))
  invisible(results_df)
}

################################################################
# SUMMARY TABLE
################################################################

summarize_results <- function(results_df) {

  estimators  <- c("horseshoe", "lasso", "ols", "doe")
  scenario_cols <- c("n", "p", "s", "beta")

  # unique scenarios
  scenarios <- unique(results_df[, scenario_cols])
  scenarios <- scenarios[order(scenarios$n, scenarios$p, scenarios$s, scenarios$beta), ]

  out_rows <- list()

  for (est in estimators) {
    bias_col  <- paste0("bias_",  est)
    ratio_col <- paste0("ratio_", est)

    for (i in seq_len(nrow(scenarios))) {
      sc  <- scenarios[i, ]
      idx <- results_df$n    == sc$n    &
             results_df$p    == sc$p    &
             results_df$s    == sc$s    &
             results_df$beta == sc$beta

      df  <- results_df[idx, ]
      bv  <- df[[bias_col]]
      rv  <- df[[ratio_col]]

      row <- data.frame(
        estimator  = est,
        n          = sc$n,
        p          = sc$p,
        s          = sc$s,
        beta       = sc$beta,
        n_reps     = sum(!is.na(bv)),
        mean_bias  = mean(bv, na.rm = TRUE),
        RMSE       = sqrt(mean(bv^2, na.rm = TRUE)),
        mean_ratio = mean(rv, na.rm = TRUE),
        sd_ratio   = sd(rv,   na.rm = TRUE),
        mean_TPR   = if (est %in% c("horseshoe","lasso")) mean(df[[paste0("TPR_",est)]], na.rm=TRUE) else NA_real_,
        mean_FPR   = if (est %in% c("horseshoe","lasso")) mean(df[[paste0("FPR_",est)]], na.rm=TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
      out_rows <- c(out_rows, list(row))
    }
  }

  do.call(rbind, out_rows)
}

################################################################
# PILOT — validate pipeline before committing to full run
# Use reduced MCMC to keep runtime short
################################################################

run_pilot <- function(n_reps = 3, nmc = 500, burn = 100, thin = 5) {
  run_simulation(
    n_reps   = n_reps,
    n_vec    = c(100, 200),
    p_vec    = c(50, 200),
    s_vec    = c(5, 20),
    beta_vec = c(0.3, 1, 1.5),
    nmc = nmc, burn = burn, thin = thin,
    output_file = "pilot_results.rds"
  )
}

################################################################
# PLOTS
#
# All plots accept the raw results_df from run_simulation().
# Call after summarize_results() or directly on results_df.
#
# Figures produced:
#   1. plot_bias_by_scenario   — mean bias across n/p/s/beta per estimator
#   2. plot_rmse_heatmap       — RMSE heatmap: estimator x (n,p) grid
#   3. plot_ratio_boxplot      — MI ratio distribution per estimator
#   4. plot_tpr_fpr            — TPR vs FPR for horseshoe and LASSO
#   5. plot_mi_scatter         — MI_hat vs MI_true scatter, all estimators
#   6. plot_bias_vs_dim        — bias as function of p, faceted by beta
#   7. plot_nselected          — distribution of n_selected for HS and LASSO
################################################################

EST_COLS <- c(
  horseshoe = "#C0392B",   # red
  lasso     = "#2980B9",   # blue
  ols       = "#27AE60",   # green
  doe       = "#8E44AD"    # purple
)
EST_LABELS <- c(
  horseshoe = "PathMI",
  lasso     = "LASSO-CV",
  ols       = "OLS",
  doe       = "DOE Linear"
)

################################################################
# FIG 1: Mean bias by n, faceted by beta strength
# One panel per beta level; lines for each estimator over n.
# Fixed p and s to a representative scenario.
################################################################

plot_bias_by_n <- function(results_df,
                            fixed_p   = 200,
                            fixed_s   = 5,
                            estimators = c("horseshoe", "lasso", "ols", "doe"),
                            save_pdf  = FALSE,
                            file      = "fig1_bias_by_n.pdf") {

  df <- results_df[results_df$p == fixed_p & results_df$s == fixed_s, ]
  if (nrow(df) == 0) stop("No rows match fixed_p and fixed_s. Adjust values.")

  beta_levels <- sort(unique(df$beta))
  n_levels    <- sort(unique(df$n))
  n_beta      <- length(beta_levels)

  if (save_pdf) pdf(file, width = 4 * n_beta, height = 4)
  par(mfrow = c(1, n_beta), mar = c(4, 4, 3, 1))

  for (b in beta_levels) {
    sub <- df[df$beta == b, ]

    means <- sapply(estimators, function(est) {
      tapply(sub[[paste0("bias_", est)]], sub$n, mean, na.rm = TRUE)
    })
    # means is n_levels x n_estimators

    ylim <- range(means, na.rm = TRUE)
    ylim <- ylim + diff(ylim) * c(-0.1, 0.1)

    plot(NA, xlim = range(n_levels), ylim = ylim,
         xlab = "n", ylab = "Mean bias (MI_hat - MI_true)",
         main = sprintf("p=%d, s=%d, beta=%.1f", fixed_p, fixed_s, b),
         xaxt = "n")
    axis(1, at = n_levels)
    abline(h = 0, lty = 2, col = "grey60")

    for (j in seq_along(estimators)) {
      est <- estimators[j]
      lines(n_levels, means[, j],
            col = EST_COLS[est], lwd = 2, type = "b", pch = 19)
    }
    legend("topright", bty = "n",
           legend = EST_LABELS[estimators],
           col    = EST_COLS[estimators],
           lwd = 2, pch = 19, cex = 0.75)
  }

  if (save_pdf) { dev.off(); cat(sprintf("Saved: %s\n", file)) }
  invisible(NULL)
}

################################################################
# FIG 2: RMSE heatmap — estimator x p, one panel per n
# Rows = estimators, columns = p values
################################################################

plot_rmse_heatmap <- function(results_df,
                               fixed_s    = 5,
                               fixed_beta = 0.7,
                               estimators = c("horseshoe", "lasso", "ols", "doe"),
                               save_pdf   = FALSE,
                               file       = "fig2_rmse_heatmap.pdf") {

  df <- results_df[results_df$s == fixed_s & results_df$beta == fixed_beta, ]
  n_levels <- sort(unique(df$n))
  p_levels <- sort(unique(df$p))

  if (save_pdf) pdf(file, width = 4 * length(n_levels), height = 4)
  par(mfrow = c(1, length(n_levels)), mar = c(5, 8, 3, 2))

  for (nv in n_levels) {
    sub <- df[df$n == nv, ]

    mat <- matrix(NA, nrow = length(estimators), ncol = length(p_levels),
                  dimnames = list(EST_LABELS[estimators], p_levels))
    for (j in seq_along(estimators)) {
      est <- estimators[j]
      for (k in seq_along(p_levels)) {
        vals <- sub[sub$p == p_levels[k], paste0("bias_", est)]
        mat[j, k] <- sqrt(mean(vals^2, na.rm = TRUE))
      }
    }

    max_val <- max(mat, na.rm = TRUE)
    image(t(mat[nrow(mat):1, ]),
          col  = hcl.colors(20, "YlOrRd", rev = FALSE),
          xaxt = "n", yaxt = "n",
          main = sprintf("RMSE | n=%d, s=%d, beta=%.1f", nv, fixed_s, fixed_beta),
          xlab = "p", ylab = "")
    axis(1, at = seq(0, 1, length.out = length(p_levels)), labels = p_levels)
    axis(2, at = seq(0, 1, length.out = length(estimators)),
         labels = rev(EST_LABELS[estimators]), las = 1)

    # cell labels
    for (j in seq_along(estimators)) {
      for (k in seq_along(p_levels)) {
        text(x = (k - 1) / (length(p_levels) - 1),
             y = 1 - (j - 1) / (length(estimators) - 1),
             labels = sprintf("%.3f", mat[j, k]), cex = 0.8)
      }
    }
  }

  if (save_pdf) { dev.off(); cat(sprintf("Saved: %s\n", file)) }
  invisible(NULL)
}

################################################################
# FIG 3: MI ratio boxplot — all estimators, faceted by p/n regime
# Shows distribution of MI_hat / MI_true across replicates.
# Reference line at 1 = perfect estimation.
################################################################

plot_ratio_boxplot <- function(results_df,
                                fixed_s    = 5,
                                fixed_beta = 0.7,
                                estimators = c("horseshoe", "lasso", "ols", "doe"),
                                save_pdf   = FALSE,
                                file       = "fig3_ratio_boxplot.pdf") {

  df <- results_df[results_df$s == fixed_s & results_df$beta == fixed_beta, ]
  df$regime <- paste0("n=", df$n, "\np=", df$p)
  regimes   <- unique(df$regime)

  if (save_pdf) pdf(file, width = max(8, 2 * length(regimes)), height = 5)
  par(mar = c(6, 4, 3, 1))

  # build list of ratio vectors: one per estimator x regime combination
  box_data  <- list()
  box_names <- c()
  box_cols  <- c()

  for (r in regimes) {
    for (est in estimators) {
      vals <- df[df$regime == r, paste0("ratio_", est)]
      vals <- vals[!is.na(vals) & is.finite(vals)]
      box_data  <- c(box_data,  list(vals))
      box_names <- c(box_names, sprintf("%s\n%s", EST_LABELS[est], r))
      box_cols  <- c(box_cols,  EST_COLS[est])
    }
  }

  boxplot(box_data,
          col    = box_cols,
          names  = box_names,
          las    = 2,
          cex.axis = 0.65,
          ylab   = "MI ratio (estimated / true)",
          main   = sprintf("MI ratio | s=%d, beta=%.1f", fixed_s, fixed_beta),
          outline = FALSE)
  abline(h = 1, lty = 2, col = "grey40", lwd = 1.5)

  if (save_pdf) { dev.off(); cat(sprintf("Saved: %s\n", file)) }
  invisible(NULL)
}

################################################################
# FIG 4: TPR vs FPR — horseshoe and LASSO only
# Each point = one (n, p, s, beta) scenario mean.
# Shows selection accuracy tradeoff.
################################################################

plot_tpr_fpr <- function(results_df,
                          save_pdf = FALSE,
                          file     = "fig4_tpr_fpr.pdf") {

  estimators <- c("horseshoe", "lasso")

  agg_rows <- list()
  scenarios <- unique(results_df[, c("n","p","s","beta")])
  for (i in seq_len(nrow(scenarios))) {
    sc  <- scenarios[i, ]
    idx <- results_df$n == sc$n & results_df$p == sc$p &
           results_df$s == sc$s & results_df$beta == sc$beta
    df  <- results_df[idx, ]
    agg_rows <- c(agg_rows, list(data.frame(
      n = sc$n, p = sc$p, s = sc$s, beta = sc$beta,
      TPR_horseshoe = mean(df$TPR_horseshoe, na.rm = TRUE),
      FPR_horseshoe = mean(df$FPR_horseshoe, na.rm = TRUE),
      TPR_lasso     = mean(df$TPR_lasso,     na.rm = TRUE),
      FPR_lasso     = mean(df$FPR_lasso,     na.rm = TRUE)
    )))
  }
  agg <- do.call(rbind, agg_rows)

  if (save_pdf) pdf(file, width = 6, height = 5)
  par(mar = c(4, 4, 3, 1))

  plot(NA, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Mean FPR", ylab = "Mean TPR",
       main = "TPR vs FPR: Horseshoe vs LASSO")
  abline(a = 0, b = 1, lty = 3, col = "grey70")

  pch_map <- c(`100` = 15, `200` = 17, `500` = 19)   # shape by n

  for (est in estimators) {
    for (nv in unique(agg$n)) {
      sub <- agg[agg$n == nv, ]
      points(sub[[paste0("FPR_", est)]],
             sub[[paste0("TPR_", est)]],
             col = EST_COLS[est],
             pch = pch_map[as.character(nv)],
             cex = 1.2)
    }
  }

  legend("bottomright", bty = "n",
         legend = c(EST_LABELS[estimators],
                    paste0("n=", names(pch_map))),
         col    = c(EST_COLS[estimators], rep("grey40", 3)),
         pch    = c(19, 19, unname(pch_map)),
         cex    = 0.8)

  if (save_pdf) { dev.off(); cat(sprintf("Saved: %s\n", file)) }
  invisible(NULL)
}

################################################################
# FIG 5: MI_hat vs MI_true scatter — all estimators
# One panel per estimator. Identity line = perfect calibration.
################################################################

plot_mi_scatter <- function(results_df,
                             estimators = c("horseshoe", "lasso", "ols", "doe"),
                             save_pdf   = FALSE,
                             file       = "fig5_mi_scatter.pdf") {

  if (save_pdf) pdf(file, width = 4 * length(estimators), height = 4)
  par(mfrow = c(1, length(estimators)), mar = c(4, 4, 3, 1))

  lims <- range(c(results_df$MI_true,
                  unlist(lapply(estimators, function(e) results_df[[paste0("MI_", e)]]))),
                na.rm = TRUE)
  lims[1] <- 0

  for (est in estimators) {
    mi_hat  <- results_df[[paste0("MI_", est)]]
    mi_true <- results_df$MI_true
    ok      <- !is.na(mi_hat) & is.finite(mi_hat)

    plot(mi_true[ok], mi_hat[ok],
         col  = adjustcolor(EST_COLS[est], alpha.f = 0.35),
         pch  = 19, cex = 0.5,
         xlim = lims, ylim = lims,
         xlab = "MI true", ylab = "MI estimated",
         main = EST_LABELS[est])
    abline(a = 0, b = 1, lty = 2, col = "grey30", lwd = 1.5)

    # add RMSE annotation
    rmse <- sqrt(mean((mi_hat[ok] - mi_true[ok])^2))
    legend("topleft", bty = "n",
           legend = sprintf("RMSE = %.4f", rmse), cex = 0.85)
  }

  if (save_pdf) { dev.off(); cat(sprintf("Saved: %s\n", file)) }
  invisible(NULL)
}

################################################################
# FIG 6: Bias vs p — shows curse of dimensionality
# Lines per estimator; panels per beta level.
# Aggregated over n and s at fixed s = 5 (sparse regime).
################################################################

plot_bias_vs_p <- function(results_df,
                            fixed_s    = 5,
                            fixed_n    = 200,
                            estimators = c("horseshoe", "lasso", "ols", "doe"),
                            save_pdf   = FALSE,
                            file       = "fig6_bias_vs_p.pdf") {

  df <- results_df[results_df$s == fixed_s & results_df$n == fixed_n, ]
  if (nrow(df) == 0) stop("No rows match fixed_s and fixed_n.")

  beta_levels <- sort(unique(df$beta))
  p_levels    <- sort(unique(df$p))

  if (save_pdf) pdf(file, width = 4 * length(beta_levels), height = 4)
  par(mfrow = c(1, length(beta_levels)), mar = c(4, 4, 3, 1))

  for (b in beta_levels) {
    sub <- df[df$beta == b, ]

    means <- sapply(estimators, function(est) {
      tapply(sub[[paste0("bias_", est)]], sub$p, mean, na.rm = TRUE)
    })

    ylim <- range(means, na.rm = TRUE)
    ylim <- ylim + diff(ylim) * c(-0.1, 0.15)

    plot(NA, xlim = range(p_levels), ylim = ylim,
         xlab = "p (dimensionality)", ylab = "Mean bias",
         main = sprintf("n=%d, s=%d, beta=%.1f", fixed_n, fixed_s, b),
         xaxt = "n")
    axis(1, at = p_levels)
    abline(h = 0, lty = 2, col = "grey60")

    for (j in seq_along(estimators)) {
      est <- estimators[j]
      lines(p_levels, means[, j],
            col = EST_COLS[est], lwd = 2, type = "b", pch = 19)
    }
    legend("topright", bty = "n",
           legend = EST_LABELS[estimators],
           col    = EST_COLS[estimators],
           lwd = 2, pch = 19, cex = 0.75)
  }

  if (save_pdf) { dev.off(); cat(sprintf("Saved: %s\n", file)) }
  invisible(NULL)
}

################################################################
# FIG 7: Distribution of n_selected — horseshoe vs LASSO
# Histograms per (n, p) combination. True s marked as vertical line.
################################################################

plot_nselected <- function(results_df,
                            fixed_beta = 0.7,
                            save_pdf   = FALSE,
                            file       = "fig7_nselected.pdf") {

  df      <- results_df[results_df$beta == fixed_beta, ]
  n_vals  <- sort(unique(df$n))
  p_vals  <- sort(unique(df$p))
  n_panel <- length(n_vals) * length(p_vals)

  if (save_pdf) pdf(file, width = 4 * length(p_vals), height = 4 * length(n_vals))
  par(mfrow = c(length(n_vals), length(p_vals)), mar = c(4, 3, 3, 1))

  for (nv in n_vals) {
    for (pv in p_vals) {
      sub     <- df[df$n == nv & df$p == pv, ]
      true_s  <- unique(sub$s)   # may be multiple — take all

      xlim <- range(c(sub$nsel_horseshoe, sub$nsel_lasso), na.rm = TRUE)

      hist(sub$nsel_horseshoe, breaks = 20,
           col  = adjustcolor(EST_COLS["horseshoe"], 0.5),
           border = NA,
           xlim = xlim,
           main = sprintf("n=%d, p=%d", nv, pv),
           xlab = "# features selected", ylab = "")
      hist(sub$nsel_lasso, breaks = 20,
           col  = adjustcolor(EST_COLS["lasso"], 0.5),
           border = NA, add = TRUE)

      # vertical lines for each true s value in this panel
      for (sv in sort(unique(sub$s))) {
        abline(v = sv, lty = 2, col = "grey30", lwd = 1.2)
      }
      legend("topright", bty = "n",
             legend = c("Horseshoe", "LASSO", "True s"),
             fill   = c(adjustcolor(EST_COLS["horseshoe"], 0.5),
                        adjustcolor(EST_COLS["lasso"], 0.5), NA),
             border = c(NA, NA, NA),
             lty    = c(NA, NA, 2),
             col    = c(NA, NA, "grey30"),
             cex    = 0.7, merge = TRUE)
    }
  }

  if (save_pdf) { dev.off(); cat(sprintf("Saved: %s\n", file)) }
  invisible(NULL)
}

################################################################
# CONVENIENCE: produce all figures at once
################################################################

plot_all <- function(results_df,
                     fixed_p    = 200,
                     fixed_s    = 5,
                     fixed_n    = 200,
                     fixed_beta = 0.7,
                     save_pdf   = TRUE,
                     prefix     = "fig") {

  plot_bias_by_n(results_df, fixed_p = fixed_p, fixed_s = fixed_s,
                 save_pdf = save_pdf, file = paste0(prefix, "1_bias_by_n.pdf"))

  plot_rmse_heatmap(results_df, fixed_s = fixed_s, fixed_beta = fixed_beta,
                    save_pdf = save_pdf, file = paste0(prefix, "2_rmse_heatmap.pdf"))

  plot_ratio_boxplot(results_df, fixed_s = fixed_s, fixed_beta = fixed_beta,
                     save_pdf = save_pdf, file = paste0(prefix, "3_ratio_boxplot.pdf"))

  plot_tpr_fpr(results_df,
               save_pdf = save_pdf, file = paste0(prefix, "4_tpr_fpr.pdf"))

  plot_mi_scatter(results_df,
                  save_pdf = save_pdf, file = paste0(prefix, "5_mi_scatter.pdf"))

  plot_bias_vs_p(results_df, fixed_s = fixed_s, fixed_n = fixed_n,
                 save_pdf = save_pdf, file = paste0(prefix, "6_bias_vs_p.pdf"))

  plot_nselected(results_df, fixed_beta = fixed_beta,
                 save_pdf = save_pdf, file = paste0(prefix, "7_nselected.pdf"))

  cat("All figures written.\n")
}


################################################################

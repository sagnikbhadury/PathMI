## Generates the figures and numbers shown in README.md from simulated data.
## Run from the repository root after installing the package:
##   Rscript tools/readme_figures.R [cores]
## No real embeddings are used anywhere in this script.
## The design mirrors the paper: four pooling rules x four projections (16
## configurations) per cohort, one shared permutation matrix, max-T over the
## projections of each pooling and over all 16 (p_max16), the numerical
## tolerance from replicate and sign-flip refits (the null covariate block is a
## specification check, reported separately). In this simulator the positive
## state also changes the within-slide dispersion of every tile coordinate (the
## shift moves the tumour centre relative to the stroma centre), which moment
## pooling captures almost perfectly; the rest of the example therefore follows
## the mean-pooling arm (the arm whose tile decomposition is exact).

suppressPackageStartupMessages(library(PathMI))
args <- commandArgs(TRUE); CORES <- if (length(args)) as.integer(args[1]) else 4L
FIG <- "man/figures"; dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
out <- file("tools/readme_numbers.txt", "w"); say <- function(...) { cat(..., "\n"); cat(..., "\n", file = out) }
png2 <- function(f, w, h) grDevices::png(file.path(FIG, f), width = w, height = h, units = "in", res = 150)
TEAL <- "#1B7F79"; ORANGE <- "#D55E00"; GREY <- "grey55"; set.seed(1)
D_PROJ <- 20L; B_NULL <- 100L; B_HALF <- 100L
DECOMPS <- c("pca", "cosine", "rbf", "poly"); POOLS <- c("mean", "norm", "moment", "attention")
KEYS <- as.vector(outer(DECOMPS, POOLS, paste, sep = "_"))          # projection_pooling, as in the paper

## ---- 1. Two simulated cohorts, 16 configurations each ---------------------
cohorts <- list(strong = simulate_cohort(400, prevalence = 0.3, D = 64, tiles_mean = 150, effect = 0.45, p_multi = 0.08, tile_sd = 0.4, seed = 11),
                weak   = simulate_cohort(400, prevalence = 0.3, D = 64, tiles_mean = 150, effect = 0.12, p_multi = 0.08, tile_sd = 0.4, seed = 12))
res <- list()
for (nm in names(cohorts)) {
  sim <- cohorts[[nm]]
  Xs <- lapply(setNames(POOLS, POOLS), function(pl) pool_cohort(sim$tiles, sim$patient, pl))
  y <- sim$y[rownames(Xs$mean)]
  say(sprintf("[%s] %d slides, %d patients, prevalence %.2f, D = %d", nm, length(sim$tiles), nrow(Xs$mean), mean(y), ncol(Xs$mean)))
  Zs <- lapply(setNames(KEYS, KEYS), function(k) project_embeddings(Xs[[sub(".*_", "", k)]], sub("_.*", "", k), D_PROJ))
  fits <- lapply(Zs, pathmi_fit, y = y)
  perms <- draw_permutations(length(y), B_NULL, seed = 42)       # one matrix, shared by all 16 configurations
  nulls <- mapply(function(Z, f) pathmi_null(Z, y, perms, observed = f, cores = CORES), Zs, fits, SIMPLIFY = FALSE)
  sel4 <- lapply(setNames(POOLS, POOLS), function(pl) select_representation(nulls[paste(DECOMPS, pl, sep = "_")]))
  sel16 <- select_representation(nulls)
  va <- verify_permutation_alignment(nulls[c(sel16$selected, setdiff(KEYS, sel16$selected)[1])], Zs, y, b = c(1L, 37L))
  for (k in KEYS) say(sprintf("[%s] %-16s estimate %.4f  null %.4f +/- %.4f  p %.3f  k* %d", nm, k, fits[[k]]$estimate, nulls[[k]]$mean, nulls[[k]]$sd, nulls[[k]]$p_value, fits[[k]]$k))
  for (pl in POOLS) say(sprintf("[%s] %s pooling: max-T over 4 projections selects %s (%.4f), p_maxT = %.3f", nm, pl, sel4[[pl]]$selected, sel4[[pl]]$estimate, sel4[[pl]]$p_maxT))
  say(sprintf("[%s] max over 16 configurations: %s (%.4f), p_max16 = %.3f | alignment refits max |diff| = %.2g", nm, sel16$selected, sel16$estimate, sel16$p_maxT, max(va$abs_diff)))
  mk <- paste(DECOMPS, "mean", sep = "_")
  rep_d <- vapply(mk, function(k) replicate_refit(Zs[[k]], y, fits[[k]]), numeric(1))
  flip_d <- vapply(mk, function(k) sign_flip_refit(Zs[[k]], y, fits[[k]]), numeric(1))
  null_d <- vapply(mk, function(k) null_block_refit(Zs[[k]], y, fits[[k]]), numeric(1))
  res[[nm]] <- list(sim = sim, Xs = Xs, y = y, Zs = Zs, fits = fits, nulls = nulls, sel4 = sel4, sel16 = sel16, rep_d = rep_d, flip_d = flip_d, null_d = null_d)
}
tol <- numerical_tolerance(replicate = c(res$strong$rep_d, res$weak$rep_d), sign_flip = c(res$strong$flip_d, res$weak$flip_d),
                           specification = c(res$strong$null_d, res$weak$null_d))
EPS <- tol$tolerance
say(sprintf("numerical tolerance: replicate max %.4f, sign-flip max %.4f -> eps_tol %.4f | specification (null block) max %.4f, %d of %d exact",
            tol$per_identity[1], tol$per_identity[2], EPS, tol$specification$max, tol$specification$exact, tol$specification$n))

for (nm in names(res)) {
  r <- res[[nm]]; key <- r$sel4$mean$selected; dc <- sub("_.*", "", key); pl <- sub(".*_", "", key); Z <- r$Zs[[key]]; X <- r$Xs[[pl]]
  r$key <- key; r$half <- half_sampling(Z, r$y, B = B_HALF, cores = CORES)
  r$usable <- usable_recoverability(X, r$y, dc, D_PROJ, K = 5, R = 5, group = NULL)
  C <- cbind(stage = as.vector(scale(r$sim$covariate[rownames(X)])), log_tiles = as.vector(scale(log(tapply(vapply(r$sim$tiles, nrow, 1), r$sim$patient, sum)[rownames(X)]))))
  r$inc <- incremental_recoverability(Z, C, r$y)
  r$curve <- learning_curve(X, r$y, dc, D_PROJ, n_grid = c(50, 100, 150, 200), B = 40, eps_tol = EPS,
                            full = list(estimate = r$fits[[key]]$estimate, se = r$half$se, null_mean = r$nulls[[key]]$mean, usable = r$usable$usable[["ridge"]]), cores = CORES)
  say(sprintf("[%s] recorded configuration %s | half-sampling s.e. %.4f, 90%% interval %.4f-%.4f", nm, key, r$half$se, r$half$interval[1], r$half$interval[2]))
  say(sprintf("[%s] usable ridge %.4f (AUC %.2f), usable PathMI readout %.4f (AUC %.2f)", nm, r$usable$usable[["ridge"]], r$usable$auc[["ridge"]], r$usable$usable[["pathmi"]], r$usable$auc[["pathmi"]]))
  say(sprintf("[%s] incremental beyond covariates: covariates %.4f, combined %.4f, incremental %.4f", nm, r$inc$covariates, r$inc$combined, r$inc$incremental))
  say(sprintf("[%s] N_evidence %s", nm, attr(r$curve, "N_evidence"))); print(r$curve)
  res[[nm]] <- r
}

## ---- 2. Tiles: exact decomposition, surprise maps, burden (fixed linear readout) ----
s1 <- res$strong; td <- tile_decomposition(s1$sim$tiles, s1$sim$patient, s1$sim$y, d = D_PROJ)
say(sprintf("tile identity max error %.2e; linear arm estimate %.4f, k* %d", td$max_identity_error, td$fit$estimate, td$fit$k))
bd <- tissue_burden(td, m_grid = c(2, 5, 10, 15, 20, 30, 50, 75, 100, 150), R = 20)
say(sprintf("burden N90 %s, N95 %s tiles per slide", bd$N[1], bd$N[2])); print(bd$curve)
for (nm in names(res)) {
  r <- res[[nm]]; key <- r$key
  rec <- recoverability_record(list(cohort = paste("simulated", nm), model = "simulated encoder (D = 64)", pooling = sub(".*_", "", key), projection = sub("_.*", "", key),
                                    d = D_PROJ, unit = "patient"),
                               r$fits[[key]], r$nulls[[key]], r$half, EPS, incremental = r$inc, usable = r$usable,
                               burden = if (nm == "strong") bd else NULL, curve = r$curve, selection = r$sel4$mean)
  print(rec); capture.output(print(rec), file = out, append = FALSE); record_to_json(rec, sprintf("tools/record_%s.json", nm)); res[[nm]]$record <- rec
}

## ---- 3. Hidden context ---------------------------------------------------
hc <- do.call(rbind, lapply(c(1, 0.75, 0.5, 0.25), function(pa) {
  est <- unlist(lapply(1:10, function(r) { s <- simulate_hidden_context(400, pi_A = pa, seed = 100 * r, mc = if (r == 1) 2e5 else 0); pathmi_fit(s$Z, s$y)$estimate }))
  tr <- simulate_hidden_context(400, pi_A = pa, seed = 100, mc = 2e5)$truth
  data.frame(pi_A = pa, mean = mean(est), lo = min(est), hi = max(est), I_A = tr[["I_A"]], I_pooled = tr[["I_pooled"]], I_cond = tr[["I_conditional"]])
})); print(hc); utils::write.csv(hc, "tools/hidden_context.csv", row.names = FALSE)
for (i in seq_len(nrow(hc))) say(sprintf("hidden context pi_A %.2f: estimate %.4f, I(Z;Y) %.4f, pi_A*I_A %.4f, I_A %.4f", hc$pi_A[i], hc$mean[i], hc$I_pooled[i], hc$I_cond[i], hc$I_A[i]))

## ---- 4. Figures ----------------------------------------------------------
par_std <- function() graphics::par(mar = c(4, 4.2, 2.6, 1), mgp = c(2.4, 0.7, 0), las = 1, cex.main = 1, font.main = 1, bty = "l")

## Figure 1: estimate, null, walk, the 16 configurations
png2("fig1_estimate.png", 11, 3.6); graphics::layout(matrix(1:3, 1), widths = c(1, 1, 1.5)); par_std()
s <- res$strong; plot_null(s$nulls[[s$key]], main = sprintf("a  Estimate against its null (%s)", gsub("_", ", ", s$key)))
plot_walk(s$fits[[s$key]], main = "b  Marginal-likelihood walk")
est <- sapply(res, function(r) sapply(KEYS, function(k) r$fits[[k]]$estimate)); thr <- sapply(res, function(r) sapply(KEYS, function(k) r$nulls[[k]]$mean + 2 * r$nulls[[k]]$sd))
graphics::par(mar = c(6.2, 4.2, 2.6, 1))
xx <- graphics::barplot(t(est), beside = TRUE, col = c(TEAL, "#9FD3CF"), border = NA, ylim = c(0, max(est) * 1.25), ylab = "Estimate (nats)",
                        main = "c  Sixteen configurations, two cohorts", names.arg = gsub("_", "\n", KEYS), cex.names = 0.55, las = 2)
graphics::segments(xx - 0.45, t(thr), xx + 0.45, t(thr), lty = 2, col = "grey30")
graphics::legend("topright", c("strong signal", "weak signal", "null mean + 2 s.d."), fill = c(TEAL, "#9FD3CF", NA), border = NA, lty = c(NA, NA, 2), bty = "n", cex = 0.8)
grDevices::dev.off()

## Figure 2: three quantities and the learning curves
png2("fig2_evidence.png", 11, 3.6); graphics::par(mfrow = c(1, 3)); par_std()
q <- sapply(res, function(r) c(represented = r$fits[[r$key]]$estimate, incremental = r$inc$incremental, usable = r$usable$usable[["ridge"]]))
graphics::barplot(q, beside = TRUE, col = c(TEAL, "#7FB8B3", ORANGE), border = NA, ylab = "nats", main = "a  Represented, incremental, usable", ylim = c(min(0, q), max(q) * 1.3))
graphics::legend("topright", c("represented I(Z;Y)", "incremental beyond covariates", "usable (held-out ridge)"), fill = c(TEAL, "#7FB8B3", ORANGE), border = NA, bty = "n", cex = 0.85)
for (nm in names(res)) { cv <- res[[nm]]$curve
  graphics::plot(cv$n, cv$estimate, ylim = range(0, cv$estimate + 1.8 * cv$se, cv$threshold), pch = 19, col = TEAL, xlab = "Patients n", ylab = "Estimate (nats)",
                 main = sprintf("%s  Learning curve, %s signal", c(strong = "b", weak = "c")[nm], nm))
  graphics::arrows(cv$n, cv$lower, cv$n, cv$estimate + 1.645 * cv$se, angle = 90, code = 3, length = 0.03, col = TEAL)
  graphics::lines(cv$n, cv$threshold, lty = 2, col = "grey30"); graphics::lines(cv$n, cv$null_mean, col = GREY)
  graphics::points(cv$n[cv$resolved], cv$lower[cv$resolved], pch = 17, col = ORANGE)
  ne <- attr(cv, "N_evidence"); graphics::mtext(sprintf("N_evidence = %s", ifelse(is.na(ne), "not reached", ne)), 3, -1.3, adj = 0.05, cex = 0.8)
  if (nm == "weak") graphics::legend("topright", c("estimate, 90% interval", "null mean", "null + eps_tol", "lower bound resolved"), col = c(TEAL, GREY, "grey30", ORANGE), lty = c(1, 1, 2, NA), pch = c(19, NA, NA, 17), bty = "n", cex = 0.8) }
grDevices::dev.off()

## Figure 3: tile maps and burden
pt <- td$patients; ppos <- pt$patient[pt$y == 1][order(-pt$p[pt$y == 1])][1]; pneg <- pt$patient[pt$y == 0][order(pt$p[pt$y == 0])][1]
pmid <- pt$patient[pt$y == 1][order(abs(pt$p[pt$y == 1] - stats::median(pt$p[pt$y == 1])))][1]
slides <- match(c(ppos, pmid, pneg), s1$sim$patient)
png2("fig3_tiles.png", 11, 6.4); graphics::layout(matrix(c(1, 2, 3, 7, 4, 5, 6, 7), 2, byrow = TRUE), widths = c(1, 1, 1, 1.25)); graphics::par(mar = c(0.5, 0.5, 2.4, 2.2), cex.main = 0.95, font.main = 1)
for (sl in slides) { tis <- as.numeric(s1$sim$tissue[[sl]] == "tumour")
  graphics::plot(range(s1$sim$coords[[sl]][, 1]) + c(-.5, .5), range(s1$sim$coords[[sl]][, 2]) + c(-.5, .5), type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "", main = sprintf("Tissue, %s", s1$sim$patient[sl]))
  xy <- s1$sim$coords[[sl]]; graphics::rect(xy[, 1] - .5, xy[, 2] - .5, xy[, 1] + .5, xy[, 2] + .5, col = ifelse(tis == 1, "#C9A3C4", "#F1E3D3"), border = NA) }
for (sl in slides) plot_tile_map(td, sl, coords = s1$sim$coords[[sl]])
graphics::par(mar = c(4.2, 4.2, 2.6, 1), mgp = c(2.4, 0.7, 0), las = 1, bty = "l", cex.main = 1, font.main = 1); plot_burden(bd, main = "Tissue-sampling burden (linear readout)")
grDevices::dev.off()

## Figure 4: hidden context
png2("fig4_hidden_context.png", 5.2, 3.8); par_std()
graphics::plot(hc$pi_A, hc$I_A, type = "b", pch = 1, ylim = c(0, max(hc$I_A, hc$hi) * 1.1), xlim = c(0.2, 1.05), col = "grey30", lty = 3,
               xlab = expression("Share of context A, " * pi[A]), ylab = "nats", main = "Pooled measurement hides a context")
graphics::lines(hc$pi_A, hc$I_cond, type = "b", pch = 0, col = "grey30", lty = 2); graphics::lines(hc$pi_A, hc$I_pooled, type = "b", pch = 15, col = "grey30")
graphics::arrows(hc$pi_A, hc$lo, hc$pi_A, hc$hi, angle = 90, code = 3, length = 0.03, col = TEAL); graphics::points(hc$pi_A, hc$mean, pch = 19, col = TEAL)
graphics::legend("topleft", c("I(Z;Y | S=A)", "pi_A I(Z;Y | S=A)", "pooled I(Z;Y)", "PathMI (10 replicates)"), pch = c(1, 0, 15, 19), lty = c(3, 2, 1, NA), col = c("grey30", "grey30", "grey30", TEAL), bty = "n", cex = 0.8)
grDevices::dev.off()
close(out); cat("done\n")

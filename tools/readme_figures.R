## Generates the figures and numbers shown in README.md from simulated data.
## Run from the repository root after installing the package:
##   Rscript tools/readme_figures.R [cores]
## No real embeddings are used anywhere in this script.

suppressPackageStartupMessages(library(PathMI))
args <- commandArgs(TRUE); CORES <- if (length(args)) as.integer(args[1]) else 4L
FIG <- "man/figures"; dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
out <- file("tools/readme_numbers.txt", "w"); say <- function(...) { cat(..., "\n"); cat(..., "\n", file = out) }
png2 <- function(f, w, h) grDevices::png(file.path(FIG, f), width = w, height = h, units = "in", res = 150)
TEAL <- "#1B7F79"; ORANGE <- "#D55E00"; GREY <- "grey55"; set.seed(1)
D_PROJ <- 20L; B_NULL <- 100L; B_HALF <- 100L
DECOMPS <- c("pca", "cosine", "rbf", "poly")

## ---- 1. Two simulated cohorts --------------------------------------------
cohorts <- list(strong = simulate_cohort(400, prevalence = 0.3, D = 64, tiles_mean = 150, effect = 0.45, p_multi = 0.08, tile_sd = 0.4, seed = 11),
                weak   = simulate_cohort(400, prevalence = 0.3, D = 64, tiles_mean = 150, effect = 0.12, p_multi = 0.08, tile_sd = 0.4, seed = 12))
res <- list()
for (nm in names(cohorts)) {
  sim <- cohorts[[nm]]; X <- pool_cohort(sim$tiles, sim$patient, "mean"); y <- sim$y[rownames(X)]
  say(sprintf("[%s] %d slides, %d patients, prevalence %.2f, D = %d", nm, length(sim$tiles), nrow(X), mean(y), ncol(X)))
  Zs <- lapply(setNames(DECOMPS, DECOMPS), function(dc) project_embeddings(X, dc, D_PROJ))
  fits <- lapply(Zs, pathmi_fit, y = y)
  perms <- draw_permutations(length(y), B_NULL, seed = 42)
  nulls <- mapply(function(Z, f) pathmi_null(Z, y, perms, observed = f, cores = CORES), Zs, fits, SIMPLIFY = FALSE)
  sel <- select_representation(nulls)
  for (dc in DECOMPS) say(sprintf("[%s] %-6s estimate %.4f  null %.4f +/- %.4f  p %.3f  k* %d", nm, dc, fits[[dc]]$estimate, nulls[[dc]]$mean, nulls[[dc]]$sd, nulls[[dc]]$p_value, fits[[dc]]$k))
  say(sprintf("[%s] max-T selects %s: %.4f, p_maxT = %.3f", nm, sel$selected, sel$estimate, sel$p_maxT))
  rep_d <- vapply(DECOMPS, function(dc) replicate_refit(Zs[[dc]], y, fits[[dc]]), numeric(1))
  flip_d <- vapply(DECOMPS, function(dc) sign_flip_refit(Zs[[dc]], y, fits[[dc]]), numeric(1))
  null_d <- vapply(DECOMPS, function(dc) null_block_refit(Zs[[dc]], y, fits[[dc]]), numeric(1))
  res[[nm]] <- list(sim = sim, X = X, y = y, Zs = Zs, fits = fits, nulls = nulls, sel = sel, rep_d = rep_d, flip_d = flip_d, null_d = null_d)
}
eps <- reproducibility_bound(replicate = c(res$strong$rep_d, res$weak$rep_d), sign_flip = c(res$strong$flip_d, res$weak$flip_d), null_block = c(res$strong$null_d, res$weak$null_d))
EPS <- eps$bound
say(sprintf("eps_rep: replicate max %.4f, sign-flip max %.4f, null-block max %.4f, bound %.4f", eps$per_identity[1], eps$per_identity[2], eps$per_identity[3], EPS))

for (nm in names(res)) {
  r <- res[[nm]]; dc <- r$sel$selected; Z <- r$Zs[[dc]]
  r$half <- half_sampling(Z, r$y, B = B_HALF, cores = CORES)
  r$usable <- usable_recoverability(r$X, r$y, dc, D_PROJ, K = 5, R = 5, group = NULL)
  C <- cbind(stage = as.vector(scale(r$sim$covariate[rownames(r$X)])), log_tiles = as.vector(scale(log(tapply(vapply(r$sim$tiles, nrow, 1), r$sim$patient, sum)[rownames(r$X)]))))
  r$inc <- incremental_recoverability(Z, C, r$y)
  r$curve <- learning_curve(r$X, r$y, dc, D_PROJ, n_grid = c(50, 100, 150, 200, 250, 300), B = 40, eps_rep = EPS,
                            full = list(estimate = r$fits[[dc]]$estimate, se = r$half$se, null_mean = r$nulls[[dc]]$mean, usable = r$usable$usable[["ridge"]]), cores = CORES)
  say(sprintf("[%s] half-sampling s.e. %.4f, 90%% interval %.4f-%.4f", nm, r$half$se, r$half$interval[1], r$half$interval[2]))
  say(sprintf("[%s] usable ridge %.4f (AUC %.2f), usable PathMI readout %.4f (AUC %.2f)", nm, r$usable$usable[["ridge"]], r$usable$auc[["ridge"]], r$usable$usable[["pathmi"]], r$usable$auc[["pathmi"]]))
  say(sprintf("[%s] attributable: covariates %.4f, combined %.4f, incremental %.4f", nm, r$inc$covariates, r$inc$combined, r$inc$incremental))
  say(sprintf("[%s] N_evidence %s", nm, attr(r$curve, "N_evidence"))); print(r$curve)
  res[[nm]] <- r
}

## ---- 2. Tiles: exact decomposition, surprise maps, burden -----------------
s1 <- res$strong; td <- tile_decomposition(s1$sim$tiles, s1$sim$patient, s1$sim$y, d = D_PROJ)
say(sprintf("tile identity max error %.2e; linear arm estimate %.4f, k* %d", td$max_identity_error, td$fit$estimate, td$fit$k))
bd <- tissue_burden(td, m_grid = c(2, 5, 10, 15, 20, 30, 50, 75, 100, 150), R = 20)
say(sprintf("burden N90 %s, N95 %s tiles per slide", bd$N[1], bd$N[2])); print(bd$curve)
for (nm in names(res)) {
  r <- res[[nm]]; dc <- r$sel$selected
  rec <- recoverability_record(list(cohort = paste("simulated", nm), model = "simulated encoder (D = 64)", pooling = "mean", projection = dc, d = D_PROJ, unit = "patient"),
                               r$fits[[dc]], r$nulls[[dc]], r$half, EPS, attributable = r$inc, usable = r$usable, burden = if (nm == "strong") bd else NULL, curve = r$curve)
  print(rec); capture.output(print(rec), file = out, append = FALSE); record_to_json(rec, sprintf("tools/record_%s.json", nm)); res[[nm]]$record <- rec
}

## ---- 3. Hidden context ---------------------------------------------------
hc <- do.call(rbind, lapply(c(1, 0.75, 0.5, 0.25), function(pa) {
  est <- unlist(parallel_or_lapply <- lapply(1:10, function(r) { s <- simulate_hidden_context(400, pi_A = pa, seed = 100 * r, mc = if (r == 1) 2e5 else 0); pathmi_fit(s$Z, s$y)$estimate }))
  tr <- simulate_hidden_context(400, pi_A = pa, seed = 100, mc = 2e5)$truth
  data.frame(pi_A = pa, mean = mean(est), lo = min(est), hi = max(est), I_A = tr[["I_A"]], I_pooled = tr[["I_pooled"]], I_cond = tr[["I_conditional"]])
})); print(hc); utils::write.csv(hc, "tools/hidden_context.csv", row.names = FALSE)
for (i in seq_len(nrow(hc))) say(sprintf("hidden context pi_A %.2f: estimate %.4f, I(Z;Y) %.4f, pi_A*I_A %.4f, I_A %.4f", hc$pi_A[i], hc$mean[i], hc$I_pooled[i], hc$I_cond[i], hc$I_A[i]))

## ---- 4. Figures ----------------------------------------------------------
par_std <- function() graphics::par(mar = c(4, 4.2, 2.6, 1), mgp = c(2.4, 0.7, 0), las = 1, cex.main = 1, font.main = 1, bty = "l")

## Figure 1: estimate, null, walk, selection
png2("fig1_estimate.png", 11, 3.6); graphics::par(mfrow = c(1, 3)); par_std()
s <- res$strong; dc <- s$sel$selected; plot_null(s$nulls[[dc]], main = sprintf("a  Estimate against its null (%s)", dc))
plot_walk(s$fits[[dc]], main = "b  Marginal-likelihood walk")
est <- sapply(res, function(r) sapply(DECOMPS, function(k) r$fits[[k]]$estimate)); nm0 <- sapply(res, function(r) sapply(DECOMPS, function(k) r$nulls[[k]]$mean))
ns0 <- sapply(res, function(r) sapply(DECOMPS, function(k) r$nulls[[k]]$sd))
xx <- graphics::barplot(t(est), beside = TRUE, col = c(TEAL, "#9FD3CF"), border = NA, ylim = c(0, max(est) * 1.25), ylab = "Estimate (nats)", main = "c  Four projections, two cohorts")
graphics::segments(xx - 0.45, t(nm0 + EPS + 2 * ns0), xx + 0.45, t(nm0 + EPS + 2 * ns0), lty = 2, col = "grey30")
graphics::legend("topright", c("strong signal", "weak signal", "null mean + 2 s.d. + eps_rep"), fill = c(TEAL, "#9FD3CF", NA), border = NA, lty = c(NA, NA, 2), bty = "n", cex = 0.85)
grDevices::dev.off()

## Figure 2: three quantities and the learning curves
png2("fig2_evidence.png", 11, 3.6); graphics::par(mfrow = c(1, 3)); par_std()
q <- sapply(res, function(r) c(represented = r$fits[[r$sel$selected]]$estimate, attributable = r$inc$incremental, usable = r$usable$usable[["ridge"]]))
graphics::barplot(q, beside = TRUE, col = c(TEAL, "#7FB8B3", ORANGE), border = NA, ylab = "nats", main = "a  Represented, attributable, usable", ylim = c(0, max(q) * 1.3))
graphics::legend("topright", c("represented I(Z;Y)", "attributable I(Z;Y|C)", "usable (held-out ridge)"), fill = c(TEAL, "#7FB8B3", ORANGE), border = NA, bty = "n", cex = 0.85)
for (nm in names(res)) { cv <- res[[nm]]$curve
  graphics::plot(cv$n, cv$estimate, ylim = range(0, cv$estimate + 1.8 * cv$se, cv$threshold), pch = 19, col = TEAL, xlab = "Patients n", ylab = "Estimate (nats)",
                 main = sprintf("%s  Learning curve, %s signal", c(strong = "b", weak = "c")[nm], nm))
  graphics::arrows(cv$n, cv$lower, cv$n, cv$estimate + 1.645 * cv$se, angle = 90, code = 3, length = 0.03, col = TEAL)
  graphics::lines(cv$n, cv$threshold, lty = 2, col = "grey30"); graphics::lines(cv$n, cv$null_mean, col = GREY)
  graphics::points(cv$n[cv$resolved], cv$lower[cv$resolved], pch = 17, col = ORANGE)
  ne <- attr(cv, "N_evidence"); graphics::mtext(sprintf("N_evidence = %s", ifelse(is.na(ne), "not reached", ne)), 3, -1.3, adj = 0.05, cex = 0.8)
  if (nm == "weak") graphics::legend("topright", c("estimate, 90% interval", "null mean", "null + eps_rep", "lower bound resolved"), col = c(TEAL, GREY, "grey30", ORANGE), lty = c(1, 1, 2, NA), pch = c(19, NA, NA, 17), bty = "n", cex = 0.8) }
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
graphics::par(mar = c(4.2, 4.2, 2.6, 1), mgp = c(2.4, 0.7, 0), las = 1, bty = "l", cex.main = 1, font.main = 1); plot_burden(bd, main = "Tissue-sampling burden")
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

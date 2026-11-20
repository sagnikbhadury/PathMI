## Base-graphics plots for the measurement objects.

#' Plot the estimate against its permutation null
#' @param null a [pathmi_null()] result.
#' @param main title.
#' @param ... passed to `hist`.
#' @export
plot_null <- function(null, main = "Permutation null", ...) {
  est <- null$observed$estimate; br <- pretty(c(0, max(null$null, est) * 1.05), 30)
  graphics::hist(null$null, breaks = br, col = "grey80", border = "white", main = main, xlab = "Estimate under permuted labels (nats)", xlim = range(br), ...)
  graphics::abline(v = est, col = "#1B7F79", lwd = 3)
  graphics::legend("top", bty = "n", legend = c(sprintf("estimate %.3f", est), sprintf("null mean %.4f", null$mean), sprintf("p = %.3f", null$p_value)), text.col = c("#1B7F79", "grey30", "grey30"))
}

#' Plot the marginal-likelihood walk of a fit
#' @param fit a [pathmi_fit()] result.
#' @param main title.
#' @export
plot_walk <- function(fit, main = "Marginal-likelihood walk") {
  w <- fit$walk; ev <- (w - max(w))
  graphics::plot(seq_along(w), ev, type = "b", pch = 19, cex = 0.5, ylim = c(min(ev), 0.12 * diff(range(ev))), xlab = "Coordinates in the model, k", ylab = "log posterior score (relative to max)", main = main)
  graphics::points(fit$k, 0, pch = 19, cex = 1.6, col = "#1B7F79"); graphics::text(fit$k, 0, sprintf("k* = %d", fit$k), pos = 3, offset = 0.8, col = "#1B7F79")
}

#' Plot a learning curve
#' @param curve a [learning_curve()] result.
#' @param main title.
#' @export
plot_learning_curve <- function(curve, main = "How many patients a claim needs") {
  m <- curve$lower - curve$threshold
  op <- graphics::par(mar = c(4.2, 4.2, 3, 4.2)); on.exit(graphics::par(op))
  graphics::plot(curve$n, m, type = "b", pch = ifelse(curve$resolved, 19, 1), col = "#1B7F79", lwd = 2, xlab = "Patients n", ylab = "Lower bound minus threshold (nats)", main = main)
  graphics::abline(h = 0)
  if (any(!is.na(curve$usable))) { sc <- max(abs(m), na.rm = TRUE) / max(abs(curve$usable), na.rm = TRUE)
    graphics::lines(curve$n, curve$usable * sc, type = "b", lty = 2, col = "#D55E00", pch = 19, cex = 0.7)
    graphics::axis(4, at = pretty(curve$usable * sc), labels = signif(pretty(curve$usable * sc) / sc, 2), col.axis = "#D55E00")
    graphics::mtext("Usable recoverability (nats)", 4, 2.6, col = "#D55E00") }
  ne <- attr(curve, "N_evidence"); graphics::mtext(sprintf("N_evidence: %s", ifelse(is.na(ne), "not reached", ne)), 3, 0.2, cex = 0.8)
}

#' Plot the specific surprise of the tiles of one slide
#'
#' Each tile is drawn as a square at its position and coloured by its specific
#' surprise (or by `values`, if given) on a scale shared across slides.
#' @param decomp a [tile_decomposition()] result.
#' @param slide index of the slide.
#' @param coords optional two-column matrix of tile positions (grid units); a
#'   square grid is used when absent.
#' @param values optional per-tile values to show instead of the surprise.
#' @param zlim colour limits (default: 0 to the 99th percentile of the surprise
#'   over all tiles).
#' @param main title.
#' @param legend draw a colour bar.
#' @export
plot_tile_map <- function(decomp, slide, coords = NULL, values = NULL, zlim = NULL, main = NULL, legend = TRUE) {
  d <- decomp$tiles[decomp$tiles$slide == slide, ]; n <- nrow(d)
  if (is.null(coords)) { w <- ceiling(sqrt(n)); coords <- cbind((seq_len(n) - 1) %% w + 1, (seq_len(n) - 1) %/% w + 1) }
  v <- if (is.null(values)) d$surprise else values
  if (is.null(zlim)) zlim <- c(0, stats::quantile(decomp$tiles$surprise, 0.99))
  pal <- grDevices::hcl.colors(100, "inferno"); col <- pal[pmax(1, pmin(100, 1 + floor(99 * (v - zlim[1]) / diff(zlim))))]
  pp <- decomp$patients$p[decomp$patients$patient == d$patient[1]]; yy <- decomp$patients$y[decomp$patients$patient == d$patient[1]]
  graphics::plot(range(coords[, 1]) + c(-0.5, 0.5), range(coords[, 2]) + c(-0.5, 0.5), type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "",
                 main = if (is.null(main)) sprintf("%s (y = %d, fitted p = %.2f)", d$patient[1], yy, pp) else main)
  graphics::rect(coords[, 1] - 0.5, coords[, 2] - 0.5, coords[, 1] + 0.5, coords[, 2] + 0.5, col = col, border = NA)
  if (legend) { u <- graphics::par("usr"); x0 <- u[2] - 0.02 * diff(u[1:2]); yb <- seq(u[3], u[3] + 0.5 * diff(u[3:4]), length.out = 101)
    graphics::rect(x0, yb[-101], x0 + 0.03 * diff(u[1:2]), yb[-1], col = pal, border = NA, xpd = TRUE)
    graphics::text(x0 + 0.035 * diff(u[1:2]), c(yb[1], yb[101]), signif(zlim, 2), adj = 0, cex = 0.7, xpd = TRUE) }
  invisible(col)
}

#' Plot the tissue-sampling burden
#' @param burden a [tissue_burden()] result.
#' @param main title.
#' @export
plot_burden <- function(burden, main = "Tissue-sampling burden") {
  cv <- burden$curve
  graphics::plot(cv$m, cv$ratio, log = "x", type = "n", ylim = c(max(-0.5, min(cv$lo)), 1.08), xlab = "Tiles per slide, m", ylab = "Plug-in estimate / whole-slide value", main = main)
  graphics::polygon(c(cv$m, rev(cv$m)), c(cv$lo, rev(cv$hi)), col = "grey88", border = NA); graphics::lines(cv$m, cv$ratio, lwd = 2, col = "#5B5EA6")
  graphics::abline(h = c(0.9, 0.95), lty = 2, col = "grey50")
  graphics::abline(v = burden$N, lty = 3, col = "#5B5EA6")
  graphics::legend("bottomright", bty = "n", legend = paste(names(burden$N), burden$N, sep = " = "), cex = 0.9)
}

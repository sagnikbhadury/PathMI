## Evidence sufficiency: patient-number learning curves, N_evidence and the
## status of a record; and the Recoverability Record itself.

#' Patient-number learning curve and N_evidence
#'
#' For each `n` in `n_grid`, draws `B` label-stratified subsamples of `n`
#' patients without replacement, refits the projection, the estimator, the
#' estimator on permuted labels (the null at that `n`) and, optionally, the
#' usable recoverability on each. The standard error of a study of `n`
#' independent patients is the spread of the subsample estimates divided by the
#' finite-population factor \eqn{(1 - n/N)^{1/2}}. A claim is resolved at `n`
#' when \eqn{\bar I(n) - z\,\mathrm{s.e.}(n) > \bar I_0(n) + \epsilon_{rep}};
#' `N_evidence` is the smallest such `n`.
#'
#' @param X pooled embeddings (rows = patients).
#' @param y binary labels.
#' @param decomp,d projection.
#' @param n_grid subsample sizes (each below the number of patients).
#' @param B subsamples per size.
#' @param eps_rep reproducibility bound (nats).
#' @param full optional list with the full-cohort `estimate`, `null_mean` and
#'   `se` (for example from [pathmi_null()] and [half_sampling()]), appended as
#'   the last point of the curve.
#' @param usable logical; also cross-fit the ridge usable recoverability.
#' @param level one-sided coverage of the lower bound (0.95 gives the 90 percent
#'   two-sided interval used in the paper).
#' @param cores parallel workers.
#' @param ... passed to [pathmi_fit()].
#' @return a data frame with one row per size and attribute `N_evidence`.
#' @export
learning_curve <- function(X, y, decomp = "pca", d = 50L, n_grid = c(50L, 100L, 150L, 200L), B = 50L, eps_rep, full = NULL,
                           usable = TRUE, level = 0.95, cores = 1L, ...) {
  X <- as.matrix(X); y <- as.numeric(y); N <- length(y); n_grid <- n_grid[n_grid < N]; z <- stats::qnorm(level)
  force(decomp); force(d); force(usable); dots <- list(...); fitp <- function(Z, yy, seed) do.call(pathmi_fit, c(list(Z, yy, seed = seed), dots))
  one <- function(job) {
    n <- job[1]; b <- job[2]; set.seed(5000L + b)
    n1 <- max(2L, round(n * mean(y))); idx <- sort(c(sample(which(y == 1), n1), sample(which(y == 0), n - n1)))
    Xs <- X[idx, , drop = FALSE]; ys <- y[idx]; Z <- project_embeddings(Xs, decomp, d)
    f <- fitp(Z, ys, 60000L + b)$estimate; set.seed(70000L + b)
    f0 <- fitp(Z, sample(ys), 70000L + b)$estimate
    u <- if (usable) usable_recoverability(Xs, ys, decomp, d, K = 5L, R = 1L, readout = "ridge", seed = 90000L + b)$usable[["ridge"]] else NA_real_
    c(n = n, b = b, I = f, I0 = f0, U = u)
  }
  jobs <- as.matrix(expand.grid(n = n_grid, b = seq_len(B)))
  dr <- as.data.frame(do.call(rbind, par_apply(split(jobs, seq_len(nrow(jobs))), one, cores)))
  cur <- do.call(rbind, lapply(split(dr, dr$n), function(s) { n <- s$n[1]; se <- stats::sd(s$I) / sqrt(1 - n / N)
    data.frame(n = n, estimate = mean(s$I), se = se, null_mean = mean(s$I0), usable = mean(s$U), source = "subsample") }))
  if (!is.null(full)) cur <- rbind(cur, data.frame(n = N, estimate = full$estimate, se = full$se, null_mean = full$null_mean,
                                                    usable = if (is.null(full$usable)) NA_real_ else full$usable, source = "full cohort"))
  cur$lower <- cur$estimate - z * cur$se; cur$threshold <- cur$null_mean + eps_rep; cur$resolved <- cur$lower > cur$threshold
  attr(cur, "N_evidence") <- if (any(cur$resolved)) min(cur$n[cur$resolved]) else NA_integer_
  attr(cur, "draws") <- dr; cur
}

#' Evidence status of a measurement
#' @param estimate,se,null_mean full-cohort estimate, its standard error and
#'   the permutation null mean.
#' @param eps_rep reproducibility bound.
#' @param z quantile for the lower bound (1.645 for a 90 percent interval).
#' @return `"resolved"` when the lower bound exceeds the null plus the bound,
#'   otherwise `"unresolved"`. An unresolved record indicates insufficient
#'   evidence under the evaluated representation, patient budget and protocol;
#'   it is not evidence that the information is absent.
#' @export
evidence_status <- function(estimate, se, null_mean, eps_rep, z = stats::qnorm(0.95)) ifelse(estimate - z * se > null_mean + eps_rep, "resolved", "unresolved")

#' Build a Recoverability Record
#'
#' Collects the measurement with its provenance, protocol, null, resolution,
#' uncertainty, attributable and usable counterparts, burdens and evidence
#' status, so that a model and a state are compared as measurements.
#'
#' @param provenance named list (model, feature key, tile size, projection,
#'   pooling, covariates, ...).
#' @param fit a [pathmi_fit()].
#' @param null a [pathmi_null()].
#' @param uncertainty a [half_sampling()] result.
#' @param eps_rep reproducibility bound.
#' @param attributable optional [incremental_recoverability()] result.
#' @param usable optional [usable_recoverability()] result.
#' @param burden optional [tissue_burden()] result.
#' @param curve optional [learning_curve()] result.
#' @return an object of class `recoverability_record`.
#' @export
recoverability_record <- function(provenance, fit, null, uncertainty, eps_rep, attributable = NULL, usable = NULL, burden = NULL, curve = NULL) {
  st <- evidence_status(fit$estimate, uncertainty$se, null$mean, eps_rep)
  structure(list(schema = "pathmi.recoverability_record/1.0", provenance = provenance,
    measurement = list(n = fit$n, label_entropy = fit$entropy, represented = fit$estimate, share_of_entropy = fit$estimate / fit$entropy, k = fit$k,
                       null_mean = null$mean, null_sd = null$sd, permutation_p = null$p_value, epsilon_rep = eps_rep, se = uncertainty$se, ci90 = uncertainty$interval,
                       attributable = if (is.null(attributable)) NULL else attributable$incremental,
                       usable = if (is.null(usable)) NULL else unname(usable$usable[1]), usable_auc = if (is.null(usable)) NULL else unname(usable$auc[1]),
                       tile_burden = if (is.null(burden)) NULL else burden$N),
    evidence = list(status = st, threshold = null$mean + eps_rep, lower_bound_90 = fit$estimate - stats::qnorm(0.95) * uncertainty$se,
                    N_evidence = if (is.null(curve)) NULL else attr(curve, "N_evidence"), learning_curve = curve),
    interpretation = if (st == "resolved") "The representation carries information about the state above the estimator's null and resolution at this patient budget."
                     else "An unresolved record indicates insufficient evidence under the evaluated representation, patient budget and protocol; it is not evidence that the molecular information is absent."),
    class = "recoverability_record")
}

#' @export
print.recoverability_record <- function(x, ...) {
  m <- x$measurement
  cat("Recoverability Record\n")
  for (nm in names(x$provenance)) cat(sprintf("  %-14s %s\n", nm, paste(x$provenance[[nm]], collapse = ", ")))
  cat(sprintf("  represented    %.4f nats (%.0f%% of H(Y) = %.3f), k* = %d, n = %d\n", m$represented, 100 * m$share_of_entropy, m$label_entropy, m$k, m$n))
  cat(sprintf("  null           %.4f +/- %.4f (p = %.3f) | eps_rep %.4f\n", m$null_mean, m$null_sd, m$permutation_p, m$epsilon_rep))
  cat(sprintf("  90%% interval   %.4f to %.4f (s.e. %.4f)\n", m$ci90[1], m$ci90[2], m$se))
  if (!is.null(m$attributable)) cat(sprintf("  attributable   %.4f\n", m$attributable))
  if (!is.null(m$usable)) cat(sprintf("  usable         %.4f (AUC %.2f)\n", m$usable, m$usable_auc))
  if (!is.null(m$tile_burden)) cat(sprintf("  tiles N90/N95  %s\n", paste(m$tile_burden, collapse = " / ")))
  cat(sprintf("  status         %s%s\n", x$evidence$status, if (is.null(x$evidence$N_evidence)) "" else sprintf(" (N_evidence %s)", ifelse(is.na(x$evidence$N_evidence), "not reached", x$evidence$N_evidence))))
  invisible(x)
}

#' Write a Recoverability Record as JSON
#' @param record a `recoverability_record`.
#' @param file path; if `NULL`, the JSON string is returned.
#' @export
record_to_json <- function(record, file = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
  r <- unclass(record); if (!is.null(r$evidence$learning_curve)) r$evidence$learning_curve <- as.data.frame(r$evidence$learning_curve)
  js <- jsonlite::toJSON(r, auto_unbox = TRUE, pretty = TRUE, digits = 6, na = "null", null = "null")
  if (is.null(file)) js else writeLines(js, file)
}

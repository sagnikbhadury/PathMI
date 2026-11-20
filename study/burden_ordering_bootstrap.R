################################################################
# burden_ordering_bootstrap.R
#
# Robustness of N_q and of the (model, state) ordering to the random
# tile draw, by resampling the 20 draws. A single draw's N_q is the
# first grid point at which a noisy ratio crosses q and is biased low;
# the estimator of N_q is the mean curve over draws. We therefore
# resample the 20 draws with replacement 2000 times, recompute the mean
# curve and N_q for every (model, state) pair, and report (i) the
# bootstrap 90% interval of N_q, (ii) Kendall's W of the ranking of the
# twelve pairs across resamples, (iii) for every pair of pairs the
# fraction of resamples in which their order agrees with the point
# estimate, and (iv) the within-model state orderings and within-state
# model orderings that hold in at least 95% of resamples.
#
# Input: burden/burden_rel_by_draw.csv. Output: burden/burden_Nq_bootstrap.csv,
#        burden/burden_ordering_bootstrap.csv, burden/burden_ordering_summary.txt
################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); OUT <- file.path(TRI, "burden")
set.seed(2026); B <- 2000L
d <- read.csv(file.path(OUT, "burden_rel_by_draw.csv"), stringsAsFactors = FALSE)
nq <- read.csv(file.path(OUT, "burden_Nq.csv"), stringsAsFactors = FALSE)
pairs <- unique(d$pair); ms <- sort(unique(d$m)); R <- max(d$draw)
# rel array: pair x m x draw
A <- array(NA_real_, c(length(pairs), length(ms), R), dimnames = list(pairs, ms, NULL))
for (i in seq_len(nrow(d))) A[d$pair[i], as.character(d$m[i]), d$draw[i]] <- d$rel[i]
Nq_of <- function(rel, q) { i <- which(rel >= q); if (length(i)) ms[min(i)] else NA_integer_ }
boot <- array(NA_integer_, c(B, length(pairs), 2), dimnames = list(NULL, pairs, c("N90", "N95")))
for (b in seq_len(B)) { idx <- sample.int(R, R, replace = TRUE); mc <- apply(A[, , idx, drop = FALSE], c(1, 2), mean)
  boot[b, , 1] <- apply(mc, 1, Nq_of, q = .90); boot[b, , 2] <- apply(mc, 1, Nq_of, q = .95) }
W_of <- function(M) { Rk <- t(apply(M, 1, rank)); k <- ncol(Rk); m <- nrow(Rk); 12 * sum((colSums(Rk) - m * (k + 1) / 2)^2) / (m^2 * (k^3 - k)) }
bs <- data.frame(pair = pairs, N90 = nq$N90[match(pairs, nq$pair)], N90_lo = apply(boot[, , 1], 2, quantile, .05), N90_hi = apply(boot[, , 1], 2, quantile, .95),
                 N95 = nq$N95[match(pairs, nq$pair)], N95_lo = apply(boot[, , 2], 2, quantile, .05), N95_hi = apply(boot[, , 2], 2, quantile, .95))
bs$model <- sub(" / .*", "", bs$pair); bs$state <- sub(".* / ", "", bs$pair)
write.csv(bs, file.path(OUT, "burden_Nq_bootstrap.csv"), row.names = FALSE)
W90 <- W_of(boot[, , 1]); W95 <- W_of(boot[, , 2])
cmb <- t(combn(pairs, 2)); pt <- setNames(nq$N90, nq$pair)
ord <- data.frame(first = cmb[, 1], second = cmb[, 2], N90_first = pt[cmb[, 1]], N90_second = pt[cmb[, 2]])
ord$agree <- apply(cmb, 1, function(pr) { s <- sign(pt[pr[1]] - pt[pr[2]]); if (s == 0) mean(boot[, pr[1], 1] == boot[, pr[2], 1]) else mean(sign(boot[, pr[1], 1] - boot[, pr[2], 1]) == s) })
ord$strict_agree <- apply(cmb, 1, function(pr) { lo <- if (pt[pr[1]] <= pt[pr[2]]) pr[1] else pr[2]; hi <- setdiff(pr, lo); mean(boot[, lo, 1] <= boot[, hi, 1]) })
write.csv(ord, file.path(OUT, "burden_ordering_bootstrap.csv"), row.names = FALSE)
lines <- c(sprintf("Kendall's W of the 12-pair ranking over %d draw-resamples: N90 %.3f | N95 %.3f", B, W90, W95),
           sprintf("pairwise N90 orders (point estimate) that hold (<=) in >= 95%% of resamples: %d of %d; in >= 90%%: %d; in 100%%: %d",
                   sum(ord$strict_agree >= .95), nrow(ord), sum(ord$strict_agree >= .90), sum(ord$strict_agree == 1)))
# within-model state orderings
for (mo in unique(bs$model)) { x <- boot[, paste(mo, c("MSI", "ERG", "ETS"), sep = " / "), 1]
  lines <- c(lines, sprintf("%-12s N90 MSI/ERG/ETS = %d/%d/%d | P(MSI <= ERG) %.2f | P(MSI <= ETS) %.2f | P(ETS <= ERG) %.2f | P(MSI < both fusions) %.2f",
    mo, pt[paste(mo, "MSI", sep = " / ")], pt[paste(mo, "ERG", sep = " / ")], pt[paste(mo, "ETS", sep = " / ")], mean(x[, 1] <= x[, 2]), mean(x[, 1] <= x[, 3]), mean(x[, 3] <= x[, 2]), mean(x[, 1] < x[, 2] & x[, 1] < x[, 3]))) }
# within-state model orderings: top tier (UNI v2, H-optimus-1) vs lower tier (Virchow2, CONCH)
for (st in c("MSI", "ERG", "ETS")) { x <- boot[, paste(c("UNI v2", "H-optimus-1", "Virchow2", "CONCH v1.5"), st, sep = " / "), 1]
  lines <- c(lines, sprintf("%-4s N90 UNI/Hopt/Virchow2/CONCH = %d/%d/%d/%d | P(UNI <= Virchow2) %.2f | P(UNI <= CONCH) %.2f | P(Hopt <= Virchow2) %.2f | P(Hopt <= CONCH) %.2f | P(UNI <= Hopt) %.2f",
    st, pt[paste("UNI v2", st, sep = " / ")], pt[paste("H-optimus-1", st, sep = " / ")], pt[paste("Virchow2", st, sep = " / ")], pt[paste("CONCH v1.5", st, sep = " / ")],
    mean(x[, 1] <= x[, 3]), mean(x[, 1] <= x[, 4]), mean(x[, 2] <= x[, 3]), mean(x[, 2] <= x[, 4]), mean(x[, 1] <= x[, 2]))) }
writeLines(lines, file.path(OUT, "burden_ordering_summary.txt")); cat(lines, sep = "\n")
cat("\nN_q with draw-bootstrap 90% intervals:\n"); print(bs[order(bs$state, bs$model), c("pair", "N90", "N90_lo", "N90_hi", "N95", "N95_lo", "N95_hi")], row.names = FALSE)

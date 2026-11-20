################################################################
# hidden_context_sim.R
#
# A two-context simulation: global recoverability can conceal a
# context in which the representation carries no information.
#
# Context S in {A, B}, unobserved by the estimator; P(S = A) = pi_A.
# Z ~ N_50(0, Sigma), Sigma_jk = 0.3^|j-k|, the same in both contexts
# (the context is not visible in the representation).
#   S = A:  logit P(Y = 1 | Z) = b0 + Z beta, beta = (1.4, -1.1, 0.8, 0, ..., 0)
#   S = B:  Y independent of Z, P(Y = 1) = P(Y = 1 | S = A)
# so H(Y) is the same in both contexts, I(Z;Y | S = A) > 0,
# I(Z;Y | S = B) = 0, and by the chain rule, since Z is independent
# of S, I(Z;Y) <= I(Z;Y | S) = pi_A I(Z;Y | S = A).
#
# Truth by Monte Carlo (10^6 draws). For pi_A in {0.25, 0.5, 0.75}
# with n = 400 patients in all (n_A = pi_A n), R = 50 replicates:
# PathMI on the pooled cohort, on the A patients and on the B patients;
# permutation null at each of the sample sizes (50 refits each).
#
# Usage (project root, R 4.4.3): Rscript TRIDENT/hidden_context_sim.R
# Output: TRIDENT/simulation_context/{context_draws,context_null,context_truth,context_summary}.csv
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(parallel) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI  <- normalizePath(dirname(this_file), winslash = "/")
ROOT <- normalizePath(file.path(TRI, ".."), winslash = "/")
setwd(ROOT); source("PathMI_SUBHO_logistic.R")
OUT <- file.path(TRI, "simulation_context"); dir.create(OUT, showWarnings = FALSE); LOG <- file.path(OUT, "hidden_context_sim.log")
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }

N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "12")); P <- 50L; RHO <- 0.3; N_TOT <- 400L; R <- 50L; N_PERM <- 50L
PI_GRID <- c(0.25, 0.5, 0.75); B0 <- -0.8; BETA <- c(1.4, -1.1, 0.8, rep(0, P - 3))
NMC <- 2000L; BURN <- 500L; THIN <- 10L
SIG <- RHO^abs(outer(seq_len(P), seq_len(P), "-")); L <- chol(SIG)
rZ <- function(n) matrix(rnorm(n * P), n, P) %*% L
hb <- function(p) { p <- pmin(pmax(p, 1e-300), 1 - 1e-16); -(p * log(p) + (1 - p) * log(1 - p)) }

# ---- truth -------------------------------------------------------------------
set.seed(1); Zmc <- rZ(1e6); pA_z <- plogis(B0 + as.vector(Zmc %*% BETA)); pA <- mean(pA_z)
I_A <- hb(pA) - mean(hb(pA_z))
truth <- do.call(rbind, lapply(PI_GRID, function(pa) data.frame(pi_A = pa, prevalence = pA, H_Y = hb(pA), I_A = I_A, I_B = 0,
  I_conditional = pa * I_A, I_pooled = hb(pA) - mean(hb(pa * pA_z + (1 - pa) * pA)))))
rm(Zmc, pA_z); write.csv(truth, file.path(OUT, "context_truth.csv"), row.names = FALSE)
logf("truth: prevalence %.3f, H(Y) %.3f, I(Z;Y|A) = %.4f; pooled I(Z;Y) = %s", pA, hb(pA), I_A, paste(sprintf("%.4f at pi_A = %.2f", truth$I_pooled, truth$pi_A), collapse = ", "))

simulate <- function(pa, seed) {
  set.seed(seed); nA <- round(pa * N_TOT); Z <- rZ(N_TOT); S <- rep(c("A", "B"), c(nA, N_TOT - nA))
  y <- ifelse(S == "A", rbinom(N_TOT, 1, plogis(B0 + as.vector(Z %*% BETA))), rbinom(N_TOT, 1, pA))
  list(Z = scale(Z), y = y, S = S)
}
fit <- function(Z, y, seed) { set.seed(seed); subho_MI_logistic(Z = Z, y = y, v_beta = 1, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)$MI }
job <- function(pa, r) {
  d <- simulate(pa, 1000L * which(PI_GRID == pa) + r); a <- d$S == "A"
  data.frame(pi_A = pa, rep = r, n_A = sum(a), n_B = sum(!a),
             I_pooled = fit(d$Z, d$y, 10L * r + 1L), I_A = fit(scale(d$Z[a, ]), d$y[a], 10L * r + 2L), I_B = fit(scale(d$Z[!a, ]), d$y[!a], 10L * r + 3L))
}
null_job <- function(n, b) { set.seed(50000L + n + b); Z <- scale(rZ(n)); y <- rbinom(n, 1, pA); data.frame(n = n, b = b, I_null = fit(Z, y, 60000L + b)) }

jobs <- expand.grid(pa = PI_GRID, r = seq_len(R)); n_sizes <- sort(unique(c(N_TOT, round(PI_GRID * N_TOT)))); njobs <- expand.grid(n = n_sizes, b = seq_len(N_PERM))
cl <- makeCluster(N_WORKERS); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(MASS); library(BayesLogit) }))
clusterExport(cl, c("jobs", "njobs", "job", "null_job", "simulate", "fit", "rZ", "L", "P", "N_TOT", "PI_GRID", "B0", "BETA", "pA", "NMC", "BURN", "THIN",
                    "sagnik.horseshoe.logistic", "lpy.x.logistic", "selection.logistic", "subho_MI_logistic", "sigmoid", "clip_eta"))
t0 <- Sys.time()
draws <- do.call(rbind, parLapply(cl, seq_len(nrow(jobs)), function(j) job(jobs$pa[j], jobs$r[j])))
nulls <- do.call(rbind, parLapply(cl, seq_len(nrow(njobs)), function(j) null_job(njobs$n[j], njobs$b[j])))
write.csv(draws, file.path(OUT, "context_draws.csv"), row.names = FALSE); write.csv(nulls, file.path(OUT, "context_null.csv"), row.names = FALSE)
nm <- tapply(nulls$I_null, nulls$n, mean); ns <- tapply(nulls$I_null, nulls$n, sd)
summ <- do.call(rbind, lapply(split(draws, draws$pi_A), function(d) { pa <- d$pi_A[1]; nA <- d$n_A[1]; nB <- d$n_B[1]
  data.frame(pi_A = pa, n_A = nA, n_B = nB,
    pooled_mean = mean(d$I_pooled), pooled_sd = sd(d$I_pooled), pooled_null = nm[[as.character(N_TOT)]], pooled_above = mean(d$I_pooled > nm[[as.character(N_TOT)]] + 2 * ns[[as.character(N_TOT)]]),
    A_mean = mean(d$I_A), A_sd = sd(d$I_A), A_null = nm[[as.character(nA)]], A_above = mean(d$I_A > nm[[as.character(nA)]] + 2 * ns[[as.character(nA)]]),
    B_mean = mean(d$I_B), B_sd = sd(d$I_B), B_null = nm[[as.character(nB)]], B_above = mean(d$I_B > nm[[as.character(nB)]] + 2 * ns[[as.character(nB)]]),
    conditional_mean = mean((nA * d$I_A + nB * d$I_B) / (nA + nB))) }))
summ <- merge(summ, truth[, c("pi_A", "I_A", "I_pooled", "I_conditional")], by = "pi_A")
write.csv(summ, file.path(OUT, "context_summary.csv"), row.names = FALSE)
logf("=== done in %.1f min ===", as.numeric(difftime(Sys.time(), t0, units = "mins")))
for (i in seq_len(nrow(summ))) with(summ[i, ], logf("pi_A %.2f: pooled %.4f (truth %.4f, above null in %.0f%%) | A %.4f (truth %.4f, %.0f%%) | B %.4f (truth 0, null %.4f, above null in %.0f%%) | conditional %.4f (truth %.4f)",
  pi_A, pooled_mean, I_pooled, 100 * pooled_above, A_mean, I_A, 100 * A_above, B_mean, B_null, 100 * B_above, conditional_mean, I_conditional))

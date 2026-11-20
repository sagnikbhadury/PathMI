################################################################
# patient_collapse.R
#
# Sensitivity of the recorded PathMI estimates to the unit of analysis.
# The rosters are slides; seven patients per cohort contribute two to
# eight slides with concordant labels, which the estimator treats as
# independent observations. Here each patient's slide-level pooled
# vectors are averaged into one row, the projection is refitted, and
# PathMI is refitted (same settings and seed as the pipeline), for the
# twelve Recoverability Record arms.
#
# Usage (project root, R 4.4.3): Rscript TRIDENT/patient_collapse.R
# Output: TRIDENT/record/patient_collapse.csv
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(parallel) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); ROOT <- normalizePath(file.path(TRI, ".."), winslash = "/")
setwd(ROOT); source("PathMI_SUBHO_logistic.R"); source(file.path(TRI, "infold_functions.R"))
POOL <- "M:/bhadury/mi_results/pool_pass"
code <- c("UNI v2" = "uni_v2", "Virchow2" = "virchow2", "H-optimus-1" = "hoptimus1", "CONCH v1.5" = "conch_v15")
rec <- read.csv(file.path(TRI, "record", "recoverability_record.csv"), stringsAsFactors = FALSE); rec <- rec[rec$record_type == "measurement", ]
one <- function(i) {
  r <- rec[i, ]; pp <- readRDS(file.path(POOL, sprintf("%s_%s", tolower(r$cohort), code[[r$model]]), "pooled.rds"))
  X <- pp$pooled[[r$pooling]]; y <- as.numeric(pp$y[[r$state]][pp$roster]); g <- substr(pp$roster, 1, 12)
  Xp <- rowsum(X, g, reorder = FALSE) / as.vector(table(g)[unique(g)]); yp <- y[match(unique(g), g)]
  Z <- unclass(apply_projection(fit_projection(Xp, r$projection, 50L), Xp))
  set.seed(42L); f <- subho_MI_logistic(Z = Z, y = yp, v_beta = 1, method.weight = "beta.binomial", nmc = 2000L, burn = 500L, thin = 10L)
  data.frame(cohort = r$cohort, state = r$state, model = r$model, n_slides = length(y), n_patients = length(yp), I_slides = r$I_hat, I_patients = f$MI, k_patients = f$n.selected)
}
cl <- makeCluster(6L); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(MASS); library(BayesLogit) }))
clusterExport(cl, c("rec", "one", "POOL", "code", "fit_projection", "apply_projection", "EIG_TOL", "sagnik.horseshoe.logistic", "lpy.x.logistic", "selection.logistic", "subho_MI_logistic", "sigmoid", "clip_eta"))
out <- do.call(rbind, parLapply(cl, seq_len(nrow(rec)), one)); out$diff <- out$I_patients - out$I_slides
write.csv(out, file.path(TRI, "record", "patient_collapse.csv"), row.names = FALSE)
print(out, digits = 3); cat(sprintf("patient minus slide: median %+.4f, range %+.4f to %+.4f\n", median(out$diff), min(out$diff), max(out$diff)))

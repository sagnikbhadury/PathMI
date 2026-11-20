################################################################
# pl_grid.R
#
# Patient-level rerun of the 192-cell recoverability grid (19 Nov 2026,
# QC finding: the original rosters were slides, and seven patients per
# cohort contributed two to eight slides that the estimator counted as
# independent observations).
#
# Starts from the per-slide pooled matrices written by pool_pass.R
# (pipeline-identical pooling; the d = 50 projections rebuilt from them
# reproduce the stored Z_scaled exactly). For each cohort x model:
#   1. MODE = "patient": average the slide-level pooled vectors of each
#      patient (one row per patient; labels are concordant within
#      patient); tile counts are summed; patient covariates taken once.
#      MODE = "slide": keep slide rows (validation: the observed fits
#      must reproduce the stored pipeline estimates exactly).
#   2. The four projections (projection_functions.R) for each of the
#      four pooling rules: 16 representations.
#   3. PathMI observed fit (set.seed(42), as the pipeline), and 200
#      permutation refits. The permutations are drawn once per cohort x
#      model (set.seed(42)) and SHARED by all 16 representations and
#      both prostate targets, so the max-T null is a joint permutation
#      (QC finding E1); each permutation fit uses set.seed(100000 + b).
#      Raw estimates and cond(H) are kept for every fit (QC finding E5).
#   4. Selection per pooling rule (and target) with the max-T p value.
#   5. Incremental fits on the selected arms: covariates only,
#      embedding only and combined; CRC: log tile count + AJCC stage;
#      PRAD: log tile count + Gleason (C0) and C0 + purity + batch
#      (C_all). The incremental block uses its own seed stream
#      (set.seed(4242) once), so the embedding-only refit is a
#      replicate refit with an independent random-number stream.
#
# Usage (project root, R 4.4.3):
#   Rscript TRIDENT/patient_level/pl_grid.R [patient|slide] [cohort_model,...|all]
# Output: TRIDENT/patient_level/<mode>/<cohort>_<model>/ (results.rds,
#   mi_summary.csv, best.csv, incremental.csv) and a progress log.
################################################################
suppressPackageStartupMessages({ library(MASS); library(BayesLogit); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) args[1] else "patient"; stopifnot(MODE %in% c("patient", "slide"))
WHICH <- if (length(args) >= 2) args[2] else "all"
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
PLD <- normalizePath(dirname(this_file), winslash = "/"); TRI <- dirname(PLD); ROOT <- dirname(TRI)
setwd(ROOT); source("PathMI_SUBHO_logistic.R"); source(file.path(TRI, "projection_functions.R"))
OUT <- file.path(PLD, MODE); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "pl_grid.log")
logf <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "%H:%M:%S "), m, "\n", file = LOG, append = TRUE) }

N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "24")); N_PERM <- as.integer(Sys.getenv("N_PERM", "200"))
V_BETA <- 1.0; NMC <- 2000L; BURN <- 500L; THIN <- 10L; SEED <- 42L; D <- 50L
POOL <- "M:/bhadury/mi_results/pool_pass"
META <- list(crc = "M:/kulkarnr/slideflow_test/slideflow_crc/annotations.csv", prad = "M:/shared_data/TCGA_PRAD/orig_annotations.csv")
MODELS <- c("uni_v2", "virchow2", "hoptimus1", "conch_v15")
STAGE <- c("Stage IA" = 1, "Stage I" = 1, "Stage IIA" = 2, "Stage IIB" = 2, "Stage IIC" = 2, "Stage II" = 2,
           "Stage IIIA" = 3, "Stage IIIB" = 3, "Stage IIIC" = 3, "Stage III" = 3, "Stage IVA" = 4, "Stage IVB" = 4, "Stage IV" = 4)
runs <- expand.grid(model = MODELS, cohort = c("crc", "prad"), stringsAsFactors = FALSE)[, 2:1]
runs$tag <- paste(runs$cohort, runs$model, sep = "_")
if (WHICH != "all") runs <- runs[runs$tag %in% strsplit(WHICH, ",")[[1]], ]
logf("=== pl_grid | mode %s | %s | %d permutations | %d workers ===", MODE, paste(runs$tag, collapse = " "), N_PERM, N_WORKERS)

# ---- build every representation -------------------------------------------------
cells <- list(); meta_runs <- list()
for (i in seq_len(nrow(runs))) {
  co <- runs$cohort[i]; mo <- runs$model[i]; tag <- runs$tag[i]
  pp <- readRDS(file.path(POOL, tag, "pooled.rds"))
  slide <- pp$roster; pat <- substr(slide, 1, 12); tc <- setNames(pp$tiles$n_tiles, pp$tiles$slide)[slide]
  unit <- if (MODE == "patient") pat else slide
  units <- unique(unit)
  collapse <- function(X) { S <- rowsum(X, unit, reorder = FALSE); S / as.vector(table(unit)[rownames(S)]) }
  raw <- lapply(pp$pooled, function(X) { if (MODE == "patient") collapse(X) else X })
  ylist <- lapply(pp$y, function(v) { v <- v[slide]; setNames(as.numeric(v[match(units, unit)]), units) })
  disc <- sapply(pp$y, function(v) sum(tapply(v[slide], pat, function(z) length(unique(z)) > 1)))
  tiles <- setNames(as.numeric(tapply(tc, unit, sum)[units]), units)
  m <- read.csv(META[[co]], stringsAsFactors = FALSE); m <- m[!is.na(m$slide) & !duplicated(m$slide), ]; rownames(m) <- m$slide
  mu <- m[slide[match(units, unit)], ]                      # one annotation row per unit (patient attributes are shared)
  if (co == "crc") {
    st <- STAGE[mu$ajcc_pathologic_stage]; st[is.na(st)] <- median(st, na.rm = TRUE)
    C0 <- cbind(log_tile_count_z = as.numeric(scale(log1p(tiles))), stage_ordinal_z = as.numeric(scale(st))); Csets <- list(C0 = C0)
  } else {
    gle <- mu$Reviewed_Gleason_sum; gle[is.na(gle)] <- median(gle, na.rm = TRUE)
    pur <- suppressWarnings(as.numeric(mu$ABSOLUTE_purity)); pur[is.na(pur)] <- median(pur, na.rm = TRUE)
    bt <- as.character(mu$Batch); keep <- names(which(table(bt) >= 8)); bt[!bt %in% keep] <- "other"
    ref <- names(which.max(table(bt))); lv <- setdiff(sort(unique(bt)), ref)
    Bm <- sapply(lv, function(l) as.numeric(bt == l)); colnames(Bm) <- paste0("batch_", lv)
    C0 <- cbind(log_tile_count_z = as.numeric(scale(log1p(tiles))), Gleason_sum_z = as.numeric(scale(gle)))
    Csets <- list(C0 = C0, C_all = cbind(C0, purity_z = as.numeric(scale(pur)), scale(Bm)))
  }
  set.seed(SEED); perms <- sapply(seq_len(N_PERM), function(b) sample.int(length(units)))   # shared label permutations
  meta_runs[[tag]] <- list(cohort = co, model = mo, units = units, ylist = ylist, Csets = Csets, perms = perms,
                           n_slides = length(slide), n_units = length(units), discordant = disc)
  for (pl in c("mean", "norm", "moment", "attention")) for (dc in c("pca", "cosine", "rbf", "poly")) {
    pr <- project(raw[[pl]], dc, D); key <- paste0(dc, "_", pl)
    cells[[paste(tag, key, sep = "|")]] <- list(tag = tag, key = key, Z = unclass(pr$Z_scaled), eig_cond = pr$eig_cond)
  }
  logf("[%s] %d slides -> %d %s rows | label-discordant patients: %s | covariate sets %s", tag, length(slide), length(units), MODE,
       paste(names(disc), disc, collapse = ", "), paste(names(Csets), sapply(Csets, ncol), collapse = ", "))
}

# ---- jobs: one per (cell, target, block of permutations); block 0 = observed ------------------
BLOCK <- 25L
jobs <- do.call(rbind, lapply(names(cells), function(cid) { tg <- names(meta_runs[[cells[[cid]]$tag]]$ylist)
  expand.grid(cid = cid, target = tg, block = 0:(N_PERM / BLOCK), stringsAsFactors = FALSE) }))
fit1 <- function(Z, y, seed) { set.seed(seed)
  r <- subho_MI_logistic(Z = Z, y = y, v_beta = V_BETA, method.weight = "beta.binomial", nmc = NMC, burn = BURN, thin = THIN)
  c(MI = r$MI, MI_raw = r$MI.raw, k = r$n.selected, cond = r$cond.H.best) }
job <- function(j) {
  cid <- jobs$cid[j]; tg <- jobs$target[j]; bl <- jobs$block[j]; c0 <- cells[[cid]]; mr <- meta_runs[[c0$tag]]; y <- mr$ylist[[tg]]
  if (bl == 0L) { f <- fit1(c0$Z, y, SEED); return(data.frame(cid = cid, target = tg, b = 0L, MI = f[["MI"]], MI_raw = f[["MI_raw"]], k = f[["k"]], cond = f[["cond"]])) }
  bs <- ((bl - 1L) * BLOCK + 1L):(bl * BLOCK)
  do.call(rbind, lapply(bs, function(b) { f <- fit1(c0$Z, y[mr$perms[, b]], 100000L + b)
    data.frame(cid = cid, target = tg, b = b, MI = f[["MI"]], MI_raw = f[["MI_raw"]], k = f[["k"]], cond = f[["cond"]]) }))
}
CK <- file.path(OUT, "checkpoint"); dir.create(CK, showWarnings = FALSE)
set.seed(11); jobs <- jobs[sample.int(nrow(jobs)), ]; jobs$chunk <- ceiling(seq_len(nrow(jobs)) / (2L * N_WORKERS))
logf("%d representations, %d jobs in %d chunks", length(cells), nrow(jobs), max(jobs$chunk))
cl <- makeCluster(N_WORKERS); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(MASS); library(BayesLogit) }))
clusterExport(cl, c("cells", "meta_runs", "jobs", "job", "fit1", "BLOCK", "SEED", "V_BETA", "NMC", "BURN", "THIN",
                    "sagnik.horseshoe.logistic", "lpy.x.logistic", "selection.logistic", "subho_MI_logistic", "sigmoid", "clip_eta"))
t0 <- Sys.time()
for (ch in sort(unique(jobs$chunk))) {
  f <- file.path(CK, sprintf("chunk_%04d.rds", ch)); if (file.exists(f)) next
  idx <- which(jobs$chunk == ch)
  o <- parLapply(cl, idx, function(j) tryCatch(job(j), error = function(e) data.frame(cid = jobs$cid[j], target = jobs$target[j], b = NA, err = conditionMessage(e))))
  saveRDS(o, f)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins")); logf("  chunk %d / %d | %.1f min | ETA %.1f min", ch, max(jobs$chunk), el, el / ch * (max(jobs$chunk) - ch))
}
draws <- do.call(rbind, lapply(list.files(CK, full.names = TRUE), function(f) do.call(rbind, lapply(readRDS(f), function(d) { if (is.null(d$err)) d else { logf("[ERR] %s %s: %s", d$cid[1], d$target[1], d$err[1]); NULL } }))))
saveRDS(draws, file.path(OUT, "all_draws.rds"))

# ---- assemble per cohort x model ----------------------------------------------------
for (tag in names(meta_runs)) {
  mr <- meta_runs[[tag]]; dd <- draws[startsWith(draws$cid, paste0(tag, "|")), ]; dir.create(file.path(OUT, tag), showWarnings = FALSE)
  res <- list(); rows <- list()
  for (cid in unique(dd$cid)) { key <- sub(".*\\|", "", cid); det <- list()
    for (tg in names(mr$ylist)) { d <- dd[dd$cid == cid & dd$target == tg, ]; o <- d[d$b == 0, ]; p <- d[d$b > 0, ]; p <- p[order(p$b), ]
      pv <- (1 + sum(p$MI >= o$MI)) / (1 + nrow(p))
      det[[tg]] <- list(target = tg, MI = o$MI, MI.raw = o$MI_raw, n.selected = o$k, cond.H = o$cond, MI.null = p$MI, p.value = pv)
      rows[[length(rows) + 1]] <- data.frame(tag = tag, key = key, decomp = sub("_.*", "", key), agg_method = sub(".*_", "", key), target = tg,
        n = length(mr$units), MI_obs = o$MI, MI_raw = o$MI_raw, k_star = o$k, cond_H = o$cond, null_mean = mean(p$MI), null_sd = sd(p$MI), p_value = pv,
        eig_cond = cells[[cid]]$eig_cond) }
    res[[key]] <- list(Z_scaled = cells[[cid]]$Z, y_list = mr$ylist, units = mr$units, result = list(details = det)) }
  saveRDS(res, file.path(OUT, tag, "results.rds"))
  summ <- do.call(rbind, rows); write.csv(summ, file.path(OUT, tag, "mi_summary.csv"), row.names = FALSE)
  # selection per pooling (and target) with joint-permutation max-T
  best <- do.call(rbind, lapply(split(summ, list(summ$agg_method, summ$target), drop = TRUE), function(s) {
    s <- s[s$eig_cond < 1e10, ]; b <- s[which.max(s$MI_obs), ]
    nm <- sapply(s$key, function(k) res[[k]]$result$details[[b$target]]$MI.null); Tn <- apply(nm, 1, max)
    data.frame(b, p_maxT = (1 + sum(Tn >= b$MI_obs)) / (1 + length(Tn)), n_candidates = nrow(s)) }))
  write.csv(best, file.path(OUT, tag, "best.csv"), row.names = FALSE)
  logf("[%s] assembled: %d cells | best: %s", tag, nrow(summ), paste(sprintf("%s/%s %s %.4f", best$target, best$agg_method, best$key, best$MI_obs), collapse = "; "))
}
saveRDS(list(meta_runs = lapply(meta_runs, function(m) m[c("cohort", "model", "units", "ylist", "Csets", "n_slides", "n_units", "discordant")])), file.path(OUT, "meta_runs.rds"))
logf("=== grid done in %.1f min ===", as.numeric(difftime(Sys.time(), t0, units = "mins")))

################################################################
# results_four_backbones_numbers.R
#
# Prints every quantity the four-model Results section quotes,
# computed from the three merged tables (192 / 48 / 48 rows) and the
# per-run rds objects (rosters, prevalences). Companion to
# figures_four_backbones.R (Figure 1) and figures_results_four.R
# (Figure 2). Usage: Rscript TRIDENT/results_four_backbones_numbers.R
################################################################
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr) })
TRI <- "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/entropy/Research - DOE Horseshoe/TRIDENT"
BB <- c("UNI v2", "Virchow2", "H-optimus-1", "CONCH v1.5")
harm <- function(e) case_when(grepl("^UNI", e) ~ BB[1], grepl("^VIRCHOW", e) ~ BB[2], grepl("^H-OPT", e) ~ BB[3], grepl("^CONCH", e) ~ BB[4], TRUE ~ NA_character_)
g <- read_excel(file.path(TRI, "mi summary table.xlsx")) %>%
  mutate(bb = factor(harm(Embeddings), BB), task = paste0(Disease, "/", target), decomp = tolower(decomp),
         z = (MI_obs - null_mean) / null_sd, sig = p_value < 0.05)
best <- read_excel(file.path(TRI, "best kpca summary.xlsx")) %>%
  mutate(bb = factor(harm(Embeddings), BB), task = paste0(Disease, "/", target), kern = sub("_.*", "", key))
inc <- read_excel(file.path(TRI, "incremental mi summary.xlsx")) %>%
  mutate(bb = factor(harm(Embeddings), BB), task = paste0(Disease, "/", target))
stopifnot(nrow(g) == 192, nrow(best) == 48, nrow(inc) == 48)

cat("=== 1. Rosters and prevalences (from the per-run rds) ===\n")
runs <- list(c("CRC", "uni_v2"), c("CRC", "virchow2"), c("CRC", "hoptimus1"), c("CRC", "conch_v15"),
             c("PRAD", "uni_v2"), c("PRAD", "virchow2"), c("PRAD", "hoptimus1"), c("PRAD", "conch_v15"))
rost <- list()
for (r in runs) {
  f <- if (r[1] == "CRC") sprintf("crc_kpca_subho_%s/kpca_subho_results_crc_%s.rds", r[2], r[2]) else sprintf("TCGA_PRAD_kpca_subho_%s/kpca_subho_results_prad_%s.rds", r[2], r[2])
  a <- readRDS(file.path(TRI, f))$pca_mean
  ids <- rownames(a$Z_scaled); rost[[paste(r, collapse = "_")]] <- ids
  if (r[1] == "CRC") cat(sprintf("%-4s %-10s n=%d  MSI-H prev %.3f  H(Y) %.3f\n", r[1], r[2], length(ids), mean(a$y), -mean(a$y) * log(mean(a$y)) - (1 - mean(a$y)) * log(1 - mean(a$y))))
  else for (t in c("ERG", "ETS")) { p <- mean(a$y_list[[t]]); cat(sprintf("%-4s %-10s n=%d  %s prev %.3f  H(Y) %.3f\n", r[1], r[2], length(ids), t, p, -p * log(p) - (1 - p) * log(1 - p))) }
}
cat("common roster CRC (all four):", length(Reduce(intersect, rost[1:4])), "| PRAD:", length(Reduce(intersect, rost[5:8])), "\n")

cat("\n=== 2. Signal: cells above null, null mean/sd, peak, by task x model ===\n")
s <- g %>% group_by(task, bb) %>% summarise(sig = sum(sig), null_mean = mean(null_mean), null_sd = mean(null_sd), peak = max(MI_obs),
  zmin_sig = min(z[sig]), zmax = max(z), pfloor = sum(p_value < 0.0051), .groups = "drop")
print(as.data.frame(s), digits = 3)
cat(sprintf("total significant: %d / 192 | at the p floor: %d\n", sum(g$sig), sum(g$p_value < 0.0051)))
cat("null mean/sd by task (over 64 cells):\n"); print(g %>% group_by(task) %>% summarise(m = mean(null_mean), s = mean(null_sd), hi = m + 2 * s) %>% as.data.frame(), digits = 3)
cat("z range among significant cells:", round(range(g$z[g$sig]), 1), "\n")
cat("peak per task (bb, decomp, pool):\n"); print(g %>% group_by(task) %>% slice_max(MI_obs, n = 1) %>% select(task, bb, decomp, agg_method, MI_obs) %>% as.data.frame(), digits = 4)

cat("\n=== 3. Conditioning: kappa of selected kernels and r(MI, kappa) ===\n")
print(best %>% group_by(task, bb) %>% summarise(kappa_mean = mean(eig_cond), kappa_min = min(eig_cond), kappa_max = max(eig_cond), MI_mean = mean(MI_obs), .groups = "drop") %>% as.data.frame(), digits = 3)
cat("r(MI, kappa) across the 16 selected cells per task (four models):\n")
print(best %>% group_by(task) %>% summarise(r = cor(MI_obs, eig_cond), r_spearman = cor(MI_obs, eig_cond, method = "spearman")) %>% as.data.frame(), digits = 3)
cat("r(MI, kappa) within task, first two models only (8 cells), for comparison with the earlier text:\n")
print(best %>% filter(bb %in% BB[1:2]) %>% group_by(task) %>% summarise(r = cor(MI_obs, eig_cond)) %>% as.data.frame(), digits = 3)
cat("r(MI, kappa) within task excluding CONCH (12 cells):\n")
print(best %>% filter(bb != "CONCH v1.5") %>% group_by(task) %>% summarise(r = cor(MI_obs, eig_cond)) %>% as.data.frame(), digits = 3)
cat("kappa on the full grid by decomposition x model (median over tasks and poolings, from best table only where selected):\n")

cat("\n=== 4. Variance partition (additive main effects), four models ===\n")
part <- function(d, f) { gm <- mean(d$MI_obs); ss <- d %>% group_by(.data[[f]]) %>% summarise(k = n(), m = mean(MI_obs), .groups = "drop"); 100 * sum(ss$k * (ss$m - gm)^2) / sum((d$MI_obs - gm)^2) }
vp <- g %>% group_by(task) %>% group_modify(~ tibble(model = part(.x, "bb"), decomposition = part(.x, "decomp"), pooling = part(.x, "agg_method"))) %>% mutate(residual = 100 - model - decomposition - pooling)
print(as.data.frame(vp), digits = 3)
cat("within-model range of the estimate over the 16 cells (max - min), by task x model:\n")
print(g %>% group_by(task, bb) %>% summarise(range = max(MI_obs) - min(MI_obs), .groups = "drop") %>% pivot_wider(names_from = bb, values_from = range) %>% as.data.frame(), digits = 3)
cat("attention best pooling (mean over decompositions) in how many of the 12 task x model combinations:\n")
pm <- g %>% group_by(task, bb, agg_method) %>% summarise(m = mean(MI_obs), .groups = "drop") %>% group_by(task, bb) %>% slice_max(m, n = 1)
print(table(pm$agg_method))

cat("\n=== 5. Kernel selection (48 arms) ===\n")
print(table(best$kern, best$bb)); print(table(best$kern, best$task))
cat(sprintf("maxT p at floor (0.005): %d / 48 | maxT p < 0.05: %d / 48\n", sum(best$p_maxT < 0.0051), sum(best$p_maxT < 0.05)))
print(best %>% filter(p_maxT >= 0.0051) %>% select(task, bb, agg_method, label, MI_obs, p_maxT) %>% as.data.frame(), digits = 3)
cat("polynomial: mean MI rank among the four decompositions per (task, model, pooling), and cells non-significant by decomposition:\n")
print(table(g$decomp[!g$sig], g$bb[!g$sig])); cat("non-significant cells total:", sum(!g$sig), "\n")
print(g %>% filter(decomp == "poly") %>% group_by(task, bb) %>% summarise(mi_min = min(MI_obs), mi_max = max(MI_obs), n_sig = sum(sig), .groups = "drop") %>% as.data.frame(), digits = 3)

cat("\n=== 6. Incremental decomposition, CRC (16 rows) ===\n")
ci <- inc %>% filter(Disease == "CRC") %>% mutate(share = covar_only / combined, synergy = combined - covar_only - emb_only) %>% arrange(bb, match(agg_method, c("mean", "norm", "moment", "attention")))
print(as.data.frame(ci %>% select(bb, agg_method, covar_only, emb_only, combined, incremental_subho_mi, share, synergy)), digits = 4)
cat("covariate share range by model:\n"); print(ci %>% group_by(bb) %>% summarise(lo = min(share), hi = max(share), covar = mean(covar_only)) %>% as.data.frame(), digits = 3)
cat("PRAD covar_only == 0 in", sum(inc$covar_only[inc$Disease == "PRAD"] == 0), "of", sum(inc$Disease == "PRAD"), "rows\n")

cat("\n=== 7. Reproducibility: replicate refit (CRC) and null-covariate identity (PRAD) ===\n")
d1 <- inc %>% filter(Disease == "CRC") %>% inner_join(best %>% select(Disease, bb, target, agg_method, MI_obs), by = c("Disease", "bb", "target", "agg_method")) %>%
  transmute(task, bb, agg_method, d = emb_only - MI_obs, src = "replicate refit")
d2 <- inc %>% filter(covar_only == 0) %>% transmute(task, bb, agg_method, d = combined - emb_only, src = "null covariate block")
d3 <- inc %>% filter(Disease == "CRC") %>% transmute(task, bb, agg_method, d = combined - emb_only, src = "CRC covariate drop (combined - emb_only, covariates carry > 0)")
disc <- bind_rows(d1, d2)
cat("nonzero discrepancies (|d| > 1e-9):\n"); print(disc %>% filter(abs(d) > 1e-9) %>% arrange(d) %>% as.data.frame(), digits = 3)
cat(sprintf("replicate refit: exact in %d / %d; null covariate: exact in %d / %d\n", sum(abs(d1$d) < 1e-9), nrow(d1), sum(abs(d2$d) < 1e-9), nrow(d2)))
cat(sprintf("bound over four models: %.4f | over the first two models: %.4f | over the top tier (UNI v2, H-optimus-1): %.4f | excluding CONCH: %.4f\n",
            max(abs(disc$d)), max(abs(disc$d[disc$bb %in% BB[1:2]])), max(abs(disc$d[disc$bb %in% BB[c(1, 3)]])), max(abs(disc$d[disc$bb != "CONCH v1.5"]))))
cat("CRC: combined below emb_only (covariates carry information but the estimate falls):\n"); print(d3 %>% filter(d < 0) %>% as.data.frame(), digits = 3)

cat("\n=== 8. Backbone contrasts within task: selected-arm means and paired sign counts vs UNI v2 (16 grid cells) ===\n")
w <- g %>% select(task, decomp, agg_method, bb, MI_obs) %>% pivot_wider(names_from = bb, values_from = MI_obs)
for (tk in unique(w$task)) { d <- w[w$task == tk, ]
  for (b in BB[-1]) cat(sprintf("%-9s %-12s higher than UNI v2 in %2d / 16 | median diff %+.4f | sign test p = %.3f\n", tk, b, sum(d[[b]] > d[["UNI v2"]]), median(d[[b]] - d[["UNI v2"]]), binom.test(sum(d[[b]] > d[["UNI v2"]]), 16)$p.value)) }
cat("selected-arm means by task x model:\n"); print(round(with(best, tapply(MI_obs, list(task, bb), mean)), 4))
cat("selected-arm peak as percent of H(Y) [H = 0.401, 0.688, 0.682]:\n")
print(round(100 * with(best, tapply(MI_obs, list(task, bb), max)) / c(0.401, 0.688, 0.682), 1))

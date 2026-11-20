################################################################
# Is there one configuration (backbone x decomposition x pooling)
# to use in general?  Tests on the three merged tables.
################################################################
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr) })
TRI <- "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/entropy/Research - DOE Horseshoe/TRIDENT"

harm <- function(e) case_when(grepl("^UNI", e) ~ "UNIv2", grepl("^VIRCHOW", e) ~ "Virchow2",
                              grepl("^H-OPT", e) ~ "Hopt1", grepl("^CONCH", e) ~ "CONCH", TRUE ~ NA_character_)
g <- read_excel(file.path(TRI, "mi summary table.xlsx")) %>%
  mutate(bb = harm(Embeddings), task = paste0(Disease, "/", target), decomp = tolower(decomp),
         z = (MI_obs - null_mean) / null_sd, sig = p_value < 0.05,
         block = paste(task, decomp, agg_method))
stopifnot(nrow(g) == 192, !anyNA(g$bb))
best <- read_excel(file.path(TRI, "best kpca summary.xlsx")) %>%
  mutate(bb = harm(Embeddings), task = paste0(Disease, "/", target), kern = sub("_.*", "", key))
inc  <- read_excel(file.path(TRI, "incremental mi summary.xlsx")) %>%
  mutate(bb = harm(Embeddings), task = paste0(Disease, "/", target))
BB <- c("UNIv2", "Virchow2", "Hopt1", "CONCH")

cat("=== 1. Cells above the permutation null (p < 0.05), by task and backbone (of 16) ===\n")
print(with(g, tapply(sig, list(task, bb), sum))[, BB])

cat("\n=== 2. Variance partition per task, additive main effects (percent), 4 backbones ===\n")
part <- function(d, f) { gm <- mean(d$MI_obs); s <- d %>% group_by(.data[[f]]) %>%
  summarise(k = n(), m = mean(MI_obs), .groups = "drop"); 100 * sum(s$k * (s$m - gm)^2) / sum((d$MI_obs - gm)^2) }
vp <- g %>% group_by(task) %>% group_modify(~ tibble(backbone = part(.x, "bb"), decomposition = part(.x, "decomp"),
  pooling = part(.x, "agg_method"))) %>% mutate(residual = 100 - backbone - decomposition - pooling)
print(as.data.frame(vp), digits = 3)

cat("\n=== 3. Does the backbone ranking depend on the task?  ANOVA, one obs per cell ===\n")
cat("Model: MI ~ (task + bb + decomp + pool)^2; residual = 3-way and higher interactions\n")
g$pool <- g$agg_method
fit <- lm(MI_obs ~ (task + bb + decomp + pool)^2, data = g)
a <- anova(fit); print(round(a[, c("Df", "Sum Sq", "F value", "Pr(>F)")], 4))

cat("\n=== 4. Friedman tests: each factor, blocks = the other three factors crossed (48 blocks) ===\n")
fr <- function(f) { d <- g %>% mutate(blk = apply(g[, setdiff(c("task","bb","decomp","pool"), f)], 1, paste, collapse = "|"))
  ft <- friedman.test(d$MI_obs, factor(d[[f]]), factor(d$blk))
  mr <- d %>% group_by(blk) %>% mutate(r = rank(-MI_obs)) %>% group_by(.data[[f]]) %>% summarise(mean_rank = mean(r), .groups = "drop")
  cat(sprintf("%-7s chi2 = %.2f, df = %d, p = %.2e | mean rank (1 = best): %s\n", f, ft$statistic, ft$parameter, ft$p.value,
              paste(sprintf("%s %.2f", mr[[1]], mr$mean_rank), collapse = ", "))) }
for (f in c("bb", "decomp", "pool")) fr(f)

cat("\n=== 5. Backbone: paired Wilcoxon signed-rank over the 48 matched (task, decomp, pool) blocks, Holm ===\n")
w <- g %>% select(block, bb, MI_obs) %>% pivot_wider(names_from = bb, values_from = MI_obs)
pw <- combn(BB, 2, function(p) { t <- wilcox.test(w[[p[1]]], w[[p[2]]], paired = TRUE, exact = FALSE)
  data.frame(a = p[1], b = p[2], a_higher = sum(w[[p[1]]] > w[[p[2]]]), n = nrow(w),
             median_diff = median(w[[p[1]]] - w[[p[2]]]), p = t$p.value) }, simplify = FALSE) %>% bind_rows()
pw$p_holm <- p.adjust(pw$p, "holm"); print(pw, digits = 3, row.names = FALSE)

cat("\n=== 6. Backbone ranking WITHIN each task (16 matched blocks per task): mean rank and Friedman p ===\n")
for (tk in unique(g$task)) {
  d <- g %>% filter(task == tk) %>% mutate(blk = paste(decomp, pool))
  ft <- friedman.test(d$MI_obs, factor(d$bb), factor(d$blk))
  mr <- d %>% group_by(blk) %>% mutate(r = rank(-MI_obs)) %>% group_by(bb) %>% summarise(mr = mean(r), .groups = "drop")
  mr <- mr[match(BB, mr$bb), ]
  cat(sprintf("%-9s Friedman p = %.1e | %s\n", tk, ft$p.value, paste(sprintf("%s %.2f", mr$bb, mr$mr), collapse = "  ")))
}
cat("Kendall's W, agreement of the three tasks on the backbone order (1 = identical rankings):\n")
R <- sapply(unique(g$task), function(tk) { d <- g %>% filter(task == tk) %>% mutate(blk = paste(decomp, pool)) %>%
  group_by(blk) %>% mutate(r = rank(-MI_obs)) %>% group_by(bb) %>% summarise(mr = mean(r), .groups = "drop")
  rank(d$mr[match(BB, d$bb)]) })
W <- function(R) { m <- ncol(R); n <- nrow(R); S <- sum((rowSums(R) - m * (n + 1) / 2)^2); 12 * S / (m^2 * (n^3 - n)) }
cat(sprintf("  backbones: W = %.2f\n", W(R)))
Rd <- sapply(unique(g$task), function(tk) { d <- g %>% filter(task == tk) %>% mutate(blk = paste(bb, pool)) %>%
  group_by(blk) %>% mutate(r = rank(-MI_obs)) %>% group_by(decomp) %>% summarise(mr = mean(r), .groups = "drop"); rank(d$mr) })
cat(sprintf("  decompositions: W = %.2f\n", W(Rd)))

cat("\n=== 7. Decomposition ranking within each task (blocks = bb x pool, 16 per task) ===\n")
for (tk in unique(g$task)) {
  d <- g %>% filter(task == tk) %>% mutate(blk = paste(bb, pool))
  ft <- friedman.test(d$MI_obs, factor(d$decomp), factor(d$blk))
  mr <- d %>% group_by(blk) %>% mutate(r = rank(-MI_obs)) %>% group_by(decomp) %>% summarise(mr = mean(r), .groups = "drop")
  cat(sprintf("%-9s Friedman p = %.1e | %s\n", tk, ft$p.value, paste(sprintf("%s %.2f", mr$decomp, mr$mr), collapse = "  ")))
}

cat("\n=== 8. Configuration level: 64 (bb, decomp, pool) configs, ranked within task (1 = best of 64) ===\n")
cfg <- g %>% group_by(task) %>% mutate(rank_in_task = rank(-MI_obs), gap_nats = max(MI_obs) - MI_obs,
                                       gap_nullsd = gap_nats / mean(null_sd)) %>% ungroup() %>%
  group_by(bb, decomp, pool) %>% summarise(mean_rank = mean(rank_in_task), worst_rank = max(rank_in_task),
    max_gap_nats = max(gap_nats), max_gap_nullsd = max(gap_nullsd), min_MI = min(MI_obs), all_sig = all(sig), .groups = "drop")
cat("Top 10 by mean rank across the three tasks:\n")
print(as.data.frame(cfg %>% arrange(mean_rank) %>% head(10)), digits = 3)
cat("Top 10 by minimax regret (smallest worst-case gap to the task's best, in nats):\n")
print(as.data.frame(cfg %>% arrange(max_gap_nats) %>% head(10)), digits = 3)
cat(sprintf("Configs within 0.008 nats (reproducibility bound) of the best in ALL three tasks: %d\n", sum(cfg$max_gap_nats <= 0.008)))
cat(sprintf("Configs within 2 null SD of the best in ALL three tasks: %d\n", sum(cfg$max_gap_nullsd <= 2)))
cat("Best configuration per task:\n")
print(as.data.frame(g %>% group_by(task) %>% slice_max(MI_obs, n = 1) %>% select(task, bb, decomp, pool, MI_obs, z)), digits = 3)

cat("\n=== 9. Selected-kernel table (48 rows): what got selected, and max-T significance, by backbone ===\n")
print(table(best$bb, best$kern)[BB, ])
print(best %>% group_by(bb) %>% summarise(arms = n(), maxT_sig = sum(p_maxT < 0.05), mean_MI = mean(MI_obs),
                                          mean_kappa = mean(eig_cond), .groups = "drop") %>% arrange(match(bb, BB)), digits = 3)
cat("Selected-arm MI by task and backbone (mean over the 4 poolings):\n")
print(round(with(best, tapply(MI_obs, list(task, bb), mean))[, BB], 4))

cat("\n=== 10. Incremental table: covariate share and additivity failures ===\n")
print(inc %>% mutate(share = ifelse(combined > 0, covar_only / combined, NA), drop = emb_only - combined) %>%
  group_by(task, bb) %>% summarise(covar = mean(covar_only), cov_share = mean(share), max_drop = max(drop), .groups = "drop") %>%
  arrange(task, match(bb, BB)) %>% as.data.frame(), digits = 3)

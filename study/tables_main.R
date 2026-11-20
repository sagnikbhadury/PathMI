################################################################
# tables_main.R  (20 Nov 2026)
#
# Generates, from the current grid, the tables that the manuscript used
# to carry as typed LaTeX: Table 1 (tab:four) and Supplementary Tables
# tab:pairwise, tab:var, tab:cond, tab:inc, tab:prad_cov. Definitions
# follow results_four_backbones_numbers.R and config_comparison_tests.R.
# Output: manuscript/table_{four,pairwise,var,cond,inc,prad_cov}.tex
################################################################
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); MS <- file.path(dirname(TRI), "manuscript")
BB <- c("UNI v2", "Virchow2", "H-optimus-1", "CONCH v1.5")
e2b <- c("UNI V2" = "UNI v2", "VIRCHOW" = "Virchow2", "H-OPTIMUS 1" = "H-optimus-1", "CONCH V1.5" = "CONCH v1.5")
g <- as.data.frame(read_excel(file.path(TRI, "mi summary table.xlsx"))) %>% mutate(bb = factor(e2b[Embeddings], BB), task = paste0(Disease, "/", target), decomp = tolower(decomp))
best <- as.data.frame(read_excel(file.path(TRI, "best kpca summary.xlsx"))) %>% mutate(bb = factor(e2b[Embeddings], BB), task = paste0(Disease, "/", target))
inc <- as.data.frame(read_excel(file.path(TRI, "incremental mi summary.xlsx"))) %>% mutate(bb = factor(e2b[Embeddings], BB), task = paste0(Disease, "/", target))
u <- read.csv(file.path(TRI, "usable", "usable_recoverability.csv"), stringsAsFactors = FALSE)
HY <- u %>% group_by(task = paste0(cohort, "/", target), bb = model) %>% summarise(n = first(n), H = first(H_Y), .groups = "drop")
TASKS <- c("CRC/MSI", "PRAD/ERG", "PRAD/ETS")
sg <- function(x, d = 3) ifelse(x < 0, sprintf("$-$%.*f", d, -x), sprintf("$+$%.*f", d, x))
pf <- function(p) ifelse(p < 1e-4, "$<0.0001$", ifelse(p < 0.001, sprintf("%.4f", p), ifelse(p < 0.01, sprintf("%.3f", p), sprintf("%.2f", p))))

# ---- Table 1 --------------------------------------------------------------------
w <- g %>% select(task, decomp, agg_method, bb, MI_obs) %>% group_by(task, decomp, agg_method) %>% mutate(rk = rank(-MI_obs)) %>% ungroup()
t1 <- g %>% group_by(task, bb) %>% summarise(sig = sum(p_value < 0.05), mx = max(MI_obs), .groups = "drop") %>%
  left_join(best %>% group_by(task, bb) %>% summarise(sel = mean(MI_obs), .groups = "drop"), by = c("task", "bb")) %>%
  left_join(w %>% group_by(task, bb) %>% summarise(mr = mean(rk), .groups = "drop"), by = c("task", "bb")) %>%
  left_join(HY %>% mutate(bb = factor(bb, BB)), by = c("task", "bb")) %>% arrange(factor(task, TASKS), bb)
L <- c(); for (tk in TASKS) { d <- t1[t1$task == tk, ]; for (i in seq_len(nrow(d))) L <- c(L, sprintf("%s & %s & %d & %d & %.3f & %.0f & %.3f & %.2f \\\\",
  if (i == 1) tk else "", d$bb[i], d$n[i], d$sig[i], d$mx[i], 100 * d$mx[i] / d$H[i], d$sel[i], d$mr[i])); if (tk != "PRAD/ETS") L <- c(L, "\\midrule") }
writeLines(c("\\begin{tabular}{llcccccc}", "\\toprule", "Task & Model & $n$ & Above null & Maximum & \\% of $H(Y)$ & Selected mean & Mean rank \\\\", "\\midrule", L, "\\bottomrule", "\\end{tabular}"),
           file.path(MS, "table_four.tex"))

# ---- pairwise -----------------------------------------------------------------------
wb <- g %>% select(task, decomp, agg_method, bb, MI_obs) %>% pivot_wider(names_from = bb, values_from = MI_obs)
pairs <- list(c("H-optimus-1", "UNI v2"), c("H-optimus-1", "Virchow2"), c("UNI v2", "Virchow2"), c("H-optimus-1", "CONCH v1.5"), c("UNI v2", "CONCH v1.5"), c("Virchow2", "CONCH v1.5"))
pw <- do.call(rbind, lapply(pairs, function(p) { d <- wb[[p[1]]] - wb[[p[2]]]
  data.frame(a = p[1], b = p[2], hi = sum(d > 0), med = median(d), p = suppressWarnings(wilcox.test(wb[[p[1]]], wb[[p[2]]], paired = TRUE, exact = FALSE)$p.value)) }))
pw$ph <- p.adjust(pw$p, "holm")
writeLines(c("\\begin{tabular}{llccc}", "\\toprule", "First & Second & First higher (of 48) & Median difference & Holm $p$ \\\\", "\\midrule",
  sprintf("%s & %s & %d & %s & %s \\\\", pw$a, pw$b, pw$hi, sg(pw$med), pf(pw$ph)), "\\bottomrule", "\\end{tabular}"), file.path(MS, "table_pairwise.tex"))

# ---- variance partition ----------------------------------------------------------
vp <- do.call(rbind, lapply(TASKS, function(tk) { d <- g[g$task == tk, ]; a <- anova(lm(MI_obs ~ bb + decomp + agg_method, data = d)); ss <- a[["Sum Sq"]]; t <- sum(ss)
  data.frame(task = tk, m = 100 * ss[1] / t, dc = 100 * ss[2] / t, pl = 100 * ss[3] / t, r = 100 * ss[4] / t) }))
writeLines(c("\\begin{tabular}{lcccc}", "\\toprule", "Task & Foundation model & Decomposition & Pooling & Residual \\\\", "\\midrule",
  sprintf("%s & %.1f & %.1f & %.1f & %.1f \\\\", vp$task, vp$m, vp$dc, vp$pl, vp$r), "\\bottomrule", "\\end{tabular}"), file.path(MS, "table_var.tex"))

# ---- conditioning --------------------------------------------------------------------
kc <- best %>% group_by(task, bb) %>% summarise(k = mean(eig_cond), .groups = "drop") %>% pivot_wider(names_from = bb, values_from = k)
rr <- best %>% group_by(task) %>% summarise(r4 = cor(MI_obs, eig_cond), r3 = cor(MI_obs[bb != "CONCH v1.5"], eig_cond[bb != "CONCH v1.5"]), .groups = "drop")
kc <- left_join(kc, rr, by = "task") %>% arrange(factor(task, TASKS))
writeLines(c("\\begin{tabular}{lcccccc}", "\\toprule", " & \\multicolumn{4}{c}{Mean $\\kappa$ of selected kernels} & \\multicolumn{2}{c}{$r(\\mathrm{MI},\\kappa)$} \\\\",
  "\\cmidrule(lr){2-5}\\cmidrule(lr){6-7}", "Task & UNI v2 & Virchow2 & H-optimus-1 & CONCH v1.5 & Four models & Without CONCH \\\\", "\\midrule",
  sprintf("%s & %.0f & %.0f & %.0f & %.0f & %s & %s \\\\", kc$task, kc[["UNI v2"]], kc[["Virchow2"]], kc[["H-optimus-1"]], kc[["CONCH v1.5"]], sg(kc$r4, 2), sg(kc$r3, 2)),
  "\\bottomrule", "\\end{tabular}"), file.path(MS, "table_cond.tex"))

# ---- incremental (CRC) -----------------------------------------------------------------
ic <- inc %>% filter(Disease == "CRC") %>% mutate(share = 100 * covar_only / combined, syn = combined - covar_only - emb_only) %>% arrange(bb, factor(agg_method, c("mean", "norm", "moment", "attention")))
writeLines(c("\\begin{tabular}{llcccccc}", "\\toprule", "Model & Pooling & Covariate & Embedding & Combined & Incremental & Share (\\%) & Synergy \\\\", "\\midrule",
  sprintf("%s & %s & %.4f & %.4f & %.4f & %.4f & %.0f & %s \\\\", ic$bb, ic$agg_method, ic$covar_only, ic$emb_only, ic$combined, ic$incremental_subho_mi, ic$share, sg(ic$syn, 4)),
  "\\bottomrule", "\\end{tabular}"), file.path(MS, "table_inc.tex"))

# ---- PRAD extended covariates (means over the four selected arms of each model and state) ---------
pe <- read.csv(file.path(TRI, "prad_incremental_batch_purity", "prad_incremental_batch_purity.csv"), stringsAsFactors = FALSE)
ps <- pe %>% group_by(target, model) %>% summarise(b = mean(covar_C0_refit), p = mean(covar_purity), bt = mean(covar_batch), a = mean(covar_all), e = mean(emb_only),
  cmb = mean(combined_all), ia = mean(incremental_all), sh = 100 * mean(covar_all) / mean(combined_all), .groups = "drop") %>% arrange(target, factor(model, BB))
f4 <- function(x) ifelse(abs(x) < 5e-5, "0", sprintf("%.4f", x)); L <- c()
for (tg in c("ERG", "ETS")) { d <- ps[ps$target == tg, ]; for (i in seq_len(nrow(d))) L <- c(L, sprintf("%s & %s & %s & %s & %s & %s & %.4f & %.4f & %.4f & %.0f \\\\",
  if (i == 1) tg else "", d$model[i], f4(d$b[i]), f4(d$p[i]), f4(d$bt[i]), f4(d$a[i]), d$e[i], d$cmb[i], d$ia[i], d$sh[i])); if (tg == "ERG") L <- c(L, "\\midrule") }
writeLines(c("\\begin{tabular}{llcccccccc}", "\\toprule", "State & Model & Base & + purity & + batch & + both & Embedding & Combined & Incremental & Share (\\%) \\\\", "\\midrule", L,
  "\\bottomrule", "\\end{tabular}"), file.path(MS, "table_prad_cov.tex"))
cat("tables written\n"); print(as.data.frame(t1), digits = 3); print(pw, digits = 3); print(vp, digits = 3); print(as.data.frame(kc), digits = 3)

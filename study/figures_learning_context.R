################################################################
# figures_learning_context.R
#
# Figure: evidence sufficiency, an unresolved record and a hidden context.
#   a  margin of the 90% lower confidence bound over the resolution
#      threshold, LCB90(n) - {Inull(n) + eps_rep}, against the number of
#      patients, one line per model, one facet per state and one for the
#      negative control; the first n above zero is N_evidence;
#   b  usable recoverability (in-fold ridge readout) against n;
#   c  three records side by side: a resolved arm, the unresolved Record
#      arm(s) and the negative control, each with Ihat and its 90%
#      interval, Inull + eps_rep and I_usable;
#   d  hidden-context simulation: PathMI on the pooled cohort, on
#      context A and on context B, against the truth, for three
#      context fractions.
# Output: figures/fig6_evidence_context.{pdf,png} (copied to manuscript/figures)
################################################################
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork); library(tidyr) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
TRI <- normalizePath(dirname(this_file), winslash = "/"); FIG <- file.path(TRI, "figures"); MS_FIG <- file.path(dirname(TRI), "manuscript", "figures")
BB <- c("UNI v2", "Virchow2", "H-optimus-1", "CONCH v1.5"); code <- c("UNI v2" = "uni_v2", "Virchow2" = "virchow2", "H-optimus-1" = "hoptimus1", "CONCH v1.5" = "conch_v15")
pal_m <- c("UNI v2" = "#0072B2", "Virchow2" = "#D55E00", "H-optimus-1" = "#009E73", "CONCH v1.5" = "#CC79A7", "Control" = "grey35")
th <- theme_bw(base_size = 7.5) + theme(panel.grid.minor = element_blank(), legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6.5),
                                      plot.tag = element_text(face = "bold", size = 9), strip.background = element_rect(fill = "grey95", colour = NA), strip.text = element_text(size = 6.8))
CTL <- NULL   # set from the learning summary below

lc <- read.csv(file.path(TRI, "learning", "learning_summary.csv"), stringsAsFactors = FALSE)
CTL <- unique(grep("^CONTROL", lc$id, value = TRUE))
ne <- read.csv(file.path(TRI, "learning", "n_evidence.csv"), stringsAsFactors = FALSE)
rec <- read.csv(file.path(TRI, "record", "recoverability_record.csv"), stringsAsFactors = FALSE)
lab <- function(id) { p <- strsplit(id, "|", fixed = TRUE)[[1]]
  if (p[1] == "CONTROL") return(c(task = "Negative control", model = "Control"))
  c(task = sprintf("%s/%s", p[1], p[3]), model = names(code)[code == p[2]]) }
L <- t(sapply(lc$id, lab)); lc$task <- L[, "task"]; lc$model <- factor(L[, "model"], c(BB, "Control"))
lc$task <- factor(lc$task, c("CRC/MSI", "PRAD/ERG", "PRAD/ETS", "Negative control"))
lc <- lc %>% filter(!(source == "subsample" & n == floor(N / 2) & !(n %in% c(50, 100, 150, 200))))   # the N/2 level serves the s.e. check only
lc$margin <- lc$LCB - lc$threshold

pa <- ggplot(lc, aes(n, margin, colour = model)) + geom_hline(yintercept = 0, linewidth = .3) +
  geom_line(linewidth = .45) + geom_point(aes(shape = resolved), size = 1.4) + scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  facet_wrap(~task, nrow = 1, scales = "free_x") + scale_colour_manual(values = pal_m, name = NULL) + scale_x_continuous(breaks = c(50, 100, 200, 300, 430)) +
  labs(x = "Patients n", y = expression(LCB[90](n) - (I[null](n) + epsilon[rep])), tag = "a") + th + theme(legend.position = "right")

pb <- ggplot(lc, aes(n, U_mean, colour = model)) + geom_hline(yintercept = 0, linewidth = .3, colour = "grey60") +
  geom_line(linewidth = .45) + geom_point(size = 1.1) +
  facet_wrap(~task, nrow = 1, scales = "free_x") + scale_colour_manual(values = pal_m, name = NULL) + scale_x_continuous(breaks = c(50, 100, 200, 300, 430)) +
  labs(x = "Patients n", y = "Usable recoverability, nats", tag = "b") + th + theme(legend.position = "none")

# c: records at the full cohort
full <- lc %>% filter(source == "full cohort") %>% mutate(lo = pmax(LCB, 0), hi = I_mean + 1.645 * se)
unres <- full$id[!full$resolved & full$id != CTL]
best_res <- full %>% filter(resolved, id != CTL) %>% slice_max(I_mean, n = 1) %>% pull(id)
show <- full %>% filter(id %in% c(best_res, unres, CTL)) %>%
  mutate(rec_lab = ifelse(id == CTL, "Control (permuted)", sprintf("%s, %s", task, model)), status = ifelse(resolved, "resolved", "unresolved"))
show$rec_lab <- factor(show$rec_lab, show$rec_lab[order(show$I_mean)])
pc <- ggplot(show, aes(rec_lab)) +
  geom_errorbar(aes(ymin = 0, ymax = threshold), width = .55, colour = "grey55", linewidth = .35) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .7, colour = "black") + geom_point(aes(y = I_mean, fill = status), shape = 21, size = 2.2) +
  geom_point(aes(y = U_mean), shape = 4, size = 1.8, colour = "#D55E00") +
  geom_text(aes(y = hi, label = status), hjust = -0.15, size = 2.2) +
  scale_fill_manual(values = c(resolved = "black", unresolved = "white"), guide = "none") +
  scale_y_continuous(limits = c(min(0, min(show$U_mean)) - 0.005, max(show$hi) * 1.45)) +
  labs(x = NULL, y = "Nats", tag = "c") + th + coord_flip()

# d: hidden-context simulation
cs <- read.csv(file.path(TRI, "simulation_context", "context_draws.csv"))
ct <- read.csv(file.path(TRI, "simulation_context", "context_truth.csv"))
cn <- read.csv(file.path(TRI, "simulation_context", "context_null.csv")) %>% group_by(n) %>% summarise(null = mean(I_null), .groups = "drop")
cd <- cs %>% pivot_longer(c(I_pooled, I_A, I_B), names_to = "set", values_to = "I") %>%
  mutate(set = factor(c(I_pooled = "Pooled (S hidden)", I_A = "Context A", I_B = "Context B")[set], c("Pooled (S hidden)", "Context A", "Context B")),
         pi_lab = sprintf("P(S = A) = %.2f", pi_A))
tr <- ct %>% transmute(pi_A, `Pooled (S hidden)` = I_pooled, `Context A` = I_A, `Context B` = I_B) %>% pivot_longer(-pi_A, names_to = "set", values_to = "truth") %>%
  mutate(set = factor(set, levels(cd$set)), pi_lab = sprintf("P(S = A) = %.2f", pi_A))
nl <- cs %>% distinct(pi_A, n_A, n_B) %>% mutate(`Pooled (S hidden)` = 400, `Context A` = n_A, `Context B` = n_B) %>% select(-n_A, -n_B) %>%
  pivot_longer(-pi_A, names_to = "set", values_to = "n") %>% left_join(cn, by = "n") %>% mutate(set = factor(set, levels(cd$set)), pi_lab = sprintf("P(S = A) = %.2f", pi_A))
pd <- ggplot(cd, aes(set, I)) + geom_boxplot(width = .55, outlier.size = .4, linewidth = .3, fill = "grey92") +
  geom_point(data = tr, aes(y = truth), shape = 23, size = 2, fill = "#E69F00") +
  geom_errorbar(data = nl, aes(y = null, ymin = null, ymax = null), width = .7, linetype = "22", colour = "#0072B2", linewidth = .35) +
  facet_wrap(~pi_lab, nrow = 1) + labs(x = NULL, y = "Recoverability, nats", tag = "d") +
  th + theme(axis.text.x = element_text(angle = 30, hjust = 1))

g <- pa / pb / (pc + pd + plot_layout(widths = c(1.3, 2))) + plot_layout(heights = c(1, 1, 1.25))
ggsave(file.path(FIG, "fig6_evidence_context.pdf"), g, width = 180, height = 150, units = "mm", device = cairo_pdf)
ggsave(file.path(FIG, "fig6_evidence_context.png"), g, width = 180, height = 150, units = "mm", dpi = 200)
file.copy(file.path(FIG, "fig6_evidence_context.pdf"), file.path(MS_FIG, "fig6_evidence_context.pdf"), overwrite = TRUE)
cat("written fig6_evidence_context\n")

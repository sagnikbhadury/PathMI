################################################################
# USAGE
#--------------------------------------------------------------#


rm(list = ls())
set.seed(2026)


source("PathMI_SUBHO_logistic.R")
source("PathMI_LASSO.R")
source("PathMI_OLS.R")
source("PathMI_DOE_linear.R")
source("simulation_MI_comparison.R")

# # Validate first:
# pilot <- run_pilot()
# head(pilot)

# Full run:
results <- run_simulation(n_reps = 50)
summary <- summarize_results(results)

# All figures saved as PDFs:
plot_all(results, save_pdf = TRUE, prefix = "sim_")

# Or individually in interactive session:
plot_mi_scatter(results)
plot_tpr_fpr(results)

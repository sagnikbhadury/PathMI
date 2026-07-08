# PathMI

An R package estimating mutual information between patient-level
histopathology tile embeddings and binary molecular targets (e.g. MSI
status, gene fusion status), via three interchangeable classifiers:

* **OLS/DOE** (`run_PathMI_OLS()`) -- unpenalized, for the low-dimensional
  regime (features << patients).
* **LASSO-CV logistic** (`run_PathMI_LASSO()`) -- for the high-dimensional
  regime.
* **SUBHO** (`pathMI_pipeline()`) -- a fully Bayesian Polya-Gamma logistic
  horseshoe estimator with a marginal-likelihood feature walk.

Each estimator ships with a permutation test and Benjamini-Hochberg FDR
correction across molecular targets.

## Installation

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/PathMI")
```

## Usage

```r
library(PathMI)

embedding_list <- readRDS("patient_embeddings.rds")  # patient -> tile x feature matrix
targets <- list(MSI = setNames(msi_df$status, msi_df$patient_id))

# High-dimensional embeddings (e.g. 2560-d): use LASSO
out <- run_PathMI_LASSO(embedding_list, targets, nfolds = 10, B = 500)
print(out$summary)

# Fully Bayesian alternative with feature selection
out_subho <- pathMI_pipeline(embedding_list, targets, n.perm = 500, v_beta = 1)
plot_MI_null(out_subho, "MSI")
plot_likelihood_surface(out_subho, "MSI")
```

## Package layout

* `R/utils.R` -- `aggregate_embeddings()`, shared across all three estimators.
* `R/ols.R` -- unpenalized OLS/DOE estimator.
* `R/lasso.R` -- LASSO-CV logistic estimator.
* `R/subho.R` -- Polya-Gamma logistic horseshoe (SUBHO) estimator and its
  diagnostics (`plot_MI_null()`, `plot_likelihood_surface()`,
  `v_beta_sensitivity()`).

See `?pathMI_pipeline`, `?run_PathMI_LASSO`, and `?run_PathMI_OLS` for full
documentation.

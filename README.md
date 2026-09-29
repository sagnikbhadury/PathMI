# PathMI

PathMI measures **molecular recoverability**: the mutual information between a patient's
pathology foundation-model embedding and a binary molecular state (microsatellite
instability, gene-fusion status), estimated for a few hundred patients and hundreds to
thousands of embedding coordinates. The estimator is Bayesian and works when the dimension
exceeds the sample size: a horseshoe prior under Polya-Gamma augmentation ranks the
coordinates, a walk along the nested models that the ranking defines scores each model by
its Laplace marginal likelihood, and the estimate is the log Bayes factor of the selected
model against the label's marginal, per observation. Around the estimate the package
provides a permutation null in which the whole estimator is refitted, and the repository
holds the results of applying it to four foundation models on TCGA colorectal and prostate
cancer, including tile-level maps of where on the slide the evidence lies.

This repository accompanies the manuscript *Measuring molecular recoverability in pathology
foundation models* (Bhadury, Kulkarni and Rao, 2026, in preparation).

## Contents

| Path | What it holds |
|---|---|
| `R/`, `man/`, `DESCRIPTION`, `NAMESPACE` | The R package: the PathMI estimator (`pathMI_pipeline()`, `subho_MI_logistic()`), two comparison estimators (`run_PathMI_LASSO()`, `run_PathMI_OLS()`), permutation tests and diagnostics. |
| `results/` | The results deposited with the manuscript: the 192-cell recoverability grid, selected projections, incremental recoverability with and without batch and purity, tissue-sampling burden, simulation results, per-run outputs for every cohort and model, and slide-level tile summaries. See [`results/README.md`](results/README.md). |
| Release [`tile-scores-v1`](https://github.com/sagnikbhadury/PathMI/releases/tag/tile-scores-v1) | The per-tile scores (4.78 million colorectal and 4.11 million prostate tiles, 173 MB and 300 MB), which exceed GitHub's file-size limit. |

## Installation

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/PathMI")
```

Imports `glmnet`, `MASS` and `BayesLogit`. Developed and run with R 4.4.3.

## Usage

The input is a list with one element per patient, each a tiles-by-features matrix of
embeddings, plus a named vector of binary labels per target.

```r
library(PathMI)

embedding_list <- readRDS("patient_embeddings.rds")  # patient -> tile x feature matrix
targets <- list(MSI = setNames(msi_df$status, msi_df$patient_id))

# PathMI: horseshoe ranking, marginal-likelihood walk, Bayes-factor estimate,
# with a permutation null that refits the whole estimator
out <- pathMI_pipeline(embedding_list, targets, n.perm = 200, v_beta = 1)
print(out$summary)
plot_MI_null(out, "MSI")                 # estimate against its permutation null
plot_likelihood_surface(out, "MSI")      # the marginal-likelihood walk and k*

# Comparison estimators
out_lasso <- run_PathMI_LASSO(embedding_list, targets, nfolds = 10, B = 200)
out_ols   <- run_PathMI_OLS(embedding_list, targets, B = 200)   # p < n only
```

`aggregate_embeddings()` pools tiles to one vector per patient (mean by default); the
manuscript's study also uses norm-weighted, moment and untrained-attention pooling and reduces
the pooled matrix to 50 coordinates by linear or kernel PCA before estimation, so that the
estimator sees an n-by-50 matrix.

### The estimator in brief

For a representation **Z** (n patients by p coordinates) and label y:

1. **Ranking.** A logistic regression with the horseshoe prior is sampled by Gibbs under
   Polya-Gamma augmentation; coordinates are ordered by the absolute posterior median of their
   coefficients.
2. **Nested model walk.** The nested models M_0 ⊂ M_1 ⊂ ... ⊂ M_K, where M_k holds the
   intercept and the k highest-ranked coordinates, are scored by the Laplace marginal
   likelihood under a working Gaussian prior of variance `v_beta` plus a uniform Beta-Bernoulli
   prior on model size; k* is the argmax.
3. **Estimate.** Î(Z; Y) = [log p(y | M_k*) − log p(y | M_0)] / n, floored at zero.
4. **Null.** Steps 1 to 3 are repeated on permuted labels; the permutation p value, null mean
   and null s.d. are reported beside the estimate.

The nested-model walk carries over a step-up selection principle introduced for Gaussian
graphical models (Bhadury, Mitra and Gaskins, 2023, arXiv:2308.02713), which is why the
function names in `R/subho.R` carry the `subho` prefix.

## Package layout

* `R/utils.R` -- `aggregate_embeddings()`, `preprocess_embeddings()`.
* `R/subho.R` -- the PathMI estimator: the Polya-Gamma horseshoe sampler, the
  marginal-likelihood walk (`selection.logistic()`), `subho_MI_logistic()`,
  `permutation_test_logistic()`, `pathMI_pipeline()`, and the diagnostics `plot_MI_null()`,
  `plot_likelihood_surface()`, `v_beta_sensitivity()`.
* `R/lasso.R` -- LASSO-CV logistic comparison estimator.
* `R/ols.R` -- unpenalised OLS/DOE comparison estimator (p < n).

See `?pathMI_pipeline`, `?run_PathMI_LASSO` and `?run_PathMI_OLS`.

## Data

Whole slide images and molecular annotation of the TCGA COAD, READ and PRAD cohorts are
available from the Genomic Data Commons (https://portal.gdc.cancer.gov). Tile embeddings were
extracted with Trident from the public UNI v2, Virchow2, H-optimus-1 and CONCH v1.5 weights.

## License

MIT; see `LICENSE`.

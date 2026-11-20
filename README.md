# PathMI

PathMI measures **molecular recoverability**: how much information a pathology foundation model's embedding of a patient's slides carries about a binary molecular state, such as microsatellite instability or an ERG fusion. The estimate is a mutual information in nats. It is built for a few hundred patients and many more embedding coordinates than patients.

Around the estimator, the package provides what is needed to report an estimate as a measurement:

- a permutation null and a resolution scale;
- patient-level uncertainty;
- what the representation adds over clinical covariates;
- what a held-out predictor can use;
- how many patients and how much tissue a claim needs;
- a Recoverability Record that records all of this with its provenance.

This repository accompanies the manuscript *Measuring molecular recoverability in pathology foundation models* (Bhadury, Kulkarni and Rao, 2026, in preparation).

**Every example below runs on simulated data.** The embeddings used in the study are not distributed with this repository.

## Installation

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/PathMI")
```

Imports `glmnet`, `MASS` and `BayesLogit`; `jsonlite` is needed to write Records as JSON. Developed with R 4.4.3.

## The estimator in brief

Take a representation **Z** (n patients by p coordinates) and a binary label y. PathMI then works in four steps:

1. **Ranking.** A logistic regression with the horseshoe prior is sampled by Gibbs under Pólya–Gamma augmentation. The coordinates are ordered by the absolute posterior median of their coefficients.
2. **Walk.** The nested models M₀ ⊂ M₁ ⊂ … ⊂ M_K are scored by their Laplace marginal likelihood with a beta-binomial prior on size. M_k holds the intercept and the k highest-ranked coordinates, with K = min(p, n − 1). The selected size k\* maximises the score.
3. **Estimate.** Î(Z;Y) = [log p(y | M_k\*) − log p(y | M₀)] / n, floored at zero. This is the log Bayes factor per patient.
4. **Null.** Steps 1–3 are refitted on permuted labels.

The nested-model walk carries over a step-up selection principle introduced for Gaussian graphical models (Bhadury, Mitra and Gaskins, 2023, arXiv:2308.02713). This is why the engine functions carry the `subho` prefix.

## A worked example

[`tools/readme_figures.R`](tools/readme_figures.R) reproduces every number and figure in this section. Running `Rscript tools/readme_figures.R 20` takes about eight minutes on 20 cores.

### 1. Simulate slides, pool them per patient, project

`simulate_cohort()` draws slides as grids of tile embeddings. Each slide has a contiguous tumour region whose share varies, as purity does. In patients with the state, the tumour tiles are shifted along a few directions, so the state leaves a trace in part of the tissue only. Each patient also carries a covariate correlated with the state, standing in for stage.

We simulate two cohorts of 400 patients. One has a strong trace and one a weak trace. Each cohort has about 460 slides, because some patients have several.

```r
library(PathMI)
sim <- simulate_cohort(400, prevalence = 0.3, D = 64, tiles_mean = 150,
                       effect = 0.45, tile_sd = 0.4, p_multi = 0.08, seed = 11)

X <- pool_cohort(sim$tiles, sim$patient, method = "mean")  # one row per patient
y <- sim$y[rownames(X)]
Zs <- lapply(c(pca = "pca", cosine = "cosine", rbf = "rbf", poly = "poly"),
             function(dc) project_embeddings(X, dc, d = 20))
```

The unit is the patient. `pool_cohort()` pools each slide and then averages a patient's slides. Treating sibling slides as separate observations would count the same patient's evidence twice in the marginal likelihood. It would also leak information between folds when the data are cross-fitted.

`pool_tiles()` also offers norm-weighted, moment (mean and s.d.) and untrained-attention pooling. `fit_projection()` and `apply_projection()` separate fitting from projecting, so that held-out patients can be projected with a training fit. They use the out-of-sample extension of kernel PCA.

### 2. Estimate, null, and selection over representations

```r
fits  <- lapply(Zs, pathmi_fit, y = y)
perms <- draw_permutations(length(y), B = 100)           # shared by all candidates
nulls <- Map(function(Z, f) pathmi_null(Z, y, perms, observed = f, cores = 20), Zs, fits)
sel   <- select_representation(nulls)                    # max-T over a joint permutation
fits$poly
#> PathMI estimate 0.0516 nats (8.2% of H(Y) = 0.631) | k* = 6 | n = 400
```

<p align="center"><img src="man/figures/fig1_estimate.png" width="100%"></p>

**a.** In the strong cohort, the selected cubic-kernel projection gives 0.052 nats. The null from 100 permutations has mean 0.0017, and no permuted estimate reaches the observed one (p = 0.010, the floor for 100 permutations).

**b.** The walk peaks at k\* = 6 coordinates and falls on either side.

**c.** All four projections recover 0.039–0.052 nats from the strong cohort, against a null band near 0.01. In the weak cohort no projection clears the band. The best, the cosine kernel, gives 0.008 nats with a max-T p of 0.059.

Because the four candidates share one set of permutations, the max-T p value accounts for having picked the best of the four.

### 3. Resolution and uncertainty

Some quantities must be identical by design: a refit with a new seed, a refit after the coordinates' signs are flipped, and a refit after appending columns that carry no information. Any difference among them is error of the measurement itself. The largest difference sets the resolution ε_rep, and contrasts are reported in units of it. Patient-level uncertainty comes from refitting on random halves of the patients, with the delete-d jackknife standard error (Shao and Wu, 1989).

```r
eps <- reproducibility_bound(
  replicate  = sapply(names(Zs), function(k) replicate_refit(Zs[[k]], y, fits[[k]])),
  sign_flip  = sapply(names(Zs), function(k) sign_flip_refit(Zs[[k]], y, fits[[k]])),
  null_block = sapply(names(Zs), function(k) null_block_refit(Zs[[k]], y, fits[[k]])))
hs <- half_sampling(Zs$poly, y, B = 100, cores = 20)
```

Across both cohorts, the replicate and sign-flip refits reproduced exactly. The null block moved the estimate by at most 0.0040 nats, so ε_rep = 0.0040. For the strong cohort, the half-sampling standard error is 0.024 and the 90 % interval runs from 0.013 to 0.090.

We do not use the bootstrap with replacement. A patient drawn twice counts twice in the marginal likelihood, and the estimate is biased upward.

### 4. Represented, attributable and usable information

These are three different quantities:

- **Represented:** the information the embedding carries, I(Z;Y).
- **Attributable:** the part not already carried by covariates, I(Z;Y | C).
- **Usable:** what a predictor fitted on other patients achieves, measured as held-out log-likelihood gain (Xu et al., 2020).

```r
n_tiles <- tapply(sapply(sim$tiles, nrow), sim$patient, sum)[rownames(X)]
C <- cbind(stage = scale(sim$covariate[rownames(X)]), log_tiles = scale(log(n_tiles)))
inc <- incremental_recoverability(Zs$poly, C, y)                 # attributable
use <- usable_recoverability(X, y, decomp = "poly", d = 20, K = 5, R = 5)  # usable
```

`usable_recoverability()` fits the standardisation, the projection and the readout inside the training patients of each fold. Supplying `group =` keeps all rows of a patient in one fold.

<p align="center"><img src="man/figures/fig2_evidence.png" width="100%"></p>

**a.** Strong cohort:
- represented: 0.052 nats;
- attributable beyond the covariates: 0.061 (covariates alone 0.078, combined 0.138);
- usable by a held-out ridge readout: 0.045 (AUC 0.68).

Weak cohort: represented 0.008, attributable 0.011 and usable −0.002 (AUC 0.51).

### 5. How many patients a claim needs

`learning_curve()` refits the estimator on label-stratified subsamples of n patients. At each n it also fits a null on permuted labels and a held-out usable readout. A claim is **resolved** at n when the lower 90 % bound exceeds the null mean plus ε_rep, and N_evidence is the first such n.

```r
lc <- learning_curve(X, y, "poly", 20, n_grid = c(50, 100, 150, 200, 250, 300), B = 40,
                     eps_rep = eps$bound, cores = 20,
                     full = list(estimate = fits$poly$estimate, se = hs$se,
                                 null_mean = nulls$poly$mean))
attr(lc, "N_evidence")
```

**b.** In the strong cohort, the estimate stays well above the null at every size, but its interval reaches the threshold until the full cohort. N_evidence is therefore 400: the claim needs every patient.

**c.** In the weak cohort, the estimate follows the null down as n grows and is never resolved.

An unresolved Record means insufficient evidence under this representation, patient budget and protocol. It is not evidence that the state leaves no trace in the tissue. Here the simulated shift is real, but too weak to resolve at 400 patients.

### 6. Where on the slide the evidence lies, and how much tissue is needed

Under mean pooling with linear PCA, a patient's fitted logit is an affine function of the mean of its tile embeddings. Every tile therefore contributes an exact share. `tile_decomposition()` returns each tile's contribution and its **specific surprise**: the KL divergence between the tile's as-if prediction and the prevalence (DeWeese and Meister, 1999). Tile contributions reproduce the patient logits to 1.3 × 10⁻¹⁵.

`tissue_burden()` recomputes the estimate from m random tiles per slide.

```r
td <- tile_decomposition(sim$tiles, sim$patient, sim$y, d = 20)
plot_tile_map(td, slide = 1, coords = sim$coords[[1]])
bd <- tissue_burden(td)
```

<p align="center"><img src="man/figures/fig3_tiles.png" width="100%"></p>

The top row shows the simulated tissue, with tumour in purple. The bottom row shows the specific surprise of each tile.

- In the confidently called positive patient (P0110, fitted p = 0.81), the high-surprise tiles sit on the tumour region.
- In the positive patient called near the prevalence (P0130, p = 0.39), most tiles carry little evidence. Its tumour trace is weaker than the noise in its tiles.
- The negative patient (P0300) shows moderate surprise everywhere, in the direction of the negative state.

Right: the plug-in estimate reaches 90 % of its whole-slide value with 15 tiles per slide and 95 % with 30. With only two to five tiles, the noisy tile means can predict worse than the prevalence. The ratio then falls below zero, and the axis is cut at −0.5.

### 7. The Recoverability Record

```r
rec <- recoverability_record(
  provenance = list(model = "simulated encoder (D = 64)", pooling = "mean",
                    projection = "poly", d = 20, unit = "patient"),
  fit = fits$poly, null = nulls$poly, uncertainty = hs, eps_rep = eps$bound,
  attributable = inc, usable = use, burden = bd, curve = lc)
rec
#> Recoverability Record
#>   represented    0.0516 nats (8% of H(Y) = 0.631), k* = 6, n = 400
#>   null           0.0017 +/- 0.0029 (p = 0.010) | eps_rep 0.0040
#>   90% interval   0.0128 to 0.0903 (s.e. 0.0236)
#>   attributable   0.0605
#>   usable         0.0445 (AUC 0.68)
#>   tiles N90/N95  15 / 30
#>   status         resolved (N_evidence 400)
record_to_json(rec, "record.json")
```

The weak cohort's Record reads `status unresolved (N_evidence not reached)`. The full JSON Records of both cohorts are in [`tools/`](tools). The 12 Records of the study and its negative control are in [`results/record`](results/record).

### 8. A pooled estimate can hide a context

Suppose a cohort mixes a context A, where the representation predicts the state, with a context B, where it does not, and the estimator cannot see the context. If the context S is independent of Z, then

I(Z;Y) ≤ I(Z;Y | S) = π_A I(Z;Y | S = A).

`simulate_hidden_context()` builds this design with Monte Carlo truths.

<p align="center"><img src="man/figures/fig4_hidden_context.png" width="55%"></p>

The information within context A stays at 0.169 nats. As the share of context A, π_A, falls from 1 to 0.25:

- the pooled information falls to 0.0095;
- the PathMI estimate (mean of 10 replicates) follows the pooled value, from 0.171 to 0.003.

A small pooled estimate is therefore compatible with a subgroup in which the state is well recovered.

## Functions

| Step | Functions |
|---|---|
| Representation | `pool_tiles()`, `pool_cohort()`, `fit_projection()`, `apply_projection()`, `project_embeddings()` |
| Estimate | `pathmi_fit()`; engine `subho_MI_logistic()`, `selection.logistic()` |
| Null and selection | `draw_permutations()`, `pathmi_null()`, `select_representation()` |
| Resolution and uncertainty | `replicate_refit()`, `sign_flip_refit()`, `null_block_refit()`, `reproducibility_bound()`, `half_sampling()` |
| Counterparts | `incremental_recoverability()`, `usable_recoverability()` |
| Evidence | `learning_curve()`, `evidence_status()` |
| Tiles | `tile_decomposition()`, `tissue_burden()` |
| Record | `recoverability_record()`, `record_to_json()` |
| Plots | `plot_null()`, `plot_walk()`, `plot_learning_curve()`, `plot_tile_map()`, `plot_burden()` |
| Simulation | `simulate_cohort()`, `simulate_hidden_context()`, `simulate_estimator_design()` |
| Comparison estimators | `doe_linear()` (Donsker–Varadhan bound, linear critic), `run_PathMI_LASSO()`, `run_PathMI_OLS()` |
| Earlier interface | `pathMI_pipeline()`, `permutation_test_logistic()`, `plot_MI_null()`, `plot_likelihood_surface()`, `v_beta_sensitivity()` |

Every function has a help page (`?pathmi_fit`, `?learning_curve`, …). The sampler is silent by default; `options(PathMI.verbose = TRUE)` prints its progress.

## Repository layout

| Path | Content |
|---|---|
| `R/`, `man/`, `tests/` | The package and its tests. |
| `tools/` | The script that produces the README figures, and its outputs (numbers, simulated Records as JSON). |
| `results/` | Patient-level result tables of the study: the 192-cell grid, Records, usable recoverability, learning curves, robustness, burden, tile summaries and simulations. It holds no embeddings. See [`results/README.md`](results/README.md). |
| `study/` | The scripts that produced the manuscript's analyses, as run. They read Trident features that are not distributed here. See [`study/README.md`](study/README.md). |

## Data

Whole-slide images and molecular annotations for TCGA COAD, READ and PRAD are available from the Genomic Data Commons (https://portal.gdc.cancer.gov). The study's tile embeddings were extracted with Trident from the public UNI v2, Virchow2, H-optimus-1 and CONCH v1.5 weights. They are not redistributed here and can be regenerated under each model's licence.

## License

MIT; see [`LICENSE.md`](LICENSE.md).

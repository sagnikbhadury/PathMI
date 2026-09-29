# Results deposited with the manuscript

Outputs of the study *Measuring molecular recoverability in pathology foundation models*
(TCGA colorectal and prostate cohorts; UNI v2, Virchow2, H-optimus-1 and CONCH v1.5
embeddings extracted with Trident at 20X). Estimates are in nats. Every file was produced by the
scripts in this repository's pipeline; the merged tables are what the figures and the numbers in
the paper are read from.

## `tables/`

| File | Content |
|---|---|
| `mi_summary_table.xlsx` | The 192-cell recoverability grid: state x model x pooling x projection, with the estimate, the permutation-null mean and s.d., the permutation p value and the condition number of the retained spectrum. |
| `best_kpca_summary.xlsx` | The 48 selected projections (one per state, model and pooling) with the max-T p value. |
| `incremental_mi_summary.xlsx` | Covariate-only, embedding-only and combined estimates and the incremental recoverability for the 48 selected configurations (base covariates: log tile count and stage or Gleason sum). |
| `prad_incremental_batch_purity/` | Prostate refits with scanning batch and ABSOLUTE purity added to the covariate set (`prad_incremental_batch_purity.csv`, `covariate_only_fits.csv`, `run.log`). |
| `sampling_burden.csv` | Plug-in estimate recomputed from m tiles per slide (m from 5 to 2500, 20 draws each) under the mean-pooling / linear-PCA arm of H-optimus-1. |
| `simulation_results.rds` | The 72-cell x 50-replicate simulation comparing PathMI, LASSO-CV, OLS and the linear Donsker-Varadhan bound. |
| `per_run/<cohort>_<model>/` | The per-run outputs of the pipeline for each cohort and model: `kpca_mi_summary_*.csv` (16 projections with nulls), `kpca_selected_decomposition_*.csv`, `best_kpca_mi_summary_*.csv`, `incremental_mi_*.csv` and the full result object `kpca_subho_results_*.rds` (fits, null draws, selected coordinates). |

Model labels: `uni_v2`, `virchow2`, `hoptimus1`, `conch_v15`. Cohort labels: `crc` (COAD + READ,
state MSI) and `prad` (TCGA-PRAD, states ERG and ETS).

## `tile_scores/`

Tile-level decomposition of the H-optimus-1 mean-pooling / linear-PCA fit.

| File | Content |
|---|---|
| `<cohort>/slide_summary_<cohort>_hoptimus1.csv` | One row per slide: label, fitted probability, posterior logit and its s.d., tile count, share of tiles called informative, per-patient difference-of-entropies term. |
| `<cohort>/fit_summary_<cohort>_hoptimus1.csv` | Selected model size, estimate, plug-in estimate with its 90 percent posterior interval, null band used by the tile rule. |
| `<cohort>/fits_<cohort>_hoptimus1.rds` | The refitted selected models (MAP coefficients, Hessian, projection parameters) needed to score new tiles. |
| `<cohort>/tile_maps_<cohort>_hoptimus1.log` | Run log with the numerical checks (projection rebuilt to machine precision, slide logits reproduced from tile averages). |

The per-tile score files, `tile_scores_crc_hoptimus1.rds` (4.78 million tiles, 173 MB) and
`tile_scores_prad_hoptimus1.rds` (4.11 million tiles, 300 MB), exceed GitHub's file-size limit
and are attached to the release **tile-scores-v1** of this repository:
<https://github.com/sagnikbhadury/PathMI/releases/tag/tile-scores-v1>. Each is a data frame with
one row per tile: slide, tile coordinates at level 0, the per-tile contribution g(z) to the slide
logit, the specific information (posterior mean) and its 5 and 95 percent posterior quantiles for
each state.

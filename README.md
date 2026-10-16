# pFFPLS simulation study

Simulation code for the paper **"Penalized function-on-function partial least-squares
regression"**. The study
compares the proposed penalized function-on-function PLS (**pFFPLS**) with four
competitors for function-on-function regression, where a response curve $Y(q)$ is
predicted from a predictor curve $X(p)$, both observed on $[0,1]$.

The simulation design is described in the paper. This README explains how to install
the dependencies, run the study and rebuild the figures.

## Methods compared

| Label | Method |
|---|---|
| **pFFPLS** | Penalized function-on-function PLS (proposed; R package `penFoFPLS`) |
| **FFPLS** | Unpenalized FFPLS with a fixed number of basis functions |
| **FFPLS_OB** | Unpenalized FFPLS with the numbers of basis functions for $X$ and $Y$ chosen by cross-validation |
| **pFFR_I** | Penalized function-on-function regression of Ivanescu et al. (2015), `refund::pffr` |
| **pFFR_RS** | Penalized function-on-function regression of Ramsay & Silverman, `fda.usc::fregre.basis.fr` |

All five methods are fitted to the same simulated dataset in every repetition.

## Installation

The code was run with R 4.2.3. Install the `penFoFPLS` package from GitHub:

```r
# install.packages("devtools")
devtools::install_github("hhroig/penFoFPLS", dependencies = TRUE)
```

and the remaining CRAN packages:

```r
install.packages(c(
  # simulation
  "fda", "fda.usc", "refund", "pls", "dplyr", "tibble", "reshape2",
  "doParallel", "foreach", "caret", "tictoc",
  # figures
  "tidyverse", "ggplot2", "ggpubr", "gridExtra", "viridis", "scales",
  "plotly", "writexl"
))
```

## Files

| File | Purpose |
|---|---|
| `run_simulation.R` | Configuration and entry point: runs the study and writes the results |
| `simulations_buried_bands.R` | Data-generating process and the per-repetition fit of all five methods (sourced by `run_simulation.R`) |
| `cv_penalties_fregre.basis.fr.R` | Cross-validation of the pFFR_RS penalties |
| `predict_fregre_fr.R` | Predictions for pFFR_RS |
| `make_plots.R` | Builds the `results_plots/` folder for a results folder |
| `compare_methods_fofr_with_ivanescu_ramsay_silverman.R` | Boxplots, error curves, fR², timings, penalties, Excel summaries, 3D surfaces (used by `make_plots.R`) |
| `plot_beta_heatmaps.R` | 2D heatmaps of the true and estimated coefficient surfaces (used by `make_plots.R`) |

All scripts are run from the project root.

## Running the study

```sh
Rscript run_simulation.R            # all 100 repetitions in one process
Rscript run_simulation.R 1 100 8    # repetitions 1-100, 8 cores
```

The arguments are `[rep_start rep_end] [n_cores] [band_sd]`. The defaults reproduce the
study reported in the paper:

| Setting | Value |
|---|---|
| Repetitions | 100 |
| Sample size | 100 training + 50 test curves per repetition |
| Grids | 100 equally spaced points for $p$ and for $q$ |
| Predictor noise | $\sigma_p = 1.5$ |
| Band width | `band_sd` = 0.035 |
| Response noise | $R^2 = 0.9$ |
| Basis | $K = L = 40$ cubic B-splines (FFPLS_OB searches 10–40) |
| Components | up to 10 (PLS methods) |
| Cross-validation | 5-fold |

The remaining settings (penalty grids, etc.) are the arguments of
`main_buried_bands_call()` in `run_simulation.R`.

**Seeds and parallel runs.** Repetition $r$ always uses `set.seed(2026 + r)`, so a
repetition gives the same result whichever process runs it. The repetitions can therefore
be split across several processes writing to the same folder.

**Computing cost.** pFFR_I is single-threaded and dominates the cost: about 35 minutes
per repetition at $K = L = 40$, with a peak of about 2 GB of RAM per process. Several
processes with a few cores each are much faster than one process with many cores. For
example, on an 8-core machine with 32 GB of RAM:

```sh
Rscript run_simulation.R  1  25 2 &
Rscript run_simulation.R 26  50 2 &
Rscript run_simulation.R 51  75 2 &
Rscript run_simulation.R 76 100 2 &
```

This takes about 17–18 hours for the 100 repetitions. With more processes than about
4 per 32 GB, memory runs out.

**Output.** Results are written to
`results_simulations/buried_bands/rep100_K40L40_sdX1.5_band0.035_seed2026/`, as one set
of `.Rds` files per repetition `<r>`:

| File | Contents |
|---|---|
| `final_models_rep_<r>_*.Rds` | IMSE of $\hat\beta$ and test IMSE of $\hat Y$, per method and number of components |
| `betas_rep_<r>_*.Rds` | Estimated and true surfaces $\beta(q,p)$ on the evaluation grid |
| `cves_rep_<r>_*.Rds` | Cross-validation errors |
| `best_lambdas_rep_<r>_*.Rds` | Selected penalties (pFFPLS, pFFR_RS) |
| `best_num_bases_FFPLS_rep_<r>_*.Rds` | Selected basis sizes (FFPLS_OB) |
| `fR2s_rep_<r>_*.Rds` | Functional $R^2$, training and test |
| `computation_times_rep_<r>_*.Rds` | Computing times |

## Figures

```sh
Rscript make_plots.R                                   # the default results folder
Rscript make_plots.R results_simulations/buried_bands/<folder>/
```

This writes `results_plots/` inside the results folder:

| Subfolder | Contents |
|---|---|
| `IMSEs_CVEs_buried_bands/` | IMSE($\hat\beta$) boxplots for each number of components $H$; training CVE and test IMSE($\hat Y$) curves |
| `2DBetaHeatmaps_buried_bands/` | Heatmaps of the true surface and of each method's mean estimate: all methods in one panel with a shared colour scale (`all_methods_ncomp<H>`), and one file per method with its own scale (`owned_scale/`) |
| `3DBeta2D_buried_bands/` | 3D views of the same surfaces (`shared_scale/`, `owned_scale/`) |
| `fR2/` | Functional $R^2$ |
| `computation_times/` | Computing times |
| `best_penalties/` | Selected penalties |
| `IMSEs_CVEs_Excel/` | Summary tables |

pFFR_I and pFFR_RS do not use components, so their results are the same for every $H$.

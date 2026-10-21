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
| `plot_example_curves.R` | Example predictor and response curves, without and with noise |

All scripts are run from the project root.

## Running the study

The scenarios reported in the paper are predefined in `run_simulation.R`. Each one has
100 repetitions and writes to its own folder under `results_simulations/buried_bands/`:

| Scenario | Setting | Results folder |
|---|---|---|
| `main` | $K = L = 40$, all five methods, measurement noise on $X$ ($\sigma_p = 1.5$) | `rep100_K40L40_sdX1.5_band0.035_seed2026/` |
| `main_noiseless` | as `main`, without measurement noise ($\sigma_p = 0$) | `rep100_K40L40_sdX0_band0.035_seed2026/` |
| `lowdim` | $K = L = 7$, without FFPLS_OB, up to 7 components, wider penalty grids | `rep100_K7L7_sdX1.5_band0.035_seed2026_noOB/` |
| `lowdim_noiseless` | as `lowdim`, $\sigma_p = 0$ | `rep100_K7L7_sdX0_band0.035_seed2026_noOB/` |
| `band0.02` | as `main`, with narrower bands (`band_sd` = 0.02) | `rep100_K40L40_sdX1.5_band0.02_seed2026/` |

```sh
Rscript run_simulation.R main              # all 100 repetitions in one process
Rscript run_simulation.R main 1 100 8      # repetitions 1-100, 8 cores
Rscript run_simulation.R lowdim_noiseless  # another scenario
```

The arguments are `<scenario> [rep_start rep_end] [n_cores]`. From R:

```r
source("run_simulation.R")
run_scenario("lowdim", rep_starts = 1, rep_end = 10, n_cores = 4)
run_scenario("main", dry_run = TRUE)   # only returns the output folder
```

Settings shared by all scenarios:

| Setting | Value |
|---|---|
| Repetitions | 100 |
| Sample size | 100 training + 50 test curves per repetition |
| Grids | 100 equally spaced points for $p$ and for $q$ |
| Band width | `band_sd` = 0.035 (except `band0.02`) |
| Response noise | $R^2 = 0.9$ |
| Components | up to 10 at $K = L = 40$, up to 7 at $K = L = 7$ (PLS methods) |
| FFPLS_OB basis search | 10–40 (scenarios with $K = L = 40$) |
| Cross-validation | 5-fold |

All other settings (penalty grids, basis sizes, etc.) are arguments of
`main_buried_bands_call()`. The scenarios list only the arguments that differ from its
defaults, and `run_scenario()` accepts any of them as an override.

**Seeds and parallel runs.** Repetition $r$ always uses `set.seed(2026 + r)` in every
scenario. Consequences:

- A repetition gives the same result whichever process runs it, so the repetitions can be
  split across several processes writing to the same folder.
- All scenarios share the same noise-free curves and responses. The noiseless scenarios
  differ from their noisy counterparts only in the noise added to $X$.

**Computing cost.** At $K = L = 40$, pFFR_I is single-threaded and dominates the cost:
about 35 minutes per repetition, with a peak of about 2 GB of RAM per process. Several
processes with a few cores each are much faster than one process with many cores. For
example, on an 8-core machine with 32 GB of RAM:

```sh
Rscript run_simulation.R main  1  25 2 &
Rscript run_simulation.R main 26  50 2 &
Rscript run_simulation.R main 51  75 2 &
Rscript run_simulation.R main 76 100 2 &
```

This takes about 17–18 hours for the 100 repetitions. With more processes than about
4 per 32 GB, memory runs out. The $K = L = 7$ scenarios are much cheaper (a few minutes
per repetition, little memory).

**Windows path length.** Windows cannot open paths longer than 260 characters.
`run_simulation.R` warns if an output folder is too deep for the figure files. In that
case, move the project to a shorter path or pass a shorter `folder_name`.

**Output.** Each results folder holds one set of `.Rds` files per repetition `<r>`:

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
Rscript make_plots.R                                   # the main scenario's folder
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

**Example curves.** To draw a few predictor curves without and with measurement noise,
and their responses, from repetition 1 of the study, run:

```sh
Rscript plot_example_curves.R [rep] [X_sd_error]
```

It writes several layout options (`option_A`–`option_E`, PNG and PDF) to
`results_simulations/buried_bands/example_curves/`. The option C grid is also saved as EPS (`option_C_grid_2x2.eps`), and option C is also
written as four separate EPS files
(`option_C_X_noisefree`, `option_C_X_observed`, `option_C_Y_noisefree`,
`option_C_Y_observed`); the two X(p) plots share one vertical scale and the two Y(q)
plots another.

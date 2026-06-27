# pFFPLS Simulation Study

Simulation code for the paper **"Penalized Function-on-Function PLS"** (3rd SMJ revision).
The study compares the proposed penalized FFPLS method against several competitors on
function-on-function regression with three different coefficient surfaces.

## Methods compared

| Label | Description |
|-------|-------------|
| **pFFPLS** | Proposed penalized Function-on-Function PLS |
| **FFPLS** | Unpenalized FFPLS (fixed basis count) |
| **FFPLS_OB** | Unpenalized FFPLS with CV-selected basis count (Setting 2 only) |
| **pFFR_I** | Ivanescu et al. penalized FFR via `refund::pffr` |
| **pFFR_RS** | Ramsay & Silverman penalized FFR via `fda.usc` |

## Coefficient functions (betas)

All three betas are defined on [0,1]×[0,1] where p is the predictor argument and
q is the response argument. They cover a spectrum of estimation difficulty:

| id | Code name | Formula | Structure | Difficulty |
|----|-----------|---------|-----------|------------|
| 1 | `cos_sin` | cos(2π p) · sin(2π q) | Separable, 4-lobe, 1 cycle per axis | Low |
| 2 | `sin_sum` | sin(3π(p + q)) | Anti-diagonal, 3 oscillation cycles | Medium |
| 3 | `cos_sum` | 0.5·cos(0.5π(p+q)) + cos(1.5π(p+q)) + 2·cos(2.5π(p+q)) | Multi-frequency anti-diagonal | High |

**cos_sin** is separable, so a single PLS component captures most variation — all
methods should perform similarly here.
**sin_sum** requires more components due to its anti-diagonal structure and serves
as an intermediate test.
**cos_sum** is the hardest case: the high-frequency dominant term is where the
roughness penalty in pFFPLS is most beneficial, and where unpenalized methods tend
to over-fit.

## Simulation settings

| Setting | K = L (bases) | FFPLS_OB | max nComp |
|---------|--------------|----------|-----------|
| 1 | 7 | No | 6 |
| 2 | 40 | Yes (CV-selected) | 8 |
| 3 | 40 | No | 8 |

Each setting is run twice: once with clean X observations (`X_sd_error = 0`) and
once with added observation noise (`X_sd_error = 0.2`). Settings 1 and 2 are
currently active; Setting 3 is commented out in `main_simulations_call.R`.

## How to run

Install the package dependency first:

```r
devtools::install_github("hhroig/penFoFPLS", dependencies = TRUE)
```

Then run the full study (simulations + plots) from the project root:

```r
source("main_simulations_call.R")
```

To re-run only the comparison/plotting step on already-saved results:

```r
library(tidyverse)
library(openxlsx)
library(reshape2)

source("compare_methods_fofr_with_ivanescu_ramsay_silverman.R", local = TRUE)
compare_methods_fun(input_folder = "results_simulations/set1_rep3_pen9_K7L7/")
```

## File structure

```
main_simulations_call.R                       # Entry point — settings and calls
simulations_fofr_v2_with_ivanescus_...R       # Simulation loop (sourced internally)
compare_methods_fofr_with_ivanescu_...R       # Plotting and summary (sourced internally)
cv_penalties_fregre.basis.fr.R                # CV wrapper for R&S method
predict_fregre_fr.R                           # Predict wrapper for R&S method
docs/                                         # Supporting documentation
```

## Saved results (per repetition and beta)

Each simulation run writes the following `.Rds` files into the setting-specific
subfolder of `results_simulations/` (e.g. `results_simulations/set1_rep3_pen9_K7L7/`):

| File pattern | Contents |
|---|---|
| `final_models_rep_*_beta_*.Rds` | IMSE on training set and IMSE on validation set, per method and nComp |
| `cves_rep_*_beta_*.Rds` | Cross-validation error used for nComp / penalty selection |
| `R2s_rep_*_beta_*.Rds` | Pointwise R²(q) vector (one value per response time point), per method and nComp |
| `fR2s_rep_*_beta_*.Rds` | Scalar functional R² pooled over all observations and time points |
| `betas_rep_*_beta_*.Rds` | Estimated and true beta surfaces on the evaluation grid |
| `best_lambdas_rep_*_beta_*.Rds` | CV-optimal penalty parameters for pFFPLS and pFFR_RS |
| `computation_times_rep_*_beta_*.Rds` | Wall-clock time per method |

## Results plots

The compare script writes everything to `results_simulations/<setting>/results_plots/`.
Subdirectories and their contents:

| Subfolder | What to look at |
|---|---|
| `IMSEs_CVEs_Excel/` | Excel summaries: mean/median/sd IMSE (training + validation) and CVE by method, beta, and nComp. **Start here for numeric comparisons.** |
| `IMSEs_CVEs_{beta}/` | Per-beta IMSE box/line plots (linear and log scale) for training and validation sets. Shows how each method degrades as nComp grows or changes. |
| `fR2/` | Functional R² plots and `summary_fR2.xlsx`. The scalar fR² captures overall fit quality; compare training vs test to detect over-fitting. |
| `R2/` | Pointwise R²(q) curves: one value per response time point. Useful for identifying where each method fits better or worse along the response domain. |
| `best_penalties/` | Selected penalty values for pFFPLS and pFFR_RS across reps and betas. Useful for checking that the penalty grid is wide enough. |
| `computation_times/` | Elapsed time per method (linear and log scale). |
| `mean_beta_{beta}/` | Mean estimated beta surface vs the true surface, one panel per method. The key visual for assessing shape recovery. |

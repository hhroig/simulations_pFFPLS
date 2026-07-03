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

See `docs/new_betas.md` for detailed motivation and shape descriptions.

## Simulation settings

| Setting | K = L (bases) | FFPLS_OB | max nComp | Status |
|---------|--------------|----------|-----------|--------|
| 1 | 7 | No | 6 | Complete (30 reps) |
| 2 | 40 | Yes (CV-selected) | 8 | Complete (30 reps) |
| 3 | 40 | No | 8 | Not run (commented out) |

Each setting is run twice: once with clean X observations (`X_sd_error = 0`) and
once with added observation noise (`X_sd_error = 0.2`), indicated by the `e` suffix
in the output folder name (e.g. `set1e_...`).

## How to run

Install the package dependency first:

```r
devtools::install_github("hhroig/penFoFPLS", dependencies = TRUE)
```

Then run the full study (simulations + comparison plots) from the project root:

```r
source("main_simulations_call.R")
```

The global parameters at the bottom of `main_simulations_call.R` control which
settings are active. Settings 1 and 2 calls are currently present; Setting 3 calls
are commented out. To re-run only the comparison and plotting step on
already-saved results:

```r
library(tidyverse)
library(openxlsx)
library(reshape2)

source("compare_methods_fofr_with_ivanescu_ramsay_silverman.R", local = TRUE)
compare_methods_fun(input_folder = "results_simulations/set1_rep30_pen100_K7L7/")
```

Replace the folder path with any completed setting folder to regenerate its plots.

## File structure

```
main_simulations_call.R                       # Entry point — settings and active calls
simulations_fofr_v2_with_ivanescus_...R       # Simulation loop (one rep per beta)
compare_methods_fofr_with_ivanescu_...R       # Aggregation, plotting, and Excel export
cv_penalties_fregre.basis.fr.R                # CV penalty grid wrapper for R&S method
cv_bases_fregre.basis.fr.R                    # CV basis-count wrapper for FFPLS_OB
predict_fregre_fr.R                           # Predict wrapper for R&S method
docs/new_betas.md                             # Motivation and formulas for the 3 betas
docs/new_fR2.md                               # Definition and rationale for fR² metric
results_simulations/                          # All output (one subfolder per setting)
```

Output folder names encode the setting parameters:
`set{setting}{e}_rep{total_reps}_pen{n_lambdas²}_K{KK}L{LL}/`

For example, `set2e_rep30_pen100_K40L40/` is Setting 2, noisy X, 30 reps, 10×10
penalty grid, 40 bases.

## Saved RDS files (per repetition and beta)

Each simulation run writes the following `.Rds` files into the setting subfolder:

| File pattern | Contents |
|---|---|
| `final_models_rep_*_beta_*.Rds` | IMSE on training set and IMSE on validation set, per method and nComp |
| `cves_rep_*_beta_*.Rds` | Cross-validation error used for nComp / penalty selection |
| `R2s_rep_*_beta_*.Rds` | Pointwise R²(q) vector per method and nComp |
| `fR2s_rep_*_beta_*.Rds` | Scalar functional R² (train and test), per method and nComp |
| `betas_rep_*_beta_*.Rds` | Estimated and true beta surfaces on the evaluation grid |
| `best_lambdas_rep_*_beta_*.Rds` | CV-optimal penalty parameters for pFFPLS and pFFR_RS |
| `computation_times_rep_*_beta_*.Rds` | Wall-clock time per method |

## Results plots

The compare script writes everything to `results_simulations/<setting>/results_plots/`.
Each subfolder is described below, roughly in order of usefulness for interpretation.

---

### `IMSEs_CVEs_Excel/`

**Start here for numeric comparisons.**  
Contains three Excel workbooks summarising all 30 reps:

| File | Contents |
|---|---|
| `summary_imse_beta.xlsx` | Mean / median / sd of beta-surface IMSE, by method, beta, and nComp |
| `summary_imse_val_y.xlsx` | Mean / median / sd of prediction error on the test Y, same breakdown |
| `summary_training_cves.xlsx` | Mean / median / sd of training CVE, same breakdown |

---

### `IMSEs_CVEs_{beta}/`

Per-beta IMSE and CVE plots (one subfolder per beta: `cos_sin`, `sin_sum`, `cos_sum`).

| File pattern | What to look at |
|---|---|
| `imse_beta_{beta}.(pdf\|png)` | Boxplots of IMSE across reps, all methods, all nComp overlaid |
| `log_imse_beta_{beta}.*` | Same on log scale — useful when methods span orders of magnitude |
| `imse_beta_{beta}_ncomp{N}.*` | IMSE at a fixed nComp — compare methods head-to-head |
| `log_imse_beta_{beta}_ncomp{N}.*` | Log scale version of the above |
| `test_imseY_{beta}.*` | Prediction IMSE on test Y, all nComp |
| `log_test_imseY_{beta}.*` | Log scale version |
| `test_PLS_imseY_{beta}.*` | Test IMSE for the PLS-family methods only (cleaner view) |
| `train_cveY_val_imseY_beta_{beta}.*` | Side-by-side panel: training CVE and validation IMSE — use this to check whether CVE is a good proxy for out-of-sample error |
| `train_PLS_cveY_{beta}.*` | Training CVE for PLS methods only |
| `train_log_cveY_{beta}.*` | Log scale training CVE |

---

### `fR2/`

Scalar functional R² (a single number per rep summarising overall goodness of fit).
Higher is better; compare train vs test to detect over-fitting.

| File pattern | What to look at |
|---|---|
| `fR2_beta_{beta}.*` | Train + test fR² by method, all nComp in one plot |
| `fR2_train_beta_{beta}.*` | Training fR² only |
| `fR2_test_beta_{beta}.*` | Test fR² only — the primary model selection diagnostic |
| `fR2_train_{beta}_ncomp{N}.*` | Training fR² at a fixed nComp |
| `fR2_test_{beta}_ncomp{N}.*` | Test fR² at a fixed nComp — compare methods at equal complexity |
| `summary_fR2.xlsx` | Numeric summary (mean / median / sd) for all the above |

---

### `R2/`

Pointwise R²(q) curves — one value per response time point q, averaged across reps.
Shows where along the response domain each method fits better or worse.

| File pattern | What to look at |
|---|---|
| `R2_beta_{beta}_nComp{N}.*` | Train + test R²(q) curves for all methods at fixed nComp |
| `R2_Training_beta_{beta}_nComp{N}.*` | Training R²(q) only |
| `R2_Test_beta_{beta}_nComp{N}.*` | Test R²(q) only |
| `R2_zoom_*` | Same plots zoomed into [0,1] on the y-axis for clearer differentiation |
| `R2_rough_{beta}_nComp{N}.*` | Highlights methods with rougher or more variable R²(q) profiles |
| `R2_rough_zoom_*` | Zoomed version of the rough comparison |

---

### `mean_beta_{beta}/`

2D heatmap (raster + contour) of the mean estimated beta surface across all reps,
one panel per method, at each nComp. Use these to assess how well each method
recovers the shape of the true coefficient surface on average.

| File pattern | Contents |
|---|---|
| `2d_mean_beta_{beta}_nComp_{N}.(pdf\|png)` | All methods + true beta side by side at nComp = N |

---

### `3DBeta2D_{beta}/`

3D perspective renders of individual beta surface estimates (one plot per method ×
nComp combination), produced with base R `persp()`. Each file shows the surface
viewed from a fixed angle, allowing direct visual comparison of shape and amplitude.
The true beta surface is included for reference.

| File pattern | Contents |
|---|---|
| `{method}_beta{beta}_ncomp{N}.(pdf\|eps\|png)` | Perspective plot for a single method at nComp = N |
| `True_beta{beta}_ncomp{N}.*` | True coefficient surface at the same nComp (for reference) |

Methods present: `FFPLS`, `FFPLS_OB` (Setting 2 only), `pFFPLS`, `pFFR_I`, `pFFR_RS`, `True`.

**Reading tip:** compare the `pFFPLS` and `FFPLS` surfaces — in `cos_sum` at high nComp,
pFFPLS should remain smoother while FFPLS starts to over-fit the high-frequency term.

---

### `best_penalties/`

2D heatmap plots of the CV-selected penalty values (λ_X, λ_Y) for pFFPLS and
pFFR_RS across all reps, one plot per beta. Use these to verify that the optimal
penalties are not hitting the boundary of the search grid (which would indicate the
grid needs extending).

| File pattern | Contents |
|---|---|
| `2d_best_lambdas_{beta}.(pdf\|png)` | Selected (λ_X, λ_Y) pairs across reps |

---

### `computation_times/`

Wall-clock time per method, useful for reporting practical cost differences.

| File pattern | Contents |
|---|---|
| `computaion_times_beta{beta}.(pdf\|png)` | Boxplot of elapsed time per method (linear scale) |
| `log_computaion_times_beta{beta}.*` | Same on log scale |
| `all_computation_times.xlsx` | Numeric summary (mean / median / sd) by method and beta |

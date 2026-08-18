# pFFPLS Simulation Study

Simulation code for the paper **"Penalized Function-on-Function PLS"** (3rd SMJ
revision). The study compares the proposed penalized FFPLS method (pFFPLS) against
several competitors on function-on-function regression: predicting an entire
response curve $Y(q)$ from an entire predictor curve $X(p)$, both observed on
$[0,1]$.

This README explains, at a level that shouldn't require re-reading the code: what's
being compared, what the coefficient surfaces mean and why there are four of them,
why the predictor curves are generated two different ways, and how the settings are
put together. For the underlying math and package-level details, see the `docs/`
folder referenced throughout.

## Methods compared

| Label | Description |
|-------|-------------|
| **pFFPLS** | Proposed penalized Function-on-Function PLS — the paper's contribution |
| **FFPLS** | Unpenalized FFPLS, fixed basis count |
| **FFPLS_OB** | Unpenalized FFPLS with a CV-selected basis count (Setting 2 only — see below) |
| **pFFR_I** | Ivanescu et al.'s penalized FFR, via `refund::pffr` |
| **pFFR_RS** | Ramsay & Silverman's penalized FFR, via `fda.usc` |

All five methods are fit to the *same* simulated data in every repetition, so
differences reflect the methods themselves, not different random draws.

## Coefficient functions (betas)

The object every method tries to recover is the coefficient surface $\beta(p,q)$ in
$Y(q) = \int X(p)\,\beta(p,q)\,dp + \varepsilon(q)$. Four surfaces are used, each
defined on $[0,1]\times[0,1]$ ($p$ = predictor argument, $q$ = response argument),
chosen to span a spectrum of shapes and estimation difficulty rather than to
resemble any particular application:

| id | Code name | Formula | Shape | Difficulty |
|----|-----------|---------|-------|------------|
| 1 | `cos_sin` | $\cos(2\pi p)\sin(2\pi q)$ | Separable, 4-lobe checkerboard, 1 cycle per axis | Low |
| 2 | `sin_sum` | $\sin(3\pi(p+q))$ | Anti-diagonal stripes, 3 cycles | Medium |
| 3 | `cos_sum` | $0.5\cos(0.5\pi(p{+}q)) + \cos(1.5\pi(p{+}q)) + 2\cos(2.5\pi(p{+}q))$ | Anti-diagonal, 3 superposed frequencies | High |
| 4 | `sin_sum_cos_diff` | $\sin(3\pi(p{+}q)) + \cos(3\pi(p{-}q))$ | Mixed diagonal **and** anti-diagonal | High (different kind) |

What each one is actually testing:

- **`cos_sin`** is a product of a function of $p$ alone and a function of $q$ alone
  — the simplest possible non-trivial shape. A single PLS component captures almost
  all of it, so this is the sanity-check case: every method should do reasonably
  well here, and any method that doesn't is oversmoothing (or otherwise misbehaving)
  on an easy target.
- **`sin_sum`** depends only on $p+q$, so it's constant along anti-diagonal lines —
  visually, parallel diagonal stripes. It can't be captured by one component the way
  `cos_sin` can, so it needs more of the model's capacity and serves as an
  intermediate-difficulty case.
- **`cos_sum`** is also a function of $p+q$ only, but it's a sum of three different
  oscillation frequencies with growing amplitude, so the highest-frequency term
  dominates and the surface is genuinely rough. This is where a roughness penalty
  should earn its keep: unpenalized methods are expected to overfit the fast
  wiggle, while pFFPLS should track it without adding spurious noise.
- **`sin_sum_cos_diff`** adds a *diagonal* wave, $\cos(3\pi(p-q))$, on top of
  `sin_sum`'s anti-diagonal one. The first three betas are all, in effect,
  functions of a single (possibly rotated) coordinate — this one genuinely isn't:
  it has independent structure along both diagonals at once, a qualitatively
  different, two-directional kind of difficulty. It's only used together with the
  `fourier_decay` predictor process described next.

## How the predictor curves $X(p)$ are generated

Two different random processes are used to generate $X$, controlled by the
`X_process` argument:

- **`uniform_bspline`** (the default, and the one used in the paper): $X$ is a
  random combination of 20 cubic B-spline basis functions with independent
  Uniform$(-1,1)$ coefficients. Every one of the 20 "directions" contributes about
  equally, so the resulting curves are moderately rough with no dominant scale.
- **`fourier_decay`**: $X$ is a sum of 10 sine/cosine harmonic pairs with
  independent Normal coefficients whose variance decays like $1/r^4$. The first
  couple of harmonics dominate, so these curves are much smoother and effectively
  low-rank — closer to what's typically simulated elsewhere in the functional-data
  literature.

The reason both exist: every method/beta combination above always shares the *same*
predictor process by default, so any performance difference is attributable purely
to $\beta$'s shape — there's no way to ask whether pFFPLS's penalty is only useful
because the predictor itself happens to be rough. The `fourier_decay` runs (an
addition beyond what's in the current paper draft — see
`docs/new_setting_smooth_X_beta4.md`) answer exactly that question.

## Simulation settings

Each setting fixes the number of basis functions, $K = L$, used to *estimate* the
curves — as opposed to the $K=20$ basis functions always used to *generate* them
(see above). Comparing a setting with fewer estimation bases than generation bases
against one with more is what lets the study probe under- and over-parameterized
regimes:

| Setting | $K=L$ (estimation bases) | Regime | FFPLS_OB run? | max # PLS components |
|---------|---------------------------|--------|----------------|----------------------|
| **1** | 7 | Under-parameterized (fewer bases than used to generate the data) | No | 6 |
| **2** | 40 | Over-parameterized (more bases than used to generate the data) | Yes — its own CV search picks the basis count | 8 |

A third setting (K=L=40, identical to Setting 2 but without the FFPLS_OB
comparison) previously existed in the code but added no information beyond what
Setting 2 already provides — it's been dropped as redundant (still present,
commented out, at the bottom of `main_simulations_call.R` for reference).

Each setting is additionally run under **four** configurations, all currently
active and complete (30 repetitions each):

| Configuration | Predictor process | $X$ noise | Betas used | Output folder |
|---|---|---|---|---|
| Clean | `uniform_bspline` | none | `cos_sin`, `sin_sum`, `cos_sum` | `set{1,2}_rep30_pen100_K{7,40}L{7,40}/` |
| Noisy | `uniform_bspline` | $\sigma=0.2$ added to each observed point | `cos_sin`, `sin_sum`, `cos_sum` | `set{1,2}e_rep30_pen100_K{7,40}L{7,40}/` |
| Clean, smooth $X$ | `fourier_decay` | none | all 4 betas, incl. `sin_sum_cos_diff` | `set{1,2}_fourier_decay_rep30_pen100_K{7,40}L{7,40}/` |
| Noisy, smooth $X$ | `fourier_decay` | $\sigma=0.2$ | all 4 betas, incl. `sin_sum_cos_diff` | `set{1,2}e_fourier_decay_rep30_pen100_K{7,40}L{7,40}/` |

i.e. $2 \text{ settings} \times 2 \text{ noise levels} \times 2 \text{ predictor
processes} = 8$ runs in total, each 30 repetitions. Only the first two rows
(`uniform_bspline`) are described in the current paper draft; the `fourier_decay`
rows are the newer, exploratory addition explained above.

For every run: $n=150$ curves are simulated (100 for training, 50 for testing),
observed on a 100-point grid. Penalty parameters for pFFPLS and pFFR\_RS are
selected by 10-fold cross-validation over a $10\times10$ log-spaced grid
($\lambda \in [10^{-2}, 10^{12}]$ for pFFPLS, $[10^{-6}, 10^{8}]$ for pFFR\_RS);
FFPLS\_OB's basis-count search uses the same 10-fold CV over 10 candidate sizes
from 9 to 40. pFFR\_I uses `pffr`'s own default (REML-based) smoothing-parameter
selection.

## How to run

Install the package dependency first. The reviewed/fixed build currently lives in
this repo's local (git-ignored) `penFoFPLS/` folder — install from there to make
sure you get it, rather than whatever is on GitHub's default branch:

```r
devtools::install("penFoFPLS")
# or, once the fix is merged upstream:
# devtools::install_github("hhroig/penFoFPLS", dependencies = TRUE)
```

Then run the full study (simulations + comparison plots) from the project root:

```r
source("main_simulations_call.R")
```

The bottom of `main_simulations_call.R` lists the active `main_simulations_call()`
calls — currently all 8 configurations from the table above, using
`global_betas = c("cos_sin", "sin_sum", "cos_sum")` for the `uniform_bspline` runs
and `global_betas_smooth_X = c(global_betas, "sin_sum_cos_diff")` for the
`fourier_decay` runs. The (redundant) Setting 3 calls are commented out at the very
end. To re-run only the comparison and plotting step on already-saved results:

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
simulations_fofr_v2_with_ivanescus_...R       # Simulation loop (one rep per beta); also defines
                                               #   the 4 betas and the 2 X-generating processes
compare_methods_fofr_with_ivanescu_...R       # Aggregation, plotting, and Excel export
cv_penalties_fregre.basis.fr.R                # CV penalty grid wrapper for R&S method
cv_bases_fregre.basis.fr.R                    # (unused/dead code — not sourced anywhere)
predict_fregre_fr.R                           # Predict wrapper for R&S method
penFoFPLS/                                    # Local (git-ignored) copy of the estimation package
results_simulations/                          # All output (one subfolder per configuration)
docs/                                         # Design notes, reviews, and investigations (git-ignored)
```

Key files under `docs/` (all git-ignored, kept locally for reference):

| File | Contents |
|---|---|
| `new_betas.md` | Motivation and shape descriptions for betas 1–3 |
| `new_setting_smooth_X_beta4.md` | Motivation for the `fourier_decay` process and beta 4 |
| `new_fR2.md` | Definition and rationale for the functional $R^2$ metric |
| `simulation_review_findings.md` | Review of implementation vs. the paper's stated methodology, with fixes |
| `bspline_boundary_basis_issue.md` | Investigation of a B-spline boundary-knot construction issue and its (non-)relevance here |
| `penFoFPLS_package_update_compatibility.md` | Compatibility check against a reviewed `penFoFPLS` package build |

Output folder names encode the configuration:
`set{setting}{e}{_fourier_decay}_rep{total_reps}_pen{n_lambdas²}_K{KK}L{LL}/`
— e.g. `set2e_fourier_decay_rep30_pen100_K40L40/` is Setting 2, noisy $X$, the
smooth-$X$ predictor process, 30 reps, a $10\times10$ penalty grid, 40 bases.

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

Per-beta IMSE and CVE plots (one subfolder per beta: `cos_sin`, `sin_sum`,
`cos_sum`, and — for `fourier_decay` runs — `sin_sum_cos_diff`).

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

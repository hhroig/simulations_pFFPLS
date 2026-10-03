# "Buried PCA signal" scenario: fits all five methods (pFFR_RS, pFFR_I,
# pFFPLS, FFPLS, FFPLS_OB) on data where X is generated from its own
# Karhunen-Loeve/Fourier representation with a non-monotonic score-variance
# profile, and the true beta loads entirely on one deliberately
# low-variance ("buried") direction -- see
# docs/new_scenario_buried_pca_signal_2026-09-11.md for the full rationale
# and simulations_buried_pca_signal.R for the data-generating mechanism.
#
# Deliberately does NOT touch main_simulations_call.R,
# simulations_fofr_v2_with_ivanescus_ramsay_silverman.R,
# pFFPLS_only_lambda_refine_sim.R, or n_sensitivity_sim.R.
#
# Output files use the exact same naming/column conventions as the main
# study (method = "pFFR_RS" / "pFFR_I" / "pFFPLS" / "FFPLS" / "FFPLS_OB"),
# so they can be loaded and analyzed with the same code used throughout
# this project.

library(pls)
library(dplyr)
library(penFoFPLS)
library(fda)
library(refund)
library(reshape2)
library(doParallel)


main_buried_pca_call <- function(
    num_betas   = c("buried_pca"),
    total_reps  = 30,
    rep_starts  = 1,

    # Data-generating (X covariance / buried-signal) parameters:
    R_fourier   = 10,   # number of Fourier sin/cos pairs used to generate X
    dominant_r  = 1:4,  # indices carrying most of X's variance
    buried_r    = 5,    # index carrying the (low-variance) true signal
    dominant_sd = 1,    # sd of each dominant direction's score
    buried_sd   = 0.5,  # sd of the buried direction's score

    # Estimation-side setup. NOT K=L=40 (the rest of the study's Setting-2
    # convention) -- deliberately reduced. Two independent problems forced
    # this: (1) at K=40, pFFR_I (Ivanescu's pffr via mgcv) took 1.85-16.9
    # HOURS per rep (confirmed from a real 5-rep partial run) because X's
    # rank (9, by construction -- see
    # docs/new_scenario_buried_pca_signal_2026-09-11.md) is so much smaller
    # than K=40 that mgcv's rank-deficiency constraint-handling explodes; at
    # K=15 the same fit takes 26-37s (confirmed twice). (2) At K=40,
    # pFFR_RS's CV-selected penalty_X was pinned at its (already-once-raised)
    # floor in all 5 completed reps -- meaning the fit was being forced to
    # operate right at the edge of numerical breakdown, which likely
    # confounds any real "PLS beats covariance-blind regression" effect with
    # a pure conditioning artifact. K=15 (much closer to X's true rank 9)
    # gives pFFR_RS room to select an interior, non-boundary penalty instead.
    KK = 15, LL = 15,
    max_nComp = 8,

    num_lambdas    = 10,
    num_lambdas_RS = 8,

    # Diagnosed directly against the real penalty grid at K=L=15 (not
    # assumed): only 1 of 64 (penalty_X, penalty_Y) combinations failed
    # ("leading minor... not positive definite"), and it was penalty_X=1e-12
    # specifically -- penalty_X=1e-10 succeeded against every penalty_Y
    # tested, so -10 is used as the floor here (the K=L=40 attempt needed
    # -6, since more of the grid failed there).
    lower_penalty_bound_RS = -10, upper_penalty_bound_RS = -2,

    # pFFPLS's penalty_X and penalty_Y grids, DECOUPLED (previously shared
    # one 10-point grid spanning [-10,12], ~2.44 decades/step -- far too
    # coarse in the region that actually matters). Two independent sweeps
    # (see docs/buried_pca_signal_results_summary_2026-09-18.md) showed:
    #   - penalty_X is the real lever: IMSE flat at ~0.00999 (matching
    #     unpenalized FFPLS almost exactly) from 1e-16 up to ~1e-5.6, then
    #     rises sharply past ~1e-4 and is badly degraded by 1e-2 (the shared
    #     grid's old floor). The 33-rep run with the widened-but-still-coarse
    #     shared floor (-10) only landed a grid point *exactly* at -10 in
    #     33% of reps -- the rest spread across [-10,-2.67], much of which is
    #     already on the rising part of the curve, which is why the mean
    #     (0.047) stayed well above the median (0.010) and above FFPLS.
    #   - penalty_Y barely matters at all once penalty_X is good: with
    #     penalty_X fixed at 1e-8, IMSE was flat (~0.00996-0.01000) across
    #     the ENTIRE range from 1e-10 to 1e12.74, only failing (NaN) at the
    #     extreme 1e14. So resolution there is cheap to sacrifice.
    # Grid reallocated accordingly: penalty_X densely covering the
    # demonstrated plateau; penalty_Y kept wide but coarse (its exact value
    # barely affects the outcome).
    #
    # First attempt used upper_penalty_bound_X = -3 (the point where the
    # sweep showed IMSE was still only mildly degraded, ~0.014). A 5-rep
    # validation at that range showed CV landing at the grid's *upper* edge
    # (-3.47/-3.00) in 4 of 5 reps, not near the true optimum at -10 -- CV
    # selects penalty_X by minimizing *prediction* error on held-out Y, not
    # *beta-recovery* error, and those are different curves: with the buried
    # signal such a small share of Y's variance, prediction error is much
    # flatter/noisier in penalty_X than the beta-recovery curve is, so CV
    # kept wandering into a region that predicts fine but recovers beta
    # poorly. Fix: narrow the upper bound to -5 -- still comfortably inside
    # the flat part of the beta-recovery curve (IMSE ~0.0099 at -5.58) but
    # excluding the region CV kept mistakenly preferring.
    lower_penalty_bound_X = -10, upper_penalty_bound_X = -5, num_lambdas_X = 16,
    lower_penalty_bound_Y = -8,  upper_penalty_bound_Y  = 12, num_lambdas_Y = 6,

    nnodesX = 100,
    nnodesY = 100,
    num_folds = 5,
    center = TRUE

) {

  KK_rs <- LL_rs <- KK
  do_opt_bases_FFPLS <- TRUE  # still run FFPLS_OB's basis-count search even at this smaller K

  argvals_X <- seq(0, 1, length.out = nnodesX) # p
  argvals_Y <- seq(0, 1, length.out = nnodesY) # q

  basisobj_X <- fda::create.bspline.basis(rangeval = range(argvals_X), nbasis = KK)
  basisobj_Y <- fda::create.bspline.basis(rangeval = range(argvals_Y), nbasis = LL)

  LL_list <- round(seq(9, LL, length.out = num_lambdas)) # basis-count search for FFPLS_OB
  KK_list <- LL_list

  # Penalty grids, same construction as main_simulations_call.R:
  lambdas_in_RS   <- 10^seq(lower_penalty_bound_RS, upper_penalty_bound_RS, length.out = num_lambdas_RS)
  penaltyvec_X_RS <- penaltyvec_Y_RS <- lambdas_in_RS

  # pFFPLS's decoupled penalty grids (see rationale above):
  penaltyvec_X <- 10^seq(lower_penalty_bound_X, upper_penalty_bound_X, length.out = num_lambdas_X)
  penaltyvec_Y <- 10^seq(lower_penalty_bound_Y, upper_penalty_bound_Y, length.out = num_lambdas_Y)

  # Own dedicated results folder, nested under results_simulations/ (already
  # entirely gitignored). "_finegrid" tags this as the decoupled/fine-X-
  # resolution pFFPLS grid, distinct from the earlier shared-coarse-grid run
  # at rep33_K15L15/ (kept on disk as a documented before/after comparison,
  # not overwritten):
  out_folder <- paste0("results_simulations/buried_pca_signal/",
                       "rep", total_reps,
                       "_K", KK, "L", LL,
                       "_finegrid",
                       "/")

  if (!dir.exists(out_folder)) {
    dir.create(out_folder, recursive = TRUE)
  }

  # Ramsay & Silverman wrappers (unchanged, existing files):
  source("cv_penalties_fregre.basis.fr.R")
  source("predict_fregre_fr.R")

  nodes_CL <- detectCores()
  cl <- makeCluster(nodes_CL)
  clusterExport(cl, c("predict_fregre_fr"))
  registerDoParallel(cl)

  source("simulations_buried_pca_signal.R", local = TRUE)

  stopCluster(cl)

}


# Run: the real study, all 5 methods, K=L=15, with pFFPLS's fine/decoupled
# penalty grid (see rationale above) -- writes to rep{total_reps}_K15L15_finegrid/.
main_buried_pca_call(
  num_betas   = c("buried_pca"),
  total_reps  = 36,
  rep_starts  = 1,
  KK = 15, LL = 15
)

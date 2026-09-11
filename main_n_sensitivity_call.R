# Standalone n-sensitivity study: fits ONLY pFFPLS and pFFR_RS (the two
# methods actually being compared) on the three reviewer-response betas
# (local_ridge, asym_step, aniso_peak), at a REDUCED training sample size,
# to test whether that narrows or flips pFFPLS's IMSE(beta) gap against
# pFFR_RS -- see n_sensitivity_sim.R for the full rationale.
#
# Deliberately does NOT touch main_simulations_call.R,
# simulations_fofr_v2_with_ivanescus_ramsay_silverman.R, or
# pFFPLS_only_lambda_refine_sim.R. Skips pFFR_I/FFPLS/FFPLS_OB entirely
# (same reasoning as the lambda-refine study: they aren't part of this
# comparison and dropping them keeps a rep to a few minutes rather than
# ~45-70). Uses the STANDARD (already-fixed) default penalty grids for both
# methods, not refined ones -- the earlier lambda-refinement study already
# ruled out hyperparameter search as an explanation, so this test is only
# about n.
#
# Output files use the exact same naming/column conventions as the main
# study (method = "pFFPLS" / "pFFR_RS"), so they can be loaded and compared
# directly against the existing n=100 results in
# results_simulations/set2_genK50_rep20_pen100_K40L40/ with the same
# analysis code -- no merging or reformatting needed. Reps here are fresh,
# independent draws (not seed-matched), since the comparison is
# distributional, not rep-by-rep paired.

library(pls)
library(dplyr)
library(penFoFPLS)
library(fda)
library(refund)
library(reshape2)
library(doParallel)


main_n_sensitivity_call <- function(
    num_betas   = c("local_ridge", "asym_step", "aniso_peak"),
    total_reps  = 10,
    rep_starts  = 1,

    n_train      = 40,  # reduced training n -- the lever being tested. The
                        # n=100 baseline run used n=100 against K=L=40 (1600
                        # coefficients), plenty of data for pFFR_RS's
                        # full-rank fit to do well via smoothness alone.
                        # n_train=40 (~ K) is a deliberately aggressive
                        # reduction, putting pFFR_RS in a genuinely
                        # high-dimensional-relative-to-n regime where
                        # pFFPLS's rank constraint should matter most, if it
                        # ever does.
    n_validation = 50,  # kept the same as the baseline so test-set IMSE(Y)
                        # stays comparably estimated

    gen_nbasis  = 50,   # matches the reviewer-response betas' generating basis

    num_lambdas    = 10,  # pFFPLS grid -- standard default (unrefined; the
                          # lambda-refinement study already showed refining
                          # this doesn't change the outcome)
    num_lambdas_RS = 8,   # pFFR_RS grid -- standard default

    lower_penalty_bound_RS = -12, upper_penalty_bound_RS = -2,
    lower_penalty_bound    = -2,  upper_penalty_bound    = 12,

    nnodesX = 100,
    nnodesY = 100,
    num_folds = 5,
    center = TRUE

) {

  # Setting 2 basis dimension, matching the existing 5-method run these
  # betas were tested under -- only n_train changes here, not K/L:
  KK <- LL <- 40
  KK_rs <- LL_rs <- KK
  max_nComp <- 8

  argvals_X <- seq(0, 1, length.out = nnodesX) # p
  argvals_Y <- seq(0, 1, length.out = nnodesY) # q

  basisobj_X <- fda::create.bspline.basis(rangeval = range(argvals_X), nbasis = KK)
  basisobj_Y <- fda::create.bspline.basis(rangeval = range(argvals_Y), nbasis = LL)

  # Penalty grids, same construction as main_simulations_call.R:
  lambdas_in_RS   <- 10^seq(lower_penalty_bound_RS, upper_penalty_bound_RS, length.out = num_lambdas_RS)
  penaltyvec_X_RS <- penaltyvec_Y_RS <- lambdas_in_RS

  lambdas_in   <- 10^seq(lower_penalty_bound, upper_penalty_bound, length.out = num_lambdas)
  penaltyvec_X <- penaltyvec_Y <- lambdas_in

  # Own dedicated results folder, nested under results_simulations/ (already
  # entirely gitignored), self-documenting on n_train so it can't be
  # confused with the n=100 baseline:
  out_folder <- paste0("results_simulations/n_sensitivity/",
                       "n", n_train,
                       "_rep", total_reps,
                       "_K", KK, "L", LL,
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

  source("n_sensitivity_sim.R", local = TRUE)

  stopCluster(cl)

}


# Run: quick test, n_train=40 (vs. the n=100 baseline), same 3
# reviewer-response betas and Setting-2 basis dims (K=L=40, gen_nbasis=50),
# 10 reps -- comparable to the n=100 numbers already in
# results_simulations/set2_genK50_rep20_pen100_K40L40/.
main_n_sensitivity_call(
  num_betas   = c("local_ridge", "asym_step", "aniso_peak"),
  total_reps  = 10,
  rep_starts  = 1,
  n_train     = 40,
  gen_nbasis  = 50
)

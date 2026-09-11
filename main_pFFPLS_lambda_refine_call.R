# Standalone pFFPLS-only study to refine the penalty (lambda) grid on the
# three reviewer-response betas (local_ridge, asym_step, aniso_peak), using
# per-beta ranges informed by where CV actually selected penalties in the
# existing 5-method, 20-rep run -- see pFFPLS_only_lambda_refine_sim.R for
# the full rationale and the ranges themselves.
#
# Deliberately does NOT touch main_simulations_call.R or
# simulations_fofr_v2_with_ivanescus_ramsay_silverman.R. It skips
# pFFR_RS/pFFR_I/FFPLS/FFPLS_OB entirely (those blocks in the main script
# never depend on pFFPLS's own block, so this is a clean removal), which is
# what makes a finer, better-targeted lambda grid affordable: a rep here
# costs a few minutes instead of the ~45-70 min a full 5-method iteration
# takes (dominated by pFFR_I's single-threaded pffr() fit).
#
# Output files use the exact same naming/column conventions as the main
# study (method = "pFFPLS" throughout), so they can be loaded and compared
# directly against the pFFR_RS/etc. numbers already in
# results_simulations/set2_genK50_rep20_pen100_K40L40/ with the same
# analysis code -- no merging or reformatting needed. Reps here are fresh,
# independent draws (not seed-matched to that run) since the comparison is
# distributional (does pFFPLS's IMSE(beta) distribution improve/shift),
# not a rep-by-rep pairing.

library(pls)
library(dplyr)
library(penFoFPLS)
library(fda)
library(refund)
library(reshape2)
library(doParallel)


main_pFFPLS_lambda_refine_call <- function(
    num_betas   = c("local_ridge", "asym_step", "aniso_peak"),
    total_reps  = 20,
    rep_starts  = 1,

    gen_nbasis  = 50,   # matches the reviewer-response betas' generating basis

    num_lambdas = 12,   # per-axis grid resolution -- higher than the shared
                        # study's 10, affordable now that pFFR_I/pFFR_RS/
                        # FFPLS/FFPLS_OB aren't being refit every rep

    nnodesX = 100,
    nnodesY = 100,
    num_folds = 5,
    center = TRUE

) {

  # Setting 2 basis dimension, matching the existing 5-method run these
  # betas were tested under:
  KK <- LL <- 40
  max_nComp <- 8

  argvals_X <- seq(0, 1, length.out = nnodesX) # p
  argvals_Y <- seq(0, 1, length.out = nnodesY) # q

  basisobj_X <- fda::create.bspline.basis(rangeval = range(argvals_X), nbasis = KK)
  basisobj_Y <- fda::create.bspline.basis(rangeval = range(argvals_Y), nbasis = LL)

  # Own dedicated results folder, nested under results_simulations/ (already
  # entirely gitignored) but kept separate from the main study's flat
  # set2_*/set1_* folders:
  out_folder <- paste0("results_simulations/pFFPLS_lambda_refine/",
                       "rep", total_reps,
                       "_K", KK, "L", LL,
                       "/")

  if (!dir.exists(out_folder)) {
    dir.create(out_folder, recursive = TRUE)
  }

  nodes_CL <- detectCores()
  cl <- makeCluster(nodes_CL)
  registerDoParallel(cl)

  source("pFFPLS_only_lambda_refine_sim.R", local = TRUE)

  stopCluster(cl)

}


# Run: pFFPLS-only, refined per-beta lambda grid, same 3 reviewer-response
# betas and Setting-2 basis dims (K=L=40, gen_nbasis=50) as the existing
# 20-rep study, so the resulting IMSE(beta) distribution is directly
# comparable to the pFFR_RS/pFFR_I/FFPLS/FFPLS_OB numbers already in
# results_simulations/set2_genK50_rep20_pen100_K40L40/.
main_pFFPLS_lambda_refine_call(
  num_betas   = c("local_ridge", "asym_step", "aniso_peak"),
  total_reps  = 20,
  rep_starts  = 1,
  gen_nbasis  = 50,
  num_lambdas = 12
)

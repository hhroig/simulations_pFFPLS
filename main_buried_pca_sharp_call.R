# "Buried PCA signal, SHARP shape" scenario: fits all five methods (pFFR_RS,
# pFFR_I, pFFPLS, FFPLS, FFPLS_OB) on data where X is generated from its own
# Fourier representation with a non-monotonic score-variance profile, and the
# true beta is a genuinely SHARP, LOCALIZED pulse (Reviewer 1's ask #1)
# loading entirely on a CLUSTER of low-variance ("buried") directions -- see
# docs/new_scenario_buried_pca_sharp_signal_2026-09-28.md for the full
# rationale and simulations_buried_pca_sharp_signal.R for the data-generating
# mechanism.
#
# Deliberately does NOT touch main_simulations_call.R,
# simulations_fofr_v2_with_ivanescus_ramsay_silverman.R,
# pFFPLS_only_lambda_refine_sim.R, n_sensitivity_sim.R,
# simulations_buried_pca_signal.R, or main_buried_pca_call.R.
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


main_buried_pca_sharp_call <- function(
    num_betas   = c("buried_pca_sharp"),
    total_reps  = 30,
    rep_starts  = 1,

    # Data-generating (X covariance / buried-bump-signal) parameters:
    R_fourier   = 10,    # number of Fourier sin/cos pairs used to generate X
    dominant_r  = 1:4,   # indices carrying most of X's variance (sin+cos)
    bump_r      = 5:10,  # CLUSTER of frequencies whose COS term is buried
                         # (low-variance) and carries the true, sharp signal
    dominant_sd = 1,     # sd of each dominant direction's score
    buried_sd   = 0.3,   # sd of each bump direction's score (~6.3% of total
                         # X variance jointly: 6*0.3^2 / (8*1^2 + 6*0.3^2) --
                         # verified numerically, see scenario doc. NOTE: an
                         # earlier design draft used buried_sd=0.45 (~13%
                         # share) before being reduced to 0.3; make sure any
                         # reported variance-share figure matches THIS value)
    Rsq         = 0.9,   # target R2(q); lowered in the SNR-mismatch variant
                         # (see docs/buried_pca_sharp_signal_results_summary_2026-09-28.md)
                         # to test whether more Y-noise, at fixed K=25 (already
                         # ~11 dims above X's true rank 14), creates the
                         # estimation-variance gap over unpenalized FFPLS that
                         # the sharp shape alone did not

    # Estimation-side setup. X's true rank here is
    # 2*length(dominant_r) + length(bump_r) = 8 + 6 = 14 (vs. rank 9 in the
    # original single-sinusoid buried_pca scenario). KK=LL=25 chosen and
    # explicitly safety-checked before use (see scenario doc): resolves the
    # bump's highest frequency (10, period 0.1) at ~2.2 knots/period
    # (comparable margin to the original scenario's K=15/freq=5 case at
    # ~2.4 knots/period), and a direct 36-combination penalty grid test
    # against real generated data found NO Cholesky failures for pFFR_RS
    # anywhere in the standard [1e-12,1e-2] range (failures only appeared at
    # extreme, unused penalties >= 1e4) -- unlike the original scenario's
    # K=40 attempt, this K does NOT need a raised penalty floor. pFFR_I
    # (Ivanescu/mgcv) measured at ~236s/fit at K=25 (vs. ~111s/fit at K=20)
    # -- both tractable; K=25 was kept for the better frequency resolution.
    KK = 25, LL = 25,
    max_nComp = 8,

    num_lambdas    = 10,

    # pFFR_RS's penalty_X/penalty_Y grids, DECOUPLED (previously shared,
    # [1e-12,1e-2] on both axes). A direct sweep of achievable IMSE(beta)
    # (not just CV-selected IMSE) found the real 30-rep run's shared grid
    # never tested penalty_Y above 1e-2, but the true optimum needs
    # penalty_Y ~1e4 paired with a low penalty_X: beta_true is smooth in q
    # (a single sin(2*pi*q) cycle, wants heavy smoothing) but sharp/
    # multi-frequency in p (wants minimal smoothing). At
    # (penalty_X, penalty_Y) = (1e-10, 1e4), mean IMSE(beta) = 0.199
    # (sd 0.009 across 10 independent draws) -- ~390x better than the
    # shared-grid run's reported pFFR_RS mean of 76.85. Still ~31x worse
    # than pFFPLS (0.0064), so the qualitative conclusion is unchanged, but
    # the OLD number was a real grid artifact, not a clean reflection of
    # the method (docs/simulation_study_buried_signal_writeup_2026-10-03.md
    # has the full sweep; see also the dedicated diagnostic scripts run to
    # find and verify this).
    #
    # penalty_Y is NOT a smooth log-sequence here: a direct safety check at
    # the realistic 80-observation CV-fold size (fda.usc::fregre.basis.fr
    # has NO error handling inside cv_penalties_fregre.basis.fr's foreach
    # loop -- a single failing combo would crash the whole CV call) found
    # penalty_Y in {1e2, 1e3, 1e5, 1e6} are genuinely failure-prone
    # (numerical pathology specific to this K=25/signal combination, not a
    # smooth function of the penalty -- e.g. 1e3 failed 8/8 draws at
    # penalty_X=1e-10 while 1e4 succeeded 8/8), while {1e-2,1e-1,1e0,1e1,1e4}
    # succeeded in all 200/200 fits tested (8 draws x 5 penalty_X values).
    # The grid below is built from that evidence directly, not a formula.
    penaltyvec_X_RS = 10^seq(-10, -2, length.out = 5),
    penaltyvec_Y_RS = c(1e-2, 1e-1, 1e0, 1e1, 1e4),

    # pFFPLS's penalty_X/penalty_Y grids. Deliberately NOT reusing the
    # tightly-tuned [-10,-5] floor from the original (single-sinusoid)
    # buried_pca scenario -- that range was fit to a different signal
    # (single frequency vs. a 6-frequency sharp bump) and assuming it
    # transfers unchanged would repeat the mistake already made once in
    # this project (see docs/buried_pca_signal_results_summary_2026-09-18.md,
    # sections 2.4-2.8: each new signal structure needed its own evidence-
    # based grid check). Starting from the main study's original wide
    # default instead; narrow later from a real CV-selection sweep once this
    # scenario's own pilot/validation results are in.
    lower_penalty_bound_X = -10, upper_penalty_bound_X = 12, num_lambdas_X = 10,
    lower_penalty_bound_Y = -10, upper_penalty_bound_Y = 12, num_lambdas_Y = 10,

    nnodesX = 100,
    nnodesY = 100,
    num_folds = 5,
    center = TRUE

) {

  KK_rs <- LL_rs <- KK
  do_opt_bases_FFPLS <- TRUE

  argvals_X <- seq(0, 1, length.out = nnodesX) # p
  argvals_Y <- seq(0, 1, length.out = nnodesY) # q

  basisobj_X <- fda::create.bspline.basis(rangeval = range(argvals_X), nbasis = KK)
  basisobj_Y <- fda::create.bspline.basis(rangeval = range(argvals_Y), nbasis = LL)

  LL_list <- round(seq(9, LL, length.out = num_lambdas)) # basis-count search for FFPLS_OB
  KK_list <- LL_list

  # penaltyvec_X_RS / penaltyvec_Y_RS are now passed in directly as
  # arguments (decoupled, evidence-based grid -- see rationale above),
  # not constructed from a shared log-sequence.

  penaltyvec_X <- 10^seq(lower_penalty_bound_X, upper_penalty_bound_X, length.out = num_lambdas_X)
  penaltyvec_Y <- 10^seq(lower_penalty_bound_Y, upper_penalty_bound_Y, length.out = num_lambdas_Y)

  # Own dedicated results folder. Rsq tagged into the path (only when it
  # departs from the 0.9 used everywhere else in this project) so the SNR-
  # mismatch variant can never silently collide with/overwrite the
  # Rsq=0.9 runs' output. "_RSfix" tags runs using the corrected, decoupled
  # pFFR_RS grid (see rationale above) -- kept distinct from the earlier
  # rep30_K25L25_Rsq0.5/ run (shared/narrow RS grid, pFFR_RS mean IMSE
  # inflated to 76.85 by the grid artifact), which stays on disk unchanged
  # as documented history, not overwritten:
  rsq_tag <- if (Rsq != 0.9) paste0("_Rsq", Rsq) else ""
  out_folder <- paste0("results_simulations/buried_pca_sharp_signal/",
                       "rep", total_reps,
                       "_K", KK, "L", LL,
                       rsq_tag,
                       "_RSfix",
                       "/")

  if (!dir.exists(out_folder)) {
    dir.create(out_folder, recursive = TRUE)
  }

  source("cv_penalties_fregre.basis.fr.R")
  source("predict_fregre_fr.R")

  nodes_CL <- detectCores()
  cl <- makeCluster(nodes_CL)
  clusterExport(cl, c("predict_fregre_fr"))
  registerDoParallel(cl)

  source("simulations_buried_pca_sharp_signal.R", local = TRUE)

  stopCluster(cl)

}


# Smoke test (1 rep, small grids, Rsq=0.9 default) and a 5-rep validation
# (Rsq=0.9, full settings) already run and confirmed clean -- see
# docs/buried_pca_sharp_signal_results_summary_2026-09-28.md. Result:
# pFFPLS still only ties FFPLS (mean IMSE 0.004978 vs 0.004974, wins 3/5
# reps) even with the sharp/localized shape -- confirms the shape lever
# alone isn't enough; proceeding to the SNR-mismatch variant below.
#
# SNR-mismatch variant: same design (K=25, sharp bump, bump_r=5:10),
# Rsq lowered from 0.9 to 0.5 (noise variance raised from ~11% to 100% of
# signal variance) to test whether more Y-noise creates the
# estimation-variance gap over unpenalized FFPLS that the sharp shape alone
# did not. 1-rep smoke test first, to confirm the noisier regime doesn't
# break anything (e.g. pFFR_RS numerics), before a real validation.
# 15-rep robustness check (rep15_K25L25_Rsq0.5/) confirmed the Rsq=0.5
# result: pFFPLS has the best mean IMSE(beta) of all five methods --
# significant win over FFPLS (14/15 reps, paired t-test p=0.025) and a
# statistical tie (slight edge) with pFFR_I (9/15, p=0.37, reversing the
# 5-rep run's apparent pFFR_I lead, confirming that was small-sample
# noise). See docs/buried_pca_sharp_signal_results_summary_2026-09-28.md,
# section 4. This is the strongest result in the project so far.
#
# 30-rep real run completed at rep30_K25L25_Rsq0.5/ (Claude ran it directly
# after an RStudio launch issue). Confirmed headline result: pFFPLS beats
# FFPLS (p=0.0086) and pFFR_I (p=0.00023) at n=30.
#
# INVESTIGATION (closed): is pFFR_RS's mean IMSE (76.85) a penalty-grid
# artifact? A direct achievable-minimum search (bypassing CV) found a real
# effect -- (penalty_X,penalty_Y)=(1e-10,1e4) achieves mean IMSE 0.199,
# ~390x better. But re-testing that SAME widened, decoupled grid (now the
# default penaltyvec_X_RS/penaltyvec_Y_RS above) through the actual 5-fold
# CV pipeline (5-rep validation, rep5_K25L25_Rsq0.5_RSfix/) found CV
# selects penalty_Y=0.01 in 5/5 reps -- never 1e4 -- giving mean IMSE 70.17,
# not meaningfully different from the original 76.85. CONCLUSION: this is
# NOT a range artifact after all; it's a genuine CV prediction-vs-recovery
# mismatch (same mechanism as pFFPLS's, documented earlier, but more severe
# here) -- full writeup in
# docs/buried_pca_sharp_signal_results_summary_2026-09-28.md section 6 and
# docs/simulation_study_buried_signal_writeup_2026-10-03.md section 6.1.
# NO RE-RUN of the 30-rep study was done (the 5-rep check's unanimous,
# zero-variation CV selection is strong enough evidence on its own) --
# rep30_K25L25_Rsq0.5/ (mean pFFR_RS IMSE = 76.853) remains the correct,
# reported result. The decoupled grid above is kept as the more thorough
# setup going forward even though it didn't change the headline number.
#
# The call below is the already-completed 5-rep validation that produced
# the 70.17 finding above -- safe to leave as the active call (cheap,
# ~15-20 min, re-running it changes nothing material). Re-run the 30-rep
# call further up in this file's history only if a fresh independent draw
# of the full study is actually needed.
main_buried_pca_sharp_call(
  num_betas   = c("buried_pca_sharp"),
  total_reps  = 5,
  rep_starts  = 1,
  KK = 25, LL = 25,
  Rsq = 0.5
)

## Standalone pFFPLS-only simulation loop, for refining the penalty (lambda)
## grid on the three reviewer-response betas (local_ridge, asym_step,
## aniso_peak) -- see docs/new_betas_reviewer1_response_2026-09-06.md.
##
## This does NOT touch simulations_fofr_v2_with_ivanescus_ramsay_silverman.R
## or main_simulations_call.R. It fits ONLY pFFPLS (skips pFFR_RS, pFFR_I,
## FFPLS, FFPLS_OB entirely -- their blocks in the main script are
## self-contained and never touched by pFFPLS's own block, so this is a
## clean removal, not a rewrite). Dropping those four cuts a ~45-70 min
## iteration (dominated by pFFR_I's single-threaded ~20-40 min pffr() fit)
## down to a few minutes, which is what makes a finer/better-targeted lambda
## grid affordable here.
##
## Per-beta penalty ranges below are informed by where CV actually selected
## penalties in the existing 20-rep, 5-method run
## (results_simulations/set2_genK50_rep20_pen100_K40L40/):
##   - local_ridge: penalty_X was pinned at the old floor (1e-2) in 98% of
##     fits -- makes sense, its defining feature is a sharp ridge in the
##     p-direction, so CV wants penalty_X as low as possible. Extended much
##     lower here (1e-8) to actually bracket the optimum instead of
##     censoring it.
##   - asym_step: penalty_Y was floor-censored ~31% of the time -- extended
##     lower (1e-6) for the same reason.
##   - aniso_peak: penalty_Y genuinely used BOTH the old floor and ceiling
##     (~23% each) -- kept wide. The ceiling here (1e12) is a hard numerical
##     limit, not a design choice: pushing pFFPLS's upper bound to 1e14
##     during the main study made cv_unique_fof_par's internal solve()
##     singular at K=40 (reciprocal condition number ~1.6e-17) and crashed
##     the run -- see docs/new_betas_reviewer1_response_2026-09-06.md,
##     "Penalty (lambda) grid design".
##
## Output files use the exact same naming/column conventions as the main
## study (final_models_rep_*, cves_rep_*, best_lambdas_rep_*, betas_rep_*,
## computation_times_rep_*, fR2s_rep_*, all with method = "pFFPLS"), so they
## can be loaded and compared directly against the existing pFFR_RS/etc.
## results with the same analysis code -- no merging or reformatting needed.

# Reuse generate_fofr_data(), redo_beta_true(), fR2_function(), and
# mean_imse_by_row() from the main simulation script WITHOUT sourcing (and
# thereby executing) its own 5-method simulation loop: parse the file and
# evaluate only the top-level expressions before the first `for` loop, which
# is exactly where those definitions live. Avoids duplicating them here (and
# risking drift) while never touching that file.
.sim_script_exprs <- parse("simulations_fofr_v2_with_ivanescus_ramsay_silverman.R")
for (.e in .sim_script_exprs) {
  if (is.call(.e) && identical(.e[[1]], as.name("for"))) break
  eval(.e, envir = .GlobalEnv)
}
rm(.sim_script_exprs, .e)
stopifnot(exists("generate_fofr_data"), exists("redo_beta_true"),
          exists("fR2_function"), exists("mean_imse_by_row"))

# Per-beta penalty ranges (log10 bounds) -- see header comment above.
lambda_ranges <- list(
  local_ridge = list(lower_X = -8, upper_X = 2, lower_Y = -2, upper_Y = 5),
  asym_step   = list(lower_X = -2, upper_X = 2, lower_Y = -6, upper_Y = 3),
  aniso_peak  = list(lower_X = -3, upper_X = 5, lower_Y = -4, upper_Y = 12)
)


for (beta.num in num_betas) {

  if (!(beta.num %in% names(lambda_ranges))) {
    stop("No refined lambda range configured for beta '", beta.num,
         "' -- add one to `lambda_ranges` in pFFPLS_only_lambda_refine_sim.R.")
  }
  rng <- lambda_ranges[[beta.num]]
  penaltyvec_X <- 10^seq(rng$lower_X, rng$upper_X, length.out = num_lambdas)
  penaltyvec_Y <- 10^seq(rng$lower_Y, rng$upper_Y, length.out = num_lambdas)

  for (rep_num in rep_starts:total_reps) {

    # Initialize savings:
    all_CVEs <- data.frame()
    all_final_res <- data.frame()
    all_beta_hats <- data.frame()
    all_fr2s <- data.frame()
    best_lambdas <- data.frame()
    computation_times <- data.frame()

    cat("Doing beta ", beta.num, " (pFFPLS-only, refined grid)\n")
    cat("     rep ", rep_num, "/", total_reps, "\n")


    # Generate data -----------------------------------------------------------
    cat("    -> generating data\n")

    gendata <- generate_fofr_data(nbasisX = gen_nbasis, nbasisY = gen_nbasis, nbeta = beta.num,
                                  nnodesX = nnodesX, nnodesY = nnodesY,
                                  n = 100, n_validation = 50,
                                  X_process = "uniform_bspline")

    X_gen <- gendata$X
    Y <- gendata$Y
    beta_true <- gendata$beta_true

    X_gen_val <- gendata$X_val
    Y_val <- gendata$Y_val

    # No predictor observation noise here -- matches the X_sd_error = 0
    # setting the 20-rep comparison run used:
    X <- X_gen
    X_val <- X_gen_val

    folds <- caret::createFolds(1:nrow(Y), k = num_folds)


    # CV - pFFPLS (refined per-beta grid) ------------------------------------

    cat("    -> starting pFFPLS (refined grid:",
        length(penaltyvec_X), "x", length(penaltyvec_Y), "penalty combos)\n")

    time_penalized <- system.time({  # Measure time
      cv_penalized <- cv_unique_fof_par(X = X,
                                        Y = Y,
                                        argvals_X = argvals_X,
                                        argvals_Y = argvals_Y,
                                        ncomp = max_nComp,
                                        center = center,
                                        folds = folds,
                                        basisobj_X = basisobj_X,
                                        basisobj_Y = basisobj_Y,
                                        penaltyvec_X = penaltyvec_X,
                                        penaltyvec_Y = penaltyvec_Y,
                                        verbose = TRUE,
                                        stripped = TRUE,
                                        maxit = 100000 )
    })

    cat("    -> pFFPLS done\n")

    best_lambdas <- rbind(
      best_lambdas,
      data.frame(lambda.penalty_X = cv_penalized$best_penalties[, "penalty_X"],
                 lambda.penalty_Y = cv_penalized$best_penalties[, "penalty_Y"],
                 method = "pFFPLS",
                 nComp = 1:max_nComp,
                 beta.num = beta.num,
                 rep_num = rep_num)
    )

    all_CVEs <- rbind(
      all_CVEs,
      data.frame(
        CVE = as.numeric(  cv_penalized$CVEs_ncomp  ),
        method = "pFFPLS",
        nComp = 1:max_nComp,
        beta.num = beta.num,
        rep_num = rep_num
      )
    )

    computation_times <- rbind(
      computation_times,
      data.frame(
        method = "pFFPLS",
        beta.num = beta.num,
        nComp = 1:max_nComp,
        rep_num = rep_num,
        elapsed_time = time_penalized[["elapsed"]]
      )
    )


    # Final models for each ncomp ---------------------------------------------

    for (nComp in 1:max_nComp) {
      cat("       ---> final model cmpt:", nComp, "/", max_nComp, "\n")

      m_final <- ffpls_bs(X = X,
                          Y = Y,
                          argvals_X = argvals_X,
                          argvals_Y = argvals_Y,
                          ncomp = nComp,
                          center = center,
                          basisobj_X = basisobj_X,
                          basisobj_Y = basisobj_Y,
                          penalty_X = cv_penalized$best_penalties[nComp, "penalty_X"],
                          penalty_Y = cv_penalized$best_penalties[nComp, "penalty_Y"],
                          verbose = FALSE,
                          stripped = F,
                          maxit = 100000)

      beta_hat <- m_final$coefficient_function[, , nComp]

      beta_df <- expand.grid(q = argvals_Y, p = argvals_X)
      beta_df$z <- as.vector(beta_hat)
      beta_df$z_true <- as.vector(beta_true)
      beta_df$method <- "pFFPLS"
      beta_df$nComp <- nComp
      beta_df$beta.num <- beta.num
      beta_df$rep_num <- rep_num

      all_beta_hats <- rbind(all_beta_hats, beta_df )

      # Compute the average IMSE for the validation set:
      mean_imse_Y_val <- mean_imse_by_row(Y_val, # true
                                          predict(object = m_final, newdata = X_val)[, , nComp], # predicted
                                          argvals_Y)

      all_final_res <- rbind(all_final_res,
                             tibble::tibble(
                               imse = penFoFPLS::imse_beta_ffpls_bs(beta_true,
                                                                    beta_hat,
                                                                    argvals_X,
                                                                    argvals_Y),
                               mean_imse_Y_val = mean_imse_Y_val,
                               nComp = nComp,
                               beta.num = beta.num,
                               rep_num = rep_num,
                               method = "pFFPLS")
      )

      # Compute fR2 for the training and validation sets:
      y_pred <- fitted.values(m_final)[ , , nComp]
      y_val_pred <- predict(object = m_final, newdata = X_val)[, , nComp]

      all_fr2s <- rbind(all_fr2s,
                        tibble::tibble(
                          fR2_train = fR2_function(Y, y_pred),
                          fR2_val   = fR2_function(Y_val, y_val_pred),
                          nComp    = nComp,
                          beta.num = beta.num,
                          rep_num  = rep_num,
                          method   = "pFFPLS")
      )

      rm(m_final, y_pred, y_val_pred, beta_df, beta_hat, mean_imse_Y_val)

    }  # end nComp loop


    # Savings:

    saveRDS(all_CVEs, file =
              paste0(out_folder, "cves_rep_", rep_num, "_beta_", beta.num, ".Rds"))

    saveRDS(all_final_res, file =
              paste0(out_folder, "final_models_rep_", rep_num, "_beta_", beta.num, ".Rds"))

    saveRDS(all_beta_hats, file =
              paste0(out_folder, "betas_rep_", rep_num, "_beta_", beta.num, ".Rds"))

    saveRDS(best_lambdas, file =
              paste0(out_folder, "best_lambdas_rep_", rep_num, "_beta_", beta.num, ".Rds"))

    saveRDS(all_fr2s, file =
              paste0(out_folder, "fR2s_rep_", rep_num, "_beta_", beta.num, ".Rds"))

    saveRDS(computation_times, file =
              paste0(out_folder, "computation_times_rep_", rep_num, "_beta_", beta.num, ".Rds"))


  } # end nested inner rep_num

} # end nested outer beta_num

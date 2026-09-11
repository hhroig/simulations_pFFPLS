## Standalone n-sensitivity study: does shrinking the training sample size
## narrow (or flip) pFFPLS's IMSE(beta) gap against pFFR_RS on the three
## reviewer-response betas (local_ridge, asym_step, aniso_peak)?
##
## Rationale (see the "what to try next" discussion): every prior test used
## n=100 against K=L=40 (1600 coefficients) -- plenty of data for a
## full-rank, globally-penalized fit (pFFR_RS) to estimate well via
## smoothness alone, leaving little room for pFFPLS's low-rank constraint to
## pay off. Shrinking n should worsen pFFR_RS's conditioning while leaving
## pFFPLS comparatively protected by its rank constraint -- the classical
## regime where PLS/PCR-type dimension reduction is expected to beat
## full-rank fits. This is the one lever not yet tested; the earlier lambda-
## grid refinement (see pFFPLS_only_lambda_refine_sim.R) ruled out
## hyperparameter search as the explanation, so this only fits pFFPLS and
## pFFR_RS -- the two methods actually being compared -- with the STANDARD
## (already-fixed) default penalty grids, not refined ones.
##
## Does NOT touch simulations_fofr_v2_with_ivanescus_ramsay_silverman.R,
## main_simulations_call.R, or pFFPLS_only_lambda_refine_sim.R. Output files
## use the same naming/column conventions as the main study (method =
## "pFFPLS" / "pFFR_RS"), so they can be loaded and compared directly against
## the existing n=100 results in
## results_simulations/set2_genK50_rep20_pen100_K40L40/ with the same
## analysis code.

# Reuse generate_fofr_data(), redo_beta_true(), fR2_function(), and
# mean_imse_by_row() from the main simulation script WITHOUT sourcing (and
# thereby executing) its own 5-method simulation loop -- same technique used
# in pFFPLS_only_lambda_refine_sim.R: parse the file, evaluate only the
# top-level expressions before the first `for` loop.
.sim_script_exprs <- parse("simulations_fofr_v2_with_ivanescus_ramsay_silverman.R")
for (.e in .sim_script_exprs) {
  if (is.call(.e) && identical(.e[[1]], as.name("for"))) break
  eval(.e, envir = .GlobalEnv)
}
rm(.sim_script_exprs, .e)
stopifnot(exists("generate_fofr_data"), exists("redo_beta_true"),
          exists("fR2_function"), exists("mean_imse_by_row"))


for (beta.num in num_betas) {

  for (rep_num in rep_starts:total_reps) {

    # Initialize savings:
    all_CVEs <- data.frame()
    all_final_res <- data.frame()
    all_beta_hats <- data.frame()
    all_fr2s <- data.frame()
    best_lambdas <- data.frame()
    computation_times <- data.frame()

    cat("Doing beta ", beta.num, " (n-sensitivity, n_train=", n_train, ")\n")
    cat("     rep ", rep_num, "/", total_reps, "\n")


    # Generate data -----------------------------------------------------------
    cat("    -> generating data (n_train=", n_train, ")\n")

    gendata <- generate_fofr_data(nbasisX = gen_nbasis, nbasisY = gen_nbasis, nbeta = beta.num,
                                  nnodesX = nnodesX, nnodesY = nnodesY,
                                  n = n_train, n_validation = n_validation,
                                  X_process = "uniform_bspline")

    X_gen <- gendata$X
    Y <- gendata$Y
    beta_true <- gendata$beta_true

    X_gen_val <- gendata$X_val
    Y_val <- gendata$Y_val

    # No predictor observation noise here -- matches the X_sd_error = 0
    # setting the n=100 comparison run used:
    X <- X_gen
    X_val <- X_gen_val

    folds <- caret::createFolds(1:nrow(Y), k = num_folds)


    # CV - Ramsay & Silverman PENALIZED --------------------------------------

    cat("    -> starting Ramsay & Silverman penalized model\n")

    time_penalized_rs <- system.time({  # Measure time
      cv_penalized_rs <- cv_penalties_fregre.basis.fr(X,
                                                      Y,
                                                      argvals_X = argvals_X,
                                                      argvals_Y = argvals_Y,
                                                      num_bases_X = KK_rs,
                                                      num_bases_Y = LL_rs,
                                                      fda_basis_func_X = fda::create.bspline.basis,
                                                      fda_basis_func_Y = fda::create.bspline.basis,
                                                      penalty_X = penaltyvec_X_RS,
                                                      penalty_Y = penaltyvec_Y_RS,
                                                      folds = folds,
                                                      verbose = TRUE,
                                                      stripped = FALSE,
                                                      Lfd_X = fda::vec2Lfd(c(0, 0), range(argvals_X)),
                                                      Lfd_Y = fda::vec2Lfd(c(0, 0), range(argvals_Y))
      )
    })

    cat("    -> Ramsay & Silverman penalized done\n")

    best_lambdas <- rbind(
      best_lambdas,
      data.frame(lambda.penalty_X = cv_penalized_rs$best_penalties['penalty_X'],
                 lambda.penalty_Y = cv_penalized_rs$best_penalties['penalty_Y'],
                 method = "pFFR_RS",
                 nComp = NA,
                 beta.num = beta.num,
                 rep_num = rep_num)
    )

    all_CVEs <- rbind(
      all_CVEs,
      data.frame(
        CVE = as.numeric(  cv_penalized_rs$CVEs  ),
        method = "pFFR_RS",
        nComp = 1:max_nComp,
        beta.num = beta.num,
        rep_num = rep_num
      )
    )

    computation_times <- rbind(
      computation_times,
      data.frame(
        method = "pFFR_RS",
        beta.num = beta.num,
        nComp = 1:max_nComp,
        rep_num = rep_num,
        elapsed_time = time_penalized_rs[["elapsed"]]
      )
    )


    # CV - pFFPLS -------------------------------------------------------------

    cat("    -> starting pFFPLS\n")

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

      ## pFFPLS -----

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

      mean_imse_Y_val <- mean_imse_by_row(Y_val,
                                          predict(object = m_final, newdata = X_val)[, , nComp],
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


      ## Ramsay & Silverman (penalized) -----

      m_final_rs <- cv_penalized_rs$final_model

      beta_hat <- t(fda::eval.bifd(argvals_Y, argvals_X, m_final_rs$beta.est))

      beta_df <- expand.grid(q = argvals_Y, p = argvals_X)
      beta_df$z <- as.vector(beta_hat)
      beta_df$z_true <- as.vector(beta_true)
      beta_df$method <- "pFFR_RS"
      beta_df$nComp <- nComp
      beta_df$beta.num <- beta.num
      beta_df$rep_num <- rep_num

      y_val_pred <- predict_fregre_fr(object = m_final_rs,
                                      new.fdataobj =
                                        fda.usc::fdata(X_val,
                                                       argvals = argvals_X,
                                                       rangeval = range(argvals_X)) )$data

      all_beta_hats <- rbind(all_beta_hats, beta_df )

      mean_imse_Y_val <- mean_imse_by_row(Y_val,
                                          y_val_pred,
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
                               method = "pFFR_RS")
      )

      y_pred <- fitted.values(m_final_rs)$data

      all_fr2s <- rbind(all_fr2s,
                        tibble::tibble(
                          fR2_train = fR2_function(Y, y_pred),
                          fR2_val   = fR2_function(Y_val, y_val_pred),
                          nComp    = nComp,
                          beta.num = beta.num,
                          rep_num  = rep_num,
                          method   = "pFFR_RS")
      )

      rm(m_final_rs, y_pred, y_val_pred, beta_df, beta_hat, mean_imse_Y_val)

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

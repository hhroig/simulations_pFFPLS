## Standalone "buried PCA signal" simulation scenario -- see
## docs/new_scenario_buried_pca_signal_2026-09-11.md for the full rationale.
##
## X is generated directly from its own Karhunen-Loeve-type (Fourier)
## representation, with a deliberately NON-monotonic score-variance profile:
## a handful of "dominant" directions carry most of X's variance, while one
## additional, still-smooth direction ("buried_r") carries only a small
## fraction of it. The true beta loads ENTIRELY on that buried direction and
## nothing else -- by Fourier orthonormality this is exact, not approximate,
## so Y's dependence on X collapses onto exactly one, low-variance score,
## with zero leakage from the dominant (Y-irrelevant) directions.
##
## This tests whether pFFPLS's SUPERVISED direction-finding (it maximizes
## Cov(Xw, Y)) recovers that buried signal better than pFFR_RS, which
## regularizes only beta's smoothness and has no mechanism tied to X's own
## covariance structure at all.
##
## Does NOT touch simulations_fofr_v2_with_ivanescus_ramsay_silverman.R,
## main_simulations_call.R, pFFPLS_only_lambda_refine_sim.R, or
## n_sensitivity_sim.R. Fits all five methods (pFFR_RS, pFFR_I, pFFPLS,
## FFPLS, FFPLS_OB), mirroring simulations_fofr_v2_with_ivanescus_ramsay_silverman.R's
## structure block-for-block, so the same downstream analysis/plotting code
## works unchanged -- only generate_fofr_data() and redo_beta_true() differ.


# Local generate_fofr_data() / redo_beta_true() overrides for this scenario --

generate_fofr_data <- function(nbeta = "buried_pca",
                               nnodesX = 100, nnodesY = 100,
                               Rsq = 0.9, n = 100, n_validation = 50,
                               R_fourier = 10,
                               dominant_r = 1:4,
                               buried_r = 5,
                               dominant_sd = 1,
                               buried_sd = 0.5) {

  if (nbeta != "buried_pca") {
    stop("generate_fofr_data() in simulations_buried_pca_signal.R only ",
         "implements nbeta = 'buried_pca'.")
  }

  TX <- 1
  TY <- 1
  K_four <- 2 * R_fourier + 1  # const + R_fourier sin/cos pairs

  # period = TX (a FULL period over the observation window [0,TX]), not
  # 2*TX as in the existing "fourier_decay" X-process. Verified numerically:
  # with period = 2*TX, sin(a*pi*p) and cos(b*pi*p) cross-terms with (a+b)
  # odd do NOT integrate to zero over just [0,TX] (only sin-sin/cos-cos
  # pairs, and sin-cos pairs with (a+b) even, do) -- e.g.
  # integral_0^1 sin(5*pi*p)*cos(4*pi*p) dp = 0.35, not 0. That breaks the
  # exact-orthogonality guarantee this scenario depends on. With
  # period = TX, ALL cross-terms (sin-sin, cos-cos, sin-cos, any indices)
  # integrate to exactly zero over [0,TX], confirmed numerically before
  # this fix was applied.
  basisX_fourier <- fda::create.fourier.basis(rangeval = c(0, TX),
                                              nbasis = K_four,
                                              period  = TX)

  n <- n + n_validation

  # Row indexing (verified numerically against fda::eval.basis, NOT just
  # assumed from the existing "fourier_decay" code's comment -- that
  # comment describes the row PATTERN correctly but this function's
  # earlier draft mistranslated it into an off-by-one formula): row 1 =
  # constant, rows 2,4,6,... = sin(2*pi*1*p), sin(2*pi*2*p), ..., rows
  # 3,5,7,... = cos(2*pi*1*p), cos(2*pi*2*p), ... (frequency 2*pi*r*p
  # because period = TX here, not r*pi*p as it would be with period = 2*TX).
  row_sin <- function(r) 2 * r
  row_cos <- function(r) 2 * r + 1

  coefs <- matrix(0, nrow = K_four, ncol = n)
  for (r in dominant_r) {
    coefs[row_sin(r), ] <- stats::rnorm(n, mean = 0, sd = dominant_sd)
    coefs[row_cos(r), ] <- stats::rnorm(n, mean = 0, sd = dominant_sd)
  }
  # Only the sin term at buried_r carries signal-relevant variance; its
  # matching cos term, row 1 (constant), and every other r are left at
  # exactly 0 -- deliberately unused, zero-variance directions:
  coefs[row_sin(buried_r), ] <- stats::rnorm(n, mean = 0, sd = buried_sd)

  X <- fda::fd(coefs, basisX_fourier)

  p <- seq(0, TX, length.out = nnodesX)
  q <- seq(0, TY, length.out = nnodesY)

  # True beta lives ENTIRELY along the buried sin(2*pi*buried_r*p)
  # direction, matching row_sin(buried_r)'s actual basis function. With
  # period = TX, ALL pairs of distinct Fourier basis functions (sin-sin,
  # cos-cos, sin-cos, any indices) integrate to exactly zero over [0,TX] --
  # verified numerically -- so Y's dependence on X should collapse exactly
  # onto the buried score, with zero leakage from the dominant,
  # high-variance, signal-irrelevant directions.
  #
  # NOTE: this is deliberately computed via direct numerical integration
  # over the EVALUATED X and beta grids (trapezoidal rule), not via
  # fda::Data2fd()+fda::inprod() on fd objects. Data2fd() applies a small
  # default roughness penalty (lambda = 3e-8/range) that turned out to
  # distort the recovered coefficients badly for this Fourier basis
  # (verified empirically: projecting a pure sin(2*pi*5*p) test function
  # recovered a coefficient vector spread across many frequencies instead
  # of a single sparse entry) -- direct numerical integration sidesteps
  # that entirely and is exact up to grid resolution.
  f <- function(q, p) sin(2 * pi * q) * sin(2 * pi * buried_r * p)

  beta <- outer(q, p, f)  # nnodesY x nnodesX

  X_eval  <- t(fda::eval.fd(evalarg = p, fdobj = X))       # n x nnodesX
  Xc_eval <- sweep(X_eval, 2, colMeans(X_eval), "-")        # center across subjects, per p

  # Trapezoidal integration weights over the (uniform) p grid:
  dp <- diff(range(p)) / (nnodesX - 1)
  w  <- rep(dp, nnodesX)
  w[1] <- w[1] / 2
  w[nnodesX] <- w[nnodesX] / 2

  # Y_clean[i, q_ind] = integral_p Xc_eval[i, p] * beta(p, q_ind) dp:
  Y_clean <- Xc_eval %*% (t(beta) * w)  # n x nnodesY
  Y       <- Y_clean

  for (q_ind in 1:length(q)) {
    # var_e chosen so that R2(q) = 1 - var_e/var(Y_clean) matches the target
    # Rsq, consistent with the R2(q) definition used throughout this study:
    var_e       <- (1 - Rsq) * stats::var(Y_clean[, q_ind])
    Y[, q_ind]  <- Y_clean[, q_ind] +
      as.matrix(stats::rnorm(n = nrow(Y_clean), mean = 0, sd = sqrt(var_e)))
  }

  X         <- t(fda::eval.fd(evalarg = p, fdobj = X))
  beta_true <- beta
  Xc        <- Xc_eval  # already computed above, no eval.fd object needed

  n_traintest <- n - n_validation

  out <- list(
    Xc       = Xc[1:n_traintest, ],
    X        = X[1:n_traintest, ],
    Y        = Y[1:n_traintest, ],
    beta_true = beta_true,
    Xc_val   = Xc[(n_traintest + 1):nrow(Xc), ],
    X_val    = X[(n_traintest + 1):nrow(X), ],
    Y_val    = Y[(n_traintest + 1):nrow(X), ]
  )

  return(out)
}


redo_beta_true <- function(q, p, nbeta = "buried_pca", buried_r = 5) {
  f <- function(q, p) sin(2 * pi * q) * sin(2 * pi * buried_r * p)
  return(outer(q, p, f))
}


# Functional R2 (scalar: pools across all observations and time points) ---

fR2_function <- function(y_true, y_pred) {
  SST <- sum(sweep(y_true, 2, colMeans(y_true))^2)
  SSR <- sum((y_true - y_pred)^2)
  return(1 - SSR / SST)
}


# Mean IMSE by row --------------------------------------

mean_imse_by_row <- function(func_true,
                             func_hat,
                             argvals){

  interval_length <- (max(argvals) - min(argvals))

  beta_diff_sq <- (func_true - func_hat)^2

  num_obs <- nrow(func_true)

  imse_by_row <- array(NA, dim = num_obs)

  for (row_argval in 1:nrow(func_true)) {

    imse_by_row[row_argval] <- (num_int_1d(argvals = argvals,
                                           f_obs = beta_diff_sq[row_argval, ]))/interval_length
  }

  mean_imse <- mean(imse_by_row)

  return(mean_imse)

}



# Repetitions -------------------------------------------------------------


for(beta.num in num_betas)       {

  for(rep_num in rep_starts:total_reps)  {


    # Initialize savings:
    all_CVEs <- data.frame()
    all_beta_hats <- data.frame()
    all_final_res <- data.frame()
    all_fr2s <- data.frame()
    best_lambdas <- data.frame()
    computation_times <- data.frame()

    if(do_opt_bases_FFPLS){
      best_num_bases <- data.frame()
    }


    cat("Doing beta ", beta.num, " (buried PCA signal)\n")
    cat("     rep ", rep_num, "/", total_reps, "\n")


    # Generate data -----------------------------------------------------------
    cat("    -> generating data\n")

    gendata <- generate_fofr_data(nbeta = beta.num,
                                  nnodesX = nnodesX, nnodesY = nnodesY,
                                  n = 100, n_validation = 50,
                                  R_fourier = R_fourier,
                                  dominant_r = dominant_r,
                                  buried_r = buried_r,
                                  dominant_sd = dominant_sd,
                                  buried_sd = buried_sd)

    X_gen <- gendata$X
    Y <- gendata$Y
    beta_true <- gendata$beta_true

    X_gen_val <- gendata$X_val
    Y_val <- gendata$Y_val

    # No predictor observation noise in this scenario -- the interesting
    # part is X's own covariance structure, not measurement error:
    X <- X_gen
    X_val <- X_gen_val

    folds <- caret::createFolds(1:nrow(Y), k = num_folds)


    # CV - Ramsay & Silverman PENALIZED ------------------------------------------------------------------

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



    # CV - Ivanescu's pffr --------------------------------

    cat("    -> starting Ivanescu's pffr\n")

    data_pffr <- list(X = X, Y = Y, argvals_X = argvals_X, argvals_Y = argvals_Y)

    time_pffr <- system.time({  # Measure time

      m_final_pffr <- pffr(
        Y ~ 0 + ff(X, xind = argvals_X,
                   splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 2)), k = c(KK, LL))),
        yind=argvals_Y,
        data=data_pffr)
    })

    cat("    -> Ivanescu's pffr model done\n")

    computation_times <- rbind(
      computation_times,
      data.frame(
        method = "pFFR_I",
        beta.num = beta.num,
        nComp = 1:max_nComp,
        rep_num = rep_num,
        elapsed_time = time_pffr[["elapsed"]]
      )
    )

    # CV - pFFPLS ------------------------------------------------------------------

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
      data.frame(lambda = cv_penalized$best_penalties,
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




    # CV - FFPLS (no penalties) ---------------------------------------------------

    cat("    -> starting FFPLS model\n")

    time_nonpen <- system.time({  # Measure time
      cv_nonpen <- cv_unique_fof_par(X = X,
                                     Y = Y,
                                     argvals_X = argvals_X,
                                     argvals_Y = argvals_Y,
                                     ncomp = max_nComp,
                                     center = center,
                                     folds = folds,
                                     basisobj_X = basisobj_X,
                                     basisobj_Y = basisobj_Y,
                                     penaltyvec_X = 0,
                                     penaltyvec_Y = 0,
                                     verbose = TRUE,
                                     stripped = TRUE,
                                     maxit = 100000 )
    })



    cat("    -> FFPLS model done\n")

    all_CVEs <- rbind(
      all_CVEs,
      data.frame(
        CVE = cv_nonpen$CVEs_ncomp,
        method = "FFPLS",
        nComp = 1:max_nComp,
        beta.num = beta.num,
        rep_num = rep_num
      )
    )

    computation_times <- rbind(
      computation_times,
      data.frame(
        method = "FFPLS",
        beta.num = beta.num,
        nComp = 1:max_nComp,
        rep_num = rep_num,
        elapsed_time = time_nonpen[["elapsed"]]
      )
    )


    # CV - FFPS basis optimization (no penalties) ---------------------------------------------------

    if (do_opt_bases_FFPLS) {

      cat("    -> starting FFPLS with optimal bases\n")

      time_nonpen_bases <- system.time({  # Measure time
        cv_nonpen_bases <- cv_bases_fof_par(X = X,
                                            Y = Y,
                                            argvals_X = argvals_X,
                                            argvals_Y = argvals_Y,
                                            ncomp = max_nComp,
                                            center = center,
                                            folds = folds,

                                            num_bases_X = KK_list,
                                            num_bases_Y = LL_list,
                                            fda_basis_func_X = fda::create.bspline.basis,
                                            fda_basis_func_Y = fda::create.bspline.basis,
                                            penalty_X = 0,
                                            penalty_Y = 0,

                                            verbose = TRUE,
                                            stripped = TRUE,
                                            maxit = 100000 )
      })


      cat("    -> FFPLS with optimal bases done\n")

      best_num_bases <- rbind(
        best_num_bases,
        data.frame(num_bases = cv_nonpen_bases$best_num_bases,
                   method = "FFPLS_OB",
                   nComp = 1:max_nComp,
                   beta.num = beta.num,
                   rep_num = rep_num)
      )

      all_CVEs <- rbind(
        all_CVEs,
        data.frame(
          CVE = cv_nonpen_bases$CVEs_ncomp,
          method = "FFPLS_OB",
          nComp = 1:max_nComp,
          beta.num = beta.num,
          rep_num = rep_num
        )
      )


      computation_times <- rbind(
        computation_times,
        data.frame(
          method = "FFPLS_OB",
          beta.num = beta.num,
          nComp = 1:max_nComp,
          rep_num = rep_num,
          elapsed_time = time_nonpen_bases[["elapsed"]]
        )
      )

    }



    # Final models for each ncomp ---------------------------------------------


    for (nComp in 1:max_nComp) {
      cat("       ---> final models cmpt:", nComp, "/", max_nComp, "\n")

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


      ## FFPLS ----

      m_final <- ffpls_bs(X = X,
                          Y = Y,
                          argvals_X = argvals_X,
                          argvals_Y = argvals_Y,
                          ncomp = nComp,
                          center = center,
                          basisobj_X = basisobj_X,
                          basisobj_Y = basisobj_Y,
                          penalty_X = 0,
                          penalty_Y = 0,
                          verbose = FALSE,
                          stripped = F,
                          maxit = 100000)

      beta_hat <- m_final$coefficient_function[, , nComp]

      beta_df <- expand.grid(q = argvals_Y, p = argvals_X)
      beta_df$z <- as.vector(beta_hat)
      beta_df$z_true <- as.vector(beta_true)
      beta_df$method <- "FFPLS"
      beta_df$nComp <- nComp
      beta_df$beta.num <- beta.num
      beta_df$rep_num <- rep_num

      all_beta_hats <-  rbind(all_beta_hats, beta_df )

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
                               method = "FFPLS")
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
                          method   = "FFPLS")
      )

      rm(m_final, y_pred, y_val_pred, beta_df, beta_hat, mean_imse_Y_val)


      ## FFPLS optimal bases ----

      if (do_opt_bases_FFPLS) {

        basisobj_X_best <- fda::create.bspline.basis(rangeval = range(argvals_X),
                                                     nbasis = cv_nonpen_bases$best_num_bases[nComp, "numbases_X"])
        basisobj_Y_best <- fda::create.bspline.basis(rangeval = range(argvals_Y),
                                                     nbasis = cv_nonpen_bases$best_num_bases[nComp, "numbases_Y"])

        m_final <- ffpls_bs(X = X,
                            Y = Y,
                            argvals_X = argvals_X,
                            argvals_Y = argvals_Y,
                            ncomp = nComp,
                            center = center,
                            basisobj_X = basisobj_X_best,
                            basisobj_Y = basisobj_Y_best,
                            penalty_X = 0,
                            penalty_Y = 0,
                            verbose = FALSE,
                            stripped = F,
                            maxit = 100000)

        beta_hat <- m_final$coefficient_function[, , nComp]

        beta_df <- expand.grid(q = argvals_Y, p = argvals_X)
        beta_df$z <- as.vector(beta_hat)
        beta_df$z_true <- as.vector(beta_true)
        beta_df$method <- "FFPLS_OB"
        beta_df$nComp <- nComp
        beta_df$beta.num <- beta.num
        beta_df$rep_num <- rep_num

        all_beta_hats <-  rbind(all_beta_hats, beta_df )

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
                                 method = "FFPLS_OB")
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
                            method   = "FFPLS_OB")
        )

        rm(m_final, y_pred, y_val_pred, beta_df, beta_hat, mean_imse_Y_val)

      }


      ## Ivanescu's pffr --------------------------------

      beta_hat <- coef(m_final_pffr)[["smterms"]][[1]]
      beta_true_rebuilt <- redo_beta_true(q = beta_hat$y, p = beta_hat$x, nbeta = beta.num, buried_r = buried_r)

      beta_df <- expand.grid(q = beta_hat$y, p = beta_hat$x)
      beta_df$z <- beta_hat$value

      # reshape to q times p matrix:
      beta_df$z <- as.vector( t(reshape2::acast(beta_df, q~p, value.var = "z"))  )

      beta_df$z_true <- as.vector(beta_true_rebuilt)
      beta_df$method <- "pFFR_I"
      beta_df$nComp <- nComp
      beta_df$beta.num <- beta.num
      beta_df$rep_num <- rep_num

      all_beta_hats <- rbind(all_beta_hats, beta_df )

      mean_imse_Y_val <- mean_imse_by_row(Y_val,
                                          predict(m_final_pffr, newdata = list(X = X_val)),
                                          argvals_Y)

      all_final_res <- rbind(all_final_res,
                             tibble::tibble(
                               imse = penFoFPLS::imse_beta_ffpls_bs(beta_true_rebuilt,
                                                                    reshape2::acast(beta_df, q~p, value.var = "z"),
                                                                    beta_hat$x,
                                                                    beta_hat$y),
                               mean_imse_Y_val = mean_imse_Y_val,
                               nComp = nComp,
                               beta.num = beta.num,
                               rep_num = rep_num,
                               method = "pFFR_I")
      )

      y_pred <- predict(m_final_pffr)
      y_val_pred <- predict(m_final_pffr, newdata = list(X = X_val))

      all_fr2s <- rbind(all_fr2s,
                        tibble::tibble(
                          fR2_train = fR2_function(Y, y_pred),
                          fR2_val   = fR2_function(Y_val, y_val_pred),
                          nComp    = nComp,
                          beta.num = beta.num,
                          rep_num  = rep_num,
                          method   = "pFFR_I")
      )

      rm(y_pred, y_val_pred, beta_df, beta_hat, mean_imse_Y_val)


      ## Ramsay & Silverman (penalized) ----

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

    if (do_opt_bases_FFPLS) {
      saveRDS(best_num_bases, file =
                paste0(out_folder, "best_num_bases_FFPLS_rep_", rep_num, "_beta_", beta.num, ".Rds"))
    }

    saveRDS(all_fr2s, file =
              paste0(out_folder, "fR2s_rep_", rep_num, "_beta_", beta.num, ".Rds"))

    saveRDS(computation_times, file =
              paste0(out_folder, "computation_times_rep_", rep_num, "_beta_", beta.num, ".Rds"))


  } # end nested inner rep_num

} # end nested outer beta_num

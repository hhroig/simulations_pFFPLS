## Standalone "buried PCA signal, SHARP/LOCALIZED shape" simulation scenario --
## see docs/new_scenario_buried_pca_sharp_signal_2026-09-28.md for the full
## rationale. Extends simulations_buried_pca_signal.R (which buries a single
## pure sinusoid in one low-variance X direction): here the true beta's
## p-shape is a genuinely SHARP, LOCALIZED pulse centered at p=0.5 (per
## Reviewer 1's ask #1: "sharp local peaks ... or rapid variation in one
## direction"), built from a CLUSTER of low-variance ("buried") Fourier
## cosine directions of X rather than one single frequency.
##
## Key construction fact (verified numerically before use): for
## fda::create.fourier.basis(period = TX = 1), cos(2*pi*r*p) is symmetric
## under reflection p -> 1-p about p=0.5 for EVERY integer r. Choosing
## coefficients c_r = w_r * (-1)^r (all terms in phase at p=0.5) makes
## sum_r c_r*cos(2*pi*r*p) a real, localized, band-limited pulse centered at
## p=0.5 with decaying sidelobes -- built ENTIRELY out of cosine directions,
## so it can be made to live EXACTLY inside a chosen low-variance direction
## cluster of X (only the COS term at each "bump" frequency is given
## nonzero variance; the matching SIN term, and every dominant/unused
## direction, stay at exactly zero -- mirroring the original scenario's
## "only sin(buried_r) is nonzero" convention, just transposed to cos and
## widened to a 6-frequency cluster).
##
## Does NOT touch simulations_fofr_v2_with_ivanescus_ramsay_silverman.R,
## main_simulations_call.R, pFFPLS_only_lambda_refine_sim.R,
## n_sensitivity_sim.R, simulations_buried_pca_signal.R, or
## main_buried_pca_call.R. Fits all five methods (pFFR_RS, pFFR_I, pFFPLS,
## FFPLS, FFPLS_OB), mirroring simulations_buried_pca_signal.R's structure
## block-for-block, so the same downstream analysis/plotting code works
## unchanged -- only generate_fofr_data() and redo_beta_true() differ.


# Local generate_fofr_data() / redo_beta_true() overrides for this scenario --

generate_fofr_data <- function(nbeta = "buried_pca_sharp",
                               nnodesX = 100, nnodesY = 100,
                               Rsq = 0.9, n = 100, n_validation = 50,
                               R_fourier = 10,
                               dominant_r = 1:4,
                               bump_r = 5:10,
                               dominant_sd = 1,
                               buried_sd = 0.3) {

  if (nbeta != "buried_pca_sharp") {
    stop("generate_fofr_data() in simulations_buried_pca_sharp_signal.R only ",
         "implements nbeta = 'buried_pca_sharp'.")
  }

  TX <- 1
  TY <- 1
  K_four <- 2 * R_fourier + 1  # const + R_fourier sin/cos pairs

  # period = TX (full period over [0,TX]) -- same fix already validated in
  # simulations_buried_pca_signal.R; required for exact cross-term
  # orthogonality of every sin/cos pair, which this scenario also depends on.
  basisX_fourier <- fda::create.fourier.basis(rangeval = c(0, TX),
                                              nbasis = K_four,
                                              period  = TX)

  n <- n + n_validation

  # Row indexing (same convention as simulations_buried_pca_signal.R,
  # verified there against fda::eval.basis): row 1 = constant, rows
  # 2,4,6,... = sin(2*pi*1*p), sin(2*pi*2*p), ..., rows 3,5,7,... =
  # cos(2*pi*1*p), cos(2*pi*2*p), ...
  row_sin <- function(r) 2 * r
  row_cos <- function(r) 2 * r + 1

  # Hann-tapered, phase-aligned bump weights: w_r tapers the ends of the
  # frequency cluster to zero (suppresses far-field ringing relative to flat
  # weights -- verified visually before committing, see
  # docs/new_scenario_buried_pca_sharp_signal_2026-09-28.md), and (-1)^r
  # aligns every term's phase at p=0.5 so they add constructively there
  # instead of canceling.
  n_bump <- length(bump_r)
  w_hann <- 0.5 * (1 - cos(2 * pi * (seq_len(n_bump) - 1) / (n_bump - 1)))
  c_bump <- w_hann * (-1)^bump_r

  coefs <- matrix(0, nrow = K_four, ncol = n)
  for (r in dominant_r) {
    coefs[row_sin(r), ] <- stats::rnorm(n, mean = 0, sd = dominant_sd)
    coefs[row_cos(r), ] <- stats::rnorm(n, mean = 0, sd = dominant_sd)
  }
  # Only the COS term at each bump_r frequency carries signal-relevant
  # variance; the matching SIN term at those same frequencies, the constant
  # (row 1), and every r outside dominant_r/bump_r are left at exactly 0 --
  # deliberately unused, zero-variance directions:
  for (r in bump_r) {
    coefs[row_cos(r), ] <- stats::rnorm(n, mean = 0, sd = buried_sd)
  }

  X <- fda::fd(coefs, basisX_fourier)

  p <- seq(0, TX, length.out = nnodesX)
  q <- seq(0, TY, length.out = nnodesY)

  # True beta's p-shape lives ENTIRELY in the span of {cos(2*pi*r*p): r in
  # bump_r} -- by the exact cross-term orthogonality of this basis (period =
  # TX), Y's dependence on X should collapse exactly onto the bump_r COS
  # scores, with zero leakage from the dominant, high-variance,
  # signal-irrelevant directions. Verified numerically before use (see
  # scenario doc): max |correlation| with any dominant score was ~0.04,
  # vs. up to ~0.66 with the bump scores.
  #
  # Computed via direct numerical (trapezoidal) integration over the
  # EVALUATED X and beta grids, not fda::Data2fd()+fda::inprod() -- same
  # reasoning as simulations_buried_pca_signal.R: Data2fd()'s default
  # roughness penalty measurably distorts recovered coefficients for a
  # basis like this one.
  beta_p <- function(pp) {
    z <- rep(0, length(pp))
    for (i in seq_along(bump_r)) {
      z <- z + c_bump[i] * cos(2 * pi * bump_r[i] * pp)
    }
    z
  }
  f <- function(q, p) sin(2 * pi * q) * beta_p(p)

  beta <- outer(q, p, f)  # nnodesY x nnodesX

  X_eval  <- t(fda::eval.fd(evalarg = p, fdobj = X))       # n x nnodesX
  Xc_eval <- sweep(X_eval, 2, colMeans(X_eval), "-")        # center across subjects, per p

  dp <- diff(range(p)) / (nnodesX - 1)
  w  <- rep(dp, nnodesX)
  w[1] <- w[1] / 2
  w[nnodesX] <- w[nnodesX] / 2

  Y_clean <- Xc_eval %*% (t(beta) * w)  # n x nnodesY
  Y       <- Y_clean

  for (q_ind in 1:length(q)) {
    var_e       <- (1 - Rsq) * stats::var(Y_clean[, q_ind])
    Y[, q_ind]  <- Y_clean[, q_ind] +
      as.matrix(stats::rnorm(n = nrow(Y_clean), mean = 0, sd = sqrt(var_e)))
  }

  X         <- t(fda::eval.fd(evalarg = p, fdobj = X))
  beta_true <- beta
  Xc        <- Xc_eval

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


redo_beta_true <- function(q, p, nbeta = "buried_pca_sharp", bump_r = 5:10) {
  n_bump <- length(bump_r)
  w_hann <- 0.5 * (1 - cos(2 * pi * (seq_len(n_bump) - 1) / (n_bump - 1)))
  c_bump <- w_hann * (-1)^bump_r
  beta_p <- function(pp) {
    z <- rep(0, length(pp))
    for (i in seq_along(bump_r)) z <- z + c_bump[i] * cos(2 * pi * bump_r[i] * pp)
    z
  }
  f <- function(q, p) sin(2 * pi * q) * beta_p(p)
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


    cat("Doing beta ", beta.num, " (buried PCA signal, sharp)\n")
    cat("     rep ", rep_num, "/", total_reps, "\n")


    # Generate data -----------------------------------------------------------
    cat("    -> generating data\n")

    gendata <- generate_fofr_data(nbeta = beta.num,
                                  nnodesX = nnodesX, nnodesY = nnodesY,
                                  Rsq = Rsq,
                                  n = 100, n_validation = 50,
                                  R_fourier = R_fourier,
                                  dominant_r = dominant_r,
                                  bump_r = bump_r,
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
      beta_true_rebuilt <- redo_beta_true(q = beta_hat$y, p = beta_hat$x, nbeta = beta.num, bump_r = bump_r)

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

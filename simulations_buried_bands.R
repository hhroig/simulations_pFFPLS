## Buried-bands simulation scenario.
##
## Data-generating process (one dataset per repetition):
##   X_i(p) = dominant part + buried part, p in [0, 1]
##     dominant part: sum_{j=1}^{10} 0.5 u_ij f_j(p), f_j(p) = (sin(0.3 j pi p) + 1) cos(0.3 j pi p),
##                    u_ij ~ N(0, 1): smooth, strongly collinear, ~97% of X's variance;
##     buried part:   b_i1 g1~(p) + b_i2 g2~(p), two Gaussian bands centred at p = 0.3 and 0.7
##                    (width band_sd), orthogonalised against the top-4 eigenfunctions of the
##                    dominant part; together they carry band_share (3%) of X's variance.
##   beta(q, p) = [h1(q) g1~(p) - h2(q) g2~(p)] / C, h_k(q) = exp(-(q - c_k)^2 / 0.01),
##     normalised to max |beta| = 1: two opposite sharp peaks, at (0.3, 0.3) and (0.7, 0.7).
##   Y_i(q) = int (X_i(p) - mean X(p)) beta(q, p) dp + eps_iq, Var(eps_q) = (1 - Rsq) Var(Y_clean(q)).
##   The methods observe X_i(p) + N(0, X_sd_error^2) noise at every grid point; Y comes from the
##   noise-free curves.
##
## This file is sourced by run_simulation.R (inside main_buried_bands_call(), with local = TRUE),
## which defines the configuration variables used below. It fits all five methods
## (pFFR_RS, pFFR_I, pFFPLS, FFPLS, FFPLS_OB) in every repetition and saves one set of .Rds
## files per repetition to `out_folder`.


# Band functions (closed form, evaluable on any grid) -----------------------

make_band_functions <- function(band_sd = 0.035, centers = c(0.3, 0.7), n_top = 4) {
  h_sd <- sqrt(0.005)   # width of the q-profiles: h_k(q) = exp(-(q - c_k)^2 / 0.01)
  fj <- function(pp) sapply(1:10, function(j) (sin(0.3 * j * pi * pp) + 1) * cos(0.3 * j * pi * pp))
  # Fine trapezoidal grid for the Gram matrix / inner products
  pf <- seq(0, 1, length.out = 4001)
  wf <- rep(pf[2] - pf[1], length(pf)); wf[c(1, length(pf))] <- wf[1] / 2
  Ff <- fj(pf)
  G  <- t(Ff) %*% (Ff * wf)                      # Gram matrix <f_j, f_k>
  # Dominant covariance in span{f_j}: C psi = 0.25 F G a, psi = F a
  # => eigen-problem 0.25 G a = lambda a, normalised so that a' G a = 1
  eg <- eigen(0.25 * G, symmetric = TRUE)
  A  <- eg$vectors[, 1:n_top, drop = FALSE]
  A  <- sweep(A, 2, sqrt(colSums(A * (G %*% A))), "/")
  psi_f <- Ff %*% A                              # eigenfunctions on the fine grid
  g_raw <- function(pp, k) exp(-(pp - centers[k])^2 / (2 * band_sd^2))
  coef_k <- lapply(1:2, function(k) colSums(psi_f * (g_raw(pf, k) * wf)))
  # bands with their projection on the top-4 dominant eigenfunctions removed:
  g_tilde <- lapply(1:2, function(k) {
    force(k)
    function(pp) as.vector(g_raw(pp, k) - fj(pp) %*% A %*% coef_k[[k]])
  })
  h <- list(function(qq)  exp(-(qq - centers[1])^2 / (2 * h_sd^2)),
            function(qq) -exp(-(qq - centers[2])^2 / (2 * h_sd^2)))
  beta_raw <- function(qq, pp) h[[1]](qq) %o% g_tilde[[1]](pp) + h[[2]](qq) %o% g_tilde[[2]](pp)
  C_norm <- max(abs(beta_raw(pf[seq(1, 4001, 4)], pf[seq(1, 4001, 4)])))
  list(fj = fj, g_tilde = g_tilde,
       beta = function(qq, pp) beta_raw(qq, pp) / C_norm,
       dom_var = 0.25 * sum(diag(G)),
       g_norm2 = sapply(g_tilde, function(g) sum(wf * g(pf)^2)))
}


# Data generation --------------------------------------------------------------

generate_fofr_data <- function(nbeta = "buried_bands",
                               nnodesX = 100, nnodesY = 100,
                               Rsq = 0.9, n = 100, n_validation = 50,
                               band_sd = 0.035, band_share = 0.03) {

  if (nbeta != "buried_bands") {
    stop("generate_fofr_data() in simulations_buried_bands.R only implements nbeta = 'buried_bands'.")
  }

  bf <- make_band_functions(band_sd = band_sd)
  n  <- n + n_validation
  p  <- seq(0, 1, length.out = nnodesX)
  q  <- seq(0, 1, length.out = nnodesY)

  # bands jointly carry `band_share` of X's total integrated variance
  sd_b <- sqrt(band_share / (1 - band_share) * bf$dom_var / sum(bf$g_norm2))

  U  <- matrix(stats::rnorm(n * 10), nrow = n)
  Bs <- matrix(stats::rnorm(n * 2, sd = sd_b), nrow = n)
  X  <- U %*% t(0.5 * bf$fj(p)) +
        Bs[, 1] %o% bf$g_tilde[[1]](p) + Bs[, 2] %o% bf$g_tilde[[2]](p)   # n x nnodesX

  beta <- bf$beta(q, p)                                                    # nnodesY x nnodesX

  Xc <- sweep(X, 2, colMeans(X), "-")
  dp <- diff(range(p)) / (nnodesX - 1)
  w  <- rep(dp, nnodesX); w[c(1, nnodesX)] <- dp / 2
  Y_clean <- Xc %*% (t(beta) * w)                                          # trapezoidal rule
  Y <- Y_clean
  for (q_ind in seq_along(q)) {
    var_e <- (1 - Rsq) * stats::var(Y_clean[, q_ind])
    Y[, q_ind] <- Y_clean[, q_ind] + stats::rnorm(n, mean = 0, sd = sqrt(var_e))
  }

  n_traintest <- n - n_validation
  list(Xc        = Xc[1:n_traintest, ],
       X         = X[1:n_traintest, ],
       Y         = Y[1:n_traintest, ],
       beta_true = beta,
       Xc_val    = Xc[(n_traintest + 1):n, ],
       X_val     = X[(n_traintest + 1):n, ],
       Y_val     = Y[(n_traintest + 1):n, ])
}


# True beta on an arbitrary grid (pFFR_I returns its estimate on its own grid)
redo_beta_true <- function(q, p, nbeta = "buried_bands", band_sd = 0.035) {
  make_band_functions(band_sd = band_sd)$beta(q, p)
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

  for(rep_num in rep_starts:rep_end)  {


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


    cat("Doing beta ", beta.num, " (buried bands)\n")
    cat("     rep ", rep_num, "/", total_reps, "\n")


    # Generate data -----------------------------------------------------------
    cat("    -> generating data\n")

    # One fixed seed per repetition, so runs with different method settings see
    # identical datasets (paired comparisons across runs):
    if (!is.null(seed_base)) set.seed(seed_base + rep_num)

    gendata <- generate_fofr_data(nbeta = beta.num,
                                  nnodesX = nnodesX, nnodesY = nnodesY,
                                  Rsq = Rsq,
                                  n = 100, n_validation = 50,
                                  band_sd = band_sd,
                                  band_share = band_share)

    X_gen <- gendata$X
    Y <- gendata$Y
    beta_true <- gendata$beta_true

    X_gen_val <- gendata$X_val
    Y_val <- gendata$Y_val

    # Measurement noise on X: one independent draw per (subject, grid point),
    # added AFTER Y was generated from the noise-free curves:
    X     <- X_gen     + matrix(stats::rnorm(length(X_gen),     mean = 0, sd = X_sd_error), nrow = nrow(X_gen))
    X_val <- X_gen_val + matrix(stats::rnorm(length(X_gen_val), mean = 0, sd = X_sd_error), nrow = nrow(X_gen_val))

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



    if (do_pffr) {
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

    }  # end if (do_pffr)
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


      if (do_pffr) {
      ## Ivanescu's pffr --------------------------------

      beta_hat <- coef(m_final_pffr)[["smterms"]][[1]]
      beta_true_rebuilt <- redo_beta_true(q = beta_hat$y, p = beta_hat$x, nbeta = beta.num, band_sd = band_sd)

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


      }  # end if (do_pffr)
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

## Run: ?fregre.basis.fr
## Then, code a CV method for this penalized approach by Ramsay and Silverman
## Adapt the same method we have for cv_unique_fof_par to this case...


# Re-do beta function for Ivanescu's --------------------------------------

redo_beta_true <- function(q, p, nbeta) {
  
  if (nbeta == 1) {
    f <- function(q, p) {
      return((q - p)^2)
    }
  }
  else if (nbeta == 2) {
    f <- function(q, p) {
      return(5 * exp(-((p - 0.75)^2 + (q - 0.75)^2)/(2 * 
                                                       0.2^2)))
    }
  }
  else if (nbeta == 3) {
    f <- function(q, p) {
      return(((p * 4 - 2)^3 - 3 * (p * 4 - 2) * ((q * 
                                                    4 - 2)^2)))
    }
  }
  else if (nbeta == 4) {
    f <- function(q, p) {
      return(5 * exp(-((p - 0.75)^2 + (q - 0.75)^2)/(2 * 0.25^2)) + 
               5 * exp(-((p - 0.1)^2 + (q - 0.1)^2)/(2 * 0.25^2)))
    }
  }
  
  return(outer(q, p, f))
}



# R2 function on q --------------------------------------

R_sqr_function <- function(y_true, y_pred) {
  
  y_mean <- colMeans(y_true)
  
  SST <- colMeans((y_true - y_mean)^2)
  SSR <- colMeans((y_true - y_pred)^2)  
  
  R2 <- 1 - SSR/SST
  
  return(R2)
  
  
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
    all_val_MSEs <- data.frame()
    all_beta_hats <- data.frame()
    all_final_res <- data.frame()
    all_r2s <- data.frame()
    best_lambdas <- data.frame()
    computation_times <- data.frame()  
    
    if(do_opt_bases_FFPLS){
      best_num_bases <- data.frame()
    }
    
    
    cat("Doing beta ", beta.num, "\n")
    cat("     rep ", rep_num, "/", total_reps, "\n")
    
    
    # Generate data -----------------------------------------------------------
    cat("    -> generating data\n")
    
    # Fix to 20 the number of basis for the generated data:
    gendata <- generate_fofr_data(nbasisX = 20, nbasisY = 20, nbeta = beta.num,
                                  nnodesX = nnodesX, nnodesY = nnodesY, 
                                  n = 100, n_validation = 50)
    
    X_gen <- gendata$X
    Y <- gendata$Y
    beta_true <- gendata$beta_true
    
    X_gen_val <- gendata$X_val
    Y_val <- gendata$Y_val
    
    if (X_sd_error > 0) {
      X <- X_gen + stats::rnorm(nrow(X_gen), mean = 0, sd = X_sd_error)
      X_val <- X_gen_val + stats::rnorm(nrow(X_gen_val), mean = 0, sd = X_sd_error)
    }else {
      X <- X_gen
      X_val <- X_gen_val
    }
    
    
    
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
      # Repeat the same unique CVE for all nComp (just for comparison purposes):
      data.frame(
        CVE = as.numeric(  cv_penalized_rs$CVEs  ),
        method = "pFFR_RS",
        nComp = 1:max_nComp,
        beta.num = beta.num,
        rep_num = rep_num
      )
    )
    
    # Record computation time
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
        Y ~ 0 + ff(X, xind=argvals_X), 
        bs.yindex = list(bs="ps",    # penalized splines bases
                         k=LL,       # number of basis for Y defined in main.R
                         m=c(2, 2)   # cubic B-splines bases with second order difference penalty
        ), 
        yind=argvals_Y, 
        data=data_pffr)
    })
    
    cat("    -> Ivanescu's pffr model done\n")
    
    # Record computation time
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
                                        stripped = FALSE,
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
    
    # Record computation time
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
                                     stripped = FALSE,
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
    
    # Record computation time
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
                                            stripped = FALSE,
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
      
      
      # Record computation time
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
      
      
      # Compute the average IMSE for the validation set:
      mean_imse_Y_val <- mean_imse_by_row(Y_val, # true
                                          predict(object = m_final, newdata = X_val)[, , nComp], # predicted
                                          argvals_Y)
      
      
      all_final_res <- rbind(all_final_res,
                             tibble(
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
      
      
      
      
      # Compute R2 for the training and validation sets:
      y_pred <- fitted.values(m_final)[ , , nComp]
      y_val_pred <- predict(object = m_final, newdata = X_val)[, , nComp]
      
      all_r2s <- rbind(all_r2s,
                       tibble(
                         R2_train = R_sqr_function(Y, y_pred),
                         R2_val = R_sqr_function(Y_val, y_val_pred),
                         q = argvals_Y,
                         nComp = nComp,
                         beta.num = beta.num,
                         rep_num = rep_num,
                         method = "pFFPLS")
      )
      
      
      rm(m_final, y_pred, y_val_pred, beta_df, beta_hat, mean_imse_Y_val) # I'll reuse the same model name
      
      
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
      
      # Compute the average IMSE for the validation set:
      mean_imse_Y_val <- mean_imse_by_row(Y_val, # true
                                          predict(object = m_final, newdata = X_val)[, , nComp], # predicted
                                          argvals_Y)
      
      all_final_res <- rbind(all_final_res,
                             tibble(
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
      
      
      # Compute R2 for the training and validation sets:
      y_pred <- fitted.values(m_final)[ , , nComp]
      y_val_pred <- predict(object = m_final, newdata = X_val)[, , nComp]
      
      all_r2s <- rbind(all_r2s,
                       tibble(
                         R2_train = R_sqr_function(Y, y_pred),
                         R2_val = R_sqr_function(Y_val, y_val_pred),
                         q = argvals_Y,
                         nComp = nComp,
                         beta.num = beta.num,
                         rep_num = rep_num,
                         method = "FFPLS")
      )
      
      
      rm(m_final, y_pred, y_val_pred, beta_df, beta_hat, mean_imse_Y_val) # I'll reuse the same model name
      
      
      ## FFPLS optimal bases ----
      
      if (do_opt_bases_FFPLS) {
        
        # B-spline basis:
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
        
        # Compute the average IMSE for the validation set:
        mean_imse_Y_val <- mean_imse_by_row(Y_val, # true
                                            predict(object = m_final, newdata = X_val)[, , nComp], # predicted
                                            argvals_Y)
        
        all_final_res <- rbind(all_final_res,
                               tibble(
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
        
        
        # Compute R2 for the training and validation sets:
        y_pred <- fitted.values(m_final)[ , , nComp]
        y_val_pred <- predict(object = m_final, newdata = X_val)[, , nComp]
        
        all_r2s <- rbind(all_r2s,
                         tibble(
                           R2_train = R_sqr_function(Y, y_pred),
                           R2_val = R_sqr_function(Y_val, y_val_pred),
                           q = argvals_Y,
                           nComp = nComp,
                           beta.num = beta.num,
                           rep_num = rep_num,
                           method = "FFPLS_OB")
        )
        
        
        rm(m_final, y_pred, y_val_pred, beta_df, beta_hat, mean_imse_Y_val) # I'll reuse the same model name
        
      }
      
      
      ## Ivanescu's pffr --------------------------------
      
      beta_hat <- coef(m_final_pffr)[["smterms"]][["ff(X,argvals_X)"]]
      beta_true_rebuilt <- redo_beta_true(q = beta_hat$y, p = beta_hat$x, nbeta = beta.num)
      
      beta_df <- expand.grid(q = beta_hat$y, p = beta_hat$x) 
      beta_df$z <- beta_hat$value
      
      # reshape to q times p matrix:
      beta_df$z <- as.vector( t(reshape2::acast(beta_df, q~p, value.var = "z"))  )
      
      beta_df$z_true <- as.vector(beta_true_rebuilt)
      beta_df$method <- "pFFR_I"
      beta_df$nComp <- nComp       # No actual components, leave it here for comparison purposes!
      beta_df$beta.num <- beta.num
      beta_df$rep_num <- rep_num
      
      all_beta_hats <- rbind(all_beta_hats, beta_df )
      
      # Compute the average IMSE for the validation set:
      mean_imse_Y_val <- mean_imse_by_row(Y_val, # true
                                          predict(m_final_pffr, newdata = list(X = X_val)), # predicted
                                          argvals_Y)
      
      all_final_res <- rbind(all_final_res,
                             tibble(
                               imse = penFoFPLS::imse_beta_ffpls_bs(beta_true_rebuilt,
                                                                    reshape2::acast(beta_df, q~p, value.var = "z"),
                                                                    beta_hat$x,
                                                                    beta_hat$y),
                               mean_imse_Y_val = mean_imse_Y_val,
                               nComp = nComp,       # No actual components, leave it here for comparison purposes!
                               beta.num = beta.num,
                               rep_num = rep_num,
                               method = "pFFR_I")
      )
      
      
      
      # Compute R2 for the training and validation sets:
      y_pred <- predict(m_final_pffr)
      y_val_pred <- predict(m_final_pffr, newdata = list(X = X_val))
      
      all_r2s <- rbind(all_r2s,
                       tibble(
                         R2_train = R_sqr_function(Y, y_pred),
                         R2_val = R_sqr_function(Y_val, y_val_pred),
                         q = argvals_Y,
                         nComp = nComp,
                         beta.num = beta.num,
                         rep_num = rep_num,
                         method = "pFFR_I")
      )
      
      
      rm(y_pred, y_val_pred, beta_df, beta_hat, mean_imse_Y_val) # I'll reuse the same model name
      
      
      
      
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
      
      # Predict validation set (re-sourced predict function with fixes):
      y_val_pred <- predict_fregre_fr(object = m_final_rs,
                                      new.fdataobj = 
                                        fda.usc::fdata(X_val,
                                                       argvals = argvals_X,
                                                       rangeval = range(argvals_X)) )$data
      
      
      all_beta_hats <- rbind(all_beta_hats, beta_df )
      
      
      # Compute the average IMSE for the validation set:
      mean_imse_Y_val <- mean_imse_by_row(Y_val, # trueX
                                          y_val_pred,
                                          argvals_Y)
      
      
      all_final_res <- rbind(all_final_res,
                             tibble(
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
      
      
      
      # Compute R2 for the training and validation sets:
      y_pred <- fitted.values(m_final_rs)$data
      
      all_r2s <- rbind(all_r2s,
                       tibble(
                         R2_train = R_sqr_function(Y, y_pred),
                         R2_val = R_sqr_function(Y_val, y_val_pred),
                         q = argvals_Y,
                         nComp = nComp,
                         beta.num = beta.num,
                         rep_num = rep_num,
                         method = "pFFR_RS")
      )
      
      
      rm(m_final_rs, y_pred, y_val_pred, beta_df, beta_hat, mean_imse_Y_val) # I'll reuse the same model name
      
      
      
    }  # end nComp loop
    
    
    # Savings:
    
    saveRDS(all_CVEs, file =
              paste0(out_folder,
                     "cves_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    saveRDS(all_final_res, file =
              paste0(out_folder,
                     "final_models_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    saveRDS(all_beta_hats, file =
              paste0(out_folder,
                     "betas_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    saveRDS(best_lambdas, file =
              paste0(out_folder,
                     "best_lambdas_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    if (do_opt_bases_FFPLS) {
      saveRDS(best_num_bases, file =
                paste0(out_folder,
                       "best_num_bases_FFPLS_rep_",
                       rep_num,
                       "_beta_",
                       beta.num,
                       ".Rds"))
    }
    
    
    saveRDS(all_r2s, file =
              paste0(out_folder,
                     "R2s_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    saveRDS(computation_times, file =
              paste0(out_folder,
                     "computation_times_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    
  } # end nested inner rep_num
  
} # end nested outer beta_num



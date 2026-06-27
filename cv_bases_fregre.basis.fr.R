#' Cross-Validation for `fda.usc::fregre.basis.fr` with bases selection
#'
#' Performs cross-validation to select the optimal number of basis functions for functional
#' predictors and responses in a `fda.usc::fregre.basis.fr` model.
#'
#' @param X Matrix of predictors. Rows represent observations, and columns represent variables.
#' @param Y Matrix of responses. Rows represent observations, and columns represent variables.
#' @param argvals_X Vector of argument values corresponding to the functional predictors. Default is `NULL`.
#' @param argvals_Y Vector of argument values corresponding to the functional responses. Default is `NULL`.
#' @param num_bases_X Integer or vector specifying the number of basis functions for the predictors. Default is 20.
#' @param num_bases_Y Integer or vector specifying the number of basis functions for the responses. Default is 20.
#' @param fda_basis_func_X Function to create the basis functions for the predictors. Defaults to `fda::create.bspline.basis`.
#' @param fda_basis_func_Y Function to create the basis functions for the responses. Defaults to `fda::create.bspline.basis`.
#' @param penalty_X Numeric penalty for smoothing the predictor basis functions. Default is 0.
#' @param penalty_Y Numeric penalty for smoothing the response basis functions. Default is 0.
#' @param folds Integer or list specifying the cross-validation folds. If an integer, it indicates the number of folds. If a list, it provides predefined fold indices. Default is 5.
#' @param verbose Logical; if `TRUE`, progress messages are displayed during cross-validation. Default is `TRUE`.
#' @param stripped Logical; if `TRUE`, the final model is not included in the output, only cross-validation results. Default is `TRUE`.
#' @param Lfd_X A linear differential operator object for the penalization of `X`. Default is `NULL`.
#' @param Lfd_Y A linear differential operator object for the penalization of `X`. Default is `NULL`.
#' @param ... Additional arguments passed to the `fregre.basis.fr` function.
#'
#' @return A list containing:
#' \describe{
#'   \item{CVEs}{Cross-validation error for the best combination of bases.}
#'   \item{MSE_fold}{Vector of mean squared errors (MSE) for each fold.}
#'   \item{best_num_bases}{Vector indicating the best number of basis functions for predictors and responses.}
#'   \item{final_model}{The final model (only included if `stripped = FALSE`).}
#'   \item{elapsed}{Elapsed time for the cross-validation process.}
#' }
#'
#' @export
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#'
#' @examples
#' # 1D example:
cv_bases_fregre.basis.fr <- function(X,
                                     Y,
                                     argvals_X = NULL,
                                     argvals_Y = NULL,
                                     num_bases_X = 20,
                                     num_bases_Y = 20,
                                     fda_basis_func_X = fda::create.bspline.basis,
                                     fda_basis_func_Y = fda::create.bspline.basis,
                                     penalty_X = 0,
                                     penalty_Y = 0,
                                     folds = 5,
                                     verbose = TRUE,
                                     stripped = TRUE,
                                     Lfd_X = NULL,
                                     Lfd_Y = NULL,
                                     ...) {
  
  tictoc::tic("Crossvalidation")
  
  
  # Get range of argvalues:
  range_X <- range(argvals_X)
  range_Y <- range(argvals_Y)
  
  # Build differential operator if NULL:
  if (is.null(Lfd_X)) {
    warning("Lfd_X is NULL. Setting penalty for X to 0.")
    Lfd_X = fda::vec2Lfd(c(0, 0), range_X)
    penalty_X = 0
  }
  if (is.null(Lfd_Y)) {
    warning("Lfd_Y is NULL. Setting penalty for Y to 0.")
    Lfd_Y = fda::vec2Lfd(c(0, 0), range_Y)
    penalty_Y = 0
  }
  
  # Initialize grid for the first component
  num_bases_grid <- expand.grid(numbases_X = num_bases_X,
                                numbases_Y = num_bases_Y)
  
  if (is.numeric(folds)) {
    
    num_folds <- folds
    folds <- caret::createFolds(1:nrow(Y), k = num_folds)
    
  }else if (is.list(folds)) {
    
    num_folds <- length(folds)
    
  }
  
  # Initialize CVEs:
  CVEs <- array(data = NA, dim = 1) # averaged
  
  MSE_fold <- matrix(data = NA,
                     nrow = 1,
                     ncol = num_folds) # MSE per component per fold
  colnames(MSE_fold) <- paste0("fold_", 1:num_folds)
  
  # Best number of bases:
  best_num_bases <- matrix(data = NA,
                           nrow = 1,
                           ncol = 2)
  
  # Loop over combination of basis number:
  
  i <- row_num_bases <- NULL
  MSE_lambda_fold <- foreach::foreach (i = 1:num_folds,
                                       .packages = c("penFoFPLS"),
                                       .combine = "cbind") %:%
    foreach::foreach(row_num_bases = 1:nrow(num_bases_grid),
                     .packages = c("penFoFPLS"),
                     .combine = 'c' ) %dopar%
    {
      
      # MSE_lambda_fold <- matrix(NA, nrow = nrow(num_bases_grid), ncol = num_folds)
      # colnames(MSE_lambda_fold) <- paste0("fold_", 1:num_folds)
      # rownames(MSE_lambda_fold) <- paste0("bases_grid_comb_", 1:nrow(num_bases_grid))
      # 
      # for (i in 1:num_folds) {
      #   for (row_num_bases in 1:nrow(num_bases_grid)) {
      
      # build train
      Y_fold_train <- fda.usc::fdata(Y[-folds[[i]], , drop = F], argvals = argvals_Y, rangeval = range_Y)
      X_fold_train <- fda.usc::fdata(X[-folds[[i]], , drop = F], argvals = argvals_X, rangeval = range_X)
      
      # build test:
      Y_fold_test <- fda.usc::fdata(Y[folds[[i]], , drop = F], argvals = argvals_Y, rangeval = range_Y)
      X_fold_test <- fda.usc::fdata(X[folds[[i]], , drop = F], argvals = argvals_X, rangeval = range_X)
      
      basisobj_X <- fda_basis_func_X(rangeval = range_X,
                                     nbasis = num_bases_grid[row_num_bases, "numbases_X"])
      
      basisobj_Y <- fda_basis_func_Y(rangeval = range_Y,
                                     nbasis = num_bases_grid[row_num_bases, "numbases_Y"])
      
      # no penalties:
      res_fpls <- fda.usc::fregre.basis.fr(x = X_fold_train,
                                           y = Y_fold_train,
                                           basis.s = basisobj_X,
                                           basis.t = basisobj_Y,
                                           lambda.s = penalty_X,
                                           lambda.t = penalty_Y,
                                           Lfdobj.s = Lfd_X,
                                           Lfdobj.t = Lfd_Y,
                                           weights = NULL )
      
      
      
      MSE_lambda_fold[row_num_bases , i] <-
        num_int_1d(argvals = argvals_Y,
                   f_obs = colMeans(
                     ((Y_fold_test - stats::predict(object = res_fpls, newdata = X_fold_test)  )^2)$data
                   ) 
        )
      
      
      
      #   } # loop row_lambda
      # } # loop fold
      
    } # nested loop parallel
  
  
  # Averaged MSE_fold:
  CVEs_bases <- rowMeans(MSE_lambda_fold)
  
  # Best penalties per component:
  sel_num_bases <- which.min(CVEs_bases)
  best_num_bases <- as.numeric(num_bases_grid[sel_num_bases, ])
  
  # Save the folds-averaged CV error:
  CVEs <- CVEs_bases[sel_num_bases]
  
  # Save MSEs per fold, for the best lambda:
  MSE_fold <- MSE_lambda_fold[sel_num_bases, ]
  
  
  names(MSE_fold) <- paste0("fold_", 1:num_folds)
  names(best_num_bases) <- colnames(num_bases_grid)
  
  
  if (stripped) {
    ret <- list(
      CVEs = CVEs,
      MSE_fold = MSE_fold,
      best_num_bases = best_num_bases,
      elapsed = tictoc::toc(quiet = !verbose)
    )
  }else {
    
    if (verbose) {
      cat("Fitting final model\n")
    }
    
    
    basisobj_X_best <- fda_basis_func_X(rangeval = range(argvals_X),
                                        nbasis = best_num_bases["numbases_X"])
    
    basisobj_Y_best <- fda_basis_func_Y(rangeval = range(argvals_Y),
                                        nbasis = best_num_bases["numbases_Y"])
    
    final_model <- fda.usc::fregre.basis.fr(x = fda.usc::fdata(X, argvals = argvals_X, rangeval = range_X),
                                            y = fda.usc::fdata(Y, argvals = argvals_Y, rangeval = range_Y),
                                            basis.s = basisobj_X_best,
                                            basis.t = basisobj_Y_best,
                                            lambda.s = penalty_X,
                                            lambda.t = penalty_Y,
                                            Lfdobj.s = Lfd_X,
                                            Lfdobj.t = Lfd_Y,
                                            weights = NULL )
    
    ret <- list(
      CVEs = CVEs,
      MSE_fold = MSE_fold,
      best_num_bases = best_num_bases,
      final_model = final_model,
      elapsed = tictoc::toc(quiet = !verbose)
    )
    
  }
  
  class(ret) <- "cv_fregre.basis.fr"
  
  return(ret)
}

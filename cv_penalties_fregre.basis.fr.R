#' Cross-Validation for `fda.usc::fregre.basis.fr` with penalties selection
#'
#' Performs cross-validation to select the optimal penalties for functional
#' predictors and responses in a `fda.usc::fregre.basis.fr` model.
#'
#' @param X Matrix of predictors. Rows represent observations, and columns represent variables.
#' @param Y Matrix of responses. Rows represent observations, and columns represent variables.
#' @param argvals_X Vector of argument values corresponding to the functional predictors. Default is `NULL`.
#' @param argvals_Y Vector of argument values corresponding to the functional responses. Default is `NULL`.
#' @param num_bases_X Integer specifying the number of basis functions for the predictors. Default is 20.
#' @param num_bases_Y Integer specifying the number of basis functions for the responses. Default is 20.
#' @param fda_basis_func_X Function to create the basis functions for the predictors. Defaults to `fda::create.bspline.basis`.
#' @param fda_basis_func_Y Function to create the basis functions for the responses. Defaults to `fda::create.bspline.basis`.
#' @param penalty_X Numeric vector of penalties for smoothing the predictor basis functions. Default is 0.
#' @param penalty_Y Numeric vector of penalties for smoothing the response basis functions. Default is 0.
#' @param folds Integer or list specifying the cross-validation folds. If an integer, it indicates the number of folds. If a list, it provides predefined fold indices. Default is 5.
#' @param verbose Logical; if `TRUE`, progress messages are displayed during cross-validation. Default is `TRUE`.
#' @param stripped Logical; if `TRUE`, the final model is not included in the output, only cross-validation results. Default is `TRUE`.
#' @param Lfd_X A linear differential operator object for the penalization of `X`. Default is `NULL`.
#' @param Lfd_Y A linear differential operator object for the penalization of `X`. Default is `NULL`.
#' @param ... Additional arguments passed to the `fregre.basis.fr` function.
#'
#' @return A list containing:
#' \describe{
#'   \item{CVEs}{Cross-validation error for the best combination of penalties}
#'   \item{MSE_fold}{Vector of mean squared errors (MSE) for each fold.}
#'   \item{best_penalties}{Vector indicating the best penalties for predictors and responses.}
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
cv_penalties_fregre.basis.fr <- function(X,
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
    warning("Lfd_X is NULL. Setting to fda::vec2Lfd(c(0, 0), range(argvals_X)).")
    Lfd_X = fda::vec2Lfd(c(0, 0), range_X)
  }
  if (is.null(Lfd_Y)) {
    warning("Lfd_Y is NULL. Setting to fda::vec2Lfd(c(0, 0), range(argvals_Y)).")
    Lfd_Y = fda::vec2Lfd(c(0, 0), range_Y)
  }
  
  # Build the basis functions:
  basisobj_X <- fda_basis_func_X(rangeval = range_X,
                                 nbasis = num_bases_X)
  
  basisobj_Y <- fda_basis_func_Y(rangeval = range_Y,
                                 nbasis = num_bases_Y)
  
  
  # Initialize grid for the first component
  penalty_grid <- expand.grid(penalty_X = penalty_X,
                              penalty_Y = penalty_Y)
  
  
  # Create folds if input is a just the number of folds:
  if (is.numeric(folds)) {
    
    num_folds <- folds
    folds <- caret::createFolds(1:nrow(Y), k = num_folds)
    
  }else if (is.list(folds)) {
    
    num_folds <- length(folds)
    
  }
  
  
  # Initialize CVEs and MSE_fold:
  CVEs <- array(data = NA, dim = 1) # averaged
  
  MSE_fold <- matrix(data = NA,
                     nrow = 1,
                     ncol = num_folds) # MSE per component per fold
  colnames(MSE_fold) <- paste0("fold_", 1:num_folds)
  
  # Best penalties:
  best_penalties <- matrix(data = NA,
                           nrow = 1,
                           ncol = 2)
  
  # Loop over combination of penalties:
  
  i <- row_lambda <- NULL
  MSE_lambda_fold <- foreach::foreach (i = 1:num_folds,
                                       .packages = c("penFoFPLS", "fda.usc"),
                                       .combine = "cbind") %:%
    foreach::foreach(row_lambda = 1:nrow(penalty_grid),
                     .packages = c("penFoFPLS", "fda.usc"),
                     .combine = 'c' ) %dopar%
    {
      
      # MSE_lambda_fold <- matrix(NA, nrow = nrow(penalty_grid), ncol = num_folds)
      # colnames(MSE_lambda_fold) <- paste0("fold_", 1:num_folds)
      # rownames(MSE_lambda_fold) <- paste0("penalty_grid_comb_", 1:nrow(penalty_grid))
      # 
      # for (i in 1:num_folds) {
      #   for (row_lambda in 1:nrow(penalty_grid)) {
      
      # build train
      Y_fold_train <- fda.usc::fdata(Y[-folds[[i]], , drop = F], argvals = argvals_Y, rangeval = range_Y)
      X_fold_train <- fda.usc::fdata(X[-folds[[i]], , drop = F], argvals = argvals_X, rangeval = range_X)
      
      # build test:
      Y_fold_test <- fda.usc::fdata(Y[folds[[i]], , drop = F], argvals = argvals_Y, rangeval = range_Y)
      X_fold_test <- fda.usc::fdata(X[folds[[i]], , drop = F], argvals = argvals_X, rangeval = range_X)
      
      
      # no penalties:
      res_fpls <- fda.usc::fregre.basis.fr(x = X_fold_train,
                                           y = Y_fold_train,
                                           basis.s = basisobj_X,
                                           basis.t = basisobj_Y,
                                           lambda.s = penalty_grid[row_lambda, "penalty_X"],
                                           lambda.t = penalty_grid[row_lambda, "penalty_Y"],
                                           Lfdobj.s = Lfd_X,
                                           Lfdobj.t = Lfd_Y,
                                           weights = NULL )
      
      # Fix class bug in fda.usc package:
      class(res_fpls) = "fregre.basis.fr"
      res_fpls$call[[1]]=="fregre.basis.fr"
      
      
      # MSE_lambda_fold[row_lambda , i] <-
      num_int_1d(argvals = argvals_Y,
                 f_obs = colMeans(
                   ((Y_fold_test - predict_fregre_fr(object = res_fpls, new.fdataobj = X_fold_test)  )^2)$data
                 ) 
      )
      
      
      
      #   } # loop row_lambda
      # } # loop fold
      
    } # nested loop parallel
  
  
  # Averaged MSE_fold:
  CVEs_penalties <- rowMeans(MSE_lambda_fold)
  
  # Best penalties per component:
  sel_penalties <- which.min(CVEs_penalties)
  best_penalties <- as.numeric(penalty_grid[sel_penalties, ])
  
  # Save the folds-averaged CV error:
  CVEs <- CVEs_penalties[sel_penalties]
  
  # Save MSEs per fold, for the best lambda:
  MSE_fold <- MSE_lambda_fold[sel_penalties, ]
  
  
  names(MSE_fold) <- paste0("fold_", 1:num_folds)
  names(best_penalties) <- colnames(penalty_grid)
  
  
  if (stripped) {
    ret <- list(
      CVEs = CVEs,
      MSE_fold = MSE_fold,
      best_penalties = best_penalties,
      elapsed = tictoc::toc(quiet = !verbose)
    )
  }else {
    
    if (verbose) {
      cat("Fitting final model\n")
    }
    
    
    
    final_model <- fda.usc::fregre.basis.fr(x = fda.usc::fdata(X, argvals = argvals_X, rangeval = range_X),
                                            y = fda.usc::fdata(Y, argvals = argvals_Y, rangeval = range_Y),
                                            basis.s = basisobj_X,
                                            basis.t = basisobj_Y,
                                            lambda.s = best_penalties["penalty_X"],
                                            lambda.t = best_penalties["penalty_Y"],
                                            Lfdobj.s = Lfd_X,
                                            Lfdobj.t = Lfd_Y,
                                            weights = NULL )
    
    # Fix class bug in fda.usc package:
    class(final_model) = "fregre.basis.fr"
    final_model$call[[1]]=="fregre.basis.fr"
    
    ret <- list(
      CVEs = CVEs,
      MSE_fold = MSE_fold,
      best_penalties = best_penalties,
      final_model = final_model,
      elapsed = tictoc::toc(quiet = !verbose)
    )
    
  }
  
  class(ret) <- "cv_fregre.basis.fr"
  
  return(ret)
}

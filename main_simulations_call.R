# Please install the R package "penFoFPLS" from GitHub:
# devtools::install_github("hhroig/penFoFPLS", dependencies = TRUE)

library(pls)
library(dplyr)
library(penFoFPLS)
library(fda)
library(refund)
library(reshape2)
library(doParallel)


main_simulations_call <- function(
    
  do_setting = 3, # settings 1, 2, or 3
  
  X_sd_error = 0, # Extra observation error for X (after data generation)
  
  center = TRUE, 
  
  num_betas = c("cos_sin", "cos_sum"),  # betas ids
  
  num_lambdas = 10, # number of lambdas to be used in the grid search (for all)
  
  lower_penalty_bound_RS = -6, # lower bound for the penalty grid search in Ramsay & Silverman
  upper_penalty_bound_RS = 8, # upper bound for the penalty grid search in Ramsay & Silverman
  
  lower_penalty_bound = -2, # lower bound for the penalty grid search in the proposed method
  upper_penalty_bound = 12, # upper bound for the penalty grid search in the proposed method
  
  nnodesX = 100, # number of nodes for the functional predictors
  nnodesY = 100, # number of nodes for the functional response
  
  total_reps  = 100, # number of repetitions
  rep_starts = 1, # starting number of the repetitions
  
  num_folds = 5 # number of folds for the cross-validation
  
) {
  
  # Grid search penalties for Ramsay & Silverman:
  lambdas_in_RS  <-  seq(lower_penalty_bound_RS, 
                         upper_penalty_bound_RS, 
                         length.out = num_lambdas)
  lambdas_in_RS  <-  10^(lambdas_in_RS)  
  penaltyvec_X_RS <- penaltyvec_Y_RS <- lambdas_in_RS
  
  
  # Grid search penalties for the proposed method:
  lambdas_in  <-  seq(lower_penalty_bound, upper_penalty_bound, length.out = num_lambdas)
  lambdas_in  <-  10^(lambdas_in)
  penaltyvec_X <- penaltyvec_Y <- lambdas_in
  
  
  # Observation nodes (for data generation):
  argvals_X <- seq(0, 1, length.out = nnodesX) # p
  argvals_Y <- seq(0, 1, length.out = nnodesY) # q
  
  
  # Number of basis depending on the setting:
  
  # # Setting 1:
  if (do_setting == 1) {
    LL <- 7 # number of basis for Y(q)
    KK <- 7 # number of basis for X(p)
    do_opt_bases_FFPLS = FALSE
    
    # number of PLS components to compute:
    max_nComp <- 6
  }
  
  
  # # Setting 2: 
  if (do_setting == 2) {
    LL <- 40 # number of basis for Y(q)
    KK <- 40 # number of basis for X(p)
    do_opt_bases_FFPLS = FALSE
    
    # number of PLS components to compute:
    max_nComp <- 8
  }
  
  
  # # Setting 3: 
  # # select number of basis using CVE for non-penalized method
  # # use the following for the penalized approach:
  if (do_setting == 3) {
    LL <- 40 # number of basis for Y(q)
    KK <- 40 # number of basis for X(p)
    
    # number of PLS components to compute:
    max_nComp <- 8
    
    LL_list <- round(seq(9, 40, length.out = num_lambdas)) # list of number of bases for Y(q)
    KK_list <- LL_list                                     # list of number of bases for X(p)
    do_opt_bases_FFPLS = TRUE
  }
  
  
  # Number of basis for Ramsay and Silverman
  # KK_rs = LL_rs = 5
  KK_rs = KK
  LL_rs = LL
  
  
  # B-spline basis:
  basisobj_X <- fda::create.bspline.basis(rangeval = range(argvals_X),
                                          nbasis = KK)
  basisobj_Y <- fda::create.bspline.basis(rangeval = range(argvals_Y),
                                          nbasis = LL)
  
  # Output folder:
  out_folder <- paste0("results_simulations/",
                       "set", do_setting,
                       ifelse(X_sd_error > 0, "e", ""),
                       "_rep", 
                       total_reps, 
                       "_pen", 
                       length(penaltyvec_X)*length(penaltyvec_Y),
                       "_K", KK, "L", LL,
                       "/")
  
  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }
  
  
  
  # Call the actual simulations ----------------------------------------------
  
  
  # Source the wrappers for R&S' method coded in fda.usc 
  source("cv_penalties_fregre.basis.fr.R")
  source("predict_fregre_fr.R")
  
  nodes_CL = detectCores()   # Detect number of cores to use
  cl = makeCluster(nodes_CL) # Specify number of threads here
  
  clusterExport(cl, c("predict_fregre_fr"))
  
  registerDoParallel(cl)
  
  source("simulations_fofr_v2_with_ivanescus_ramsay_silverman.R", local = TRUE)
  
  stopCluster(cl)
  
  
  # Plot comparisons --------------------------------------------------------
  
  
  source("compare_methods_fofr_with_ivanescu_ramsay_silverman.R", local = TRUE)
  
  compare_methods_fun(input_folder = out_folder)
  
  
}



# Run the simulations -----------------------------------------------------


global_num_lambdas = 3
global_total_reps = 3
global_start_reps = 1
global_betas = c("cos_sin", "sin_sum", "cos_sum")


# Setting 1:
main_simulations_call(
  do_setting = 1, 
  X_sd_error = 0, 
  
  num_betas = global_betas,  # betas ids
  
  num_lambdas = global_num_lambdas, 
  total_reps  = global_total_reps,
  rep_starts = global_start_reps
)

main_simulations_call(
  do_setting = 1, 
  X_sd_error = 0.2, 
  
  num_betas = global_betas,  # betas ids
  
  num_lambdas = global_num_lambdas, 
  total_reps  = global_total_reps,
  rep_starts = global_start_reps
)


# 
# # Setting 2:
# main_simulations_call(
#   do_setting = 2, 
#   X_sd_error = 0, 
#   
#   num_betas = global_betas,  # betas ids
#   
#   num_lambdas = global_num_lambdas, 
#   total_reps  = global_total_reps,
#   rep_starts = global_start_reps
# )
# 
# main_simulations_call(
#   do_setting = 2, 
#   X_sd_error = 0.2, 
#   
#   num_betas = global_betas,  # betas ids
#   
#   num_lambdas = global_num_lambdas, 
#   total_reps  = global_total_reps,
#   rep_starts = global_start_reps
# )



# Setting 3:
main_simulations_call(
  do_setting = 3, 
  X_sd_error = 0, 
  
  num_betas = global_betas,  # betas ids
  
  num_lambdas = global_num_lambdas, 
  total_reps  = global_total_reps,
  rep_starts = global_start_reps
)

main_simulations_call(
  do_setting = 3, 
  X_sd_error = 0.2, 
  
  num_betas = global_betas,  # betas ids
  
  num_lambdas = global_num_lambdas, 
  total_reps  = global_total_reps,
  rep_starts = global_start_reps
)
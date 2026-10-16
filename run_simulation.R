# Buried-bands simulation study: configuration and entry point.
#
# Fits all five methods -- pFFPLS (proposed), FFPLS, FFPLS_OB, pFFR_I (refund::pffr) and
# pFFR_RS (fda.usc::fregre.basis.fr) -- on data from simulations_buried_bands.R and saves one
# set of .Rds files per repetition to
#   <results_dir>/rep<total_reps>_K<KK>L<LL>_sdX<X_sd_error>_band<band_sd>_seed<seed_base>/
# The defaults reproduce the reference study:
#   results_simulations/buried_bands/rep100_K40L40_sdX1.5_band0.035_seed2026/
#
# Usage, from the project root:
#   Rscript run_simulation.R [rep_start rep_end] [n_cores] [band_sd]
#   e.g. the full study:  Rscript run_simulation.R 1 100 8
# Repetition r always uses set.seed(seed_base + r), so the repetitions can be split across
# several concurrent processes writing to the same folder, e.g. four processes:
#   Rscript run_simulation.R  1  25 2  &  Rscript run_simulation.R 26  50 2  &
#   Rscript run_simulation.R 51  75 2  &  Rscript run_simulation.R 76 100 2
# pFFR_I is single-threaded and dominates the cost (~35 min/rep at K = L = 40, ~2 GB RAM
# each), so several processes with few cores each is much faster than one process.
#
# Requires the penFoFPLS package:
#   devtools::install_github("hhroig/penFoFPLS", dependencies = TRUE)

library(pls)
library(dplyr)
library(penFoFPLS)
library(fda)
library(refund)
library(reshape2)
library(doParallel)


main_buried_bands_call <- function(
    num_betas   = c("buried_bands"),
    total_reps  = 100,
    rep_starts  = 1,
    rep_end     = total_reps, # last repetition run by THIS process (see header)

    # Data-generating parameters:
    band_sd     = 0.035, # width of the two Gaussian bands in p
    band_share  = 0.03,  # share of X's total integrated variance carried by the bands
    X_sd_error  = 1.5,   # sd of the white measurement noise added to X (sigma_p)
    Rsq         = 0.9,   # target R^2(q) for the response noise

    # Estimation: K = L = 40 cubic B-splines for X and Y (all methods except FFPLS_OB).
    KK = 40, LL = 40,
    max_nComp = 10,
    num_lambdas = 10,    # number of candidate basis sizes per axis for FFPLS_OB

    # pFFR_RS penalty grid (lambda_X x lambda_Y):
    penaltyvec_X_RS = 10^seq(-10, -2, length.out = 5),
    penaltyvec_Y_RS = c(1e-2, 1e-1, 1e0, 1e1, 1e4),

    # pFFPLS penalty grid: log10 bounds and number of points per axis:
    lower_penalty_bound_X = -6, upper_penalty_bound_X = 4, num_lambdas_X = 8,
    lower_penalty_bound_Y = -6, upper_penalty_bound_Y = 4, num_lambdas_Y = 8,

    nnodesX = 100,
    nnodesY = 100,
    num_folds = 5,
    center = TRUE,
    n_cores = parallel::detectCores(),
    do_pffr = TRUE,      # FALSE skips pFFR_I (output folder tagged "_noPFFR")
    seed_base = 2026,    # repetition r uses set.seed(seed_base + r)
    results_dir = "results_simulations/buried_bands/"
) {

  KK_rs <- LL_rs <- KK
  do_opt_bases_FFPLS <- TRUE

  argvals_X <- seq(0, 1, length.out = nnodesX) # p
  argvals_Y <- seq(0, 1, length.out = nnodesY) # q

  basisobj_X <- fda::create.bspline.basis(rangeval = range(argvals_X), nbasis = KK)
  basisobj_Y <- fda::create.bspline.basis(rangeval = range(argvals_Y), nbasis = LL)

  # FFPLS_OB basis-size search; never below max_nComp (a PLS fit cannot extract more
  # components than its basis dimension):
  LL_list <- round(seq(max(10, max_nComp), LL, length.out = num_lambdas))
  KK_list <- LL_list

  penaltyvec_X <- 10^seq(lower_penalty_bound_X, upper_penalty_bound_X, length.out = num_lambdas_X)
  penaltyvec_Y <- 10^seq(lower_penalty_bound_Y, upper_penalty_bound_Y, length.out = num_lambdas_Y)

  # Tag the folder whenever pFFPLS's penalty grid departs from the default, so runs with
  # different grids never overwrite each other:
  pgrid_tag <- if (!identical(c(lower_penalty_bound_X, upper_penalty_bound_X, num_lambdas_X,
                                lower_penalty_bound_Y, upper_penalty_bound_Y, num_lambdas_Y),
                              c(-6, 4, 8, -6, 4, 8))) {
    paste0("_pgX", lower_penalty_bound_X, "_", upper_penalty_bound_X, "n", num_lambdas_X,
           "_Y", lower_penalty_bound_Y, "_", upper_penalty_bound_Y, "n", num_lambdas_Y)
  } else ""

  out_folder <- paste0(results_dir,
                       "rep", total_reps,
                       "_K", KK, "L", LL,
                       "_sdX", X_sd_error,
                       "_band", band_sd,
                       if (!is.null(seed_base)) paste0("_seed", seed_base) else "",
                       if (Rsq != 0.9) paste0("_Rsq", Rsq) else "",
                       if (!do_pffr) "_noPFFR" else "",
                       pgrid_tag,
                       "/")

  if (!dir.exists(out_folder)) {
    dir.create(out_folder, recursive = TRUE)
  }

  source("cv_penalties_fregre.basis.fr.R")
  source("predict_fregre_fr.R")

  cl <- makeCluster(n_cores)
  clusterExport(cl, c("predict_fregre_fr"))
  registerDoParallel(cl)

  source("simulations_buried_bands.R", local = TRUE)

  stopCluster(cl)

  invisible(out_folder)
}


# Run from the command line (does nothing when this file is source()d):
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  main_buried_bands_call(
    rep_starts = if (length(args) >= 1) as.integer(args[1]) else 1,
    rep_end    = if (length(args) >= 2) as.integer(args[2]) else 100,
    n_cores    = if (length(args) >= 3) as.integer(args[3]) else parallel::detectCores(),
    band_sd    = if (length(args) >= 4) as.numeric(args[4]) else 0.035
  )
}

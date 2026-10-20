# Buried-bands simulation study: configuration and entry point.
#
# Fits pFFPLS (proposed), FFPLS, FFPLS_OB, pFFR_I (refund::pffr) and pFFR_RS
# (fda.usc::fregre.basis.fr) on data from simulations_buried_bands.R and saves one set of
# .Rds files per repetition to results_simulations/buried_bands/<folder>/.
#
# The scenarios of the paper are predefined in `scenarios` below (100 repetitions each):
#   main              K = L = 40, all five methods, sigma_p = 1.5
#                     -> rep100_K40L40_sdX1.5_band0.035_seed2026/
#   main_noiseless    as main, without measurement noise on X (sigma_p = 0)
#                     -> rep100_K40L40_sdX0_band0.035_seed2026/
#   lowdim            K = L = 7, without FFPLS_OB, up to 7 components, wider penalty grids
#                     -> rep100_K7L7_sdX1.5_band0.035_seed2026_noOB/
#   lowdim_noiseless  as lowdim, sigma_p = 0
#                     -> rep100_K7L7_sdX0_band0.035_seed2026_noOB/
#   band0.02          as main, with narrower bands (band_sd = 0.02)
#                     -> rep100_K40L40_sdX1.5_band0.02_seed2026/
# Every scenario uses the same seeds, so all scenarios share the same clean curves and
# responses (the noiseless scenarios differ only in the noise added to X).
#
# Usage, from the project root:
#   Rscript run_simulation.R <scenario> [rep_start rep_end] [n_cores]
#   e.g.  Rscript run_simulation.R main 1 100 8
# Repetition r always uses set.seed(2026 + r), so the repetitions can be split across
# several concurrent processes writing to the same folder, e.g. four processes:
#   Rscript run_simulation.R main  1  25 2  &  Rscript run_simulation.R main 26  50 2  &
#   Rscript run_simulation.R main 51  75 2  &  Rscript run_simulation.R main 76 100 2
# At K = L = 40, pFFR_I is single-threaded and dominates the cost (~35 min/rep, ~2 GB RAM
# per process), so several processes with few cores each is much faster than one process.
# From R: source("run_simulation.R"); run_scenario("lowdim", rep_starts = 1, rep_end = 10).
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
    do_opt_bases = TRUE, # FALSE skips FFPLS_OB (output folder tagged "_noOB"), e.g. for small K = L
    seed_base = 2026,    # repetition r uses set.seed(seed_base + r)
    results_dir = "results_simulations/buried_bands/",
    folder_name = NULL,  # output folder name; NULL builds it from the settings
    dry_run = FALSE      # TRUE: only return the output folder, run nothing
) {

  KK_rs <- LL_rs <- KK
  do_opt_bases_FFPLS <- do_opt_bases

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

  auto_name <- paste0("rep", total_reps,
                       "_K", KK, "L", LL,
                       "_sdX", X_sd_error,
                       "_band", band_sd,
                       if (!is.null(seed_base)) paste0("_seed", seed_base) else "",
                       if (Rsq != 0.9) paste0("_Rsq", Rsq) else "",
                       if (!do_pffr) "_noPFFR" else "",
                       if (!do_opt_bases) "_noOB" else "",
                       pgrid_tag)
  out_folder <- paste0(results_dir, if (is.null(folder_name)) auto_name else folder_name, "/")

  # Windows cannot open paths longer than 260 characters; the longest file written later
  # (a results_plots/ figure) adds about 80 characters to the folder path.
  if (nchar(normalizePath(out_folder, mustWork = FALSE)) + 80 > 259) {
    warning("Output folder path is long (", nchar(normalizePath(out_folder, mustWork = FALSE)),
            " characters); saving results or figures may fail on Windows. Use a shorter `folder_name`.")
  }

  if (dry_run) return(out_folder)

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


# Scenarios of the paper ------------------------------------------------------------
# Each entry lists only the arguments that differ from the defaults of
# main_buried_bands_call(); total_reps = 100 and seed_base = 2026 for all of them.

lowdim_settings <- list(
  KK = 7, LL = 7, max_nComp = 7, do_opt_bases = FALSE,
  # wider grids: at K = L = 7 the default grids put the CV selections on their edges
  penaltyvec_X_RS = 10^seq(-12, 4, by = 2), penaltyvec_Y_RS = 10^seq(-14, 2, by = 2),
  lower_penalty_bound_X = -12, upper_penalty_bound_X = 4, num_lambdas_X = 17,
  lower_penalty_bound_Y = -8,  upper_penalty_bound_Y = 4, num_lambdas_Y = 13
)

scenarios <- list(
  main             = list(),
  main_noiseless   = list(X_sd_error = 0),
  lowdim           = c(lowdim_settings, list(folder_name = "rep100_K7L7_sdX1.5_band0.035_seed2026_noOB")),
  lowdim_noiseless = c(lowdim_settings, list(X_sd_error = 0,
                                             folder_name = "rep100_K7L7_sdX0_band0.035_seed2026_noOB")),
  band0.02         = list(band_sd = 0.02)
)

# Run one scenario; `...` overrides its settings (e.g. rep_starts, rep_end, n_cores).
run_scenario <- function(name, ...) {
  if (!name %in% names(scenarios)) {
    stop("Unknown scenario '", name, "'. Available: ", paste(names(scenarios), collapse = ", "))
  }
  args <- modifyList(scenarios[[name]], list(...))
  do.call(main_buried_bands_call, args)
}


# Run from the command line (does nothing when this file is source()d):
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: Rscript run_simulation.R <scenario> [rep_start rep_end] [n_cores]")
  run_scenario(
    args[1],
    rep_starts = if (length(args) >= 2) as.integer(args[2]) else 1,
    rep_end    = if (length(args) >= 3) as.integer(args[3]) else 100,
    n_cores    = if (length(args) >= 4) as.integer(args[4]) else parallel::detectCores()
  )
}

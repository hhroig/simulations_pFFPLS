# Builds the results_plots/ folder for one results folder written by run_simulation.R:
#   - compare_methods_fun(): IMSE / CVE / test IMSE(Y) boxplots and curves, fR2, computing
#     times, selected penalties, Excel summaries, and 3D (persp) mean surfaces;
#   - plot_beta_heatmaps(): 2D heatmaps of the true and mean estimated surfaces.
#
# Usage, from the project root:
#   Rscript make_plots.R [results_folder]
# Default: the reference study,
#   results_simulations/buried_bands/rep100_K40L40_sdX1.5_band0.035_seed2026/

source("compare_methods_fofr_with_ivanescu_ramsay_silverman.R", local = TRUE)
source("plot_beta_heatmaps.R", local = TRUE)

args <- commandArgs(trailingOnly = TRUE)
input_folder <- if (length(args) >= 1) args[1] else
  "results_simulations/buried_bands/rep100_K40L40_sdX1.5_band0.035_seed2026/"
if (!grepl("/$", input_folder)) input_folder <- paste0(input_folder, "/")

compare_methods_fun(input_folder = input_folder)
plot_beta_heatmaps(input_folder)

cat("\n=== DONE: results_plots/ written under", input_folder, "===\n")

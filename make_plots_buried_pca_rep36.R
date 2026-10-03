# Generates the standard results_plots/ folder (IMSE/CVE boxplots, fR2
# plots, computation-time plots, penalty plots, 3D beta surfaces, Excel
# summaries) for the final buried-PCA-signal 36-rep run, using the existing
# compare_methods_fun() unchanged -- does not modify that file or any other
# existing script.

source("compare_methods_fofr_with_ivanescu_ramsay_silverman.R", local = TRUE)

compare_methods_fun(
  input_folder = "results_simulations/buried_pca_signal/rep36_K15L15_finegrid/"
)

cat("\n=== DONE: results_plots/ written under rep36_K15L15_finegrid/ ===\n")

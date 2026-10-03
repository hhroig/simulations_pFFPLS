# Generates the standard results_plots/ folder (IMSE/CVE boxplots, fR2
# plots, computation-time plots, penalty plots, 3D beta surfaces, Excel
# summaries) for the final buried-PCA-SHARP-signal 30-rep run (Rsq=0.5,
# K=L=25), using the existing compare_methods_fun() unchanged -- does not
# modify that file or any other existing script.

source("compare_methods_fofr_with_ivanescu_ramsay_silverman.R", local = TRUE)

compare_methods_fun(
  input_folder = "results_simulations/buried_pca_sharp_signal/rep30_K25L25_Rsq0.5/"
)

cat("\n=== DONE: results_plots/ written under rep30_K25L25_Rsq0.5/ ===\n")

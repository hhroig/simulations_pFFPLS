# 3D (persp) views of the coefficient surfaces estimated in single repetitions, one
# folder per repetition. Complements the 3D views of the mean surfaces made by
# compare_methods_fun() (results_plots/3DBeta2D_<beta>/): a mean over repetitions
# smooths out rep-to-rep variability, so these show what one fit actually looks like.
#
# Output, under <input_folder>/results_plots/3DBeta2D_<beta>_rep<r>/:
#   shared_scale/<method>_beta<beta>_ncomp<H>.{pdf,eps,png}
#                -- one z-axis range for every method, every H and the true beta of
#                   repetition r (compare magnitudes directly)
#   owned_scale/<method>_beta<beta>_ncomp<H>.{pdf,eps,png}
#                -- each surface with its own z-axis range (see each surface's shape)
#   True_beta<beta>_ncomp<H>.* in both subfolders, repeated for every H so that files
#   pair up by name. pFFR_I and pFFR_RS do not depend on H; they are repeated too.
# File names, sizes, viewing angles and styling follow 3DBeta2D_<beta>/.
#
# Usage, from the project root:
#   Rscript plot_beta_reps_3D.R [results_folder] [reps]
# reps is an R vector without spaces, e.g. 1:5 (default) or 1,7,12. Seeds are paired
# across scenarios (2026 + r), so rep r shows the same data set in every scenario.

plot_beta_reps_3D <- function(input_folder, reps = 1:5, n_comps = NULL,
                              theta = 40, phi = 25, formats = c("pdf", "eps", "png")) {

  if (!grepl("/$", input_folder)) input_folder <- paste0(input_folder, "/")

  # Surface as a matrix with rows = p and columns = q, as persp(x = p, y = q, z) expects:
  surface_matrix <- function(p, q, z) {
    o  <- order(p, q)
    pp <- sort(unique(p))
    qq <- sort(unique(q))
    stopifnot(length(z) == length(pp) * length(qq))
    list(p = pp, q = qq, z = matrix(z[o], nrow = length(pp), ncol = length(qq), byrow = TRUE))
  }

  save_persp <- function(file_stem, s, zlim) {
    draw <- function() {
      persp(x = s$p, y = s$q, z = s$z, col = "white",
            xlab = "p", ylab = "q", zlab = "z", zlim = zlim,
            theta = theta, phi = phi, expand = 0.5, shade = 0.5, ticktype = "detailed")
    }
    if ("pdf" %in% formats) {
      pdf(paste0(file_stem, ".pdf"), width = 7, height = 5); draw(); dev.off()
    }
    if ("eps" %in% formats) {
      postscript(paste0(file_stem, ".eps"), width = 7, height = 5,
                 horizontal = FALSE, paper = "special")
      draw(); dev.off()
    }
    if ("png" %in% formats) {
      png(paste0(file_stem, ".png"), width = 800, height = 600, res = 100); draw(); dev.off()
    }
  }

  for (r in reps) {

    file <- list.files(input_folder, pattern = paste0("^betas_rep_", r, "_.*\\.Rds$"),
                       full.names = TRUE)
    if (length(file) != 1) {
      warning("Skipping rep ", r, ": expected one betas_rep_", r, "_*.Rds file in ",
              input_folder, ", found ", length(file), ".")
      next
    }
    b <- readRDS(file)
    H_all <- if (is.null(n_comps)) sort(unique(b$nComp)) else n_comps

    for (beta_num in unique(b$beta.num)) {

      bb <- b[b$beta.num == beta_num, ]

      out_dir     <- paste0(input_folder, "results_plots/3DBeta2D_", beta_num, "_rep", r, "/")
      path_shared <- paste0(out_dir, "shared_scale/")
      path_owned  <- paste0(out_dir, "owned_scale/")
      dir.create(path_shared, recursive = TRUE, showWarnings = FALSE)
      dir.create(path_owned,  recursive = TRUE, showWarnings = FALSE)

      # True beta, on the full-resolution grid (taken from the pFFPLS rows):
      tr <- bb[bb$method == "pFFPLS" & bb$nComp == min(bb$nComp), ]
      s_true <- surface_matrix(tr$p, tr$q, tr$z_true)

      # Shared scale: every method and every H of this repetition, plus the truth.
      zs_scale <- range(bb$z, s_true$z, na.rm = TRUE)

      for (H in H_all) {
        for (m in unique(bb$method)) {
          d <- bb[bb$method == m & bb$nComp == H, ]
          if (nrow(d) == 0) next
          s <- surface_matrix(d$p, d$q, d$z)
          stem <- paste0(m, "_beta", beta_num, "_ncomp", H)
          save_persp(paste0(path_shared, stem), s, zs_scale)
          save_persp(paste0(path_owned, stem), s, range(s$z, na.rm = TRUE))
        }
        stem <- paste0("True_beta", beta_num, "_ncomp", H)
        save_persp(paste0(path_shared, stem), s_true, zs_scale)
        save_persp(paste0(path_owned, stem), s_true, range(s_true$z, na.rm = TRUE))
      }

      cat("Saved rep", r, "surfaces to", out_dir, " shared zlim =",
          paste(signif(zs_scale, 3), collapse = " to "), "\n")
    }
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  input_folder <- if (length(args) >= 1) args[1] else
    "results_simulations/buried_bands/rep100_K40L40_sdX1.5_band0.035_seed2026/"
  reps <- if (length(args) >= 2) as.integer(eval(parse(text = paste0("c(", args[2], ")")))) else 1:5
  plot_beta_reps_3D(input_folder, reps = reps)
}

# 2D heatmaps of the true coefficient surface beta(q, p) and of each method's mean
# estimated surface (averaged over repetitions), one set of figures per number of
# components H. Complements the 3D (persp) surfaces made by compare_methods_fun().
#
# Output, under <input_folder>/results_plots/2DBetaHeatmaps_<beta>/:
#   all_methods_ncomp<H>.{png,pdf}  -- true beta + every method in one panel, one shared
#                                      colour scale (compare magnitudes directly)
#   owned_scale/<method>_ncomp<H>.{png,pdf} -- one surface per file, its own colour scale
#                                      (see each surface's shape, e.g. pFFR_RS's flat estimate)
# pFFR_I and pFFR_RS do not depend on H; they are repeated in every panel.
#
# Colour: diverging blue (negative) - light gray (0) - red (positive), with limits
# symmetric around 0 so equal colours mean equal magnitudes on both sides.

library(ggplot2)

plot_beta_heatmaps <- function(input_folder, n_comps = NULL, width = 10, height = 6.5) {

  files <- list.files(input_folder, pattern = "^betas_rep_.*\\.Rds$", full.names = TRUE)
  if (length(files) == 0) stop("No betas_rep_*.Rds files in ", input_folder)

  # Mean surface over repetitions. Every file has the same row layout
  # (method, nComp, q, p), so the estimates can be summed row by row; fall back to a
  # keyed merge if a file ever differs.
  first <- readRDS(files[1])
  keys  <- first[, c("method", "nComp", "beta.num", "q", "p")]
  z_sum <- first$z
  for (f in files[-1]) {
    b <- readRDS(f)
    if (identical(b[, c("method", "nComp", "q", "p")], keys[, c("method", "nComp", "q", "p")])) {
      z_sum <- z_sum + b$z
    } else {
      m <- merge(cbind(keys, row = seq_len(nrow(keys))), b[, c("method", "nComp", "q", "p", "z")],
                 by = c("method", "nComp", "q", "p"))
      z_sum[m$row] <- z_sum[m$row] + m$z
    }
  }
  est <- cbind(keys, z = z_sum / length(files))

  # True beta, on the full-resolution grid (taken from the pFFPLS rows):
  truth <- first[first$method == "pFFPLS" & first$nComp == 1, c("q", "p", "z_true", "beta.num")]
  truth <- data.frame(method = "True beta", q = truth$q, p = truth$p, z = truth$z_true,
                      beta.num = truth$beta.num)

  method_order <- c("True beta", "pFFPLS", "FFPLS", "FFPLS_OB", "pFFR_I", "pFFR_RS")
  if (is.null(n_comps)) n_comps <- sort(unique(est$nComp))

  for (beta_num in unique(est$beta.num)) {

    out_dir   <- file.path(input_folder, "results_plots", paste0("2DBetaHeatmaps_", beta_num))
    owned_dir <- file.path(out_dir, "owned_scale")
    dir.create(owned_dir, recursive = TRUE, showWarnings = FALSE)

    for (H in n_comps) {

      d <- rbind(truth[truth$beta.num == beta_num, c("method", "q", "p", "z")],
                 est[est$beta.num == beta_num & est$nComp == H, c("method", "q", "p", "z")])
      d$method <- factor(d$method, levels = intersect(method_order, unique(d$method)))

      # tile size per method (pFFR_I uses its own, coarser grid)
      res <- do.call(rbind, lapply(split(d, d$method), function(x)
        data.frame(method = x$method[1], w = min(diff(sort(unique(x$p)))), h = min(diff(sort(unique(x$q)))))))
      d <- merge(d, res, by = "method")

      heat <- function(dd, lim, title = NULL) {
        ggplot(dd, aes(x = p, y = q, fill = z)) +
          geom_tile(aes(width = w, height = h)) +
          scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948",
                               midpoint = 0, limits = c(-lim, lim),
                               name = expression(beta(q, p))) +
          scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
          scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
          coord_fixed() +
          labs(x = "p (predictor domain)", y = "q (response domain)", title = title) +
          theme_bw(base_size = 12) +
          theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"),
                panel.spacing = grid::unit(1.2, "lines"))
      }

      # one panel, shared symmetric scale
      lim_all <- max(abs(d$z))
      g <- heat(d, lim_all, title = paste0("Mean estimated surfaces, PLS methods with H = ", H)) +
        facet_wrap(~ method, nrow = 2)
      for (ext in c("png", "pdf")) {
        ggsave(file.path(out_dir, paste0("all_methods_ncomp", H, ".", ext)), g,
               width = width, height = height, dpi = 300)
      }

      # one file per surface, own symmetric scale
      for (m in levels(d$method)) {
        dm <- d[d$method == m, ]
        gm <- heat(dm, max(abs(dm$z)),
                   title = if (m == "True beta") "True beta" else paste0(m, " (mean estimate, H = ", H, ")"))
        tag <- if (m == "True beta") "True_beta" else m
        for (ext in c("png", "pdf")) {
          ggsave(file.path(owned_dir, paste0(tag, "_ncomp", H, ".", ext)), gm,
                 width = 5.5, height = 4.5, dpi = 300)
        }
      }
    }
  }

  invisible(est)
}

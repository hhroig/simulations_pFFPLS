# Example curves for the paper: a few predictor curves X_i(p) without and with measurement
# noise, and their responses Y_i(q). Several layout options are written so that one can be
# chosen for the paper.
#
# The data are exactly those of one repetition of the study (repetition r uses
# set.seed(2026 + r)): the noise-free curves and responses come from generate_fofr_data()
# in simulations_buried_bands.R, and the noisy curves add the same measurement-noise draw
# the simulation adds (sd X_sd_error).
#
# The response is generated from the noise-free curves. In the figures, "noise-free" x is the
# curve before measurement noise and "observed" x includes it; "noise-free" y is the signal
# int (x_i(p) - mean x(p)) beta(q, p) dp and "observed" y adds the response noise.
#
# Usage, from the project root:
#   Rscript plot_example_curves.R [rep] [X_sd_error]
# Output: results_simulations/buried_bands/example_curves/  (option_*.png / .pdf)

library(ggplot2)

# Load only the function definitions of simulations_buried_bands.R (sourcing the whole file
# would run the simulation loop):
for (e in parse("simulations_buried_bands.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))) eval(e, envir = globalenv())
}

example_curve_data <- function(rep = 1, X_sd_error = 1.5, n_subjects = 3, seed_base = 2026,
                               band_sd = 0.035, band_share = 0.03) {
  p <- seq(0, 1, length.out = 100); q <- p
  set.seed(seed_base + rep)
  d <- generate_fofr_data(band_sd = band_sd, band_share = band_share)      # as in the simulation
  X_noisy <- d$X + matrix(stats::rnorm(length(d$X), 0, X_sd_error), nrow(d$X))

  # Split the noise-free curves into their dominant and buried parts, by redrawing the same
  # scores generate_fofr_data() drew (checked against d$X below):
  set.seed(seed_base + rep)
  bf <- make_band_functions(band_sd = band_sd)
  n  <- 150
  sd_b <- sqrt(band_share / (1 - band_share) * bf$dom_var / sum(bf$g_norm2))
  U  <- matrix(stats::rnorm(n * 10), nrow = n)
  Bs <- matrix(stats::rnorm(n * 2, sd = sd_b), nrow = n)
  dominant <- (U %*% t(0.5 * bf$fj(p)))[1:100, ]
  bands    <- (Bs[, 1] %o% bf$g_tilde[[1]](p) + Bs[, 2] %o% bf$g_tilde[[2]](p))[1:100, ]
  stopifnot(max(abs(dominant + bands - d$X)) < 1e-10)

  # noise-free response signal (as in generate_fofr_data(): trapezoidal rule, centred X)
  w <- rep(p[2] - p[1], 100); w[c(1, 100)] <- w[1] / 2
  Y_signal <- d$Xc %*% (t(d$beta_true) * w)

  # three subjects whose first band amplitude is low, median and high, so that the
  # effect of the buried bands on y is visible:
  b1 <- Bs[1:100, 1]
  ids <- order(b1)[round(seq(5, 95, length.out = n_subjects))]

  long <- function(M, domain, grid, part) {
    data.frame(subject = factor(rep(paste("Subject", seq_along(ids)), each = length(grid))),
               t = rep(grid, length(ids)), value = as.vector(t(M[ids, , drop = FALSE])),
               domain = domain, part = part)
  }
  rbind(long(d$X,      "x", p, "noise-free"),
        long(X_noisy,  "x", p, "with noise"),
        long(dominant, "x", p, "dominant part"),
        long(bands,    "x", p, "buried bands"),
        long(Y_signal, "y", q, "noise-free"),
        long(d$Y,      "y", q, "with noise"))
}

plot_example_curves <- function(rep = 1, X_sd_error = 1.5,
                                out_dir = "results_simulations/buried_bands/example_curves/") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dd <- example_curve_data(rep = rep, X_sd_error = X_sd_error)
  dd$obs <- factor(ifelse(dd$part == "with noise", "observed", "noise-free"), levels = c("noise-free", "observed"))
  # colours chosen by lightness so the subjects stay distinct in black-and-white print
  # (greyscale ~0.15 / 0.47 / 0.66): dark navy, vermillion, mustard
  cols <- c("Subject 1" = "#08306B", "Subject 2" = "#D55E00", "Subject 3" = "#E6AB02")
  base <- theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 12),
          legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(0, 0, 0, 0))
  save <- function(g, name, w, h) {
    for (ext in c("png", "pdf")) ggsave(file.path(out_dir, paste0(name, ".", ext)), g, width = w, height = h, dpi = 300, bg = "white")
  }
  combine <- function(..., nrow = 1) {
    ggpubr::ggarrange(..., nrow = nrow, ncol = length(list(...)) / nrow,
                      common.legend = TRUE, legend = "bottom", align = "hv")
  }
  # one panel: noise-free curves solid and thick; observed curves in the given style
  panel <- function(domain, parts, title, style = c("dashed", "translucent", "points", "solid")) {
    style <- match.arg(style)
    x  <- dd[dd$domain == domain & dd$part %in% parts, ]
    clean <- x[x$obs == "noise-free", ]; noisy <- x[x$obs == "observed", ]
    g  <- ggplot(x, aes(t, value, colour = subject))
    if (nrow(noisy) > 0) {
      g <- g + switch(style,
        dashed      = geom_line(data = noisy, aes(linetype = obs), linewidth = 0.5, alpha = 0.7),
        translucent = geom_line(data = noisy, aes(linetype = obs), linewidth = 0.5, alpha = 0.5),
        points      = geom_point(data = noisy, size = 0.8, alpha = 0.6),
        solid       = geom_line(data = noisy, linewidth = 0.6))
    }
    if (nrow(clean) > 0) {
      g <- g + if (style %in% c("dashed", "translucent")) geom_line(data = clean, aes(linetype = obs), linewidth = 0.9)
               else geom_line(data = clean, linewidth = 0.9)
    }
    g <- g + scale_colour_manual(values = cols, name = NULL) +
      labs(x = if (domain == "x") "p" else "q", y = NULL, title = title) + base
    if (style %in% c("dashed", "translucent")) {
      g <- g + scale_linetype_manual(values = c("noise-free" = "solid", "observed" = if (style == "dashed") "22" else "solid"),
                                     name = NULL, drop = FALSE) +
        guides(linetype = guide_legend(override.aes = list(alpha = c(1, if (style == "dashed") 0.7 else 0.5),
                                                           linewidth = c(0.9, 0.5))))
    }
    g
  }
  tx <- "Predictor curves  X(p)"; ty <- "Response curves  Y(q)"
  both <- c("noise-free", "with noise")

  # A: noise-free solid, observed dashed and translucent; x and y side by side
  save(combine(panel("x", both, tx, "dashed"), panel("y", both, ty, "dashed")), "option_A_solid_dashed", 9, 4.3)

  # B: as A, but the observed curves are thin translucent solid lines (easier to read than
  # dashes on a very jagged curve)
  save(combine(panel("x", both, tx, "translucent"), panel("y", both, ty, "translucent")),
       "option_B_solid_translucent", 9, 4.3)

  # C: 2 x 2 grid: top row noise-free, bottom row observed
  save(combine(panel("x", "noise-free", paste(tx, "- noise-free"), "solid"),
               panel("y", "noise-free", paste(ty, "- noise-free"), "solid"),
               panel("x", "with noise", paste(tx, "- observed"), "solid"),
               panel("y", "with noise", paste(ty, "- observed"), "solid"), nrow = 2),
       "option_C_grid_2x2", 9, 6.5)

  # D: decomposition: x (noise-free and observed), its buried bands alone, and y
  pb <- panel("x", "buried bands", "Buried bands of X(p) (the part that drives Y)", "solid") +
    geom_vline(xintercept = c(0.3, 0.7), colour = "grey55", linetype = "dotted")
  save(combine(panel("x", both, tx, "translucent"), pb, panel("y", both, ty, "translucent")),
       "option_D_decomposition", 13, 4.3)

  # E: observed curves as points, noise-free as lines
  save(ggpubr::annotate_figure(combine(panel("x", both, tx, "points"), panel("y", both, ty, "points")),
                               bottom = ggpubr::text_grob("Lines: noise-free; points: observed.", size = 10)),
       "option_E_points", 9, 4.5)

  invisible(dd)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  plot_example_curves(rep        = if (length(args) >= 1) as.integer(args[1]) else 1,
                      X_sd_error = if (length(args) >= 2) as.numeric(args[2]) else 1.5)
}

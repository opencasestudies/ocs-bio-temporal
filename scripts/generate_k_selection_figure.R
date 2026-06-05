#!/usr/bin/env Rscript

# Generate the K-selection figure used in the Data Visualization section.
# The figure is built only from the lightweight source table committed with the
# case study, so it can be regenerated without running the full CoGAPS sweep.

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
fig_dir <- file.path(root, "data", "processed", "figures")
source_path <- file.path(fig_dir, "source_tables", "figure_01_k_selection_source.csv")
out_png <- file.path(fig_dir, "figure_01_k_selection.png")
out_svg <- file.path(fig_dir, "figure_01_k_selection.svg")

k_summary <- utils::read.csv(source_path, stringsAsFactors = FALSE)
k_summary <- k_summary[order(k_summary$K), ]

as_bool <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  tolower(as.character(x)) %in% c("true", "1", "yes")
}

k_summary$selected_for_case_study <- as_bool(k_summary$selected_for_case_study)
k_summary$on_stability_plateau <- as_bool(k_summary$on_stability_plateau)

draw_figure <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  layout(matrix(c(1, 2), nrow = 2), heights = c(1.15, 1.0))
  par(
    family = "sans",
    oma = c(0.5, 0.5, 3.4, 0.5),
    mar = c(4.8, 5.1, 3.3, 1.7),
    cex.axis = 1.05,
    cex.lab = 1.12
  )

  point_colors <- rep("#b8b8b8", nrow(k_summary))
  point_colors[k_summary$on_stability_plateau & k_summary$max_eta_timepoint < 0.12] <- "#4c72b0"
  point_colors[k_summary$selected_for_case_study] <- "#c44e52"

  plot(
    k_summary$stability_core,
    k_summary$max_eta_timepoint,
    pch = 21,
    bg = point_colors,
    col = "white",
    cex = 1.75,
    lwd = 1.1,
    xlab = "Seed-to-seed stability core",
    ylab = "Maximum temporal effect size",
    main = "A. Stability and temporal resolution",
    xlim = c(min(k_summary$stability_core) - 0.015, 1.005),
    ylim = c(-0.02, max(k_summary$max_eta_timepoint) + 0.08),
    las = 1
  )
  grid(col = "#e6edf3", lwd = 1)
  plateau_threshold <- min(k_summary$stability_core[k_summary$on_stability_plateau])
  abline(v = plateau_threshold, lty = 2, col = "#6b7280", lwd = 1.2)
  text(
    k_summary$stability_core + 0.003,
    k_summary$max_eta_timepoint + 0.006,
    labels = paste0("K=", k_summary$K),
    cex = 0.86,
    col = "#263238"
  )

  selected <- k_summary[k_summary$selected_for_case_study, ][1, ]
  compact <- k_summary[k_summary$K == 5, ][1, ]
  arrows(
    selected$stability_core - 0.025,
    selected$max_eta_timepoint - 0.12,
    selected$stability_core - 0.003,
    selected$max_eta_timepoint - 0.01,
    length = 0.08,
    col = "#9f2f20",
    lwd = 1.3
  )
  text(
    selected$stability_core - 0.04,
    selected$max_eta_timepoint - 0.14,
    "selected K=10\nstable + temporal",
    col = "#9f2f20",
    font = 2,
    cex = 0.95,
    adj = c(0, 0.5)
  )
  arrows(
    compact$stability_core - 0.02,
    compact$max_eta_timepoint + 0.12,
    compact$stability_core - 0.003,
    compact$max_eta_timepoint + 0.012,
    length = 0.08,
    col = "#2f5f9f",
    lwd = 1.2
  )
  text(
    compact$stability_core - 0.055,
    compact$max_eta_timepoint + 0.13,
    "K=5\nstable, compact",
    col = "#2f5f9f",
    cex = 0.95,
    adj = c(0, 0.5)
  )
  legend(
    "bottomleft",
    legend = c("selected K=10", "stable but compact", "other tested K"),
    pt.bg = c("#c44e52", "#4c72b0", "#b8b8b8"),
    pch = 21,
    col = "white",
    bty = "n",
    cex = 0.9
  )

  par(mar = c(6.0, 5.1, 3.3, 1.7))
  bar_colors <- ifelse(k_summary$selected_for_case_study, "#c44e52", "#55a868")
  bars <- barplot(
    k_summary$goal_aligned_score,
    names.arg = paste0("K=", k_summary$K),
    col = bar_colors,
    border = "white",
    las = 2,
    ylim = c(0, max(k_summary$goal_aligned_score) * 1.18),
    ylab = "Goal-aligned score",
    main = "B. Revised K-selection ranking"
  )
  grid(nx = NA, ny = NULL, col = "#e6edf3", lwd = 1)
  selected_index <- which(k_summary$selected_for_case_study)[1]
  text(
    bars[selected_index],
    k_summary$goal_aligned_score[selected_index] + 0.035,
    labels = "selected",
    col = "#263238",
    cex = 0.95
  )

  mtext(
    "K selection: stability, temporal resolution, and selected rank",
    outer = TRUE,
    side = 3,
    line = 1.3,
    cex = 1.28,
    font = 2,
    col = "#12263a"
  )
}

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

grDevices::png(out_png, width = 1600, height = 1900, res = 190)
draw_figure()
grDevices::dev.off()

grDevices::svg(out_svg, width = 8.5, height = 10.2)
draw_figure()
grDevices::dev.off()

message("Wrote ", out_png)
message("Wrote ", out_svg)

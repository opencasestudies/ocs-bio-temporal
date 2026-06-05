#!/usr/bin/env Rscript

# Generate the temporal small-multiples figure used in the Data Visualization
# section. The figure is built from the lightweight pattern-by-time source table,
# so it can be regenerated without rerunning CoGAPS.

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
if (length(script_path) == 0 || is.na(script_path)) {
  script_path <- normalizePath("scripts/generate_pattern_time_small_multiples.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
fig_dir <- file.path(root, "data", "processed", "figures")
source_path <- file.path(fig_dir, "source_tables", "figure_03_pattern_usage_by_time_matrix.csv")
out_png <- file.path(fig_dir, "figure_03_pattern_usage_by_time_small_multiples.png")
out_svg <- file.path(fig_dir, "figure_03_pattern_usage_by_time_small_multiples.svg")

suppressPackageStartupMessages({
  library(ggplot2)
})

timepoint_order <- c("D0", "D2", "D4", "D6", "D8", "D10", "D14/15", "D28")

pattern_time_matrix <- utils::read.csv(
  source_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pattern_time_long <- stats::reshape(
  pattern_time_matrix,
  varying = timepoint_order,
  v.names = "mean_usage",
  timevar = "timepoint",
  times = timepoint_order,
  direction = "long"
)
rownames(pattern_time_long) <- NULL
pattern_time_long$timepoint <- factor(pattern_time_long$timepoint, levels = timepoint_order)
pattern_time_long$pattern <- factor(pattern_time_long$pattern, levels = pattern_time_matrix$pattern)
pattern_time_long$day_index <- as.integer(pattern_time_long$timepoint)
pattern_time_long$is_ifn_candidate <- pattern_time_long$pattern == "Pattern3"

pattern_plot <- ggplot(
  pattern_time_long,
  aes(x = day_index, y = mean_usage, group = pattern)
) +
  annotate(
    "rect",
    xmin = which(timepoint_order == "D10") - 0.35,
    xmax = which(timepoint_order == "D14/15") + 0.35,
    ymin = -Inf,
    ymax = Inf,
    fill = "#fee2e2",
    alpha = 0.75
  ) +
  geom_line(aes(color = is_ifn_candidate), linewidth = 0.9) +
  geom_point(aes(color = is_ifn_candidate), size = 2.0) +
  facet_wrap(~ pattern, ncol = 5) +
  scale_color_manual(values = c(`FALSE` = "#4f6f8f", `TRUE` = "#c23b22"), guide = "none") +
  scale_x_continuous(
    breaks = seq_along(timepoint_order),
    labels = timepoint_order,
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  scale_y_continuous(limits = c(0, 0.31), expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Pattern usage trajectories across experimental infection days",
    subtitle = "Mean selected-model CoGAPS pattern scores in the balanced discovery set",
    x = "Experimental day",
    y = "Mean pattern usage"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", color = "#12263a", size = 20),
    plot.subtitle = element_text(color = "#425f7a", size = 14),
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1.1, "lines")
  )

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

grDevices::png(out_png, width = 1800, height = 1160, res = 180)
print(pattern_plot)
grDevices::dev.off()

grDevices::svg(out_svg, width = 10, height = 6.45)
print(pattern_plot)
grDevices::dev.off()

message("Wrote ", out_png)
message("Wrote ", out_svg)

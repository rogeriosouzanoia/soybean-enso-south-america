#!/usr/bin/env Rscript

# Figure 4
# Regional climate regimes associated with ENSO and soybean yield.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else {
  override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
  normalizePath(if (nzchar(override)) override else default_project_dir,
                mustWork = TRUE)
}
data_dir <- if (length(args) >= 2L) args[[2]] else
  file.path(project_dir, "outputs", "figure_data", "figure_4")

find_one <- function(root, candidates) {
  direct <- file.path(root, candidates)
  if (any(file.exists(direct))) {
    return(normalizePath(direct[file.exists(direct)][[1]]))
  }
  all_files <- list.files(root, recursive = TRUE, full.names = TRUE)
  hits <- all_files[basename(all_files) %in% candidates]
  if (length(hits) == 0L) {
    stop("Could not find any of: ", paste(candidates, collapse = ", "),
         call. = FALSE)
  }
  normalizePath(hits[order(nchar(hits), hits)][[1]], mustWork = TRUE)
}

paths <- list(
  scores = find_one(data_dir, c("Figure_4_local_climate_PCA_scores.csv.gz")),
  loadings = find_one(data_dir, c("Figure_4_macroregion_PCA_loadings.csv")),
  variance = find_one(data_dir, c("Figure_4_macroregion_PCA_variance.csv")),
  productivity = find_one(data_dir, c(
    "Figure_4_productivity_supplementary_vectors.csv"
  )),
  centroids = find_one(data_dir, c("Figure_4_ENSO_phase_centroids.csv")),
  tests = find_one(data_dir, c("Figure_4_ENSO_phase_separation_tests.csv"))
)

scores <- read_csv(paths$scores, show_col_types = FALSE)
loadings <- read_csv(paths$loadings, show_col_types = FALSE)
variance <- read_csv(paths$variance, show_col_types = FALSE)
productivity <- read_csv(paths$productivity, show_col_types = FALSE)
centroids <- read_csv(paths$centroids, show_col_types = FALSE)
tests <- read_csv(paths$tests, show_col_types = FALSE)

required_scores <- c(
  "macroregion", "ENSO_phase", "PC1", "PC2", "analysis_weight",
  "productivity_anomaly_percent"
)
missing_scores <- setdiff(required_scores, names(scores))
if (length(missing_scores) > 0L) {
  stop("PCA score file is missing: ", paste(missing_scores, collapse = ", "),
       call. = FALSE)
}

macroregion_order <- c(
  "Argentina", "Southern Brazil", "Central Brazil",
  "MATOPIBA and northeastern frontier"
)
phase_order <- c("El Nino", "Neutral", "La Nina")
phase_labels <- c(
  "El Nino" = "El Niño", "Neutral" = "Neutral", "La Nina" = "La Niña"
)
phase_colours <- c(
  "El Nino" = "#D6604D", "Neutral" = "#8A8A8A", "La Nina" = "#4477AA"
)

scores <- scores |>
  mutate(
    macroregion = factor(macroregion, levels = macroregion_order),
    ENSO_phase = factor(ENSO_phase, levels = phase_order)
  )
centroids <- centroids |>
  mutate(
    macroregion = factor(macroregion, levels = macroregion_order),
    ENSO_phase = factor(ENSO_phase, levels = phase_order)
  )

theme_rogerio <- function(base_size = 13, base_family = "serif") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(colour = "black"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(colour = "black", linewidth = 0.55),
      axis.ticks = element_line(colour = "black", linewidth = 0.45),
      axis.ticks.length = unit(2.2, "mm"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      plot.title = element_text(
        size = base_size + 1, face = "bold", hjust = 0,
        margin = margin(b = 7)
      ),
      plot.tag = element_text(face = "bold", size = base_size + 2),
      plot.tag.position = c(0.005, 0.995),
      plot.margin = margin(12, 15, 10, 12)
    )
}

weighted_quantile <- function(x, w, probabilities) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  ordering <- order(x)
  x <- x[ordering]
  cumulative <- cumsum(w[ordering]) / sum(w)
  vapply(probabilities, function(probability) {
    x[which(cumulative >= probability)[1]]
  }, numeric(1))
}

format_p <- function(p) {
  if (!is.finite(p)) return("P unavailable")
  if (p < 0.001) return("Pperm < 0.001")
  paste0("Pperm = ", formatC(p, format = "f", digits = 3))
}

prepare_vectors <- function(region_scores, region_loadings, region_productivity) {
  x_limits <- weighted_quantile(
    region_scores$PC1, region_scores$analysis_weight, c(0.002, 0.998)
  )
  y_limits <- weighted_quantile(
    region_scores$PC2, region_scores$analysis_weight, c(0.002, 0.998)
  )
  x_padding <- diff(x_limits) * 0.06
  y_padding <- diff(y_limits) * 0.06
  x_limits <- x_limits + c(-x_padding, x_padding)
  y_limits <- y_limits + c(-y_padding, y_padding)

  selected <- region_loadings |>
    arrange(loading_rank) |>
    slice_head(n = 6L)
  maximum_x <- max(abs(selected$PC1_loading), na.rm = TRUE)
  maximum_y <- max(abs(selected$PC2_loading), na.rm = TRUE)
  selected <- selected |>
    mutate(
      arrow_x = PC1_loading / maximum_x * 0.31 * diff(x_limits),
      arrow_y = PC2_loading / maximum_y * 0.31 * diff(y_limits),
      label_x = arrow_x * 1.08,
      label_y = arrow_y * 1.08,
      vector_type = "Climate variable"
    )

  productivity_vector <- region_productivity |>
    transmute(
      variable_label = "Productivity anomaly",
      arrow_x = PC1_correlation * 0.43 * diff(x_limits),
      arrow_y = PC2_correlation * 0.43 * diff(y_limits),
      label_x = arrow_x * 1.08,
      label_y = arrow_y * 1.08,
      vector_type = "Productivity"
    )

  list(
    limits = list(x = x_limits, y = y_limits),
    climate = selected,
    productivity = productivity_vector
  )
}

make_region_plot <- function(selected_macroregion) {
  region_scores <- scores |> filter(macroregion == selected_macroregion)
  region_loadings <- loadings |> filter(macroregion == selected_macroregion)
  region_productivity <- productivity |>
    filter(macroregion == selected_macroregion)
  region_centroids <- centroids |> filter(macroregion == selected_macroregion)
  region_variance <- variance |>
    filter(macroregion == selected_macroregion, component %in% c("PC1", "PC2"))
  region_test <- tests |> filter(macroregion == selected_macroregion)

  vectors <- prepare_vectors(region_scores, region_loadings, region_productivity)
  pc1_percent <- 100 * region_variance$variance_fraction[
    region_variance$component == "PC1"
  ]
  pc2_percent <- 100 * region_variance$variance_fraction[
    region_variance$component == "PC2"
  ]
  annotation <- paste0(
    "Climate space: ", sprintf("%.1f", pc1_percent + pc2_percent), "%   ",
    format_p(region_test$permutation_p_value[[1]])
  )

  ggplot(region_scores, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, colour = "#D3D3D3", linewidth = 0.45) +
    geom_vline(xintercept = 0, colour = "#D3D3D3", linewidth = 0.45) +
    geom_point(
      aes(colour = ENSO_phase, size = analysis_weight),
      alpha = 0.055, stroke = 0, show.legend = c(colour = TRUE, size = FALSE)
    ) +
    geom_segment(
      data = vectors$climate,
      aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
      inherit.aes = FALSE, colour = "#666666", linewidth = 0.55,
      alpha = 0.72, arrow = arrow(length = unit(1.8, "mm"), type = "closed")
    ) +
    geom_text_repel(
      data = vectors$climate,
      aes(label = variable_label, x = label_x, y = label_y),
      inherit.aes = FALSE, colour = "#4D4D4D", size = 3.25,
      family = "serif", segment.colour = "#B5B5B5", segment.size = 0.3,
      min.segment.length = 0, box.padding = 0.25, point.padding = 0.1,
      max.overlaps = Inf, seed = 20260907
    ) +
    geom_segment(
      data = vectors$productivity,
      aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
      inherit.aes = FALSE, colour = "black", linewidth = 0.9,
      arrow = arrow(length = unit(2.2, "mm"), type = "closed")
    ) +
    geom_text_repel(
      data = vectors$productivity,
      aes(label = variable_label, x = label_x, y = label_y),
      inherit.aes = FALSE, colour = "black", fontface = "bold", size = 3.4,
      family = "serif", segment.colour = "black", segment.size = 0.35,
      min.segment.length = 0, box.padding = 0.3, max.overlaps = Inf,
      seed = 20260907
    ) +
    geom_errorbar(
      data = region_centroids,
      aes(x = PC1, ymin = PC2_lower, ymax = PC2_upper, colour = ENSO_phase),
      inherit.aes = FALSE, width = 0, linewidth = 0.8, alpha = 0.95
    ) +
    geom_errorbarh(
      data = region_centroids,
      aes(y = PC2, xmin = PC1_lower, xmax = PC1_upper, colour = ENSO_phase),
      inherit.aes = FALSE, height = 0, linewidth = 0.8, alpha = 0.95
    ) +
    geom_point(
      data = region_centroids,
      aes(x = PC1, y = PC2, fill = ENSO_phase),
      inherit.aes = FALSE, shape = 21, size = 4.3, stroke = 0.7,
      colour = "white", show.legend = FALSE
    ) +
    annotate(
      "text", x = vectors$limits$x[1], y = vectors$limits$y[2],
      label = annotation, hjust = 0, vjust = 1, family = "serif",
      size = 3.6, colour = "#333333"
    ) +
    scale_colour_manual(
      values = phase_colours, labels = phase_labels, drop = FALSE,
      guide = guide_legend(
        nrow = 1, byrow = TRUE,
        override.aes = list(alpha = 1, size = 3.4, stroke = 0)
      )
    ) +
    scale_fill_manual(values = phase_colours, labels = phase_labels, drop = FALSE) +
    scale_size_continuous(range = c(0.18, 0.85), guide = "none") +
    coord_cartesian(
      xlim = vectors$limits$x, ylim = vectors$limits$y, clip = "off"
    ) +
    labs(
      title = selected_macroregion,
      x = paste0("PC1 (", sprintf("%.1f", pc1_percent), "%)"),
      y = paste0("PC2 (", sprintf("%.1f", pc2_percent), "%)")
    ) +
    theme_rogerio()
}

plots <- lapply(macroregion_order, make_region_plot)

f4 <- wrap_plots(plots, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "top",
    legend.justification = "center",
    plot.tag = element_text(face = "bold", family = "serif", size = 15)
  )

# Run this final line in RStudio to display Figure 4.
f4

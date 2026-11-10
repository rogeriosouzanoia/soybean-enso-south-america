#!/usr/bin/env Rscript

# Figure 4
# Regional climate regimes associated with ENSO and soybean yield.
# FINAL LABEL REVISION: full variable names, white label boxes, and lower-right
# PCA summary statistics. If Y/RWS/HWS appear, an older script was executed.

packages <- c(
  "dplyr", "readr", "ggplot2", "ggrepel", "patchwork", "scales",
  "sf", "rnaturalearth", "systemfonts"
)
missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(grid)
  library(sf)
  library(rnaturalearth)
  library(systemfonts)
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

geometry_file <- find_one(project_dir, c(
  "ENSO_phase_specific_spatial_layers.gpkg",
  "ENSO_phase_specific_spatial_layers(1).gpkg"
))

scores <- read_csv(paths$scores, show_col_types = FALSE)
loadings <- read_csv(paths$loadings, show_col_types = FALSE)
variance <- read_csv(paths$variance, show_col_types = FALSE)
productivity <- read_csv(paths$productivity, show_col_types = FALSE)
centroids <- read_csv(paths$centroids, show_col_types = FALSE)
tests <- read_csv(paths$tests, show_col_types = FALSE)

required_scores <- c(
  "country", "spatial_unit_id", "macroregion", "ENSO_phase",
  "PC1", "PC2", "analysis_weight",
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

# Canonical climate-index palette. Keep these exact values in every figure in
# which indices are encoded by colour. CWS is ochre-brown to communicate crop
# water deficit; blue is reserved for precipitation-related conditions.
index_colours <- c(
  "Temperature" = "#D73027",
  "VPD" = "#E66101",
  "Precipitation" = "#56B4E9",
  "Heavy rainfall" = "#5E3C99",
  "CWS" = "#A6761D",
  "Hot-dry conditions" = "#8C2D04",
  "Hot-wet conditions" = "#008837"
)
# Light fills keep the locator map subordinate to the PCA panels. Darker
# companions are used only for region names and arrows.
macroregion_colours <- c(
  "Argentina" = "#F2C2B8",
  "Southern Brazil" = "#B8CEE5",
  "Central Brazil" = "#C9E1EC",
  "MATOPIBA and northeastern frontier" = "#D8C9E6"
)
macroregion_label_colours <- c(
  "Argentina" = "#B34A3A",
  "Southern Brazil" = "#315F96",
  "Central Brazil" = "#4F899F",
  "MATOPIBA and northeastern frontier" = "#76519A"
)

available_fonts <- unique(systemfonts::system_fonts()$family)
figure_font_family <- dplyr::case_when(
  "Calibri" %in% available_fonts ~ "Calibri",
  "Carlito" %in% available_fonts ~ "Carlito",
  "Arial" %in% available_fonts ~ "Arial",
  TRUE ~ "sans"
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

theme_rogerio <- function(base_size = 12.4, base_family = figure_font_family) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(colour = "black"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(colour = "black", linewidth = 0.90),
      axis.ticks = element_line(colour = "black", linewidth = 0.85),
      axis.ticks.length = unit(2.8, "mm"),
      axis.title = element_text(
        size = base_size + 0.2, face = "bold", colour = "black"
      ),
      axis.text = element_text(
        size = base_size - 0.5, face = "bold", colour = "black"
      ),
      legend.position = "top",
      legend.text = element_text(size = base_size - 0.5, face = "bold"),
      legend.title = element_text(size = base_size, face = "bold"),
      plot.title = element_text(
        size = base_size + 1, face = "bold", hjust = 0,
        margin = margin(b = 7)
      ),
      plot.tag = element_text(face = "bold", size = base_size + 2),
      plot.tag.position = c(0.005, 0.995),
      plot.margin = margin(7, 8, 6, 7)
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

standardize_id <- function(x, country) {
  digits <- gsub("[^0-9]", "", as.character(x))
  widths <- ifelse(country == "Brazil", 7L, 5L)
  mapply(function(value, width) {
    if (is.na(value) || value == "") return(NA_character_)
    sprintf(paste0("%0", width, "d"), as.integer(value))
  }, digits, widths, USE.NAMES = FALSE)
}

standardize_geometry <- function(x) {
  active <- attr(x, "sf_column")
  if (is.null(active) || !active %in% names(x)) {
    stop("Active sf geometry column was not found.", call. = FALSE)
  }
  if (active != "geometry") names(x)[names(x) == active] <- "geometry"
  st_geometry(x) <- "geometry"
  x
}

short_vector_label <- function(x) {
  dplyr::recode(
    x,
    "VPD" = "Vapor pressure\ndeficit",
    "Vapor pressure deficit" = "VPD",
    "Hot and dry" = "Hot-dry conditions",
    "Hot and dry conditions" = "Hot-dry\nconditions",
    "Hot and wet" = "Hot-wet conditions",
    "Hot and wet conditions" = "Hot-wet\nconditions",
    "Cold and wet" = "Cold-wet\nconditions",
    "Cold and wet conditions" = "Cold-wet\nconditions",
    "Consecutive dry days" = "Consecutive\ndry days",
    "Hot water stress" = "CWS",
    "Reproductive water satisfaction" = "CWS",
    "Reproductive-stage water adequacy" = "CWS",
    "Crop water stress" = "CWS",
    .default = x
  )
}

geometry_layer <- st_read(
  geometry_file, layer = "local_phase_effects", quiet = TRUE
) |>
  standardize_geometry() |>
  st_make_valid() |>
  mutate(spatial_unit_id_join = standardize_id(spatial_unit_id, country))

unit_geometry <- geometry_layer |>
  arrange(country, spatial_unit_id_join) |>
  group_by(country, spatial_unit_id_join) |>
  slice(1L) |>
  ungroup() |>
  select(country, spatial_unit_id_join, geometry) |>
  st_transform(3857) |>
  st_simplify(dTolerance = 150, preserveTopology = TRUE) |>
  st_transform(4326)

macroregion_sf <- scores |>
  distinct(country, spatial_unit_id, macroregion) |>
  mutate(
    spatial_unit_id_join = standardize_id(spatial_unit_id, country),
    macroregion = factor(macroregion, levels = macroregion_order)
  ) |>
  inner_join(unit_geometry, by = c("country", "spatial_unit_id_join")) |>
  st_as_sf()

# Dissolved region polygons provide clean macroregional outlines and stable
# positions for labels. point_on_surface keeps every label inside its region.
macroregion_boundaries <- macroregion_sf |>
  group_by(macroregion) |>
  summarise(do_union = TRUE, .groups = "drop") |>
  st_make_valid()

macroregion_labels <- macroregion_boundaries |>
  st_point_on_surface() |>
  mutate(
    region_label = recode(
      as.character(macroregion),
      "Argentina" = "Argentina",
      "Southern Brazil" = "South Brazil",
      "Central Brazil" = "Central Brazil",
      "MATOPIBA and northeastern frontier" = "MATOPIBA–NE"
    )
  )

# External labels avoid crowding the small locator map. Label positions are
# fixed for a stable publication layout; target coordinates come from each
# dissolved macroregion and therefore remain data driven.
macroregion_callout_xy <- st_coordinates(macroregion_labels)
macroregion_callouts <- macroregion_labels |>
  st_drop_geometry() |>
  mutate(
    target_x = macroregion_callout_xy[, 1],
    target_y = macroregion_callout_xy[, 2],
    label_x = case_when(
      as.character(macroregion) == "Argentina" ~ -73.0,
      TRUE ~ -35.5
    ),
    label_y = case_when(
      as.character(macroregion) == "Argentina" ~ -36.0,
      as.character(macroregion) == "Southern Brazil" ~ -28.0,
      as.character(macroregion) == "Central Brazil" ~ -17.0,
      TRUE ~ -7.0
    ),
    label_hjust = if_else(as.character(macroregion) == "Argentina", 0, 1)
  )

south_america <- rnaturalearth::ne_countries(
  continent = "South America", scale = "medium", returnclass = "sf"
) |>
  st_make_valid()

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
    slice_head(n = 4L)
  maximum_x <- max(abs(selected$PC1_loading), na.rm = TRUE)
  maximum_y <- max(abs(selected$PC2_loading), na.rm = TRUE)
  selected <- selected |>
    mutate(
      variable_label = short_vector_label(variable_label),
      arrow_x = PC1_loading / maximum_x * 0.31 * diff(x_limits),
      arrow_y = PC2_loading / maximum_y * 0.31 * diff(y_limits),
      label_x = arrow_x * 1.08,
      label_y = arrow_y * 1.08,
      vector_type = "Climate variable"
    )

  productivity_vector <- region_productivity |>
    transmute(
      variable_label = "Yield",
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

make_region_plot <- function(selected_macroregion, panel_letter) {
  region_scores <- scores |>
    filter(macroregion == selected_macroregion) |>
    mutate(draw_order = if_else(ENSO_phase == "Neutral", 1L, 2L)) |>
    arrange(draw_order)
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
    sprintf("%.1f", pc1_percent + pc2_percent), "%  |  ",
    format_p(region_test$permutation_p_value[[1]])
  )

  ggplot(region_scores, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, colour = "#D3D3D3", linewidth = 0.45) +
    geom_vline(xintercept = 0, colour = "#D3D3D3", linewidth = 0.45) +
    geom_point(
      aes(colour = ENSO_phase, size = analysis_weight),
      alpha = 0.20, stroke = 0, show.legend = c(colour = TRUE, size = FALSE)
    ) +
    geom_segment(
      data = vectors$climate,
      aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
      inherit.aes = FALSE, colour = "#666666", linewidth = 0.55,
      alpha = 0.72, arrow = arrow(length = unit(1.8, "mm"), type = "closed")
    ) +
    geom_label_repel(
      data = vectors$climate,
      aes(label = variable_label, x = label_x, y = label_y),
      inherit.aes = FALSE, colour = "#4D4D4D", size = 3.15,
      family = figure_font_family, fill = scales::alpha("white", 0.88),
      linewidth = 0, label.padding = unit(0.10, "lines"),
      label.r = unit(0.04, "lines"), lineheight = 0.92,
      segment.colour = "#AFAFAF", segment.size = 0.3,
      min.segment.length = 0, box.padding = 0.32, point.padding = 0.12,
      force = 1.4, max.overlaps = Inf, seed = 20260907
    ) +
    geom_segment(
      data = vectors$productivity,
      aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
      inherit.aes = FALSE, colour = "black", linewidth = 0.9,
      arrow = arrow(length = unit(2.2, "mm"), type = "closed")
    ) +
    geom_label_repel(
      data = vectors$productivity,
      aes(label = variable_label, x = label_x, y = label_y),
      inherit.aes = FALSE, colour = "black", fontface = "bold", size = 3.3,
      family = figure_font_family, fill = scales::alpha("white", 0.90),
      linewidth = 0, label.padding = unit(0.11, "lines"),
      label.r = unit(0.04, "lines"), segment.colour = "black",
      segment.size = 0.35, min.segment.length = 0, box.padding = 0.34,
      force = 1.5, max.overlaps = Inf,
      seed = 20260907
    ) +
    geom_errorbar(
      data = region_centroids,
      aes(x = PC1, ymin = PC2_lower, ymax = PC2_upper, colour = ENSO_phase),
      inherit.aes = FALSE, width = 0, linewidth = 0.8, alpha = 0.95
    ) +
    geom_errorbar(
      data = region_centroids,
      aes(y = PC2, xmin = PC1_lower, xmax = PC1_upper, colour = ENSO_phase),
      inherit.aes = FALSE, width = 0, orientation = "y",
      linewidth = 0.8, alpha = 0.95
    ) +
    geom_point(
      data = region_centroids,
      aes(x = PC1, y = PC2, fill = ENSO_phase),
      inherit.aes = FALSE, shape = 21, size = 4.3, stroke = 0.7,
      colour = "white", show.legend = FALSE
    ) +
    annotate(
      "label", x = vectors$limits$x[2], y = vectors$limits$y[1],
      label = annotation, hjust = 1, vjust = 0, family = figure_font_family,
      size = 3.15, colour = "#333333",
      fill = scales::alpha("white", 0.90), linewidth = 0,
      label.padding = unit(0.12, "lines"), label.r = unit(0.04, "lines")
    ) +
    scale_colour_manual(
      values = phase_colours, labels = phase_labels, drop = FALSE,
      guide = guide_legend(
        nrow = 1, byrow = TRUE,
        override.aes = list(alpha = 1, size = 3.4, stroke = 0)
      )
    ) +
    scale_fill_manual(values = phase_colours, labels = phase_labels, drop = FALSE) +
    scale_size_continuous(range = c(0.75, 2.70), guide = "none") +
    coord_cartesian(
      xlim = vectors$limits$x, ylim = vectors$limits$y, clip = "off"
    ) +
    labs(
      title = paste0(
        panel_letter, ") ",
        recode(
          selected_macroregion,
          "MATOPIBA and northeastern frontier" = "MATOPIBA Brazil"
        )
      ),
      x = paste0("PC1 (", sprintf("%.1f", pc1_percent), "%)"),
      y = paste0("PC2 (", sprintf("%.1f", pc2_percent), "%)")
    ) +
    theme_rogerio() +
    theme(
      plot.title = element_text(
        size = 13.2, face = "bold", hjust = 0,
        margin = margin(b = 8, l = 0)
      ),
      legend.position = "none"
    )
}

plots <- Map(
  function(selected_macroregion, panel_letter) {
    make_region_plot(selected_macroregion, panel_letter)
  },
  macroregion_order, letters[1:4]
)

left_grid <- wrap_plots(plots, ncol = 2)

panel_e <- ggplot() +
  geom_sf(
    data = south_america, fill = "grey98", colour = "black", linewidth = 0.36
  ) +
  geom_sf(
    data = macroregion_sf, aes(fill = macroregion),
    colour = NA, linewidth = 0
  ) +
  geom_sf(
    data = macroregion_boundaries,
    fill = NA, colour = "white", linewidth = 0.72
  ) +
  geom_segment(
    data = macroregion_callouts,
    aes(
      x = label_x, y = label_y, xend = target_x, yend = target_y,
      colour = macroregion
    ),
    inherit.aes = FALSE, linewidth = 0.58,
    arrow = arrow(length = unit(1.35, "mm"), type = "closed"),
    show.legend = FALSE
  ) +
  geom_label(
    data = macroregion_callouts,
    aes(
      x = label_x, y = label_y, label = region_label,
      colour = macroregion, hjust = label_hjust
    ),
    inherit.aes = FALSE, family = figure_font_family,
    fontface = "bold", size = 2.7,
    fill = scales::alpha("white", 0.94), linewidth = 0.35,
    label.padding = unit(0.65, "mm"), lineheight = 0.92,
    show.legend = FALSE
  ) +
  geom_sf(data = south_america, fill = NA, colour = "black", linewidth = 0.42) +
  scale_fill_manual(
    values = macroregion_colours,
    labels = c(
      "Argentina", "Southern Brazil", "Central Brazil",
      "MATOPIBA and Northeast"
    ),
    drop = FALSE
  ) +
  scale_colour_manual(values = macroregion_label_colours, drop = FALSE) +
  coord_sf(xlim = c(-75, -34), ylim = c(-42, 1), expand = FALSE) +
  labs(fill = NULL) +
  theme_void(base_family = figure_font_family) +
  theme(
    legend.position = "none",
    plot.background = element_rect(
      fill = scales::alpha("white", 0.97), colour = "black", linewidth = 0.55
    ),
    plot.margin = margin(3, 3, 3, 3)
  )

enso_legend <- ggplot(
  data.frame(
    ENSO_phase = factor(phase_order, levels = phase_order),
    x = seq_along(phase_order), y = 1
  ),
  aes(x, y, colour = ENSO_phase)
) +
  geom_point(size = 3.2, alpha = 0) +
  scale_colour_manual(
    values = phase_colours, labels = phase_labels, drop = FALSE
  ) +
  guides(colour = guide_legend(
    title = "ENSO phase",
    nrow = 1, byrow = TRUE,
    title.position = "left", title.vjust = 0.5,
    override.aes = list(alpha = 1, size = 4.3)
  )) +
  theme_void(base_family = figure_font_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12.5),
    legend.text = element_text(face = "bold", size = 12.0),
    legend.spacing.x = unit(3.5, "mm"),
    legend.key.width = unit(5.0, "mm")
  )

# Keep only a narrow right-hand gutter. Cropping the locator to the soybean
# domain makes the map readable without shrinking the four PCA panels.
main_canvas <- (left_grid | plot_spacer()) +
  plot_layout(widths = c(1, 0.05))

main_layout <- main_canvas +
  inset_element(
    panel_e,
    left = 0.755, bottom = 0.595, right = 0.995, top = 0.995,
    align_to = "full", on_top = TRUE, clip = TRUE
  )

f4 <- main_layout / enso_legend +
  plot_layout(heights = c(1, 0.075)) &
  theme(
    plot.tag = element_text(
      face = "bold", family = figure_font_family, size = 13
    )
  )

if (interactive() && .Platform$OS.type == "windows") {
  grDevices::windows(width = 15.5, height = 10.5, record = FALSE)
  print(f4)
} else {
  print(f4)
}

# The assembled figure remains available as the object f4.

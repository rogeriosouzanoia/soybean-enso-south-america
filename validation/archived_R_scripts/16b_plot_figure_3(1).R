#!/usr/bin/env Rscript

# Figure 3
# Climate mechanisms connecting ENSO to soybean yield.

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1]] else {
  override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
  if (nzchar(override)) override else default_project_dir
}
figure_data_dir <- if (length(args) >= 2L) args[[2]] else
  file.path(project_dir, "outputs", "figure_data", "figure_3")

packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "sf", "ggplot2", "scales",
  "patchwork", "rnaturalearth", "systemfonts", "ggrepel"
)
missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(sf); library(ggplot2); library(scales); library(patchwork)
  library(rnaturalearth); library(ggrepel)
})

country_colours <- c(Argentina = "#D6604D", Brazil = "#477DB3")
driver_colours <- c(
  "VPD" = "#D55E00",
  "Temperature" = "#B2182B",
  "Reproductive water satisfaction" = "#2166AC",
  "Precipitation" = "#56B4E9",
  "Heavy rainfall" = "#5E3C99"
)
mechanism_colours <- c(
  "Thermal and atmospheric demand" = "#D55E00",
  "Crop water balance" = "#2166AC",
  "Water availability" = "#56B4E9",
  "Hot and dry compound" = "#B2182B",
  "Wet compound extremes" = "#5E3C99"
)
model_colours <- c(
  "ENSO" = "#777777",
  "Climate" = "#D55E00",
  "Water balance" = "#2166AC",
  "Climate and water balance" = "#6A51A3",
  "ENSO, climate and water balance" = "#111111"
)

available_fonts <- unique(systemfonts::system_fonts()$family)
figure_font_family <- case_when(
  "Times New Roman" %in% available_fonts ~ "Times New Roman",
  "Liberation Serif" %in% available_fonts ~ "Liberation Serif",
  TRUE ~ "serif"
)
message("Figure font: ", figure_font_family)

theme_rogerio <- function(base_size = 14) {
  theme_minimal(base_size = base_size, base_family = figure_font_family) +
    theme(
      panel.grid = element_blank(),
      axis.line.x = element_line(colour = "black", linewidth = 0.48),
      axis.line.y = element_line(colour = "black", linewidth = 0.48),
      axis.ticks = element_line(colour = "black", linewidth = 0.42),
      axis.ticks.length = grid::unit(1.8, "mm"),
      axis.text = element_text(colour = "black", size = rel(0.88)),
      axis.title = element_text(colour = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.box = "horizontal",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      panel.spacing = grid::unit(0.8, "lines"),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(5, 7, 5, 7),
      plot.tag = element_text(face = "bold", size = rel(1.12)),
      plot.tag.position = c(0.01, 0.99)
    )
}

theme_map <- function(base_size = 11.5) {
  theme_rogerio(base_size) +
    theme(
      legend.key.width = grid::unit(7, "mm"),
      legend.key.height = grid::unit(3.5, "mm")
    )
}

find_one <- function(root, candidates) {
  direct <- file.path(root, candidates)
  if (any(file.exists(direct))) return(direct[file.exists(direct)][[1]])
  all_files <- list.files(root, recursive = TRUE, full.names = TRUE)
  hits <- all_files[basename(all_files) %in% candidates]
  if (length(hits) == 0L) {
    stop("Missing required file: ", paste(candidates, collapse = " or "),
         call. = FALSE)
  }
  hits[order(nchar(hits), hits)][[1]]
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

paths <- list(
  map = find_one(figure_data_dir, "Figure_3_local_dominant_driver_map.csv"),
  mechanisms = find_one(figure_data_dir, "Figure_3_country_mechanism_summary.csv"),
  enso_climate = find_one(figure_data_dir, "Figure_3_ENSO_to_climate_evidence.csv"),
  models = find_one(figure_data_dir, "Figure_3_macroregion_model_performance.csv")
)
geometry_file <- find_one(project_dir, c(
  "ENSO_phase_specific_spatial_layers.gpkg",
  "ENSO_phase_specific_spatial_layers(1).gpkg"
))

map_data <- read_csv(
  paths$map,
  col_types = cols(
    spatial_unit_id = col_character(),
    spatial_unit_id_join = col_character(), .default = col_guess()
  ), show_col_types = FALSE
)
mechanisms <- read_csv(paths$mechanisms, show_col_types = FALSE)
enso_climate <- read_csv(paths$enso_climate, show_col_types = FALSE)
models <- read_csv(paths$models, show_col_types = FALSE)

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
  select(country, spatial_unit_id_join, geometry)

map_sf <- map_data |>
  mutate(
    spatial_unit_id_join = standardize_id(spatial_unit_id_join, country),
    evidence_class = factor(
      evidence_class,
      levels = c("Uncertain", "Moderate", "Robust", "Very robust")
    )
  ) |>
  inner_join(unit_geometry, by = c("country", "spatial_unit_id_join")) |>
  st_as_sf()

if (nrow(map_sf) == 0L) {
  stop("The spatial unit identifiers did not join to the geometry.", call. = FALSE)
}

south_america <- rnaturalearth::ne_countries(
  continent = "South America", scale = "medium", returnclass = "sf"
) |>
  st_make_valid()
map_bbox <- st_bbox(c(xmin = -82, ymin = -56, xmax = -34, ymax = 13),
                    crs = st_crs(4326))

panel_a <- ggplot() +
  geom_sf(data = south_america, fill = "grey98", colour = "black", linewidth = 0.36) +
  geom_sf(
    data = map_sf,
    aes(fill = driver_label, alpha = evidence_class),
    colour = NA, linewidth = 0
  ) +
  geom_sf(data = south_america, fill = NA, colour = "black", linewidth = 0.42) +
  scale_fill_manual(values = driver_colours, drop = FALSE) +
  scale_alpha_manual(
    values = c("Uncertain" = 0.24, "Moderate" = 0.48,
               "Robust" = 0.76, "Very robust" = 0.92),
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE, order = 1),
    alpha = guide_legend(nrow = 1, order = 2,
                         override.aes = list(fill = "grey35"))
  ) +
  coord_sf(
    xlim = c(map_bbox[["xmin"]], map_bbox[["xmax"]]),
    ylim = c(map_bbox[["ymin"]], map_bbox[["ymax"]]),
    expand = FALSE
  ) +
  labs(x = "Longitude", y = "Latitude", tag = "a") +
  theme_map(11.2) +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.box = "vertical"
  )

mechanism_plot <- mechanisms |>
  mutate(
    country = factor(country, levels = c("Argentina", "Brazil")),
    mechanism_group = factor(
      mechanism_group,
      levels = rev(c(
        "Thermal and atmospheric demand", "Crop water balance",
        "Water availability", "Hot and dry compound", "Wet compound extremes"
      ))
    ),
    aligned_size = if_else(aligned_score_percent > 0, aligned_score_percent, 0)
  )

panel_b <- ggplot(
  mechanism_plot,
  aes(
    x = ENSO_climate_skill_percent,
    y = productivity_importance_percent,
    colour = mechanism_group,
    size = aligned_size
  )
) +
  geom_hline(yintercept = 0, colour = "grey68", linewidth = 0.42) +
  geom_vline(xintercept = 0, colour = "grey68", linewidth = 0.42) +
  geom_point(alpha = 0.84, stroke = 0) +
  ggrepel::geom_text_repel(
    aes(label = mechanism_group),
    family = figure_font_family, size = 3.05, colour = "black",
    box.padding = 0.28, point.padding = 0.18,
    min.segment.length = 0, segment.colour = "grey68",
    show.legend = FALSE, max.overlaps = Inf
  ) +
  facet_wrap(~ country, nrow = 1, scales = "free") +
  scale_colour_manual(values = mechanism_colours, guide = "none") +
  scale_size_area(
    max_size = 9,
    breaks = c(25, 50, 75),
    labels = label_percent(scale = 1),
    name = "Share of aligned score"
  ) +
  scale_x_continuous(
    labels = label_percent(scale = 1), breaks = pretty_breaks(5),
    expand = expansion(mult = c(0.10, 0.18))
  ) +
  scale_y_continuous(
    labels = label_percent(scale = 1), breaks = pretty_breaks(5),
    expand = expansion(mult = c(0.08, 0.14))
  ) +
  labs(
    x = "ENSO contribution to climate prediction",
    y = "Climate contribution to yield prediction",
    tag = "b"
  ) +
  theme_rogerio(10.8) +
  theme(
    legend.position = "top", legend.title = element_text(face = "bold"),
    panel.spacing.x = grid::unit(1.2, "lines")
  )

driver_order <- enso_climate |>
  group_by(driver_label) |>
  summarise(maximum = max(ENSO_error_reduction_percent, na.rm = TRUE),
            .groups = "drop") |>
  arrange(maximum) |>
  pull(driver_label)

enso_climate_plot <- enso_climate |>
  mutate(
    driver_label = factor(driver_label, levels = driver_order),
    evidence = if_else(
      passes_circular_shift_test & positive_crossvalidated_skill,
      "Supported", "Uncertain or unsupported"
    )
  )

panel_c <- ggplot(
  enso_climate_plot,
  aes(x = ENSO_error_reduction_percent, y = driver_label, colour = country)
) +
  geom_vline(xintercept = 0, colour = "grey62", linewidth = 0.45) +
  geom_segment(
    aes(x = 0, xend = ENSO_error_reduction_percent,
        y = driver_label, yend = driver_label),
    colour = "grey82", linewidth = 0.55
  ) +
  geom_point(aes(shape = evidence), size = 2.75, stroke = 0.85) +
  facet_wrap(~ country, nrow = 1) +
  scale_colour_manual(values = country_colours, guide = "none") +
  scale_shape_manual(values = c(
    "Supported" = 16,
    "Uncertain or unsupported" = 1
  )) +
  scale_x_continuous(
    labels = label_percent(scale = 1), breaks = pretty_breaks(5),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    x = "Out of sample ENSO contribution to climate prediction",
    y = NULL, tag = "c"
  ) +
  theme_rogerio(10.6) +
  theme(
    legend.position = "top", legend.justification = "left",
    panel.spacing.x = grid::unit(1.25, "lines")
  )

model_labels <- c(
  "ENSO" = "ENSO",
  "climate" = "Climate",
  "water_balance" = "Water balance",
  "climate_plus_water" = "Climate and water balance",
  "ENSO_plus_climate_water" = "ENSO, climate and water balance"
)
macroregion_levels <- c(
  "Argentina", "Southern Brazil", "Central Brazil",
  "MATOPIBA and northeastern frontier"
)
models_plot <- models |>
  filter(model %in% names(model_labels)) |>
  mutate(
    model_label = recode(model, !!!model_labels),
    model_label = factor(model_label, levels = unname(model_labels)),
    macroregion = recode(
      macroregion,
      "MATOPIBA and northeastern frontier" =
        "MATOPIBA and\nnortheastern frontier"
    ),
    macroregion = factor(
      macroregion,
      levels = rev(c(
        "Argentina", "Southern Brazil", "Central Brazil",
        "MATOPIBA and\nnortheastern frontier"
      ))
    )
  )

panel_d <- ggplot(
  models_plot,
  aes(
    x = weighted_error_reduction_percent,
    y = macroregion,
    colour = model_label,
    shape = model_label
  )
) +
  geom_vline(xintercept = 0, colour = "grey62", linewidth = 0.45) +
  geom_errorbarh(
    aes(xmin = lower_region_percent, xmax = upper_region_percent),
    position = position_dodge(width = 0.62), height = 0,
    linewidth = 0.48, alpha = 0.40
  ) +
  geom_point(position = position_dodge(width = 0.62), size = 2.7) +
  scale_colour_manual(values = model_colours, drop = FALSE) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8), drop = FALSE) +
  scale_x_continuous(
    labels = label_percent(scale = 1), breaks = pretty_breaks(6),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    x = "Out of sample error reduction relative to the local null",
    y = NULL, tag = "d"
  ) +
  theme_rogerio(10.8) +
  theme(legend.position = "top", legend.justification = "left")

f3 <- (
  (panel_a | panel_b) /
  panel_c /
  panel_d
) +
  plot_layout(heights = c(1.42, 1.02, 0.78), widths = c(0.88, 1.12),
              guides = "keep")

# Run this final line in RStudio to display Figure 3.
f3

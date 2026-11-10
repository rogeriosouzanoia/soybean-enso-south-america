#!/usr/bin/env Rscript

# Figure 2
# Spatial distribution, regional strength and scale dependence of the
# out of sample predictive contribution of ENSO to soybean yield.

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1]] else {
  override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
  if (nzchar(override)) override else default_project_dir
}
figure_data_dir <- if (length(args) >= 2L) args[[2]] else
  file.path(project_dir, "outputs", "figure_data", "figure_2")

packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "sf", "ggplot2", "scales",
  "patchwork", "rnaturalearth", "systemfonts"
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
  library(rnaturalearth)
})

country_colours <- c(Argentina = "#D6604D", Brazil = "#477DB3")
macroregion_colours <- c(
  "Argentina" = "#C95340",
  "Southern Brazil" = "#315F91",
  "Central Brazil" = "#5D91BD",
  "MATOPIBA and northeastern frontier" = "#8EB6D2"
)
macroregion_levels <- c(
  "Argentina", "Southern Brazil", "Central Brazil",
  "MATOPIBA and northeastern frontier"
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
      plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.05)),
      legend.title = element_text(face = "bold", size = rel(0.88)),
      legend.key.width = grid::unit(8.5, "mm"),
      legend.key.height = grid::unit(3.4, "mm")
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

required_files <- c(
  local = "Figure_2_local_ENSO_metrics.csv",
  relationship = "Figure_2_local_skill_consistency_relationship.csv",
  macroregion = "Figure_2_macroregion_field_summary.csv",
  scale = "Figure_2_country_scale_comparison.csv"
)
paths <- lapply(required_files, function(x) find_one(figure_data_dir, x))
geometry_file <- find_one(project_dir, c(
  "ENSO_phase_specific_spatial_layers.gpkg",
  "ENSO_phase_specific_spatial_layers(1).gpkg"
))

local_metrics <- read_csv(
  paths$local,
  col_types = cols(
    spatial_unit_id = col_character(),
    spatial_unit_id_join = col_character(), .default = col_guess()
  ), show_col_types = FALSE
)
relationship <- read_csv(
  paths$relationship,
  col_types = cols(
    spatial_unit_id = col_character(),
    spatial_unit_id_join = col_character(), .default = col_guess()
  ), show_col_types = FALSE
)
macroregion <- read_csv(paths$macroregion, show_col_types = FALSE)
scale_comparison <- read_csv(paths$scale, show_col_types = FALSE)

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

local_sf <- local_metrics |>
  mutate(spatial_unit_id_join = standardize_id(spatial_unit_id_join, country)) |>
  inner_join(unit_geometry, by = c("country", "spatial_unit_id_join")) |>
  st_as_sf()

if (nrow(local_sf) == 0L) {
  stop("The spatial unit identifiers did not join to the geometry.", call. = FALSE)
}

south_america <- rnaturalearth::ne_countries(
  continent = "South America", scale = "medium", returnclass = "sf"
) |>
  st_make_valid()
map_bbox <- st_bbox(c(xmin = -82, ymin = -56, xmax = -34, ymax = 13),
                    crs = st_crs(4326))

predictive_levels <- c(
  "No positive skill", "0 to 2%", "2 to 5%", "5 to 10%",
  "10 to 20%", "Above 20%"
)
signal_levels <- c(
  "Below 1%", "1 to 2%", "2 to 5%", "5 to 10%", "10 to 20%",
  "Above 20%"
)
predictive_colours <- c(
  "No positive skill" = "#F1F1F1", "0 to 2%" = "#E1EEF5",
  "2 to 5%" = "#B9D9E8", "5 to 10%" = "#77B5D5",
  "10 to 20%" = "#2F8FBF", "Above 20%" = "#0B559F"
)
signal_colours <- c(
  "Below 1%" = "#F2F0F7", "1 to 2%" = "#DADAEB",
  "2 to 5%" = "#BCBDDC", "5 to 10%" = "#9E9AC8",
  "10 to 20%" = "#756BB1", "Above 20%" = "#54278F"
)

local_sf <- local_sf |>
  mutate(
    predictive_skill_class = factor(predictive_skill_class,
                                    levels = predictive_levels),
    signal_energy_class = factor(signal_energy_class, levels = signal_levels)
  )

make_class_map <- function(variable, palette, breaks, panel_title,
                           legend_title, legend_rows = 2L) {
  ggplot() +
    geom_sf(data = south_america, fill = "grey98", colour = "black",
            linewidth = 0.36) +
    geom_sf(data = local_sf, aes(fill = .data[[variable]]),
            colour = NA, linewidth = 0) +
    geom_sf(data = south_america, fill = NA, colour = "black", linewidth = 0.42) +
    scale_fill_manual(values = palette, breaks = breaks, drop = FALSE) +
    guides(fill = guide_legend(
      title = legend_title, nrow = legend_rows, byrow = TRUE,
      title.position = "top", title.hjust = 0
    )) +
    coord_sf(
      xlim = c(map_bbox[["xmin"]], map_bbox[["xmax"]]),
      ylim = c(map_bbox[["ymin"]], map_bbox[["ymax"]]),
      expand = FALSE
    ) +
    labs(x = "Longitude", y = "Latitude", title = panel_title) +
    theme_map(10.8) +
    theme(legend.position = "top", legend.justification = "left")
}

panel_a <- make_class_map(
  "predictive_skill_class", predictive_colours, predictive_levels,
  panel_title = "Local predictive contribution",
  legend_title = "Out of sample error reduction"
) +
  labs(tag = "a")

panel_b <- make_class_map(
  "signal_energy_class", signal_colours, signal_levels,
  panel_title = "Local signal magnitude",
  legend_title = "ENSO signal energy"
) +
  labs(tag = "b")

macroregion_long <- bind_rows(
  macroregion |>
    transmute(
      macroregion,
      metric = "Predictive contribution",
      estimate = cv_r_squared_percent,
      lower = cv_r2_lower_percent,
      upper = cv_r2_upper_percent
    ),
  macroregion |>
    transmute(
      macroregion,
      metric = "Signal energy",
      estimate = signal_energy_percent,
      lower = signal_lower_percent,
      upper = signal_upper_percent
    )
) |>
  mutate(
    macroregion = factor(macroregion, levels = rev(macroregion_levels)),
    macroregion_label = recode(
      as.character(macroregion),
      "MATOPIBA and northeastern frontier" =
        "MATOPIBA and\nnortheastern frontier"
    ),
    metric = factor(metric, levels = c("Predictive contribution", "Signal energy"))
  )

metric_colours <- c(
  "Predictive contribution" = "#1B6CA8",
  "Signal energy" = "#7A5195"
)

panel_c <- ggplot(
  macroregion_long,
  aes(
    x = estimate,
    y = factor(macroregion_label, levels = c(
      "MATOPIBA and\nnortheastern frontier", "Central Brazil",
      "Southern Brazil", "Argentina"
    )),
    colour = metric, shape = metric
  )
) +
  geom_vline(xintercept = 0, colour = "grey58", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 position = position_dodge(width = 0.42), height = 0,
                 linewidth = 0.72) +
  geom_point(position = position_dodge(width = 0.42), size = 3.0) +
  scale_colour_manual(values = metric_colours) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_continuous(labels = label_percent(scale = 1),
                     breaks = pretty_breaks(5), expand = expansion(mult = c(0.03, 0.08))) +
  labs(x = "Contribution or signal magnitude", y = NULL, tag = "c") +
  theme_rogerio(11.5) +
  theme(legend.position = "top", legend.justification = "left")

relationship_plot <- relationship |>
  filter(
    is.finite(cv_enso_signal_energy_percent),
    is.finite(sign_consistency),
    is.finite(reference_area_ha), reference_area_ha > 0
  ) |>
  mutate(
    macroregion = factor(macroregion, levels = macroregion_levels),
    phase = factor(phase, levels = c("El Niño", "La Niña"))
  )

signal_cap <- quantile(
  relationship_plot$cv_enso_signal_energy_percent, 0.99, na.rm = TRUE
)
relationship_plot <- relationship_plot |>
  mutate(
    signal_display = pmin(cv_enso_signal_energy_percent, signal_cap),
    signal_group = cut(
      cv_enso_signal_energy_percent,
      breaks = c(-Inf, 1, 2, 5, 10, 20, Inf),
      labels = c("<1", "1 to 2", "2 to 5", "5 to 10", "10 to 20", ">20"),
      right = FALSE
    )
  )

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

relationship_summary <- relationship_plot |>
  group_by(phase, macroregion, signal_group) |>
  summarise(
    consistency = weighted_mean_safe(
      sign_consistency,
      reference_area_ha * pmax(evaluated_events, 1)
    ),
    productive_area_ha = sum(reference_area_ha, na.rm = TRUE),
    spatial_units = n(),
    .groups = "drop"
  ) |>
  filter(is.finite(consistency)) |>
  mutate(
    signal_group = factor(
      signal_group,
      levels = c("<1", "1 to 2", "2 to 5", "5 to 10", "10 to 20", ">20")
    )
  )

panel_e <- ggplot(
  relationship_summary,
  aes(
    x = signal_group, y = consistency,
    colour = macroregion, group = macroregion
  )
) +
  geom_hline(yintercept = 0.5, colour = "grey72", linewidth = 0.45) +
  geom_line(linewidth = 0.78, alpha = 0.62) +
  geom_point(
    aes(size = productive_area_ha), alpha = 0.78,
    fill = "white", shape = 21, stroke = 0.9
  ) +
  facet_wrap(~ phase, nrow = 1) +
  scale_colour_manual(values = macroregion_colours, drop = FALSE) +
  scale_size_area(max_size = 5.0, guide = "none") +
  scale_y_continuous(
    labels = label_percent(accuracy = 1), limits = c(0, 1),
    breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.03))
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    x = "Local ENSO signal energy (%)",
    y = "Area weighted sign consistency", tag = "e"
  ) +
  theme_rogerio(11) +
  theme(
    legend.position = "top", legend.justification = "left",
    panel.spacing.x = grid::unit(1.4, "lines")
  )

scale_long <- scale_comparison |>
  transmute(
    country,
    `Spatial field` = spatial_field_cv_r2_percent,
    `National aggregate` = national_cv_r2_percent,
    `Null 95th percentile` = null_95_percentile_cv_r2_percent,
    spatial_lower = spatial_field_lower_percent,
    spatial_upper = spatial_field_upper_percent,
    null_threshold = null_95_percentile_cv_r2_percent,
    randomization_p_value
  ) |>
  pivot_longer(
    cols = c(`Spatial field`, `National aggregate`, `Null 95th percentile`),
    names_to = "scale", values_to = "estimate"
  ) |>
  mutate(scale = factor(
    scale,
    levels = c("Spatial field", "National aggregate", "Null 95th percentile")
  ))

panel_d <- ggplot(scale_long, aes(x = estimate, y = country, colour = country)) +
  geom_vline(xintercept = 0, colour = "grey58", linewidth = 0.45) +
  geom_line(
    data = filter(scale_long, scale != "Null 95th percentile"),
    aes(group = country), colour = "grey76", linewidth = 0.75
  ) +
  geom_errorbarh(
    data = filter(scale_long, scale == "Spatial field"),
    aes(xmin = spatial_lower, xmax = spatial_upper), height = 0,
    linewidth = 0.7
  ) +
  geom_point(
    aes(shape = scale), size = 3.2, fill = "white", stroke = 1.05
  ) +
  geom_text(
    data = filter(scale_long, scale == "Spatial field"),
    aes(
      x = spatial_upper,
      label = paste0("  p = ", formatC(randomization_p_value, format = "f", digits = 2))
    ),
    hjust = 0, vjust = -0.75, size = 3.4, show.legend = FALSE,
    family = figure_font_family
  ) +
  scale_colour_manual(values = country_colours, guide = "none") +
  scale_shape_manual(values = c(
    "Spatial field" = 16, "National aggregate" = 18,
    "Null 95th percentile" = 4
  )) +
  scale_x_continuous(labels = label_percent(scale = 1),
                     breaks = pretty_breaks(5), expand = expansion(mult = c(0.04, 0.18))) +
  labs(x = "Out of sample predictive contribution", y = NULL, tag = "d") +
  guides(shape = guide_legend(override.aes = list(colour = "black"))) +
  theme_rogerio(11.5) +
  theme(legend.position = "top", legend.justification = "left")

f2 <- (
  (panel_a | panel_b) /
  (panel_c | panel_d) /
  panel_e
) +
  plot_layout(heights = c(1.45, 0.72, 0.88), guides = "keep")

# Run this final line in RStudio to display Figure 2.
f2

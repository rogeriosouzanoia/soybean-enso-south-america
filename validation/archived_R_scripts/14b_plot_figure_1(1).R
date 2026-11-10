#!/usr/bin/env Rscript

# Figure 1
# Mean ENSO responses, temporal consistency, intensity dependence and
# regional spatial cancellation in South American soybean yield.

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1]] else {
  override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
  if (nzchar(override)) override else default_project_dir
}

figure_data_dir <- if (length(args) >= 2L) args[[2]] else
  file.path(project_dir, "outputs", "figure_data", "figure_1")

packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "purrr", "sf", "ggplot2",
  "scales", "patchwork", "rnaturalearth", "systemfonts"
)
missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(purrr); library(sf); library(ggplot2); library(scales)
  library(patchwork); library(rnaturalearth)
})

bootstrap_replicates <- 2000L
bootstrap_seed <- 1402L

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

theme_map <- function(base_size = 12) {
  theme_rogerio(base_size) +
    theme(
      legend.key.width = grid::unit(10, "mm"),
      legend.key.height = grid::unit(3.2, "mm")
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

bootstrap_mean <- function(x, B = bootstrap_replicates) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) return(tibble(mean = NA_real_, lower = NA_real_, upper = NA_real_))
  if (n < 3L) return(tibble(mean = mean(x), lower = NA_real_, upper = NA_real_))
  draws <- replicate(B, mean(sample(x, n, replace = TRUE)))
  tibble(
    mean = mean(x),
    lower = unname(quantile(draws, 0.025)),
    upper = unname(quantile(draws, 0.975))
  )
}

phase_response_file <- find_one(figure_data_dir, "Figure_1_local_phase_response.csv")
consistency_file <- find_one(figure_data_dir, "Figure_1_crossvalidated_consistency.csv")
event_file <- find_one(figure_data_dir, "Figure_1_event_spatial_extent.csv")
cancellation_file <- find_one(figure_data_dir, "Figure_1_regional_cancellation.csv")
geometry_file <- find_one(project_dir, c(
  "ENSO_phase_specific_spatial_layers.gpkg",
  "ENSO_phase_specific_spatial_layers(1).gpkg"
))

phase_response <- read_csv(
  phase_response_file,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
)
consistency <- read_csv(
  consistency_file,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
)
event_extent <- read_csv(event_file, show_col_types = FALSE)
cancellation <- read_csv(cancellation_file, show_col_types = FALSE)

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

south_america <- rnaturalearth::ne_countries(
  continent = "South America", scale = "medium", returnclass = "sf"
) |>
  st_make_valid()

phase_sf <- phase_response |>
  mutate(spatial_unit_id_join = standardize_id(spatial_unit_id, country)) |>
  inner_join(unit_geometry, by = c("country", "spatial_unit_id_join")) |>
  st_as_sf()

consistency_sf <- consistency |>
  mutate(spatial_unit_id_join = standardize_id(spatial_unit_id, country)) |>
  inner_join(unit_geometry, by = c("country", "spatial_unit_id_join")) |>
  st_as_sf()

if (nrow(phase_sf) == 0L || nrow(consistency_sf) == 0L) {
  stop("The spatial unit identifiers did not join to the geometry.", call. = FALSE)
}

map_bbox <- st_bbox(c(xmin = -82, ymin = -56, xmax = -34, ymax = 13),
                    crs = st_crs(4326))
response_limit <- max(abs(quantile(
  phase_sf$mean_yield_anomaly_t_ha, c(0.01, 0.99), na.rm = TRUE
)))
response_limit <- ceiling(response_limit * 10) / 10

make_response_map <- function(phase_name, tag) {
  ggplot() +
    geom_sf(data = south_america, fill = "grey98", colour = "black", linewidth = 0.38) +
    geom_sf(
      data = filter(phase_sf, phase == phase_name),
      aes(fill = mean_yield_anomaly_t_ha), colour = NA, linewidth = 0
    ) +
    geom_sf(data = south_america, fill = NA, colour = "black", linewidth = 0.42) +
    scale_fill_gradient2(
      low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
      midpoint = 0, limits = c(-response_limit, response_limit),
      oob = squish,
      breaks = pretty(c(-response_limit, response_limit), n = 5),
      name = expression("Mean yield anomaly (t ha"^{-1}*")")
    ) +
    coord_sf(
      xlim = c(map_bbox[["xmin"]], map_bbox[["xmax"]]),
      ylim = c(map_bbox[["ymin"]], map_bbox[["ymax"]]),
      expand = FALSE
    ) +
    labs(title = phase_name, x = "Longitude", y = "Latitude", tag = tag) +
    theme_map(11.5) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

make_consistency_map <- function(phase_name, tag) {
  ggplot() +
    geom_sf(data = south_america, fill = "grey98", colour = "black", linewidth = 0.38) +
    geom_sf(
      data = filter(consistency_sf, phase == phase_name),
      aes(fill = 100 * sign_consistency), colour = NA, linewidth = 0
    ) +
    geom_sf(data = south_america, fill = NA, colour = "black", linewidth = 0.42) +
    scale_fill_gradientn(
      colours = c("#F2F2F2", "#D9E8F2", "#85B6D1", "#2B78A6", "#084C7F"),
      values = rescale(c(0, 40, 55, 70, 100)), limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100), labels = label_percent(scale = 1),
      oob = squish, name = "Out of sample sign consistency"
    ) +
    coord_sf(
      xlim = c(map_bbox[["xmin"]], map_bbox[["xmax"]]),
      ylim = c(map_bbox[["ymin"]], map_bbox[["ymax"]]),
      expand = FALSE
    ) +
    labs(title = phase_name, x = "Longitude", y = "Latitude", tag = tag) +
    theme_map(11.5) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

map_el_nino <- make_response_map("El Niño", "a")
map_la_nina <- make_response_map("La Niña", "b")
consistency_el_nino <- make_consistency_map("El Niño", "c")
consistency_la_nina <- make_consistency_map("La Niña", "d")

# Strong and Very strong events are pooled. This avoids presenting a separate
# uncertainty interval for the single Very strong La Niña event.
event_plot_data <- event_extent |>
  mutate(
    intensity_plot = case_when(
      str_to_lower(intensity) %in% c("strong", "very strong") ~ "Strong",
      str_to_lower(intensity) == "moderate" ~ "Moderate",
      str_to_lower(intensity) %in% c("weak", "very weak") ~ "Weak",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(intensity_plot))

set.seed(bootstrap_seed)
gain_summary <- event_plot_data |>
  group_by(macroregion, phase, intensity_plot) |>
  group_modify(~ bootstrap_mean(.x$gain_area_percent) |>
                 mutate(n_events = n_distinct(.x$event_id))) |>
  ungroup() |>
  mutate(
    macroregion_plot = recode(
      macroregion,
      "MATOPIBA and northeastern frontier" = "MATOPIBA and Northeast"
    ),
    phase = factor(phase, levels = c("El Niño", "La Niña")),
    intensity_plot = factor(intensity_plot, levels = c("Weak", "Moderate", "Strong")),
    category = paste0(intensity_plot, " ", phase, " (n = ", n_events, ")"),
    category_order = as.integer(phase) * 10 + as.integer(intensity_plot)
  ) |>
  arrange(macroregion, desc(category_order)) |>
  mutate(category = factor(category, levels = unique(category)))

intensity_bars <- bind_rows(
  gain_summary |> transmute(macroregion_plot, category, direction = "Gain", value = mean),
  gain_summary |> transmute(macroregion_plot, category, direction = "Loss", value = 100 - mean)
) |>
  mutate(direction = factor(direction, levels = c("Loss", "Gain")))

intensity_plot <- ggplot(
  intensity_bars,
  aes(x = value, y = category, fill = direction)
) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25, alpha = 0.58) +
  geom_segment(
    data = gain_summary,
    aes(x = lower, xend = upper, y = category, yend = category),
    inherit.aes = FALSE, colour = "grey25", linewidth = 0.48, alpha = 0.65
  ) +
  geom_point(
    data = gain_summary,
    aes(x = mean, y = category), inherit.aes = FALSE,
    shape = 21, size = 1.8, stroke = 0.45, fill = "white",
    colour = "grey20", alpha = 0.78
  ) +
  facet_wrap(~ macroregion_plot, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("Loss" = "#C95D4B", "Gain" = "#3B78A8")) +
  scale_x_continuous(
    limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
    labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = "Productive area with gain or loss",
    y = NULL,
    tag = "e"
  ) +
  theme_rogerio(10.5) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 8.2),
    strip.text = element_text(size = 10.5, face = "bold")
  )

cancellation_plot <- cancellation |>
  mutate(
    macroregion = factor(
      macroregion,
      levels = rev(c(
        "Argentina", "Southern Brazil", "Central Brazil",
        "MATOPIBA and northeastern frontier"
      ))
    ),
    phase = factor(phase, levels = c("El Niño", "La Niña"))
  ) |>
  ggplot(aes(x = mean, y = macroregion, colour = phase)) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    position = position_dodge(width = 0.42), width = 0,
    linewidth = 0.7, orientation = "y", alpha = 0.62
  ) +
  geom_point(position = position_dodge(width = 0.42), size = 2.8, alpha = 0.76) +
  scale_colour_manual(values = c("El Niño" = "#C73E4D", "La Niña" = "#2166AC")) +
  scale_x_continuous(
    limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
    labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Spatial cancellation", y = NULL, tag = "f") +
  theme_rogerio(10.5) +
  theme(axis.text.y = element_text(size = 8.8))

mean_maps <- (map_el_nino + map_la_nina) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

consistency_maps <- (consistency_el_nino + consistency_la_nina) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

bottom_row <- intensity_plot + cancellation_plot +
  plot_layout(widths = c(1.58, 1))

f1 <- mean_maps / consistency_maps / bottom_row +
  plot_layout(heights = c(1.05, 1.05, 1.12))

# Run this final line in RStudio to display Figure 1.
f1

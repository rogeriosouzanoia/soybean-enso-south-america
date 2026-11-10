#!/usr/bin/env Rscript

# Figure 5
# Spatial reliability, predictability and cancellation of ENSO yield effects.

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1]] else {
  override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
  if (nzchar(override)) override else default_project_dir
}
figure_data_dir <- if (length(args) >= 2L) args[[2]] else
  file.path(project_dir, "outputs", "figure_data", "figure_5")

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

available_fonts <- unique(systemfonts::system_fonts()$family)
figure_font_family <- dplyr::case_when(
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
      legend.title = element_text(face = "bold", size = rel(0.9)),
      legend.text = element_text(size = rel(0.84)),
      legend.box = "vertical",
      legend.box.just = "center",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      panel.spacing = grid::unit(0.8, "lines"),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(6, 8, 6, 8),
      plot.tag = element_text(face = "bold", size = rel(1.12)),
      plot.tag.position = c(0.01, 0.99)
    )
}

theme_map <- function(base_size = 11.5) {
  theme_rogerio(base_size) +
    theme(
      legend.key.width = grid::unit(8, "mm"),
      legend.key.height = grid::unit(3.5, "mm"),
      plot.title = element_text(face = "bold", hjust = 0.5,
                                size = rel(1.12), margin = margin(b = 4))
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

local_file <- find_one(
  figure_data_dir, "Figure_5_local_ENSO_evidence_classification.csv"
)
area_file <- find_one(
  figure_data_dir, "Figure_5_regional_evidence_area_summary.csv"
)
cancellation_file <- find_one(
  figure_data_dir, "Figure_5_regional_cancellation.csv"
)
geometry_file <- find_one(project_dir, c(
  "ENSO_phase_specific_spatial_layers.gpkg",
  "ENSO_phase_specific_spatial_layers(1).gpkg"
))

local_evidence <- read_csv(
  local_file,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
)
area_summary <- read_csv(area_file, show_col_types = FALSE)
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

class_order <- c(
  "Consistent and predictable",
  "Predictable but variable",
  "Consistent but low predictive skill",
  "Uncertain local response"
)
class_labels <- c(
  "Consistent and predictable" = "Consistent and predictable",
  "Predictable but variable" = "Predictable, variable sign",
  "Consistent but low predictive skill" = "Consistent, low predictive skill",
  "Uncertain local response" = "Uncertain response"
)
class_alpha <- c(
  "Consistent and predictable" = 0.96,
  "Predictable but variable" = 0.58,
  "Consistent but low predictive skill" = 0.40,
  "Uncertain local response" = 0.16
)
class_colours <- c(
  "Consistent and predictable" = "#2F5D50",
  "Predictable but variable" = "#78909C",
  "Consistent but low predictive skill" = "#C8A24B",
  "Uncertain local response" = "#D9D9D9"
)
sign_colours <- c("Loss" = "#C64E43", "Gain" = "#3F78A8", "No change" = "#BDBDBD")
phase_colours <- c("El Nino" = "#C64E43", "La Nina" = "#3F78A8")
phase_labels <- c("El Nino" = "El Niño", "La Nina" = "La Niña")
macroregion_order <- c(
  "Argentina", "Southern Brazil", "Central Brazil",
  "MATOPIBA and northeastern frontier"
)

local_sf <- local_evidence |>
  mutate(
    spatial_unit_id_join = standardize_id(spatial_unit_id, country),
    evidence_class = factor(evidence_class, levels = class_order)
  ) |>
  inner_join(unit_geometry, by = c("country", "spatial_unit_id_join")) |>
  st_as_sf()

if (nrow(local_sf) != nrow(local_evidence)) {
  warning("Some local evidence rows did not match the spatial geometry.")
}

map_bbox <- st_bbox(c(xmin = -82, ymin = -56, xmax = -34, ymax = 13),
                    crs = st_crs(4326))

make_evidence_map <- function(phase_code, phase_title) {
  selected <- local_sf |> filter(as.character(ENSO_phase) == phase_code)
  ggplot() +
    geom_sf(
      data = south_america, fill = "grey98", colour = "black", linewidth = 0.36
    ) +
    geom_sf(
      data = selected,
      aes(fill = response_sign, alpha = evidence_class),
      colour = NA, linewidth = 0
    ) +
    geom_sf(data = south_america, fill = NA, colour = "black", linewidth = 0.42) +
    scale_fill_manual(
      values = sign_colours, breaks = c("Loss", "Gain"),
      name = "Mean yield response"
    ) +
    scale_alpha_manual(
      values = class_alpha, breaks = class_order, labels = class_labels,
      name = "Evidence strength",
      guide = guide_legend(override.aes = list(fill = "#333333"), nrow = 2)
    ) +
    coord_sf(
      xlim = c(map_bbox[["xmin"]], map_bbox[["xmax"]]),
      ylim = c(map_bbox[["ymin"]], map_bbox[["ymax"]]),
      expand = FALSE
    ) +
    labs(title = phase_title, x = "Longitude", y = "Latitude") +
    theme_map() +
    theme(legend.position = "top")
}

panel_a <- make_evidence_map("El Nino", "El Niño")
panel_b <- make_evidence_map("La Nina", "La Niña")

area_plot_data <- area_summary |>
  mutate(
    macroregion = factor(macroregion, levels = rev(macroregion_order)),
    ENSO_phase = factor(ENSO_phase, levels = c("El Nino", "La Nina")),
    evidence_class = factor(evidence_class, levels = rev(class_order))
  ) |>
  complete(
    macroregion, ENSO_phase, evidence_class,
    fill = list(productive_area_percent = 0)
  )

panel_c <- ggplot(
  area_plot_data,
  aes(x = productive_area_percent, y = macroregion, fill = evidence_class)
) +
  geom_col(width = 0.64, colour = "white", linewidth = 0.18) +
  facet_wrap(
    vars(ENSO_phase), ncol = 1,
    labeller = as_labeller(phase_labels)
  ) +
  scale_x_continuous(
    limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
    labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.01))
  ) +
  scale_fill_manual(
    values = class_colours, breaks = class_order, labels = class_labels,
    name = "Local evidence class"
  ) +
  labs(x = "Share of soybean productive area", y = NULL) +
  theme_rogerio(12.5) +
  theme(
    legend.position = "top",
    legend.key.size = grid::unit(3.8, "mm"),
    strip.text = element_text(hjust = 0)
  )

cancellation_plot_data <- cancellation |>
  mutate(
    macroregion = factor(macroregion, levels = rev(macroregion_order)),
    ENSO_phase = factor(ENSO_phase, levels = c("El Nino", "La Nina"))
  )

panel_d <- ggplot(
  cancellation_plot_data,
  aes(
    x = cancellation_mean_percent, y = macroregion,
    xmin = cancellation_lower_percent, xmax = cancellation_upper_percent,
    colour = ENSO_phase
  )
) +
  geom_vline(xintercept = 0, colour = "#BFBFBF", linewidth = 0.45) +
  geom_errorbarh(
    aes(y = macroregion), height = 0, linewidth = 0.85,
    position = position_dodge(width = 0.42)
  ) +
  geom_point(size = 3.3, position = position_dodge(width = 0.42)) +
  scale_x_continuous(
    limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
    labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    values = phase_colours, labels = phase_labels,
    name = "ENSO phase"
  ) +
  labs(x = "Spatial cancellation of local responses", y = NULL) +
  theme_rogerio(12.5) +
  theme(legend.position = "top")

top_row <- (panel_a | panel_b) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

bottom_row <- panel_c | panel_d

f5 <- (top_row / bottom_row) +
  plot_layout(heights = c(1.58, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.justification = "center",
    plot.tag = element_text(
      face = "bold", family = figure_font_family, size = 15
    )
  )

# Run this final line in RStudio to display Figure 5.
f5

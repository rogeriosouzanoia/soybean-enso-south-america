# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 13_synthesize_ENSO_soybean_figures.R
# PURPOSE: Build a result driven figure sequence from local ENSO responses to
#          spatial cancellation and climate mechanisms
# VERSION: 3.2
# ==============================================================================

# 1. SETTINGS ------------------------------------------------------------------

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(override)) override else default_project_dir

figure_dpi <- 600
phase_limit <- 15

# 2. PACKAGES ------------------------------------------------------------------

packages <- c(
  "dplyr", "readr", "tidyr", "stringr", "sf", "ggplot2", "scales",
  "rnaturalearth", "patchwork", "ragg", "systemfonts"
)
missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
  library(sf); library(ggplot2); library(scales)
  library(rnaturalearth); library(patchwork)
})

installed_font_families <- unique(systemfonts::system_fonts()$family)
figure_font_family <- if ("Calibri" %in% installed_font_families) {
  "Calibri"
} else if ("Arial" %in% installed_font_families) {
  "Arial"
} else {
  "sans"
}

message("Figure font: ", figure_font_family)

# 3. PATHS ---------------------------------------------------------------------

field_dir <- file.path(project_dir, "enso_multivariate_field_attribution")
local_dir <- file.path(project_dir, "enso_spatial_attribution")
climate_dir <- file.path(project_dir, "climate_mechanisms")
attribution_dir <- file.path(climate_dir, "spatial_attribution")
map_dir <- file.path(project_dir, "figures", "enso_spatial_responses")
output_dir <- file.path(project_dir, "figures", "final_story")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  geometry = file.path(map_dir, "ENSO_phase_specific_spatial_layers.gpkg"),
  phase_local = file.path(map_dir, "ENSO_local_phase_specific_effects.csv"),
  phase_region = file.path(map_dir, "ENSO_state_province_phase_effects.csv"),
  local_enso = file.path(local_dir, "ENSO_spatially_resolved_local_summary.csv"),
  field_country = file.path(field_dir, "ENSO_multivariate_field_country_summary.csv"),
  field_region = file.path(field_dir, "ENSO_multivariate_field_region_summary.csv"),
  climate_country = file.path(attribution_dir, "ENSO_climate_attribution_country_summary.csv"),
  climate_region = file.path(attribution_dir, "ENSO_climate_attribution_region_summary.csv"),
  mechanism_importance = file.path(attribution_dir, "ENSO_climate_variable_permutation_importance.csv"),
  enso_mechanism = file.path(attribution_dir, "ENSO_to_mechanism_group_crossvalidated_skill.csv"),
  compounds = file.path(climate_dir, "ENSO_compound_event_frequency_by_phase.csv")
)
missing_files <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing_files) > 0) {
  stop("Missing required outputs: ", paste(missing_files, collapse = "; "),
       call. = FALSE)
}

# 4. THEME AND HELPERS ---------------------------------------------------------

theme_rogerio <- function(base_size = 14) {
  theme_classic(base_size = base_size, base_family = figure_font_family) +
    theme(
      axis.line.x = element_line(colour = "black", linewidth = 0.48),
      axis.line.y = element_line(colour = "black", linewidth = 0.48),
      axis.ticks = element_line(colour = "black", linewidth = 0.42),
      axis.ticks.length = grid::unit(1.7, "mm"),
      axis.text = element_text(colour = "black", size = rel(0.88)),
      axis.title = element_text(colour = "black"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "top", legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black", size = rel(0.95)),
      panel.spacing = grid::unit(1.15, "lines"),
      plot.margin = margin(7, 8, 7, 8),
      plot.tag = element_text(face = "bold", size = rel(1.15)),
      plot.tag.position = c(0.01, 0.99)
    )
}

theme_map <- function(base_size = 14) {
  theme_rogerio(base_size) +
    theme(
      legend.key.height = grid::unit(3.0, "mm"),
      legend.key.width = grid::unit(9.0, "mm")
    )
}

tag_panel <- function(plot, tag) {
  plot + labs(tag = tag) +
    theme(plot.tag = element_text(face = "bold", colour = "black"))
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

standardize_id <- function(x, country) {
  digits <- gsub("[^0-9]", "", as.character(x))
  widths <- ifelse(country == "Brazil", 7L, 5L)
  mapply(
    function(value, width) {
      if (is.na(value) || value == "") return(NA_character_)
      sprintf(paste0("%0", width, "d"), as.integer(value))
    }, digits, widths, USE.NAMES = FALSE
  )
}

save_figure <- function(plot, name, width = 12, height = 8.2) {
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "in", bg = "white",
         device = grDevices::cairo_pdf)
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = figure_dpi,
         bg = "white", device = ragg::agg_png)
  ggsave(file.path(output_dir, paste0(name, ".tiff")), plot,
         width = width, height = height, units = "in", dpi = figure_dpi,
         compression = "lzw", bg = "white", device = ragg::agg_tiff)
}

phase_label <- function(x) {
  ifelse(grepl("El", x, ignore.case = TRUE),
         "El Niño versus neutral", "La Niña versus neutral")
}

mechanism_labels <- c(
  crop_water_balance = "Crop water balance",
  water_availability = "Water availability",
  thermal_atmospheric_demand = "Thermal and atmospheric demand",
  hot_dry_compound = "Hot and dry conditions",
  wet_compound_extremes = "Wet compound extremes"
)

# 5. DATA ----------------------------------------------------------------------

phase_local <- read_csv(
  paths$phase_local,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
)
phase_region <- read_csv(paths$phase_region, show_col_types = FALSE)
local_enso <- read_csv(
  paths$local_enso,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
)
field_country <- read_csv(paths$field_country, show_col_types = FALSE)
field_region <- read_csv(paths$field_region, show_col_types = FALSE)
climate_country <- read_csv(paths$climate_country, show_col_types = FALSE)
climate_region <- read_csv(paths$climate_region, show_col_types = FALSE)
importance <- read_csv(paths$mechanism_importance, show_col_types = FALSE)
enso_mechanism <- read_csv(paths$enso_mechanism, show_col_types = FALSE)
compounds <- read_csv(paths$compounds, show_col_types = FALSE)

geometry_layer <- st_read(
  paths$geometry, layer = "local_phase_effects", quiet = TRUE
) |>
  standardize_geometry() |>
  st_make_valid() |>
  mutate(
    spatial_unit_id_join = standardize_id(spatial_unit_id, country)
  )

unit_geometry <- geometry_layer |>
  arrange(country, spatial_unit_id_join) |>
  group_by(country, spatial_unit_id_join) |>
  slice(1L) |>
  ungroup() |>
  select(country, spatial_unit_id_join, geometry)

# Do not derive state or province outlines from the local polygons. Small
# topological gaps in the municipal source recreate thousands of internal
# boundaries after union and visually contaminate high resolution maps.
# Clean administrative level 1 outlines can be added later from a dedicated
# state and province geometry, but local polygons must never supply this layer.

south_america <- rnaturalearth::ne_countries(
  continent = "South America", scale = "medium", returnclass = "sf"
) |>
  st_make_valid()

phase_map <- phase_local |>
  mutate(
    spatial_unit_id_join = standardize_id(spatial_unit_id, country),
    phase_panel = factor(
      phase_label(enso_phase),
      levels = c("El Niño versus neutral", "La Niña versus neutral")
    )
  ) |>
  inner_join(unit_geometry, by = c("country", "spatial_unit_id_join")) |>
  st_as_sf()

phase_levels <- levels(phase_map$phase_panel)
sa_phase <- bind_rows(lapply(phase_levels, function(label) {
  south_america |> mutate(phase_panel = factor(label, levels = phase_levels))
}))

# 6. FIGURE 1: THE CENTRAL PARADOX --------------------------------------------

map_panel <- ggplot() +
  geom_sf(data = sa_phase, fill = "grey98", colour = "black", linewidth = 0.42) +
  # A hairline in the same colour as the fill prevents antialiasing seams
  # between adjacent municipalities and departments in high resolution output.
  geom_sf(
    data = phase_map,
    aes(fill = effect, colour = effect),
    linewidth = 0.075,
    linejoin = "mitre"
  ) +
  geom_sf(data = sa_phase, fill = NA, colour = "black", linewidth = 0.48) +
  facet_wrap(~ phase_panel, nrow = 1, axes = "all", axis.labels = "all") +
  scale_fill_gradient2(
    low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0,
    limits = c(-phase_limit, phase_limit), oob = squish,
    breaks = c(-15, -7.5, 0, 7.5, 15),
    labels = c("−15", "−7.5", "0", "7.5", "15"), na.value = "transparent"
  ) +
  scale_colour_gradient2(
    low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0,
    limits = c(-phase_limit, phase_limit), oob = squish,
    na.value = "transparent", guide = "none"
  ) +
  coord_sf(xlim = c(-82, -32), ylim = c(-56, 13), expand = FALSE) +
  guides(fill = guide_colourbar(barwidth = grid::unit(6.5, "cm"))) +
  labs(x = "Longitude", y = "Latitude") +
  theme_map(12) +
  theme(legend.position = "top")

local_distribution <- phase_local |>
  filter(is.finite(effect)) |>
  mutate(
    phase_panel = factor(
      phase_label(enso_phase),
      levels = c("El Niño versus neutral", "La Niña versus neutral")
    ),
    country = factor(country, levels = c("Brazil", "Argentina"))
  )

distribution_panel <- ggplot(
  local_distribution,
  aes(x = effect, y = country, fill = country)
) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.75,
              colour = NA, width = 0.85) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.2,
               fill = "white", colour = "black") +
  facet_wrap(~ phase_panel, nrow = 1, axes = "all", axis.labels = "all") +
  coord_cartesian(xlim = c(-30, 30)) +
  scale_fill_manual(values = c("Brazil" = "#4C78A8", "Argentina" = "#D6604D")) +
  labs(x = "Local productivity response, percentage points", y = NULL) +
  theme_rogerio(11) +
  theme(legend.position = "none")

cancellation_panel <- field_country |>
  select(
    country, spatial_field_cv_r_squared_percent,
    national_net_cv_r_squared_percent, cancellation_index
  ) |>
  pivot_longer(
    c(spatial_field_cv_r_squared_percent, national_net_cv_r_squared_percent),
    names_to = "scale", values_to = "skill"
  ) |>
  mutate(
    scale = recode(
      scale,
      spatial_field_cv_r_squared_percent = "Spatial field",
      national_net_cv_r_squared_percent = "National aggregate"
    )
  ) |>
  ggplot(aes(x = scale, y = skill, fill = country)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(
    aes(label = sprintf("%.1f%%", skill)),
    position = position_dodge(width = 0.72), vjust = -0.4,
    family = figure_font_family, size = 3.1
  ) +
  scale_fill_manual(values = c("Brazil" = "#4C78A8", "Argentina" = "#D6604D")) +
  scale_y_continuous(labels = label_number(suffix = "%", accuracy = 1),
                     expand = expansion(mult = c(0.08, 0.15))) +
  labs(x = NULL, y = "Cross validated predictive skill") +
  theme_rogerio(11) +
  theme(legend.position = "top")

figure_1 <- tag_panel(map_panel, "a") /
  (tag_panel(distribution_panel, "b") | tag_panel(cancellation_panel, "c")) +
  plot_layout(heights = c(2.85, 1)) &
  theme(plot.margin = margin(7, 7, 7, 7))

save_figure(figure_1, "Figure_1_spatial_contrast_and_national_cancellation", 13, 9.2)

# 7. FIGURE 2: HOW MUCH ENSO INFORMATION IS RETAINED --------------------------

local_skill_map <- local_enso |>
  mutate(
    spatial_unit_id_join = standardize_id(spatial_unit_id, country),
    skill = pmax(0, cv_r_squared_percent),
    skill_class = cut(
      skill, c(-Inf, 0, 2, 5, 10, 20, Inf),
      labels = c("No positive skill", "0 to 2%", "2 to 5%", "5 to 10%",
                 "10 to 20%", "Above 20%"), include.lowest = TRUE
    )
  ) |>
  inner_join(unit_geometry, by = c("country", "spatial_unit_id_join")) |>
  st_as_sf()

map_skill <- ggplot() +
  geom_sf(data = south_america, fill = "grey98", colour = "black", linewidth = 0.42) +
  # Matching fill and stroke removes polygon seams without showing local borders.
  geom_sf(
    data = local_skill_map,
    aes(fill = skill_class, colour = skill_class),
    linewidth = 0.075,
    linejoin = "mitre"
  ) +
  geom_sf(data = south_america, fill = NA, colour = "black", linewidth = 0.48) +
  scale_fill_manual(
    values = c(
      "No positive skill" = "#F2F2F2", "0 to 2%" = "#E1EEF6",
      "2 to 5%" = "#BDDBEA", "5 to 10%" = "#74B3D3",
      "10 to 20%" = "#2B8CBE", "Above 20%" = "#08589E"
    ), drop = FALSE, na.value = "transparent"
  ) +
  scale_colour_manual(
    values = c(
      "No positive skill" = "#F2F2F2", "0 to 2%" = "#E1EEF6",
      "2 to 5%" = "#BDDBEA", "5 to 10%" = "#74B3D3",
      "10 to 20%" = "#2B8CBE", "Above 20%" = "#08589E"
    ), drop = FALSE, na.value = "transparent", guide = "none"
  ) +
  coord_sf(xlim = c(-90, -32), ylim = c(-56, 13), expand = FALSE) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_map(11)

country_metrics <- field_country |>
  transmute(
    country,
    `Spatial predictive skill` = spatial_field_cv_r_squared_percent,
    `ENSO signal energy` = spatial_signal_energy_percent,
    `National predictive skill` = national_net_cv_r_squared_percent,
    `Spatial cancellation` = 100 * cancellation_index
  ) |>
  pivot_longer(-country, names_to = "metric", values_to = "value") |>
  mutate(
    metric = factor(
      metric,
      levels = c("Spatial predictive skill", "ENSO signal energy",
                 "National predictive skill", "Spatial cancellation")
    )
  )

metric_panel <- ggplot(country_metrics, aes(x = value, y = metric, colour = country)) +
  geom_segment(aes(x = 0, xend = value, yend = metric),
               colour = "grey78", linewidth = 0.65) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.1f%%", value)), nudge_y = 0.18,
            family = figure_font_family, size = 3.1) +
  facet_wrap(~ country, ncol = 1) +
  scale_colour_manual(values = c("Brazil" = "#4C78A8", "Argentina" = "#D6604D")) +
  scale_x_continuous(labels = label_number(suffix = "%", accuracy = 1)) +
  labs(x = "Magnitude", y = NULL) +
  theme_rogerio(11) +
  theme(legend.position = "none")

metric_argentina <- metric_panel %+% filter(country_metrics, country == "Argentina") +
  labs(title = "Argentina") + theme(strip.text = element_blank())
metric_brazil <- metric_panel %+% filter(country_metrics, country == "Brazil") +
  labs(title = "Brazil") + theme(strip.text = element_blank())

figure_2 <- tag_panel(map_skill, "a") |
  (tag_panel(metric_argentina, "b") / tag_panel(metric_brazil, "c")) +
  plot_layout(widths = c(2.1, 1))

save_figure(figure_2, "Figure_2_ENSO_information_and_spatial_cancellation", 13, 8.2)

# 8. FIGURE 3: FROM ENSO TO CLIMATE TO PRODUCTIVITY ---------------------------

pathway_data <- importance |>
  select(country, mechanism_group, normalized_positive_importance) |>
  left_join(
    enso_mechanism |>
      select(country, mechanism_group, mean_positive_ENSO_error_reduction),
    by = c("country", "mechanism_group")
  ) |>
  mutate(
    mechanism = recode(mechanism_group, !!!mechanism_labels),
    `ENSO to climate` = 100 * mean_positive_ENSO_error_reduction,
    `Climate to productivity` = 100 * normalized_positive_importance
  ) |>
  select(country, mechanism, `ENSO to climate`, `Climate to productivity`) |>
  pivot_longer(c(`ENSO to climate`, `Climate to productivity`),
               names_to = "link", values_to = "strength") |>
  mutate(
    mechanism = factor(
      mechanism,
      levels = rev(c(
        "Thermal and atmospheric demand", "Crop water balance",
        "Water availability", "Hot and dry conditions",
        "Wet compound extremes"
      ))
    )
  )

make_pathway_panel <- function(link_name, x_limits) {
  point_data <- pathway_data |> filter(link == link_name)
  range_data <- point_data |>
    group_by(mechanism) |>
    summarise(
      xmin = min(strength, na.rm = TRUE),
      xmax = max(strength, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot(point_data, aes(x = strength, y = mechanism, colour = country)) +
    geom_segment(
      data = range_data,
      aes(x = xmin, xend = xmax, y = mechanism, yend = mechanism),
      inherit.aes = FALSE, colour = "grey78", linewidth = 0.7
    ) +
    geom_point(size = 2.9) +
    scale_colour_manual(
      values = c("Brazil" = "#4C78A8", "Argentina" = "#D6604D")
    ) +
    scale_x_continuous(labels = label_number(suffix = "%", accuracy = 1)) +
    coord_cartesian(xlim = x_limits) +
    labs(x = "Predictive contribution (%)", y = NULL, title = link_name) +
    theme_rogerio(11)
}

climate_to_yield_panel <- make_pathway_panel(
  "Climate to productivity", c(0, 90)
)

enso_to_climate_panel <- make_pathway_panel(
  "ENSO to climate", c(0, 15)
)

climate_total <- climate_country |>
  filter(model %in% c("ENSO", "climate", "water_balance", "climate_plus_water")) |>
  mutate(
    model = recode(
      model,
      ENSO = "ENSO only", climate = "Climate",
      water_balance = "Crop water balance",
      climate_plus_water = "Climate and water balance"
    ),
    explained = 100 * out_of_sample_r_squared
  )

total_panel <- ggplot(climate_total, aes(x = explained, y = model, colour = country)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.5), size = 2.8) +
  scale_colour_manual(values = c("Brazil" = "#4C78A8", "Argentina" = "#D6604D")) +
  scale_x_continuous(labels = label_number(suffix = "%", accuracy = 1)) +
  labs(x = "Local out of sample explained variability", y = NULL) +
  theme_rogerio(11)

figure_3 <- (tag_panel(climate_to_yield_panel, "a") |
             tag_panel(enso_to_climate_panel, "b")) /
  tag_panel(total_panel, "c") +
  plot_layout(heights = c(1.15, 0.82), guides = "collect") &
  theme(legend.position = "top")

save_figure(figure_3, "Figure_3_ENSO_climate_productivity_pathways", 12.5, 7.6)

# 9. FIGURE 4: REGIONAL SYNTHESIS ---------------------------------------------

region_story <- field_region |>
  select(country, parent_region_name,
         enso_field_skill = cv_r_squared_percent,
         enso_signal = signal_energy_percent) |>
  full_join(
    climate_region |>
      filter(model == "climate_plus_water") |>
      transmute(country, parent_region_name,
                climate_skill = 100 * out_of_sample_r_squared),
    by = c("country", "parent_region_name")
  ) |>
  left_join(
    phase_region |>
      transmute(
        country, parent_region_name,
        phase = phase_label(enso_phase), effect
      ) |>
      pivot_wider(names_from = phase, values_from = effect),
    by = c("country", "parent_region_name")
  ) |>
  rename(
    `ENSO field skill` = enso_field_skill,
    `ENSO signal magnitude` = enso_signal,
    `Climate predictive skill` = climate_skill,
    `El Niño response` = `El Niño versus neutral`,
    `La Niña response` = `La Niña versus neutral`
  ) |>
  pivot_longer(
    c(`ENSO field skill`, `ENSO signal magnitude`, `Climate predictive skill`,
      `El Niño response`, `La Niña response`),
    names_to = "metric", values_to = "value"
  ) |>
  mutate(
    region_id = paste(country, parent_region_name, sep = "___"),
    metric = factor(
      metric,
      levels = c("ENSO field skill", "ENSO signal magnitude",
                 "Climate predictive skill", "El Niño response", "La Niña response")
    ),
    metric_short = factor(
      metric,
      levels = levels(metric),
      labels = c(
        "ENSO field\nskill", "ENSO signal\nmagnitude",
        "Climate predictive\nskill", "El Niño\nresponse",
        "La Niña\nresponse"
      )
    )
  )

regional_matrix_base <- ggplot(region_story, aes(x = metric_short, y = region_id)) +
  geom_point(
    aes(size = pmin(abs(value), 30), fill = value),
    shape = 21, colour = "grey35", stroke = 0.25
  ) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
    limits = c(-30, 30), oob = squish
  ) +
  scale_size_continuous(range = c(1.8, 8), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_rogerio(10.5) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.95))

matrix_argentina <- regional_matrix_base %+%
  filter(region_story, country == "Argentina") +
  labs(title = "Argentina")
matrix_brazil <- regional_matrix_base %+%
  filter(region_story, country == "Brazil") +
  labs(title = "Brazil")

regional_phase_contrast <- phase_region |>
  mutate(phase = phase_label(enso_phase)) |>
  select(country, parent_region_name, phase, effect) |>
  pivot_wider(names_from = phase, values_from = effect) |>
  ggplot(aes(x = `El Niño versus neutral`, y = `La Niña versus neutral`,
             shape = country)) +
  geom_hline(yintercept = 0, colour = "grey72", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey72", linewidth = 0.4) +
  geom_abline(slope = -1, intercept = 0, linetype = 3,
              colour = "grey65", linewidth = 0.4) +
  geom_point(size = 2.6, colour = "grey20", alpha = 0.9) +
  scale_shape_manual(values = c("Argentina" = 16, "Brazil" = 17)) +
  scale_x_continuous(limits = c(-20, 25), breaks = seq(-20, 20, 10)) +
  scale_y_continuous(limits = c(-20, 25), breaks = seq(-20, 20, 10)) +
  labs(x = "El Niño response (percentage points)",
       y = "La Niña response (percentage points)") +
  theme_rogerio(10.5)

figure_4 <- (tag_panel(matrix_argentina, "a") |
             tag_panel(matrix_brazil, "b")) /
  tag_panel(regional_phase_contrast, "c") +
  plot_layout(heights = c(2.35, 0.78), guides = "collect") &
  theme(legend.position = "top")

save_figure(figure_4, "Figure_4_regional_evidence_matrix", 13, 8.7)

# 10. FIGURE 5: COMPOUND EVENTS -----------------------------------------------

compound_labels <- c(
  hot_dry_fraction = "Hot and dry",
  hot_wet_fraction = "Hot and wet",
  cold_wet_fraction = "Cold and wet",
  cold_dry_fraction = "Cold and dry"
)

compound_plot <- compounds |>
  filter(climate_metric %in% names(compound_labels)) |>
  mutate(
    event = recode(climate_metric, !!!compound_labels),
    phase = recode(enso_phase, `El Nino` = "El Niño", `La Nina` = "La Niña"),
    phase = factor(phase, levels = c("El Niño", "La Niña")),
    window = recode(
      climate_window,
      growing_season = "Growing season",
      reproductive = "Reproductive period",
      .default = str_to_sentence(str_replace_all(climate_window, "_", " "))
    ),
    window = factor(window, levels = c("Growing season", "Reproductive period")),
    region_phase = paste(country, phase, sep = ": ")
  )

compound_duplicates <- compound_plot |>
  count(country, window, event, phase, name = "n") |>
  filter(n > 1L)

if (nrow(compound_duplicates) > 0) {
  stop(
    "Compound event estimates are not unique by country, climate window, ",
    "event, and ENSO phase. Inspect the source table before plotting.",
    call. = FALSE
  )
}

compound_base <- ggplot(
  compound_plot,
  aes(
    x = phase_difference_percentage_points, y = event,
    colour = phase, shape = window,
    group = interaction(phase, window)
  )
) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.45) +
  geom_errorbarh(
    aes(xmin = ci_lower_percentage_points, xmax = ci_upper_percentage_points),
    height = 0, linewidth = 0.55,
    position = position_dodge(width = 0.48)
  ) +
  geom_point(size = 2.7, position = position_dodge(width = 0.48)) +
  scale_colour_manual(values = c("El Niño" = "#B2182B", "La Niña" = "#2166AC")) +
  scale_shape_manual(values = c("Growing season" = 16, "Reproductive period" = 17)) +
  labs(x = "Difference from neutral years (percentage points)", y = NULL) +
  theme_rogerio(11)

compound_argentina <- compound_base %+%
  filter(compound_plot, country == "Argentina") + labs(title = "Argentina")
compound_brazil <- compound_base %+%
  filter(compound_plot, country == "Brazil") + labs(title = "Brazil")

compound_xlim <- range(
  c(
    compound_plot$ci_lower_percentage_points,
    compound_plot$ci_upper_percentage_points
  ),
  na.rm = TRUE
)
compound_padding <- diff(compound_xlim) * 0.04
compound_xlim <- compound_xlim + c(-compound_padding, compound_padding)

compound_argentina <- compound_argentina + coord_cartesian(xlim = compound_xlim)
compound_brazil <- compound_brazil + coord_cartesian(xlim = compound_xlim)

compound_separation <- compound_plot |>
  select(country, window, event, phase, phase_difference_percentage_points) |>
  pivot_wider(
    id_cols = c(country, window, event),
    names_from = phase,
    values_from = phase_difference_percentage_points
  ) |>
  mutate(separation = `El Niño` - `La Niña`) |>
  ggplot(aes(x = separation, y = event, shape = window)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.45) +
  geom_point(
    size = 2.7, colour = "grey20",
    position = position_dodge(width = 0.45)
  ) +
  facet_wrap(~country, nrow = 1, axes = "all", axis.labels = "all") +
  scale_shape_manual(values = c("Growing season" = 16, "Reproductive period" = 17)) +
  labs(x = "El Niño minus La Niña (percentage points)", y = NULL) +
  theme_rogerio(11)

figure_5 <- (tag_panel(compound_argentina, "a") |
             tag_panel(compound_brazil, "b")) /
  tag_panel(compound_separation, "c") +
  plot_layout(heights = c(1.35, 0.72), guides = "collect") &
  theme(legend.position = "top")

save_figure(figure_5, "Figure_5_compound_climate_events_by_ENSO_phase", 12.5, 7.4)

# 11. TABLES AND QUALITY CONTROL ----------------------------------------------

country_story <- field_country |>
  select(
    country, spatial_field_cv_r_squared_percent,
    field_r2_ci_lower_percent, field_r2_ci_upper_percent,
    spatial_signal_energy_percent, national_net_cv_r_squared_percent,
    cancellation_index
  ) |>
  left_join(
    climate_country |>
      filter(model %in% c("ENSO", "climate", "water_balance", "climate_plus_water")) |>
      transmute(country, model,
                local_out_of_sample_r_squared_percent = 100 * out_of_sample_r_squared) |>
      pivot_wider(names_from = model, values_from = local_out_of_sample_r_squared_percent,
                  names_prefix = "local_"),
    by = "country"
  )

regional_story_table <- region_story |>
  select(country, parent_region_name, metric, value) |>
  pivot_wider(names_from = metric, values_from = value) |>
  arrange(country, desc(`Climate predictive skill`))

write_csv(country_story,
          file.path(output_dir, "Table_1_country_story_metrics.csv"), na = "")
write_csv(regional_story_table,
          file.path(output_dir, "Table_2_regional_story_metrics.csv"), na = "")
write_csv(pathway_data,
          file.path(output_dir, "Table_3_ENSO_climate_pathways.csv"), na = "")

qc <- tibble::tribble(
  ~test, ~value,
  "Mapped spatial units", n_distinct(phase_map$country, phase_map$spatial_unit_id_join),
  "Phase specific local estimates", nrow(phase_map),
  "Country summaries", nrow(field_country),
  "Regional summaries", n_distinct(region_story$country, region_story$parent_region_name),
  "Mechanism pathway rows", nrow(pathway_data),
  "Compound event rows", nrow(compound_plot)
)
write_csv(qc, file.path(output_dir, "QC_final_story_figures.csv"), na = "")

if (qc$value[qc$test == "Mapped spatial units"] != 2146) {
  stop("The spatial figures did not retain all 2146 units.", call. = FALSE)
}

message("Script 13 version 3.2 completed successfully.")
message("Story driven figures saved in: ", output_dir)

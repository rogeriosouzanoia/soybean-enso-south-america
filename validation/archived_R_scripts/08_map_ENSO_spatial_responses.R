# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 08_map_ENSO_spatial_responses.R
# PURPOSE: Estimate and map phase-specific ENSO responses in one South America
#          frame, using the author's publication theme
# VERSION: 3.1
# ==============================================================================


# 1. USER SETTINGS -------------------------------------------------------------

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")

project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) project_dir_override else default_project_dir

brazil_geometry_year <- 2022L
minimum_brazil_match_rate <- 0.99
minimum_argentina_match_rate <- 0.97
minimum_units_for_regional_pooling <- 5L
posterior_probability_threshold <- 0.95
temporal_sign_stability_threshold <- 0.80
number_temporal_folds <- 5L
colour_limit_percentage_points <- 15
figure_width_in <- 12
figure_height_in <- 8.2
figure_dpi <- 600


# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c(
  "dplyr", "readr", "tibble", "purrr", "sf", "ggplot2", "scales",
  "geobr", "geoAr", "rnaturalearth", "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(sf)
  library(ggplot2)
  library(scales)
  library(geobr)
  library(geoAr)
  library(rnaturalearth)
  library(patchwork)
})


# 3. PATHS ---------------------------------------------------------------------

merged_dir <- file.path(project_dir, "data_merged")
variance_dir <- file.path(project_dir, "enso_variance_models")
output_dir <- file.path(project_dir, "figures", "enso_spatial_responses")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panel_file <- file.path(merged_dir, "soybean_ENSO_spatial_panel_1974_2024.csv")
spatial_variance_file <- file.path(variance_dir, "ENSO_spatial_unit_models.csv")
country_variance_file <- file.path(variance_dir, "ENSO_spatial_panel_fixed_effects.csv")
if (!file.exists(panel_file)) stop("Spatial panel not found: ", panel_file, call. = FALSE)
if (!file.exists(spatial_variance_file)) stop("Spatial variance file not found: ", spatial_variance_file, call. = FALSE)
if (!file.exists(country_variance_file)) stop("Country variance file not found: ", country_variance_file, call. = FALSE)


# 4. AUTHOR THEME --------------------------------------------------------------

theme_rogerio <- function(base_size = 14) {
  theme_minimal(base_size = base_size, base_family = "Times New Roman") +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.45),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "top",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = grid::unit(0.8, "lines"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
}


# 5. HELPERS -------------------------------------------------------------------

check_required_columns <- function(data, required, dataset_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(dataset_name, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

standardize_numeric_code <- function(x, width) {
  x <- gsub("[^0-9]", "", as.character(x))
  x[x == ""] <- NA_character_
  ifelse(is.na(x), NA_character_, sprintf(paste0("%0", width, "d"), as.integer(x)))
}

first_existing_column <- function(data, candidates, dataset_name) {
  found <- candidates[candidates %in% names(data)]
  if (length(found) == 0) {
    stop("No identifier column found in ", dataset_name, ". Columns: ",
         paste(names(data), collapse = ", "), call. = FALSE)
  }
  found[1]
}

standardize_geometry_column <- function(x) {
  active <- attr(x, "sf_column")
  if (!identical(active, "geometry")) names(x)[names(x) == active] <- "geometry"
  st_geometry(x) <- "geometry"
  x
}

read_geobr_safely <- function(reader, dataset_name, ...) {
  result <- tryCatch(reader(..., cache = TRUE), error = function(e) {
    message("Cached ", dataset_name, " could not be read: ", conditionMessage(e))
    NULL
  })
  if (is.null(result) || nrow(result) == 0) {
    message("Downloading a fresh copy of ", dataset_name, "...")
    result <- tryCatch(reader(..., cache = FALSE), error = function(e) {
      stop("Fresh download of ", dataset_name, " failed: ", conditionMessage(e), call. = FALSE)
    })
  }
  if (is.null(result) || nrow(result) == 0) stop(dataset_name, " returned no features.")
  if (!inherits(result, "sf")) result <- st_as_sf(result)
  result
}

save_figure <- function(plot, stem, width = figure_width_in, height = figure_height_in) {
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot, width = width,
         height = height, units = "in", dpi = figure_dpi, bg = "white")
  ggsave(file.path(output_dir, paste0(stem, ".tiff")), plot, width = width,
         height = height, units = "in", dpi = figure_dpi,
         compression = "lzw", bg = "white")
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot, width = width,
         height = height, units = "in",
         device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf",
         bg = "white")
}

phase_panel_label <- function(x) {
  ifelse(x == "El Nino", "El Niño versus neutral", "La Niña versus neutral")
}

make_year_folds <- function(year) {
  ranked_years <- sort(unique(year))
  fold_lookup <- tibble(
    analysis_year = ranked_years,
    fold = pmin(
      number_temporal_folds,
      ceiling(seq_along(ranked_years) * number_temporal_folds / length(ranked_years))
    )
  )
  fold_lookup$fold[match(year, fold_lookup$analysis_year)]
}

welch_contrast <- function(y_phase, y_neutral) {
  n_phase <- sum(is.finite(y_phase))
  n_neutral <- sum(is.finite(y_neutral))
  if (n_phase < 3 || n_neutral < 3) {
    return(tibble(effect = NA_real_, standard_error = NA_real_,
                  ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_))
  }
  y_phase <- y_phase[is.finite(y_phase)]
  y_neutral <- y_neutral[is.finite(y_neutral)]
  effect <- mean(y_phase) - mean(y_neutral)
  se <- sqrt(var(y_phase) / n_phase + var(y_neutral) / n_neutral)
  if (!is.finite(se) || se <= 0) {
    return(tibble(effect = effect, standard_error = NA_real_,
                  ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_))
  }
  numerator <- (var(y_phase) / n_phase + var(y_neutral) / n_neutral)^2
  denominator <- (var(y_phase) / n_phase)^2 / (n_phase - 1) +
    (var(y_neutral) / n_neutral)^2 / (n_neutral - 1)
  df <- numerator / denominator
  critical <- qt(0.975, df)
  tibble(
    effect = effect,
    standard_error = se,
    ci_lower = effect - critical * se,
    ci_upper = effect + critical * se,
    p_value = 2 * pt(-abs(effect / se), df)
  )
}

calculate_phase_contrasts <- function(data) {
  folds <- make_year_folds(data$analysis_year)
  map_dfr(c("El Nino", "La Nina"), function(target_phase) {
    full <- welch_contrast(
      data$yield_anomaly_percent[data$enso_phase == target_phase],
      data$yield_anomaly_percent[data$enso_phase == "Neutral"]
    )
    fold_effects <- map_dbl(seq_len(number_temporal_folds), function(fold_id) {
      keep <- folds != fold_id
      result <- welch_contrast(
        data$yield_anomaly_percent[keep & data$enso_phase == target_phase],
        data$yield_anomaly_percent[keep & data$enso_phase == "Neutral"]
      )
      result$effect
    })
    target_sign <- sign(full$effect)
    sign_stability <- if (is.finite(target_sign) && target_sign != 0) {
      mean(sign(fold_effects[is.finite(fold_effects)]) == target_sign)
    } else NA_real_
    full |>
      mutate(
        enso_phase = target_phase,
        phase_years = sum(data$enso_phase == target_phase, na.rm = TRUE),
        neutral_years = sum(data$enso_phase == "Neutral", na.rm = TRUE),
        temporal_sign_stability = sign_stability
      )
  })
}

estimate_tau_squared <- function(effect, variance) {
  valid <- is.finite(effect) & is.finite(variance) & variance > 0
  effect <- effect[valid]
  variance <- variance[valid]
  k <- length(effect)
  if (k < minimum_units_for_regional_pooling) return(0)
  weight <- 1 / variance
  fixed_mean <- sum(weight * effect) / sum(weight)
  q <- sum(weight * (effect - fixed_mean)^2)
  c_value <- sum(weight) - sum(weight^2) / sum(weight)
  max(0, (q - (k - 1)) / c_value)
}

shrink_within_region <- function(data) {
  valid <- is.finite(data$effect) & is.finite(data$standard_error) &
    data$standard_error > 0
  result <- data |>
    mutate(
      between_unit_variance = NA_real_,
      effect_eb = effect,
      effect_eb_standard_error = standard_error,
      shrinkage_fraction = 0
    )
  if (sum(valid) < minimum_units_for_regional_pooling) return(result)
  variance <- data$standard_error[valid]^2
  tau_squared <- estimate_tau_squared(data$effect[valid], variance)
  weights <- 1 / (variance + tau_squared)
  regional_mean <- sum(weights * data$effect[valid]) / sum(weights)
  regional_mean_variance <- 1 / sum(weights)
  prior_variance <- tau_squared + regional_mean_variance
  posterior_variance <- 1 / (1 / variance + 1 / prior_variance)
  posterior_mean <- posterior_variance *
    (data$effect[valid] / variance + regional_mean / prior_variance)
  result$between_unit_variance[valid] <- tau_squared
  result$effect_eb[valid] <- posterior_mean
  result$effect_eb_standard_error[valid] <- sqrt(posterior_variance)
  result$shrinkage_fraction[valid] <- 1 - posterior_variance / variance
  result
}


# 6. PHASE-SPECIFIC LOCAL AND REGIONAL RESPONSES -------------------------------

message("Estimating separate El Niño and La Niña contrasts...")

panel <- read_csv(
  panel_file,
  col_types = cols(
    spatial_unit_id = col_character(), parent_region_code = col_character(),
    analysis_year = col_integer(), .default = col_guess()
  ),
  show_col_types = FALSE
)

spatial_variance <- read_csv(
  spatial_variance_file,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
)

country_variance <- read_csv(
  country_variance_file,
  show_col_types = FALSE
)

check_required_columns(
  spatial_variance,
  c("country", "spatial_unit_id", "spatial_unit_name", "parent_region_name",
    "r_squared", "variance_explained_percent", "cv_r_squared"),
  "Spatial variance dataset"
)

check_required_columns(
  country_variance,
  c("country", "area_weighted_variance_explained_percent"),
  "Country variance dataset"
)

check_required_columns(
  panel,
  c("country", "spatial_unit_id", "spatial_unit_name", "parent_region_code",
    "parent_region_name", "spatial_level", "analysis_year",
    "yield_anomaly_percent", "oni_harvest_year_phase",
    "production_t", "selected_trend_t_ha", "harvested_area_ha"),
  "Spatial ENSO panel"
)

panel <- panel |>
  mutate(
    enso_phase = case_when(
      grepl("El", as.character(oni_harvest_year_phase), ignore.case = TRUE) ~ "El Nino",
      grepl("La", as.character(oni_harvest_year_phase), ignore.case = TRUE) ~ "La Nina",
      TRUE ~ "Neutral"
    )
  )

local_effects <- panel |>
  filter(is.finite(yield_anomaly_percent)) |>
  group_by(country, spatial_unit_id, spatial_unit_name, parent_region_code,
           parent_region_name, spatial_level) |>
  group_modify(~ calculate_phase_contrasts(.x)) |>
  ungroup() |>
  group_by(country, parent_region_code, parent_region_name, enso_phase) |>
  group_modify(~ shrink_within_region(.x)) |>
  ungroup() |>
  mutate(
    effect_eb_ci_lower = effect_eb - 1.96 * effect_eb_standard_error,
    effect_eb_ci_upper = effect_eb + 1.96 * effect_eb_standard_error,
    posterior_probability_positive = pnorm(effect_eb / effect_eb_standard_error),
    posterior_probability_negative = pnorm(-effect_eb / effect_eb_standard_error),
    # Phase evidence is based on each unit's own estimate. The EB estimate is
    # retained for audit, but is not used to declare phase-specific robustness.
    # This avoids pseudoreplication when many units share the same ENSO years.
    robust_phase_class = case_when(
      effect > 0 & p_value < 0.05 &
        temporal_sign_stability >= temporal_sign_stability_threshold ~ "Higher than neutral",
      effect < 0 & p_value < 0.05 &
        temporal_sign_stability >= temporal_sign_stability_threshold ~ "Lower than neutral",
      TRUE ~ "Limited or uncertain"
    ),
    phase_panel = factor(
      phase_panel_label(enso_phase),
      levels = c("El Niño versus neutral", "La Niña versus neutral")
    )
  )

region_year <- panel |>
  filter(is.finite(production_t), is.finite(selected_trend_t_ha),
         is.finite(harvested_area_ha), harvested_area_ha > 0) |>
  group_by(country, parent_region_code, parent_region_name, analysis_year, enso_phase) |>
  summarise(
    observed_production_t = sum(production_t, na.rm = TRUE),
    expected_production_t = sum(selected_trend_t_ha * harvested_area_ha, na.rm = TRUE),
    yield_anomaly_percent = 100 *
      (observed_production_t / expected_production_t - 1),
    .groups = "drop"
  )

regional_effects <- region_year |>
  group_by(country, parent_region_code, parent_region_name) |>
  group_modify(~ calculate_phase_contrasts(.x)) |>
  ungroup() |>
  mutate(
    phase_panel = factor(
      phase_panel_label(enso_phase),
      levels = c("El Niño versus neutral", "La Niña versus neutral")
    ),
    evidence_class = case_when(
      p_value < 0.05 & temporal_sign_stability >= temporal_sign_stability_threshold ~
        "Supported and temporally stable",
      temporal_sign_stability >= temporal_sign_stability_threshold ~
        "Temporally stable",
      TRUE ~ "Limited evidence"
    )
  )

write_csv(local_effects, file.path(output_dir, "ENSO_local_phase_specific_effects.csv"))
write_csv(regional_effects, file.path(output_dir, "ENSO_state_province_phase_effects.csv"))


# 7. OFFICIAL GEOMETRIES -------------------------------------------------------

message("Reading official Brazilian municipality geometries...")
brazil_municipalities <- read_geobr_safely(
  geobr::read_municipality, "Brazilian municipality geometry",
  year = brazil_geometry_year, code_muni = "all", simplified = TRUE,
  showProgress = TRUE
) |>
  st_make_valid() |>
  standardize_geometry_column() |>
  st_transform(4326)

brazil_code_column <- first_existing_column(
  brazil_municipalities, c("code_muni", "code_municipality", "code_muni_7"),
  "geobr municipality geometry"
)
brazil_municipalities <- brazil_municipalities |>
  mutate(
    spatial_unit_id_join = standardize_numeric_code(.data[[brazil_code_column]], 7),
    region_code_join = substr(spatial_unit_id_join, 1, 2)
  )
# Use the official state layer directly. Dissolving thousands of municipality
# polygons can preserve slivers and internal line fragments in the state layer.
brazil_regions <- read_geobr_safely(
  geobr::read_state, "Brazilian state geometry",
  year = 2020L, code_state = "all", simplified = TRUE,
  showProgress = TRUE
) |>
  st_make_valid() |>
  standardize_geometry_column() |>
  st_transform(4326)
brazil_country <- brazil_municipalities |>
  summarise(country = "Brazil") |>
  st_make_valid()

message("Reading official Argentine department geometries...")
argentina_departments <- geoAr::get_geo(
  geo = "ARGENTINA", level = "departamento", simplified = TRUE
) |>
  st_as_sf() |>
  st_make_valid() |>
  standardize_geometry_column() |>
  st_transform(4326)

arg_department_column <- first_existing_column(
  argentina_departments,
  c("coddepto_censo", "codigo_departamento", "id_departamento", "coddepto", "INDEC", "indec"),
  "geoAr department geometry"
)
arg_province_code_column <- first_existing_column(
  argentina_departments,
  c("codprov_censo", "codigo_provincia", "id_provincia", "codprov", "provincia_id"),
  "geoAr department geometry"
)
argentina_departments <- argentina_departments |>
  mutate(
    department_code_raw = gsub("[^0-9]", "", as.character(.data[[arg_department_column]])),
    region_code_join = standardize_numeric_code(.data[[arg_province_code_column]], 2),
    department_code_3 = standardize_numeric_code(department_code_raw, 3),
    spatial_unit_id_join = if_else(
      nchar(department_code_raw) <= 3,
      paste0(region_code_join, department_code_3),
      standardize_numeric_code(department_code_raw, 5)
    )
  )
argentina_regions <- geoAr::get_geo(
  geo = "ARGENTINA", level = "provincia", simplified = TRUE
) |>
  st_as_sf() |>
  st_make_valid() |>
  standardize_geometry_column() |>
  st_transform(4326)
argentina_country <- argentina_departments |>
  summarise(country = "Argentina") |>
  st_make_valid()

# All South American countries are shown for geographical context. Analytical
# fills remain restricted to Brazil and Argentina.
south_america_countries <- rnaturalearth::ne_countries(
  scale = "medium", continent = "South America", returnclass = "sf"
) |>
  st_make_valid() |>
  st_transform(4326)


# 8. JOINS AND AUDIT -----------------------------------------------------------

brazil_effects <- local_effects |>
  filter(country == "Brazil") |>
  mutate(spatial_unit_id_join = standardize_numeric_code(spatial_unit_id, 7))
argentina_effects <- local_effects |>
  filter(country == "Argentina") |>
  mutate(spatial_unit_id_join = standardize_numeric_code(spatial_unit_id, 5))

brazil_match <- brazil_effects |>
  distinct(country, spatial_unit_id, spatial_unit_name, spatial_unit_id_join) |>
  mutate(geometry_matched = spatial_unit_id_join %in% brazil_municipalities$spatial_unit_id_join)
argentina_match <- argentina_effects |>
  distinct(country, spatial_unit_id, spatial_unit_name, spatial_unit_id_join) |>
  mutate(geometry_matched = spatial_unit_id_join %in% argentina_departments$spatial_unit_id_join)
join_audit <- bind_rows(brazil_match, argentina_match)

join_summary <- join_audit |>
  group_by(country) |>
  summarise(
    analytical_units = n(), matched_units = sum(geometry_matched),
    unmatched_units = sum(!geometry_matched), match_rate = mean(geometry_matched),
    match_rate_percent = 100 * match_rate, .groups = "drop"
  )

write_csv(join_audit, file.path(output_dir, "spatial_join_audit_all_units.csv"))
write_csv(filter(join_audit, !geometry_matched),
          file.path(output_dir, "spatial_join_unmatched_units.csv"))
write_csv(join_summary, file.path(output_dir, "QC_spatial_join_summary.csv"))

brazil_rate <- join_summary$match_rate[join_summary$country == "Brazil"]
argentina_rate <- join_summary$match_rate[join_summary$country == "Argentina"]
if (length(brazil_rate) != 1 || brazil_rate < minimum_brazil_match_rate) {
  stop("Brazil geometry match rate below the required threshold.", call. = FALSE)
}
if (length(argentina_rate) != 1 || argentina_rate < minimum_argentina_match_rate) {
  stop("Argentina geometry match rate below the required threshold. Inspect join audit files.",
       call. = FALSE)
}

brazil_map <- brazil_municipalities |>
  inner_join(brazil_effects, by = "spatial_unit_id_join") |>
  select(country, spatial_unit_id, spatial_unit_name, parent_region_name,
         enso_phase, phase_panel, effect, standard_error, p_value,
         temporal_sign_stability,
         robust_phase_class, geometry)
argentina_map <- argentina_departments |>
  inner_join(argentina_effects, by = "spatial_unit_id_join") |>
  select(country, spatial_unit_id, spatial_unit_name, parent_region_name,
         enso_phase, phase_panel, effect, standard_error, p_value,
         temporal_sign_stability,
         robust_phase_class, geometry)

response_map <- bind_rows(brazil_map, argentina_map) |>
  st_as_sf() |>
  st_make_valid()
regional_outlines <- bind_rows(
  brazil_regions |> select(geometry),
  argentina_regions |> select(geometry)
) |>
  st_as_sf() |>
  st_make_valid()

phase_levels <- levels(local_effects$phase_panel)
south_america_phase <- bind_rows(lapply(phase_levels, function(label) {
  south_america_countries |> mutate(phase_panel = factor(label, levels = phase_levels))
}))
regional_outlines_phase <- bind_rows(lapply(phase_levels, function(label) {
  regional_outlines |> mutate(phase_panel = factor(label, levels = phase_levels))
}))

brazil_variance <- spatial_variance |>
  filter(country == "Brazil") |>
  mutate(spatial_unit_id_join = standardize_numeric_code(spatial_unit_id, 7))
argentina_variance <- spatial_variance |>
  filter(country == "Argentina") |>
  mutate(spatial_unit_id_join = standardize_numeric_code(spatial_unit_id, 5))

variance_map <- bind_rows(
  brazil_municipalities |>
    inner_join(brazil_variance, by = "spatial_unit_id_join") |>
    select(country, spatial_unit_id, spatial_unit_name, parent_region_name,
           variance_explained_percent, cv_r_squared, geometry),
  argentina_departments |>
    inner_join(argentina_variance, by = "spatial_unit_id_join") |>
    select(country, spatial_unit_id, spatial_unit_name, parent_region_name,
           variance_explained_percent, cv_r_squared, geometry)
) |>
  st_as_sf() |>
  st_make_valid() |>
  mutate(
    variance_class = cut(
      variance_explained_percent,
      breaks = c(-Inf, 2, 5, 10, 20, Inf),
      labels = c("0 to 2%", "2 to 5%", "5 to 10%", "10 to 20%", "Above 20%"),
      right = FALSE
    )
  )

st_write(
  response_map,
  file.path(output_dir, "ENSO_phase_specific_spatial_layers.gpkg"),
  layer = "local_phase_effects", delete_dsn = TRUE, quiet = TRUE
)


# 9. FIGURE 1: VARIABILITY EXPLAINED BY ENSO ----------------------------------

message("Generating the ENSO variance explained map...")

country_bar_data <- country_variance |>
  filter(country %in% c("Brazil", "Argentina")) |>
  transmute(
    country,
    variance_explained_percent = area_weighted_variance_explained_percent
  )

country_bar <- ggplot(
  country_bar_data,
  aes(y = reorder(country, variance_explained_percent),
      x = variance_explained_percent, fill = country)
) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", variance_explained_percent)),
            hjust = -0.25, family = "Times New Roman", size = 3.2) +
  scale_fill_manual(values = c("Argentina" = "#2B8CBE", "Brazil" = "#7B3294")) +
  scale_x_continuous(limits = c(0, 9), breaks = c(0, 4, 8),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "National estimates", x = "Explained variability (%)", y = NULL) +
  theme_rogerio(base_size = 9) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 9.5, hjust = 0),
    plot.background = element_rect(fill = "white", colour = "grey70"),
    panel.border = element_rect(fill = NA, colour = "grey70", linewidth = 0.3),
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(4, 7, 4, 5)
  )

variance_base <- ggplot() +
  geom_sf(data = south_america_countries, fill = "grey97", colour = "black", linewidth = 0.38) +
  geom_sf(data = variance_map, aes(fill = variance_class), colour = NA) +
  geom_sf(data = regional_outlines, fill = NA, colour = "grey62", linewidth = 0.22) +
  geom_sf(data = south_america_countries, fill = NA, colour = "black", linewidth = 0.42) +
  scale_fill_manual(
    values = c(
      "0 to 2%" = "#F2F2F2",
      "2 to 5%" = "#D6E6F4",
      "5 to 10%" = "#9ECAE1",
      "10 to 20%" = "#4292C6",
      "Above 20%" = "#08519C"
    ),
    drop = FALSE,
    na.value = "transparent"
  ) +
  coord_sf(xlim = c(-90, -32), ylim = c(-56, 13), expand = FALSE) +
  labs(
    title = "Share of local soybean productivity variability explained by ENSO",
    subtitle = "Colour shows local R² for DJF ONI. National bars show area weighted partial R².",
    x = "Longitude", y = "Latitude"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_rogerio() +
  theme(legend.key.width = grid::unit(1.15, "cm"))

figure_1 <- variance_base +
  inset_element(country_bar, left = 0.025, bottom = 0.08, right = 0.34, top = 0.25,
                align_to = "panel")

save_figure(figure_1, "Figure_1_South_America_ENSO_variance_explained")


# 10. FIGURE 2: PHASE-SPECIFIC PRODUCTIVITY RESPONSE --------------------------

message("Generating the clean phase-specific response map...")

figure_2 <- ggplot() +
  geom_sf(data = south_america_phase, fill = "grey97", colour = "black", linewidth = 0.38) +
  geom_sf(data = response_map, aes(fill = effect), colour = NA) +
  geom_sf(data = regional_outlines_phase, fill = NA, colour = "grey62", linewidth = 0.22) +
  geom_sf(data = south_america_phase, fill = NA, colour = "black", linewidth = 0.42) +
  facet_wrap(~ phase_panel, nrow = 1) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0,
    limits = c(-colour_limit_percentage_points, colour_limit_percentage_points),
    oob = scales::squish,
    breaks = c(-15, -7.5, 0, 7.5, 15),
    labels = c("−15", "−7.5", "0", "7.5", "15"),
    na.value = "transparent"
  ) +
  coord_sf(xlim = c(-82, -32), ylim = c(-56, 13), expand = FALSE) +
  labs(
    title = "Soybean productivity response to ENSO phases",
    subtitle = paste(
      "Difference from neutral years in percentage points.",
      "Blue indicates higher productivity and red indicates lower productivity."
    ),
    x = "Longitude", y = "Latitude"
  ) +
  guides(fill = guide_colourbar(barwidth = grid::unit(7.5, "cm"))) +
  theme_rogerio()

save_figure(figure_2, "Figure_2_South_America_ENSO_phase_effects")


# 11. REGIONAL TABLE AND DIAGNOSTIC FIGURE ------------------------------------

regional_variance_table <- spatial_variance |>
  group_by(country, parent_region_name) |>
  summarise(
    spatial_units = n(),
    median_local_variance_explained_percent = median(variance_explained_percent, na.rm = TRUE),
    area_unweighted_mean_local_variance_explained_percent = mean(variance_explained_percent, na.rm = TRUE),
    median_cv_r_squared = median(cv_r_squared, na.rm = TRUE),
    .groups = "drop"
  )

regional_results_table <- regional_effects |>
  select(country, parent_region_name, enso_phase, phase_years, neutral_years,
         effect, standard_error, ci_lower, ci_upper, p_value,
         temporal_sign_stability, evidence_class) |>
  left_join(regional_variance_table, by = c("country", "parent_region_name")) |>
  arrange(country, parent_region_name, enso_phase)

write_csv(
  regional_results_table,
  file.path(output_dir, "Table_state_province_ENSO_variance_and_phase_effects.csv")
)

regional_plot_data <- regional_effects |>
  filter(is.finite(effect), is.finite(ci_lower), is.finite(ci_upper)) |>
  mutate(
    region_label = parent_region_name,
    region_phase_id = paste(region_label, country, sep = "___"),
    region_phase_id = reorder(region_phase_id, effect)
  )

figure_3 <- ggplot(regional_plot_data, aes(x = effect, y = region_phase_id)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0,
                 linewidth = 0.5, colour = "grey40") +
  geom_point(aes(colour = evidence_class), size = 2.3) +
  facet_wrap(~ phase_panel, nrow = 1, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_colour_manual(values = c(
    "Supported and temporally stable" = "#542788",
    "Temporally stable" = "#2B8CBE",
    "Limited evidence" = "#969696"
  )) +
  labs(
    title = "State and province responses to ENSO phases",
    subtitle = "Points show differences from neutral years. Lines show 95% confidence intervals.",
    x = "Productivity anomaly difference from neutral years, percentage points",
    y = NULL
  ) +
  theme_rogerio(base_size = 11)

save_figure(figure_3, "Figure_3_state_province_ENSO_phase_effects", 12, 10)


# 12. QUALITY CONTROL ----------------------------------------------------------

qc_phase_maps <- local_effects |>
  group_by(country, enso_phase, phase_panel) |>
  summarise(
    spatial_units = n(),
    finite_effects = sum(is.finite(effect_eb)),
    higher_than_neutral = sum(robust_phase_class == "Higher than neutral", na.rm = TRUE),
    lower_than_neutral = sum(robust_phase_class == "Lower than neutral", na.rm = TRUE),
    limited_or_uncertain = sum(robust_phase_class == "Limited or uncertain", na.rm = TRUE),
    median_effect_percentage_points = median(effect, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(qc_phase_maps, file.path(output_dir, "QC_ENSO_phase_maps.csv"))

rs_phase_audit <- regional_effects |>
  filter(country == "Brazil", parent_region_name == "RS") |>
  mutate(
    supported_at_005 = p_value < 0.05,
    interpretation = case_when(
      supported_at_005 & effect > 0 ~ "Higher than neutral, supported",
      supported_at_005 & effect < 0 ~ "Lower than neutral, supported",
      TRUE ~ "Direction uncertain"
    )
  )
write_csv(
  rs_phase_audit,
  file.path(output_dir, "QC_Rio_Grande_do_Sul_ENSO_phase_effects.csv")
)

message("Script 08 completed successfully.")
message("All South American countries are shown. Analytical fills cover Brazil and Argentina.")
message("Effects are reported separately as El Niño versus neutral and La Niña versus neutral.")
message("Outputs saved in: ", output_dir)

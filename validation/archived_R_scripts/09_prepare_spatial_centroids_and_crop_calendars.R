# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 09_prepare_spatial_centroids_and_crop_calendars.R
# PURPOSE: Prepare spatial sampling points and regional soybean calendars for
#          the climate mechanism analysis
# VERSION: 1.2
# ==============================================================================


# 1. USER SETTINGS -------------------------------------------------------------

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")

project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) project_dir_override else default_project_dir

# NASA POWER daily meteorology is consistently available from 1981 onward.
climate_start_year <- 1981L
climate_end_year <- 2024L

# The broad download window surrounds all regional crop calendars. Precise
# phenological windows are applied only after the weather data are downloaded.
download_start_month_previous_year <- 8L
download_end_month_harvest_year <- 7L

minimum_geometry_match_rate <- 0.99
recent_reference_start_year <- 2015L
recent_reference_end_year <- 2024L


# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c("dplyr", "readr", "tibble", "tidyr", "stringr", "sf")
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
  library(tidyr)
  library(stringr)
  library(sf)
})


# 3. PATHS ---------------------------------------------------------------------

merged_dir <- file.path(project_dir, "data_merged")
map_dir <- file.path(project_dir, "figures", "enso_spatial_responses")
output_dir <- file.path(project_dir, "climate_mechanisms")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panel_file <- file.path(merged_dir, "soybean_ENSO_spatial_panel_1974_2024.csv")
geometry_file <- file.path(map_dir, "ENSO_phase_specific_spatial_layers.gpkg")

if (!file.exists(panel_file)) {
  stop("Spatial soybean panel not found: ", panel_file, call. = FALSE)
}
if (!file.exists(geometry_file)) {
  stop(
    "Validated spatial layer not found: ", geometry_file,
    ". Run the final version of script 08 first.", call. = FALSE
  )
}


# 4. HELPERS -------------------------------------------------------------------

check_required_columns <- function(data, required, dataset_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(dataset_name, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

standardize_numeric_code <- function(x, width) {
  value <- gsub("[^0-9]", "", as.character(x))
  value[value == ""] <- NA_character_
  ifelse(
    is.na(value), NA_character_,
    sprintf(paste0("%0", width, "d"), as.integer(value))
  )
}

standardize_unit_id <- function(country, spatial_unit_id) {
  ifelse(
    country == "Brazil",
    standardize_numeric_code(spatial_unit_id, 7L),
    standardize_numeric_code(spatial_unit_id, 5L)
  )
}

standardize_region_code <- function(x) standardize_numeric_code(x, 2L)

derive_parent_region_code <- function(country, parent_region_code, parent_region_name) {
  standardized_code <- standardize_region_code(parent_region_code)

  brazil_state_lookup <- c(
    AC = "12", AL = "27", AP = "16", AM = "13", BA = "29", CE = "23",
    DF = "53", ES = "32", GO = "52", MA = "21", MT = "51", MS = "50",
    MG = "31", PA = "15", PB = "25", PR = "41", PE = "26", PI = "22",
    RJ = "33", RN = "24", RS = "43", RO = "11", RR = "14", SC = "42",
    SP = "35", SE = "28", TO = "17"
  )

  brazil_state_names <- c(
    "ACRE" = "12", "ALAGOAS" = "27", "AMAPA" = "16",
    "AMAZONAS" = "13", "BAHIA" = "29", "CEARA" = "23",
    "DISTRITO FEDERAL" = "53", "ESPIRITO SANTO" = "32",
    "GOIAS" = "52", "MARANHAO" = "21", "MATO GROSSO" = "51",
    "MATO GROSSO DO SUL" = "50", "MINAS GERAIS" = "31",
    "PARA" = "15", "PARAIBA" = "25", "PARANA" = "41",
    "PERNAMBUCO" = "26", "PIAUI" = "22", "RIO DE JANEIRO" = "33",
    "RIO GRANDE DO NORTE" = "24", "RIO GRANDE DO SUL" = "43",
    "RONDONIA" = "11", "RORAIMA" = "14", "SANTA CATARINA" = "42",
    "SAO PAULO" = "35", "SERGIPE" = "28", "TOCANTINS" = "17"
  )

  region_name_clean <- toupper(trimws(as.character(parent_region_name)))
  region_name_ascii <- iconv(
    region_name_clean,
    from = "UTF-8", to = "ASCII//TRANSLIT"
  )
  region_name_ascii[is.na(region_name_ascii)] <- region_name_clean[is.na(region_name_ascii)]

  brazil_code_from_name <- unname(brazil_state_lookup[region_name_ascii])
  missing_from_abbreviation <- is.na(brazil_code_from_name)
  brazil_code_from_name[missing_from_abbreviation] <- unname(
    brazil_state_names[region_name_ascii[missing_from_abbreviation]]
  )

  use_name <- country == "Brazil" & is.na(standardized_code)
  standardized_code[use_name] <- brazil_code_from_name[use_name]
  standardized_code
}

month_sequence <- function(start_month, end_month) {
  if (start_month <= end_month) {
    seq.int(start_month, end_month)
  } else {
    c(seq.int(start_month, 12L), seq.int(1L, end_month))
  }
}

month_label <- function(x) paste(month.abb[x], collapse = "-")

first_layer_name <- function(gpkg_file) {
  layers <- st_layers(gpkg_file)$name
  preferred <- c("local_phase_effects", "local_effects", "phase_effects")
  selected <- preferred[preferred %in% layers]
  if (length(selected) > 0) return(selected[1])
  if (length(layers) == 0) stop("The GeoPackage contains no spatial layers.", call. = FALSE)
  layers[1]
}

standardize_geometry_column <- function(x) {
  if (!inherits(x, "sf")) {
    stop("The spatial object is not an sf object.", call. = FALSE)
  }
  active_geometry <- attr(x, "sf_column")
  if (is.null(active_geometry) || !active_geometry %in% names(x)) {
    stop("No active geometry column was found in the spatial object.", call. = FALSE)
  }
  if (!identical(active_geometry, "geometry")) {
    names(x)[names(x) == active_geometry] <- "geometry"
    st_geometry(x) <- "geometry"
  }
  x
}


# 5. READ AND STANDARDIZE INPUTS ----------------------------------------------

message("Reading the validated soybean panel and spatial layer...")

panel <- read_csv(panel_file, show_col_types = FALSE)
check_required_columns(
  panel,
  c(
    "country", "spatial_unit_id", "spatial_unit_name", "parent_region_code",
    "parent_region_name", "spatial_level", "analysis_year", "harvested_area_ha",
    "production_t"
  ),
  "soybean ENSO spatial panel"
)

panel_units <- panel |>
  mutate(
    country = as.character(country),
    spatial_unit_id_join = standardize_unit_id(country, spatial_unit_id),
    parent_region_code = derive_parent_region_code(
      country, parent_region_code, parent_region_name
    )
  ) |>
  filter(country %in% c("Brazil", "Argentina")) |>
  group_by(
    country, spatial_unit_id_join, spatial_unit_name, parent_region_code,
    parent_region_name, spatial_level
  ) |>
  summarise(
    first_soybean_year = min(analysis_year, na.rm = TRUE),
    last_soybean_year = max(analysis_year, na.rm = TRUE),
    soybean_years = n_distinct(analysis_year),
    mean_harvested_area_ha = mean(harvested_area_ha, na.rm = TRUE),
    recent_mean_harvested_area_ha = {
      values <- harvested_area_ha[
        analysis_year >= recent_reference_start_year &
          analysis_year <= recent_reference_end_year
      ]
      if (length(values) == 0 || all(!is.finite(values))) NA_real_ else mean(values, na.rm = TRUE)
    },
    recent_mean_production_t = {
      values <- production_t[
        analysis_year >= recent_reference_start_year &
          analysis_year <= recent_reference_end_year
      ]
      if (length(values) == 0 || all(!is.finite(values))) NA_real_ else mean(values, na.rm = TRUE)
    },
    .groups = "drop"
  ) |>
  mutate(
    reference_harvested_area_ha = coalesce(
      recent_mean_harvested_area_ha,
      mean_harvested_area_ha
    )
  )

geometry_layer <- first_layer_name(geometry_file)
geometry_raw <- st_read(geometry_file, layer = geometry_layer, quiet = TRUE) |>
  standardize_geometry_column() |>
  st_make_valid() |>
  standardize_geometry_column()

check_required_columns(
  geometry_raw,
  c("country", "spatial_unit_id", "spatial_unit_name"),
  "validated spatial layer"
)

# The script 08 layer contains one repeated geometry for each ENSO phase.
# Keeping the first record removes this intentional duplication without a
# computationally expensive or geometrically risky polygon union.
geometry_units <- geometry_raw |>
  mutate(
    country = as.character(country),
    spatial_unit_id_join = standardize_unit_id(country, spatial_unit_id)
  ) |>
  arrange(country, spatial_unit_id_join) |>
  group_by(country, spatial_unit_id_join) |>
  slice(1L) |>
  ungroup() |>
  select(country, spatial_unit_id_join, geometry) |>
  standardize_geometry_column()


# 6. GEOMETRY JOIN AND INTERNAL SAMPLING POINTS --------------------------------

join_audit <- panel_units |>
  st_drop_geometry() |>
  mutate(
    geometry_matched = paste(country, spatial_unit_id_join) %in%
      paste(geometry_units$country, geometry_units$spatial_unit_id_join)
  )

join_summary <- join_audit |>
  group_by(country) |>
  summarise(
    analytical_units = n(),
    matched_units = sum(geometry_matched),
    unmatched_units = sum(!geometry_matched),
    match_rate = mean(geometry_matched),
    .groups = "drop"
  )

write_csv(join_summary, file.path(output_dir, "QC_centroid_geometry_join_summary.csv"))
write_csv(
  filter(join_audit, !geometry_matched),
  file.path(output_dir, "QC_centroid_unmatched_units.csv")
)

if (any(join_summary$match_rate < minimum_geometry_match_rate)) {
  stop(
    "Geometry match rate was below ", 100 * minimum_geometry_match_rate,
    "%. Inspect QC_centroid_unmatched_units.csv.", call. = FALSE
  )
}

units_sf <- geometry_units |>
  inner_join(
    panel_units,
    by = c("country", "spatial_unit_id_join")
  ) |>
  st_transform(6933)

# st_point_on_surface guarantees an internal representative point. This avoids
# offshore or neighbouring-unit coordinates for concave and multipart polygons.
sampling_points <- suppressWarnings(st_point_on_surface(units_sf)) |>
  st_transform(4326)

coordinates <- st_coordinates(sampling_points)

centroids <- sampling_points |>
  st_drop_geometry() |>
  mutate(
    longitude = coordinates[, "X"],
    latitude = coordinates[, "Y"],
    climate_location_id = paste0(
      if_else(country == "Brazil", "BRA", "ARG"), "_", spatial_unit_id_join
    )
  ) |>
  select(
    climate_location_id, country, spatial_unit_id = spatial_unit_id_join,
    spatial_unit_name, parent_region_code, parent_region_name, spatial_level,
    longitude, latitude, first_soybean_year, last_soybean_year, soybean_years,
    reference_harvested_area_ha, recent_mean_production_t
  ) |>
  arrange(country, parent_region_name, spatial_unit_name)

if (any(!is.finite(centroids$longitude)) || any(!is.finite(centroids$latitude))) {
  stop("Non-finite sampling coordinates were produced.", call. = FALSE)
}
if (any(centroids$longitude < -90 | centroids$longitude > -30) ||
    any(centroids$latitude < -60 | centroids$latitude > 15)) {
  stop("A sampling coordinate falls outside the South America validation box.", call. = FALSE)
}
if (anyDuplicated(centroids$climate_location_id)) {
  stop("Duplicated climate location identifiers were produced.", call. = FALSE)
}


# 7. REGIONAL SOYBEAN CALENDARS ------------------------------------------------

# These calendars define transparent regional analysis windows. They are broad
# agronomic priors rather than estimated windows, which prevents selecting a
# climate window because it happens to maximize association with productivity.
# A future curated file named soybean_regional_crop_calendars_user.csv in the
# output directory can replace these defaults without changing the script.

brazil_calendar <- tribble(
  ~parent_region_code, ~calendar_group,
  "11", "Brazil_Central", "12", "Brazil_North", "13", "Brazil_North",
  "14", "Brazil_North", "15", "Brazil_North", "16", "Brazil_North",
  "17", "Brazil_Central", "21", "Brazil_Northeast", "22", "Brazil_Northeast",
  "23", "Brazil_Northeast", "24", "Brazil_Northeast", "25", "Brazil_Northeast",
  "26", "Brazil_Northeast", "27", "Brazil_Northeast", "28", "Brazil_Northeast",
  "29", "Brazil_Northeast", "31", "Brazil_Southeast", "32", "Brazil_Southeast",
  "33", "Brazil_Southeast", "35", "Brazil_Southeast", "41", "Brazil_South",
  "42", "Brazil_South", "43", "Brazil_South", "50", "Brazil_Central",
  "51", "Brazil_Central", "52", "Brazil_Central", "53", "Brazil_Central"
) |>
  mutate(country = "Brazil", .before = 1)

argentina_calendar <- tribble(
  ~parent_region_code, ~calendar_group,
  "06", "Argentina_Central", "10", "Argentina_North", "14", "Argentina_Central",
  "18", "Argentina_Northeast", "22", "Argentina_Northeast",
  "30", "Argentina_Central", "34", "Argentina_Northeast", "38", "Argentina_North",
  "42", "Argentina_Central", "54", "Argentina_Northeast", "66", "Argentina_North",
  "74", "Argentina_Central", "82", "Argentina_Central", "86", "Argentina_North",
  "90", "Argentina_North"
) |>
  mutate(country = "Argentina", .before = 1)

calendar_definitions <- tribble(
  ~calendar_group, ~sowing_start_month, ~sowing_end_month, ~vegetative_start_month, ~vegetative_end_month, ~reproductive_start_month, ~reproductive_end_month, ~harvest_start_month, ~harvest_end_month,
  "Brazil_North",          11L,  1L, 12L,  2L,  2L,  4L,  3L,  6L,
  "Brazil_Northeast",      11L,  1L, 12L,  2L,  2L,  4L,  3L,  6L,
  "Brazil_Central",         9L, 11L, 10L,  1L, 12L,  3L,  1L,  4L,
  "Brazil_Southeast",      10L, 12L, 11L,  1L,  1L,  3L,  2L,  5L,
  "Brazil_South",          10L, 12L, 11L,  1L,  1L,  3L,  2L,  5L,
  "Argentina_North",       11L,  1L, 12L,  2L,  2L,  4L,  4L,  6L,
  "Argentina_Northeast",   10L,  1L, 11L,  2L,  1L,  4L,  3L,  6L,
  "Argentina_Central",     10L, 12L, 11L,  1L,  1L,  3L,  3L,  5L
) |>
  mutate(
    growing_season_start_month = sowing_start_month,
    growing_season_end_month = harvest_end_month,
    sowing_window = paste(month.abb[sowing_start_month], month.abb[sowing_end_month], sep = "-"),
    vegetative_window = paste(month.abb[vegetative_start_month], month.abb[vegetative_end_month], sep = "-"),
    reproductive_window = paste(month.abb[reproductive_start_month], month.abb[reproductive_end_month], sep = "-"),
    harvest_window = paste(month.abb[harvest_start_month], month.abb[harvest_end_month], sep = "-")
  )

regional_calendars_default <- bind_rows(brazil_calendar, argentina_calendar) |>
  left_join(calendar_definitions, by = "calendar_group") |>
  mutate(
    calendar_source = "Prespecified regional agronomic window for sensitivity analysis",
    calendar_status = "provisional"
  )

user_calendar_file <- file.path(output_dir, "soybean_regional_crop_calendars_user.csv")
if (file.exists(user_calendar_file)) {
  message("Using user-curated regional crop calendars: ", user_calendar_file)
  regional_calendars <- read_csv(user_calendar_file, show_col_types = FALSE)
  check_required_columns(
    regional_calendars,
    c(
      "country", "parent_region_code", "calendar_group",
      "sowing_start_month", "sowing_end_month",
      "vegetative_start_month", "vegetative_end_month",
      "reproductive_start_month", "reproductive_end_month",
      "harvest_start_month", "harvest_end_month"
    ),
    "user crop calendar"
  )
  if (!"calendar_source" %in% names(regional_calendars)) {
    regional_calendars$calendar_source <- "User-curated regional crop calendar"
  }
  regional_calendars <- regional_calendars |>
    mutate(
      parent_region_code = standardize_region_code(parent_region_code),
      growing_season_start_month = sowing_start_month,
      growing_season_end_month = harvest_end_month,
      sowing_window = paste(month.abb[sowing_start_month], month.abb[sowing_end_month], sep = "-"),
      vegetative_window = paste(month.abb[vegetative_start_month], month.abb[vegetative_end_month], sep = "-"),
      reproductive_window = paste(month.abb[reproductive_start_month], month.abb[reproductive_end_month], sep = "-"),
      harvest_window = paste(month.abb[harvest_start_month], month.abb[harvest_end_month], sep = "-"),
      calendar_source = coalesce(as.character(calendar_source), "User-curated regional crop calendar"),
      calendar_status = "user_curated"
    )
} else {
  regional_calendars <- regional_calendars_default
  write_csv(regional_calendars_default, user_calendar_file)
  message(
    "A provisional editable calendar was created. Review it before the final ",
    "climate inference: ", user_calendar_file
  )
}

# Preserve any user edits while adding default rows for regions introduced by
# a revised code lookup. This also makes reruns safe when the editable calendar
# was created by an earlier script version.
regional_calendars <- bind_rows(
  regional_calendars |> mutate(calendar_priority = 1L),
  regional_calendars_default |> mutate(calendar_priority = 2L)
) |>
  arrange(calendar_priority) |>
  distinct(country, parent_region_code, .keep_all = TRUE) |>
  select(-calendar_priority)

month_columns <- c(
  "sowing_start_month", "sowing_end_month",
  "vegetative_start_month", "vegetative_end_month",
  "reproductive_start_month", "reproductive_end_month",
  "harvest_start_month", "harvest_end_month"
)
if (any(!unlist(regional_calendars[month_columns]) %in% 1:12)) {
  stop("All crop-calendar month fields must be integers from 1 to 12.", call. = FALSE)
}


# 8. ATTACH CALENDARS TO SPATIAL UNITS ----------------------------------------

unit_calendar <- centroids |>
  mutate(parent_region_code = standardize_region_code(parent_region_code)) |>
  left_join(
    regional_calendars,
    by = c("country", "parent_region_code")
  ) |>
  mutate(
    download_start_month_previous_year = download_start_month_previous_year,
    download_end_month_harvest_year = download_end_month_harvest_year,
    climate_start_year = climate_start_year,
    climate_end_year = climate_end_year
  )

missing_calendars <- unit_calendar |>
  filter(is.na(calendar_group)) |>
  distinct(
    country, parent_region_code, parent_region_name,
    spatial_unit_id, spatial_unit_name
  )

write_csv(missing_calendars, file.path(output_dir, "QC_units_without_crop_calendar.csv"))
if (nrow(missing_calendars) > 0) {
  stop(
    nrow(missing_calendars),
    " spatial units have no regional calendar. Update ", user_calendar_file,
    " and rerun.", call. = FALSE
  )
}

# One request covers the entire daily record for each point. Script 10 will use
# this table as a restartable download manifest.
request_manifest <- unit_calendar |>
  transmute(
    climate_location_id,
    country,
    spatial_unit_id,
    spatial_unit_name,
    parent_region_code,
    parent_region_name,
    latitude = round(latitude, 5),
    longitude = round(longitude, 5),
    start_date = sprintf("%d%02d01", climate_start_year - 1L, download_start_month_previous_year),
    end_date = sprintf("%d%02d31", climate_end_year, download_end_month_harvest_year),
    reference_harvested_area_ha,
    calendar_group
  ) |>
  arrange(country, climate_location_id)

calendar_summary <- unit_calendar |>
  group_by(
    country, parent_region_code, parent_region_name, calendar_group,
    sowing_window, vegetative_window, reproductive_window, harvest_window,
    calendar_source, calendar_status
  ) |>
  summarise(
    spatial_units = n(),
    reference_harvested_area_ha = sum(reference_harvested_area_ha, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(country, parent_region_name)

coordinate_summary <- unit_calendar |>
  group_by(country) |>
  summarise(
    spatial_units = n(),
    minimum_longitude = min(longitude),
    maximum_longitude = max(longitude),
    minimum_latitude = min(latitude),
    maximum_latitude = max(latitude),
    climate_start_year = first(climate_start_year),
    climate_end_year = first(climate_end_year),
    .groups = "drop"
  )


# 9. SAVE OUTPUTS --------------------------------------------------------------

message("Saving centroid and crop-calendar outputs...")

write_csv(centroids, file.path(output_dir, "soybean_spatial_sampling_points.csv"))
write_csv(regional_calendars, file.path(output_dir, "soybean_regional_crop_calendars.csv"))
write_csv(unit_calendar, file.path(output_dir, "soybean_spatial_units_with_calendars.csv"))
write_csv(request_manifest, file.path(output_dir, "NASA_POWER_soybean_request_manifest.csv"))
write_csv(calendar_summary, file.path(output_dir, "QC_crop_calendar_summary.csv"))
write_csv(coordinate_summary, file.path(output_dir, "QC_spatial_sampling_summary.csv"))

sampling_points_output <- sampling_points |>
  standardize_geometry_column() |>
  select(country, spatial_unit_id = spatial_unit_id_join, geometry) |>
  standardize_geometry_column()

st_write(
  sampling_points_output,
  file.path(output_dir, "soybean_spatial_sampling_points.gpkg"),
  layer = "sampling_points",
  delete_dsn = TRUE,
  quiet = TRUE
)


# 10. FINAL VALIDATION ---------------------------------------------------------

if (nrow(request_manifest) != nrow(centroids)) {
  stop("Request manifest and centroid table have different row counts.", call. = FALSE)
}
if (anyDuplicated(request_manifest$climate_location_id)) {
  stop("Request manifest contains duplicated location identifiers.", call. = FALSE)
}
if (any(is.na(request_manifest$calendar_group))) {
  stop("Request manifest contains missing crop calendars.", call. = FALSE)
}

message("")
message("Script 09 completed successfully.")
message("Spatial sampling units: ", nrow(request_manifest))
message("Brazilian municipalities: ", sum(request_manifest$country == "Brazil"))
message("Argentine departments: ", sum(request_manifest$country == "Argentina"))
message("Climate period prepared: ", climate_start_year, " to ", climate_end_year)
message("Output directory: ", normalizePath(output_dir, winslash = "/"))
message("")
message("IMPORTANT: Review soybean_regional_crop_calendars_user.csv before final inference.")
message("The broad Aug-to-Jul weather download does not depend on the provisional calendar.")

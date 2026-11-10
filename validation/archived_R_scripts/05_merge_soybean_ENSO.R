# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 05_merge_soybean_ENSO.R
# PURPOSE: Merge detrended soybean productivity with ENSO indices and generate
#          descriptive annual, phase composite and crop failure datasets
# VERSION: 1.1
# IMPORTANT: This script prepares descriptive evidence. Formal variance
#            attribution is reserved for script 06.
# ==============================================================================


# 1. USER SETTINGS -------------------------------------------------------------

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")

project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) {
  project_dir_override
} else {
  default_project_dir
}

analysis_start_year <- 1974L
analysis_end_year <- 2024L

# Thresholds used to describe negative productivity shocks.
break_thresholds_percent <- c(-5, -10, -15, -20)

# Reproducible bootstrap confidence intervals based on harvest years.
bootstrap_repetitions <- 5000L
bootstrap_seed <- 20260904L


# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c("dplyr", "readr", "tidyr", "tibble")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})


# 3. PATHS ---------------------------------------------------------------------

detrended_dir <- file.path(project_dir, "data_detrended")
enso_dir <- file.path(project_dir, "climate_indices", "processed")
output_dir <- file.path(project_dir, "data_merged")
metadata_dir <- file.path(project_dir, "metadata")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

soybean_file <- file.path(
  detrended_dir,
  "soybean_Brazil_Argentina_detrended_selected.csv"
)

enso_file <- file.path(
  enso_dir,
  paste0(
    "ENSO_indices_by_soybean_harvest_year_",
    analysis_start_year,
    "_",
    analysis_end_year,
    ".csv"
  )
)

if (!file.exists(soybean_file)) {
  stop("Detrended soybean file not found: ", soybean_file, call. = FALSE)
}

if (!file.exists(enso_file)) {
  stop("ENSO harvest year file not found: ", enso_file, call. = FALSE)
}


# 4. HELPERS -------------------------------------------------------------------

check_required_columns <- function(data, required, dataset_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      dataset_name,
      " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

safe_ratio <- function(numerator, denominator) {
  if_else(
    is.finite(denominator) & denominator > 0,
    numerator / denominator,
    NA_real_
  )
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

bootstrap_mean_interval <- function(
  x,
  repetitions = bootstrap_repetitions,
  confidence = 0.95
) {
  x <- x[is.finite(x)]

  if (length(x) < 2) {
    return(tibble(
      bootstrap_mean = if (length(x) == 1) x else NA_real_,
      bootstrap_ci_lower = NA_real_,
      bootstrap_ci_upper = NA_real_
    ))
  }

  bootstrap_means <- replicate(
    repetitions,
    mean(sample(x, size = length(x), replace = TRUE))
  )

  alpha <- (1 - confidence) / 2

  tibble(
    bootstrap_mean = mean(x),
    bootstrap_ci_lower = unname(quantile(
      bootstrap_means,
      probs = alpha,
      na.rm = TRUE
    )),
    bootstrap_ci_upper = unname(quantile(
      bootstrap_means,
      probs = 1 - alpha,
      na.rm = TRUE
    ))
  )
}


# 5. READ INPUTS ---------------------------------------------------------------

message("Reading detrended soybean and ENSO data...")

soybean <- read_csv(
  soybean_file,
  col_types = cols(
    spatial_unit_id = col_character(),
    parent_region_code = col_character(),
    analysis_year = col_integer(),
    harvested_area_ha = col_double(),
    production_t = col_double(),
    yield_harvested_t_ha = col_double(),
    selected_trend_t_ha = col_double(),
    yield_anomaly_t_ha = col_double(),
    yield_anomaly_percent = col_double(),
    .default = col_guess()
  ),
  na = c("", "NA", "N/A"),
  show_col_types = FALSE
)

enso <- read_csv(
  enso_file,
  col_types = cols(
    harvest_year = col_integer(),
    .default = col_guess()
  ),
  na = c("", "NA", "N/A"),
  show_col_types = FALSE
)

check_required_columns(
  soybean,
  c(
    "country", "spatial_unit_id", "spatial_unit_name",
    "parent_region_name", "spatial_level", "analysis_year",
    "harvested_area_ha", "production_t", "yield_harvested_t_ha",
    "selected_trend_method", "temporal_trend_detected",
    "selected_trend_t_ha", "yield_anomaly_t_ha",
    "yield_anomaly_percent", "yield_anomaly_standardized"
  ),
  "Detrended soybean dataset"
)

check_required_columns(
  enso,
  c(
    "harvest_year", "oni_djf_c", "oni_harvest_year_phase",
    "oni_harvest_year_strength", "nino34_djf_c",
    "nino34_growing_season_oct_apr_c",
    "nino34_reproductive_jan_mar_c",
    "meiv2_reproductive_jan_mar"
  ),
  "ENSO harvest year dataset"
)


# 6. INPUT QUALITY CONTROL -----------------------------------------------------

if (anyDuplicated(
  soybean[c("country", "spatial_unit_id", "analysis_year")]
) > 0) {
  stop("Duplicated soybean spatial unit and year keys were found.", call. = FALSE)
}

if (anyDuplicated(enso["harvest_year"]) > 0) {
  stop("Duplicated ENSO harvest years were found.", call. = FALSE)
}

if (any(!is.finite(soybean$selected_trend_t_ha)) ||
    any(soybean$selected_trend_t_ha <= 0, na.rm = TRUE)) {
  stop("Invalid technological trend values were found.", call. = FALSE)
}

soybean_years <- sort(unique(soybean$analysis_year))
missing_enso_years <- setdiff(soybean_years, enso$harvest_year)

if (length(missing_enso_years) > 0) {
  stop(
    "ENSO data are unavailable for soybean years: ",
    paste(missing_enso_years, collapse = ", "),
    call. = FALSE
  )
}


# 7. MERGE SOYBEAN AND ENSO ----------------------------------------------------

soybean_enso_panel <- soybean |>
  filter(
    analysis_year >= analysis_start_year,
    analysis_year <= analysis_end_year
  ) |>
  left_join(
    enso,
    by = c("analysis_year" = "harvest_year"),
    relationship = "many-to-one"
  ) |>
  mutate(
    oni_harvest_year_phase = factor(
      oni_harvest_year_phase,
      levels = c("La Nina", "Neutral", "El Nino")
    ),
    oni_harvest_year_strength = factor(
      oni_harvest_year_strength,
      levels = c("Neutral", "Weak", "Moderate", "Strong", "Very strong")
    ),
    expected_production_t = selected_trend_t_ha * harvested_area_ha,
    production_anomaly_t = production_t - expected_production_t,
    production_anomaly_percent = 100 * safe_ratio(
      production_anomaly_t,
      expected_production_t
    ),
    break_5_percent = yield_anomaly_percent <= -5,
    break_10_percent = yield_anomaly_percent <= -10,
    break_15_percent = yield_anomaly_percent <= -15,
    break_20_percent = yield_anomaly_percent <= -20,
    break_one_standard_deviation = yield_anomaly_standardized <= -1
  ) |>
  arrange(country, spatial_unit_id, analysis_year)

if (any(is.na(soybean_enso_panel$oni_djf_c))) {
  stop("Missing ONI values were generated after the merge.", call. = FALSE)
}

if (any(is.na(soybean_enso_panel$oni_harvest_year_phase))) {
  stop("Missing ENSO phases were generated after the merge.", call. = FALSE)
}

if (nrow(soybean_enso_panel) != nrow(soybean)) {
  stop("The merge changed the number of soybean observations.", call. = FALSE)
}


# 8. COUNTRY YEAR AGGREGATION --------------------------------------------------

# Expected production is calculated from the local technological trend and the
# observed harvested area. Aggregation therefore measures the productivity
# shock while retaining the spatial distribution of soybean in each year.
country_year <- soybean_enso_panel |>
  group_by(
    country,
    analysis_year,
    oni_harvest_year_phase,
    oni_harvest_year_strength
  ) |>
  summarise(
    active_spatial_units = n_distinct(spatial_unit_id),
    harvested_area_ha = sum(harvested_area_ha, na.rm = TRUE),
    observed_production_t = sum(production_t, na.rm = TRUE),
    expected_production_t = sum(expected_production_t, na.rm = TRUE),
    observed_yield_t_ha = safe_ratio(
      observed_production_t,
      harvested_area_ha
    ),
    expected_yield_t_ha = safe_ratio(
      expected_production_t,
      harvested_area_ha
    ),
    yield_anomaly_t_ha = observed_yield_t_ha - expected_yield_t_ha,
    yield_anomaly_percent = 100 * safe_ratio(
      observed_production_t - expected_production_t,
      expected_production_t
    ),
    spatial_units_break_5_percent = sum(break_5_percent, na.rm = TRUE),
    spatial_units_break_10_percent = sum(break_10_percent, na.rm = TRUE),
    spatial_units_break_15_percent = sum(break_15_percent, na.rm = TRUE),
    spatial_units_break_20_percent = sum(break_20_percent, na.rm = TRUE),
    spatial_exposure_break_10_percent = mean(
      break_10_percent,
      na.rm = TRUE
    ),
    spatial_exposure_break_20_percent = mean(
      break_20_percent,
      na.rm = TRUE
    ),
    across(
      c(
        oni_son_previous_c, oni_ond_previous_c, oni_ndj_c,
        oni_djf_c, oni_jfm_c, oni_fma_c, oni_mam_c,
        nino34_son_previous_c, nino34_ond_previous_c, nino34_ndj_c,
        nino34_djf_c, nino34_jfm_c, nino34_fma_c, nino34_mam_c,
        nino34_growing_season_oct_apr_c,
        nino34_reproductive_jan_mar_c,
        meiv2_growing_season_oct_apr,
        meiv2_reproductive_jan_mar
      ),
      ~ first(.x)
    ),
    .groups = "drop"
  ) |>
  arrange(country, analysis_year)


# 9. SOUTH AMERICA YEAR AGGREGATION -------------------------------------------

south_america_year <- soybean_enso_panel |>
  group_by(
    analysis_year,
    oni_harvest_year_phase,
    oni_harvest_year_strength
  ) |>
  summarise(
    countries = n_distinct(country),
    active_spatial_units = n_distinct(interaction(country, spatial_unit_id)),
    harvested_area_ha = sum(harvested_area_ha, na.rm = TRUE),
    observed_production_t = sum(production_t, na.rm = TRUE),
    expected_production_t = sum(expected_production_t, na.rm = TRUE),
    observed_yield_t_ha = safe_ratio(
      observed_production_t,
      harvested_area_ha
    ),
    expected_yield_t_ha = safe_ratio(
      expected_production_t,
      harvested_area_ha
    ),
    yield_anomaly_t_ha = observed_yield_t_ha - expected_yield_t_ha,
    yield_anomaly_percent = 100 * safe_ratio(
      observed_production_t - expected_production_t,
      expected_production_t
    ),
    spatial_exposure_break_10_percent = mean(
      break_10_percent,
      na.rm = TRUE
    ),
    spatial_exposure_break_20_percent = mean(
      break_20_percent,
      na.rm = TRUE
    ),
    across(
      c(
        oni_son_previous_c, oni_ond_previous_c, oni_ndj_c,
        oni_djf_c, oni_jfm_c, oni_fma_c, oni_mam_c,
        nino34_son_previous_c, nino34_ond_previous_c, nino34_ndj_c,
        nino34_djf_c, nino34_jfm_c, nino34_fma_c, nino34_mam_c,
        nino34_growing_season_oct_apr_c,
        nino34_reproductive_jan_mar_c,
        meiv2_growing_season_oct_apr,
        meiv2_reproductive_jan_mar
      ),
      ~ first(.x)
    ),
    .groups = "drop"
  ) |>
  mutate(country = "Brazil and Argentina", .before = 1) |>
  arrange(analysis_year)


# 10. ENSO PHASE COMPOSITES ----------------------------------------------------

annual_all_geographies <- bind_rows(
  country_year,
  south_america_year
)

set.seed(bootstrap_seed)

phase_composite_basic <- annual_all_geographies |>
  group_by(country, oni_harvest_year_phase) |>
  summarise(
    harvest_years = n(),
    first_year = min(analysis_year),
    last_year = max(analysis_year),
    mean_yield_anomaly_percent = mean(yield_anomaly_percent, na.rm = TRUE),
    median_yield_anomaly_percent = median(yield_anomaly_percent, na.rm = TRUE),
    sd_yield_anomaly_percent = sd(yield_anomaly_percent, na.rm = TRUE),
    minimum_yield_anomaly_percent = min(yield_anomaly_percent, na.rm = TRUE),
    maximum_yield_anomaly_percent = max(yield_anomaly_percent, na.rm = TRUE),
    mean_yield_anomaly_t_ha = mean(yield_anomaly_t_ha, na.rm = TRUE),
    mean_spatial_exposure_break_10_percent = 100 * mean(
      spatial_exposure_break_10_percent,
      na.rm = TRUE
    ),
    mean_spatial_exposure_break_20_percent = 100 * mean(
      spatial_exposure_break_20_percent,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

phase_composite_bootstrap <- annual_all_geographies |>
  group_by(country, oni_harvest_year_phase) |>
  group_modify(~ bootstrap_mean_interval(.x$yield_anomaly_percent)) |>
  ungroup()

phase_composite <- phase_composite_basic |>
  left_join(
    phase_composite_bootstrap,
    by = c("country", "oni_harvest_year_phase")
  ) |>
  arrange(country, oni_harvest_year_phase)


# 11. REGIONAL BREAK PROBABILITY ----------------------------------------------

# Each observation in this table is an independent harvest year within a
# geography. Probability is the empirical fraction of years below a threshold.
regional_break_probability <- annual_all_geographies |>
  select(
    country,
    analysis_year,
    oni_harvest_year_phase,
    yield_anomaly_percent
  ) |>
  crossing(break_threshold_percent = break_thresholds_percent) |>
  group_by(country, oni_harvest_year_phase, break_threshold_percent) |>
  summarise(
    harvest_years = n_distinct(analysis_year),
    break_years = sum(
      yield_anomaly_percent <= break_threshold_percent,
      na.rm = TRUE
    ),
    empirical_break_probability = break_years / harvest_years,
    empirical_break_probability_percent = 100 * empirical_break_probability,
    .groups = "drop"
  ) |>
  arrange(country, oni_harvest_year_phase, break_threshold_percent)


# 12. SPATIAL EXPOSURE BY ENSO PHASE ------------------------------------------

# This table describes the fraction of available unit year observations with a
# local productivity anomaly below each threshold. It is spatial exposure and
# must not be interpreted as the probability of a national production failure.
spatial_exposure_by_phase <- soybean_enso_panel |>
  group_by(country, oni_harvest_year_phase) |>
  summarise(
    spatial_unit_year_observations = n(),
    spatial_units = n_distinct(spatial_unit_id),
    harvest_years = n_distinct(analysis_year),
    exposure_break_5_percent = mean(break_5_percent, na.rm = TRUE),
    exposure_break_10_percent = mean(break_10_percent, na.rm = TRUE),
    exposure_break_15_percent = mean(break_15_percent, na.rm = TRUE),
    exposure_break_20_percent = mean(break_20_percent, na.rm = TRUE),
    exposure_break_one_standard_deviation = mean(
      break_one_standard_deviation,
      na.rm = TRUE
    ),
    mean_local_anomaly_percent = mean(yield_anomaly_percent, na.rm = TRUE),
    median_local_anomaly_percent = median(yield_anomaly_percent, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    across(
      starts_with("exposure_"),
      ~ 100 * .x,
      .names = "{.col}_percent"
    )
  )


# 13. SPATIAL UNIT PHASE SUMMARY ----------------------------------------------

spatial_unit_phase_summary <- soybean_enso_panel |>
  group_by(
    country,
    spatial_unit_id,
    spatial_unit_name,
    parent_region_name,
    spatial_level,
    oni_harvest_year_phase
  ) |>
  summarise(
    observations = n(),
    first_year = min(analysis_year),
    last_year = max(analysis_year),
    mean_yield_anomaly_t_ha = mean(yield_anomaly_t_ha, na.rm = TRUE),
    mean_yield_anomaly_percent = mean(yield_anomaly_percent, na.rm = TRUE),
    median_yield_anomaly_percent = median(yield_anomaly_percent, na.rm = TRUE),
    sd_yield_anomaly_percent = sd(yield_anomaly_percent, na.rm = TRUE),
    break_10_percent_observations = sum(break_10_percent, na.rm = TRUE),
    break_20_percent_observations = sum(break_20_percent, na.rm = TRUE),
    probability_break_10_percent = mean(break_10_percent, na.rm = TRUE),
    probability_break_20_percent = mean(break_20_percent, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(country, spatial_unit_id, oni_harvest_year_phase)


# 14. EXPLORATORY SEASONAL CORRELATIONS ---------------------------------------

seasonal_index_columns <- c(
  "oni_son_previous_c", "oni_ond_previous_c", "oni_ndj_c",
  "oni_djf_c", "oni_jfm_c", "oni_fma_c", "oni_mam_c",
  "nino34_son_previous_c", "nino34_ond_previous_c", "nino34_ndj_c",
  "nino34_djf_c", "nino34_jfm_c", "nino34_fma_c", "nino34_mam_c",
  "nino34_growing_season_oct_apr_c",
  "nino34_reproductive_jan_mar_c",
  "meiv2_growing_season_oct_apr",
  "meiv2_reproductive_jan_mar"
)

exploratory_seasonal_correlations <- bind_rows(lapply(
  split(annual_all_geographies, annual_all_geographies$country),
  function(country_data) {
    bind_rows(lapply(seasonal_index_columns, function(index_name) {
      valid <- complete.cases(
        country_data$yield_anomaly_percent,
        country_data[[index_name]]
      )
      tibble(
        country = country_data$country[1],
        enso_index_window = index_name,
        harvest_years = sum(valid),
        pearson_correlation = if (sum(valid) >= 10) {
          cor(
            country_data$yield_anomaly_percent[valid],
            country_data[[index_name]][valid]
          )
        } else {
          NA_real_
        },
        spearman_correlation = if (sum(valid) >= 10) {
          cor(
            country_data$yield_anomaly_percent[valid],
            country_data[[index_name]][valid],
            method = "spearman"
          )
        } else {
          NA_real_
        }
      )
    }))
  }
)) |>
  arrange(country, desc(abs(pearson_correlation)))


# 15. MERGE QUALITY CONTROL ----------------------------------------------------

qc_merge_summary <- soybean_enso_panel |>
  group_by(country) |>
  summarise(
    observations = n(),
    spatial_units = n_distinct(spatial_unit_id),
    first_year = min(analysis_year),
    last_year = max(analysis_year),
    missing_oni = sum(is.na(oni_djf_c)),
    missing_nino34 = sum(is.na(nino34_djf_c)),
    missing_meiv2 = sum(is.na(meiv2_reproductive_jan_mar)),
    invalid_expected_production = sum(
      !is.finite(expected_production_t) | expected_production_t <= 0
    ),
    duplicated_keys = 0L,
    .groups = "drop"
  )

qc_year_coverage <- soybean_enso_panel |>
  group_by(country, analysis_year) |>
  summarise(
    spatial_units = n_distinct(spatial_unit_id),
    harvested_area_ha = sum(harvested_area_ha, na.rm = TRUE),
    observed_production_t = sum(production_t, na.rm = TRUE),
    expected_production_t = sum(expected_production_t, na.rm = TRUE),
    oni_harvest_year_phase = first(oni_harvest_year_phase),
    oni_djf_c = first(oni_djf_c),
    .groups = "drop"
  ) |>
  arrange(country, analysis_year)

known_enso_years <- tibble(
  analysis_year = c(1983L, 1989L, 1998L, 2010L, 2012L, 2016L, 2023L, 2024L),
  expected_phase = c(
    "El Nino", "La Nina", "El Nino", "El Nino",
    "La Nina", "El Nino", "La Nina", "El Nino"
  )
)

qc_known_enso_years <- enso |>
  filter(harvest_year %in% known_enso_years$analysis_year) |>
  select(
    analysis_year = harvest_year,
    observed_phase = oni_harvest_year_phase,
    observed_strength = oni_harvest_year_strength,
    oni_djf_c
  ) |>
  left_join(known_enso_years, by = "analysis_year") |>
  mutate(phase_matches_expectation = observed_phase == expected_phase)

if (!all(qc_known_enso_years$phase_matches_expectation)) {
  stop("Known ENSO years failed the phase alignment check.", call. = FALSE)
}


# 16. METADATA -----------------------------------------------------------------

analysis_metadata <- tibble(
  parameter = c(
    "analysis_start_year",
    "analysis_end_year",
    "break_thresholds_percent",
    "bootstrap_repetitions",
    "bootstrap_seed",
    "primary_enso_index",
    "primary_enso_window",
    "aggregation_weight"
  ),
  value = c(
    as.character(analysis_start_year),
    as.character(analysis_end_year),
    paste(break_thresholds_percent, collapse = ", "),
    as.character(bootstrap_repetitions),
    as.character(bootstrap_seed),
    "ONI",
    "DJF aligned with harvest year",
    "Harvested area through expected production"
  ),
  description = c(
    "First soybean harvest year",
    "Last soybean harvest year",
    "Regional and local productivity break thresholds",
    "Number of bootstrap resamples of harvest years",
    "Seed used for reproducible bootstrap intervals",
    "Primary ENSO index for classification and later attribution",
    "Primary seasonal ENSO window",
    "Weight used to aggregate local productivity anomalies"
  )
)


# 17. SAVE OUTPUTS -------------------------------------------------------------

message("Saving merged and descriptive outputs...")

write_csv(
  soybean_enso_panel,
  file.path(output_dir, "soybean_ENSO_spatial_panel_1974_2024.csv"),
  na = ""
)
write_csv(
  country_year,
  file.path(output_dir, "soybean_ENSO_country_year_1974_2024.csv"),
  na = ""
)
write_csv(
  south_america_year,
  file.path(output_dir, "soybean_ENSO_Brazil_Argentina_year_1974_2024.csv"),
  na = ""
)
write_csv(
  phase_composite,
  file.path(output_dir, "descriptive_ENSO_phase_composite.csv"),
  na = ""
)
write_csv(
  regional_break_probability,
  file.path(output_dir, "descriptive_regional_break_probability.csv"),
  na = ""
)
write_csv(
  spatial_exposure_by_phase,
  file.path(output_dir, "descriptive_spatial_exposure_by_ENSO_phase.csv"),
  na = ""
)
write_csv(
  spatial_unit_phase_summary,
  file.path(output_dir, "descriptive_spatial_unit_phase_summary.csv"),
  na = ""
)
write_csv(
  exploratory_seasonal_correlations,
  file.path(output_dir, "exploratory_ENSO_seasonal_correlations.csv"),
  na = ""
)
write_csv(
  qc_merge_summary,
  file.path(output_dir, "QC_merge_summary.csv"),
  na = ""
)
write_csv(
  qc_year_coverage,
  file.path(output_dir, "QC_merge_year_coverage.csv"),
  na = ""
)
write_csv(
  qc_known_enso_years,
  file.path(output_dir, "QC_known_ENSO_years.csv"),
  na = ""
)
write_csv(
  analysis_metadata,
  file.path(metadata_dir, "merge_soybean_ENSO_parameters.csv"),
  na = ""
)
capture.output(
  sessionInfo(),
  file = file.path(metadata_dir, "merge_soybean_ENSO_R_session_info.txt")
)


# 18. FINAL REPORT -------------------------------------------------------------

message("")
message("Soybean and ENSO merge completed successfully.")
message("Spatial observations: ", nrow(soybean_enso_panel))
message("Country year observations: ", nrow(country_year))
message("Combined annual observations: ", nrow(south_america_year))
message("Output directory: ", normalizePath(output_dir, winslash = "/"))
message("")

print(qc_merge_summary)
print(phase_composite)
print(regional_break_probability)

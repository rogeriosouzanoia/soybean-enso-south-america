#!/usr/bin/env Rscript

# Figure 3 data preparation
# Seven climate indices connecting ENSO and soybean yield.
#
# Regional dominant drivers are derived from standardized hierarchical ridge
# coefficients. Robustness describes stability across outer validation folds.
# It is not interpreted as a local significance test.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
})

options(dplyr.summarise.inform = FALSE)

args <- commandArgs(trailingOnly = TRUE)
default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else {
  override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
  normalizePath(if (nzchar(override)) override else default_project_dir,
                mustWork = TRUE)
}
output_dir <- if (length(args) >= 2L) args[[2]] else
  file.path(project_dir, "outputs", "figure_data", "figure_3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

find_one <- function(root, candidates) {
  direct <- file.path(root, candidates)
  if (any(file.exists(direct))) return(normalizePath(direct[file.exists(direct)][[1]]))
  all_files <- list.files(root, recursive = TRUE, full.names = TRUE)
  hits <- all_files[basename(all_files) %in% candidates]
  if (length(hits) == 0L) {
    stop("Could not find any of: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  normalizePath(hits[order(nchar(hits), hits)][[1]], mustWork = TRUE)
}

standardize_id <- function(x, country) {
  digits <- gsub("[^0-9]", "", as.character(x))
  widths <- ifelse(country == "Brazil", 7L, 5L)
  mapply(function(value, width) {
    if (is.na(value) || value == "") return(NA_character_)
    sprintf(paste0("%0", width, "d"), as.integer(value))
  }, digits, widths, USE.NAMES = FALSE)
}

classify_macroregion <- function(country, parent_region) {
  region <- str_to_upper(str_squish(as.character(parent_region)))
  case_when(
    country == "Argentina" ~ "Argentina",
    country == "Brazil" & region %in% c("PR", "SC", "RS") ~ "Southern Brazil",
    country == "Brazil" & region %in% c(
      "TO", "MA", "PI", "BA", "AL", "CE", "PB", "PE", "RN", "SE"
    ) ~ "MATOPIBA and northeastern frontier",
    country == "Brazil" ~ "Central Brazil",
    TRUE ~ "Outside focal regions"
  )
}

predictor_dictionary <- tribble(
  ~predictor, ~driver_label,
  "climate_z_mean_temperature_c", "Temperature",
  "climate_z_mean_vpd_kpa", "VPD",
  "climate_z_precipitation_total_mm", "Precipitation",
  "climate_z_heavy_rain_fraction", "Heavy rainfall",
  "water_z_crop_water_stress", "CWS",
  "climate_z_hot_dry_fraction", "Hot-dry conditions",
  "climate_z_hot_wet_fraction", "Hot-wet conditions"
)

paths <- list(
  coefficients = find_one(project_dir, c(
    "ENSO_local_climate_regional_coefficients.csv.gz",
    "checkpoint_12_v6_coefficients.csv.gz"
  )),
  importance = find_one(project_dir, c(
    "ENSO_climate_variable_permutation_importance.csv",
    "ENSO_climate_variable_permutation_importance(3).csv"
  )),
  mediated = find_one(project_dir, c(
    "ENSO_climate_mediated_predictive_scores.csv",
    "ENSO_climate_mediated_predictive_scores(3).csv"
  )),
  enso_skill = find_one(project_dir, c(
    "ENSO_to_climate_crossvalidated_skill.csv",
    "ENSO_to_climate_crossvalidated_skill(3).csv"
  )),
  circular = find_one(project_dir, c(
    "ENSO_climate_circular_shift_tests.csv",
    "ENSO_climate_circular_shift_tests(3).csv"
  )),
  regional_models = find_one(project_dir, c(
    "ENSO_climate_attribution_region_summary.csv",
    "ENSO_climate_attribution_region_summary(3).csv"
  )),
  country_models = find_one(project_dir, c(
    "ENSO_climate_water_model_comparison.csv",
    "ENSO_climate_water_model_comparison(2).csv"
  )),
  local_units = find_one(project_dir, c(
    "Figure_2_local_ENSO_metrics.csv"
  ))
)

coefficients <- read_csv(paths$coefficients, show_col_types = FALSE)
importance <- read_csv(paths$importance, show_col_types = FALSE)
mediated <- read_csv(paths$mediated, show_col_types = FALSE)
enso_skill <- read_csv(paths$enso_skill, show_col_types = FALSE)
circular <- read_csv(paths$circular, show_col_types = FALSE)
regional_models <- read_csv(paths$regional_models, show_col_types = FALSE)
country_models <- read_csv(paths$country_models, show_col_types = FALSE)
local_units <- read_csv(
  paths$local_units,
  col_types = cols(
    spatial_unit_id = col_character(),
    spatial_unit_id_join = col_character(), .default = col_guess()
  ), show_col_types = FALSE
)

required_coefficient_columns <- c(
  "country", "model", "outer_fold", "parent_region_name", "predictor",
  "effective_coefficient"
)
missing_coefficient_columns <- setdiff(required_coefficient_columns, names(coefficients))
if (length(missing_coefficient_columns) > 0L) {
  stop(
    "Coefficient file is missing: ",
    paste(missing_coefficient_columns, collapse = ", "), call. = FALSE
  )
}

# Standardized coefficients are comparable within the fitted model. The
# dominant driver is selected independently in each outer fold and stability is
# then summarized across folds.
coefficient_fold <- coefficients |>
  filter(model == "climate_plus_water", predictor %in% predictor_dictionary$predictor) |>
  left_join(predictor_dictionary, by = "predictor") |>
  group_by(country, parent_region_name, outer_fold) |>
  mutate(
    absolute_coefficient = abs(effective_coefficient),
    fold_maximum = max(absolute_coefficient, na.rm = TRUE),
    dominant_in_fold = absolute_coefficient == fold_maximum
  ) |>
  ungroup()

regional_driver_summary <- coefficient_fold |>
  group_by(
    country, parent_region_name, predictor, driver_label
  ) |>
  summarise(
    folds = n_distinct(outer_fold),
    median_coefficient = median(effective_coefficient, na.rm = TRUE),
    median_absolute_coefficient = median(absolute_coefficient, na.rm = TRUE),
    lower_coefficient = quantile(effective_coefficient, 0.10, na.rm = TRUE,
                                 names = FALSE, type = 8),
    upper_coefficient = quantile(effective_coefficient, 0.90, na.rm = TRUE,
                                 names = FALSE, type = 8),
    positive_fold_fraction = mean(effective_coefficient > 0, na.rm = TRUE),
    dominance_fraction = mean(dominant_in_fold, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    sign_stability = pmax(positive_fold_fraction, 1 - positive_fold_fraction),
    robust_sign = sign_stability >= 0.80,
    direction = case_when(
      median_coefficient > 0 ~ "Positive",
      median_coefficient < 0 ~ "Negative",
      TRUE ~ "Neutral"
    )
  )

dominant_regional_driver <- regional_driver_summary |>
  group_by(country, parent_region_name) |>
  arrange(desc(dominance_fraction), desc(median_absolute_coefficient), driver_label) |>
  slice(1L) |>
  ungroup() |>
  mutate(
    robust_dominance = dominance_fraction >= 0.60 & robust_sign,
    evidence_class = case_when(
      robust_dominance & sign_stability >= 1.00 ~ "Very robust",
      robust_dominance ~ "Robust",
      dominance_fraction >= 0.40 ~ "Moderate",
      TRUE ~ "Uncertain"
    )
  )

local_driver_map <- local_units |>
  transmute(
    country,
    spatial_unit_id = as.character(spatial_unit_id),
    spatial_unit_id_join = standardize_id(spatial_unit_id_join, country),
    spatial_unit_name,
    parent_region_name,
    macroregion,
    reference_area_ha,
    longitude,
    latitude
  ) |>
  left_join(
    dominant_regional_driver,
    by = c("country", "parent_region_name")
  )

if (mean(!is.na(local_driver_map$driver_label)) < 0.95) {
  warning("Less than 95% of local units received a regional dominant driver.")
}

country_index_summary <- importance |>
  filter(predictor %in% predictor_dictionary$predictor) |>
  left_join(predictor_dictionary, by = "predictor") |>
  select(
    country, predictor, driver_label,
    importance_relative_to_null_error,
    positive_importance_probability,
    normalized_positive_importance
  ) |>
  left_join(
    mediated |>
      filter(predictor %in% predictor_dictionary$predictor) |>
      select(
        country, predictor, enso_error_reduction_relative_to_null,
        positive_ENSO_climate_skill,
        ENSO_aligned_predictive_score,
        normalized_ENSO_aligned_predictive_score
      ),
    by = c("country", "predictor")
  ) |>
  mutate(
    yield_importance_percent = 100 * importance_relative_to_null_error,
    ENSO_index_skill_percent = 100 * enso_error_reduction_relative_to_null,
    aligned_score_percent = 100 * normalized_ENSO_aligned_predictive_score
  )

enso_climate_evidence <- enso_skill |>
  left_join(predictor_dictionary, by = "predictor") |>
  left_join(
    circular |>
      select(
        country, predictor, observed_partial_r_squared,
        circular_shift_p_value, null_95th_percentile
      ),
    by = c("country", "predictor")
  ) |>
  mutate(
    ENSO_error_reduction_percent = 100 * enso_error_reduction_relative_to_null,
    partial_r_squared_percent = 100 * observed_partial_r_squared,
    null_95th_percentile_percent = 100 * null_95th_percentile,
    passes_circular_shift_test = circular_shift_p_value <= 0.05,
    positive_crossvalidated_skill = enso_error_reduction_relative_to_null > 0
  )

regional_model_summary <- regional_models |>
  mutate(
    model = replace_na(model, "null"),
    macroregion = classify_macroregion(country, parent_region_name),
    error_reduction_percent = 100 * error_reduction_relative_to_local_null,
    row_observations = observations
  ) |>
  filter(model %in% c(
    "ENSO", "climate", "water_balance", "climate_plus_water",
    "ENSO_plus_climate", "ENSO_plus_climate_water"
  )) |>
  group_by(macroregion, model) |>
  summarise(
    regions = n_distinct(parent_region_name),
    observations = sum(row_observations, na.rm = TRUE),
    weighted_error_reduction_percent = weighted.mean(
      error_reduction_percent, pmax(row_observations, 1), na.rm = TRUE
    ),
    median_error_reduction_percent = median(error_reduction_percent, na.rm = TRUE),
    lower_region_percent = quantile(error_reduction_percent, 0.10, na.rm = TRUE,
                                    names = FALSE, type = 8),
    upper_region_percent = quantile(error_reduction_percent, 0.90, na.rm = TRUE,
                                    names = FALSE, type = 8),
    .groups = "drop"
  )

country_model_summary <- country_models |>
  select(country, starts_with("out_of_sample_r_squared_"),
         ENSO_gain_beyond_combined_mechanisms) |>
  pivot_longer(
    cols = starts_with("out_of_sample_r_squared_"),
    names_prefix = "out_of_sample_r_squared_",
    names_to = "model", values_to = "out_of_sample_r_squared"
  ) |>
  mutate(out_of_sample_r_squared_percent = 100 * out_of_sample_r_squared)

qc_summary <- tribble(
  ~test, ~value,
  "Regional coefficient rows", nrow(regional_driver_summary),
  "Parent regions", n_distinct(interaction(
    dominant_regional_driver$country,
    dominant_regional_driver$parent_region_name
  )),
  "Outer folds", n_distinct(coefficients$outer_fold),
  "Local productive units", nrow(local_driver_map),
  "Local units with mapped driver", sum(!is.na(local_driver_map$driver_label)),
  "Fraction of local units with mapped driver",
    mean(!is.na(local_driver_map$driver_label)),
  "Robust regional dominance fraction",
    mean(dominant_regional_driver$robust_dominance, na.rm = TRUE),
  "ENSO climate links passing circular shift test",
    sum(enso_climate_evidence$passes_circular_shift_test, na.rm = TRUE)
)

methods <- tribble(
  ~item, ~value,
  "figure_question", "Which of seven harmonized climate indices connect ENSO to soybean yield?",
  "dominant_driver_scale", "State or province within the hierarchical ridge model",
  "map_support", "Dominant regional driver displayed only across local productive areas",
  "coefficient_comparability", "All climate predictors were standardized before fitting",
  "dominant_driver_definition", "Highest fold frequency as the largest absolute standardized coefficient",
  "sign_stability", "Fraction of outer folds agreeing with the dominant coefficient sign",
  "robust_dominance", "Dominance in at least 60% of folds and sign agreement in at least 80%",
  "uncertain_areas", "Retained and explicitly marked as uncertain",
  "yield_link", "Held out MSE increase after permutation of one harmonized climate index",
  "ENSO_climate_link", "Out of sample error reduction when ENSO predicts the same climate index",
  "temporal_test", "Exact nonzero circular shifts of the ENSO sequence",
  "causal_language", "All quantities are predictive associations and are not causal mediation estimates"
)

write_csv(
  dominant_regional_driver,
  file.path(output_dir, "Figure_3_dominant_regional_climate_driver.csv"), na = ""
)
write_csv(
  local_driver_map,
  file.path(output_dir, "Figure_3_local_dominant_driver_map.csv"), na = ""
)
write_csv(
  regional_driver_summary,
  file.path(output_dir, "Figure_3_regional_driver_stability.csv"), na = ""
)
write_csv(
  country_index_summary,
  file.path(output_dir, "Figure_3_country_index_summary.csv"), na = ""
)
write_csv(
  enso_climate_evidence,
  file.path(output_dir, "Figure_3_ENSO_to_climate_evidence.csv"), na = ""
)
write_csv(
  regional_model_summary,
  file.path(output_dir, "Figure_3_macroregion_model_performance.csv"), na = ""
)
write_csv(
  country_model_summary,
  file.path(output_dir, "Figure_3_country_model_performance.csv"), na = ""
)
write_csv(
  qc_summary,
  file.path(output_dir, "QC_Figure_3_data_preparation.csv"), na = ""
)
write_csv(
  methods,
  file.path(output_dir, "README_Figure_3_methods.csv"), na = ""
)

message("Figure 3 data written to: ", normalizePath(output_dir))

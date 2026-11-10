#!/usr/bin/env Rscript

# Figure 2 data preparation
# ENSO contribution to spatial soybean productivity variability
#
# Predictive contribution is defined as out of sample error reduction relative
# to the null model. Signal energy is retained as a distinct measure of ENSO
# footprint magnitude and is never labelled as causal explained variance.

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
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else
  Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
output_dir <- if (length(args) >= 2L) args[[2]] else
  file.path(project_dir, "outputs", "figure_data", "figure_2")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_repetitions <- 1999L
bootstrap_block_length <- 3L
bootstrap_seed <- 1501L

find_one <- function(root, candidates) {
  all_files <- list.files(root, recursive = TRUE, full.names = TRUE)
  hits <- all_files[basename(all_files) %in% candidates]
  if (length(hits) == 0L) {
    stop("Could not find any of: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  hits <- hits[order(nchar(hits), hits)]
  normalizePath(hits[[1]], mustWork = TRUE)
}

paths <- list(
  local = find_one(project_dir, c(
    "ENSO_spatially_resolved_local_summary.csv",
    "ENSO_spatially_resolved_local_summary(1).csv"
  )),
  oof = find_one(project_dir, c(
    "ENSO_multivariate_field_oof_predictions.csv"
  )),
  country = find_one(project_dir, c(
    "ENSO_multivariate_field_country_summary.csv"
  )),
  phase = find_one(project_dir, c(
    "ENSO_multivariate_field_phase_summary.csv"
  )),
  randomization = find_one(project_dir, c(
    "ENSO_multivariate_field_circular_shift_test.csv"
  )),
  area = find_one(project_dir, c(
    "soybean_spatial_units_with_calendars.csv"
  )),
  consistency = find_one(project_dir, c(
    "Figure_1_crossvalidated_consistency.csv"
  )),
  cancellation = find_one(project_dir, c(
    "Figure_1_regional_cancellation.csv"
  ))
)

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

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  if (!any(ok) || sum(w[ok]) <= 0) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_sse <- function(y, pred, w) {
  ok <- is.finite(y) & is.finite(pred) & is.finite(w) & w >= 0
  if (!any(ok) || sum(w[ok]) <= 0) return(NA_real_)
  sum(w[ok] * (y[ok] - pred[ok])^2)
}

predictive_r2 <- function(y, pred, null_pred, w) {
  denominator <- weighted_sse(y, null_pred, w)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  1 - weighted_sse(y, pred, w) / denominator
}

sample_blocks <- function(n, block_length) {
  starts <- sample.int(n, ceiling(n / block_length), replace = TRUE)
  ids <- unlist(lapply(starts, function(s) {
    ((s - 1L + seq_len(block_length) - 1L) %% n) + 1L
  }))
  ids[seq_len(n)]
}

summarise_field <- function(x) {
  observed_energy <- weighted_mean_safe(
    (x$yield_anomaly_percent - x$null_prediction_percent)^2,
    x$score_weight
  )
  signal_energy <- weighted_mean_safe(
    x$enso_component_percent^2, x$score_weight
  )
  signed <- weighted_mean_safe(x$enso_component_percent, x$score_weight)
  absolute <- weighted_mean_safe(abs(x$enso_component_percent), x$score_weight)
  tibble(
    observations = nrow(x),
    years = n_distinct(x$analysis_year),
    regions = n_distinct(x$region_key),
    cv_r_squared = predictive_r2(
      x$yield_anomaly_percent, x$model_prediction_percent,
      x$null_prediction_percent, x$score_weight
    ),
    signal_energy_fraction = if_else(
      observed_energy > 0, signal_energy / observed_energy, NA_real_
    ),
    mean_signed_component_percent = signed,
    mean_absolute_component_percent = absolute,
    cancellation_index = if_else(
      absolute > 0, 1 - abs(signed) / absolute, NA_real_
    )
  )
}

bootstrap_field <- function(x) {
  yearly <- x |>
    group_by(analysis_year) |>
    summarise(
      sse_model = sum(score_weight * squared_error_model, na.rm = TRUE),
      sse_null = sum(score_weight * squared_error_null, na.rm = TRUE),
      signal = sum(score_weight * enso_component_percent^2, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(analysis_year)

  if (nrow(yearly) < 8L) {
    return(tibble(
      cv_r2_lower = NA_real_, cv_r2_upper = NA_real_,
      signal_lower = NA_real_, signal_upper = NA_real_
    ))
  }

  draws <- replicate(bootstrap_repetitions, {
    z <- yearly[sample_blocks(nrow(yearly), bootstrap_block_length), ]
    c(
      r2 = 1 - sum(z$sse_model) / sum(z$sse_null),
      energy = sum(z$signal) / sum(z$sse_null)
    )
  })

  tibble(
    cv_r2_lower = quantile(draws["r2", ], 0.025, na.rm = TRUE),
    cv_r2_upper = quantile(draws["r2", ], 0.975, na.rm = TRUE),
    signal_lower = quantile(draws["energy", ], 0.025, na.rm = TRUE),
    signal_upper = quantile(draws["energy", ], 0.975, na.rm = TRUE)
  )
}

message("Reading spatial and field attribution outputs...")

local_raw <- read_csv(
  paths$local,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
)

area <- read_csv(
  paths$area,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
) |>
  transmute(
    country,
    spatial_unit_id_join = standardize_id(spatial_unit_id, country),
    reference_area_ha = as.numeric(reference_harvested_area_ha),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  ) |>
  distinct(country, spatial_unit_id_join, .keep_all = TRUE)

consistency <- read_csv(
  paths$consistency,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
) |>
  mutate(
    spatial_unit_id_join = standardize_id(spatial_unit_id, country),
    phase = if_else(str_detect(phase, "El"), "El Niño", "La Niña")
  ) |>
  select(
    country, spatial_unit_id_join, phase, evaluated_events,
    sign_consistency, dominant_expected_sign
  )

local_metrics <- local_raw |>
  mutate(
    spatial_unit_id_join = standardize_id(spatial_unit_id, country),
    macroregion = classify_macroregion(country, parent_region_name),
    positive_predictive_skill = cv_r_squared > 0,
    cv_r_squared_display_percent = pmax(0, cv_r_squared_percent),
    predictive_skill_class = cut(
      cv_r_squared_display_percent,
      breaks = c(-Inf, 0, 2, 5, 10, 20, Inf),
      labels = c("No positive skill", "0 to 2%", "2 to 5%", "5 to 10%",
                 "10 to 20%", "Above 20%"),
      right = TRUE
    ),
    signal_energy_class = cut(
      cv_enso_signal_energy_percent,
      breaks = c(-Inf, 1, 2, 5, 10, 20, Inf),
      labels = c("Below 1%", "1 to 2%", "2 to 5%", "5 to 10%",
                 "10 to 20%", "Above 20%"),
      right = TRUE
    )
  ) |>
  left_join(area, by = c("country", "spatial_unit_id_join"))

local_consistency_relationship <- local_metrics |>
  select(
    country, spatial_unit_id, spatial_unit_name, parent_region_name,
    spatial_level, spatial_unit_id_join, macroregion, observations,
    cv_r_squared_percent, cv_r_squared_display_percent,
    cv_enso_signal_energy_percent, mean_absolute_component_percent,
    cancellation_index, reference_area_ha, longitude, latitude
  ) |>
  inner_join(consistency, by = c("country", "spatial_unit_id_join"))

local_macroregion_summary <- local_metrics |>
  filter(is.finite(reference_area_ha), reference_area_ha > 0) |>
  group_by(macroregion) |>
  summarise(
    country = first(country),
    spatial_units = n_distinct(spatial_unit_id_join),
    productive_area_ha = sum(reference_area_ha),
    area_with_positive_skill_ha = sum(
      reference_area_ha[positive_predictive_skill %in% TRUE], na.rm = TRUE
    ),
    area_with_positive_skill_percent = 100 * area_with_positive_skill_ha /
      productive_area_ha,
    area_weighted_cv_r2_percent = weighted_mean_safe(
      cv_r_squared_percent, reference_area_ha
    ),
    area_weighted_positive_cv_r2_percent = weighted_mean_safe(
      cv_r_squared_display_percent, reference_area_ha
    ),
    area_weighted_signal_energy_percent = weighted_mean_safe(
      cv_enso_signal_energy_percent, reference_area_ha
    ),
    area_weighted_absolute_component_percent = weighted_mean_safe(
      mean_absolute_component_percent, reference_area_ha
    ),
    .groups = "drop"
  )

oof <- read_csv(paths$oof, show_col_types = FALSE) |>
  mutate(
    macroregion = classify_macroregion(country, parent_region_name),
    phase = case_when(
      str_detect(str_to_lower(phase_simple), "nino") ~ "El Niño",
      str_detect(str_to_lower(phase_simple), "nina") ~ "La Niña",
      TRUE ~ "Neutral"
    )
  )

set.seed(bootstrap_seed)
macroregion_field_summary <- oof |>
  group_by(macroregion) |>
  group_modify(~ bind_cols(summarise_field(.x), bootstrap_field(.x))) |>
  ungroup() |>
  mutate(
    cv_r_squared_percent = 100 * cv_r_squared,
    cv_r2_lower_percent = 100 * cv_r2_lower,
    cv_r2_upper_percent = 100 * cv_r2_upper,
    signal_energy_percent = 100 * signal_energy_fraction,
    signal_lower_percent = 100 * signal_lower,
    signal_upper_percent = 100 * signal_upper,
    cancellation_percent = 100 * cancellation_index,
    positive_predictive_skill = cv_r_squared > 0
  )

set.seed(bootstrap_seed + 1L)
macroregion_phase_summary <- oof |>
  filter(phase != "Neutral") |>
  group_by(macroregion, phase) |>
  group_modify(~ bind_cols(summarise_field(.x), bootstrap_field(.x))) |>
  ungroup() |>
  mutate(
    cv_r_squared_percent = 100 * cv_r_squared,
    cv_r2_lower_percent = 100 * cv_r2_lower,
    cv_r2_upper_percent = 100 * cv_r2_upper,
    signal_energy_percent = 100 * signal_energy_fraction,
    signal_lower_percent = 100 * signal_lower,
    signal_upper_percent = 100 * signal_upper,
    cancellation_percent = 100 * cancellation_index
  )

country_summary <- read_csv(paths$country, show_col_types = FALSE) |>
  select(
    country, observations, years, regions,
    spatial_field_cv_r_squared_percent,
    field_r2_ci_lower_percent, field_r2_ci_upper_percent,
    spatial_signal_energy_percent,
    signal_energy_ci_lower_percent, signal_energy_ci_upper_percent,
    national_net_cv_r_squared_percent,
    mean_absolute_component_percent,
    national_mean_absolute_component_percent,
    cancellation_index
  )

randomization <- read_csv(paths$randomization, show_col_types = FALSE) |>
  transmute(
    country,
    randomization_p_value,
    shifts_tested,
    null_95_percentile_cv_r2_percent =
      100 * null_95_percentile_spatial_field_cv_r_squared
  )

country_scale_comparison <- country_summary |>
  left_join(randomization, by = "country") |>
  transmute(
    country,
    spatial_field_cv_r2_percent = spatial_field_cv_r_squared_percent,
    spatial_field_lower_percent = field_r2_ci_lower_percent,
    spatial_field_upper_percent = field_r2_ci_upper_percent,
    national_cv_r2_percent = national_net_cv_r_squared_percent,
    scale_difference_percentage_points =
      national_net_cv_r_squared_percent - spatial_field_cv_r_squared_percent,
    spatial_signal_energy_percent,
    signal_energy_lower_percent = signal_energy_ci_lower_percent,
    signal_energy_upper_percent = signal_energy_ci_upper_percent,
    field_mean_absolute_component_percent = mean_absolute_component_percent,
    national_mean_absolute_component_percent,
    field_cancellation_percent = 100 * cancellation_index,
    randomization_p_value,
    shifts_tested,
    null_95_percentile_cv_r2_percent
  )

phase_country_summary <- read_csv(paths$phase, show_col_types = FALSE) |>
  mutate(
    phase = case_when(
      str_detect(str_to_lower(phase_simple), "nino") ~ "El Niño",
      str_detect(str_to_lower(phase_simple), "nina") ~ "La Niña",
      TRUE ~ "Neutral"
    )
  ) |>
  filter(phase != "Neutral") |>
  transmute(
    country, phase, observations, years, regions,
    spatial_field_cv_r_squared_percent,
    spatial_signal_energy_percent,
    mean_signed_component_percent,
    mean_absolute_component_percent,
    cancellation_percent = 100 * cancellation_index
  )

regional_cancellation <- read_csv(paths$cancellation, show_col_types = FALSE) |>
  rename(
    cancellation_mean_percent = mean,
    cancellation_lower_percent = lower,
    cancellation_upper_percent = upper
  )

qc_input_alignment <- tribble(
  ~dataset, ~rows, ~countries, ~spatial_units, ~years,
  "local_metrics", nrow(local_metrics), n_distinct(local_metrics$country),
  n_distinct(local_metrics$spatial_unit_id_join), NA_integer_,
  "local_consistency_relationship", nrow(local_consistency_relationship),
  n_distinct(local_consistency_relationship$country),
  n_distinct(local_consistency_relationship$spatial_unit_id_join), NA_integer_,
  "regional_oof", nrow(oof), n_distinct(oof$country),
  n_distinct(oof$region_key), n_distinct(oof$analysis_year)
)

qc_metric_ranges <- bind_rows(
  local_metrics |>
    summarise(
      level = "local",
      minimum_cv_r2_percent = min(cv_r_squared_percent, na.rm = TRUE),
      median_cv_r2_percent = median(cv_r_squared_percent, na.rm = TRUE),
      maximum_cv_r2_percent = max(cv_r_squared_percent, na.rm = TRUE),
      minimum_signal_energy_percent = min(cv_enso_signal_energy_percent, na.rm = TRUE),
      median_signal_energy_percent = median(cv_enso_signal_energy_percent, na.rm = TRUE),
      maximum_signal_energy_percent = max(cv_enso_signal_energy_percent, na.rm = TRUE)
    ),
  macroregion_field_summary |>
    summarise(
      level = "macroregion",
      minimum_cv_r2_percent = min(cv_r_squared_percent, na.rm = TRUE),
      median_cv_r2_percent = median(cv_r_squared_percent, na.rm = TRUE),
      maximum_cv_r2_percent = max(cv_r_squared_percent, na.rm = TRUE),
      minimum_signal_energy_percent = min(signal_energy_percent, na.rm = TRUE),
      median_signal_energy_percent = median(signal_energy_percent, na.rm = TRUE),
      maximum_signal_energy_percent = max(signal_energy_percent, na.rm = TRUE)
    )
)

methods <- tribble(
  ~item, ~value,
  "primary_estimand", "out of sample error reduction relative to the null model",
  "predictive_contribution", "100 times one minus weighted model SSE divided by weighted null SSE",
  "negative_skill", "retained numerically and displayed separately from positive skill",
  "map_display", "negative predictive skill is classified as no positive skill",
  "signal_energy", "squared ENSO component divided by null residual energy",
  "signal_energy_interpretation", "ENSO footprint magnitude, not causal explained variance",
  "macroregion_uncertainty", "circular moving block bootstrap over harvest years",
  "bootstrap_block_length_years", as.character(bootstrap_block_length),
  "bootstrap_repetitions", as.character(bootstrap_repetitions),
  "randomization_test", "circular shifts of the ENSO time series at country scale",
  "region_definition", paste(
    "Argentina; Southern Brazil; Central Brazil including MG, SP, PA, RO, AM and RR;",
    "MATOPIBA and northeastern frontier"
  )
)

outputs <- list(
  "Figure_2_local_ENSO_metrics.csv" = local_metrics,
  "Figure_2_local_skill_consistency_relationship.csv" = local_consistency_relationship,
  "Figure_2_local_macroregion_summary.csv" = local_macroregion_summary,
  "Figure_2_macroregion_field_summary.csv" = macroregion_field_summary,
  "Figure_2_macroregion_phase_summary.csv" = macroregion_phase_summary,
  "Figure_2_country_scale_comparison.csv" = country_scale_comparison,
  "Figure_2_country_phase_summary.csv" = phase_country_summary,
  "Figure_2_regional_cancellation.csv" = regional_cancellation,
  "QC_Figure_2_input_alignment.csv" = qc_input_alignment,
  "QC_Figure_2_metric_ranges.csv" = qc_metric_ranges,
  "README_Figure_2_methods.csv" = methods
)

iwalk(outputs, ~ write_csv(.x, file.path(output_dir, .y), na = ""))

if (nrow(local_consistency_relationship) < 0.95 * 2 * nrow(local_metrics)) {
  warning("Less than 95% of local metrics joined to both phase consistency records.")
}

message("Script 15a completed.")
message("Prepared ", nrow(local_metrics), " local ENSO records, ",
        nrow(macroregion_field_summary), " macroregion records and ",
        nrow(country_scale_comparison), " country scale comparisons.")
message("Outputs: ", normalizePath(output_dir, mustWork = TRUE))

#!/usr/bin/env Rscript

# Figure 1 data preparation
# Quantifying ENSO impacts on South American soybean productivity
#
# This script prepares, but does not plot, the five analytical components of
# Figure 1. The expected local ENSO sign is cross fitted: the observation from
# the event being evaluated is never used to define its expected sign.

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
  file.path(project_dir, "outputs", "figure_data", "figure_1")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

minimum_total_training_years <- 15L
minimum_phase_training_events <- 4L
bootstrap_replicates <- 2000L
bootstrap_seed <- 1401L

find_one <- function(root, candidates) {
  all_files <- list.files(root, recursive = TRUE, full.names = TRUE)
  hits <- all_files[basename(all_files) %in% candidates]
  if (length(hits) == 0L) {
    stop("Could not find any of: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  hits <- hits[order(nchar(hits), hits)]
  normalizePath(hits[[1]], mustWork = TRUE)
}

yield_file <- find_one(project_dir, c(
  "soybean_Brazil_Argentina_analysis_ready_yield.csv",
  "soybean_Brazil_Argentina_analysis_ready_yield(1).csv",
  "soybean_Brazil_Argentina_analysis_ready_yield(2).csv"
))

enso_file <- find_one(project_dir, c(
  "ENSO_indices_by_soybean_harvest_year_1974_2024.csv"
))

area_file <- find_one(project_dir, c(
  "soybean_spatial_units_with_calendars.csv"
))

message("Reading annual soybean observations...")

yield_raw <- read_csv(yield_file, show_col_types = FALSE, progress = FALSE) |>
  mutate(
    country = as.character(country),
    spatial_unit_id = as.character(spatial_unit_id),
    analysis_year = as.integer(analysis_year),
    yield_t_ha = as.numeric(yield_harvested_t_ha),
    analysis_ready_yield = as.logical(analysis_ready_yield),
    unit_key = paste(country, spatial_unit_id, sep = "::")
  ) |>
  filter(
    analysis_ready_yield %in% TRUE,
    is.finite(yield_t_ha),
    yield_t_ha >= 0,
    is.finite(analysis_year)
  ) |>
  distinct(unit_key, analysis_year, .keep_all = TRUE)

enso <- read_csv(enso_file, show_col_types = FALSE) |>
  transmute(
    analysis_year = as.integer(harvest_year),
    oni_djf_c = as.numeric(oni_djf_c),
    phase = case_when(
      str_detect(str_to_lower(oni_harvest_year_phase), "nino") ~ "El Niño",
      str_detect(str_to_lower(oni_harvest_year_phase), "nina") ~ "La Niña",
      TRUE ~ "Neutral"
    ),
    intensity = str_to_sentence(as.character(oni_harvest_year_strength)),
    event_id = paste0("ENSO_", analysis_year)
  ) |>
  mutate(
    intensity = case_when(
      phase == "Neutral" ~ "Neutral",
      is.na(intensity) | intensity == "" ~ "Unclassified",
      TRUE ~ intensity
    )
  )

area_reference <- read_csv(area_file, show_col_types = FALSE, progress = FALSE) |>
  transmute(
    unit_key = paste(as.character(country), as.character(spatial_unit_id), sep = "::"),
    reference_area_ha = as.numeric(reference_harvested_area_ha),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude),
    calendar_group = as.character(calendar_group)
  ) |>
  distinct(unit_key, .keep_all = TRUE)

analysis <- yield_raw |>
  inner_join(enso, by = "analysis_year") |>
  left_join(area_reference, by = "unit_key") |>
  mutate(
    parent_region = coalesce(as.character(parent_region_name),
                             as.character(parent_region_code)),
    parent_region = str_to_upper(str_squish(parent_region)),
    macroregion = case_when(
      country == "Argentina" ~ "Argentina",
      country == "Brazil" & parent_region %in% c("PR", "SC", "RS") ~ "Southern Brazil",
      country == "Brazil" & parent_region %in% c(
        "MT", "MS", "GO", "DF", "MG", "SP", "PA", "RO", "AM", "RR"
      ) ~ "Central Brazil",
      country == "Brazil" & parent_region %in% c(
        "TO", "MA", "PI", "BA", "AL", "CE", "PB", "PE", "RN", "SE"
      ) ~ "MATOPIBA and northeastern frontier",
      country == "Brazil" ~ "Central Brazil",
      TRUE ~ "Outside focal regions"
    ),
    focal_region = macroregion %in% c(
      "Argentina", "Southern Brazil", "Central Brazil",
      "MATOPIBA and northeastern frontier"
    ),
    area_ha = case_when(
      is.finite(reference_area_ha) & reference_area_ha > 0 ~ reference_area_ha,
      is.finite(harvested_area_ha) & harvested_area_ha > 0 ~ as.numeric(harvested_area_ha),
      TRUE ~ NA_real_
    )
  )

if (nrow(analysis) == 0L) stop("No aligned soybean and ENSO observations.", call. = FALSE)

# Cross fitted linear detrending using OLS sufficient statistics. All rows from
# the held out event are removed before the trend for that event is estimated.
message("Computing event held out yield anomalies...")

unit_totals <- analysis |>
  group_by(unit_key) |>
  summarise(
    n_all = n(), sx_all = sum(analysis_year), sy_all = sum(yield_t_ha),
    sxx_all = sum(analysis_year^2), sxy_all = sum(analysis_year * yield_t_ha),
    .groups = "drop"
  )

event_totals <- analysis |>
  group_by(unit_key, event_id) |>
  summarise(
    n_event = n(), sx_event = sum(analysis_year), sy_event = sum(yield_t_ha),
    sxx_event = sum(analysis_year^2),
    sxy_event = sum(analysis_year * yield_t_ha), .groups = "drop"
  )

crossfit <- analysis |>
  left_join(unit_totals, by = "unit_key") |>
  left_join(event_totals, by = c("unit_key", "event_id")) |>
  mutate(
    n_train = n_all - n_event,
    sx = sx_all - sx_event, sy = sy_all - sy_event,
    sxx = sxx_all - sxx_event, sxy = sxy_all - sxy_event,
    trend_denominator = n_train * sxx - sx^2,
    trend_slope = if_else(
      n_train >= minimum_total_training_years & abs(trend_denominator) > 1e-10,
      (n_train * sxy - sx * sy) / trend_denominator, NA_real_
    ),
    trend_intercept = if_else(
      is.finite(trend_slope), (sy - trend_slope * sx) / n_train, NA_real_
    ),
    expected_yield_t_ha = trend_intercept + trend_slope * analysis_year,
    yield_anomaly_t_ha = yield_t_ha - expected_yield_t_ha,
    yield_anomaly_percent = 100 * yield_anomaly_t_ha / expected_yield_t_ha
  ) |>
  filter(is.finite(yield_anomaly_t_ha), is.finite(yield_anomaly_percent))

# Define the expected phase sign without the held out event. Cross fitted
# anomalies are used, and the evaluated row is subtracted from the phase sum.
phase_totals <- crossfit |>
  filter(phase != "Neutral") |>
  group_by(unit_key, phase) |>
  summarise(
    phase_sum = sum(yield_anomaly_t_ha),
    phase_n = n(),
    phase_event_n = n_distinct(event_id),
    .groups = "drop"
  )

consistency <- crossfit |>
  filter(phase != "Neutral") |>
  left_join(phase_totals, by = c("unit_key", "phase")) |>
  mutate(
    training_phase_n = phase_n - 1L,
    training_phase_event_n = phase_event_n - 1L,
    expected_phase_anomaly_t_ha = if_else(
      training_phase_event_n >= minimum_phase_training_events,
      (phase_sum - yield_anomaly_t_ha) / training_phase_n,
      NA_real_
    ),
    expected_sign = sign(expected_phase_anomaly_t_ha),
    observed_sign = sign(yield_anomaly_t_ha),
    sign_consistent = if_else(
      is.finite(expected_phase_anomaly_t_ha) & expected_sign != 0 & observed_sign != 0,
      expected_sign == observed_sign, NA
    )
  )

local_phase_response <- crossfit |>
  filter(phase != "Neutral") |>
  group_by(
    country, spatial_unit_id, spatial_unit_name, parent_region_name,
    spatial_level, unit_key, macroregion, focal_region, phase
  ) |>
  summarise(
    observations = n(), events = n_distinct(event_id),
    mean_yield_anomaly_t_ha = mean(yield_anomaly_t_ha),
    median_yield_anomaly_t_ha = median(yield_anomaly_t_ha),
    mean_yield_anomaly_percent = mean(yield_anomaly_percent),
    reference_area_ha = first(area_ha),
    longitude = first(longitude), latitude = first(latitude),
    .groups = "drop"
  )

crossvalidated_consistency <- consistency |>
  group_by(
    country, spatial_unit_id, spatial_unit_name, parent_region_name,
    spatial_level, unit_key, macroregion, focal_region, phase
  ) |>
  summarise(
    evaluated_events = sum(!is.na(sign_consistent)),
    consistent_events = sum(sign_consistent %in% TRUE),
    sign_consistency = mean(sign_consistent, na.rm = TRUE),
    expected_positive_fraction = mean(expected_sign > 0, na.rm = TRUE),
    expected_negative_fraction = mean(expected_sign < 0, na.rm = TRUE),
    dominant_expected_sign = case_when(
      expected_positive_fraction > expected_negative_fraction ~ "Gain",
      expected_negative_fraction > expected_positive_fraction ~ "Loss",
      TRUE ~ "Indeterminate"
    ),
    reference_area_ha = first(area_ha),
    longitude = first(longitude), latitude = first(latitude),
    .groups = "drop"
  ) |>
  mutate(sign_consistency = if_else(is.nan(sign_consistency), NA_real_, sign_consistency))

event_spatial_extent <- consistency |>
  filter(focal_region, is.finite(area_ha), area_ha > 0) |>
  group_by(macroregion, phase, intensity, event_id, analysis_year) |>
  summarise(
    n_spatial_units = n_distinct(unit_key),
    observed_area_ha = sum(area_ha),
    gain_area_ha = sum(area_ha[observed_sign > 0], na.rm = TRUE),
    loss_area_ha = sum(area_ha[observed_sign < 0], na.rm = TRUE),
    expected_gain_area_ha = sum(area_ha[expected_sign > 0], na.rm = TRUE),
    expected_loss_area_ha = sum(area_ha[expected_sign < 0], na.rm = TRUE),
    consistent_area_ha = sum(area_ha[sign_consistent %in% TRUE], na.rm = TRUE),
    weighted_anomaly_t_ha = weighted.mean(yield_anomaly_t_ha, area_ha),
    weighted_absolute_anomaly_t_ha = weighted.mean(abs(yield_anomaly_t_ha), area_ha),
    .groups = "drop"
  ) |>
  mutate(
    gain_area_percent = 100 * gain_area_ha / observed_area_ha,
    loss_area_percent = 100 * loss_area_ha / observed_area_ha,
    consistent_area_percent = 100 * consistent_area_ha / observed_area_ha,
    cancellation_percent = 100 * (1 - abs(weighted_anomaly_t_ha) /
                                      weighted_absolute_anomaly_t_ha)
  )

bootstrap_event_mean <- function(x, B = bootstrap_replicates) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) return(tibble(mean = NA_real_, lower = NA_real_, upper = NA_real_))
  if (n < 3L) return(tibble(mean = mean(x), lower = NA_real_, upper = NA_real_))
  draws <- replicate(B, mean(sample(x, n, replace = TRUE)))
  tibble(mean = mean(x), lower = unname(quantile(draws, 0.025)),
         upper = unname(quantile(draws, 0.975)))
}

set.seed(bootstrap_seed)

intensity_summary <- event_spatial_extent |>
  select(macroregion, phase, intensity, event_id,
         gain_area_percent, loss_area_percent, consistent_area_percent) |>
  pivot_longer(ends_with("_percent"), names_to = "metric", values_to = "value") |>
  group_by(macroregion, phase, intensity, metric) |>
  group_modify(~ bootstrap_event_mean(.x$value) |>
                 mutate(n_events = n_distinct(.x$event_id))) |>
  ungroup() |>
  mutate(
    intensity_order = match(str_to_lower(intensity),
                            c("very weak", "weak", "moderate", "strong", "very strong")),
    category_label = paste0(intensity, " ", phase, " (n = ", n_events, ")")
  ) |>
  arrange(macroregion, phase, intensity_order, metric)

set.seed(bootstrap_seed + 1L)

regional_cancellation <- event_spatial_extent |>
  group_by(macroregion, phase) |>
  group_modify(~ bootstrap_event_mean(.x$cancellation_percent) |>
                 mutate(n_events = n_distinct(.x$event_id))) |>
  ungroup()

region_coverage_qc <- analysis |>
  filter(is.finite(area_ha), area_ha > 0) |>
  distinct(unit_key, country, macroregion, focal_region, area_ha) |>
  group_by(country, macroregion, focal_region) |>
  summarise(n_spatial_units = n(), productive_area_ha = sum(area_ha), .groups = "drop") |>
  group_by(country) |>
  mutate(country_area_percent = 100 * productive_area_ha / sum(productive_area_ha)) |>
  ungroup()

sample_size_qc <- crossfit |>
  count(country, phase, intensity, name = "observations") |>
  left_join(
    crossfit |> distinct(country, phase, intensity, event_id) |>
      count(country, phase, intensity, name = "events"),
    by = c("country", "phase", "intensity")
  )

method_qc <- tribble(
  ~item, ~value,
  "yield_metric", "harvested area yield in t ha-1",
  "detrending", "unit specific linear trend estimated after excluding the evaluated ENSO event",
  "expected_sign", "mean phase anomaly from other events for the same spatial unit",
  "minimum_total_training_years", as.character(minimum_total_training_years),
  "minimum_phase_training_events", as.character(minimum_phase_training_events),
  "area_weight", "reference harvested area, with annual harvested area fallback",
  "uncertainty", "nonparametric bootstrap over complete ENSO events",
  "bootstrap_replicates", as.character(bootstrap_replicates),
  "bootstrap_seed", as.character(bootstrap_seed),
  "focal_regions", paste(c("Argentina", "Southern Brazil", "Central Brazil",
                            "MATOPIBA and northeastern frontier"), collapse = "; ")
)

outputs <- list(
  "Figure_1_local_phase_response.csv" = local_phase_response,
  "Figure_1_crossvalidated_consistency.csv" = crossvalidated_consistency,
  "Figure_1_event_spatial_extent.csv" = event_spatial_extent,
  "Figure_1_intensity_summary.csv" = intensity_summary,
  "Figure_1_regional_cancellation.csv" = regional_cancellation,
  "QC_Figure_1_region_coverage.csv" = region_coverage_qc,
  "QC_Figure_1_sample_sizes.csv" = sample_size_qc,
  "README_Figure_1_methods.csv" = method_qc
)

iwalk(outputs, ~ write_csv(.x, file.path(output_dir, .y), na = ""))

message("Script 14a completed.")
message("Prepared ", nrow(local_phase_response), " local phase response records and ",
        nrow(event_spatial_extent), " region event records.")
message("Outputs: ", normalizePath(output_dir, mustWork = TRUE))

# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 12_spatial_climate_mediation_ENSO.R
# PURPOSE: Quantify spatially resolved climate prediction, identify the most
#          informative indices, and estimate the ENSO-aligned climate signal
# VERSION: 3.1
# ==============================================================================


# 1. USER SETTINGS -------------------------------------------------------------

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")

project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) project_dir_override else default_project_dir

analysis_start_year <- 1985L
analysis_end_year <- 2024L
outer_folds <- 5L
inner_folds <- 3L
ridge_candidates <- c(0.1, 1, 10, 100, 1000, 10000)
common_slope_penalty_multiplier <- 0.25
permutation_repetitions <- 20L
minimum_local_validation_rows <- 10L
weight_cap_probability <- 0.99
random_seed <- 20260907L
reuse_existing_checkpoints <- TRUE

climate_metrics <- c(
  "mean_temperature_c",
  "precipitation_total_mm",
  "mean_vpd_kpa",
  "hot_dry_fraction",
  "hot_wet_fraction",
  "heavy_rain_fraction"
)

water_balance_metrics <- c(
  "reproductive_requirement_satisfaction"
)

# One harmonized seven-index set is used in every model, contribution estimate,
# map and ordination. No post-hoc mechanism groups are constructed.
index_predictors <- c(
  "climate_z_mean_temperature_c",
  "climate_z_mean_vpd_kpa",
  "climate_z_precipitation_total_mm",
  "climate_z_heavy_rain_fraction",
  "water_z_crop_water_stress",
  "climate_z_hot_dry_fraction",
  "climate_z_hot_wet_fraction"
)


# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c("dplyr", "readr", "tibble", "tidyr", "purrr", "stringr")
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
  library(purrr)
  library(stringr)
})


# 3. PATHS ---------------------------------------------------------------------

climate_dir <- file.path(project_dir, "climate_mechanisms")
merged_dir <- file.path(project_dir, "data_merged")
output_dir <- file.path(climate_dir, "spatial_attribution")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

grid_manifest_file <- file.path(climate_dir, "NASA_POWER_grid_manifest.csv")
grid_climate_file <- file.path(climate_dir, "soybean_climate_indices_grid_year.csv.gz")
soybean_panel_file <- file.path(merged_dir, "soybean_ENSO_spatial_panel_1974_2024.csv")
water_balance_file <- file.path(
  climate_dir, "soybean_water_balance_primary_AWC100.csv"
)

required_files <- c(
  grid_manifest_file, grid_climate_file, soybean_panel_file, water_balance_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Required inputs are missing: ", paste(missing_files, collapse = "; "),
    call. = FALSE
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

standardize_unit_id <- function(country, x) {
  ifelse(
    country == "Brazil",
    standardize_numeric_code(x, 7L),
    standardize_numeric_code(x, 5L)
  )
}

safe_weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  weighted.mean(x[keep], w[keep])
}

safe_weighted_mse <- function(observed, predicted, weights) {
  keep <- is.finite(observed) & is.finite(predicted) &
    is.finite(weights) & weights > 0
  if (sum(keep) < 3) return(NA_real_)
  weighted.mean((observed[keep] - predicted[keep])^2, weights[keep])
}

safe_r_squared <- function(observed, predicted, weights) {
  keep <- is.finite(observed) & is.finite(predicted) &
    is.finite(weights) & weights > 0
  if (sum(keep) < 3) return(NA_real_)
  observed <- observed[keep]
  predicted <- predicted[keep]
  weights <- weights[keep]
  centre <- weighted.mean(observed, weights)
  sse <- sum(weights * (observed - predicted)^2)
  sst <- sum(weights * (observed - centre)^2)
  if (!is.finite(sst) || sst <= 0) return(NA_real_)
  1 - sse / sst
}

group_weighted_mean <- function(x, weights, groups) {
  indices <- split(seq_along(x), as.character(groups), drop = TRUE)
  vapply(
    indices,
    function(index) safe_weighted_mean(x[index], weights[index]),
    numeric(1)
  )
}

make_phase_stratified_folds <- function(years, phases, number_folds) {
  tibble(analysis_year = as.integer(years), phase = as.character(phases)) |>
    distinct(analysis_year, .keep_all = TRUE) |>
    arrange(analysis_year) |>
    group_by(phase) |>
    mutate(fold = ((row_number() - 1L) %% number_folds) + 1L) |>
    ungroup() |>
    select(analysis_year, fold)
}

build_hierarchical_matrix <- function(predictor_matrix, regions, region_levels) {
  predictor_matrix <- as.matrix(predictor_matrix)
  storage.mode(predictor_matrix) <- "double"
  region_matrix <- vapply(
    region_levels,
    function(region_value) as.numeric(as.character(regions) == region_value),
    numeric(length(regions))
  )
  if (length(region_levels) == 1L) {
    region_matrix <- matrix(region_matrix, ncol = 1L)
  }
  colnames(region_matrix) <- paste0("region_", make.names(region_levels))

  deviations <- do.call(
    cbind,
    lapply(seq_len(ncol(predictor_matrix)), function(predictor_index) {
      result <- region_matrix * predictor_matrix[, predictor_index]
      colnames(result) <- paste0(
        colnames(predictor_matrix)[predictor_index], "__", colnames(region_matrix)
      )
      result
    })
  )
  cbind(predictor_matrix, deviations)
}

solve_weighted_ridge <- function(x, y, weights, penalty) {
  keep <- is.finite(y) & is.finite(weights) & weights > 0 &
    apply(x, 1L, function(value) all(is.finite(value)))
  if (sum(keep) < 3) stop("Insufficient finite rows for ridge estimation")
  x <- x[keep, , drop = FALSE]
  y <- y[keep]
  weights <- weights[keep]
  weights <- weights / mean(weights)
  x_weighted <- x * sqrt(weights)
  y_weighted <- y * sqrt(weights)
  system_matrix <- crossprod(x_weighted) + diag(penalty, ncol(x_weighted))
  rhs <- crossprod(x_weighted, y_weighted)
  decomposition <- eigen(system_matrix, symmetric = TRUE)
  tolerance <- max(decomposition$values) * .Machine$double.eps *
    max(dim(system_matrix))
  inverse_values <- ifelse(
    decomposition$values > tolerance,
    1 / decomposition$values,
    0
  )
  as.numeric(
    decomposition$vectors %*%
      (inverse_values * crossprod(decomposition$vectors, rhs))
  )
}

fit_hierarchical_ridge <- function(
  training, predictors, response_column, lambda, weight_column = "model_weight"
) {
  if (length(predictors) == 0) stop("At least one predictor is required")
  response <- training[[response_column]]
  weights <- training[[weight_column]]
  units <- as.character(training$spatial_unit_id)
  regions <- as.character(training$parent_region_name)

  response_means <- group_weighted_mean(response, weights, units)
  global_response_mean <- safe_weighted_mean(response, weights)
  response_baseline <- unname(response_means[units])
  response_baseline[!is.finite(response_baseline)] <- global_response_mean

  predictor_means <- lapply(predictors, function(predictor) {
    group_weighted_mean(training[[predictor]], weights, units)
  })
  names(predictor_means) <- predictors
  global_predictor_means <- vapply(
    predictors,
    function(predictor) safe_weighted_mean(training[[predictor]], weights),
    numeric(1)
  )

  centred_predictors <- sapply(seq_along(predictors), function(index) {
    predictor <- predictors[index]
    baseline <- unname(predictor_means[[index]][units])
    baseline[!is.finite(baseline)] <- global_predictor_means[index]
    training[[predictor]] - baseline
  })
  centred_predictors <- matrix(
    centred_predictors, nrow = nrow(training), ncol = length(predictors),
    dimnames = list(NULL, predictors)
  )

  region_levels <- sort(unique(regions[!is.na(regions)]))
  design <- build_hierarchical_matrix(
    centred_predictors, regions, region_levels
  )
  common_penalty <- rep(
    lambda * common_slope_penalty_multiplier,
    length(predictors)
  )
  deviation_penalty <- rep(lambda, length(predictors) * length(region_levels))
  penalty <- c(common_penalty, deviation_penalty)
  coefficients <- solve_weighted_ridge(
    design,
    response - response_baseline,
    weights,
    penalty
  )

  list(
    predictors = predictors,
    response_column = response_column,
    response_means = response_means,
    global_response_mean = global_response_mean,
    predictor_means = predictor_means,
    global_predictor_means = global_predictor_means,
    region_levels = region_levels,
    coefficients = coefficients,
    lambda = lambda
  )
}

predict_hierarchical_ridge <- function(
  fit, new_data, permute_predictors = NULL, permutation_seed = NULL,
  return_signal = FALSE
) {
  units <- as.character(new_data$spatial_unit_id)
  regions <- as.character(new_data$parent_region_name)
  response_baseline <- unname(fit$response_means[units])
  response_baseline[!is.finite(response_baseline)] <- fit$global_response_mean

  centred_predictors <- sapply(seq_along(fit$predictors), function(index) {
    predictor <- fit$predictors[index]
    baseline <- unname(fit$predictor_means[[index]][units])
    baseline[!is.finite(baseline)] <- fit$global_predictor_means[index]
    new_data[[predictor]] - baseline
  })
  centred_predictors <- matrix(
    centred_predictors, nrow = nrow(new_data), ncol = length(fit$predictors),
    dimnames = list(NULL, fit$predictors)
  )

  if (!is.null(permute_predictors)) {
    permute_predictors <- as.character(permute_predictors)
    unknown_predictors <- setdiff(permute_predictors, fit$predictors)
    if (length(unknown_predictors) > 0) {
      stop(
        "Unknown permutation predictors: ",
        paste(unknown_predictors, collapse = ", ")
      )
    }
    if (!is.null(permutation_seed)) set.seed(permutation_seed)
    column_indices <- match(permute_predictors, fit$predictors)
    region_indices <- split(seq_len(nrow(new_data)), regions, drop = TRUE)
    for (index in region_indices) {
      permutation_index <- sample(index, length(index), replace = FALSE)
      centred_predictors[index, column_indices] <- centred_predictors[
        permutation_index, column_indices, drop = FALSE
      ]
    }
  }

  design <- build_hierarchical_matrix(
    centred_predictors, regions, fit$region_levels
  )
  signal <- as.numeric(design %*% fit$coefficients)
  prediction <- response_baseline + signal
  prediction[!is.finite(prediction)] <- NA_real_
  if (return_signal) signal else prediction
}

predict_unit_baseline <- function(training, testing, response_column) {
  means <- group_weighted_mean(
    training[[response_column]], training$model_weight,
    training$spatial_unit_id
  )
  global_mean <- safe_weighted_mean(
    training[[response_column]], training$model_weight
  )
  prediction <- unname(means[as.character(testing$spatial_unit_id)])
  prediction[!is.finite(prediction)] <- global_mean
  prediction
}

tune_lambda <- function(training, predictors, response_column) {
  # Hyperparameters are tuned on the state or province by year field. This
  # preserves temporal independence while avoiding repeated multiplication of
  # thousands of local rows that share the same regional shrinkage parameter.
  tuning_source <- training |>
    group_by(country, parent_region_name, analysis_year, enso_phase_simple) |>
    summarise(
      across(
        all_of(unique(c(response_column, predictors))),
        ~ safe_weighted_mean(.x, model_weight)
      ),
      model_weight = sum(model_weight, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(spatial_unit_id = as.character(parent_region_name))

  folds <- make_phase_stratified_folds(
    tuning_source$analysis_year,
    tuning_source$enso_phase_simple,
    inner_folds
  ) |>
    rename(inner_fold = fold)
  tuning_data <- tuning_source |> left_join(folds, by = "analysis_year")
  scores <- map_dfr(ridge_candidates, function(lambda) {
    fold_scores <- map_dbl(seq_len(inner_folds), function(fold_id) {
      inner_training <- filter(tuning_data, inner_fold != fold_id)
      inner_testing <- filter(tuning_data, inner_fold == fold_id)
      fit <- tryCatch(
        fit_hierarchical_ridge(
          inner_training, predictors, response_column, lambda
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NA_real_)
      prediction <- predict_hierarchical_ridge(fit, inner_testing)
      safe_weighted_mse(
        inner_testing[[response_column]], prediction, inner_testing$model_weight
      )
    })
    tibble(
      lambda = lambda,
      validation_mse = mean(fold_scores, na.rm = TRUE),
      valid_folds = sum(is.finite(fold_scores))
    )
  })
  valid <- scores |> filter(valid_folds == inner_folds, is.finite(validation_mse))
  if (nrow(valid) == 0) stop("No ridge candidate completed inner validation")
  valid |> arrange(validation_mse, desc(lambda)) |> slice(1L) |> pull(lambda)
}

extract_effective_coefficients <- function(fit, country, model, outer_fold) {
  predictor_count <- length(fit$predictors)
  region_count <- length(fit$region_levels)
  common <- fit$coefficients[seq_len(predictor_count)]
  deviations <- matrix(
    fit$coefficients[-seq_len(predictor_count)],
    nrow = region_count,
    ncol = predictor_count,
    byrow = FALSE
  )
  colnames(deviations) <- fit$predictors
  tibble(
    country = country,
    model = model,
    outer_fold = outer_fold,
    parent_region_name = rep(fit$region_levels, each = predictor_count),
    predictor = rep(fit$predictors, times = region_count),
    common_coefficient = rep(common, times = region_count),
    regional_deviation = as.numeric(t(deviations)),
    effective_coefficient = common_coefficient + regional_deviation,
    selected_lambda = fit$lambda
  )
}

fit_weighted_formula_mse <- function(formula, data) {
  frame <- model.frame(formula, data = data, na.action = na.pass)
  y <- model.response(frame)
  x <- model.matrix(terms(frame), frame)
  w <- data$model_weight
  keep <- is.finite(y) & is.finite(w) & w > 0 &
    apply(x, 1L, function(value) all(is.finite(value)))
  x <- x[keep, , drop = FALSE]
  y <- y[keep]
  w <- w[keep] / mean(w[keep])
  xw <- x * sqrt(w)
  yw <- y * sqrt(w)
  decomposition <- svd(xw)
  tolerance <- max(dim(xw)) * max(decomposition$d) * .Machine$double.eps
  retained <- decomposition$d > tolerance
  coefficients <- decomposition$v[, retained, drop = FALSE] %*%
    (crossprod(decomposition$u[, retained, drop = FALSE], yw) /
       decomposition$d[retained])
  fitted <- as.numeric(x %*% coefficients)
  weighted.mean((y - fitted)^2, w)
}


# 5. READ AND ASSEMBLE LOCAL PANEL --------------------------------------------

message("Reading spatial soybean and climate data...")

grid_manifest <- read_csv(grid_manifest_file, show_col_types = FALSE)
grid_climate <- read_csv(grid_climate_file, show_col_types = FALSE)
soybean <- read_csv(soybean_panel_file, show_col_types = FALSE)
water_balance <- read_csv(water_balance_file, show_col_types = FALSE)

check_required_columns(
  grid_manifest,
  c("country", "spatial_unit_id", "parent_region_name", "grid_id", "calendar_group"),
  "NASA POWER grid manifest"
)
check_required_columns(
  grid_climate,
  c(
    "grid_id", "country", "calendar_group", "climate_window", "harvest_year",
    "oni_djf_c", "enso_phase_simple", "climate_metric",
    "climate_anomaly_standardized"
  ),
  "Grid climate indices"
)
check_required_columns(
  soybean,
  c(
    "country", "spatial_unit_id", "parent_region_name", "analysis_year",
    "yield_anomaly_percent", "harvested_area_ha", "oni_harvest_year_phase"
  ),
  "Soybean ENSO panel"
)
check_required_columns(
  water_balance,
  c("country", "spatial_unit_id", "harvest_year", "awc_mm", water_balance_metrics),
  "AWC100 soybean water balance"
)

grid_manifest <- grid_manifest |>
  transmute(
    country = as.character(country),
    spatial_unit_id = standardize_unit_id(country, spatial_unit_id),
    parent_region_name = as.character(parent_region_name),
    grid_id = as.character(grid_id),
    calendar_group = as.character(calendar_group)
  ) |>
  distinct(country, spatial_unit_id, .keep_all = TRUE)

climate_wide <- grid_climate |>
  filter(
    climate_window == "reproductive",
    climate_metric %in% climate_metrics,
    harvest_year >= analysis_start_year,
    harvest_year <= analysis_end_year
  ) |>
  select(
    grid_id, country, calendar_group, harvest_year, oni_djf_c,
    enso_phase_simple, climate_metric, climate_anomaly_standardized
  ) |>
  pivot_wider(
    names_from = climate_metric,
    values_from = climate_anomaly_standardized,
    names_prefix = "climate_z_"
  )

climate_predictors <- paste0("climate_z_", climate_metrics)
water_predictors <- "water_z_crop_water_stress"

water_prepared <- water_balance |>
  filter(
    awc_mm == 100,
    harvest_year >= analysis_start_year,
    harvest_year <= analysis_end_year
  ) |>
  transmute(
    country = as.character(country),
    spatial_unit_id = standardize_unit_id(country, spatial_unit_id),
    analysis_year = as.integer(harvest_year),
    crop_water_stress_raw = -as.numeric(reproductive_requirement_satisfaction)
  ) |>
  group_by(country) |>
  mutate(
    water_z_crop_water_stress = as.numeric(scale(crop_water_stress_raw))
  ) |>
  ungroup() |>
  select(country, spatial_unit_id, analysis_year, all_of(water_predictors))

if (any(duplicated(water_prepared[, c(
  "country", "spatial_unit_id", "analysis_year"
)]))) {
  stop("Duplicated AWC100 water-balance unit-year rows were detected", call. = FALSE)
}

climate_availability <- climate_wide |>
  mutate(
    raw_climate_available = TRUE,
    complete_raw_climate = if_all(
      all_of(climate_predictors),
      ~ is.finite(.x)
    )
  ) |>
  group_by(country, grid_id, calendar_group) |>
  summarise(
    raw_climate_available = any(raw_climate_available),
    complete_raw_climate = all(complete_raw_climate),
    .groups = "drop"
  )

missing_predictor_qc <- climate_wide |>
  pivot_longer(
    cols = all_of(climate_predictors),
    names_to = "predictor",
    values_to = "value"
  ) |>
  group_by(country, predictor) |>
  summarise(
    rows = n(),
    nonfinite_rows = sum(!is.finite(value)),
    affected_grid_calendar_combinations = n_distinct(
      paste(grid_id[!is.finite(value)], calendar_group[!is.finite(value)])
    ),
    .groups = "drop"
  )

write_csv(
  missing_predictor_qc,
  file.path(output_dir, "QC_spatial_climate_missing_predictors.csv"),
  na = ""
)

input_alignment_qc <- water_prepared |>
  distinct(country, spatial_unit_id) |>
  left_join(
    grid_manifest |>
      transmute(country, spatial_unit_id, grid_id, calendar_group),
    by = c("country", "spatial_unit_id")
  ) |>
  left_join(
    climate_availability,
    by = c("country", "grid_id", "calendar_group")
  ) |>
  mutate(
    manifest_available = !is.na(grid_id),
    raw_climate_available = coalesce(raw_climate_available, FALSE),
    complete_raw_climate = coalesce(complete_raw_climate, FALSE)
  ) |>
  group_by(country) |>
  summarise(
    expected_spatial_units = n(),
    manifest_spatial_units = sum(manifest_available),
    raw_climate_spatial_units = sum(raw_climate_available),
    raw_climate_unit_coverage =
      raw_climate_spatial_units / expected_spatial_units,
    complete_raw_climate_spatial_units = sum(complete_raw_climate),
    complete_raw_climate_unit_coverage =
      complete_raw_climate_spatial_units / expected_spatial_units,
    .groups = "drop"
  )

write_csv(
  input_alignment_qc,
  file.path(output_dir, "QC_spatial_climate_input_alignment.csv"),
  na = ""
)

if (any(input_alignment_qc$complete_raw_climate_unit_coverage < 0.97)) {
  failed <- input_alignment_qc |>
    filter(complete_raw_climate_unit_coverage < 0.97) |>
    transmute(
      label = paste0(
        country, ": ", complete_raw_climate_spatial_units, "/",
        expected_spatial_units, " units (",
        round(100 * complete_raw_climate_unit_coverage, 1), "%)"
      )
    ) |>
    pull(label)
  stop(
    "Finite climate-predictor coverage is below 97%: ",
    paste(failed, collapse = "; "),
    ". Inspect QC_spatial_climate_missing_predictors.csv and rerun script 11 ",
    "with the current spatial manifest before script 12.",
    call. = FALSE
  )
}

soybean_prepared <- soybean |>
  transmute(
    country = as.character(country),
    spatial_unit_id = standardize_unit_id(country, spatial_unit_id),
    soybean_region = as.character(parent_region_name),
    analysis_year = as.integer(analysis_year),
    yield_anomaly_percent = as.numeric(yield_anomaly_percent),
    harvested_area_ha = as.numeric(harvested_area_ha),
    oni_harvest_year_phase = as.character(oni_harvest_year_phase)
  )

all_soybean_units <- soybean_prepared |>
  distinct(country, spatial_unit_id, .keep_all = TRUE) |>
  select(country, spatial_unit_id, soybean_region)

analysis_soybean <- soybean_prepared |>
  filter(
    analysis_year >= analysis_start_year,
    analysis_year <= analysis_end_year,
    is.finite(yield_anomaly_percent),
    is.finite(harvested_area_ha), harvested_area_ha > 0
  )

local_panel <- analysis_soybean |>
  inner_join(grid_manifest, by = c("country", "spatial_unit_id")) |>
  inner_join(
    climate_wide,
    by = c(
      "country", "grid_id", "calendar_group",
      "analysis_year" = "harvest_year"
    )
  ) |>
  inner_join(
    water_prepared,
    by = c("country", "spatial_unit_id", "analysis_year")
  ) |>
  mutate(
    parent_region_name = coalesce(parent_region_name, soybean_region),
    enso_phase_simple = case_when(
      as.character(enso_phase_simple) == "El Nino" ~ "El Nino",
      as.character(enso_phase_simple) == "La Nina" ~ "La Nina",
      TRUE ~ "Neutral"
    ),
    enso_warm = pmax(oni_djf_c, 0),
    enso_cold = pmin(oni_djf_c, 0)
  ) |>
  select(-soybean_region)

all_mechanism_predictors <- index_predictors
complete_predictors <- complete.cases(local_panel[, all_mechanism_predictors])
local_panel <- local_panel[complete_predictors, , drop = FALSE]

unit_inclusion_audit <- all_soybean_units |>
  left_join(
    analysis_soybean |>
      group_by(country, spatial_unit_id) |>
      summarise(
        valid_productivity_rows_1985_2024 = n(),
        .groups = "drop"
      ),
    by = c("country", "spatial_unit_id")
  ) |>
  left_join(
    grid_manifest |>
      transmute(
        country, spatial_unit_id,
        manifest_available = TRUE
      ),
    by = c("country", "spatial_unit_id")
  ) |>
  left_join(
    local_panel |>
      group_by(country, spatial_unit_id) |>
      summarise(
        analytical_rows = n(),
        climate_years = n_distinct(analysis_year),
        .groups = "drop"
      ),
    by = c("country", "spatial_unit_id")
  ) |>
  mutate(
    valid_productivity_rows_1985_2024 = coalesce(
      valid_productivity_rows_1985_2024, 0L
    ),
    manifest_available = coalesce(manifest_available, FALSE),
    analytical_rows = coalesce(analytical_rows, 0L),
    climate_years = coalesce(climate_years, 0L),
    included = analytical_rows > 0,
    exclusion_reason = case_when(
      included ~ "Included",
      valid_productivity_rows_1985_2024 == 0 ~
        "No valid positive-area productivity observation from 1985 to 2024",
      !manifest_available ~ "No NASA POWER spatial manifest match",
      TRUE ~ "No complete climate and AWC100 water-balance match"
    )
  )

unit_inclusion_summary <- unit_inclusion_audit |>
  count(country, included, exclusion_reason, name = "spatial_units") |>
  group_by(country) |>
  mutate(
    country_units = sum(spatial_units),
    fraction_of_country_units = spatial_units / country_units
  ) |>
  ungroup()

local_panel <- local_panel |>
  group_by(country) |>
  mutate(
    area_weight_cap = quantile(
      harvested_area_ha, weight_cap_probability,
      na.rm = TRUE, names = FALSE, type = 8
    ),
    model_weight = pmin(harvested_area_ha, area_weight_cap),
    model_weight = model_weight / mean(model_weight, na.rm = TRUE)
  ) |>
  ungroup()

if (nrow(local_panel) == 0) stop("The local climate and soybean panels did not match")
if (any(!complete.cases(local_panel[, c(
  "yield_anomaly_percent", "oni_djf_c", "enso_warm", "enso_cold",
  "model_weight", all_mechanism_predictors
)]))) {
  stop("The assembled local panel contains missing analytical values", call. = FALSE)
}


# 6. NESTED TEMPORAL VALIDATION ------------------------------------------------

message("Running nested temporal validation at municipal and departmental scale...")

model_predictor_sets <- list(
  ENSO = c("enso_warm", "enso_cold"),
  climate = climate_predictors,
  water_balance = water_predictors,
  climate_plus_water = all_mechanism_predictors,
  ENSO_plus_climate = c("enso_warm", "enso_cold", climate_predictors),
  ENSO_plus_climate_water = c(
    "enso_warm", "enso_cold", all_mechanism_predictors
  )
)

validation_checkpoint_files <- c(
  oof = file.path(output_dir, "checkpoint_12_v6_oof_predictions.csv.gz"),
  tuning = file.path(output_dir, "checkpoint_12_v6_hyperparameters.csv"),
  coefficients = file.path(output_dir, "checkpoint_12_v6_coefficients.csv.gz"),
  permutations = file.path(output_dir, "checkpoint_12_v6_index_permutations.csv.gz")
)

if (reuse_existing_checkpoints && all(file.exists(validation_checkpoint_files))) {
  message("Reusing completed nested-validation checkpoint...")
  oof_predictions <- read_csv(validation_checkpoint_files[["oof"]], show_col_types = FALSE)
  selected_hyperparameters <- read_csv(
    validation_checkpoint_files[["tuning"]], show_col_types = FALSE
  )
  spatial_coefficients <- read_csv(
    validation_checkpoint_files[["coefficients"]], show_col_types = FALSE
  )
  permutation_distribution <- read_csv(
    validation_checkpoint_files[["permutations"]], show_col_types = FALSE
  )
} else {
oof_results <- list()
tuning_results <- list()
coefficient_results <- list()
importance_results <- list()
result_position <- 1L
tuning_position <- 1L
coefficient_position <- 1L
importance_position <- 1L

for (selected_country in c("Argentina", "Brazil")) {
  message("Processing ", selected_country, "...")
  country_data <- local_panel |> filter(country == selected_country)
  fold_table <- make_phase_stratified_folds(
    country_data$analysis_year,
    country_data$enso_phase_simple,
    outer_folds
  )
  country_data <- country_data |> left_join(fold_table, by = "analysis_year")

  for (fold_id in seq_len(outer_folds)) {
    message("  Outer fold ", fold_id, " of ", outer_folds)
    training <- country_data |> filter(fold != fold_id)
    testing <- country_data |> filter(fold == fold_id)
    null_prediction <- predict_unit_baseline(
      training, testing, "yield_anomaly_percent"
    )

    fold_predictions <- testing |>
      transmute(
        country, spatial_unit_id, parent_region_name, analysis_year,
        enso_phase_simple, observed = yield_anomaly_percent,
        model_weight, harvested_area_ha, outer_fold = fold_id,
        prediction_null = null_prediction
      )

    fitted_models <- list()
    for (model_name in names(model_predictor_sets)) {
      predictors <- model_predictor_sets[[model_name]]
      selected_lambda <- tune_lambda(
        training, predictors, "yield_anomaly_percent"
      )
      fit <- fit_hierarchical_ridge(
        training, predictors, "yield_anomaly_percent", selected_lambda
      )
      prediction <- predict_hierarchical_ridge(fit, testing)
      prediction_column <- paste0("prediction_", model_name)
      fold_predictions[[prediction_column]] <- prediction
      fitted_models[[model_name]] <- fit

      tuning_results[[tuning_position]] <- tibble(
        country = selected_country,
        outer_fold = fold_id,
        model = model_name,
        selected_lambda = selected_lambda,
        training_rows = nrow(training),
        testing_rows = nrow(testing)
      )
      tuning_position <- tuning_position + 1L

      coefficient_results[[coefficient_position]] <-
        extract_effective_coefficients(
          fit, selected_country, model_name, fold_id
        )
      coefficient_position <- coefficient_position + 1L
    }

    climate_fit <- fitted_models[["climate_plus_water"]]
    baseline_climate_mse <- safe_weighted_mse(
      testing$yield_anomaly_percent,
      fold_predictions$prediction_climate_plus_water,
      testing$model_weight
    )
    null_mse <- safe_weighted_mse(
      testing$yield_anomaly_percent,
      fold_predictions$prediction_null,
      testing$model_weight
    )

    for (index_predictor in all_mechanism_predictors) {
      for (permutation_id in seq_len(permutation_repetitions)) {
        permutation_prediction <- predict_hierarchical_ridge(
          climate_fit,
          testing,
          permute_predictors = index_predictor,
          permutation_seed = random_seed + 100000L * fold_id +
            1000L * match(selected_country, c("Argentina", "Brazil")) +
            10L * match(index_predictor, all_mechanism_predictors) +
            permutation_id
        )
        permutation_mse <- safe_weighted_mse(
          testing$yield_anomaly_percent,
          permutation_prediction,
          testing$model_weight
        )
        importance_results[[importance_position]] <- tibble(
          country = selected_country,
          outer_fold = fold_id,
          predictor = index_predictor,
          permutation_id = permutation_id,
          baseline_climate_mse = baseline_climate_mse,
          permutation_mse = permutation_mse,
          null_mse = null_mse,
          mse_increase = permutation_mse - baseline_climate_mse,
          importance_relative_to_null_error =
            (permutation_mse - baseline_climate_mse) / null_mse
        )
        importance_position <- importance_position + 1L
      }
    }

    oof_results[[result_position]] <- fold_predictions
    result_position <- result_position + 1L
  }
}

oof_predictions <- bind_rows(oof_results)
selected_hyperparameters <- bind_rows(tuning_results)
spatial_coefficients <- bind_rows(coefficient_results)
permutation_distribution <- bind_rows(importance_results)

write_csv(oof_predictions, validation_checkpoint_files[["oof"]], na = "")
write_csv(selected_hyperparameters, validation_checkpoint_files[["tuning"]], na = "")
write_csv(spatial_coefficients, validation_checkpoint_files[["coefficients"]], na = "")
write_csv(permutation_distribution, validation_checkpoint_files[["permutations"]], na = "")
}


# 7. COUNTRY, REGION, AND LOCAL ATTRIBUTION -----------------------------------

prediction_columns <- c(
  null = "prediction_null",
  ENSO = "prediction_ENSO",
  climate = "prediction_climate",
  water_balance = "prediction_water_balance",
  climate_plus_water = "prediction_climate_plus_water",
  ENSO_plus_climate = "prediction_ENSO_plus_climate",
  ENSO_plus_climate_water = "prediction_ENSO_plus_climate_water"
)

summarise_predictions <- function(data) {
  null_mse <- safe_weighted_mse(
    data$observed, data$prediction_null, data$model_weight
  )
  map_dfr(names(prediction_columns), function(model_name) {
    prediction <- data[[prediction_columns[[model_name]]]]
    model_mse <- safe_weighted_mse(data$observed, prediction, data$model_weight)
    tibble(
      model = model_name,
      observations = sum(is.finite(data$observed) & is.finite(prediction)),
      years = n_distinct(data$analysis_year[is.finite(prediction)]),
      validation_mse = model_mse,
      out_of_sample_r_squared = safe_r_squared(
        data$observed, prediction, data$model_weight
      ),
      error_reduction_relative_to_local_null = if (
        model_name == "null"
      ) 0 else 1 - model_mse / null_mse
    )
  })
}

country_attribution <- oof_predictions |>
  group_by(country) |>
  group_modify(~ summarise_predictions(.x)) |>
  ungroup()

regional_attribution <- oof_predictions |>
  group_by(country, parent_region_name) |>
  group_modify(~ summarise_predictions(.x)) |>
  ungroup()

local_attribution <- oof_predictions |>
  group_by(country, spatial_unit_id, parent_region_name) |>
  group_modify(~ summarise_predictions(.x)) |>
  ungroup() |>
  filter(observations >= minimum_local_validation_rows)

national_annual_oof <- oof_predictions |>
  group_by(country, analysis_year, enso_phase_simple) |>
  summarise(
    observed = safe_weighted_mean(observed, harvested_area_ha),
    across(
      all_of(unname(prediction_columns)),
      ~ safe_weighted_mean(.x, harvested_area_ha)
    ),
    harvested_area_ha = sum(harvested_area_ha, na.rm = TRUE),
    spatial_units = n_distinct(spatial_unit_id),
    .groups = "drop"
  )

national_net_attribution <- national_annual_oof |>
  group_by(country) |>
  group_modify(function(data, keys) {
    data$model_weight <- data$harvested_area_ha
    summarise_predictions(data)
  }) |>
  ungroup()

country_model_comparison <- country_attribution |>
  select(country, model, validation_mse, out_of_sample_r_squared) |>
  pivot_wider(
    names_from = model,
    values_from = c(validation_mse, out_of_sample_r_squared)
  ) |>
  mutate(
    water_gain_over_raw_climate =
      1 - validation_mse_climate_plus_water / validation_mse_climate,
    raw_climate_gain_over_water =
      1 - validation_mse_climate_plus_water / validation_mse_water_balance,
    ENSO_gain_beyond_combined_mechanisms =
      1 - validation_mse_ENSO_plus_climate_water /
        validation_mse_climate_plus_water
  )

permutation_importance <- permutation_distribution |>
  group_by(country, predictor) |>
  summarise(
    permutations = n(),
    mean_mse_increase = mean(mse_increase, na.rm = TRUE),
    median_mse_increase = median(mse_increase, na.rm = TRUE),
    importance_relative_to_null_error = mean(
      importance_relative_to_null_error, na.rm = TRUE
    ),
    positive_importance_probability = mean(mse_increase > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(country) |>
  mutate(
    positive_importance = pmax(importance_relative_to_null_error, 0),
    importance_denominator = sum(positive_importance, na.rm = TRUE),
    normalized_positive_importance = if_else(
      importance_denominator > 0,
      positive_importance / importance_denominator,
      0
    )
  ) |>
  ungroup() |>
  select(-importance_denominator)


# 8. CROSS-VALIDATED ENSO TO CLIMATE LINKS -----------------------------------

message("Estimating the cross-validated ENSO signal in each harmonized index...")

enso_climate_checkpoint <- file.path(
  output_dir, "checkpoint_12_v6_ENSO_to_indices_oof.csv.gz"
)

if (reuse_existing_checkpoints && file.exists(enso_climate_checkpoint)) {
  message("Reusing completed ENSO-to-climate checkpoint...")
  enso_climate_oof <- read_csv(enso_climate_checkpoint, show_col_types = FALSE)
} else {
enso_climate_oof_results <- list()
enso_climate_position <- 1L

for (selected_country in c("Argentina", "Brazil")) {
  country_data <- local_panel |> filter(country == selected_country)
  fold_table <- make_phase_stratified_folds(
    country_data$analysis_year,
    country_data$enso_phase_simple,
    outer_folds
  )
  country_data <- country_data |> left_join(fold_table, by = "analysis_year")

  for (climate_predictor in all_mechanism_predictors) {
    for (fold_id in seq_len(outer_folds)) {
      training <- country_data |> filter(fold != fold_id)
      testing <- country_data |> filter(fold == fold_id)
      selected_lambda <- tune_lambda(
        training, c("enso_warm", "enso_cold"), climate_predictor
      )
      fit <- fit_hierarchical_ridge(
        training,
        c("enso_warm", "enso_cold"),
        climate_predictor,
        selected_lambda
      )
      enso_prediction <- predict_hierarchical_ridge(fit, testing)
      null_prediction <- predict_unit_baseline(
        training, testing, climate_predictor
      )
      enso_climate_oof_results[[enso_climate_position]] <- tibble(
        country = selected_country,
        predictor = climate_predictor,
        spatial_unit_id = testing$spatial_unit_id,
        parent_region_name = testing$parent_region_name,
        analysis_year = testing$analysis_year,
        observed_climate = testing[[climate_predictor]],
        prediction_null = null_prediction,
        prediction_ENSO = enso_prediction,
        model_weight = testing$model_weight,
        outer_fold = fold_id,
        selected_lambda = selected_lambda
      )
      enso_climate_position <- enso_climate_position + 1L
    }
  }
}

enso_climate_oof <- bind_rows(enso_climate_oof_results)
write_csv(enso_climate_oof, enso_climate_checkpoint, na = "")
}

enso_to_climate_skill <- enso_climate_oof |>
  group_by(country, predictor) |>
  summarise(
    observations = n(),
    years = n_distinct(analysis_year),
    null_mse = safe_weighted_mse(
      observed_climate, prediction_null, model_weight
    ),
    enso_mse = safe_weighted_mse(
      observed_climate, prediction_ENSO, model_weight
    ),
    enso_climate_out_of_sample_r_squared = safe_r_squared(
      observed_climate, prediction_ENSO, model_weight
    ),
    enso_error_reduction_relative_to_null = 1 - enso_mse / null_mse,
    .groups = "drop"
  )

mediated_predictive_scores <- permutation_importance |>
  left_join(
    enso_to_climate_skill,
    by = c("country", "predictor")
  ) |>
  mutate(
    positive_ENSO_climate_skill = pmax(
      enso_error_reduction_relative_to_null, 0
    ),
    ENSO_aligned_predictive_score = positive_importance *
      positive_ENSO_climate_skill
  ) |>
  group_by(country) |>
  mutate(
    score_denominator = sum(ENSO_aligned_predictive_score, na.rm = TRUE),
    normalized_ENSO_aligned_predictive_score = if_else(
      score_denominator > 0,
      ENSO_aligned_predictive_score / score_denominator,
      0
    )
  ) |>
  ungroup() |>
  select(-score_denominator)


# 9. CIRCULAR SHIFT TESTS FOR ENSO TO CLIMATE ALIGNMENT -----------------------

message("Running exact nonzero circular shifts for ENSO to climate alignment...")

state_year_panel <- local_panel |>
  group_by(
    country, parent_region_name, analysis_year,
    enso_warm, enso_cold, enso_phase_simple
  ) |>
  summarise(
    across(
      all_of(all_mechanism_predictors),
      ~ safe_weighted_mean(.x, model_weight)
    ),
    model_weight = sum(model_weight, na.rm = TRUE),
    .groups = "drop"
  )

circular_shift_results <- map_dfr(c("Argentina", "Brazil"), function(selected_country) {
  country_data <- state_year_panel |> filter(country == selected_country)
  years <- sort(unique(country_data$analysis_year))
  oni_table <- country_data |>
    distinct(analysis_year, enso_warm, enso_cold) |>
    arrange(analysis_year)

  map_dfr(all_mechanism_predictors, function(climate_predictor) {
    null_formula <- as.formula(
      paste(climate_predictor, "~ factor(parent_region_name)")
    )
    full_formula <- as.formula(paste(
      climate_predictor,
      "~ factor(parent_region_name) +",
      "enso_warm:factor(parent_region_name) +",
      "enso_cold:factor(parent_region_name)"
    ))
    null_mse <- fit_weighted_formula_mse(null_formula, country_data)
    observed_full_mse <- fit_weighted_formula_mse(full_formula, country_data)
    observed_partial_r2 <- 1 - observed_full_mse / null_mse

    null_values <- map_dbl(seq_len(length(years) - 1L), function(shift) {
      shifted_oni <- oni_table |>
        transmute(
          analysis_year,
          enso_warm_shifted = dplyr::lead(enso_warm, shift),
          enso_cold_shifted = dplyr::lead(enso_cold, shift)
        )
      shifted_oni$enso_warm_shifted <- c(
        oni_table$enso_warm[(shift + 1L):length(years)],
        oni_table$enso_warm[seq_len(shift)]
      )
      shifted_oni$enso_cold_shifted <- c(
        oni_table$enso_cold[(shift + 1L):length(years)],
        oni_table$enso_cold[seq_len(shift)]
      )
      shifted_data <- country_data |>
        select(-enso_warm, -enso_cold) |>
        left_join(shifted_oni, by = "analysis_year") |>
        rename(
          enso_warm = enso_warm_shifted,
          enso_cold = enso_cold_shifted
        )
      shifted_mse <- fit_weighted_formula_mse(full_formula, shifted_data)
      1 - shifted_mse / null_mse
    })

    tibble(
      country = selected_country,
      predictor = climate_predictor,
      observed_partial_r_squared = observed_partial_r2,
      nonzero_circular_shifts = length(null_values),
      circular_shift_p_value =
        (1 + sum(null_values >= observed_partial_r2, na.rm = TRUE)) /
        (1 + sum(is.finite(null_values))),
      null_median = median(null_values, na.rm = TRUE),
      null_95th_percentile = quantile(
        null_values, 0.95, na.rm = TRUE, names = FALSE, type = 8
      )
    )
  })
})


# 10. QUALITY CONTROL AND METHODOLOGY -----------------------------------------

qc_summary <- tibble(
  test = c(
    "Local panel rows",
    "Countries",
    "Spatial units",
    "Regions",
    "Harvest years",
    "Missing analytical values",
    "Duplicated country unit year rows",
    "OOF prediction rows",
    "Minimum model prediction coverage",
    "Local units retained in validation",
    "Permutation evaluations",
    "ENSO climate validation rows"
  ),
  value = c(
    nrow(local_panel),
    n_distinct(local_panel$country),
    n_distinct(interaction(local_panel$country, local_panel$spatial_unit_id)),
    n_distinct(interaction(local_panel$country, local_panel$parent_region_name)),
    n_distinct(local_panel$analysis_year),
    sum(!complete.cases(local_panel[, c(
      "yield_anomaly_percent", "model_weight", "enso_warm", "enso_cold",
      all_mechanism_predictors
    )])),
    sum(duplicated(local_panel[, c("country", "spatial_unit_id", "analysis_year")])),
    nrow(oof_predictions),
    min(vapply(
      prediction_columns,
      function(column) mean(is.finite(oof_predictions[[column]])),
      numeric(1)
    )),
    n_distinct(interaction(local_attribution$country, local_attribution$spatial_unit_id)),
    nrow(permutation_distribution),
    nrow(enso_climate_oof)
  )
)

methodology <- tribble(
  ~component, ~definition,
  "Observation scale", "Municipality-year in Brazil and department-year in Argentina",
  "Outcome", "Detrended soybean productivity anomaly in percent",
  "Weather window", "Reproductive crop-calendar window from script 11",
  "Local baseline", "Spatial-unit mean estimated only from training years",
  "Hierarchy", "National climate slope plus ridge-shrunk state or province deviations",
  "Weight", "Harvested area capped at the country-specific 99th percentile",
  "Outer validation", "Five phase-stratified folds holding out complete harvest years",
  "Inner validation", "Three phase-stratified folds selecting the regional ridge penalty",
  "ENSO model", "Separate positive and negative DJF ONI components",
  "Index set", "Seven harmonized reproductive-window indices: temperature, VPD, precipitation, heavy rainfall, CWS, hot-dry conditions, and hot-wet conditions",
  "CWS orientation", "CWS is the sign-reversed and standardized AWC100 reproductive requirement satisfaction index; higher values indicate greater crop water stress",
  "Combined model", "The same seven indices are used jointly over the common 1985 to 2024 period",
  "AWC sensitivity", "AWC100 is primary; AWC75 and AWC125 remain available for a separate robustness analysis",
  "Permutation importance", "Increase in held-out MSE after within-region permutation of each of the seven indices separately",
  "ENSO to index skill", "Held-out error reduction for predicting each of the seven indices from ONI",
  "ENSO-aligned score", "Positive yield importance multiplied by positive ENSO-to-climate skill",
  "Circular shift", "All nonzero temporal rotations of ONI against the regional climate field",
  "Interpretation", "Predictive attribution and mechanistic consistency, not causal identification"
)


# 11. SAVE OUTPUTS -------------------------------------------------------------

message("Saving spatial climate attribution outputs...")

write_csv(
  oof_predictions,
  file.path(output_dir, "ENSO_local_climate_model_oof_predictions.csv.gz"),
  na = ""
)
write_csv(
  selected_hyperparameters,
  file.path(output_dir, "ENSO_local_climate_selected_hyperparameters.csv"),
  na = ""
)
write_csv(
  spatial_coefficients,
  file.path(output_dir, "ENSO_local_climate_regional_coefficients.csv.gz"),
  na = ""
)
write_csv(
  country_attribution,
  file.path(output_dir, "ENSO_climate_attribution_country_summary.csv"),
  na = ""
)
write_csv(
  national_net_attribution,
  file.path(output_dir, "ENSO_climate_attribution_national_net_summary.csv"),
  na = ""
)
write_csv(
  country_model_comparison,
  file.path(output_dir, "ENSO_climate_water_model_comparison.csv"),
  na = ""
)
write_csv(
  regional_attribution,
  file.path(output_dir, "ENSO_climate_attribution_region_summary.csv"),
  na = ""
)
write_csv(
  local_attribution,
  file.path(output_dir, "ENSO_climate_attribution_local_summary.csv.gz"),
  na = ""
)
write_csv(
  national_annual_oof,
  file.path(output_dir, "ENSO_climate_attribution_national_annual_oof.csv"),
  na = ""
)
write_csv(
  permutation_distribution,
  file.path(output_dir, "ENSO_climate_permutation_distribution.csv.gz"),
  na = ""
)
write_csv(
  permutation_importance,
  file.path(output_dir, "ENSO_climate_variable_permutation_importance.csv"),
  na = ""
)
write_csv(
  enso_climate_oof,
  file.path(output_dir, "ENSO_to_climate_oof_predictions.csv.gz"),
  na = ""
)
write_csv(
  enso_to_climate_skill,
  file.path(output_dir, "ENSO_to_climate_crossvalidated_skill.csv"),
  na = ""
)
write_csv(
  mediated_predictive_scores,
  file.path(output_dir, "ENSO_climate_mediated_predictive_scores.csv"),
  na = ""
)
write_csv(
  circular_shift_results,
  file.path(output_dir, "ENSO_climate_circular_shift_tests.csv"),
  na = ""
)
write_csv(
  qc_summary,
  file.path(output_dir, "QC_spatial_climate_attribution.csv"),
  na = ""
)
write_csv(
  methodology,
  file.path(output_dir, "README_spatial_climate_attribution_metrics.csv"),
  na = ""
)
write_csv(
  unit_inclusion_audit,
  file.path(output_dir, "QC_spatial_climate_unit_inclusion.csv"),
  na = ""
)
write_csv(
  unit_inclusion_summary,
  file.path(output_dir, "QC_spatial_climate_unit_inclusion_summary.csv"),
  na = ""
)


# 12. FINAL VALIDATION ---------------------------------------------------------

if (qc_summary$value[qc_summary$test == "Missing analytical values"] > 0) {
  stop("Analytical values are missing after panel assembly", call. = FALSE)
}
if (qc_summary$value[qc_summary$test == "Duplicated country unit year rows"] > 0) {
  stop("Duplicated spatial unit-year rows were detected", call. = FALSE)
}
if (qc_summary$value[
  qc_summary$test == "Minimum model prediction coverage"
] < 0.99) {
  stop("At least one model produced less than 99% OOF prediction coverage", call. = FALSE)
}
if (any(!is.finite(country_attribution$out_of_sample_r_squared))) {
  stop("A country attribution metric is non-finite", call. = FALSE)
}
if (any(!is.finite(enso_to_climate_skill$enso_error_reduction_relative_to_null))) {
  stop("An ENSO-to-climate validation metric is non-finite", call. = FALSE)
}

message("")
message("Script 12 completed successfully.")
message("Local panel rows: ", nrow(local_panel))
message(
  "Spatial units: ",
  n_distinct(interaction(local_panel$country, local_panel$spatial_unit_id))
)
message("Harvest years: ", min(local_panel$analysis_year), " to ", max(local_panel$analysis_year))
message("Output directory: ", normalizePath(output_dir, winslash = "/"))

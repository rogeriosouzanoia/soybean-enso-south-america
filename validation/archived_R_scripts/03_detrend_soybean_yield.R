# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 03_detrend_soybean_yield.R
# PURPOSE: Remove technological trends from municipal and departmental soybean
#          productivity series in Brazil and Argentina
# VERSION: 1.1
# IMPORTANT: The quality controlled source file is never overwritten
# ==============================================================================


# 1. USER SETTINGS -------------------------------------------------------------

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")

project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) {
  project_dir_override
} else {
  default_project_dir
}

minimum_valid_years <- 20L
number_cv_folds <- 5L

# A simpler method is retained when its cross validated RMSE is no more than
# 2 percent above the best RMSE. This reduces unnecessary trend flexibility.
selection_tolerance <- 0.02

# The constant method is the null reference representing no temporal trend.
# LOESS span controls smoothness. Spline df defines the effective degrees of
# freedom used for the technological trend.
candidate_methods <- c(
  "constant",
  "linear",
  "quadratic",
  "spline_df4",
  "spline_df6",
  "loess_span100",
  "loess_span075",
  "loess_span050"
)

# Lower values indicate simpler models and determine preference inside the
# selection tolerance.
method_complexity <- c(
  constant = 0L,
  linear = 1L,
  quadratic = 2L,
  spline_df4 = 3L,
  loess_span100 = 4L,
  spline_df6 = 5L,
  loess_span075 = 6L,
  loess_span050 = 7L
)

# A fitted technological trend at or below this value is considered invalid.
minimum_valid_trend_t_ha <- 0.05

# Minimum fraction of held out observations that must receive predictions.
minimum_cv_prediction_coverage <- 0.80


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

input_dir <- file.path(project_dir, "data_qc")
output_dir <- file.path(project_dir, "data_detrended")
metadata_dir <- file.path(project_dir, "metadata")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(
  input_dir,
  "soybean_Brazil_Argentina_analysis_ready_yield.csv"
)

if (!file.exists(input_file)) {
  stop("Analysis ready productivity file not found: ", input_file, call. = FALSE)
}


# 4. HELPER FUNCTIONS ----------------------------------------------------------

check_required_columns <- function(data, required) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      "The input dataset is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

rmse <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  if (!any(valid)) return(NA_real_)
  sqrt(mean((observed[valid] - predicted[valid])^2))
}

mae <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  if (!any(valid)) return(NA_real_)
  mean(abs(observed[valid] - predicted[valid]))
}

r_squared <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  if (sum(valid) < 3) return(NA_real_)
  denominator <- sum((observed[valid] - mean(observed[valid]))^2)
  if (denominator <= 0) return(NA_real_)
  1 - sum((observed[valid] - predicted[valid])^2) / denominator
}

fit_and_predict <- function(method, train, new_year) {
  train <- train |>
    filter(is.finite(analysis_year), is.finite(yield_harvested_t_ha)) |>
    arrange(analysis_year)

  if (nrow(train) < 8 || dplyr::n_distinct(train$analysis_year) < 8) {
    return(rep(NA_real_, length(new_year)))
  }

  year_center <- mean(train$analysis_year)
  train <- train |>
    mutate(year_centered = analysis_year - year_center)
  new_data <- tibble(
    analysis_year = new_year,
    year_centered = new_year - year_center
  )

  prediction <- tryCatch({
    if (method == "constant") {
      rep(mean(train$yield_harvested_t_ha, na.rm = TRUE), length(new_year))

    } else if (method == "linear") {
      model <- lm(yield_harvested_t_ha ~ year_centered, data = train)
      predict(model, newdata = new_data)

    } else if (method == "quadratic") {
      model <- lm(
        yield_harvested_t_ha ~ year_centered + I(year_centered^2),
        data = train
      )
      predict(model, newdata = new_data)

    } else if (method == "spline_df4") {
      model <- smooth.spline(
        x = train$analysis_year,
        y = train$yield_harvested_t_ha,
        df = 4
      )
      predict(model, x = new_year)$y

    } else if (method == "spline_df6") {
      model <- smooth.spline(
        x = train$analysis_year,
        y = train$yield_harvested_t_ha,
        df = 6
      )
      predict(model, x = new_year)$y

    } else if (method == "loess_span100") {
      model <- loess(
        yield_harvested_t_ha ~ analysis_year,
        data = train,
        span = 1.00,
        degree = 1,
        control = loess.control(surface = "direct")
      )
      predict(model, newdata = new_data)

    } else if (method == "loess_span075") {
      model <- loess(
        yield_harvested_t_ha ~ analysis_year,
        data = train,
        span = 0.75,
        degree = 1,
        control = loess.control(surface = "direct")
      )
      predict(model, newdata = new_data)

    } else if (method == "loess_span050") {
      model <- loess(
        yield_harvested_t_ha ~ analysis_year,
        data = train,
        span = 0.50,
        degree = 1,
        control = loess.control(surface = "direct")
      )
      predict(model, newdata = new_data)

    } else {
      stop("Unknown detrending method: ", method)
    }
  }, error = function(e) {
    rep(NA_real_, length(new_year))
  })

  as.numeric(prediction)
}


# 5. TEMPORAL CROSS VALIDATION -------------------------------------------------

# Observations are sorted by year and assigned systematically to folds. Each
# fold therefore samples the full historical period while every prediction is
# generated without using the corresponding observation in model fitting.
cross_validate_method <- function(unit_data, method, k = 5L) {
  unit_data <- unit_data |>
    arrange(analysis_year) |>
    mutate(cv_fold = ((row_number() - 1L) %% k) + 1L)

  predictions <- rep(NA_real_, nrow(unit_data))

  for (fold in seq_len(k)) {
    test_index <- which(unit_data$cv_fold == fold)
    train_data <- unit_data[-test_index, , drop = FALSE]
    predictions[test_index] <- fit_and_predict(
      method = method,
      train = train_data,
      new_year = unit_data$analysis_year[test_index]
    )
  }

  valid_prediction <- is.finite(predictions)
  prediction_coverage <- mean(valid_prediction)

  tibble(
    method = method,
    cv_rmse_t_ha = rmse(unit_data$yield_harvested_t_ha, predictions),
    cv_mae_t_ha = mae(unit_data$yield_harvested_t_ha, predictions),
    cv_r_squared = r_squared(unit_data$yield_harvested_t_ha, predictions),
    cv_prediction_coverage = prediction_coverage,
    cv_predictions = sum(valid_prediction)
  )
}


# 6. FIT AND SELECT TREND FOR ONE SPATIAL UNIT --------------------------------

process_spatial_unit <- function(unit_data) {
  unit_data <- unit_data |>
    arrange(analysis_year)

  model_metrics <- bind_rows(lapply(
    candidate_methods,
    function(method) cross_validate_method(unit_data, method, number_cv_folds)
  )) |>
    mutate(
      complexity = unname(method_complexity[method]),
      cv_eligible = is.finite(cv_rmse_t_ha) &
        cv_prediction_coverage >= minimum_cv_prediction_coverage
    )

  final_predictions <- lapply(candidate_methods, function(method) {
    fit_and_predict(method, unit_data, unit_data$analysis_year)
  })
  names(final_predictions) <- candidate_methods

  final_prediction_table <- as_tibble(final_predictions) |>
    rename_with(~ paste0("trend_", .x, "_t_ha"))

  final_validity <- tibble(
    method = candidate_methods,
    final_fit_valid = vapply(
      final_predictions,
      function(x) {
        all(is.finite(x)) && all(x > minimum_valid_trend_t_ha)
      },
      logical(1)
    )
  )

  model_metrics <- model_metrics |>
    left_join(final_validity, by = "method") |>
    mutate(selection_eligible = cv_eligible & final_fit_valid)

  eligible_metrics <- model_metrics |>
    filter(selection_eligible)

  if (nrow(eligible_metrics) == 0) {
    selected_method <- NA_character_
    selection_status <- "no_valid_candidate"
    selected_prediction <- rep(NA_real_, nrow(unit_data))
  } else {
    best_rmse <- min(eligible_metrics$cv_rmse_t_ha, na.rm = TRUE)
    selected_method <- eligible_metrics |>
      filter(cv_rmse_t_ha <= best_rmse * (1 + selection_tolerance)) |>
      arrange(complexity, cv_rmse_t_ha) |>
      slice(1) |>
      pull(method)
    selection_status <- "selected"
    selected_prediction <- final_predictions[[selected_method]]
  }

  observed <- unit_data$yield_harvested_t_ha
  residual <- observed - selected_prediction
  residual_sd <- sd(residual, na.rm = TRUE)

  standardized_residual <- if (
    is.finite(residual_sd) && residual_sd > 0
  ) {
    residual / residual_sd
  } else {
    rep(NA_real_, length(residual))
  }

  result_data <- bind_cols(unit_data, final_prediction_table) |>
    mutate(
      selected_trend_method = selected_method,
      temporal_trend_detected = !is.na(selected_method) &
        selected_method != "constant",
      selected_trend_t_ha = selected_prediction,
      yield_anomaly_t_ha = residual,
      yield_anomaly_fraction = if_else(
        is.finite(selected_prediction) & selected_prediction > 0,
        residual / selected_prediction,
        NA_real_
      ),
      yield_anomaly_percent = 100 * yield_anomaly_fraction,
      yield_ratio_observed_to_trend = if_else(
        is.finite(selected_prediction) & selected_prediction > 0,
        observed / selected_prediction,
        NA_real_
      ),
      yield_anomaly_standardized = standardized_residual
    )

  selected_metrics <- model_metrics |>
    mutate(
      selected = !is.na(selected_method) & method == selected_method
    )

  unit_summary <- tibble(
    selected_trend_method = selected_method,
    temporal_trend_detected = !is.na(selected_method) &
      selected_method != "constant",
    selection_status = selection_status,
    observations = nrow(unit_data),
    first_year = min(unit_data$analysis_year),
    last_year = max(unit_data$analysis_year),
    observed_mean_t_ha = mean(observed, na.rm = TRUE),
    observed_sd_t_ha = sd(observed, na.rm = TRUE),
    residual_mean_t_ha = mean(residual, na.rm = TRUE),
    residual_sd_t_ha = residual_sd,
    selected_cv_rmse_t_ha = selected_metrics |>
      filter(selected) |>
      pull(cv_rmse_t_ha) |>
      dplyr::first(default = NA_real_),
    selected_cv_mae_t_ha = selected_metrics |>
      filter(selected) |>
      pull(cv_mae_t_ha) |>
      dplyr::first(default = NA_real_),
    selected_cv_r_squared = selected_metrics |>
      filter(selected) |>
      pull(cv_r_squared) |>
      dplyr::first(default = NA_real_)
  )

  list(
    data = result_data,
    metrics = selected_metrics,
    summary = unit_summary
  )
}


# 7. READ AND VALIDATE INPUT ---------------------------------------------------

message("Reading analysis ready productivity data...")

soybean <- read_csv(
  input_file,
  col_types = cols(
    spatial_unit_id = col_character(),
    parent_region_code = col_character(),
    analysis_year = col_integer(),
    yield_harvested_t_ha = col_double(),
    .default = col_guess()
  ),
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
)

required_columns <- c(
  "country", "spatial_unit_id", "spatial_unit_name", "parent_region_name",
  "spatial_level", "analysis_year", "harvested_area_ha", "production_t",
  "yield_harvested_t_ha", "analysis_ready_yield"
)
check_required_columns(soybean, required_columns)

soybean <- soybean |>
  filter(
    analysis_ready_yield,
    is.finite(analysis_year),
    is.finite(yield_harvested_t_ha)
  ) |>
  arrange(country, spatial_unit_id, analysis_year)

if (anyDuplicated(soybean[c("country", "spatial_unit_id", "analysis_year")]) > 0) {
  stop("Duplicated spatial unit and year keys were found.", call. = FALSE)
}

coverage_check <- soybean |>
  count(country, spatial_unit_id, name = "valid_years")

if (any(coverage_check$valid_years < minimum_valid_years)) {
  stop("At least one spatial unit has fewer than the required valid years.", call. = FALSE)
}


# 8. PROCESS ALL SPATIAL UNITS -------------------------------------------------

unit_list <- soybean |>
  group_by(country, spatial_unit_id) |>
  group_split(.keep = TRUE)

number_units <- length(unit_list)
message("Spatial units to process: ", number_units)
message("This step can take several minutes.")

processed <- vector("list", number_units)

for (i in seq_along(unit_list)) {
  if (i == 1L || i %% 100L == 0L || i == number_units) {
    message("Processing spatial unit ", i, " of ", number_units, "...")
  }
  processed[[i]] <- process_spatial_unit(unit_list[[i]])
}


# 9. COMBINE RESULTS -----------------------------------------------------------

detrended_all_methods <- bind_rows(lapply(
  processed,
  function(result) result$data
))

model_selection <- bind_rows(lapply(seq_along(processed), function(i) {
  identity_columns <- unit_list[[i]] |>
    slice(1) |>
    select(
      country, spatial_unit_id, spatial_unit_name, parent_region_name,
      spatial_level
    )
  bind_cols(identity_columns[rep(1, nrow(processed[[i]]$metrics)), ],
            processed[[i]]$metrics)
}))

unit_summary <- bind_rows(lapply(seq_along(processed), function(i) {
  identity_columns <- unit_list[[i]] |>
    slice(1) |>
    select(
      country, spatial_unit_id, spatial_unit_name, parent_region_name,
      spatial_level
    )
  bind_cols(identity_columns, processed[[i]]$summary)
}))

detrended_selected <- detrended_all_methods |>
  select(
    country, spatial_unit_id, spatial_unit_name, parent_region_code,
    parent_region_name, spatial_level, season, analysis_year,
    planted_area_ha, harvested_area_ha, production_t,
    yield_harvested_t_ha, selected_trend_method, temporal_trend_detected,
    selected_trend_t_ha,
    yield_anomaly_t_ha, yield_anomaly_fraction, yield_anomaly_percent,
    yield_ratio_observed_to_trend, yield_anomaly_standardized
  )

detrending_failures <- unit_summary |>
  filter(selection_status != "selected")


# 10. FINAL QUALITY CONTROL ----------------------------------------------------

if (nrow(detrended_selected) != nrow(soybean)) {
  stop("The detrended output has a different number of rows from the input.", call. = FALSE)
}

if (anyDuplicated(
  detrended_selected[c("country", "spatial_unit_id", "analysis_year")]
) > 0) {
  stop("Duplicated keys were generated during detrending.", call. = FALSE)
}

if (nrow(detrending_failures) == 0 && any(
  !is.finite(detrended_selected$selected_trend_t_ha)
)) {
  stop("Missing selected trend values were generated.", call. = FALSE)
}

method_frequency <- unit_summary |>
  count(country, selected_trend_method, selection_status, name = "spatial_units") |>
  group_by(country) |>
  mutate(percent_spatial_units = 100 * spatial_units / sum(spatial_units)) |>
  ungroup()

country_diagnostic <- detrended_selected |>
  group_by(country) |>
  summarise(
    observations = n(),
    spatial_units = n_distinct(spatial_unit_id),
    first_year = min(analysis_year),
    last_year = max(analysis_year),
    mean_observed_yield_t_ha = mean(yield_harvested_t_ha, na.rm = TRUE),
    mean_absolute_anomaly_t_ha = mean(yield_anomaly_t_ha, na.rm = TRUE),
    sd_absolute_anomaly_t_ha = sd(yield_anomaly_t_ha, na.rm = TRUE),
    mean_relative_anomaly_percent = mean(yield_anomaly_percent, na.rm = TRUE),
    failed_spatial_units = n_distinct(
      spatial_unit_id[!is.finite(selected_trend_t_ha)]
    ),
    .groups = "drop"
  )


# 11. METADATA -----------------------------------------------------------------

parameter_metadata <- tibble(
  parameter = c(
    "minimum_valid_years",
    "number_cv_folds",
    "selection_tolerance",
    "minimum_valid_trend_t_ha",
    "minimum_cv_prediction_coverage",
    "candidate_methods"
  ),
  value = c(
    as.character(minimum_valid_years),
    as.character(number_cv_folds),
    as.character(selection_tolerance),
    as.character(minimum_valid_trend_t_ha),
    as.character(minimum_cv_prediction_coverage),
    paste(candidate_methods, collapse = ", ")
  ),
  description = c(
    "Minimum valid annual observations required for each spatial unit",
    "Number of systematic temporal cross validation folds",
    "Maximum proportional RMSE difference for selecting a simpler method",
    "Minimum biologically plausible fitted technological trend",
    "Minimum share of observations with cross validated predictions",
    "Null reference and technological trend methods evaluated for every spatial unit"
  )
)


# 12. SAVE OUTPUTS -------------------------------------------------------------

message("Saving detrended productivity outputs...")

write_csv(
  detrended_selected,
  file.path(output_dir, "soybean_Brazil_Argentina_detrended_selected.csv"),
  na = ""
)

write_csv(
  detrended_all_methods,
  file.path(output_dir, "soybean_Brazil_Argentina_detrended_all_methods.csv"),
  na = ""
)

write_csv(
  model_selection,
  file.path(output_dir, "QC_detrending_model_selection.csv"),
  na = ""
)

write_csv(
  unit_summary,
  file.path(output_dir, "QC_detrending_spatial_unit_summary.csv"),
  na = ""
)

write_csv(
  detrending_failures,
  file.path(output_dir, "QC_detrending_failures.csv"),
  na = ""
)

write_csv(
  method_frequency,
  file.path(output_dir, "QC_detrending_method_frequency.csv"),
  na = ""
)

write_csv(
  country_diagnostic,
  file.path(output_dir, "QC_detrending_country_diagnostic.csv"),
  na = ""
)

write_csv(
  parameter_metadata,
  file.path(metadata_dir, "detrending_analysis_parameters.csv"),
  na = ""
)

capture.output(
  sessionInfo(),
  file = file.path(metadata_dir, "detrending_R_session_info.txt")
)


# 13. FINAL REPORT -------------------------------------------------------------

message("")
message("Detrending completed successfully.")
message("Input observations: ", nrow(soybean))
message("Spatial units processed: ", number_units)
message("Spatial units without a valid trend: ", nrow(detrending_failures))
message("Output directory: ", normalizePath(output_dir, winslash = "/"))
message("")

print(method_frequency)
print(country_diagnostic)

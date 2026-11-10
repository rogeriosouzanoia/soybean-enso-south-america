# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 06_quantify_ENSO_yield_variability.R
# PURPOSE: Quantify the association between ENSO and interannual soybean
#          productivity variability at national, continental and local levels
# VERSION: 1.0
#
# INTERPRETATION:
#   The script estimates statistical variance explained by ENSO. These results
#   quantify association and predictive contribution, not causal attribution.
#   ONI in DJF is the preregistered primary exposure. Other windows are treated
#   as sensitivity analyses and are corrected for multiple window screening.
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

# Five contiguous temporal folds are used for out of sample validation.
cv_folds <- 5L

# Random year permutations preserve the spatial dependence among all units
# observed in the same harvest year.
permutation_repetitions <- 1999L
permutation_seed <- 20260904L

# Minimum number of observations for unit level inference.
minimum_unit_years <- 20L

# Primary ENSO exposure defined before screening alternative windows.
primary_index <- "oni_djf_c"


# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c("dplyr", "readr", "tidyr", "tibble", "purrr")

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
  library(purrr)
})


# 3. PATHS ---------------------------------------------------------------------

merged_dir <- file.path(project_dir, "data_merged")
output_dir <- file.path(project_dir, "enso_variance_models")
metadata_dir <- file.path(project_dir, "metadata")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

panel_file <- file.path(
  merged_dir,
  "soybean_ENSO_spatial_panel_1974_2024.csv"
)

country_year_file <- file.path(
  merged_dir,
  "soybean_ENSO_country_year_1974_2024.csv"
)

combined_year_file <- file.path(
  merged_dir,
  "soybean_ENSO_Brazil_Argentina_year_1974_2024.csv"
)

required_files <- c(panel_file, country_year_file, combined_year_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Required script 05 outputs were not found: ",
    paste(missing_files, collapse = "; "),
    call. = FALSE
  )
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

safe_r_squared <- function(observed, predicted, reference = NULL) {
  valid <- is.finite(observed) & is.finite(predicted)
  observed <- observed[valid]
  predicted <- predicted[valid]

  if (length(observed) < 3) return(NA_real_)

  if (is.null(reference)) {
    denominator <- sum((observed - mean(observed))^2)
  } else {
    reference <- reference[valid]
    denominator <- sum((observed - reference)^2)
  }

  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  1 - sum((observed - predicted)^2) / denominator
}

adjusted_r_squared <- function(model) {
  if (is.null(model)) NA_real_ else unname(summary(model)$adj.r.squared)
}

model_r_squared <- function(model) {
  if (is.null(model)) NA_real_ else unname(summary(model)$r.squared)
}

rmse <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  if (sum(valid) == 0) return(NA_real_)
  sqrt(mean((observed[valid] - predicted[valid])^2))
}

mae <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  if (sum(valid) == 0) return(NA_real_)
  mean(abs(observed[valid] - predicted[valid]))
}

make_temporal_folds <- function(year, k = cv_folds) {
  unique_years <- sort(unique(year))
  year_groups <- split(
    unique_years,
    cut(seq_along(unique_years), breaks = k, labels = FALSE)
  )
  fold_lookup <- tibble(
    analysis_year = unique_years,
    cv_fold = unlist(map2(year_groups, seq_along(year_groups), ~ rep(.y, length(.x))))
  )
  fold_lookup$cv_fold[match(year, fold_lookup$analysis_year)]
}

predict_temporal_cv <- function(data, model_type, index_name) {
  data <- data |>
    arrange(analysis_year) |>
    mutate(cv_fold = make_temporal_folds(analysis_year))

  predictions <- rep(NA_real_, nrow(data))
  null_predictions <- rep(NA_real_, nrow(data))

  for (fold in sort(unique(data$cv_fold))) {
    train <- data[data$cv_fold != fold, , drop = FALSE]
    test <- data[data$cv_fold == fold, , drop = FALSE]
    test_rows <- which(data$cv_fold == fold)

    null_predictions[test_rows] <- mean(
      train$yield_anomaly_percent,
      na.rm = TRUE
    )

    model <- tryCatch({
      if (model_type == "linear") {
        lm(
          yield_anomaly_percent ~ enso_value,
          data = transform(train, enso_value = train[[index_name]])
        )
      } else if (model_type == "quadratic") {
        lm(
          yield_anomaly_percent ~ enso_value + I(enso_value^2),
          data = transform(train, enso_value = train[[index_name]])
        )
      } else if (model_type == "phase") {
        lm(yield_anomaly_percent ~ oni_harvest_year_phase, data = train)
      } else {
        stop("Unknown model type: ", model_type)
      }
    }, error = function(e) NULL)

    if (!is.null(model)) {
      prediction_data <- test
      prediction_data$enso_value <- test[[index_name]]
      predictions[test_rows] <- tryCatch(
        as.numeric(predict(model, newdata = prediction_data)),
        error = function(e) rep(NA_real_, nrow(test))
      )
    }
  }

  tibble(
    cv_r_squared = safe_r_squared(
      data$yield_anomaly_percent,
      predictions,
      null_predictions
    ),
    cv_rmse_percent = rmse(data$yield_anomaly_percent, predictions),
    cv_mae_percent = mae(data$yield_anomaly_percent, predictions),
    cv_complete_predictions = sum(is.finite(predictions)),
    cv_total_observations = nrow(data)
  )
}

fit_annual_model <- function(data, model_type, index_name = primary_index) {
  model_data <- data |>
    transmute(
      analysis_year,
      yield_anomaly_percent,
      oni_harvest_year_phase,
      enso_value = .data[[index_name]]
    ) |>
    filter(
      is.finite(yield_anomaly_percent),
      if (model_type == "phase") !is.na(oni_harvest_year_phase) else is.finite(enso_value)
    )

  model <- if (model_type == "linear") {
    lm(yield_anomaly_percent ~ enso_value, data = model_data)
  } else if (model_type == "quadratic") {
    lm(yield_anomaly_percent ~ enso_value + I(enso_value^2), data = model_data)
  } else if (model_type == "phase") {
    lm(yield_anomaly_percent ~ oni_harvest_year_phase, data = model_data)
  } else {
    stop("Unknown model type: ", model_type)
  }

  fitted_values <- fitted(model)
  model_anova <- anova(model)
  total_ss <- sum((model_data$yield_anomaly_percent - mean(model_data$yield_anomaly_percent))^2)
  residual_ss <- sum(residuals(model)^2)

  model_p_value <- if (nrow(model_anova) >= 2) {
    model_anova$`Pr(>F)`[1]
  } else {
    NA_real_
  }

  coefficients <- coef(model)
  linear_slope <- if ("enso_value" %in% names(coefficients)) {
    unname(coefficients["enso_value"])
  } else {
    NA_real_
  }
  quadratic_term <- if ("I(enso_value^2)" %in% names(coefficients)) {
    unname(coefficients["I(enso_value^2)"])
  } else {
    NA_real_
  }

  cv <- predict_temporal_cv(data, model_type, index_name)

  tibble(
    model_type = model_type,
    enso_index_window = if (model_type == "phase") "ONI phase" else index_name,
    observations = nrow(model_data),
    r_squared = model_r_squared(model),
    adjusted_r_squared = adjusted_r_squared(model),
    variance_explained_percent = 100 * model_r_squared(model),
    adjusted_variance_explained_percent = 100 * adjusted_r_squared(model),
    residual_sum_squares = residual_ss,
    total_sum_squares = total_ss,
    model_p_value = model_p_value,
    linear_slope_percent_per_index_unit = linear_slope,
    quadratic_term = quadratic_term,
    cv_r_squared = cv$cv_r_squared,
    cv_variance_explained_percent = 100 * cv$cv_r_squared,
    cv_rmse_percent = cv$cv_rmse_percent,
    cv_mae_percent = cv$cv_mae_percent,
    cv_complete_predictions = cv$cv_complete_predictions
  )
}

permute_primary_annual <- function(data, repetitions, seed) {
  model_data <- data |>
    transmute(
      analysis_year,
      y = yield_anomaly_percent,
      x = .data[[primary_index]]
    ) |>
    filter(is.finite(y), is.finite(x)) |>
    arrange(analysis_year)

  observed_r2 <- summary(lm(y ~ x, data = model_data))$r.squared
  observed_slope <- coef(lm(y ~ x, data = model_data))["x"]

  set.seed(seed)
  permuted_r2 <- numeric(repetitions)
  permuted_slope <- numeric(repetitions)

  for (i in seq_len(repetitions)) {
    permuted_x <- sample(model_data$x, replace = FALSE)
    permuted_model <- lm(model_data$y ~ permuted_x)
    permuted_r2[i] <- summary(permuted_model)$r.squared
    permuted_slope[i] <- coef(permuted_model)[2]
  }

  # Circular shifts retain the temporal ordering and autocorrelation of ENSO.
  circular_r2 <- map_dbl(seq_len(nrow(model_data) - 1L), function(shift) {
    shifted_x <- c(
      tail(model_data$x, shift),
      head(model_data$x, -shift)
    )
    summary(lm(model_data$y ~ shifted_x))$r.squared
  })

  tibble(
    observed_r_squared = observed_r2,
    observed_slope = observed_slope,
    random_permutation_p_value =
      (1 + sum(permuted_r2 >= observed_r2)) / (repetitions + 1),
    random_permutation_slope_two_sided_p_value =
      (1 + sum(abs(permuted_slope) >= abs(observed_slope))) /
      (repetitions + 1),
    circular_shift_p_value =
      (1 + sum(circular_r2 >= observed_r2)) /
      (length(circular_r2) + 1),
    permutations = repetitions,
    circular_shifts = length(circular_r2)
  )
}

screen_annual_window <- function(data, index_name) {
  model_data <- data |>
    transmute(
      analysis_year,
      y = yield_anomaly_percent,
      x = .data[[index_name]]
    ) |>
    filter(is.finite(y), is.finite(x))

  if (nrow(model_data) < 20 || sd(model_data$x) == 0) {
    return(tibble(
      enso_index_window = index_name,
      observations = nrow(model_data),
      slope = NA_real_,
      r_squared = NA_real_,
      adjusted_r_squared = NA_real_,
      nominal_p_value = NA_real_,
      cv_r_squared = NA_real_
    ))
  }

  model <- lm(y ~ x, data = model_data)
  cv <- predict_temporal_cv(data, "linear", index_name)

  tibble(
    enso_index_window = index_name,
    observations = nrow(model_data),
    slope = unname(coef(model)["x"]),
    r_squared = summary(model)$r.squared,
    adjusted_r_squared = summary(model)$adj.r.squared,
    nominal_p_value = coef(summary(model))["x", "Pr(>|t|)"],
    cv_r_squared = cv$cv_r_squared
  )
}

max_window_permutation_p <- function(data, index_columns, repetitions, seed) {
  complete_columns <- index_columns[vapply(
    data[index_columns],
    function(x) sum(is.finite(x)) >= 20,
    logical(1)
  )]

  y <- data$yield_anomaly_percent
  observed_correlations <- vapply(complete_columns, function(index_name) {
    cor(y, data[[index_name]], use = "complete.obs")
  }, numeric(1))

  observed_max <- max(abs(observed_correlations), na.rm = TRUE)

  set.seed(seed)
  permutation_maxima <- replicate(repetitions, {
    permuted_y <- sample(y, replace = FALSE)
    max(abs(vapply(complete_columns, function(index_name) {
      cor(permuted_y, data[[index_name]], use = "complete.obs")
    }, numeric(1))), na.rm = TRUE)
  })

  tibble(
    observed_maximum_absolute_correlation = observed_max,
    familywise_permutation_p_value =
      (1 + sum(permutation_maxima >= observed_max)) / (repetitions + 1),
    screened_windows = length(complete_columns),
    permutations = repetitions
  )
}

weighted_fixed_effect_linear <- function(data, index_name = primary_index) {
  model_data <- data |>
    transmute(
      country,
      spatial_unit_id,
      analysis_year,
      y = yield_anomaly_percent,
      x = .data[[index_name]],
      weight = harvested_area_ha
    ) |>
    filter(
      is.finite(y), is.finite(x), is.finite(weight), weight > 0
    ) |>
    mutate(unit_key = interaction(country, spatial_unit_id, drop = TRUE)) |>
    group_by(unit_key) |>
    mutate(
      unit_weight_total = sum(weight),
      y_within = y - weighted.mean(y, weight),
      x_within = x - weighted.mean(x, weight)
    ) |>
    ungroup()

  denominator <- sum(model_data$weight * model_data$x_within^2)
  slope <- sum(
    model_data$weight * model_data$x_within * model_data$y_within
  ) / denominator

  residual_null <- model_data$y_within
  residual_full <- model_data$y_within - slope * model_data$x_within
  sse_null <- sum(model_data$weight * residual_null^2)
  sse_full <- sum(model_data$weight * residual_full^2)
  partial_r2 <- 1 - sse_full / sse_null

  list(
    data = model_data,
    summary = tibble(
      enso_index_window = index_name,
      observations = nrow(model_data),
      spatial_units = n_distinct(model_data$unit_key),
      harvest_years = n_distinct(model_data$analysis_year),
      area_weighted_within_unit_slope = slope,
      area_weighted_partial_r_squared = partial_r2,
      area_weighted_variance_explained_percent = 100 * partial_r2,
      null_sum_squares = sse_null,
      full_sum_squares = sse_full
    )
  )
}

permute_panel_years <- function(fixed_effect_result, repetitions, seed) {
  data <- fixed_effect_result$data
  observed <- fixed_effect_result$summary$area_weighted_partial_r_squared
  year_x <- data |>
    distinct(analysis_year, x) |>
    arrange(analysis_year)

  year_position <- match(data$analysis_year, year_x$analysis_year)
  unit_index <- as.integer(data$unit_key)
  unit_weight <- as.numeric(rowsum(
    data$weight,
    unit_index,
    reorder = FALSE
  ))

  set.seed(seed)
  permuted_r2 <- numeric(repetitions)

  for (i in seq_len(repetitions)) {
    permuted_year_x <- sample(year_x$x, replace = FALSE)
    x_permuted <- permuted_year_x[year_position]
    unit_weighted_x <- as.numeric(rowsum(
      data$weight * x_permuted,
      unit_index,
      reorder = FALSE
    )) / unit_weight
    x_permuted_within <- x_permuted - unit_weighted_x[unit_index]

    denominator <- sum(
      data$weight * x_permuted_within^2
    )
    slope <- sum(
      data$weight * x_permuted_within * data$y_within
    ) / denominator
    residual <- data$y_within - slope * x_permuted_within
    sse <- sum(data$weight * residual^2)
    permuted_r2[i] <- 1 - sse / fixed_effect_result$summary$null_sum_squares
  }

  tibble(
    year_cluster_permutation_p_value =
      (1 + sum(permuted_r2 >= observed)) / (repetitions + 1),
    permutations = repetitions
  )
}

fit_unit_model <- function(unit_data, index_name = primary_index) {
  model_data <- unit_data |>
    transmute(
      analysis_year,
      y = yield_anomaly_percent,
      x = .data[[index_name]]
    ) |>
    filter(is.finite(y), is.finite(x))

  if (nrow(model_data) < minimum_unit_years || sd(model_data$x) == 0) {
    return(tibble(
      observations = nrow(model_data),
      slope_percent_per_oni_c = NA_real_,
      r_squared = NA_real_,
      adjusted_r_squared = NA_real_,
      variance_explained_percent = NA_real_,
      p_value = NA_real_,
      cv_r_squared = NA_real_
    ))
  }

  model <- lm(y ~ x, data = model_data)
  cv_input <- unit_data |>
    rename(yield_anomaly_percent_original = yield_anomaly_percent) |>
    mutate(yield_anomaly_percent = yield_anomaly_percent_original)
  cv <- predict_temporal_cv(cv_input, "linear", index_name)

  tibble(
    observations = nrow(model_data),
    slope_percent_per_oni_c = unname(coef(model)["x"]),
    r_squared = summary(model)$r.squared,
    adjusted_r_squared = summary(model)$adj.r.squared,
    variance_explained_percent = 100 * summary(model)$r.squared,
    p_value = coef(summary(model))["x", "Pr(>|t|)"],
    cv_r_squared = cv$cv_r_squared
  )
}


# 5. READ INPUTS ---------------------------------------------------------------

message("Reading merged soybean and ENSO datasets...")

panel <- read_csv(
  panel_file,
  col_types = cols(
    spatial_unit_id = col_character(),
    analysis_year = col_integer(),
    harvested_area_ha = col_double(),
    yield_anomaly_percent = col_double(),
    .default = col_guess()
  ),
  na = c("", "NA", "N/A"),
  show_col_types = FALSE
)

country_year <- read_csv(
  country_year_file,
  col_types = cols(analysis_year = col_integer(), .default = col_guess()),
  na = c("", "NA", "N/A"),
  show_col_types = FALSE
)

combined_year <- read_csv(
  combined_year_file,
  col_types = cols(analysis_year = col_integer(), .default = col_guess()),
  na = c("", "NA", "N/A"),
  show_col_types = FALSE
)

combined_year$country <- "Brazil and Argentina"
annual <- bind_rows(country_year, combined_year) |>
  filter(
    analysis_year >= analysis_start_year,
    analysis_year <= analysis_end_year
  ) |>
  mutate(
    oni_harvest_year_phase = factor(
      oni_harvest_year_phase,
      levels = c("La Nina", "Neutral", "El Nino")
    )
  ) |>
  arrange(country, analysis_year)

index_columns <- c(
  "oni_son_previous_c", "oni_ond_previous_c", "oni_ndj_c",
  "oni_djf_c", "oni_jfm_c", "oni_fma_c", "oni_mam_c",
  "nino34_son_previous_c", "nino34_ond_previous_c", "nino34_ndj_c",
  "nino34_djf_c", "nino34_jfm_c", "nino34_fma_c", "nino34_mam_c",
  "nino34_growing_season_oct_apr_c",
  "nino34_reproductive_jan_mar_c",
  "meiv2_growing_season_oct_apr",
  "meiv2_reproductive_jan_mar"
)

common_required <- c(
  "country", "analysis_year", "yield_anomaly_percent",
  "oni_harvest_year_phase", index_columns
)

check_required_columns(annual, common_required, "Annual merged dataset")
check_required_columns(
  panel,
  c(
    "country", "spatial_unit_id", "spatial_unit_name",
    "parent_region_name", "analysis_year", "harvested_area_ha",
    "yield_anomaly_percent", "oni_harvest_year_phase", primary_index
  ),
  "Spatial panel"
)


# 6. INPUT QUALITY CONTROL -----------------------------------------------------

if (anyDuplicated(annual[c("country", "analysis_year")]) > 0) {
  stop("Duplicated country and year keys were found.", call. = FALSE)
}

if (anyDuplicated(panel[c("country", "spatial_unit_id", "analysis_year")]) > 0) {
  stop("Duplicated spatial unit and year keys were found.", call. = FALSE)
}

annual_year_counts <- annual |>
  count(country, name = "years")

if (any(annual_year_counts$years != 51L)) {
  stop("Each annual geography must contain exactly 51 harvest years.", call. = FALSE)
}

if (any(is.na(annual[[primary_index]]))) {
  stop("The primary ONI DJF index contains missing values.", call. = FALSE)
}


# 7. PRIMARY ANNUAL MODELS -----------------------------------------------------

message("Fitting primary annual models...")

annual_models <- annual |>
  group_by(country) |>
  group_modify(function(.x, .y) {
    bind_rows(
      fit_annual_model(.x, "linear", primary_index),
      fit_annual_model(.x, "quadratic", primary_index),
      fit_annual_model(.x, "phase", primary_index)
    )
  }) |>
  ungroup() |>
  arrange(country, model_type)

annual_permutation_tests <- annual |>
  group_by(country) |>
  group_modify(function(.x, .y) {
    permute_primary_annual(
      .x,
      repetitions = permutation_repetitions,
      seed = permutation_seed + match(.y$country, sort(unique(annual$country)))
    )
  }) |>
  ungroup()

# Test whether the quadratic term improves the linear specification.
annual_nonlinearity_tests <- annual |>
  group_by(country) |>
  group_modify(function(.x, .y) {
    model_data <- .x |>
      transmute(
        y = yield_anomaly_percent,
        x = .data[[primary_index]]
      ) |>
      filter(is.finite(y), is.finite(x))
    linear_model <- lm(y ~ x, data = model_data)
    quadratic_model <- lm(y ~ x + I(x^2), data = model_data)
    comparison <- anova(linear_model, quadratic_model)
    tibble(
      linear_r_squared = summary(linear_model)$r.squared,
      quadratic_r_squared = summary(quadratic_model)$r.squared,
      incremental_r_squared =
        summary(quadratic_model)$r.squared - summary(linear_model)$r.squared,
      quadratic_increment_p_value = comparison$`Pr(>F)`[2]
    )
  }) |>
  ungroup()


# 8. SEASONAL WINDOW SENSITIVITY ----------------------------------------------

message("Screening ENSO seasonal windows...")

window_sensitivity <- annual |>
  group_by(country) |>
  group_modify(function(country_data, group_key) {
    bind_rows(map(
      index_columns,
      ~ screen_annual_window(country_data, .x)
    )) |>
      mutate(
        p_value_fdr = p.adjust(nominal_p_value, method = "BH"),
        rank_by_r_squared = rank(-r_squared, ties.method = "min"),
        is_primary_window = enso_index_window == primary_index
      )
  }) |>
  ungroup() |>
  arrange(country, rank_by_r_squared, enso_index_window)

window_familywise_tests <- annual |>
  group_by(country) |>
  group_modify(function(.x, .y) {
    max_window_permutation_p(
      .x,
      index_columns,
      repetitions = permutation_repetitions,
      seed = permutation_seed + 100L +
        match(.y$country, sort(unique(annual$country)))
    )
  }) |>
  ungroup()


# 9. AREA WEIGHTED SPATIAL PANEL ----------------------------------------------

message("Fitting area weighted spatial fixed effect models...")

panel_results <- list()
panel_permutation_results <- list()

for (geography in c("Argentina", "Brazil", "Brazil and Argentina")) {
  panel_subset <- if (geography == "Brazil and Argentina") {
    panel
  } else {
    panel |> filter(country == geography)
  }

  fitted_panel <- weighted_fixed_effect_linear(panel_subset, primary_index)
  fitted_panel$summary$country <- geography
  panel_results[[geography]] <- fitted_panel$summary

  panel_test <- permute_panel_years(
    fitted_panel,
    repetitions = permutation_repetitions,
    seed = permutation_seed + 200L + match(
      geography,
      c("Argentina", "Brazil", "Brazil and Argentina")
    )
  )
  panel_test$country <- geography
  panel_permutation_results[[geography]] <- panel_test
}

panel_fixed_effects <- bind_rows(panel_results) |>
  relocate(country)

panel_permutation_tests <- bind_rows(panel_permutation_results) |>
  relocate(country)


# 10. SPATIAL UNIT MODELS ------------------------------------------------------

message("Fitting unit level models for 2,146 spatial units...")

unit_models <- panel |>
  group_by(
    country,
    spatial_unit_id,
    spatial_unit_name,
    parent_region_name
  ) |>
  group_modify(~ fit_unit_model(.x, primary_index)) |>
  ungroup() |>
  group_by(country) |>
  mutate(
    p_value_fdr = p.adjust(p_value, method = "BH"),
    significant_nominal_005 = is.finite(p_value) & p_value < 0.05,
    significant_fdr_005 = is.finite(p_value_fdr) & p_value_fdr < 0.05,
    response_direction = case_when(
      !is.finite(slope_percent_per_oni_c) ~ "Unavailable",
      slope_percent_per_oni_c > 0 ~ "Positive",
      slope_percent_per_oni_c < 0 ~ "Negative",
      TRUE ~ "Zero"
    )
  ) |>
  ungroup() |>
  arrange(country, spatial_unit_id)

unit_model_summary <- unit_models |>
  group_by(country) |>
  summarise(
    spatial_units = n(),
    successfully_fitted = sum(is.finite(r_squared)),
    median_r_squared = median(r_squared, na.rm = TRUE),
    mean_r_squared = mean(r_squared, na.rm = TRUE),
    median_cv_r_squared = median(cv_r_squared, na.rm = TRUE),
    positive_slopes = sum(response_direction == "Positive"),
    negative_slopes = sum(response_direction == "Negative"),
    nominally_significant_units = sum(significant_nominal_005),
    fdr_significant_units = sum(significant_fdr_005),
    percent_fdr_significant = 100 * mean(significant_fdr_005),
    .groups = "drop"
  )


# 11. CONSOLIDATED INTERPRETATION TABLE ---------------------------------------

primary_linear <- annual_models |>
  filter(model_type == "linear") |>
  select(
    country,
    annual_r_squared = r_squared,
    annual_adjusted_r_squared = adjusted_r_squared,
    annual_variance_explained_percent = variance_explained_percent,
    annual_cv_r_squared = cv_r_squared,
    annual_cv_variance_explained_percent = cv_variance_explained_percent,
    annual_slope_percent_per_oni_c = linear_slope_percent_per_index_unit,
    annual_nominal_p_value = model_p_value
  )

primary_results_summary <- primary_linear |>
  left_join(annual_permutation_tests, by = "country") |>
  left_join(
    panel_fixed_effects |>
      select(
        country,
        panel_partial_r_squared = area_weighted_partial_r_squared,
        panel_variance_explained_percent =
          area_weighted_variance_explained_percent,
        panel_slope_percent_per_oni_c =
          area_weighted_within_unit_slope
      ),
    by = "country"
  ) |>
  left_join(
    panel_permutation_tests |>
      select(country, panel_year_permutation_p_value =
        year_cluster_permutation_p_value),
    by = "country"
  ) |>
  left_join(annual_nonlinearity_tests, by = "country") |>
  mutate(
    evidence_interpretation = case_when(
      random_permutation_p_value < 0.05 & annual_cv_r_squared > 0 ~
        "Supported association with positive out of sample skill",
      random_permutation_p_value < 0.05 & annual_cv_r_squared <= 0 ~
        "In sample association without positive out of sample skill",
      TRUE ~ "Limited evidence for a general annual association"
    )
  ) |>
  arrange(country)


# 12. QUALITY CONTROL ----------------------------------------------------------

qc_model_summary <- tibble(
  check = c(
    "Annual geographies",
    "Annual observations",
    "Spatial panel observations",
    "Spatial units",
    "Primary annual models",
    "Seasonal window models",
    "Unit models",
    "Missing annual primary index",
    "Duplicated annual keys",
    "Duplicated panel keys"
  ),
  value = c(
    n_distinct(annual$country),
    nrow(annual),
    nrow(panel),
    n_distinct(interaction(panel$country, panel$spatial_unit_id)),
    nrow(annual_models),
    nrow(window_sensitivity),
    nrow(unit_models),
    sum(is.na(annual[[primary_index]])),
    sum(duplicated(annual[c("country", "analysis_year")])),
    sum(duplicated(panel[c("country", "spatial_unit_id", "analysis_year")]))
  )
)

if (nrow(annual_models) != 9L) {
  stop("Unexpected number of primary annual models.", call. = FALSE)
}

if (nrow(window_sensitivity) != 3L * length(index_columns)) {
  stop("Unexpected number of seasonal window sensitivity models.", call. = FALSE)
}

if (any(!is.finite(primary_results_summary$annual_r_squared))) {
  stop("At least one primary annual model failed.", call. = FALSE)
}

if (any(primary_results_summary$annual_r_squared < 0 |
        primary_results_summary$annual_r_squared > 1)) {
  stop("Annual R squared values outside the valid range were generated.", call. = FALSE)
}

if (any(panel_fixed_effects$area_weighted_partial_r_squared < 0 |
        panel_fixed_effects$area_weighted_partial_r_squared > 1)) {
  stop("Panel partial R squared values outside the valid range were generated.", call. = FALSE)
}


# 13. WRITE OUTPUTS ------------------------------------------------------------

message("Writing script 06 outputs...")

write_csv(
  primary_results_summary,
  file.path(output_dir, "ENSO_primary_variance_explained_summary.csv")
)

write_csv(
  annual_models,
  file.path(output_dir, "ENSO_annual_model_comparison.csv")
)

write_csv(
  annual_permutation_tests,
  file.path(output_dir, "ENSO_annual_permutation_tests.csv")
)

write_csv(
  annual_nonlinearity_tests,
  file.path(output_dir, "ENSO_annual_nonlinearity_tests.csv")
)

write_csv(
  window_sensitivity,
  file.path(output_dir, "ENSO_seasonal_window_sensitivity.csv")
)

write_csv(
  window_familywise_tests,
  file.path(output_dir, "ENSO_window_familywise_permutation_tests.csv")
)

write_csv(
  panel_fixed_effects,
  file.path(output_dir, "ENSO_spatial_panel_fixed_effects.csv")
)

write_csv(
  panel_permutation_tests,
  file.path(output_dir, "ENSO_spatial_panel_permutation_tests.csv")
)

write_csv(
  unit_models,
  file.path(output_dir, "ENSO_spatial_unit_models.csv")
)

write_csv(
  unit_model_summary,
  file.path(output_dir, "ENSO_spatial_unit_model_summary.csv")
)

write_csv(
  qc_model_summary,
  file.path(output_dir, "QC_ENSO_variance_models.csv")
)


# 14. COMPLETION MESSAGE -------------------------------------------------------

message("")
message("Script 06 completed successfully.")
message("Primary exposure: ONI DJF")
message("Output directory: ", output_dir)
message("")
message("Primary annual results:")
print(
  primary_results_summary |>
    select(
      country,
      annual_variance_explained_percent,
      annual_cv_variance_explained_percent,
      random_permutation_p_value,
      panel_variance_explained_percent,
      panel_year_permutation_p_value,
      evidence_interpretation
    )
)

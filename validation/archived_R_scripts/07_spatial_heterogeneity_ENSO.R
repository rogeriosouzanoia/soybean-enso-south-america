# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 07_spatial_heterogeneity_ENSO.R
# PURPOSE: Quantify municipal, departmental, state and provincial differences
#          in soybean productivity responses to ENSO
# VERSION: 1.1
#
# METHODS:
#   1. Local OLS response to the preregistered DJF ONI exposure
#   2. Region centered empirical Bayes shrinkage of local slopes
#   3. Temporal cross validation for every spatial unit
#   4. State and province annual models weighted by expected production
#   5. Country level heterogeneity tests and regional summaries
#
# INTERPRETATION:
#   Local slopes describe statistical sensitivity to ENSO. Shrinkage estimates
#   are the primary mapping values because isolated OLS estimates can be noisy.
#   The analysis quantifies association and does not establish causality.
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
primary_index <- "oni_djf_c"
minimum_unit_years <- 20L
minimum_region_years <- 20L
cv_folds <- 5L
permutation_repetitions <- 1999L
permutation_seed <- 20260904L


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
model_dir <- file.path(project_dir, "enso_variance_models")
output_dir <- file.path(project_dir, "enso_spatial_heterogeneity")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panel_file <- file.path(
  merged_dir,
  "soybean_ENSO_spatial_panel_1974_2024.csv"
)

script06_unit_file <- file.path(
  model_dir,
  "ENSO_spatial_unit_models.csv"
)

if (!file.exists(panel_file)) {
  stop("Spatial panel not found: ", panel_file, call. = FALSE)
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
  ifelse(
    is.finite(denominator) & denominator > 0,
    numerator / denominator,
    NA_real_
  )
}

make_temporal_folds <- function(year, k = cv_folds) {
  unique_years <- sort(unique(year))
  fold_number <- cut(
    seq_along(unique_years),
    breaks = k,
    labels = FALSE
  )
  fold_number[match(year, unique_years)]
}

temporal_cv_linear <- function(data, index_name = primary_index) {
  model_data <- data |>
    transmute(
      analysis_year,
      y = yield_anomaly_percent,
      x = .data[[index_name]]
    ) |>
    filter(is.finite(y), is.finite(x)) |>
    arrange(analysis_year) |>
    mutate(fold = make_temporal_folds(analysis_year))

  prediction <- rep(NA_real_, nrow(model_data))
  null_prediction <- rep(NA_real_, nrow(model_data))

  for (current_fold in sort(unique(model_data$fold))) {
    train <- model_data[model_data$fold != current_fold, , drop = FALSE]
    test <- model_data[model_data$fold == current_fold, , drop = FALSE]
    test_rows <- which(model_data$fold == current_fold)

    model <- lm(y ~ x, data = train)
    prediction[test_rows] <- predict(model, newdata = test)
    null_prediction[test_rows] <- mean(train$y)
  }

  null_sse <- sum((model_data$y - null_prediction)^2)
  model_sse <- sum((model_data$y - prediction)^2)

  tibble(
    cv_r_squared = if (null_sse > 0) 1 - model_sse / null_sse else NA_real_,
    cv_rmse_percent = sqrt(mean((model_data$y - prediction)^2)),
    cv_mae_percent = mean(abs(model_data$y - prediction))
  )
}

fit_local_linear <- function(data, index_name = primary_index) {
  model_data <- data |>
    transmute(
      analysis_year,
      y = yield_anomaly_percent,
      x = .data[[index_name]]
    ) |>
    filter(is.finite(y), is.finite(x))

  if (nrow(model_data) < minimum_unit_years || sd(model_data$x) == 0) {
    return(tibble(
      observations = nrow(model_data),
      slope_ols = NA_real_,
      slope_standard_error = NA_real_,
      slope_ci_lower_ols = NA_real_,
      slope_ci_upper_ols = NA_real_,
      r_squared = NA_real_,
      adjusted_r_squared = NA_real_,
      p_value = NA_real_,
      cv_r_squared = NA_real_,
      cv_rmse_percent = NA_real_,
      cv_mae_percent = NA_real_
    ))
  }

  model <- lm(y ~ x, data = model_data)
  coefficient_table <- coef(summary(model))
  slope <- unname(coef(model)["x"])
  slope_se <- coefficient_table["x", "Std. Error"]
  slope_df <- df.residual(model)
  critical_t <- qt(0.975, df = slope_df)
  cv <- temporal_cv_linear(data, index_name)

  tibble(
    observations = nrow(model_data),
    slope_ols = slope,
    slope_standard_error = slope_se,
    slope_ci_lower_ols = slope - critical_t * slope_se,
    slope_ci_upper_ols = slope + critical_t * slope_se,
    r_squared = summary(model)$r.squared,
    adjusted_r_squared = summary(model)$adj.r.squared,
    p_value = coefficient_table["x", "Pr(>|t|)"],
    cv_r_squared = cv$cv_r_squared,
    cv_rmse_percent = cv$cv_rmse_percent,
    cv_mae_percent = cv$cv_mae_percent
  )
}

estimate_empirical_bayes <- function(data) {
  valid <- data |>
    filter(
      is.finite(slope_ols),
      is.finite(slope_standard_error),
      slope_standard_error > 0
    )

  if (nrow(valid) < 2) {
    return(data |>
      mutate(
        country_mean_slope_random_effects = NA_real_,
        between_unit_variance_tau_squared = NA_real_,
        slope_eb = NA_real_,
        slope_eb_standard_error = NA_real_,
        slope_eb_ci_lower = NA_real_,
        slope_eb_ci_upper = NA_real_,
        posterior_probability_positive = NA_real_,
        posterior_probability_negative = NA_real_
      ))
  }

  sampling_variance <- valid$slope_standard_error^2
  fixed_weight <- 1 / sampling_variance
  fixed_mean <- sum(fixed_weight * valid$slope_ols) / sum(fixed_weight)
  q_statistic <- sum(fixed_weight * (valid$slope_ols - fixed_mean)^2)
  c_statistic <- sum(fixed_weight) -
    sum(fixed_weight^2) / sum(fixed_weight)
  tau_squared <- max(0, (q_statistic - (nrow(valid) - 1)) / c_statistic)

  random_weight <- 1 / (sampling_variance + tau_squared)
  random_mean <- sum(random_weight * valid$slope_ols) / sum(random_weight)
  random_mean_variance <- 1 / sum(random_weight)

  if (tau_squared > 0) {
    conditional_posterior_variance <- 1 / (
      1 / sampling_variance + 1 / tau_squared
    )
    prior_fraction <- conditional_posterior_variance / tau_squared
    posterior_mean <- conditional_posterior_variance * (
      valid$slope_ols / sampling_variance + random_mean / tau_squared
    )
    # Include uncertainty in the estimated regional prior mean.
    posterior_variance <- conditional_posterior_variance +
      prior_fraction^2 * random_mean_variance
  } else {
    posterior_variance <- rep(random_mean_variance, nrow(valid))
    posterior_mean <- rep(random_mean, nrow(valid))
  }

  valid <- valid |>
    mutate(
      country_mean_slope_random_effects = random_mean,
      between_unit_variance_tau_squared = tau_squared,
      slope_eb = posterior_mean,
      slope_eb_standard_error = sqrt(posterior_variance),
      slope_eb_ci_lower = slope_eb - 1.96 * slope_eb_standard_error,
      slope_eb_ci_upper = slope_eb + 1.96 * slope_eb_standard_error,
      posterior_probability_positive = if_else(
        slope_eb_standard_error > 0,
        pnorm(0, mean = slope_eb, sd = slope_eb_standard_error,
              lower.tail = FALSE),
        as.numeric(slope_eb > 0)
      ),
      posterior_probability_negative = 1 - posterior_probability_positive
    )

  data |>
    select(-any_of(c(
      "country_mean_slope_random_effects",
      "between_unit_variance_tau_squared",
      "slope_eb", "slope_eb_standard_error",
      "slope_eb_ci_lower", "slope_eb_ci_upper",
      "posterior_probability_positive",
      "posterior_probability_negative"
    ))) |>
    left_join(
      valid |>
        select(
          spatial_unit_id,
          country_mean_slope_random_effects,
          between_unit_variance_tau_squared,
          slope_eb,
          slope_eb_standard_error,
          slope_eb_ci_lower,
          slope_eb_ci_upper,
          posterior_probability_positive,
          posterior_probability_negative
        ),
      by = "spatial_unit_id"
    )
}

estimate_region_aware_empirical_bayes <- function(data) {
  # At least five units are required to estimate a regional prior. For smaller
  # groups, the local estimate is retained and its sampling uncertainty is not
  # artificially reduced.
  valid_units <- sum(
    is.finite(data$slope_ols) &
      is.finite(data$slope_standard_error) &
      data$slope_standard_error > 0
  )

  if (valid_units >= 5) {
    result <- estimate_empirical_bayes(data) |>
      rename(
        shrinkage_prior_mean_slope = country_mean_slope_random_effects,
        within_region_tau_squared = between_unit_variance_tau_squared
      ) |>
      mutate(
        shrinkage_reference_level = "State or province",
        regional_units_used_for_prior = valid_units
      )
    return(result)
  }

  data |>
    mutate(
      shrinkage_prior_mean_slope = NA_real_,
      within_region_tau_squared = NA_real_,
      slope_eb = slope_ols,
      slope_eb_standard_error = slope_standard_error,
      slope_eb_ci_lower = slope_ci_lower_ols,
      slope_eb_ci_upper = slope_ci_upper_ols,
      posterior_probability_positive = if_else(
        is.finite(slope_ols) & is.finite(slope_standard_error) &
          slope_standard_error > 0,
        pnorm(
          0,
          mean = slope_ols,
          sd = slope_standard_error,
          lower.tail = FALSE
        ),
        NA_real_
      ),
      posterior_probability_negative = 1 - posterior_probability_positive,
      shrinkage_reference_level = "No pooling, fewer than five regional units",
      regional_units_used_for_prior = valid_units
    )
}

heterogeneity_statistics <- function(data) {
  valid <- data |>
    filter(
      is.finite(slope_ols),
      is.finite(slope_standard_error),
      slope_standard_error > 0
    )

  sampling_variance <- valid$slope_standard_error^2
  weight <- 1 / sampling_variance
  fixed_mean <- sum(weight * valid$slope_ols) / sum(weight)
  q_statistic <- sum(weight * (valid$slope_ols - fixed_mean)^2)
  degrees_freedom <- nrow(valid) - 1
  c_statistic <- sum(weight) - sum(weight^2) / sum(weight)
  tau_squared <- max(0, (q_statistic - degrees_freedom) / c_statistic)
  i_squared <- max(0, 100 * (q_statistic - degrees_freedom) / q_statistic)
  random_weight <- 1 / (sampling_variance + tau_squared)
  random_mean <- sum(random_weight * valid$slope_ols) / sum(random_weight)
  random_mean_se <- sqrt(1 / sum(random_weight))

  tibble(
    spatial_units = nrow(valid),
    fixed_effect_mean_slope = fixed_mean,
    random_effect_mean_slope = random_mean,
    random_effect_mean_ci_lower = random_mean - 1.96 * random_mean_se,
    random_effect_mean_ci_upper = random_mean + 1.96 * random_mean_se,
    q_statistic = q_statistic,
    q_degrees_freedom = degrees_freedom,
    heterogeneity_p_value = pchisq(
      q_statistic,
      df = degrees_freedom,
      lower.tail = FALSE
    ),
    i_squared_percent = i_squared,
    tau_squared = tau_squared
  )
}

fit_region_model <- function(data) {
  model_data <- data |>
    filter(
      is.finite(yield_anomaly_percent),
      is.finite(.data[[primary_index]])
    )

  if (nrow(model_data) < minimum_region_years) {
    return(tibble(
      harvest_years = nrow(model_data),
      slope_percent_per_oni_c = NA_real_,
      slope_standard_error = NA_real_,
      r_squared = NA_real_,
      adjusted_r_squared = NA_real_,
      p_value = NA_real_,
      cv_r_squared = NA_real_,
      phase_r_squared = NA_real_,
      phase_p_value = NA_real_
    ))
  }

  linear_model <- lm(
    yield_anomaly_percent ~ oni_djf_c,
    data = model_data
  )
  phase_model <- lm(
    yield_anomaly_percent ~ oni_harvest_year_phase,
    data = model_data
  )
  coefficient_table <- coef(summary(linear_model))
  phase_anova <- anova(phase_model)
  cv <- temporal_cv_linear(model_data, primary_index)

  tibble(
    harvest_years = nrow(model_data),
    slope_percent_per_oni_c = unname(coef(linear_model)["oni_djf_c"]),
    slope_standard_error = coefficient_table["oni_djf_c", "Std. Error"],
    r_squared = summary(linear_model)$r.squared,
    adjusted_r_squared = summary(linear_model)$adj.r.squared,
    p_value = coefficient_table["oni_djf_c", "Pr(>|t|)"],
    cv_r_squared = cv$cv_r_squared,
    phase_r_squared = summary(phase_model)$r.squared,
    phase_p_value = phase_anova$`Pr(>F)`[1]
  )
}

permutation_p_value_region <- function(data, observed_slope, repetitions, seed) {
  model_data <- data |>
    transmute(
      y = yield_anomaly_percent,
      x = .data[[primary_index]]
    ) |>
    filter(is.finite(y), is.finite(x))

  set.seed(seed)
  permuted_slopes <- replicate(repetitions, {
    permuted_x <- sample(model_data$x, replace = FALSE)
    unname(coef(lm(model_data$y ~ permuted_x))[2])
  })

  (1 + sum(abs(permuted_slopes) >= abs(observed_slope))) /
    (repetitions + 1)
}


# 5. READ AND VALIDATE THE PANEL ----------------------------------------------

message("Reading the spatial panel...")

panel <- read_csv(
  panel_file,
  col_types = cols(
    spatial_unit_id = col_character(),
    parent_region_code = col_character(),
    analysis_year = col_integer(),
    harvested_area_ha = col_double(),
    production_t = col_double(),
    expected_production_t = col_double(),
    yield_anomaly_percent = col_double(),
    oni_djf_c = col_double(),
    .default = col_guess()
  ),
  na = c("", "NA", "N/A"),
  show_col_types = FALSE
) |>
  filter(
    analysis_year >= analysis_start_year,
    analysis_year <= analysis_end_year
  ) |>
  mutate(
    oni_harvest_year_phase = factor(
      oni_harvest_year_phase,
      levels = c("La Nina", "Neutral", "El Nino")
    )
  )

check_required_columns(
  panel,
  c(
    "country", "spatial_unit_id", "spatial_unit_name",
    "parent_region_code", "parent_region_name", "spatial_level",
    "analysis_year", "harvested_area_ha", "production_t",
    "expected_production_t", "yield_anomaly_percent",
    "oni_djf_c", "oni_harvest_year_phase"
  ),
  "Spatial panel"
)

if (anyDuplicated(panel[c(
  "country", "spatial_unit_id", "analysis_year"
)]) > 0) {
  stop("Duplicated spatial unit and year keys were found.", call. = FALSE)
}

if (any(is.na(panel[[primary_index]]))) {
  stop("Missing values were found in the primary ONI exposure.", call. = FALSE)
}


# 6. MUNICIPAL AND DEPARTMENTAL MODELS ----------------------------------------

message("Estimating local slopes and temporal validation...")

local_models_ols <- panel |>
  group_by(
    country,
    spatial_unit_id,
    spatial_unit_name,
    parent_region_code,
    parent_region_name,
    spatial_level
  ) |>
  group_modify(~ fit_local_linear(.x, primary_index)) |>
  ungroup()

local_models <- local_models_ols |>
  group_by(
    country,
    parent_region_code,
    parent_region_name
  ) |>
  group_modify(~ estimate_region_aware_empirical_bayes(.x)) |>
  ungroup() |>
  group_by(country) |>
  mutate(
    p_value_fdr = p.adjust(p_value, method = "BH"),
    shrinkage_fraction = safe_ratio(slope_ols - slope_eb, slope_ols),
    response_direction_eb = case_when(
      !is.finite(slope_eb) ~ "Unavailable",
      slope_eb > 0 ~ "Positive",
      slope_eb < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    response_evidence_eb = case_when(
      posterior_probability_positive >= 0.95 ~ "Positive, high probability",
      posterior_probability_negative >= 0.95 ~ "Negative, high probability",
      TRUE ~ "Uncertain"
    ),
    positive_out_of_sample_skill = is.finite(cv_r_squared) & cv_r_squared > 0,
    robust_local_pattern = case_when(
      response_evidence_eb == "Positive, high probability" &
        positive_out_of_sample_skill ~ "Robust positive",
      response_evidence_eb == "Negative, high probability" &
        positive_out_of_sample_skill ~ "Robust negative",
      response_evidence_eb != "Uncertain" ~
        "High posterior probability, limited predictive skill",
      TRUE ~ "Uncertain"
    )
  ) |>
  ungroup() |>
  arrange(country, parent_region_name, spatial_unit_name)


# 7. COUNTRY HETEROGENEITY -----------------------------------------------------

country_heterogeneity <- local_models |>
  group_by(country) |>
  group_modify(~ heterogeneity_statistics(.x)) |>
  ungroup()


# 8. STATE AND PROVINCE ANNUAL SERIES -----------------------------------------

# Expected production is the preferred weight because it measures how much
# production would have occurred without the annual productivity anomaly.
region_year <- panel |>
  group_by(
    country,
    parent_region_code,
    parent_region_name,
    analysis_year,
    oni_harvest_year_phase
  ) |>
  summarise(
    spatial_units = n_distinct(spatial_unit_id),
    harvested_area_ha = sum(harvested_area_ha, na.rm = TRUE),
    observed_production_t = sum(production_t, na.rm = TRUE),
    expected_production_t = sum(expected_production_t, na.rm = TRUE),
    yield_anomaly_percent = 100 * safe_ratio(
      observed_production_t - expected_production_t,
      expected_production_t
    ),
    oni_djf_c = first(oni_djf_c),
    .groups = "drop"
  ) |>
  filter(
    is.finite(yield_anomaly_percent),
    expected_production_t > 0
  ) |>
  arrange(country, parent_region_name, analysis_year)


# 9. STATE AND PROVINCE MODELS -------------------------------------------------

message("Estimating state and province responses...")

region_models <- region_year |>
  group_by(country, parent_region_code, parent_region_name) |>
  group_modify(function(region_data, region_key) {
    fitted <- fit_region_model(region_data)
    region_seed <- permutation_seed + sum(utf8ToInt(
      paste(region_key$country, region_key$parent_region_name)
    ))
    fitted$permutation_p_value <- if (
      is.finite(fitted$slope_percent_per_oni_c)
    ) {
      permutation_p_value_region(
        region_data,
        fitted$slope_percent_per_oni_c,
        permutation_repetitions,
        region_seed
      )
    } else {
      NA_real_
    }
    fitted
  }) |>
  ungroup() |>
  group_by(country) |>
  mutate(
    p_value_fdr = p.adjust(p_value, method = "BH"),
    permutation_p_value_fdr = p.adjust(
      permutation_p_value,
      method = "BH"
    ),
    response_direction = case_when(
      !is.finite(slope_percent_per_oni_c) ~ "Unavailable",
      slope_percent_per_oni_c > 0 ~ "Positive",
      slope_percent_per_oni_c < 0 ~ "Negative",
      TRUE ~ "Zero"
    )
  ) |>
  ungroup() |>
  arrange(country, desc(r_squared))


# 10. WITHIN REGION MUNICIPAL PATTERNS ----------------------------------------

within_region_summary <- local_models |>
  group_by(
    country,
    parent_region_code,
    parent_region_name
  ) |>
  summarise(
    spatial_units = n(),
    median_years = median(observations),
    mean_ols_slope = mean(slope_ols, na.rm = TRUE),
    median_ols_slope = median(slope_ols, na.rm = TRUE),
    mean_eb_slope = mean(slope_eb, na.rm = TRUE),
    median_eb_slope = median(slope_eb, na.rm = TRUE),
    minimum_eb_slope = min(slope_eb, na.rm = TRUE),
    maximum_eb_slope = max(slope_eb, na.rm = TRUE),
    positive_units = sum(slope_eb > 0, na.rm = TRUE),
    negative_units = sum(slope_eb < 0, na.rm = TRUE),
    percent_positive_units = 100 * mean(slope_eb > 0, na.rm = TRUE),
    percent_negative_units = 100 * mean(slope_eb < 0, na.rm = TRUE),
    high_probability_positive_units = sum(
      posterior_probability_positive >= 0.95,
      na.rm = TRUE
    ),
    high_probability_negative_units = sum(
      posterior_probability_negative >= 0.95,
      na.rm = TRUE
    ),
    positive_cv_units = sum(cv_r_squared > 0, na.rm = TRUE),
    percent_positive_cv_units = 100 * mean(cv_r_squared > 0, na.rm = TRUE),
    robust_positive_units = sum(robust_local_pattern == "Robust positive"),
    robust_negative_units = sum(robust_local_pattern == "Robust negative"),
    .groups = "drop"
  ) |>
  left_join(
    region_models |>
      select(
        country,
        parent_region_code,
        parent_region_name,
        regional_slope_percent_per_oni_c = slope_percent_per_oni_c,
        regional_r_squared = r_squared,
        regional_cv_r_squared = cv_r_squared,
        regional_permutation_p_value = permutation_p_value,
        regional_permutation_p_value_fdr = permutation_p_value_fdr
      ),
    by = c("country", "parent_region_code", "parent_region_name")
  ) |>
  arrange(country, desc(abs(regional_slope_percent_per_oni_c)))


# 11. COUNTRY SUMMARY ----------------------------------------------------------

country_spatial_summary <- local_models |>
  group_by(country) |>
  summarise(
    spatial_units = n(),
    positive_eb_slopes = sum(slope_eb > 0, na.rm = TRUE),
    negative_eb_slopes = sum(slope_eb < 0, na.rm = TRUE),
    percent_positive_eb_slopes = 100 * mean(slope_eb > 0, na.rm = TRUE),
    percent_negative_eb_slopes = 100 * mean(slope_eb < 0, na.rm = TRUE),
    high_probability_positive = sum(
      posterior_probability_positive >= 0.95,
      na.rm = TRUE
    ),
    high_probability_negative = sum(
      posterior_probability_negative >= 0.95,
      na.rm = TRUE
    ),
    positive_out_of_sample_skill = sum(cv_r_squared > 0, na.rm = TRUE),
    percent_positive_out_of_sample_skill = 100 * mean(
      cv_r_squared > 0,
      na.rm = TRUE
    ),
    robust_positive = sum(robust_local_pattern == "Robust positive"),
    robust_negative = sum(robust_local_pattern == "Robust negative"),
    median_eb_slope = median(slope_eb, na.rm = TRUE),
    eb_slope_10th_percentile = quantile(slope_eb, 0.10, na.rm = TRUE),
    eb_slope_90th_percentile = quantile(slope_eb, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(country_heterogeneity, by = "country")


# 12. CONSISTENCY WITH SCRIPT 06 ----------------------------------------------

script06_comparison <- tibble()

if (file.exists(script06_unit_file)) {
  script06_units <- read_csv(
    script06_unit_file,
    col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
    show_col_types = FALSE
  )

  script06_comparison <- local_models |>
    select(
      country,
      spatial_unit_id,
      slope_script07 = slope_ols,
      r_squared_script07 = r_squared,
      cv_r_squared_script07 = cv_r_squared
    ) |>
    left_join(
      script06_units |>
        select(
          country,
          spatial_unit_id,
          slope_script06 = slope_percent_per_oni_c,
          r_squared_script06 = r_squared,
          cv_r_squared_script06 = cv_r_squared
        ),
      by = c("country", "spatial_unit_id")
    ) |>
    summarise(
      compared_units = n(),
      maximum_absolute_slope_difference = max(
        abs(slope_script07 - slope_script06),
        na.rm = TRUE
      ),
      maximum_absolute_r_squared_difference = max(
        abs(r_squared_script07 - r_squared_script06),
        na.rm = TRUE
      ),
      maximum_absolute_cv_r_squared_difference = max(
        abs(cv_r_squared_script07 - cv_r_squared_script06),
        na.rm = TRUE
      )
    )
}


# 13. QUALITY CONTROL ----------------------------------------------------------

qc_spatial_heterogeneity <- tibble(
  check = c(
    "Panel observations",
    "Spatial units",
    "Countries",
    "Local models",
    "Local model failures",
    "Regions",
    "Region models",
    "Missing primary ONI",
    "Duplicated panel keys",
    "Empirical Bayes estimates"
  ),
  value = c(
    nrow(panel),
    n_distinct(interaction(panel$country, panel$spatial_unit_id)),
    n_distinct(panel$country),
    nrow(local_models),
    sum(!is.finite(local_models$slope_ols)),
    n_distinct(interaction(panel$country, panel$parent_region_name)),
    nrow(region_models),
    sum(is.na(panel[[primary_index]])),
    sum(duplicated(panel[c(
      "country", "spatial_unit_id", "analysis_year"
    )])),
    sum(is.finite(local_models$slope_eb))
  )
)

if (nrow(local_models) != 2146L) {
  stop("The expected 2,146 spatial unit models were not generated.", call. = FALSE)
}

if (any(!is.finite(local_models$slope_eb))) {
  stop("At least one empirical Bayes slope is unavailable.", call. = FALSE)
}

if (any(local_models$posterior_probability_positive < 0 |
        local_models$posterior_probability_positive > 1)) {
  stop("Invalid posterior probabilities were generated.", call. = FALSE)
}

if (any(country_heterogeneity$i_squared_percent < 0 |
        country_heterogeneity$i_squared_percent > 100)) {
  stop("Invalid heterogeneity statistics were generated.", call. = FALSE)
}


# 14. WRITE OUTPUTS ------------------------------------------------------------

message("Writing spatial heterogeneity outputs...")

write_csv(
  local_models,
  file.path(output_dir, "ENSO_local_responses_empirical_bayes.csv")
)

write_csv(
  region_year,
  file.path(output_dir, "ENSO_state_province_annual_series.csv")
)

write_csv(
  region_models,
  file.path(output_dir, "ENSO_state_province_models.csv")
)

write_csv(
  within_region_summary,
  file.path(output_dir, "ENSO_within_state_province_summary.csv")
)

write_csv(
  country_heterogeneity,
  file.path(output_dir, "ENSO_country_spatial_heterogeneity_tests.csv")
)

write_csv(
  country_spatial_summary,
  file.path(output_dir, "ENSO_country_spatial_response_summary.csv")
)

write_csv(
  script06_comparison,
  file.path(output_dir, "QC_script06_script07_consistency.csv")
)

write_csv(
  qc_spatial_heterogeneity,
  file.path(output_dir, "QC_ENSO_spatial_heterogeneity.csv")
)


# 15. COMPLETION MESSAGE -------------------------------------------------------

message("")
message("Script 07 completed successfully.")
message("Output directory: ", output_dir)
message("")
message("Country spatial response summary:")
print(country_spatial_summary)
message("")
message("The primary mapping variable is slope_eb.")
message("Use robust_local_pattern to distinguish predictive local responses.")

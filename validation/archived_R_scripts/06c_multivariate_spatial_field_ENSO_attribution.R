# ==============================================================================
# PROJECT: ENSO contribution to soybean productivity in South America
# SCRIPT: 06c_multivariate_spatial_field_ENSO_attribution.R
# VERSION: 1.1
# PURPOSE: Independent spatial field attribution using nested, event blocked,
#          multivariate ridge regression at state and province level.
# NOTE: Estimates predictive association. It does not establish causality.
# ==============================================================================

# 1. SETTINGS ------------------------------------------------------------------
default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) project_dir_override else default_project_dir

analysis_start_year <- 1974L
analysis_end_year <- 2024L
minimum_local_unit_years <- 20L
minimum_region_years <- 25L
minimum_training_years_per_region <- 10L
outer_folds <- 5L
inner_folds <- 4L
cv_gap_years <- 1L
reference_area_start_year <- 2015L
reference_area_end_year <- 2024L
lambda_grid <- c(0, 10^seq(-4, 4, length.out = 41))
candidate_families <- c("symmetric", "asymmetric")
bootstrap_repetitions <- 1999L
bootstrap_block_length <- 3L
bootstrap_seed <- 20260906L
run_circular_shift_test <- TRUE
circular_shift_repetitions <- 49L
circular_shift_seed <- 20260907L

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
panel_file <- file.path(
  project_dir, "data_merged", "soybean_ENSO_spatial_panel_1974_2024.csv"
)
output_dir <- file.path(project_dir, "enso_multivariate_field_attribution")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(panel_file)) stop("Required file not found: ", panel_file, call. = FALSE)

# 4. GENERAL HELPERS -----------------------------------------------------------
check_columns <- function(x, required) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

normalize_spatial_id <- function(country, spatial_id) {
  digits <- gsub("[^0-9]", "", as.character(spatial_id))
  result <- digits
  brazil <- country == "Brazil" & !is.na(digits) & nzchar(digits)
  argentina <- country == "Argentina" & !is.na(digits) & nzchar(digits)
  result[brazil] <- sprintf("%07d", as.integer(digits[brazil]))
  result[argentina] <- sprintf("%05d", as.integer(digits[argentina]))
  result
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

# Entire consecutive ENSO episodes remain in one fold. Balancing is conducted
# within El Nino, La Nina and neutral episodes.
create_event_folds <- function(year_table, k, fold_name = "cv_fold") {
  tagged <- year_table |>
    arrange(analysis_year) |>
    mutate(
      phase_simple = as.character(official_enso_phase),
      episode_id = cumsum(
        phase_simple != lag(phase_simple, default = first(phase_simple))
      ) + 1L
    )
  episodes <- tagged |>
    count(episode_id, phase_simple, name = "episode_years") |>
    group_by(phase_simple) |>
    group_modify(~ {
      z <- arrange(.x, desc(episode_years), episode_id)
      load <- rep(0L, k)
      z$.fold <- NA_integer_
      for (i in seq_len(nrow(z))) {
        f <- which.min(load)
        z$.fold[i] <- f
        load[f] <- load[f] + z$episode_years[i]
      }
      z
    }) |>
    ungroup()
  result <- tagged |>
    left_join(select(episodes, episode_id, .fold), by = "episode_id") |>
    select(analysis_year, phase_simple, episode_id, .fold)
  names(result)[names(result) == ".fold"] <- fold_name
  result
}

exclude_gap_years <- function(test_years, gap) {
  unique(unlist(lapply(test_years, function(y) seq(y - gap, y + gap))))
}

# 5. RIDGE HELPERS -------------------------------------------------------------
# The intercept is never penalized. ENSO intensity is scaled using training
# data only. Zero ONI remains the physical boundary between warm and cold.
make_design <- function(oni, family, training_scale) {
  if (!is.finite(training_scale) || training_scale <= 0) training_scale <- 1
  if (family == "symmetric") {
    cbind(intercept = 1, oni = oni / training_scale)
  } else if (family == "asymmetric") {
    cbind(
      intercept = 1,
      el_nino = pmax(oni, 0) / training_scale,
      la_nina = pmax(-oni, 0) / training_scale
    )
  } else {
    stop("Unknown model family: ", family, call. = FALSE)
  }
}

ridge_fit_one_region <- function(train, family, lambda, training_scale) {
  ok <- is.finite(train$yield_anomaly_percent) & is.finite(train$oni_model)
  z <- train[ok, , drop = FALSE]
  if (nrow(z) < minimum_training_years_per_region) return(NULL)
  X <- make_design(z$oni_model, family, training_scale)
  y <- z$yield_anomaly_percent
  penalty <- diag(ncol(X))
  penalty[1, 1] <- 0
  A <- crossprod(X) + lambda * penalty
  b <- as.numeric(crossprod(X, y))
  # SVD gives the Moore Penrose solution when a shifted training sample lacks
  # one ENSO phase and the asymmetric design becomes rank deficient. For full
  # rank systems it is numerically equivalent to the ordinary ridge solution.
  decomposition <- svd(A)
  tolerance <- max(dim(A)) * max(decomposition$d) * sqrt(.Machine$double.eps)
  inverse_d <- rep(0, length(decomposition$d))
  retained <- decomposition$d > tolerance
  inverse_d[retained] <- 1 / decomposition$d[retained]
  beta <- decomposition$v %*%
    (inverse_d * as.numeric(crossprod(decomposition$u, b)))
  list(beta = as.numeric(beta), null_mean = mean(y), n_training_years = nrow(z))
}

predict_one_region <- function(fit, oni, family, training_scale) {
  if (is.null(fit)) return(rep(NA_real_, length(oni)))
  as.numeric(make_design(oni, family, training_scale) %*% fit$beta)
}

predict_neutral_one_region <- function(fit, n, family, training_scale) {
  predict_one_region(fit, rep(0, n), family, training_scale)
}

fit_field <- function(train, family, lambda) {
  training_scale <- sd(unique(train$oni_model[is.finite(train$oni_model)]))
  if (!is.finite(training_scale) || training_scale <= 0) training_scale <- 1
  fits <- train |>
    group_split(region_key) |>
    set_names(vapply(train |> group_split(region_key),
                     function(x) first(x$region_key), character(1))) |>
    lapply(ridge_fit_one_region, family = family, lambda = lambda,
           training_scale = training_scale)
  list(fits = fits, family = family, lambda = lambda, training_scale = training_scale)
}

predict_field <- function(object, newdata) {
  rows <- newdata |>
    group_split(region_key) |>
    lapply(function(z) {
      key <- first(z$region_key)
      fit <- object$fits[[key]]
      z$model_prediction_percent <- predict_one_region(
        fit, z$oni_model, object$family, object$training_scale
      )
      z$neutral_prediction_percent <- predict_neutral_one_region(
        fit, nrow(z), object$family, object$training_scale
      )
      z$null_prediction_percent <- if (is.null(fit)) NA_real_ else fit$null_mean
      z
    })
  bind_rows(rows)
}

score_predictions <- function(predictions) {
  predictive_r2(
    predictions$yield_anomaly_percent,
    predictions$model_prediction_percent,
    predictions$null_prediction_percent,
    predictions$score_weight
  )
}

# Inner validation chooses both functional family and lambda without accessing
# the outer test years.
tune_field_model <- function(outer_train, inner_k = inner_folds) {
  inner_years <- outer_train |>
    distinct(analysis_year, official_enso_phase) |>
    create_event_folds(inner_k, fold_name = "inner_fold")
  d <- outer_train |> left_join(select(inner_years, analysis_year, inner_fold),
                                by = "analysis_year")
  candidates <- crossing(family = candidate_families, lambda = lambda_grid)
  candidate_scores <- vector("list", nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    family_i <- candidates$family[i]
    lambda_i <- candidates$lambda[i]
    fold_sse_model <- 0
    fold_sse_null <- 0
    scored_folds <- 0L
    for (fold in sort(unique(d$inner_fold))) {
      validation_years <- unique(d$analysis_year[d$inner_fold == fold])
      excluded <- exclude_gap_years(validation_years, cv_gap_years)
      train_i <- filter(d, !analysis_year %in% excluded)
      valid_i <- filter(d, inner_fold == fold)
      object <- fit_field(train_i, family_i, lambda_i)
      pred_i <- predict_field(object, valid_i)
      sse_m <- weighted_sse(
        pred_i$yield_anomaly_percent, pred_i$model_prediction_percent,
        pred_i$score_weight
      )
      sse_n <- weighted_sse(
        pred_i$yield_anomaly_percent, pred_i$null_prediction_percent,
        pred_i$score_weight
      )
      if (is.finite(sse_m) && is.finite(sse_n) && sse_n > 0) {
        fold_sse_model <- fold_sse_model + sse_m
        fold_sse_null <- fold_sse_null + sse_n
        scored_folds <- scored_folds + 1L
      }
    }
    candidate_scores[[i]] <- tibble(
      family = family_i,
      lambda = lambda_i,
      inner_cv_r_squared = if_else(
        fold_sse_null > 0, 1 - fold_sse_model / fold_sse_null, NA_real_
      ),
      scored_inner_folds = scored_folds
    )
  }
  scores <- bind_rows(candidate_scores) |>
    filter(is.finite(inner_cv_r_squared), scored_inner_folds == inner_k) |>
    arrange(desc(inner_cv_r_squared), desc(lambda), family)
  if (!nrow(scores)) stop("No ridge candidate completed inner validation.", call. = FALSE)
  # A simpler symmetric response is preferred when it is within 0.25 percentage
  # points of the best inner validation score. Within family, stronger shrinkage
  # is preferred when tied.
  best_score <- max(scores$inner_cv_r_squared)
  eligible <- scores |>
    filter(inner_cv_r_squared >= best_score - 0.0025) |>
    mutate(family_rank = if_else(family == "symmetric", 1L, 2L)) |>
    arrange(family_rank, desc(lambda), desc(inner_cv_r_squared))
  list(selected = slice(eligible, 1), scores = scores)
}

run_nested_country <- function(country_data, save_tuning = TRUE) {
  predictions <- list()
  tuning <- list()
  z <- 0L
  for (fold in sort(unique(country_data$outer_fold))) {
    test_years <- unique(country_data$analysis_year[country_data$outer_fold == fold])
    excluded <- exclude_gap_years(test_years, cv_gap_years)
    train <- filter(country_data, !analysis_year %in% excluded)
    test <- filter(country_data, outer_fold == fold)
    tuned <- tune_field_model(train)
    family_selected <- tuned$selected$family[[1]]
    lambda_selected <- tuned$selected$lambda[[1]]
    object <- fit_field(train, family_selected, lambda_selected)
    pred <- predict_field(object, test) |>
      mutate(
        outer_fold = fold,
        selected_family = family_selected,
        selected_lambda = lambda_selected,
        enso_component_percent = model_prediction_percent - neutral_prediction_percent,
        squared_error_model = (yield_anomaly_percent - model_prediction_percent)^2,
        squared_error_null = (yield_anomaly_percent - null_prediction_percent)^2
      )
    z <- z + 1L
    predictions[[z]] <- pred
    if (save_tuning) {
      tuning[[z]] <- tuned$scores |>
        mutate(
          country = first(country_data$country),
          outer_fold = fold,
          selected = family == family_selected & lambda == lambda_selected,
          outer_test_years = n_distinct(test$analysis_year),
          outer_training_years = n_distinct(train$analysis_year)
        ) |>
        relocate(country, outer_fold, selected, outer_test_years,
                 outer_training_years)
    }
  }
  list(predictions = bind_rows(predictions), tuning = bind_rows(tuning))
}

# 6. DATA PREPARATION ----------------------------------------------------------
panel <- read_csv(panel_file, show_col_types = FALSE)
check_columns(panel, c(
  "country", "spatial_unit_id", "spatial_unit_name", "parent_region_code",
  "parent_region_name", "spatial_level", "analysis_year",
  "yield_anomaly_percent", "harvested_area_ha", "oni_djf_c"
))

panel <- panel |>
  filter(
    country %in% c("Brazil", "Argentina"),
    between(analysis_year, analysis_start_year, analysis_end_year),
    is.finite(yield_anomaly_percent), is.finite(harvested_area_ha),
    harvested_area_ha > 0, is.finite(oni_djf_c)
  ) |>
  mutate(
    spatial_unit_id = normalize_spatial_id(country, spatial_unit_id),
    derived_region_code = substr(spatial_unit_id, 1, 2),
    parent_region_code = case_when(
      grepl("^[0-9]{2}$", derived_region_code) ~ derived_region_code,
      !is.na(parent_region_code) & nzchar(trimws(as.character(parent_region_code))) ~
        as.character(parent_region_code),
      TRUE ~ as.character(parent_region_name)
    ),
    region_key = paste(country, parent_region_code, sep = "__"),
    unit_key = paste(country, spatial_unit_id, sep = "__")
  ) |>
  group_by(country, unit_key) |>
  filter(n_distinct(analysis_year) >= minimum_local_unit_years) |>
  ungroup()

if ("oni_harvest_year_phase" %in% names(panel)) {
  panel <- panel |> mutate(official_enso_phase = as.character(oni_harvest_year_phase))
} else {
  phase_candidates <- c(
    file.path(project_dir, "data_merged", "soybean_ENSO_country_year_1974_2024.csv"),
    file.path(project_dir, "enso_indices", "ENSO_indices_by_soybean_harvest_year_1974_2024.csv"),
    file.path(project_dir, "data_processed", "ENSO_indices_by_soybean_harvest_year_1974_2024.csv")
  )
  phase_file <- phase_candidates[file.exists(phase_candidates)][1]
  if (is.na(phase_file)) stop("Official ENSO phases were not found.", call. = FALSE)
  phase_lookup <- read_csv(phase_file, show_col_types = FALSE)
  if ("country" %in% names(phase_lookup)) {
    phase_lookup <- phase_lookup |> filter(country %in% c("Brazil", "Argentina"))
  }
  phase_lookup <- phase_lookup |>
    transmute(analysis_year, official_enso_phase = as.character(oni_harvest_year_phase)) |>
    distinct(analysis_year, .keep_all = TRUE)
  panel <- panel |> left_join(phase_lookup, by = "analysis_year")
}

panel <- panel |>
  mutate(
    official_enso_phase = case_when(
      grepl("El", official_enso_phase, ignore.case = TRUE) ~ "El Nino",
      grepl("La", official_enso_phase, ignore.case = TRUE) ~ "La Nina",
      grepl("Neu", official_enso_phase, ignore.case = TRUE) ~ "Neutral",
      TRUE ~ NA_character_
    )
  )
if (anyNA(panel$official_enso_phase)) stop("Missing official ENSO phases.", call. = FALSE)

duplicates <- panel |> count(country, unit_key, analysis_year) |> filter(n > 1)
if (nrow(duplicates)) stop("Duplicated local unit years were found.", call. = FALSE)

# State and province anomalies are weighted by harvested area within each year.
# The outcome remains productivity anomaly in percentage points.
region_annual <- panel |>
  group_by(country, region_key, parent_region_code, parent_region_name, analysis_year) |>
  summarise(
    yield_anomaly_percent = weighted_mean_safe(yield_anomaly_percent, harvested_area_ha),
    harvested_area_ha = sum(harvested_area_ha),
    local_units = n_distinct(unit_key),
    oni_djf_c = first(oni_djf_c),
    oni_values = n_distinct(oni_djf_c),
    official_enso_phase = first(official_enso_phase),
    phase_values = n_distinct(official_enso_phase),
    .groups = "drop"
  )
if (any(region_annual$oni_values != 1) || any(region_annual$phase_values != 1)) {
  stop("ENSO information is inconsistent within a region and year.", call. = FALSE)
}

eligible_regions <- region_annual |>
  group_by(country, region_key) |>
  summarise(valid_years = n_distinct(analysis_year), .groups = "drop") |>
  filter(valid_years >= minimum_region_years)
region_annual <- region_annual |>
  semi_join(eligible_regions, by = c("country", "region_key"))

reference_region_area <- region_annual |>
  group_by(country, region_key) |>
  summarise(
    reference_area_ha = {
      recent <- harvested_area_ha[between(
        analysis_year, reference_area_start_year, reference_area_end_year
      )]
      if (length(recent) && any(is.finite(recent))) mean(recent, na.rm = TRUE)
      else mean(harvested_area_ha, na.rm = TRUE)
    },
    .groups = "drop"
  ) |>
  group_by(country) |>
  mutate(reference_weight_raw = reference_area_ha / sum(reference_area_ha)) |>
  ungroup()

region_annual <- region_annual |>
  left_join(reference_region_area, by = c("country", "region_key")) |>
  group_by(country, analysis_year) |>
  mutate(
    reference_area_coverage = sum(reference_weight_raw),
    score_weight = reference_weight_raw / sum(reference_weight_raw)
  ) |>
  ungroup()

year_table <- region_annual |>
  group_by(analysis_year) |>
  summarise(
    oni_djf_c = first(oni_djf_c),
    oni_values = n_distinct(oni_djf_c),
    official_enso_phase = first(official_enso_phase),
    phase_values = n_distinct(official_enso_phase),
    .groups = "drop"
  )
if (any(year_table$oni_values != 1) || any(year_table$phase_values != 1)) {
  stop("ENSO is not unique within harvest year.", call. = FALSE)
}
expected_counts <- tibble(
  official_enso_phase = c("El Nino", "La Nina", "Neutral"),
  expected_years = c(17L, 19L, 15L)
)
phase_counts <- year_table |>
  count(official_enso_phase, name = "harvest_years") |>
  left_join(expected_counts, by = "official_enso_phase") |>
  mutate(matches_expected = harvest_years == expected_years)
if (!all(phase_counts$matches_expected)) stop("Unexpected ENSO phase counts.", call. = FALSE)

outer_fold_table <- create_event_folds(year_table, outer_folds, "outer_fold")
region_annual <- region_annual |>
  left_join(select(outer_fold_table, analysis_year, phase_simple, episode_id,
                   outer_fold), by = "analysis_year") |>
  mutate(oni_model = oni_djf_c)

# 7. OBSERVED NESTED OUT OF SAMPLE ANALYSIS -----------------------------------
message("Running observed nested spatial field models...")
country_results <- list()
for (country_name in sort(unique(region_annual$country))) {
  message("Processing ", country_name, "...")
  country_results[[country_name]] <- run_nested_country(
    filter(region_annual, country == country_name), save_tuning = TRUE
  )
}
oof <- bind_rows(lapply(country_results, `[[`, "predictions"))
tuning_results <- bind_rows(lapply(country_results, `[[`, "tuning"))
if (any(!is.finite(oof$model_prediction_percent))) {
  stop("Missing outer fold predictions. Review regional training coverage.", call. = FALSE)
}

# 8. FIELD, NATIONAL, REGIONAL AND PHASE METRICS ------------------------------
summarise_field_metrics <- function(x) {
  observed_energy <- weighted_mean_safe(
    (x$yield_anomaly_percent - x$null_prediction_percent)^2, x$score_weight
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
    spatial_field_cv_r_squared = predictive_r2(
      x$yield_anomaly_percent, x$model_prediction_percent,
      x$null_prediction_percent, x$score_weight
    ),
    spatial_signal_energy_fraction = if_else(
      observed_energy > 0, signal_energy / observed_energy, NA_real_
    ),
    mean_signed_component_percent = signed,
    mean_absolute_component_percent = absolute,
    rms_component_percent = sqrt(weighted_mean_safe(
      x$enso_component_percent^2, x$score_weight
    )),
    cancellation_index = if_else(
      absolute > 0, 1 - abs(signed) / absolute, NA_real_
    )
  )
}

bootstrap_field <- function(x) {
  yearly <- x |>
    group_by(analysis_year) |>
    summarise(
      sse_model = sum(score_weight * squared_error_model),
      sse_null = sum(score_weight * squared_error_null),
      signal = sum(score_weight * enso_component_percent^2),
      weight = sum(score_weight),
      .groups = "drop"
    ) |>
    arrange(analysis_year)
  draws <- replicate(bootstrap_repetitions, {
    z <- yearly[sample_blocks(nrow(yearly), bootstrap_block_length), ]
    c(
      r2 = 1 - sum(z$sse_model) / sum(z$sse_null),
      energy = sum(z$signal) / sum(z$sse_null)
    )
  })
  tibble(
    field_r2_ci_lower = quantile(draws["r2", ], 0.025, na.rm = TRUE),
    field_r2_ci_upper = quantile(draws["r2", ], 0.975, na.rm = TRUE),
    signal_energy_ci_lower = quantile(draws["energy", ], 0.025, na.rm = TRUE),
    signal_energy_ci_upper = quantile(draws["energy", ], 0.975, na.rm = TRUE)
  )
}

set.seed(bootstrap_seed)
country_summary <- oof |>
  group_by(country) |>
  group_modify(~ bind_cols(summarise_field_metrics(.x), bootstrap_field(.x))) |>
  ungroup() |>
  mutate(
    spatial_field_cv_r_squared_percent = 100 * spatial_field_cv_r_squared,
    spatial_signal_energy_percent = 100 * spatial_signal_energy_fraction,
    field_r2_ci_lower_percent = 100 * field_r2_ci_lower,
    field_r2_ci_upper_percent = 100 * field_r2_ci_upper,
    signal_energy_ci_lower_percent = 100 * signal_energy_ci_lower,
    signal_energy_ci_upper_percent = 100 * signal_energy_ci_upper
  )

national_year <- oof |>
  group_by(country, analysis_year, phase_simple, episode_id, outer_fold) |>
  summarise(
    observed_national_percent = weighted_mean_safe(yield_anomaly_percent, score_weight),
    predicted_national_percent = weighted_mean_safe(model_prediction_percent, score_weight),
    null_national_percent = weighted_mean_safe(null_prediction_percent, score_weight),
    neutral_national_percent = weighted_mean_safe(neutral_prediction_percent, score_weight),
    enso_national_component_percent = predicted_national_percent - neutral_national_percent,
    reference_area_coverage = first(reference_area_coverage),
    .groups = "drop"
  )
national_summary <- national_year |>
  group_by(country) |>
  summarise(
    national_net_cv_r_squared = predictive_r2(
      observed_national_percent, predicted_national_percent,
      null_national_percent, rep(1, n())
    ),
    national_net_cv_r_squared_percent = 100 * national_net_cv_r_squared,
    national_mean_absolute_component_percent = mean(
      abs(enso_national_component_percent), na.rm = TRUE
    ),
    .groups = "drop"
  )
country_summary <- country_summary |>
  left_join(national_summary, by = "country") |>
  mutate(
    interpretation = paste(
      "Spatial field R2 measures out of sample error reduction before national",
      "cancellation. Signal energy measures footprint magnitude and is not causal variance."
    )
  )

region_summary <- oof |>
  group_by(country, region_key, parent_region_code, parent_region_name) |>
  summarise(
    observations = n(),
    years = n_distinct(analysis_year),
    cv_r_squared = predictive_r2(
      yield_anomaly_percent, model_prediction_percent,
      null_prediction_percent, rep(1, n())
    ),
    cv_r_squared_percent = 100 * cv_r_squared,
    signal_energy_fraction = sum(enso_component_percent^2) /
      sum((yield_anomaly_percent - null_prediction_percent)^2),
    signal_energy_percent = 100 * signal_energy_fraction,
    mean_signed_component_percent = mean(enso_component_percent),
    mean_absolute_component_percent = mean(abs(enso_component_percent)),
    positive_temporal_skill = cv_r_squared > 0,
    .groups = "drop"
  ) |>
  left_join(reference_region_area, by = c("country", "region_key"))

phase_summary <- oof |>
  filter(phase_simple != "Neutral") |>
  group_by(country, phase_simple) |>
  group_modify(~ summarise_field_metrics(.x)) |>
  ungroup() |>
  mutate(
    spatial_field_cv_r_squared_percent = 100 * spatial_field_cv_r_squared,
    spatial_signal_energy_percent = 100 * spatial_signal_energy_fraction
  )

selected_hyperparameters <- tuning_results |>
  filter(selected) |>
  select(country, outer_fold, selected_family = family,
         selected_lambda = lambda, inner_cv_r_squared,
         outer_test_years, outer_training_years)

# Save the observed analysis before the computationally intensive randomization
# test. These files remain available if an external interruption occurs later.
write_csv(oof, file.path(output_dir, "ENSO_multivariate_field_oof_predictions.csv"))
write_csv(country_summary,
          file.path(output_dir, "ENSO_multivariate_field_country_summary.csv"))
write_csv(region_summary,
          file.path(output_dir, "ENSO_multivariate_field_region_summary.csv"))
write_csv(phase_summary,
          file.path(output_dir, "ENSO_multivariate_field_phase_summary.csv"))
write_csv(selected_hyperparameters,
          file.path(output_dir, "ENSO_multivariate_field_selected_hyperparameters.csv"))

# 9. CIRCULAR SHIFT TEST WITH THE SAME NESTED OOF PROCEDURE -------------------
shift_distribution <- tibble()
shift_summary <- tibble()
if (run_circular_shift_test) {
  set.seed(circular_shift_seed)
  possible <- seq_len(nrow(year_table) - 1L)
  shifts <- if (circular_shift_repetitions >= length(possible)) {
    possible
  } else {
    sample(possible, circular_shift_repetitions)
  }
  observed_lookup <- country_summary |>
    select(country, observed_field_r2 = spatial_field_cv_r_squared)
  checkpoint_file <- file.path(
    output_dir, "QC_checkpoint_06c_v1_1_circular_shifts.csv"
  )
  completed_shifts <- if (file.exists(checkpoint_file)) {
    read_csv(checkpoint_file, show_col_types = FALSE) |>
      filter(
        country %in% unique(region_annual$country),
        shift_years %in% shifts,
        is.finite(shifted_spatial_field_cv_r_squared)
      ) |>
      distinct(country, shift_years, .keep_all = TRUE)
  } else {
    tibble()
  }
  annual_oni <- year_table |> arrange(analysis_year) |>
    select(analysis_year, oni_djf_c)
  for (shift_index in seq_along(shifts)) {
    shift <- shifts[shift_index]
    already_complete <- nrow(completed_shifts) > 0 &&
      sum(completed_shifts$shift_years == shift) ==
      n_distinct(region_annual$country)
    if (already_complete) {
      message("Circular shift ", shift_index, " of ", length(shifts),
              " already completed. Resuming...")
      next
    }
    message("Circular shift ", shift_index, " of ", length(shifts), "...")
    shifted_values <- annual_oni$oni_djf_c[
      ((seq_len(nrow(annual_oni)) - 1L - shift) %% nrow(annual_oni)) + 1L
    ]
    lookup <- tibble(
      analysis_year = annual_oni$analysis_year,
      oni_shifted = shifted_values
    )
    shifted_data <- region_annual |>
      select(-oni_model) |>
      left_join(lookup, by = "analysis_year") |>
      mutate(oni_model = oni_shifted) |>
      select(-oni_shifted)
    current_rows <- list()
    current_index <- 0L
    for (country_name in sort(unique(shifted_data$country))) {
      if (nrow(completed_shifts) > 0 && any(
        completed_shifts$country == country_name &
          completed_shifts$shift_years == shift
      )) next
      result <- run_nested_country(
        filter(shifted_data, country == country_name), save_tuning = FALSE
      )
      metric <- summarise_field_metrics(result$predictions)
      current_index <- current_index + 1L
      current_rows[[current_index]] <- tibble(
        country = country_name,
        shift_years = shift,
        shifted_spatial_field_cv_r_squared = metric$spatial_field_cv_r_squared
      )
    }
    completed_shifts <- bind_rows(completed_shifts, bind_rows(current_rows)) |>
      distinct(country, shift_years, .keep_all = TRUE) |>
      arrange(country, shift_years)
    write_csv(completed_shifts, checkpoint_file)
  }
  shift_distribution <- completed_shifts |>
    left_join(observed_lookup, by = "country")
  shift_summary <- shift_distribution |>
    group_by(country) |>
    summarise(
      observed_spatial_field_cv_r_squared = first(observed_field_r2),
      shifts_tested = n(),
      randomization_p_value =
        (1 + sum(shifted_spatial_field_cv_r_squared >=
                   observed_spatial_field_cv_r_squared)) / (n() + 1),
      null_median_spatial_field_cv_r_squared = median(
        shifted_spatial_field_cv_r_squared
      ),
      null_95_percentile_spatial_field_cv_r_squared = quantile(
        shifted_spatial_field_cv_r_squared, 0.95
      ),
      .groups = "drop"
    )
}

# 10. COMPARISON WITH SCRIPT 06B ----------------------------------------------
comparison_06b <- tibble()
file_06b <- file.path(
  project_dir, "enso_spatial_attribution",
  "ENSO_spatially_resolved_country_summary.csv"
)
if (file.exists(file_06b)) {
  old <- read_csv(file_06b, show_col_types = FALSE) |>
    filter(
      weighting == "fixed 2015 to 2024 reference area",
      model == "state_specific_symmetric"
    ) |>
    select(
      country,
      script06b_spatial_cv_r2_percent = cv_r_squared_percent,
      script06b_signal_energy_percent = cv_enso_signal_energy_percent
    )
  comparison_06b <- country_summary |>
    select(
      country,
      script06c_spatial_field_cv_r2_percent = spatial_field_cv_r_squared_percent,
      script06c_signal_energy_percent = spatial_signal_energy_percent,
      script06c_national_net_cv_r2_percent = national_net_cv_r_squared_percent
    ) |>
    left_join(old, by = "country") |>
    mutate(
      same_sign_spatial_skill = sign(script06c_spatial_field_cv_r2_percent) ==
        sign(script06b_spatial_cv_r2_percent)
    )
}

# 11. QUALITY CONTROL AND OUTPUTS ---------------------------------------------
qc <- tribble(
  ~check, ~value, ~status,
  "Retained state and province observations", nrow(region_annual), "Information",
  "Retained states and provinces", n_distinct(region_annual$region_key), "Information",
  "Missing outer fold predictions", sum(!is.finite(oof$model_prediction_percent)),
  if_else(all(is.finite(oof$model_prediction_percent)), "Pass", "Fail"),
  "Minimum region years", min(eligible_regions$valid_years), "Information",
  "Minimum reference area coverage", min(region_annual$reference_area_coverage),
  "Information",
  "Outer folds", n_distinct(oof$outer_fold),
  if_else(n_distinct(oof$outer_fold) == outer_folds, "Pass", "Fail"),
  "Countries", n_distinct(oof$country),
  if_else(n_distinct(oof$country) == 2, "Pass", "Fail")
)

dictionary <- tribble(
  ~metric, ~definition,
  "spatial_field_cv_r_squared_percent",
  "Out of sample reduction in weighted state and province prediction error versus a region mean null model.",
  "national_net_cv_r_squared_percent",
  "Out of sample reduction in national aggregate prediction error after spatial cancellation.",
  "spatial_signal_energy_percent",
  "Squared out of sample ENSO component divided by null residual energy. A footprint magnitude, not causal variance.",
  "cancellation_index",
  "One minus absolute spatial mean response divided by mean absolute local response. Higher values indicate stronger cancellation.",
  "randomization_p_value",
  "Probability under circularly shifted ENSO series, with the complete nested prediction procedure repeated."
)

write_csv(region_annual, file.path(output_dir, "ENSO_state_province_annual_field.csv"))
write_csv(oof, file.path(output_dir, "ENSO_multivariate_field_oof_predictions.csv"))
write_csv(selected_hyperparameters,
          file.path(output_dir, "ENSO_multivariate_field_selected_hyperparameters.csv"))
write_csv(tuning_results,
          file.path(output_dir, "ENSO_multivariate_field_inner_tuning.csv"))
write_csv(country_summary,
          file.path(output_dir, "ENSO_multivariate_field_country_summary.csv"))
write_csv(national_year,
          file.path(output_dir, "ENSO_multivariate_field_national_annual.csv"))
write_csv(region_summary,
          file.path(output_dir, "ENSO_multivariate_field_region_summary.csv"))
write_csv(phase_summary,
          file.path(output_dir, "ENSO_multivariate_field_phase_summary.csv"))
write_csv(shift_summary,
          file.path(output_dir, "ENSO_multivariate_field_circular_shift_test.csv"))
write_csv(shift_distribution,
          file.path(output_dir, "ENSO_multivariate_field_circular_shift_null.csv"))
write_csv(comparison_06b,
          file.path(output_dir, "QC_script06b_script06c_comparison.csv"))
write_csv(phase_counts, file.path(output_dir, "QC_ENSO_phase_counts_06c.csv"))
write_csv(outer_fold_table, file.path(output_dir, "QC_ENSO_event_folds_06c.csv"))
write_csv(qc, file.path(output_dir, "QC_multivariate_field_attribution.csv"))
write_csv(dictionary, file.path(output_dir, "README_multivariate_field_metrics.csv"))

message("Script 06c completed successfully.")
message("Outputs saved in: ", output_dir)
print(country_summary, n = Inf, width = Inf)
print(selected_hyperparameters, n = Inf, width = Inf)
print(shift_summary, n = Inf, width = Inf)

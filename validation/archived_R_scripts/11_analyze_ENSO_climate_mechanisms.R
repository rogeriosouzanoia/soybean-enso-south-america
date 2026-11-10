# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 11_analyze_ENSO_climate_mechanisms.R
# PURPOSE: Quantify ENSO effects on crop-season climate, simple extremes,
#          compound events, and climate-mediated soybean productivity responses
# VERSION: 1.1
# ==============================================================================


# 1. USER SETTINGS -------------------------------------------------------------

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")

project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) project_dir_override else default_project_dir

analysis_start_harvest_year <- 1982L
analysis_end_harvest_year <- 2024L
climatology_start_year <- 1981L
climatology_end_year <- 2010L
wet_day_threshold_mm <- 1
heavy_rain_absolute_threshold_mm <- 20
temperature_upper_probability <- 0.90
temperature_lower_probability <- 0.10
moisture_lower_probability <- 0.20
moisture_upper_probability <- 0.80
heavy_rain_probability <- 0.95
calendar_window_half_width_days <- 2L
rolling_precipitation_days <- 30L
minimum_window_daily_coverage <- 0.95
minimum_years_for_contrast <- 8L
number_temporal_folds <- 5L
bootstrap_repetitions <- 2000L
bootstrap_seed <- 20260906L


# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c(
  "dplyr", "readr", "tibble", "tidyr", "purrr", "stringr",
  "lubridate", "zoo"
)
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
  library(lubridate)
  library(zoo)
})


# 3. PATHS ---------------------------------------------------------------------

climate_dir <- file.path(project_dir, "climate_mechanisms")
weather_dir <- file.path(climate_dir, "NASA_POWER_daily_grid")
merged_dir <- file.path(project_dir, "data_merged")
enso_dir <- file.path(project_dir, "climate_indices", "processed")
output_dir <- climate_dir

grid_manifest_file <- file.path(climate_dir, "NASA_POWER_grid_manifest.csv")
file_index_file <- file.path(climate_dir, "NASA_POWER_daily_file_index.csv")
unit_calendar_file <- file.path(climate_dir, "soybean_spatial_units_with_calendars.csv")
soybean_panel_file <- file.path(merged_dir, "soybean_ENSO_spatial_panel_1974_2024.csv")
enso_file <- file.path(
  enso_dir,
  "ENSO_indices_by_soybean_harvest_year_1974_2024.csv"
)

required_files <- c(
  grid_manifest_file, file_index_file, unit_calendar_file,
  soybean_panel_file, enso_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Required input files are missing: ", paste(missing_files, collapse = "; "), call. = FALSE)
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
  ifelse(is.na(value), NA_character_, sprintf(paste0("%0", width, "d"), as.integer(value)))
}

standardize_unit_id <- function(country, x) {
  ifelse(
    country == "Brazil",
    standardize_numeric_code(x, 7L),
    standardize_numeric_code(x, 5L)
  )
}

month_in_window <- function(month_value, start_month, end_month) {
  if (start_month <= end_month) {
    month_value >= start_month & month_value <= end_month
  } else {
    month_value >= start_month | month_value <= end_month
  }
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NA_real_)
  as.numeric(quantile(x, probability, na.rm = TRUE, names = FALSE, type = 8))
}

maximum_consecutive_true <- function(x) {
  x[is.na(x)] <- FALSE
  if (!any(x)) return(0L)
  max(rle(x)$lengths[rle(x)$values])
}

safe_weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  weighted.mean(x[keep], w[keep])
}

safe_r_squared <- function(observed, predicted, weights = NULL) {
  keep <- is.finite(observed) & is.finite(predicted)
  if (!is.null(weights)) keep <- keep & is.finite(weights) & weights > 0
  observed <- observed[keep]
  predicted <- predicted[keep]
  if (length(observed) < 3) return(NA_real_)
  if (is.null(weights)) weights <- rep(1, length(keep)) else weights <- weights[keep]
  weights <- weights / sum(weights)
  baseline <- sum(weights * observed)
  sse <- sum(weights * (observed - predicted)^2)
  sst <- sum(weights * (observed - baseline)^2)
  if (!is.finite(sst) || sst <= 0) return(NA_real_)
  1 - sse / sst
}

make_temporal_folds <- function(years, phases) {
  year_table <- tibble(
    analysis_year = sort(unique(years))
  ) |>
    left_join(
      tibble(analysis_year = years, phase = phases) |>
        distinct(analysis_year, .keep_all = TRUE),
      by = "analysis_year"
    ) |>
    arrange(analysis_year) |>
    mutate(
      phase_sequence = ave(
        analysis_year, phase,
        FUN = function(x) seq_along(x)
      ),
      outer_fold = ((phase_sequence - 1L) %% number_temporal_folds) + 1L
    )
  year_table
}

fit_predict_lm <- function(formula, training, testing, weights_column) {
  # Weighted least squares solved through an SVD Moore-Penrose inverse.
  # This gives the same fitted projection as ordinary least squares when the
  # design matrix has full rank and remains defined under exact collinearity.
  # It therefore prevents predict.lm from discarding otherwise usable rows.
  result <- tryCatch({
    training_frame <- model.frame(
      formula, data = training, na.action = na.pass,
      drop.unused.levels = FALSE
    )
    formula_terms <- terms(training_frame)
    training_matrix <- model.matrix(formula_terms, training_frame)
    training_response <- model.response(training_frame)

    predictor_terms <- delete.response(formula_terms)
    training_factor_levels <- lapply(
      training_frame[vapply(training_frame, is.factor, logical(1))],
      levels
    )
    testing_frame <- model.frame(
      predictor_terms, data = testing, na.action = na.pass,
      xlev = training_factor_levels,
      drop.unused.levels = FALSE
    )
    testing_matrix <- model.matrix(
      predictor_terms,
      testing_frame,
      contrasts.arg = attr(training_matrix, "contrasts")
    )

    missing_columns <- setdiff(colnames(training_matrix), colnames(testing_matrix))
    if (length(missing_columns) > 0) {
      zero_columns <- matrix(
        0, nrow = nrow(testing_matrix), ncol = length(missing_columns),
        dimnames = list(NULL, missing_columns)
      )
      testing_matrix <- cbind(testing_matrix, zero_columns)
    }
    testing_matrix <- testing_matrix[, colnames(training_matrix), drop = FALSE]

    training_weights <- training[[weights_column]]
    valid_training <- is.finite(training_response) &
      is.finite(training_weights) & training_weights > 0 &
      apply(training_matrix, 1L, function(x) all(is.finite(x)))
    if (sum(valid_training) < 3) stop("Insufficient finite training rows")

    weighted_matrix <- training_matrix[valid_training, , drop = FALSE] *
      sqrt(training_weights[valid_training])
    weighted_response <- training_response[valid_training] *
      sqrt(training_weights[valid_training])
    decomposition <- svd(weighted_matrix, nu = min(dim(weighted_matrix)))
    tolerance <- max(dim(weighted_matrix)) * max(decomposition$d) *
      .Machine$double.eps
    retained <- decomposition$d > tolerance
    if (!any(retained)) stop("No estimable model dimension")

    coefficients <- decomposition$v[, retained, drop = FALSE] %*%
      ((crossprod(
        decomposition$u[, retained, drop = FALSE], weighted_response
      )) / decomposition$d[retained])
    prediction <- as.numeric(testing_matrix %*% coefficients)
    valid_testing <- apply(testing_matrix, 1L, function(x) all(is.finite(x)))
    prediction[!valid_testing | !is.finite(prediction)] <- NA_real_
    prediction
  }, error = function(e) rep(NA_real_, nrow(testing)))

  if (length(result) != nrow(testing)) return(rep(NA_real_, nrow(testing)))
  result
}

safe_weighted_mse <- function(observed, predicted, weights) {
  keep <- is.finite(observed) & is.finite(predicted) &
    is.finite(weights) & weights > 0
  if (sum(keep) < 3) return(NA_real_)
  weighted.mean((observed[keep] - predicted[keep])^2, weights[keep])
}

bootstrap_phase_difference <- function(data, target_phase, value_column, weight_column) {
  target <- data |>
    filter(enso_phase_simple == target_phase)
  neutral <- data |>
    filter(enso_phase_simple == "Neutral")
  if (n_distinct(target$harvest_year) < 3 || n_distinct(neutral$harvest_year) < 3) {
    return(tibble(ci_lower = NA_real_, ci_upper = NA_real_))
  }
  target_values <- target[[value_column]]
  target_weights <- target[[weight_column]]
  neutral_values <- neutral[[value_column]]
  neutral_weights <- neutral[[weight_column]]
  simulated <- replicate(bootstrap_repetitions, {
    target_indices <- sample.int(
      length(target_values), length(target_values), replace = TRUE
    )
    neutral_indices <- sample.int(
      length(neutral_values), length(neutral_values), replace = TRUE
    )
    safe_weighted_mean(
      target_values[target_indices], target_weights[target_indices]
    ) - safe_weighted_mean(
      neutral_values[neutral_indices], neutral_weights[neutral_indices]
    )
  })
  tibble(
    ci_lower = quantile(simulated, 0.025, na.rm = TRUE, names = FALSE),
    ci_upper = quantile(simulated, 0.975, na.rm = TRUE, names = FALSE)
  )
}


# 5. READ AND VALIDATE INPUTS --------------------------------------------------

message("Reading spatial, climate, ENSO, and soybean inputs...")

grid_manifest <- read_csv(
  grid_manifest_file,
  col_types = cols(
    climate_location_id = col_character(), country = col_character(),
    spatial_unit_id = col_character(), parent_region_code = col_character(),
    grid_id = col_character(), .default = col_guess()
  ),
  show_col_types = FALSE
)

file_index <- read_csv(file_index_file, show_col_types = FALSE)
unit_calendars <- read_csv(
  unit_calendar_file,
  col_types = cols(
    climate_location_id = col_character(), country = col_character(),
    spatial_unit_id = col_character(), parent_region_code = col_character(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)
soybean <- read_csv(soybean_panel_file, show_col_types = FALSE)
enso <- read_csv(enso_file, show_col_types = FALSE)

check_required_columns(
  grid_manifest,
  c(
    "climate_location_id", "country", "spatial_unit_id", "grid_id",
    "power_latitude", "power_longitude", "reference_harvested_area_ha",
    "calendar_group"
  ), "NASA POWER grid manifest"
)
check_required_columns(
  unit_calendars,
  c(
    "climate_location_id", "country", "spatial_unit_id", "parent_region_code",
    "parent_region_name", "calendar_group", "growing_season_start_month",
    "growing_season_end_month", "vegetative_start_month",
    "vegetative_end_month", "reproductive_start_month",
    "reproductive_end_month"
  ), "spatial unit crop calendars"
)
check_required_columns(
  soybean,
  c(
    "country", "spatial_unit_id", "analysis_year", "yield_anomaly_percent",
    "harvested_area_ha", "oni_djf_c", "oni_harvest_year_phase"
  ), "soybean ENSO panel"
)
check_required_columns(
  enso,
  c("harvest_year", "oni_harvest_year_phase", "oni_djf_c"),
  "ENSO harvest-year table"
)

if (!all(file_index$complete)) {
  stop("NASA POWER file index contains incomplete cells. Finish script 10 first.", call. = FALSE)
}

missing_grid_files <- setdiff(
  unique(as.character(grid_manifest$grid_id)),
  unique(as.character(file_index$grid_id[file_index$complete]))
)
if (length(missing_grid_files) > 0) {
  stop(
    length(missing_grid_files),
    " current spatial-grid cells are absent from the completed NASA POWER index. ",
    "Rerun script 10 with the current spatial manifest before script 11.",
    call. = FALSE
  )
}

calendar_table <- unit_calendars |>
  distinct(
    country, calendar_group,
    growing_season_start_month, growing_season_end_month,
    vegetative_start_month, vegetative_end_month,
    reproductive_start_month, reproductive_end_month
  )

grid_calendar <- grid_manifest |>
  distinct(grid_id, country, calendar_group) |>
  left_join(
    calendar_table,
    by = c("country", "calendar_group")
  )

if (any(is.na(grid_calendar$reproductive_start_month))) {
  stop("At least one grid and calendar combination has incomplete calendar fields.", call. = FALSE)
}

enso_analysis <- enso |>
  transmute(
    harvest_year = as.integer(harvest_year),
    oni_djf_c = as.numeric(oni_djf_c),
    enso_phase_simple = case_when(
      as.character(oni_harvest_year_phase) == "El Nino" ~ "El Nino",
      as.character(oni_harvest_year_phase) == "La Nina" ~ "La Nina",
      TRUE ~ "Neutral"
    )
  ) |>
  filter(
    harvest_year >= analysis_start_harvest_year,
    harvest_year <= analysis_end_harvest_year
  )


# 6. DAILY CLIMATOLOGICAL THRESHOLDS AND CROP-WINDOW INDICES -------------------

message("Calculating local climate thresholds and crop-window indices...")

grid_ids <- sort(unique(grid_manifest$grid_id))
threshold_results <- vector("list", length(grid_ids))
index_results <- vector("list", length(grid_ids))
daily_qc_results <- vector("list", length(grid_ids))

for (grid_position in seq_along(grid_ids)) {
  grid_id <- grid_ids[[grid_position]]
  if (grid_position == 1L || grid_position %% 20L == 0L || grid_position == length(grid_ids)) {
    message("Climate grid ", grid_position, " of ", length(grid_ids), "...")
  }

  weather_file <- file.path(weather_dir, paste0(grid_id, ".csv.gz"))
  if (!file.exists(weather_file)) stop("Weather file missing: ", weather_file, call. = FALSE)

  daily <- read_csv(
    weather_file,
    col_types = cols(
      grid_id = col_character(), power_latitude = col_double(),
      power_longitude = col_double(), date = col_date(), .default = col_double()
    ),
    progress = FALSE
  ) |>
    arrange(date) |>
    mutate(
      year = year(date),
      month = month(date),
      harvest_year = if_else(month >= 8L, year + 1L, year),
      climate_day = yday(as.Date(sprintf("2000-%02d-%02d", month, day(date)))),
      precipitation_30d_mm = zoo::rollsumr(
        PRECTOTCORR,
        k = rolling_precipitation_days,
        fill = NA_real_, na.rm = FALSE
      ),
      vapour_pressure_deficit_kpa = pmax(
        0,
        0.6108 * exp(17.27 * T2M / (T2M + 237.3)) -
          0.6108 * exp(17.27 * T2MDEW / (T2MDEW + 237.3))
      )
    )

  baseline <- daily |>
    filter(year >= climatology_start_year, year <= climatology_end_year) |>
    select(climate_day, T2M_MAX, T2M_MIN, PRECTOTCORR, precipitation_30d_mm)

  expanded_baseline <- map_dfr(
    seq.int(-calendar_window_half_width_days, calendar_window_half_width_days),
    function(day_shift) {
      baseline |>
        mutate(target_climate_day = ((climate_day - 1L - day_shift) %% 366L) + 1L)
    }
  )

  thresholds <- expanded_baseline |>
    group_by(target_climate_day) |>
    summarise(
      tmax_p90_c = safe_quantile(T2M_MAX, temperature_upper_probability),
      tmin_p10_c = safe_quantile(T2M_MIN, temperature_lower_probability),
      precipitation_30d_p20_mm = safe_quantile(
        precipitation_30d_mm, moisture_lower_probability
      ),
      precipitation_30d_p80_mm = safe_quantile(
        precipitation_30d_mm, moisture_upper_probability
      ),
      wet_day_precipitation_p95_mm = safe_quantile(
        PRECTOTCORR[PRECTOTCORR >= wet_day_threshold_mm],
        heavy_rain_probability
      ),
      .groups = "drop"
    ) |>
    mutate(grid_id = grid_id, .before = 1)

  daily_classified <- daily |>
    left_join(
      thresholds,
      by = c("climate_day" = "target_climate_day")
    ) |>
    mutate(
      hot_day = T2M_MAX > tmax_p90_c,
      cold_day = T2M_MIN < tmin_p10_c,
      dry_condition = precipitation_30d_mm < precipitation_30d_p20_mm,
      wet_condition = precipitation_30d_mm > precipitation_30d_p80_mm,
      heavy_rain_day = PRECTOTCORR > wet_day_precipitation_p95_mm,
      hot_dry_day = hot_day & dry_condition,
      hot_wet_day = hot_day & wet_condition,
      cold_dry_day = cold_day & dry_condition,
      cold_wet_day = cold_day & wet_condition
    )

  calendar_rows <- grid_calendar |>
    filter(.data$grid_id == .env$grid_id)

  grid_indices <- map_dfr(seq_len(nrow(calendar_rows)), function(calendar_position) {
    calendar <- calendar_rows[calendar_position, ]

    window_definitions <- tribble(
      ~climate_window, ~start_month, ~end_month,
      "growing_season", calendar$growing_season_start_month[[1]], calendar$growing_season_end_month[[1]],
      "vegetative", calendar$vegetative_start_month[[1]], calendar$vegetative_end_month[[1]],
      "reproductive", calendar$reproductive_start_month[[1]], calendar$reproductive_end_month[[1]]
    )

    map_dfr(seq_len(nrow(window_definitions)), function(window_position) {
      window <- window_definitions[window_position, ]
      selected <- daily_classified |>
        filter(
          month_in_window(month, window$start_month[[1]], window$end_month[[1]]),
          harvest_year >= analysis_start_harvest_year,
          harvest_year <= analysis_end_harvest_year
        )

      selected |>
        group_by(harvest_year) |>
        summarise(
          days_expected = n(),
          days_with_core_climate = sum(
            is.finite(T2M) & is.finite(T2M_MAX) &
              is.finite(T2M_MIN) & is.finite(PRECTOTCORR)
          ),
          daily_coverage = days_with_core_climate / days_expected,
          mean_temperature_c = mean(T2M, na.rm = TRUE),
          mean_maximum_temperature_c = mean(T2M_MAX, na.rm = TRUE),
          maximum_temperature_c = max(T2M_MAX, na.rm = TRUE),
          mean_minimum_temperature_c = mean(T2M_MIN, na.rm = TRUE),
          minimum_temperature_c = min(T2M_MIN, na.rm = TRUE),
          mean_vpd_kpa = mean(vapour_pressure_deficit_kpa, na.rm = TRUE),
          precipitation_total_mm = sum(PRECTOTCORR, na.rm = TRUE),
          rainy_days = sum(PRECTOTCORR >= wet_day_threshold_mm, na.rm = TRUE),
          heavy_rain_days_20mm = sum(
            PRECTOTCORR >= heavy_rain_absolute_threshold_mm, na.rm = TRUE
          ),
          percentile_heavy_rain_days = sum(heavy_rain_day, na.rm = TRUE),
          maximum_one_day_precipitation_mm = max(PRECTOTCORR, na.rm = TRUE),
          maximum_consecutive_dry_days = maximum_consecutive_true(
            PRECTOTCORR < wet_day_threshold_mm
          ),
          hot_days = sum(hot_day, na.rm = TRUE),
          cold_days = sum(cold_day, na.rm = TRUE),
          dry_condition_days = sum(dry_condition, na.rm = TRUE),
          wet_condition_days = sum(wet_condition, na.rm = TRUE),
          hot_dry_days = sum(hot_dry_day, na.rm = TRUE),
          hot_wet_days = sum(hot_wet_day, na.rm = TRUE),
          cold_dry_days = sum(cold_dry_day, na.rm = TRUE),
          cold_wet_days = sum(cold_wet_day, na.rm = TRUE),
          .groups = "drop"
        ) |>
        mutate(
          hot_day_fraction = hot_days / days_expected,
          cold_day_fraction = cold_days / days_expected,
          dry_condition_fraction = dry_condition_days / days_expected,
          wet_condition_fraction = wet_condition_days / days_expected,
          heavy_rain_fraction = percentile_heavy_rain_days / days_expected,
          hot_dry_fraction = hot_dry_days / days_expected,
          hot_wet_fraction = hot_wet_days / days_expected,
          cold_dry_fraction = cold_dry_days / days_expected,
          cold_wet_fraction = cold_wet_days / days_expected,
          grid_id = grid_id,
          country = calendar$country[[1]],
          calendar_group = calendar$calendar_group[[1]],
          climate_window = window$climate_window[[1]],
          .before = 1
        )
    })
  })

  threshold_results[[grid_position]] <- thresholds
  index_results[[grid_position]] <- grid_indices
  daily_qc_results[[grid_position]] <- tibble(
    grid_id = grid_id,
    daily_rows = nrow(daily),
    first_date = min(daily$date),
    last_date = max(daily$date),
    threshold_days = nrow(thresholds),
    grid_calendar_combinations = nrow(calendar_rows),
    crop_window_rows = nrow(grid_indices),
    minimum_window_coverage = min(grid_indices$daily_coverage, na.rm = TRUE)
  )
}

climate_thresholds <- bind_rows(threshold_results)
grid_year_indices_raw <- bind_rows(index_results) |>
  left_join(enso_analysis, by = "harvest_year")
daily_processing_qc <- bind_rows(daily_qc_results)

if (any(grid_year_indices_raw$daily_coverage < minimum_window_daily_coverage)) {
  stop("At least one crop window has insufficient daily coverage.", call. = FALSE)
}
if (any(is.na(grid_year_indices_raw$enso_phase_simple))) {
  stop("ENSO phase is missing after joining climate indices by harvest year.", call. = FALSE)
}


# 7. DETREND AND STANDARDIZE CLIMATE INDICES ----------------------------------

metric_columns <- c(
  "mean_temperature_c", "mean_maximum_temperature_c",
  "maximum_temperature_c", "mean_minimum_temperature_c",
  "minimum_temperature_c", "mean_vpd_kpa", "precipitation_total_mm",
  "rainy_days", "heavy_rain_days_20mm", "percentile_heavy_rain_days",
  "maximum_one_day_precipitation_mm", "maximum_consecutive_dry_days",
  "hot_day_fraction", "cold_day_fraction", "dry_condition_fraction",
  "wet_condition_fraction", "heavy_rain_fraction", "hot_dry_fraction",
  "hot_wet_fraction", "cold_dry_fraction", "cold_wet_fraction"
)

message("Removing long-term climate trends and standardizing local indices...")

grid_year_indices_long <- grid_year_indices_raw |>
  pivot_longer(
    cols = all_of(metric_columns),
    names_to = "climate_metric",
    values_to = "climate_value"
  ) |>
  group_by(grid_id, country, calendar_group, climate_window, climate_metric) |>
  group_modify(function(data, keys) {
    valid <- is.finite(data$climate_value) & is.finite(data$harvest_year)
    anomaly <- rep(NA_real_, nrow(data))
    if (sum(valid) >= minimum_years_for_contrast) {
      fit <- lm(climate_value ~ harvest_year, data = data[valid, ])
      anomaly[valid] <- residuals(fit)
    }
    anomaly_sd <- sd(anomaly, na.rm = TRUE)
    data |>
      mutate(
        climate_anomaly = anomaly,
        climate_anomaly_standardized = case_when(
          !valid ~ NA_real_,
          is.finite(anomaly_sd) && anomaly_sd > 0 ~ anomaly / anomaly_sd,
          TRUE ~ 0
        )
      )
  }) |>
  ungroup()

grid_year_indices <- grid_year_indices_long |>
  select(
    grid_id, country, calendar_group,
    climate_window, harvest_year, oni_djf_c, enso_phase_simple,
    climate_metric, climate_value, climate_anomaly,
    climate_anomaly_standardized
  )


# 8. GRID-LEVEL ENSO PHASE CONTRASTS ------------------------------------------

message("Estimating local ENSO phase contrasts in climate...")

grid_phase_contrasts <- grid_year_indices |>
  filter(enso_phase_simple %in% c("El Nino", "La Nina", "Neutral")) |>
  group_by(
    grid_id, country, calendar_group, climate_window, climate_metric
  ) |>
  group_modify(function(data, keys) {
    map_dfr(c("El Nino", "La Nina"), function(target_phase) {
      target <- data$climate_anomaly[data$enso_phase_simple == target_phase]
      neutral <- data$climate_anomaly[data$enso_phase_simple == "Neutral"]
      target_z <- data$climate_anomaly_standardized[data$enso_phase_simple == target_phase]
      neutral_z <- data$climate_anomaly_standardized[data$enso_phase_simple == "Neutral"]
      target <- target[is.finite(target)]
      neutral <- neutral[is.finite(neutral)]
      target_z <- target_z[is.finite(target_z)]
      neutral_z <- neutral_z[is.finite(neutral_z)]
      test <- if (length(target) >= 3 && length(neutral) >= 3) {
        tryCatch(t.test(target, neutral), error = function(e) NULL)
      } else NULL
      tibble(
        enso_phase = target_phase,
        phase_years = length(target),
        neutral_years = length(neutral),
        phase_difference = mean(target, na.rm = TRUE) - mean(neutral, na.rm = TRUE),
        standardized_phase_difference = mean(target_z, na.rm = TRUE) - mean(neutral_z, na.rm = TRUE),
        standard_error = if (is.null(test)) NA_real_ else unname(test$stderr),
        ci_lower = if (is.null(test)) NA_real_ else unname(test$conf.int[1]),
        ci_upper = if (is.null(test)) NA_real_ else unname(test$conf.int[2]),
        p_value = if (is.null(test)) NA_real_ else test$p.value
      )
    })
  }) |>
  ungroup() |>
  group_by(country, climate_window, climate_metric, enso_phase) |>
  mutate(fdr_adjusted_p_value = p.adjust(p_value, method = "BH")) |>
  ungroup()


# 9. COUNTRY AND STATE/PROVINCE CLIMATE COMPOSITES ----------------------------

grid_weights <- grid_manifest |>
  group_by(
    grid_id, country, parent_region_code, parent_region_name, calendar_group
  ) |>
  summarise(
    climate_weight = sum(reference_harvested_area_ha, na.rm = TRUE),
    spatial_units = n(),
    .groups = "drop"
  )

weighted_climate <- grid_year_indices |>
  left_join(
    grid_weights,
    by = c("grid_id", "country", "calendar_group")
  )

country_year_climate <- weighted_climate |>
  group_by(
    country, climate_window, climate_metric, harvest_year,
    enso_phase_simple, oni_djf_c
  ) |>
  summarise(
    climate_value = safe_weighted_mean(climate_value, climate_weight),
    climate_anomaly = safe_weighted_mean(climate_anomaly, climate_weight),
    climate_anomaly_standardized = safe_weighted_mean(
      climate_anomaly_standardized, climate_weight
    ),
    climate_weight = sum(climate_weight, na.rm = TRUE),
    grid_cells = n_distinct(grid_id),
    .groups = "drop"
  )

state_year_climate <- weighted_climate |>
  group_by(
    country, parent_region_code, parent_region_name, climate_window,
    climate_metric, harvest_year, enso_phase_simple, oni_djf_c
  ) |>
  summarise(
    climate_value = safe_weighted_mean(climate_value, climate_weight),
    climate_anomaly = safe_weighted_mean(climate_anomaly, climate_weight),
    climate_anomaly_standardized = safe_weighted_mean(
      climate_anomaly_standardized, climate_weight
    ),
    climate_weight = sum(climate_weight, na.rm = TRUE),
    grid_cells = n_distinct(grid_id),
    .groups = "drop"
  )

set.seed(bootstrap_seed)

country_phase_summary <- country_year_climate |>
  group_by(country, climate_window, climate_metric) |>
  group_modify(function(data, keys) {
    map_dfr(c("El Nino", "La Nina"), function(target_phase) {
      target <- filter(data, enso_phase_simple == target_phase)
      neutral <- filter(data, enso_phase_simple == "Neutral")
      difference <- safe_weighted_mean(target$climate_anomaly, target$climate_weight) -
        safe_weighted_mean(neutral$climate_anomaly, neutral$climate_weight)
      interval <- bootstrap_phase_difference(
        data, target_phase, "climate_anomaly", "climate_weight"
      )
      tibble(
        enso_phase = target_phase,
        phase_years = n_distinct(target$harvest_year),
        neutral_years = n_distinct(neutral$harvest_year),
        phase_difference = difference,
        ci_lower = interval$ci_lower,
        ci_upper = interval$ci_upper
      )
    })
  }) |>
  ungroup()

state_phase_summary <- state_year_climate |>
  group_by(
    country, parent_region_code, parent_region_name,
    climate_window, climate_metric
  ) |>
  group_modify(function(data, keys) {
    map_dfr(c("El Nino", "La Nina"), function(target_phase) {
      target <- filter(data, enso_phase_simple == target_phase)
      neutral <- filter(data, enso_phase_simple == "Neutral")
      tibble(
        enso_phase = target_phase,
        phase_years = n_distinct(target$harvest_year),
        neutral_years = n_distinct(neutral$harvest_year),
        phase_difference = mean(target$climate_anomaly, na.rm = TRUE) -
          mean(neutral$climate_anomaly, na.rm = TRUE),
        standardized_phase_difference = mean(
          target$climate_anomaly_standardized, na.rm = TRUE
        ) - mean(neutral$climate_anomaly_standardized, na.rm = TRUE)
      )
    })
  }) |>
  ungroup()

compound_metrics <- c(
  "hot_dry_fraction", "hot_wet_fraction",
  "cold_dry_fraction", "cold_wet_fraction"
)

compound_frequency_summary <- country_phase_summary |>
  filter(
    climate_window %in% c("growing_season", "reproductive"),
    climate_metric %in% compound_metrics
  ) |>
  mutate(
    phase_difference_percentage_points = 100 * phase_difference,
    ci_lower_percentage_points = 100 * ci_lower,
    ci_upper_percentage_points = 100 * ci_upper
  ) |>
  group_by(country, climate_window, enso_phase) |>
  mutate(
    absolute_change_rank = rank(
      -abs(phase_difference_percentage_points), ties.method = "first"
    ),
    most_responsive_compound = absolute_change_rank == 1
  ) |>
  ungroup()


# 10. CLIMATE AND PRODUCTIVITY PANEL ------------------------------------------

message("Linking climate mechanisms to soybean productivity...")

primary_climate_metrics <- c(
  "mean_temperature_c", "precipitation_total_mm", "mean_vpd_kpa",
  "hot_dry_fraction", "hot_wet_fraction", "cold_wet_fraction",
  "heavy_rain_fraction", "maximum_consecutive_dry_days"
)

unit_climate <- grid_manifest |>
  select(
    climate_location_id, country, spatial_unit_id, parent_region_code,
    parent_region_name, grid_id, calendar_group
  ) |>
  inner_join(
    grid_year_indices |>
      filter(
        climate_window == "reproductive",
        climate_metric %in% primary_climate_metrics
      ) |>
      select(
        grid_id, country, calendar_group, harvest_year, climate_metric,
        climate_anomaly_standardized
      ),
    by = c("grid_id", "country", "calendar_group")
  ) |>
  select(-climate_location_id) |>
  pivot_wider(
    names_from = climate_metric,
    values_from = climate_anomaly_standardized,
    names_prefix = "climate_z_"
  )

yield_climate_panel <- soybean |>
  mutate(
    country = as.character(country),
    spatial_unit_id = standardize_unit_id(country, spatial_unit_id),
    parent_region_code = standardize_numeric_code(parent_region_code, 2L),
    analysis_year = as.integer(analysis_year)
  ) |>
  filter(
    analysis_year >= analysis_start_harvest_year,
    analysis_year <= analysis_end_harvest_year
  ) |>
  inner_join(
    unit_climate,
    by = c("country", "spatial_unit_id", "analysis_year" = "harvest_year"),
    suffix = c("", "_climate")
  )

if (nrow(yield_climate_panel) == 0) {
  stop("The climate and soybean panels did not match.", call. = FALSE)
}

state_yield_climate <- yield_climate_panel |>
  mutate(
    model_region = interaction(country, parent_region_name, drop = TRUE),
    model_weight = pmax(harvested_area_ha, 0)
  ) |>
  group_by(
    country, parent_region_name, model_region, analysis_year,
    oni_djf_c, oni_harvest_year_phase
  ) |>
  summarise(
    yield_anomaly_percent = safe_weighted_mean(
      yield_anomaly_percent, model_weight
    ),
    across(
      starts_with("climate_z_"),
      ~ safe_weighted_mean(.x, model_weight)
    ),
    model_weight = sum(model_weight, na.rm = TRUE),
    spatial_units = n_distinct(spatial_unit_id),
    .groups = "drop"
  ) |>
  mutate(
    enso_phase_simple = case_when(
      as.character(oni_harvest_year_phase) == "El Nino" ~ "El Nino",
      as.character(oni_harvest_year_phase) == "La Nina" ~ "La Nina",
      TRUE ~ "Neutral"
    )
  )


# 11. EVENT-BLOCKED CLIMATE MEDIATION MODELS ----------------------------------

model_predictors <- paste0("climate_z_", primary_climate_metrics)
model_predictors <- model_predictors[model_predictors %in% names(state_yield_climate)]

climate_formula <- as.formula(paste(
  "yield_anomaly_percent ~ factor(model_region) +",
  paste(model_predictors, collapse = " + ")
))
null_formula <- yield_anomaly_percent ~ factor(model_region)
enso_formula <- yield_anomaly_percent ~ factor(model_region) +
  oni_djf_c:factor(model_region)
full_formula <- update(
  climate_formula,
  paste(". ~ . + oni_djf_c:factor(model_region)")
)

model_comparison <- map_dfr(c("Argentina", "Brazil"), function(selected_country) {
  country_data <- state_yield_climate |>
    filter(country == selected_country) |>
    mutate(model_region = droplevels(factor(model_region)))

  folds <- make_temporal_folds(
    country_data$analysis_year,
    country_data$enso_phase_simple
  )
  country_data <- country_data |>
    left_join(folds, by = "analysis_year")

  predictions <- map_dfr(seq_len(number_temporal_folds), function(fold_id) {
    training <- filter(country_data, outer_fold != fold_id)
    testing <- filter(country_data, outer_fold == fold_id)
    testing |>
      transmute(
        country,
        model_region,
        analysis_year,
        observed = yield_anomaly_percent,
        model_weight,
        outer_fold,
        prediction_null = fit_predict_lm(
          null_formula, training, testing, "model_weight"
        ),
        prediction_enso = fit_predict_lm(
          enso_formula, training, testing, "model_weight"
        ),
        prediction_climate = fit_predict_lm(
          climate_formula, training, testing, "model_weight"
        ),
        prediction_full = fit_predict_lm(
          full_formula, training, testing, "model_weight"
        )
      )
  })

  model_names <- c("null", "ENSO", "climate", "ENSO_plus_climate")
  prediction_names <- paste0("prediction_", c("null", "enso", "climate", "full"))
  validation_rows <- map_dfr(seq_along(model_names), function(model_index) {
    prediction <- predictions[[prediction_names[model_index]]]
    keep <- is.finite(predictions$observed) & is.finite(prediction) &
      is.finite(predictions$model_weight) & predictions$model_weight > 0
    tibble(
      country = selected_country,
      model = model_names[model_index],
      validation_observations = sum(keep),
      validation_years = n_distinct(predictions$analysis_year[keep]),
      prediction_coverage = mean(is.finite(prediction)),
      validation_mse = safe_weighted_mse(
        predictions$observed, prediction, predictions$model_weight
      ),
      out_of_sample_r_squared = safe_r_squared(
        predictions$observed, prediction, predictions$model_weight
      )
    )
  })

  baseline_error <- validation_rows$validation_mse[validation_rows$model == "null"]

  validation_rows |>
    mutate(
    observations = nrow(predictions),
    years = n_distinct(predictions$analysis_year),
    regions = n_distinct(predictions$model_region),
    error_reduction_relative_to_null = if_else(
      model == "null",
      0,
      1 - validation_mse / baseline_error
    )
  )
}) |>
  group_by(country) |>
  mutate(
    incremental_climate_skill_beyond_ENSO = if_else(
      model == "ENSO_plus_climate",
      out_of_sample_r_squared - out_of_sample_r_squared[model == "ENSO"],
      NA_real_
    ),
    incremental_ENSO_skill_beyond_climate = if_else(
      model == "ENSO_plus_climate",
      out_of_sample_r_squared - out_of_sample_r_squared[model == "climate"],
      NA_real_
    )
  ) |>
  ungroup()


# 12. QUALITY CONTROL ----------------------------------------------------------

qc_summary <- tibble(
  test = c(
    "NASA POWER grid cells processed",
    "Daily threshold rows",
    "Grid crop-window year rows",
    "Minimum crop-window daily coverage",
    "Harvest years represented",
    "Climate and yield panel rows",
    "State and province model rows",
    "Missing primary climate predictors",
    "Duplicated grid-year-window-metric rows"
  ),
  value = c(
    length(grid_ids),
    nrow(climate_thresholds),
    nrow(grid_year_indices_raw),
    min(grid_year_indices_raw$daily_coverage, na.rm = TRUE),
    n_distinct(grid_year_indices_raw$harvest_year),
    nrow(yield_climate_panel),
    nrow(state_yield_climate),
    sum(!complete.cases(state_yield_climate[, model_predictors])),
    sum(duplicated(grid_year_indices[, c(
      "grid_id", "calendar_group", "climate_window",
      "harvest_year", "climate_metric"
    )]))
  )
)

methodology <- tribble(
  ~component, ~definition,
  "Reference climatology", "1981 to 2010",
  "Calendar-day thresholds", "Local five-day moving window centred on each day",
  "Hot day", "Daily maximum temperature above local calendar-day 90th percentile",
  "Cold day", "Daily minimum temperature below local calendar-day 10th percentile",
  "Dry condition", "Previous 30-day precipitation below local calendar-day 20th percentile",
  "Wet condition", "Previous 30-day precipitation above local calendar-day 80th percentile",
  "Heavy rain", "Daily precipitation above local wet-day 95th percentile",
  "Hot and dry", "Hot day coinciding with dry condition",
  "Hot and wet", "Hot day coinciding with wet condition",
  "Cold and dry", "Cold day coinciding with dry condition",
  "Cold and wet", "Cold day coinciding with wet condition",
  "Trend adjustment", "Linear trend removed separately from every grid, calendar, window, and metric",
  "ENSO contrast", "Detrended climate anomaly in El Nino or La Nina minus neutral years",
  "Yield mechanism validation", "Five event-stratified temporal folds at state and province level"
)


# 13. SAVE OUTPUTS -------------------------------------------------------------

message("Saving climate mechanism outputs...")

write_csv(
  climate_thresholds,
  file.path(output_dir, "NASA_POWER_local_daily_extreme_thresholds.csv.gz"),
  na = ""
)
write_csv(
  grid_year_indices,
  file.path(output_dir, "soybean_climate_indices_grid_year.csv.gz"),
  na = ""
)
write_csv(
  country_year_climate,
  file.path(output_dir, "soybean_climate_indices_country_year.csv"),
  na = ""
)
write_csv(
  state_year_climate,
  file.path(output_dir, "soybean_climate_indices_state_province_year.csv.gz"),
  na = ""
)
write_csv(
  grid_phase_contrasts,
  file.path(output_dir, "ENSO_climate_phase_contrasts_grid.csv.gz"),
  na = ""
)
write_csv(
  country_phase_summary,
  file.path(output_dir, "ENSO_climate_phase_summary_country.csv"),
  na = ""
)
write_csv(
  state_phase_summary,
  file.path(output_dir, "ENSO_climate_phase_summary_state_province.csv"),
  na = ""
)
write_csv(
  compound_frequency_summary,
  file.path(output_dir, "ENSO_compound_event_frequency_by_phase.csv"),
  na = ""
)
write_csv(
  state_yield_climate,
  file.path(output_dir, "soybean_yield_climate_state_province_panel.csv.gz"),
  na = ""
)
write_csv(
  model_comparison,
  file.path(output_dir, "ENSO_climate_yield_model_comparison.csv"),
  na = ""
)
write_csv(
  daily_processing_qc,
  file.path(output_dir, "QC_climate_daily_processing.csv"),
  na = ""
)
write_csv(
  qc_summary,
  file.path(output_dir, "QC_ENSO_climate_mechanisms.csv"),
  na = ""
)
write_csv(
  methodology,
  file.path(output_dir, "README_ENSO_climate_mechanism_metrics.csv"),
  na = ""
)


# 14. FINAL VALIDATION ---------------------------------------------------------

if (qc_summary$value[qc_summary$test == "Missing primary climate predictors"] > 0) {
  stop("Primary climate predictors contain missing values. Inspect QC outputs.", call. = FALSE)
}
if (qc_summary$value[qc_summary$test == "Duplicated grid-year-window-metric rows"] > 0) {
  stop("Duplicated grid-year-window-metric rows were produced.", call. = FALSE)
}
invalid_models <- model_comparison |>
  filter(
    !is.finite(out_of_sample_r_squared) |
      validation_observations < 3 |
      validation_years < 3 |
      prediction_coverage < 0.95
  )
write_csv(
  invalid_models,
  file.path(output_dir, "QC_climate_yield_invalid_models.csv"),
  na = ""
)
if (nrow(invalid_models) > 0) {
  stop(
    "One or more climate-yield models lacked valid predictions for at least ",
    "95% of observations. Inspect QC_climate_yield_invalid_models.csv.",
    call. = FALSE
  )
}

message("")
message("Script 11 completed successfully.")
message("Climate grids analysed: ", length(grid_ids))
message("Complete harvest years: ", analysis_start_harvest_year, " to ", analysis_end_harvest_year)
message("Local climate metrics: ", length(metric_columns))
message("Compound event types: 4")
message("Output directory: ", normalizePath(output_dir, winslash = "/"))

# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 06b_spatially_resolved_ENSO_attribution.R
# VERSION: 3.0
# PURPOSE: Estimate ENSO predictive contribution without cancelling opposite
#          local responses. Statistical association, not causal attribution.
# ==============================================================================

# 1. SETTINGS ------------------------------------------------------------------
default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) project_dir_override else default_project_dir

analysis_start_year <- 1974L
analysis_end_year <- 2024L
minimum_unit_years <- 20L
cv_folds <- 5L
cv_gap_years <- 1L
reference_area_start_year <- 2015L
reference_area_end_year <- 2024L
bootstrap_repetitions <- 1999L
bootstrap_block_length <- 3L
bootstrap_seed <- 20260904L
run_circular_shift_test <- TRUE
circular_shift_repetitions <- 49L
circular_shift_seed <- 20260905L

# 2. PACKAGES ------------------------------------------------------------------
required_packages <- c("dplyr", "readr", "tidyr", "tibble", "purrr", "fixest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(tibble)
  library(purrr); library(fixest)
})

# 3. PATHS ---------------------------------------------------------------------
panel_file <- file.path(
  project_dir, "data_merged", "soybean_ENSO_spatial_panel_1974_2024.csv"
)
output_dir <- file.path(project_dir, "enso_spatial_attribution")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(panel_file)) stop("Required file not found: ", panel_file, call. = FALSE)

# 4. HELPERS -------------------------------------------------------------------
check_columns <- function(x, required) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  if (!any(ok) || sum(w[ok]) <= 0) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_sse <- function(y, pred, w) {
  ok <- is.finite(y) & is.finite(pred) & is.finite(w) & w >= 0
  if (!any(ok)) return(NA_real_)
  sum(w[ok] * (y[ok] - pred[ok])^2)
}

predictive_r2 <- function(y, pred, null_pred, w) {
  denominator <- weighted_sse(y, null_pred, w)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  1 - weighted_sse(y, pred, w) / denominator
}

fit_model_safe <- function(formula, data) {
  tryCatch(
    feols(formula, data = data, vcov = ~analysis_year, warn = FALSE,
          notes = FALSE),
    error = function(e) stop("Fixed effects model failed: ", conditionMessage(e),
                             call. = FALSE)
  )
}

model_specs <- tribble(
  ~model, ~family, ~structure,
  "common_symmetric", "symmetric ONI", "common",
  "state_specific_symmetric", "symmetric ONI", "state specific",
  "common_asymmetric", "asymmetric ONI", "common",
  "state_specific_asymmetric", "asymmetric ONI", "state specific"
)

# ENSO slopes vary among states or provinces. Unit fixed effects absorb every
# time invariant municipal or departmental difference. Standard errors are
# clustered by harvest year because ENSO is shared by all units in a year.
model_formula <- function(model_name, data) {
  if (n_distinct(data$unit_key[!is.na(data$unit_key)]) < 2) {
    stop("Fewer than two local units are available for the mixed model.",
         call. = FALSE)
  }
  rhs <- switch(
    model_name,
    null = "1",
    common_symmetric = "oni_djf_c",
    state_specific_symmetric = "0 + state_key:oni_djf_c",
    common_asymmetric = "el_nino_intensity + la_nina_intensity",
    state_specific_asymmetric = paste(
      "0 + state_key:el_nino_intensity + state_key:la_nina_intensity"
    ),
    stop("Unknown model: ", model_name, call. = FALSE)
  )
  as.formula(paste("yield_anomaly_percent ~", rhs, "| unit_key"))
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

neutral_data <- function(x) {
  x$oni_djf_c <- 0
  x$el_nino_intensity <- 0
  x$la_nina_intensity <- 0
  x
}

# Entire consecutive ENSO episodes are assigned to one fold. Allocation is
# balanced separately within El Nino, La Nina and neutral episodes.
create_event_folds <- function(year_table, k) {
  tagged <- year_table |>
    arrange(analysis_year) |>
    mutate(
      phase_simple = official_enso_phase,
      episode_id = cumsum(phase_simple != lag(phase_simple,
                                               default = first(phase_simple))) + 1L
    )
  episodes <- tagged |>
    count(episode_id, phase_simple, name = "episode_years") |>
    group_by(phase_simple) |>
    group_modify(~ {
      z <- arrange(.x, desc(episode_years), episode_id)
      load <- rep(0L, k); z$cv_fold <- NA_integer_
      for (i in seq_len(nrow(z))) {
        f <- which.min(load); z$cv_fold[i] <- f
        load[f] <- load[f] + z$episode_years[i]
      }
      z
    }) |>
    ungroup()
  tagged |>
    left_join(select(episodes, episode_id, cv_fold), by = "episode_id") |>
    select(analysis_year, phase_simple, episode_id, cv_fold)
}

# 5. DATA PREPARATION ----------------------------------------------------------
panel <- read_csv(panel_file, show_col_types = FALSE)
check_columns(panel, c(
  "country", "spatial_unit_id", "spatial_unit_name", "parent_region_code",
  "parent_region_name", "spatial_level", "analysis_year",
  "yield_anomaly_percent", "harvested_area_ha", "oni_djf_c"
))

panel <- panel |>
  filter(
    between(analysis_year, analysis_start_year, analysis_end_year),
    is.finite(yield_anomaly_percent), is.finite(oni_djf_c),
    is.finite(harvested_area_ha), harvested_area_ha > 0
  ) |>
  mutate(
    spatial_unit_id_original = as.character(spatial_unit_id),
    parent_region_code_original = as.character(parent_region_code),
    spatial_unit_id = normalize_spatial_id(country, spatial_unit_id),
    derived_region_code = substr(spatial_unit_id, 1, 2),
    parent_region_code = case_when(
      country %in% c("Brazil", "Argentina") &
        grepl("^[0-9]{2}$", derived_region_code) ~ derived_region_code,
      !is.na(parent_region_code_original) &
        nzchar(trimws(parent_region_code_original)) ~ parent_region_code_original,
      TRUE ~ parent_region_name
    ),
    state_key = paste(country, parent_region_code, sep = "__"),
    unit_key = paste(country, spatial_unit_id, sep = "__"),
    el_nino_intensity = pmax(oni_djf_c, 0),
    la_nina_intensity = pmax(-oni_djf_c, 0)
  )

# Recover the persistent NOAA event classification created by script 04.
# The continuous models still use ONI for every year, including neutral years.
if ("oni_harvest_year_phase" %in% names(panel)) {
  panel <- panel |>
    mutate(official_enso_phase = as.character(oni_harvest_year_phase))
} else {
  phase_candidates <- c(
    file.path(project_dir, "data_merged",
              "soybean_ENSO_country_year_1974_2024.csv"),
    file.path(project_dir, "enso_indices",
              "ENSO_indices_by_soybean_harvest_year_1974_2024.csv"),
    file.path(project_dir, "data_processed",
              "ENSO_indices_by_soybean_harvest_year_1974_2024.csv")
  )
  phase_file <- phase_candidates[file.exists(phase_candidates)][1]
  if (is.na(phase_file)) {
    stop(
      "Official ENSO phases were not found. Rerun scripts 04 and 05 first.",
      call. = FALSE
    )
  }
  phase_lookup <- read_csv(phase_file, show_col_types = FALSE)
  if ("country" %in% names(phase_lookup)) {
    phase_lookup <- phase_lookup |>
      filter(country %in% c("Brazil", "Argentina"))
  }
  phase_lookup <- phase_lookup |>
    transmute(
      analysis_year,
      official_enso_phase = as.character(oni_harvest_year_phase)
    ) |>
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
  ) |>
  group_by(country, unit_key) |>
  filter(n_distinct(analysis_year) >= minimum_unit_years) |>
  ungroup()

duplicates <- panel |> count(country, unit_key, analysis_year) |> filter(n > 1)
if (nrow(duplicates)) stop("Duplicated unit years were found.", call. = FALSE)

reference_area <- panel |>
  group_by(country, unit_key) |>
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

panel <- panel |>
  left_join(reference_area, by = c("country", "unit_key")) |>
  group_by(country, analysis_year) |>
  mutate(
    reference_area_coverage = sum(reference_weight_raw),
    weight_fixed = reference_weight_raw / sum(reference_weight_raw),
    weight_dynamic = harvested_area_ha / sum(harvested_area_ha)
  ) |>
  ungroup()

year_table <- panel |>
  group_by(analysis_year) |>
  summarise(
    oni_djf_c = first(oni_djf_c), oni_values = n_distinct(oni_djf_c),
    el_nino_intensity = first(el_nino_intensity),
    la_nina_intensity = first(la_nina_intensity),
    official_enso_phase = first(official_enso_phase),
    phase_values = n_distinct(official_enso_phase), .groups = "drop"
  )
if (any(year_table$oni_values != 1)) stop("ONI is not unique within year.")
if (any(year_table$phase_values != 1) || anyNA(year_table$official_enso_phase)) {
  stop("Official ENSO phase is missing or inconsistent within year.", call. = FALSE)
}
phase_counts <- count(year_table, official_enso_phase, name = "harvest_years")
if (!all(c("El Nino", "La Nina", "Neutral") %in% phase_counts$official_enso_phase)) {
  stop("The official ENSO classification does not contain all three phases.",
       call. = FALSE)
}
expected_phase_counts <- tibble(
  official_enso_phase = c("El Nino", "La Nina", "Neutral"),
  expected_harvest_years = c(17L, 19L, 15L)
)
phase_counts <- phase_counts |>
  left_join(expected_phase_counts, by = "official_enso_phase") |>
  mutate(matches_expected = harvest_years == expected_harvest_years)
if (!all(phase_counts$matches_expected)) {
  stop("Official phase counts differ from 17 El Nino, 19 La Nina and 15 neutral years.",
       call. = FALSE)
}
fold_table <- create_event_folds(year_table, cv_folds)
panel <- left_join(panel, fold_table, by = "analysis_year")

hierarchy_qc <- panel |>
  group_by(country) |>
  summarise(
    states_or_provinces = n_distinct(state_key),
    spatial_units = n_distinct(unit_key),
    missing_region_keys = sum(is.na(state_key) | grepl("__NA$", state_key)),
    .groups = "drop"
  )
if (any(hierarchy_qc$spatial_units < 2)) {
  stop("A country has fewer than two valid spatial units.", call. = FALSE)
}

# 6. EVENT BLOCKED OUT OF SAMPLE PREDICTIONS ----------------------------------
fit_fold <- function(country_data, fold, spec) {
  test_years <- unique(country_data$analysis_year[country_data$cv_fold == fold])
  excluded <- unique(unlist(lapply(test_years, function(y) {
    seq(y - cv_gap_years, y + cv_gap_years)
  })))
  train <- filter(country_data, !analysis_year %in% excluded)
  test <- filter(country_data, cv_fold == fold)
  training_units <- n_distinct(train$unit_key[!is.na(train$unit_key)])
  training_regions <- n_distinct(train$state_key[!is.na(train$state_key)])
  if (!nrow(train) || training_units < 2) {
    stop(
      "Cross-validation fold ", fold, " for ", first(country_data$country),
      " has fewer than two spatial units in training.", call. = FALSE
    )
  }
  null_model <- fit_model_safe(model_formula("null", train), train)
  enso_model <- fit_model_safe(model_formula(spec$model[[1]], train), train)
  pred_null <- as.numeric(predict(null_model, test))
  pred_enso <- as.numeric(predict(enso_model, test))
  pred_neutral <- as.numeric(predict(
    enso_model, neutral_data(test)
  ))
  predictions <- test |>
    mutate(
      model = spec$model[[1]], family = spec$family[[1]],
      structure = spec$structure[[1]],
      null_prediction_percent = pred_null,
      enso_prediction_percent = pred_enso,
      neutral_prediction_percent = pred_neutral,
      enso_component_percent = pred_enso - pred_neutral,
      squared_error_null = (yield_anomaly_percent - pred_null)^2,
      squared_error_enso = (yield_anomaly_percent - pred_enso)^2
    ) |>
    select(
      country, spatial_unit_id, spatial_unit_name, parent_region_code,
      parent_region_name, spatial_level, state_key, unit_key, analysis_year,
      phase_simple, episode_id, cv_fold, oni_djf_c, el_nino_intensity,
      la_nina_intensity, harvested_area_ha, reference_area_ha,
      reference_area_coverage, weight_fixed, weight_dynamic,
      yield_anomaly_percent, model, family, structure,
      null_prediction_percent, enso_prediction_percent,
      neutral_prediction_percent, enso_component_percent,
      squared_error_null, squared_error_enso
    )
  diagnostic <- tibble(
    country = first(test$country), model = spec$model[[1]], cv_fold = fold,
    test_episodes = n_distinct(test$episode_id),
    test_years = n_distinct(test$analysis_year),
    training_years = n_distinct(train$analysis_year),
    training_regions = training_regions,
    training_spatial_units = training_units,
    model_singular = FALSE,
    convergence_message = ""
  )
  list(predictions = predictions, diagnostic = diagnostic)
}

predictions <- list(); diagnostics <- list(); z <- 0L
for (country_name in sort(unique(panel$country))) {
  country_data <- filter(panel, country == country_name)
  message("Processing ", country_name, "...")
  for (j in seq_len(nrow(model_specs))) {
    for (fold in sort(unique(country_data$cv_fold))) {
      z <- z + 1L
      result <- fit_fold(country_data, fold, model_specs[j, ])
      predictions[[z]] <- result$predictions
      diagnostics[[z]] <- result$diagnostic
    }
  }
}
oof <- bind_rows(predictions)
fit_qc <- bind_rows(diagnostics)
if (any(!is.finite(oof$enso_prediction_percent))) stop("Missing OOF predictions.")

# 7. METRICS -------------------------------------------------------------------
summarise_metrics <- function(data, weight_name) {
  w <- data[[weight_name]]
  signed <- weighted_mean_safe(data$enso_component_percent, w)
  absolute <- weighted_mean_safe(abs(data$enso_component_percent), w)
  r2 <- predictive_r2(data$yield_anomaly_percent, data$enso_prediction_percent,
                      data$null_prediction_percent, w)
  observed_energy <- weighted_mean_safe(
    (data$yield_anomaly_percent - data$null_prediction_percent)^2, w
  )
  enso_signal_energy <- weighted_mean_safe(
    data$enso_component_percent^2, w
  )
  energy_fraction <- if_else(
    is.finite(observed_energy) & observed_energy > 0,
    enso_signal_energy / observed_energy, NA_real_
  )
  tibble(
    observations = nrow(data), years = n_distinct(data$analysis_year),
    spatial_units = if ("unit_key" %in% names(data))
      n_distinct(data$unit_key) else 1L,
    cv_r_squared = r2,
    cv_r_squared_percent = 100 * r2,
    cv_enso_signal_energy_fraction = energy_fraction,
    cv_enso_signal_energy_percent = 100 * energy_fraction,
    mean_signed_component_percent = signed,
    mean_absolute_component_percent = absolute,
    rms_component_percent = sqrt(weighted_mean_safe(
      data$enso_component_percent^2, w
    )),
    cancellation_index = if_else(
      is.finite(absolute) & absolute > 0, 1 - abs(signed) / absolute, NA_real_
    )
  )
}

sample_blocks <- function(n, block_length) {
  starts <- sample.int(n, ceiling(n / block_length), replace = TRUE)
  ids <- unlist(lapply(starts, function(s) {
    ((s - 1L + seq_len(block_length) - 1L) %% n) + 1L
  }))
  ids[seq_len(n)]
}

bootstrap_oof <- function(data, weight_name) {
  yearly <- data |>
    mutate(.w = .data[[weight_name]]) |>
    group_by(analysis_year) |>
    summarise(
      sse_model = sum(.w * squared_error_enso),
      sse_null = sum(.w * squared_error_null),
      signal_energy = sum(.w * enso_component_percent^2),
      absolute = sum(.w * abs(enso_component_percent)),
      weight = sum(.w), .groups = "drop"
    ) |>
    arrange(analysis_year)
  draws <- replicate(bootstrap_repetitions, {
    x <- yearly[sample_blocks(nrow(yearly), bootstrap_block_length), ]
    c(r2 = 1 - sum(x$sse_model) / sum(x$sse_null),
      energy = sum(x$signal_energy) / sum(x$sse_null),
      absolute = sum(x$absolute) / sum(x$weight))
  })
  tibble(
    cv_r2_ci_lower = quantile(draws["r2", ], 0.025, na.rm = TRUE),
    cv_r2_ci_upper = quantile(draws["r2", ], 0.975, na.rm = TRUE),
    energy_fraction_ci_lower = quantile(draws["energy", ], 0.025, na.rm = TRUE),
    energy_fraction_ci_upper = quantile(draws["energy", ], 0.975, na.rm = TRUE),
    absolute_ci_lower = quantile(draws["absolute", ], 0.025, na.rm = TRUE),
    absolute_ci_upper = quantile(draws["absolute", ], 0.975, na.rm = TRUE)
  )
}

set.seed(bootstrap_seed)
country_summary <- bind_rows(lapply(c("weight_fixed", "weight_dynamic"), function(w) {
  oof |>
    group_by(country, model, family, structure) |>
    group_modify(~ bind_cols(summarise_metrics(.x, w), bootstrap_oof(.x, w))) |>
    ungroup() |>
    mutate(
      weighting = if_else(w == "weight_fixed",
        "fixed 2015 to 2024 reference area", "annual harvested area sensitivity"),
      cv_r2_ci_lower_percent = 100 * cv_r2_ci_lower,
      cv_r2_ci_upper_percent = 100 * cv_r2_ci_upper,
      energy_fraction_ci_lower_percent = 100 * energy_fraction_ci_lower,
      energy_fraction_ci_upper_percent = 100 * energy_fraction_ci_upper
    )
}))

common <- country_summary |>
  filter(structure == "common") |>
  select(country, family, weighting, common_cv_r2_percent = cv_r_squared_percent)
country_summary <- country_summary |>
  left_join(common, by = c("country", "family", "weighting")) |>
  mutate(spatial_heterogeneity_gain_percentage_points =
           cv_r_squared_percent - common_cv_r2_percent)

phase_summary <- oof |>
  filter(model == "state_specific_asymmetric", phase_simple != "Neutral") |>
  group_by(country, phase_simple) |>
  group_modify(~ summarise_metrics(.x, "weight_fixed")) |>
  ungroup()

regional_summary <- oof |>
  filter(model == "state_specific_asymmetric") |>
  group_by(country, parent_region_code, parent_region_name, analysis_year) |>
  mutate(weight_region = reference_area_ha / sum(reference_area_ha)) |>
  ungroup() |>
  group_by(country, parent_region_code, parent_region_name) |>
  group_modify(~ summarise_metrics(.x, "weight_region")) |>
  ungroup()

local_summary <- oof |>
  filter(model == "state_specific_asymmetric") |>
  group_by(country, spatial_unit_id, spatial_unit_name, parent_region_code,
           parent_region_name, spatial_level, unit_key) |>
  group_modify(~ { x <- .x; x$weight_local <- 1; summarise_metrics(x, "weight_local") }) |>
  ungroup()

# 8. FULL PRIMARY MODEL AND STATE SPECIFIC COEFFICIENTS -----------------------
primary_spec <- filter(model_specs, model == "state_specific_asymmetric")
full_models <- list(); full_summary <- list(); coefficient_list <- list()
for (country_name in sort(unique(panel$country))) {
  d <- filter(panel, country == country_name)
  fit <- fit_model_safe(model_formula("state_specific_asymmetric", d), d)
  full_models[[country_name]] <- fit
  full_summary[[country_name]] <- tibble(
    country = country_name,
    model = "state_specific_asymmetric",
    singular_fit = FALSE,
    convergence_message = ""
  )
  metadata <- distinct(d, country, unit_key, spatial_unit_id, spatial_unit_name,
                       parent_region_code, parent_region_name, spatial_level,
                       state_key)
  states <- d |>
    group_by(state_key) |>
    slice(1) |>
    ungroup() |>
    select(state_key, unit_key)
  neutral_state <- states |>
    mutate(oni_djf_c = 0, el_nino_intensity = 0, la_nina_intensity = 0)
  el_state <- neutral_state |> mutate(oni_djf_c = 1, el_nino_intensity = 1)
  la_state <- neutral_state |> mutate(oni_djf_c = -1, la_nina_intensity = 1)
  state_effects <- states |>
    mutate(
      state_intercept_percent = as.numeric(predict(
        fit, neutral_state
      )),
      el_nino_effect_per_degree_percent = as.numeric(predict(
        fit, el_state
      )) - state_intercept_percent,
      la_nina_effect_per_degree_percent = as.numeric(predict(
        fit, la_state
      )) - state_intercept_percent
    ) |>
    select(-unit_key)
  unit_intercepts <- enframe(fixef(fit)$unit_key,
                             name = "unit_key", value = "unit_fixed_effect")
  coefficient_list[[country_name]] <- metadata |>
    left_join(state_effects, by = "state_key") |>
    left_join(unit_intercepts, by = "unit_key") |>
    mutate(intercept_percent = unit_fixed_effect) |>
    relocate(country, spatial_unit_id, spatial_unit_name, parent_region_code,
             parent_region_name, spatial_level, state_key, unit_key,
             intercept_percent, el_nino_effect_per_degree_percent,
             la_nina_effect_per_degree_percent)
}
full_summary <- bind_rows(full_summary)
local_coefficients <- bind_rows(coefficient_list)
state_coefficients <- local_coefficients |>
  distinct(country, parent_region_code, parent_region_name, state_key,
           state_intercept_percent, el_nino_effect_per_degree_percent,
           la_nina_effect_per_degree_percent)

# 9. CIRCULAR SHIFT TEST -------------------------------------------------------
circular_shift_summary <- tibble(); circular_shift_distribution <- tibble()
if (run_circular_shift_test) {
  set.seed(circular_shift_seed)
  possible <- seq_len(nrow(year_table) - 1L)
  shifts <- if (circular_shift_repetitions >= length(possible)) possible else
    sample(possible, circular_shift_repetitions)
  rows <- list(); z <- 0L
  for (country_name in sort(unique(panel$country))) {
    d <- filter(panel, country == country_name)
    null_fit <- fit_model_safe(model_formula("null", d), d)
    null_pred <- as.numeric(predict(null_fit, d))
    observed_pred <- as.numeric(predict(full_models[[country_name]], d))
    observed_r2 <- predictive_r2(d$yield_anomaly_percent, observed_pred,
                                 null_pred, d$weight_fixed)
    annual <- distinct(d, analysis_year, oni_djf_c) |> arrange(analysis_year)
    for (shift in shifts) {
      shifted <- annual$oni_djf_c[
        ((seq_len(nrow(annual)) - 1L - shift) %% nrow(annual)) + 1L
      ]
      lookup <- tibble(analysis_year = annual$analysis_year, shifted_oni = shifted)
      dp <- d |>
        select(-oni_djf_c, -el_nino_intensity, -la_nina_intensity) |>
        left_join(lookup, by = "analysis_year") |>
        mutate(
          oni_djf_c = shifted_oni,
          el_nino_intensity = pmax(oni_djf_c, 0),
          la_nina_intensity = pmax(-oni_djf_c, 0)
        ) |>
        select(-shifted_oni)
      perm_fit <- fit_model_safe(
        model_formula("state_specific_asymmetric", dp), dp
      )
      perm_pred <- as.numeric(predict(perm_fit, dp))
      z <- z + 1L
      rows[[z]] <- tibble(
        country = country_name, shift_years = shift,
        observed_r_squared = observed_r2,
        shifted_r_squared = predictive_r2(
          dp$yield_anomaly_percent, perm_pred, null_pred, dp$weight_fixed
        )
      )
    }
  }
  circular_shift_distribution <- bind_rows(rows)
  circular_shift_summary <- circular_shift_distribution |>
    group_by(country) |>
    summarise(
      observed_r_squared = first(observed_r_squared), shifts_tested = n(),
      randomization_p_value =
        (1 + sum(shifted_r_squared >= observed_r_squared)) / (n() + 1),
      null_median_r_squared = median(shifted_r_squared),
      null_95_percentile_r_squared = quantile(shifted_r_squared, 0.95),
      .groups = "drop"
    )
}

# 10. QUALITY CONTROL AND OUTPUTS ---------------------------------------------
weight_qc <- panel |>
  group_by(country, analysis_year) |>
  summarise(
    spatial_units = n_distinct(unit_key),
    reference_area_coverage = first(reference_area_coverage),
    fixed_weight_sum = sum(weight_fixed), dynamic_weight_sum = sum(weight_dynamic),
    .groups = "drop"
  )
qc <- tribble(
  ~check, ~value, ~status,
  "Retained observations", nrow(panel), "Information",
  "Retained spatial units", n_distinct(panel$unit_key), "Information",
  "Missing OOF predictions", sum(!is.finite(oof$enso_prediction_percent)),
  if_else(all(is.finite(oof$enso_prediction_percent)), "Pass", "Fail"),
  "Singular CV models", sum(fit_qc$model_singular), "Review",
  "CV convergence messages", sum(nzchar(fit_qc$convergence_message)), "Review",
  "Minimum reference area coverage", min(weight_qc$reference_area_coverage),
  "Information"
)
dictionary <- tribble(
  ~field, ~meaning,
  "cv_r_squared_percent", "Out of sample prediction error reduction versus null",
  "cv_enso_signal_energy_percent",
  "Cross validated ENSO signal energy relative to observed anomaly energy; this is a spatial footprint magnitude, not predictive R squared",
  "spatial_heterogeneity_gain_percentage_points",
  "Gain from allowing state and local responses to differ",
  "mean_absolute_component_percent", "Absolute local ENSO effect without cancellation",
  "mean_signed_component_percent", "Net area weighted ENSO effect",
  "cancellation_index", "Zero means no cancellation and one means complete cancellation"
)

write_csv(country_summary, file.path(output_dir, "ENSO_spatially_resolved_country_summary.csv"))
write_csv(phase_summary, file.path(output_dir, "ENSO_spatially_resolved_phase_summary.csv"))
write_csv(regional_summary, file.path(output_dir, "ENSO_spatially_resolved_regional_summary.csv"))
write_csv(local_summary, file.path(output_dir, "ENSO_spatially_resolved_local_summary.csv"))
write_csv(local_coefficients, file.path(output_dir, "ENSO_spatially_resolved_local_coefficients.csv"))
write_csv(state_coefficients,
          file.path(output_dir, "ENSO_state_province_ONI_coefficients.csv"))
write_csv(oof, file.path(output_dir, "ENSO_spatially_resolved_oof_predictions.csv"))
write_csv(full_summary,
          file.path(output_dir, "ENSO_state_specific_asymmetric_model_summary.csv"))
write_csv(circular_shift_summary, file.path(output_dir, "ENSO_circular_shift_test.csv"))
write_csv(circular_shift_distribution,
          file.path(output_dir, "ENSO_circular_shift_null_distribution.csv"))
write_csv(fold_table, file.path(output_dir, "QC_ENSO_event_folds.csv"))
write_csv(phase_counts, file.path(output_dir, "QC_ENSO_official_phase_counts.csv"))
write_csv(hierarchy_qc, file.path(output_dir, "QC_spatial_attribution_hierarchy.csv"))
write_csv(fit_qc, file.path(output_dir, "QC_spatial_attribution_model_fits.csv"))
write_csv(weight_qc, file.path(output_dir, "QC_spatial_attribution_weights.csv"))
write_csv(qc, file.path(output_dir, "QC_spatial_attribution_summary.csv"))
write_csv(dictionary, file.path(output_dir, "README_spatial_attribution_metrics.csv"))

message("\nSpatial attribution completed: ", output_dir)
print(
  country_summary |>
    filter(model == "state_specific_asymmetric",
           weighting == "fixed 2015 to 2024 reference area") |>
    select(country, cv_r_squared_percent, cv_r2_ci_lower_percent,
           cv_r2_ci_upper_percent, cv_enso_signal_energy_percent,
           energy_fraction_ci_lower_percent, energy_fraction_ci_upper_percent,
           spatial_heterogeneity_gain_percentage_points,
           mean_absolute_component_percent, cancellation_index),
  n = Inf
)
print(circular_shift_summary, n = Inf)

#!/usr/bin/env Rscript

# Figure 5 data preparation
# Spatial synthesis of reliable and predictable ENSO effects on productivity.

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
  file.path(project_dir, "outputs", "figure_data", "figure_5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

consistency_alpha <- 0.10

find_one <- function(root, candidates) {
  direct <- file.path(root, candidates)
  if (any(file.exists(direct))) {
    return(normalizePath(direct[file.exists(direct)][[1]]))
  }
  all_files <- list.files(root, recursive = TRUE, full.names = TRUE)
  hits <- all_files[basename(all_files) %in% candidates]
  if (length(hits) == 0L) {
    stop("Could not find any of: ", paste(candidates, collapse = ", "),
         call. = FALSE)
  }
  normalizePath(hits[order(nchar(hits), hits)][[1]], mustWork = TRUE)
}

check_required_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(label, " is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
}

standardize_numeric_code <- function(x, width) {
  value <- gsub("[^0-9]", "", as.character(x))
  value[value == ""] <- NA_character_
  ifelse(is.na(value), NA_character_,
         sprintf(paste0("%0", width, "d"), as.integer(value)))
}

standardize_unit_id <- function(country, x) {
  ifelse(country == "Brazil", standardize_numeric_code(x, 7L),
         standardize_numeric_code(x, 5L))
}

standardize_phase <- function(x) {
  value <- iconv(str_to_upper(str_squish(as.character(x))),
                 from = "", to = "ASCII//TRANSLIT")
  value <- gsub("[^A-Z]", "", value)
  case_when(
    value %in% c("ELNINO", "ELNIO") ~ "El Nino",
    value %in% c("LANINA", "LANIA") ~ "La Nina",
    TRUE ~ NA_character_
  )
}

safe_weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  weighted.mean(x[keep], w[keep])
}

exact_consistency_p <- function(successes, trials) {
  ifelse(
    is.finite(successes) & is.finite(trials) & trials > 0,
    pbinom(successes - 1L, size = trials, prob = 0.5, lower.tail = FALSE),
    NA_real_
  )
}

paths <- list(
  response = find_one(project_dir, c("Figure_1_local_phase_response.csv")),
  consistency = find_one(project_dir, c("Figure_1_crossvalidated_consistency.csv")),
  local_skill = find_one(project_dir, c("Figure_2_local_ENSO_metrics.csv")),
  cancellation = find_one(project_dir, c(
    "Figure_2_regional_cancellation.csv",
    "Figure_1_regional_cancellation.csv"
  )),
  driver = find_one(project_dir, c("Figure_3_local_dominant_driver_map.csv"))
)

response <- read_csv(
  paths$response,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
)
consistency <- read_csv(
  paths$consistency,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
)
local_skill <- read_csv(
  paths$local_skill,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
)
cancellation <- read_csv(paths$cancellation, show_col_types = FALSE)
driver <- read_csv(
  paths$driver,
  col_types = cols(spatial_unit_id = col_character(), .default = col_guess()),
  show_col_types = FALSE
)

check_required_columns(
  response,
  c("country", "spatial_unit_id", "spatial_unit_name", "parent_region_name",
    "macroregion", "phase", "mean_yield_anomaly_t_ha",
    "mean_yield_anomaly_percent", "events", "reference_area_ha", "longitude",
    "latitude"),
  "Local phase response"
)
check_required_columns(
  consistency,
  c("country", "spatial_unit_id", "phase", "evaluated_events",
    "consistent_events", "sign_consistency", "dominant_expected_sign"),
  "Cross validated consistency"
)
check_required_columns(
  local_skill,
  c("country", "spatial_unit_id", "cv_r_squared", "cv_r_squared_percent",
    "cv_enso_signal_energy_percent", "mean_absolute_component_percent",
    "reference_area_ha", "longitude", "latitude", "macroregion"),
  "Local ENSO metrics"
)
check_required_columns(
  cancellation,
  c("macroregion", "phase"),
  "Regional cancellation"
)
check_required_columns(
  driver,
  c("country", "spatial_unit_id", "driver_label", "mechanism_group",
    "robust_dominance", "evidence_class"),
  "Local dominant climate driver"
)

response_prepared <- response |>
  transmute(
    country = as.character(country),
    spatial_unit_id = standardize_unit_id(country, spatial_unit_id),
    spatial_unit_name,
    parent_region_name,
    macroregion,
    ENSO_phase = standardize_phase(phase),
    phase_events = as.integer(events),
    mean_productivity_anomaly_t_ha = as.numeric(mean_yield_anomaly_t_ha),
    mean_productivity_anomaly_percent = as.numeric(mean_yield_anomaly_percent),
    response_sign = case_when(
      mean_yield_anomaly_t_ha > 0 ~ "Gain",
      mean_yield_anomaly_t_ha < 0 ~ "Loss",
      TRUE ~ "No change"
    ),
    reference_area_ha = as.numeric(reference_area_ha),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  ) |>
  filter(!is.na(ENSO_phase))

consistency_prepared <- consistency |>
  transmute(
    country = as.character(country),
    spatial_unit_id = standardize_unit_id(country, spatial_unit_id),
    ENSO_phase = standardize_phase(phase),
    evaluated_events = as.integer(evaluated_events),
    consistent_events = as.integer(consistent_events),
    sign_consistency = as.numeric(sign_consistency),
    dominant_expected_sign,
    consistency_p_value = exact_consistency_p(
      consistent_events, evaluated_events
    ),
    statistically_consistent =
      consistency_p_value <= consistency_alpha & sign_consistency > 0.5
  ) |>
  filter(!is.na(ENSO_phase))

skill_prepared <- local_skill |>
  transmute(
    country = as.character(country),
    spatial_unit_id = standardize_unit_id(country, spatial_unit_id),
    cv_r_squared = as.numeric(cv_r_squared),
    cv_r_squared_percent = as.numeric(cv_r_squared_percent),
    ENSO_signal_energy_percent = as.numeric(cv_enso_signal_energy_percent),
    mean_absolute_ENSO_component_percent = as.numeric(mean_absolute_component_percent),
    positive_predictive_skill = is.finite(cv_r_squared) & cv_r_squared > 0
  )

driver_prepared <- driver |>
  transmute(
    country = as.character(country),
    spatial_unit_id = standardize_unit_id(country, spatial_unit_id),
    driver_label,
    mechanism_group,
    robust_driver = as.logical(robust_dominance),
    driver_evidence_class = evidence_class
  ) |>
  distinct(country, spatial_unit_id, .keep_all = TRUE)

local_synthesis <- response_prepared |>
  inner_join(
    consistency_prepared,
    by = c("country", "spatial_unit_id", "ENSO_phase")
  ) |>
  inner_join(skill_prepared, by = c("country", "spatial_unit_id")) |>
  left_join(driver_prepared, by = c("country", "spatial_unit_id")) |>
  mutate(
    evidence_class = case_when(
      statistically_consistent & positive_predictive_skill ~
        "Consistent and predictable",
      !statistically_consistent & positive_predictive_skill ~
        "Predictable but variable",
      statistically_consistent & !positive_predictive_skill ~
        "Consistent but low predictive skill",
      TRUE ~ "Uncertain local response"
    ),
    signed_evidence_class = paste(response_sign, evidence_class, sep = ": "),
    expected_sign_matches_mean = dominant_expected_sign == response_sign
  )

if (anyDuplicated(local_synthesis[c("country", "spatial_unit_id", "ENSO_phase")])) {
  stop("Duplicated local phase rows were detected after integration.",
       call. = FALSE)
}

expected_rows <- response_prepared |>
  count(country, spatial_unit_id, ENSO_phase)
join_qc <- response_prepared |>
  distinct(country, spatial_unit_id, ENSO_phase) |>
  left_join(
    local_synthesis |>
      transmute(country, spatial_unit_id, ENSO_phase, matched = TRUE),
    by = c("country", "spatial_unit_id", "ENSO_phase")
  ) |>
  summarise(
    response_rows = n(),
    matched_rows = sum(matched, na.rm = TRUE),
    match_fraction = mean(coalesce(matched, FALSE))
  )

if (join_qc$match_fraction < 0.95) {
  stop("Less than 95% of local phase responses matched skill and consistency data.",
       call. = FALSE)
}

regional_area_summary <- local_synthesis |>
  group_by(macroregion, ENSO_phase, evidence_class) |>
  summarise(
    spatial_units = n_distinct(spatial_unit_id),
    productive_area_ha = sum(reference_area_ha, na.rm = TRUE),
    .groups = "drop_last"
  ) |>
  mutate(
    total_productive_area_ha = sum(productive_area_ha),
    productive_area_fraction = productive_area_ha / total_productive_area_ha,
    productive_area_percent = 100 * productive_area_fraction
  ) |>
  ungroup()

regional_sign_summary <- local_synthesis |>
  group_by(macroregion, ENSO_phase, response_sign) |>
  summarise(
    spatial_units = n_distinct(spatial_unit_id),
    productive_area_ha = sum(reference_area_ha, na.rm = TRUE),
    .groups = "drop_last"
  ) |>
  mutate(
    total_productive_area_ha = sum(productive_area_ha),
    productive_area_percent = 100 * productive_area_ha / total_productive_area_ha
  ) |>
  ungroup()

regional_evidence_summary <- local_synthesis |>
  group_by(macroregion, ENSO_phase) |>
  summarise(
    spatial_units = n_distinct(spatial_unit_id),
    productive_area_ha = sum(reference_area_ha, na.rm = TRUE),
    area_weighted_productivity_anomaly_t_ha = safe_weighted_mean(
      mean_productivity_anomaly_t_ha, reference_area_ha
    ),
    area_weighted_productivity_anomaly_percent = safe_weighted_mean(
      mean_productivity_anomaly_percent, reference_area_ha
    ),
    area_weighted_sign_consistency_percent = 100 * safe_weighted_mean(
      sign_consistency, reference_area_ha
    ),
    area_with_statistically_consistent_signal_percent = 100 * sum(
      reference_area_ha[statistically_consistent], na.rm = TRUE
    ) / sum(reference_area_ha, na.rm = TRUE),
    area_with_positive_predictive_skill_percent = 100 * sum(
      reference_area_ha[positive_predictive_skill], na.rm = TRUE
    ) / sum(reference_area_ha, na.rm = TRUE),
    area_with_consistent_predictable_signal_percent = 100 * sum(
      reference_area_ha[statistically_consistent & positive_predictive_skill],
      na.rm = TRUE
    ) / sum(reference_area_ha, na.rm = TRUE),
    area_weighted_cv_r_squared_percent = safe_weighted_mean(
      cv_r_squared_percent, reference_area_ha
    ),
    area_weighted_ENSO_signal_energy_percent = safe_weighted_mean(
      ENSO_signal_energy_percent, reference_area_ha
    ),
    .groups = "drop"
  )

if (!"cancellation_mean_percent" %in% names(cancellation) &&
    "mean" %in% names(cancellation)) {
  cancellation <- cancellation |>
    rename(
      cancellation_mean_percent = mean,
      cancellation_lower_percent = lower,
      cancellation_upper_percent = upper
    )
}
check_required_columns(
  cancellation,
  c("cancellation_mean_percent", "cancellation_lower_percent",
    "cancellation_upper_percent", "n_events"),
  "Regional cancellation estimates"
)

cancellation_prepared <- cancellation |>
  transmute(
    macroregion,
    ENSO_phase = standardize_phase(phase),
    cancellation_mean_percent = as.numeric(cancellation_mean_percent),
    cancellation_lower_percent = as.numeric(cancellation_lower_percent),
    cancellation_upper_percent = as.numeric(cancellation_upper_percent),
    events = as.integer(n_events)
  ) |>
  filter(!is.na(ENSO_phase))

regional_evidence_summary <- regional_evidence_summary |>
  left_join(cancellation_prepared, by = c("macroregion", "ENSO_phase"))

class_order <- c(
  "Consistent and predictable",
  "Predictable but variable",
  "Consistent but low predictive skill",
  "Uncertain local response"
)

local_synthesis <- local_synthesis |>
  mutate(
    macroregion = factor(
      macroregion,
      levels = c(
        "Argentina", "Southern Brazil", "Central Brazil",
        "MATOPIBA and northeastern frontier"
      )
    ),
    ENSO_phase = factor(ENSO_phase, levels = c("El Nino", "La Nina")),
    evidence_class = factor(evidence_class, levels = class_order)
  ) |>
  arrange(ENSO_phase, macroregion, spatial_unit_id)

qc_class_counts <- local_synthesis |>
  count(macroregion, ENSO_phase, evidence_class, name = "spatial_units",
        .drop = FALSE)

qc_summary <- bind_rows(
  tibble(test = names(join_qc), value = as.character(unlist(join_qc[1, ]))),
  tribble(
    ~test, ~value,
    "Integrated local phase rows", as.character(nrow(local_synthesis)),
    "Unique productive spatial units", as.character(n_distinct(
      interaction(local_synthesis$country, local_synthesis$spatial_unit_id)
    )),
    "Macroregions", as.character(n_distinct(local_synthesis$macroregion)),
    "Phase specific classifications", as.character(n_distinct(
      interaction(local_synthesis$ENSO_phase, local_synthesis$evidence_class)
    )),
    "Mean agreement between expected and observed response sign",
      as.character(mean(local_synthesis$expected_sign_matches_mean, na.rm = TRUE)),
    "Missing coordinates", as.character(sum(
      !is.finite(local_synthesis$longitude) |
        !is.finite(local_synthesis$latitude)
    ))
  )
)

methods <- tribble(
  ~item, ~value,
  "figure_question", "Where are ENSO productivity effects locally reliable and predictively useful?",
  "spatial_unit", "Soybean productive municipality or department represented by its productive area",
  "response_sign", "Sign of the mean phase specific productivity anomaly in tonnes per hectare",
  "consistency", "Fraction of held out ENSO events matching the expected local response sign",
  "consistency_test", "One sided exact binomial test against a probability of 0.5",
  "consistency_alpha", as.character(consistency_alpha),
  "predictive_skill", "Positive local cross validated R squared from ENSO predictors",
  "consistent_and_predictable", "Sign consistency P value at or below 0.10 and positive cross validated R squared",
  "regional_weighting", "Reference soybean harvested area in hectares",
  "cancellation", "One minus the absolute regional mean response divided by the regional mean absolute local response",
  "interpretation", "Classes describe out of sample reliability and predictive association, not causality"
)

write_csv(
  local_synthesis,
  file.path(output_dir, "Figure_5_local_ENSO_evidence_classification.csv"), na = ""
)
write_csv(
  regional_area_summary,
  file.path(output_dir, "Figure_5_regional_evidence_area_summary.csv"), na = ""
)
write_csv(
  regional_sign_summary,
  file.path(output_dir, "Figure_5_regional_response_sign_area.csv"), na = ""
)
write_csv(
  regional_evidence_summary,
  file.path(output_dir, "Figure_5_regional_integrated_summary.csv"), na = ""
)
write_csv(
  cancellation_prepared,
  file.path(output_dir, "Figure_5_regional_cancellation.csv"), na = ""
)
write_csv(
  qc_class_counts,
  file.path(output_dir, "QC_Figure_5_class_counts.csv"), na = ""
)
write_csv(
  qc_summary,
  file.path(output_dir, "QC_Figure_5_data_integration.csv"), na = ""
)
write_csv(
  methods,
  file.path(output_dir, "README_Figure_5_methods.csv"), na = ""
)

message("Figure 5 data written to: ", normalizePath(output_dir))

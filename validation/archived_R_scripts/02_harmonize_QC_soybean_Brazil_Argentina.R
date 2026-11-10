# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 02_harmonize_QC_soybean_Brazil_Argentina.R
# VERSION: 2.0
# PURPOSE: Harmonize and quality control soybean data from Brazil and Argentina
# IMPORTANT: Source files are never overwritten
# ==============================================================================

# 1. SETTINGS ------------------------------------------------------------------

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")

project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) {
  project_dir_override
} else {
  default_project_dir
}

minimum_valid_years <- 20L
yield_review_lower_t_ha <- 0.05
yield_review_upper_t_ha <- 8.00
yield_exclusion_upper_t_ha <- 8.00
minimum_planted_area_coverage <- 0.99
yield_difference_absolute_t_ha <- 0.01
yield_difference_relative <- 0.01

# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c("dplyr", "readr", "tidyr")
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
})

# 3. PATHS ---------------------------------------------------------------------

input_dir <- file.path(project_dir, "data_processed")
output_dir <- file.path(project_dir, "data_qc")
metadata_dir <- file.path(project_dir, "metadata")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

brazil_file <- file.path(input_dir, "soybean_Brazil_municipality_historical.csv")
argentina_file <- file.path(input_dir, "soybean_Argentina_department_historical.csv")
if (!file.exists(brazil_file)) stop("Brazilian input file not found: ", brazil_file, call. = FALSE)
if (!file.exists(argentina_file)) stop("Argentinian input file not found: ", argentina_file, call. = FALSE)

# 4. HELPERS -------------------------------------------------------------------

check_required_columns <- function(data, required, dataset_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(dataset_name, " is missing required columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

safe_min_integer <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_integer_ else min(x)
}

safe_max_integer <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_integer_ else max(x)
}

safe_ratio <- function(numerator, denominator) {
  if_else(!is.na(denominator) & denominator > 0,
          numerator / denominator, NA_real_)
}

safe_sum_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}

# 5. READ DATA -----------------------------------------------------------------

message("Reading Brazil and Argentina data...")

brazil_raw <- read_csv(
  brazil_file,
  col_types = cols(
    country = col_character(), municipality_code = col_character(),
    municipality = col_character(), state = col_character(), year = col_integer(),
    planted_area_ha = col_double(), harvested_area_ha = col_double(),
    production_t = col_double(), yield_kg_ha_reported = col_double(),
    yield_t_ha = col_double(), yield_t_ha_calculated = col_double(),
    source = col_character(), spatial_level = col_character()
  ),
  na = c("", "NA", "N/A")
)

argentina_raw <- read_csv(
  argentina_file,
  col_types = cols(
    country = col_character(), province_code = col_character(),
    province = col_character(), department_code = col_character(),
    department = col_character(), campaign = col_character(),
    campaign_start_year = col_integer(), harvest_year = col_integer(),
    planted_area_ha = col_double(), harvested_area_ha = col_double(),
    production_t = col_double(), yield_kg_ha_reported = col_double(),
    yield_t_ha = col_double(), yield_t_ha_calculated = col_double(),
    source = col_character(), spatial_level = col_character()
  ),
  na = c("", "NA", "N/A")
)

check_required_columns(
  brazil_raw,
  c("municipality_code", "municipality", "state", "year", "planted_area_ha",
    "harvested_area_ha", "production_t", "yield_t_ha", "yield_t_ha_calculated"),
  "Brazilian dataset"
)
check_required_columns(
  argentina_raw,
  c("province_code", "province", "department_code", "department", "campaign",
    "harvest_year", "planted_area_ha", "harvested_area_ha", "production_t",
    "yield_t_ha", "yield_t_ha_calculated"),
  "Argentinian dataset"
)

# 6. DUPLICATES ----------------------------------------------------------------

brazil_duplicates <- brazil_raw |>
  count(municipality_code, year, name = "duplicate_count") |>
  filter(duplicate_count > 1)
argentina_duplicates <- argentina_raw |>
  count(department_code, harvest_year, name = "duplicate_count") |>
  filter(duplicate_count > 1)

write_csv(brazil_duplicates, file.path(output_dir, "QC_duplicates_Brazil.csv"), na = "")
write_csv(argentina_duplicates, file.path(output_dir, "QC_duplicates_Argentina.csv"), na = "")
if (nrow(brazil_duplicates) > 0 || nrow(argentina_duplicates) > 0) {
  stop("Duplicated spatial unit and year combinations were found. Review the duplicate files.", call. = FALSE)
}

# 7. HARMONIZATION -------------------------------------------------------------

brazil <- brazil_raw |>
  transmute(
    country = "Brazil", spatial_unit_id = municipality_code,
    spatial_unit_name = municipality, parent_region_code = NA_character_,
    parent_region_name = state, spatial_level = "municipality",
    season = as.character(year), analysis_year = year, planted_area_ha,
    harvested_area_ha, production_t, yield_reported_t_ha = yield_t_ha,
    yield_source_calculated_t_ha = yield_t_ha_calculated, source,
    source_row = row_number()
  )

argentina <- argentina_raw |>
  transmute(
    country = "Argentina", spatial_unit_id = department_code,
    spatial_unit_name = department, parent_region_code = province_code,
    parent_region_name = province, spatial_level = "department", season = campaign,
    analysis_year = harvest_year, planted_area_ha, harvested_area_ha, production_t,
    yield_reported_t_ha = yield_t_ha,
    yield_source_calculated_t_ha = yield_t_ha_calculated, source,
    source_row = row_number()
  )

common_start_year <- max(min(brazil$analysis_year, na.rm = TRUE),
                         min(argentina$analysis_year, na.rm = TRUE))
common_end_year <- min(max(brazil$analysis_year, na.rm = TRUE),
                       max(argentina$analysis_year, na.rm = TRUE))
common_period_tag <- paste0(common_start_year, "_", common_end_year)
message("Common period detected: ", common_start_year, " to ", common_end_year)

# 8. QUALITY CONTROL FLAGS -----------------------------------------------------

harmonized <- bind_rows(brazil, argentina) |>
  mutate(
    in_common_period = analysis_year >= common_start_year & analysis_year <= common_end_year,
    yield_harvested_t_ha = safe_ratio(production_t, harvested_area_ha),
    output_planted_t_ha = safe_ratio(production_t, planted_area_ha),
    unharvested_area_ha = if_else(!is.na(planted_area_ha) & !is.na(harvested_area_ha),
                                  planted_area_ha - harvested_area_ha, NA_real_),
    unharvested_area_fraction = safe_ratio(unharvested_area_ha, planted_area_ha),
    yield_absolute_difference_t_ha = abs(yield_reported_t_ha - yield_harvested_t_ha),
    yield_relative_difference = safe_ratio(yield_absolute_difference_t_ha,
                                           yield_reported_t_ha),
    flag_missing_production = is.na(production_t),
    flag_missing_harvested_area = is.na(harvested_area_ha),
    flag_negative_production = !is.na(production_t) & production_t < 0,
    flag_negative_harvested_area = !is.na(harvested_area_ha) & harvested_area_ha < 0,
    flag_negative_planted_area = !is.na(planted_area_ha) & planted_area_ha < 0,
    flag_production_without_harvested_area = !is.na(production_t) & production_t > 0 &
      (is.na(harvested_area_ha) | harvested_area_ha <= 0),
    flag_harvested_area_above_planted = !is.na(planted_area_ha) &
      !is.na(harvested_area_ha) & harvested_area_ha > planted_area_ha,
    flag_zero_harvest_complete_loss = !is.na(planted_area_ha) & planted_area_ha > 0 &
      !is.na(harvested_area_ha) & harvested_area_ha == 0 &
      !is.na(production_t) & production_t == 0,
    flag_zero_production_with_harvested_area = !is.na(harvested_area_ha) &
      harvested_area_ha > 0 & !is.na(production_t) & production_t == 0,
    flag_yield_difference = !is.na(yield_absolute_difference_t_ha) &
      yield_absolute_difference_t_ha > yield_difference_absolute_t_ha &
      !is.na(yield_relative_difference) &
      yield_relative_difference > yield_difference_relative,
    flag_yield_low_review = !is.na(yield_harvested_t_ha) &
      yield_harvested_t_ha < yield_review_lower_t_ha,
    flag_yield_high_review = !is.na(yield_harvested_t_ha) &
      yield_harvested_t_ha >= yield_review_upper_t_ha,
    flag_extreme_yield_exclude = !is.na(yield_harvested_t_ha) &
      yield_harvested_t_ha > yield_exclusion_upper_t_ha,
    flag_structural_error = flag_negative_production | flag_negative_harvested_area |
      flag_negative_planted_area | flag_production_without_harvested_area |
      flag_harvested_area_above_planted,
    flag_any_review = flag_structural_error | flag_zero_harvest_complete_loss |
      flag_zero_production_with_harvested_area | flag_yield_difference |
      flag_yield_low_review | flag_yield_high_review,
    valid_yield_observation = in_common_period & !flag_structural_error &
      !flag_extreme_yield_exclude & !is.na(harvested_area_ha) & harvested_area_ha > 0 &
      !is.na(production_t) & production_t >= 0 & !is.na(yield_harvested_t_ha),
    valid_production_observation = in_common_period & !flag_structural_error &
      !flag_extreme_yield_exclude & !is.na(production_t) & production_t >= 0,
    valid_annual_aggregation = valid_production_observation
  )

# 9. COVERAGE AND ELIGIBILITY --------------------------------------------------

unit_coverage <- harmonized |>
  filter(in_common_period) |>
  group_by(country, spatial_unit_id, spatial_unit_name, parent_region_code,
           parent_region_name, spatial_level) |>
  summarise(
    first_year_with_yield = safe_min_integer(analysis_year[valid_yield_observation]),
    last_year_with_yield = safe_max_integer(analysis_year[valid_yield_observation]),
    years_with_valid_yield = sum(valid_yield_observation, na.rm = TRUE),
    years_with_production_record = sum(valid_production_observation, na.rm = TRUE),
    years_with_planted_area = sum(!is.na(planted_area_ha) & planted_area_ha > 0,
                                  na.rm = TRUE),
    excluded_extreme_records = sum(flag_extreme_yield_exclude, na.rm = TRUE),
    total_production_t = sum(production_t[valid_annual_aggregation], na.rm = TRUE),
    mean_harvested_yield_t_ha = if_else(
      any(valid_yield_observation),
      mean(yield_harvested_t_ha[valid_yield_observation], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  ) |>
  mutate(eligible_minimum_years = years_with_valid_yield >= minimum_valid_years)

harmonized <- harmonized |>
  left_join(unit_coverage |>
              select(country, spatial_unit_id, eligible_minimum_years),
            by = c("country", "spatial_unit_id")) |>
  mutate(
    eligible_minimum_years = replace_na(eligible_minimum_years, FALSE),
    analysis_ready_yield = valid_yield_observation & eligible_minimum_years,
    analysis_ready_production = valid_production_observation & eligible_minimum_years
  )

# 10. COUNTRY AGGREGATES -------------------------------------------------------

country_annual <- harmonized |>
  filter(in_common_period) |>
  group_by(country, analysis_year) |>
  summarise(
    planted_area_observed_sum_ha = safe_sum_or_na(planted_area_ha[valid_annual_aggregation]),
    harvested_area_total_ha = sum(
      harvested_area_ha[valid_annual_aggregation],
      na.rm = TRUE
    ),
    production_total_t = sum(
      production_t[valid_annual_aggregation],
      na.rm = TRUE
    ),
    production_with_planted_area_t = sum(
      if_else(valid_annual_aggregation & !is.na(planted_area_ha), production_t, 0),
      na.rm = TRUE
    ),
    raw_harvested_area_ha = sum(harvested_area_ha, na.rm = TRUE),
    raw_production_t = sum(production_t, na.rm = TRUE),
    excluded_extreme_records = sum(flag_extreme_yield_exclude, na.rm = TRUE),
    excluded_extreme_production_t = sum(production_t[flag_extreme_yield_exclude],
                                        na.rm = TRUE),
    active_spatial_units = n_distinct(spatial_unit_id[
      valid_annual_aggregation & !is.na(production_t) & production_t > 0
    ]),
    complete_harvest_failures = sum(flag_zero_harvest_complete_loss, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    planted_area_data_coverage = safe_ratio(
      production_with_planted_area_t,
      production_total_t
    ),
    planted_area_coverage_complete = !is.na(planted_area_data_coverage) &
      planted_area_data_coverage >= minimum_planted_area_coverage,
    planted_area_ha = if_else(planted_area_coverage_complete,
                              planted_area_observed_sum_ha, NA_real_),
    yield_harvested_t_ha = safe_ratio(
      production_total_t,
      harvested_area_total_ha
    ),
    output_planted_t_ha = if_else(
      planted_area_coverage_complete,
      safe_ratio(production_total_t, planted_area_observed_sum_ha),
      NA_real_
    )
  ) |>
  transmute(
    country,
    analysis_year,
    planted_area_ha,
    planted_area_observed_sum_ha,
    harvested_area_ha = harvested_area_total_ha,
    production_t = production_total_t,
    yield_harvested_t_ha,
    output_planted_t_ha,
    production_with_planted_area_t,
    planted_area_data_coverage,
    planted_area_coverage_complete,
    active_spatial_units,
    complete_harvest_failures,
    excluded_extreme_records,
    excluded_extreme_production_t,
    raw_harvested_area_ha,
    raw_production_t
  )

# 11. COMBINED AGGREGATES ------------------------------------------------------

combined_annual <- country_annual |>
  group_by(analysis_year) |>
  summarise(
    planted_area_observed_sum_ha = safe_sum_or_na(planted_area_observed_sum_ha),
    harvested_area_ha = sum(harvested_area_ha, na.rm = TRUE),
    production_t = sum(production_t, na.rm = TRUE),
    production_with_planted_area_t = sum(production_with_planted_area_t, na.rm = TRUE),
    raw_harvested_area_ha = sum(raw_harvested_area_ha, na.rm = TRUE),
    raw_production_t = sum(raw_production_t, na.rm = TRUE),
    excluded_extreme_records = sum(excluded_extreme_records, na.rm = TRUE),
    excluded_extreme_production_t = sum(excluded_extreme_production_t, na.rm = TRUE),
    active_spatial_units = sum(active_spatial_units, na.rm = TRUE),
    complete_harvest_failures = sum(complete_harvest_failures, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    planted_area_data_coverage = safe_ratio(production_with_planted_area_t, production_t),
    planted_area_coverage_complete = !is.na(planted_area_data_coverage) &
      planted_area_data_coverage >= minimum_planted_area_coverage,
    planted_area_ha = if_else(planted_area_coverage_complete,
                              planted_area_observed_sum_ha, NA_real_),
    yield_harvested_t_ha = safe_ratio(production_t, harvested_area_ha),
    output_planted_t_ha = if_else(
      planted_area_coverage_complete,
      safe_ratio(production_t, planted_area_observed_sum_ha), NA_real_
    ),
    country = "Brazil and Argentina", .before = 1
  ) |>
  select(country, analysis_year, planted_area_ha, planted_area_observed_sum_ha,
         harvested_area_ha, production_t, yield_harvested_t_ha,
         output_planted_t_ha, production_with_planted_area_t,
         planted_area_data_coverage, planted_area_coverage_complete,
         active_spatial_units, complete_harvest_failures, excluded_extreme_records,
         excluded_extreme_production_t, raw_harvested_area_ha, raw_production_t)

# 12. QC SUMMARIES -------------------------------------------------------------

qc_summary <- harmonized |>
  group_by(country) |>
  summarise(
    first_year = min(analysis_year, na.rm = TRUE),
    last_year = max(analysis_year, na.rm = TRUE), total_rows = n(),
    spatial_units = n_distinct(spatial_unit_id), duplicate_keys = 0L,
    missing_production = sum(flag_missing_production, na.rm = TRUE),
    missing_harvested_area = sum(flag_missing_harvested_area, na.rm = TRUE),
    negative_production = sum(flag_negative_production, na.rm = TRUE),
    negative_harvested_area = sum(flag_negative_harvested_area, na.rm = TRUE),
    production_without_harvested_area = sum(flag_production_without_harvested_area,
                                             na.rm = TRUE),
    harvested_area_above_planted = sum(flag_harvested_area_above_planted,
                                        na.rm = TRUE),
    complete_harvest_failures = sum(flag_zero_harvest_complete_loss, na.rm = TRUE),
    zero_production_with_harvested_area = sum(flag_zero_production_with_harvested_area,
                                              na.rm = TRUE),
    yield_difference_flags = sum(flag_yield_difference, na.rm = TRUE),
    low_yield_review_flags = sum(flag_yield_low_review, na.rm = TRUE),
    high_yield_review_flags = sum(flag_yield_high_review, na.rm = TRUE),
    excluded_extreme_yield_records = sum(flag_extreme_yield_exclude, na.rm = TRUE),
    valid_yield_rows_common_period = sum(valid_yield_observation, na.rm = TRUE),
    eligible_spatial_units = n_distinct(spatial_unit_id[eligible_minimum_years]),
    analysis_ready_yield_rows = sum(analysis_ready_yield, na.rm = TRUE),
    analysis_ready_production_rows = sum(analysis_ready_production, na.rm = TRUE),
    .groups = "drop"
  )

qc_review_rows <- harmonized |>
  filter(flag_any_review) |>
  arrange(country, analysis_year, spatial_unit_id)
qc_excluded_extremes <- harmonized |>
  filter(flag_extreme_yield_exclude) |>
  arrange(country, analysis_year, spatial_unit_id)
planted_area_coverage_qc <- bind_rows(
  country_annual |> select(country, analysis_year, planted_area_data_coverage,
                           planted_area_coverage_complete, planted_area_ha,
                           planted_area_observed_sum_ha, production_t),
  combined_annual |> select(country, analysis_year, planted_area_data_coverage,
                            planted_area_coverage_complete, planted_area_ha,
                            planted_area_observed_sum_ha, production_t)
)

period_metadata <- tibble(
  parameter = c("common_start_year", "common_end_year", "minimum_valid_years",
                "yield_review_lower_t_ha", "yield_review_upper_t_ha",
                "yield_exclusion_upper_t_ha", "minimum_planted_area_coverage",
                "yield_difference_absolute_t_ha", "yield_difference_relative"),
  value = as.character(c(common_start_year, common_end_year, minimum_valid_years,
                         yield_review_lower_t_ha, yield_review_upper_t_ha,
                         yield_exclusion_upper_t_ha, minimum_planted_area_coverage,
                         yield_difference_absolute_t_ha, yield_difference_relative)),
  description = c(
    "First year available in both countries",
    "Last year available in both countries",
    "Minimum valid productivity observations required for regional analysis",
    "Productivity below this value is flagged for review",
    "Productivity at or above this value is flagged for review",
    "Productivity above this value is excluded from analysis ready outputs",
    "Minimum production weighted coverage required for planted area metrics",
    "Absolute reported versus calculated productivity difference threshold",
    "Relative reported versus calculated productivity difference threshold"
  )
)

# 13. INTERNAL VALIDATION ------------------------------------------------------

if (anyDuplicated(harmonized[c("country", "spatial_unit_id", "analysis_year")]) > 0) {
  stop("Duplicated keys remain after harmonization.", call. = FALSE)
}
# Planted area indicators must be unavailable whenever coverage is below the
# selected threshold. The validation is based on actual annual coverage because
# the first complete year can change when source data are updated.
if (any(
  !combined_annual$planted_area_coverage_complete &
    (!is.na(combined_annual$planted_area_ha) |
       !is.na(combined_annual$output_planted_t_ha)),
  na.rm = TRUE
)) {
  stop(
    "Combined planted area indicators were available despite incomplete coverage.",
    call. = FALSE
  )
}

if (any(
  combined_annual$planted_area_coverage_complete &
    (is.na(combined_annual$planted_area_ha) |
       is.na(combined_annual$output_planted_t_ha)),
  na.rm = TRUE
)) {
  stop(
    "Combined planted area indicators were missing despite complete coverage.",
    call. = FALSE
  )
}
country_sum_check <- country_annual |>
  group_by(analysis_year) |>
  summarise(expected_production_t = sum(production_t, na.rm = TRUE), .groups = "drop") |>
  left_join(combined_annual |>
              select(analysis_year, observed_production_t = production_t),
            by = "analysis_year") |>
  mutate(difference_t = observed_production_t - expected_production_t)
if (any(abs(country_sum_check$difference_t) > 1e-6, na.rm = TRUE)) {
  stop("Combined production does not equal the country sum.", call. = FALSE)
}

# 14. SAVE OUTPUTS -------------------------------------------------------------

message("Saving harmonized and quality control outputs...")
write_csv(harmonized,
          file.path(output_dir, "soybean_Brazil_Argentina_harmonized_all.csv"), na = "")
write_csv(filter(harmonized, analysis_ready_yield),
          file.path(output_dir, "soybean_Brazil_Argentina_analysis_ready_yield.csv"), na = "")
write_csv(filter(harmonized, analysis_ready_production),
          file.path(output_dir, "soybean_Brazil_Argentina_analysis_ready_production.csv"), na = "")
write_csv(unit_coverage, file.path(output_dir, "QC_spatial_unit_coverage.csv"), na = "")
write_csv(qc_review_rows, file.path(output_dir, "QC_observations_for_review.csv"), na = "")
write_csv(qc_excluded_extremes, file.path(output_dir, "QC_excluded_extreme_yield.csv"), na = "")
write_csv(planted_area_coverage_qc,
          file.path(output_dir, "QC_planted_area_coverage.csv"), na = "")
write_csv(country_annual,
          file.path(output_dir, paste0("soybean_country_annual_totals_",
                                      common_period_tag, ".csv")), na = "")
write_csv(combined_annual,
          file.path(output_dir, paste0("soybean_Brazil_Argentina_combined_annual_",
                                      common_period_tag, ".csv")), na = "")
write_csv(qc_summary, file.path(metadata_dir, "QC_summary_Brazil_Argentina.csv"), na = "")
write_csv(period_metadata, file.path(metadata_dir, "QC_analysis_parameters.csv"), na = "")
capture.output(sessionInfo(), file = file.path(metadata_dir, "QC_R_session_info.txt"))

# 15. FINAL REPORT -------------------------------------------------------------

message("")
message("Quality control and harmonization completed successfully.")
message("Common analysis period: ", common_start_year, " to ", common_end_year)
message("Minimum valid years: ", minimum_valid_years)
message("Extreme productivity records excluded: ",
        sum(harmonized$flag_extreme_yield_exclude, na.rm = TRUE))
first_complete_planted_year <- combined_annual |>
  filter(planted_area_coverage_complete) |>
  summarise(first_year = safe_min_integer(analysis_year)) |>
  pull(first_year)
message("First combined year with adequate planted area coverage: ",
        first_complete_planted_year)
message("Output directory: ", normalizePath(output_dir, winslash = "/"))
message("")
print(qc_summary)

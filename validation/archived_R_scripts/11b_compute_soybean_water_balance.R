# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 11b_compute_soybean_water_balance.R
# PURPOSE: Calculate FAO 56 ET0, stage specific ETc and daily soybean water balance
# VERSION: 1.3
# ==============================================================================

# 1. USER SETTINGS -------------------------------------------------------------

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) project_dir_override else default_project_dir

# NASA POWER solar radiation begins in 1984. Harvest year 1985 is the first
# crop season with complete Aug to Jul coverage for FAO 56 Penman Monteith ET0.
analysis_start_harvest_year <- 1985L
analysis_end_harvest_year <- 2024L
soil_water_capacities_mm <- c(75, 100, 125)
primary_capacity_mm <- 100
depletion_fraction <- 0.50
kc_initial <- 0.40
kc_midseason <- 1.15
kc_end <- 0.50
initial_stage_fraction <- 0.15
development_stage_fraction <- 0.25
midseason_stage_fraction <- 0.40
minimum_daily_coverage <- 0.95

# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c("dplyr", "readr", "tibble", "tidyr", "purrr", "lubridate")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) install.packages(missing_packages, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr)
  library(purrr); library(lubridate)
})

# 3. PATHS ---------------------------------------------------------------------

climate_dir <- file.path(project_dir, "climate_mechanisms")
weather_dir <- file.path(climate_dir, "NASA_POWER_daily_grid")
et0_dir <- file.path(climate_dir, "NASA_POWER_daily_ET0_inputs")
grid_manifest_file <- file.path(climate_dir, "NASA_POWER_grid_manifest.csv")
calendar_file <- file.path(climate_dir, "soybean_spatial_units_with_calendars.csv")
output_file <- file.path(climate_dir, "soybean_water_balance_spatial_unit_year.csv")
primary_file <- file.path(climate_dir, "soybean_water_balance_primary_AWC100.csv")
qc_file <- file.path(climate_dir, "QC_soybean_water_balance.csv")
metadata_file <- file.path(climate_dir, "README_soybean_water_balance_metrics.csv")

required_files <- c(grid_manifest_file, calendar_file)
if (any(!file.exists(required_files))) stop("Run scripts 09, 10 and 10b first.", call. = FALSE)

# 4. HELPERS -------------------------------------------------------------------

month_date <- function(harvest_year, month_value, day = 1L) {
  year_value <- ifelse(month_value >= 8L, harvest_year - 1L, harvest_year)
  as.Date(sprintf("%04d-%02d-%02d", year_value, month_value, day))
}

month_end_date <- function(harvest_year, month_value) {
  start <- month_date(harvest_year, month_value, 1L)
  ceiling_date(start, "month") - days(1)
}

window_midpoint_date <- function(harvest_year, start_month, end_month) {
  start <- month_date(harvest_year, start_month, 1L)
  end <- month_end_date(harvest_year, end_month)
  start + floor(as.numeric(end - start) / 2)
}

interpolate <- function(date, start, end, from, to) {
  if (end <= start) return(rep(to, length(date)))
  fraction <- pmin(1, pmax(0, as.numeric(date - start) / as.numeric(end - start)))
  from + fraction * (to - from)
}

calculate_kc <- function(dates) {
  n_days <- length(dates)
  relative_day <- seq_len(n_days)
  initial_end <- max(1L, round(n_days * initial_stage_fraction))
  development_end <- max(initial_end + 1L, round(
    n_days * (initial_stage_fraction + development_stage_fraction)
  ))
  midseason_end <- max(development_end + 1L, round(
    n_days * (
      initial_stage_fraction + development_stage_fraction +
        midseason_stage_fraction
    )
  ))
  development_end <- min(development_end, n_days)
  midseason_end <- min(midseason_end, n_days)

  kc <- rep(kc_initial, n_days)
  development <- relative_day > initial_end & relative_day <= development_end
  if (any(development)) {
    kc[development] <- seq(
      kc_initial, kc_midseason, length.out = sum(development) + 1L
    )[-1L]
  }
  middle <- relative_day > development_end & relative_day <= midseason_end
  kc[middle] <- kc_midseason
  late <- relative_day > midseason_end
  if (any(late)) {
    kc[late] <- seq(kc_midseason, kc_end, length.out = sum(late) + 1L)[-1L]
  }
  tibble(
    date = dates, kc = kc,
    crop_stage = case_when(
      relative_day <= initial_end ~ "initial",
      relative_day <= development_end ~ "vegetative",
      relative_day <= midseason_end ~ "reproductive",
      TRUE ~ "maturation"
    )
  )
}

extraterrestrial_radiation <- function(date, latitude_degree) {
  j <- yday(date)
  phi <- latitude_degree * pi / 180
  dr <- 1 + 0.033 * cos(2 * pi * j / 365)
  delta <- 0.409 * sin(2 * pi * j / 365 - 1.39)
  argument <- pmin(1, pmax(-1, -tan(phi) * tan(delta)))
  ws <- acos(argument)
  (24 * 60 / pi) * 0.0820 * dr *
    (ws * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(ws))
}

calculate_fao56_et0 <- function(x, latitude) {
  tmean <- x$T2M
  es_tmax <- 0.6108 * exp(17.27 * x$T2M_MAX / (x$T2M_MAX + 237.3))
  es_tmin <- 0.6108 * exp(17.27 * x$T2M_MIN / (x$T2M_MIN + 237.3))
  es <- (es_tmax + es_tmin) / 2
  ea <- 0.6108 * exp(17.27 * x$T2MDEW / (x$T2MDEW + 237.3))
  delta <- 4098 * (0.6108 * exp(17.27 * tmean / (tmean + 237.3))) /
    (tmean + 237.3)^2
  gamma <- 0.000665 * x$PS
  ra <- extraterrestrial_radiation(x$date, latitude)
  rso <- 0.75 * ra
  rs <- pmax(0, x$ALLSKY_SFC_SW_DWN)
  rns <- (1 - 0.23) * rs
  sigma <- 4.903e-9
  rnl <- sigma * (((x$T2M_MAX + 273.16)^4 + (x$T2M_MIN + 273.16)^4) / 2) *
    (0.34 - 0.14 * sqrt(pmax(0, ea))) *
    (1.35 * pmin(1, rs / pmax(rso, 0.01)) - 0.35)
  rn <- rns - rnl
  numerator <- 0.408 * delta * rn +
    gamma * (900 / (tmean + 273)) * x$WS2M * pmax(0, es - ea)
  denominator <- delta + gamma * (1 + 0.34 * x$WS2M)
  pmax(0, numerator / denominator)
}

run_bucket <- function(precipitation, etc, capacity, p) {
  n <- length(etc)
  storage <- eta <- ks <- drainage <- rep(NA_real_, n)
  previous_storage <- capacity
  for (i in seq_len(n)) {
    if (!is.finite(precipitation[i]) || !is.finite(etc[i])) next
    water_before_et <- min(capacity, previous_storage + max(0, precipitation[i]))
    drainage[i] <- max(0, previous_storage + max(0, precipitation[i]) - capacity)
    depletion <- capacity - water_before_et
    raw <- p * capacity
    ks[i] <- if (depletion <= raw) 1 else
      max(0, min(1, (capacity - depletion) / ((1 - p) * capacity)))
    eta[i] <- min(water_before_et, etc[i] * ks[i])
    storage[i] <- max(0, water_before_et - eta[i])
    previous_storage <- storage[i]
  }
  tibble(storage_mm = storage, eta_mm = eta, ks = ks, drainage_mm = drainage)
}

summarize_balance <- function(daily, capacity) {
  valid <- is.finite(daily$etc_mm) & is.finite(daily$eta_mm) &
    is.finite(daily$PRECTOTCORR)
  reproductive <- daily$crop_stage == "reproductive" & valid
  deficit <- pmax(0, daily$etc_mm - daily$eta_mm)
  stress <- valid & daily$ks < 0.999
  longest <- if (any(stress)) max(rle(stress)$lengths[rle(stress)$values]) else 0L
  tibble(
    awc_mm = capacity,
    daily_coverage = mean(valid),
    precipitation_mm = sum(daily$PRECTOTCORR[valid], na.rm = TRUE),
    et0_mm = sum(daily$et0_mm[valid], na.rm = TRUE),
    etc_mm = sum(daily$etc_mm[valid], na.rm = TRUE),
    eta_mm = sum(daily$eta_mm[valid], na.rm = TRUE),
    crop_aridity_index = sum(daily$PRECTOTCORR[valid], na.rm = TRUE) /
      sum(daily$etc_mm[valid], na.rm = TRUE),
    water_requirement_satisfaction = sum(daily$eta_mm[valid], na.rm = TRUE) /
      sum(daily$etc_mm[valid], na.rm = TRUE),
    cumulative_et_deficit_mm = sum(deficit[valid], na.rm = TRUE),
    reproductive_et_deficit_mm = sum(deficit[reproductive], na.rm = TRUE),
    reproductive_requirement_satisfaction = sum(daily$eta_mm[reproductive], na.rm = TRUE) /
      sum(daily$etc_mm[reproductive], na.rm = TRUE),
    water_stress_days = sum(stress),
    longest_water_stress_spell_days = longest,
    minimum_soil_water_storage_mm = min(daily$storage_mm[valid], na.rm = TRUE),
    minimum_relative_soil_water_storage = min(daily$storage_mm[valid] / capacity, na.rm = TRUE),
    drainage_mm = sum(daily$drainage_mm[valid], na.rm = TRUE),
    hot_water_stress_days = sum(stress & daily$T2M_MAX >= 30, na.rm = TRUE)
  )
}

# 5. INPUT TABLES --------------------------------------------------------------

grid_manifest <- read_csv(grid_manifest_file, show_col_types = FALSE) |>
  select(climate_location_id, grid_id, power_latitude, power_longitude) |>
  distinct()
calendars <- read_csv(calendar_file, show_col_types = FALSE) |>
  inner_join(grid_manifest, by = "climate_location_id")

calendar_columns <- c(
  "sowing_start_month", "sowing_end_month", "vegetative_start_month",
  "vegetative_end_month", "reproductive_start_month", "reproductive_end_month",
  "harvest_start_month", "harvest_end_month", "growing_season_start_month",
  "growing_season_end_month"
)
if (any(!calendar_columns %in% names(calendars))) stop("Crop calendar columns are missing", call. = FALSE)

analysis_groups <- calendars |>
  distinct(
    grid_id, power_latitude, power_longitude, calendar_group,
    across(all_of(calendar_columns))
  ) |>
  arrange(grid_id, calendar_group)

if (anyDuplicated(
  analysis_groups |> select(grid_id, calendar_group)
)) {
  stop(
    "A NASA POWER grid and calendar combination has multiple definitions.",
    call. = FALSE
  )
}

# 6. DAILY WATER BALANCE -------------------------------------------------------

results <- vector("list", nrow(analysis_groups) *
                         (analysis_end_harvest_year - analysis_start_harvest_year + 1L))
result_counter <- 0L

for (g in seq_len(nrow(analysis_groups))) {
  group <- analysis_groups[g, ]
  id <- group$grid_id[[1]]
  message("Water balance group ", g, " of ", nrow(analysis_groups), " [", id, "]...")
  base_path <- file.path(weather_dir, paste0(id, ".csv.gz"))
  et0_path <- file.path(et0_dir, paste0(id, "_ET0.csv.gz"))
  if (!file.exists(base_path) || !file.exists(et0_path)) stop("Missing weather for ", id, call. = FALSE)
  weather <- read_csv(base_path, show_col_types = FALSE, progress = FALSE) |>
    mutate(date = as.Date(date)) |>
    inner_join(
      read_csv(et0_path, show_col_types = FALSE, progress = FALSE) |>
        mutate(date = as.Date(date)) |>
        select(date, ALLSKY_SFC_SW_DWN, PS),
      by = "date"
    ) |>
    arrange(date)
  weather$et0_mm <- calculate_fao56_et0(weather, group$power_latitude[[1]])

  for (harvest_year in analysis_start_harvest_year:analysis_end_harvest_year) {
    # Regional windows describe spatial variation in farm operations. Their
    # midpoints define one representative crop cycle for each region.
    start_date <- window_midpoint_date(
      harvest_year,
      group$sowing_start_month[[1]], group$sowing_end_month[[1]]
    )
    end_date <- window_midpoint_date(
      harvest_year,
      group$harvest_start_month[[1]], group$harvest_end_month[[1]]
    )
    if (end_date <= start_date) {
      stop("Invalid representative crop cycle for ", group$calendar_group[[1]], call. = FALSE)
    }
    season <- weather |> filter(date >= start_date, date <= end_date)
    if (nrow(season) == 0) next
    stage <- calculate_kc(season$date)
    season <- bind_cols(season, stage |> select(-date)) |>
      mutate(etc_mm = kc * et0_mm)

    for (capacity in soil_water_capacities_mm) {
      daily <- bind_cols(
        season,
        run_bucket(season$PRECTOTCORR, season$etc_mm, capacity, depletion_fraction)
      )
      result_counter <- result_counter + 1L
      results[[result_counter]] <- summarize_balance(daily, capacity) |>
        mutate(
          grid_id = id, calendar_group = group$calendar_group[[1]],
          harvest_year = harvest_year, season_start_date = start_date,
          season_end_date = end_date, .before = 1
        )
    }
  }
}

grid_year <- bind_rows(results[seq_len(result_counter)])

# 7. ASSIGN GRID RESULTS TO SPATIAL UNITS -------------------------------------

unit_keys <- calendars |>
  select(
    climate_location_id, country, spatial_unit_id, spatial_unit_name,
    parent_region_code, parent_region_name, grid_id, calendar_group,
    reference_harvested_area_ha
  ) |>
  distinct()

unit_year <- unit_keys |>
  inner_join(grid_year, by = c("grid_id", "calendar_group")) |>
  arrange(country, spatial_unit_id, harvest_year, awc_mm)

duplicate_unit_year <- unit_year |>
  count(country, climate_location_id, harvest_year, awc_mm, name = "records") |>
  filter(records > 1L)

expected_rows <- n_distinct(unit_keys$climate_location_id) *
  (analysis_end_harvest_year - analysis_start_harvest_year + 1L) *
  length(soil_water_capacities_mm)

if (nrow(duplicate_unit_year) > 0) {
  write_csv(
    duplicate_unit_year,
    file.path(climate_dir, "QC_duplicated_soybean_water_balance_records.csv")
  )
  stop(
    nrow(duplicate_unit_year),
    " duplicated spatial unit, year and AWC combinations were produced.",
    call. = FALSE
  )
}
if (nrow(unit_year) != expected_rows) {
  stop(
    "Unexpected water balance row count. Expected ", expected_rows,
    " and produced ", nrow(unit_year), ".", call. = FALSE
  )
}

write_csv(unit_year, output_file, na = "")
write_csv(filter(unit_year, awc_mm == primary_capacity_mm), primary_file, na = "")

# 8. QUALITY CONTROL AND DOCUMENTATION ----------------------------------------

qc <- unit_year |>
  group_by(country, awc_mm) |>
  summarise(
    spatial_units = n_distinct(spatial_unit_id), years = n_distinct(harvest_year),
    rows = n(), minimum_daily_coverage = min(daily_coverage, na.rm = TRUE),
    mean_water_requirement_satisfaction = mean(water_requirement_satisfaction, na.rm = TRUE),
    mean_reproductive_deficit_mm = mean(reproductive_et_deficit_mm, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(qc, qc_file)

write_csv(tribble(
  ~metric, ~definition,
  "et0_mm", "FAO 56 Penman Monteith reference evapotranspiration",
  "etc_mm", "Reference evapotranspiration multiplied by the daily stage specific crop coefficient",
  "eta_mm", "Actual crop evapotranspiration constrained by soil water availability",
  "crop_aridity_index", "Crop season precipitation divided by crop evapotranspiration demand",
  "water_requirement_satisfaction", "Actual divided by potential crop evapotranspiration",
  "cumulative_et_deficit_mm", "Crop season sum of ETc minus ETa",
  "reproductive_et_deficit_mm", "ETc minus ETa accumulated during the reproductive window",
  "water_stress_days", "Days with FAO 56 water stress coefficient below one",
  "hot_water_stress_days", "Water stress days with maximum temperature at least 30 degrees Celsius",
  "awc_mm", "Maximum soil water storage. Primary value 100 mm, sensitivity values 75 and 125 mm",
  "representative_crop_cycle", "Midpoint of regional sowing window to midpoint of regional harvest window",
  "kc_stage_fractions", "Initial 15 percent, development 25 percent, midseason 40 percent, late season 20 percent"
), metadata_file)

if (any(unit_year$daily_coverage < minimum_daily_coverage, na.rm = TRUE)) {
  stop("At least one crop season has less than 95% daily coverage. Inspect ", qc_file, call. = FALSE)
}
crop_cycle_days <- as.numeric(
  as.Date(unit_year$season_end_date) - as.Date(unit_year$season_start_date)
) + 1
if (any(crop_cycle_days < 90 | crop_cycle_days > 190, na.rm = TRUE)) {
  stop("A representative crop cycle is outside 90 to 190 days.", call. = FALSE)
}
if (any(!is.finite(unit_year$water_requirement_satisfaction))) {
  stop("Non finite water balance metrics were produced.", call. = FALSE)
}
if (any(unit_year$water_requirement_satisfaction < 0 |
        unit_year$water_requirement_satisfaction > 1.000001, na.rm = TRUE)) {
  stop("Water requirement satisfaction is outside zero to one.", call. = FALSE)
}
message("Script 11b completed: ", nrow(unit_year), " spatial unit year capacity records.")

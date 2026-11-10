# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 10_download_NASA_POWER_soybean_regions.R
# PURPOSE: Download restartable daily NASA POWER weather data for the soybean
#          spatial sampling grid
# VERSION: 1.1
# ==============================================================================


# 1. USER SETTINGS -------------------------------------------------------------

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")

project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) project_dir_override else default_project_dir

# NASA POWER daily data begin on 1 January 1981. The August to December 1980
# portion of harvest year 1981 is unavailable, so complete crop-year climate
# summaries will begin with harvest year 1982.
power_first_date <- as.Date("1981-01-01")
power_last_date <- as.Date("2024-07-31")

# Current NASA POWER documentation describes a 0.5 degree global grid and asks
# users to submit only one request per grid cell. All analytical units assigned
# to the same cell therefore reuse one downloaded time series.
power_grid_resolution_degree <- 0.5

power_parameters <- c(
  "T2M", "T2M_MAX", "T2M_MIN", "T2MDEW", "RH2M",
  "PRECTOTCORR", "WS2M"
)

maximum_attempts <- 5L
request_timeout_seconds <- 300
pause_between_requests_seconds <- 0.35
checkpoint_every_requests <- 5L
minimum_daily_coverage <- 0.98

# Use Inf for the complete download. A small integer can be used for a test run.
maximum_grid_cells_this_run <- Inf


# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c(
  "dplyr", "readr", "tibble", "tidyr", "purrr", "stringr",
  "lubridate", "httr", "jsonlite"
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
  library(httr)
  library(jsonlite)
})


# 3. PATHS ---------------------------------------------------------------------

input_dir <- file.path(project_dir, "climate_mechanisms")
output_dir <- input_dir
weather_dir <- file.path(output_dir, "NASA_POWER_daily_grid")
dir.create(weather_dir, recursive = TRUE, showWarnings = FALSE)

manifest_file <- file.path(input_dir, "NASA_POWER_soybean_request_manifest.csv")
grid_manifest_file <- file.path(output_dir, "NASA_POWER_grid_manifest.csv")
file_index_file <- file.path(output_dir, "NASA_POWER_daily_file_index.csv")
status_file <- file.path(output_dir, "QC_NASA_POWER_download_status.csv")
variable_coverage_file <- file.path(output_dir, "QC_NASA_POWER_variable_coverage.csv")
year_coverage_file <- file.path(output_dir, "QC_NASA_POWER_year_coverage.csv")
metadata_file <- file.path(output_dir, "README_NASA_POWER_daily_data.csv")

if (!file.exists(manifest_file)) {
  stop("NASA POWER request manifest not found: ", manifest_file, call. = FALSE)
}


# 4. HELPERS -------------------------------------------------------------------

check_required_columns <- function(data, required, dataset_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(dataset_name, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

format_power_date <- function(x) format(as.Date(x), "%Y%m%d")

safe_file_name <- function(grid_id) paste0(grid_id, ".csv.gz")

expected_dates <- seq.Date(power_first_date, power_last_date, by = "day")
expected_daily_rows <- length(expected_dates)

read_cached_weather <- function(path) {
  tryCatch(
    read_csv(
      path,
      col_types = cols(
        grid_id = col_character(),
        power_latitude = col_double(),
        power_longitude = col_double(),
        date = col_date(),
        .default = col_double()
      ),
      progress = FALSE
    ),
    error = function(e) NULL
  )
}

validate_weather_table <- function(data) {
  if (is.null(data) || nrow(data) == 0) return(FALSE)
  required <- c(
    "grid_id", "power_latitude", "power_longitude", "date",
    power_parameters
  )
  if (!all(required %in% names(data))) return(FALSE)
  if (anyDuplicated(data$date)) return(FALSE)
  if (min(data$date, na.rm = TRUE) > power_first_date) return(FALSE)
  if (max(data$date, na.rm = TRUE) < power_last_date) return(FALSE)
  nrow(data) >= floor(minimum_daily_coverage * expected_daily_rows)
}

parse_power_json <- function(response, grid_id, latitude, longitude) {
  payload <- jsonlite::fromJSON(
    rawToChar(response$content),
    simplifyVector = TRUE
  )

  if (is.null(payload$properties$parameter)) {
    api_message <- paste(unlist(payload$messages), collapse = " | ")
    stop("NASA POWER returned no parameter data. ", api_message, call. = FALSE)
  }

  parameter_data <- payload$properties$parameter
  missing_parameters <- setdiff(power_parameters, names(parameter_data))
  if (length(missing_parameters) > 0) {
    stop(
      "NASA POWER response omitted: ",
      paste(missing_parameters, collapse = ", "), call. = FALSE
    )
  }

  all_date_keys <- sort(unique(unlist(lapply(parameter_data, names))))
  if (length(all_date_keys) == 0) {
    stop("NASA POWER response contained no daily date keys.", call. = FALSE)
  }

  result <- tibble(
    grid_id = grid_id,
    power_latitude = latitude,
    power_longitude = longitude,
    date = as.Date(all_date_keys, format = "%Y%m%d")
  )

  for (parameter in power_parameters) {
    values <- parameter_data[[parameter]]
    result[[parameter]] <- as.numeric(values[match(all_date_keys, names(values))])
    result[[parameter]][result[[parameter]] <= -999] <- NA_real_
  }

  result |>
    filter(date >= power_first_date, date <= power_last_date) |>
    arrange(date)
}

download_power_grid <- function(grid_id, latitude, longitude) {
  endpoint <- "https://power.larc.nasa.gov/api/temporal/daily/point"
  query <- list(
    parameters = paste(power_parameters, collapse = ","),
    community = "AG",
    longitude = format(longitude, scientific = FALSE, trim = TRUE),
    latitude = format(latitude, scientific = FALSE, trim = TRUE),
    start = format_power_date(power_first_date),
    end = format_power_date(power_last_date),
    format = "JSON",
    `time-standard` = "LST"
  )

  last_error <- "Unknown download error"

  for (attempt in seq_len(maximum_attempts)) {
    response <- tryCatch(
      httr::GET(
        endpoint,
        query = query,
        httr::timeout(request_timeout_seconds),
        httr::user_agent("Soybean ENSO climate research; NASA POWER API")
      ),
      error = function(e) e
    )

    if (inherits(response, "response") && httr::status_code(response) == 200L) {
      parsed <- tryCatch(
        parse_power_json(response, grid_id, latitude, longitude),
        error = function(e) e
      )
      if (is.data.frame(parsed) && validate_weather_table(parsed)) return(parsed)
      if (inherits(parsed, "error")) {
        last_error <- conditionMessage(parsed)
      } else {
        last_error <- "Downloaded table failed the completeness validation"
      }
    } else if (inherits(response, "response")) {
      last_error <- paste("HTTP", httr::status_code(response))
    } else if (inherits(response, "error")) {
      last_error <- conditionMessage(response)
    }

    if (attempt < maximum_attempts) {
      waiting_time <- min(60, 2^(attempt + 1)) + runif(1, 0, 1)
      message(
        "Attempt ", attempt, " failed for ", grid_id, ": ", last_error,
        ". Retrying in ", round(waiting_time, 1), " seconds..."
      )
      Sys.sleep(waiting_time)
    }
  }

  stop(last_error, call. = FALSE)
}

summarize_grid_file <- function(data, grid_id, file_path, source_status) {
  coverage <- vapply(
    power_parameters,
    function(parameter) mean(is.finite(data[[parameter]])),
    numeric(1)
  )
  tibble(
    grid_id = grid_id,
    file_name = basename(file_path),
    file_path = normalizePath(file_path, winslash = "/", mustWork = FALSE),
    status = source_status,
    rows = nrow(data),
    first_date = min(data$date, na.rm = TRUE),
    last_date = max(data$date, na.rm = TRUE),
    minimum_variable_coverage = min(coverage),
    downloaded_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    error_message = NA_character_
  )
}


# 5. CREATE THE UNIQUE NASA POWER GRID ----------------------------------------

message("Preparing unique 0.5 degree NASA POWER grid cells...")

request_manifest <- read_csv(
  manifest_file,
  col_types = cols(
    climate_location_id = col_character(),
    country = col_character(),
    spatial_unit_id = col_character(),
    parent_region_code = col_character(),
    latitude = col_double(),
    longitude = col_double(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

check_required_columns(
  request_manifest,
  c(
    "climate_location_id", "country", "spatial_unit_id", "spatial_unit_name",
    "parent_region_code", "parent_region_name", "latitude", "longitude",
    "reference_harvested_area_ha", "calendar_group"
  ),
  "NASA POWER soybean request manifest"
)

if (any(!is.finite(request_manifest$latitude)) ||
    any(!is.finite(request_manifest$longitude))) {
  stop("The request manifest contains invalid coordinates.", call. = FALSE)
}

grid_lookup <- request_manifest |>
  transmute(
    power_latitude = round(latitude / power_grid_resolution_degree) *
      power_grid_resolution_degree,
    power_longitude = round(longitude / power_grid_resolution_degree) *
      power_grid_resolution_degree
  ) |>
  distinct() |>
  arrange(power_latitude, power_longitude) |>
  mutate(grid_id = sprintf("POWER_GRID_%04d", row_number()))

grid_manifest <- request_manifest |>
  mutate(
    power_latitude = round(latitude / power_grid_resolution_degree) *
      power_grid_resolution_degree,
    power_longitude = round(longitude / power_grid_resolution_degree) *
      power_grid_resolution_degree
  ) |>
  left_join(
    grid_lookup,
    by = c("power_latitude", "power_longitude")
  ) |>
  mutate(
    weather_file = file.path(
      "NASA_POWER_daily_grid", safe_file_name(grid_id)
    ),
    centroid_to_grid_distance_degree = sqrt(
      (latitude - power_latitude)^2 +
        ((longitude - power_longitude) * cos(latitude * pi / 180))^2
    )
  ) |>
  arrange(country, climate_location_id)

write_csv(grid_manifest, grid_manifest_file)

grid_requests <- grid_manifest |>
  group_by(grid_id, power_latitude, power_longitude, weather_file) |>
  summarise(
    spatial_units = n(),
    countries = paste(sort(unique(country)), collapse = "; "),
    reference_harvested_area_ha = sum(reference_harvested_area_ha, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(grid_id)

if (anyDuplicated(grid_manifest$climate_location_id)) {
  stop("Duplicated analytical units were produced in the grid manifest.", call. = FALSE)
}
if (anyDuplicated(grid_requests$grid_id)) {
  stop("Duplicated NASA POWER grid identifiers were produced.", call. = FALSE)
}

message(
  nrow(request_manifest), " spatial units were assigned to ",
  nrow(grid_requests), " unique NASA POWER grid cells."
)


# 6. DOWNLOAD OR REUSE EACH GRID CELL -----------------------------------------

if (file.exists(status_file)) {
  previous_status <- read_csv(status_file, show_col_types = FALSE) |>
    mutate(
      grid_id = as.character(grid_id),
      file_name = as.character(file_name),
      file_path = as.character(file_path),
      status = as.character(status),
      rows = as.integer(rows),
      first_date = as.Date(first_date),
      last_date = as.Date(last_date),
      minimum_variable_coverage = as.numeric(minimum_variable_coverage),
      downloaded_at_utc = as.character(downloaded_at_utc),
      error_message = as.character(error_message)
    )
} else {
  previous_status <- tibble()
}

run_requests <- if (is.finite(maximum_grid_cells_this_run)) {
  grid_requests |>
    slice_head(n = min(nrow(grid_requests), as.integer(maximum_grid_cells_this_run)))
} else {
  grid_requests
}

status_rows <- vector("list", nrow(run_requests))

for (i in seq_len(nrow(run_requests))) {
  request <- run_requests[i, ]
  grid_id <- request$grid_id[[1]]
  file_path <- file.path(weather_dir, safe_file_name(grid_id))

  message(
    "NASA POWER grid ", i, " of ", nrow(run_requests),
    " [", grid_id, "]..."
  )

  cached <- if (file.exists(file_path)) read_cached_weather(file_path) else NULL

  if (validate_weather_table(cached)) {
    status_rows[[i]] <- summarize_grid_file(
      cached, grid_id, file_path, "cached_valid"
    )
  } else {
    if (file.exists(file_path)) {
      message("Existing cache is incomplete or corrupted and will be replaced.")
    }

    downloaded <- tryCatch(
      download_power_grid(
        grid_id,
        request$power_latitude[[1]],
        request$power_longitude[[1]]
      ),
      error = function(e) e
    )

    if (inherits(downloaded, "error")) {
      status_rows[[i]] <- tibble(
        grid_id = grid_id,
        file_name = basename(file_path),
        file_path = normalizePath(file_path, winslash = "/", mustWork = FALSE),
        status = "failed",
        rows = NA_integer_,
        first_date = as.Date(NA),
        last_date = as.Date(NA),
        minimum_variable_coverage = NA_real_,
        downloaded_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        error_message = conditionMessage(downloaded)
      )
    } else {
      write_csv(downloaded, file_path, na = "")
      status_rows[[i]] <- summarize_grid_file(
        downloaded, grid_id, file_path, "downloaded"
      )
      Sys.sleep(pause_between_requests_seconds)
    }
  }

  if (i %% checkpoint_every_requests == 0L || i == nrow(run_requests)) {
    current_status <- bind_rows(status_rows[seq_len(i)])
    combined_status <- bind_rows(
      previous_status |> filter(!grid_id %in% current_status$grid_id),
      current_status
    ) |>
      arrange(grid_id)
    write_csv(combined_status, status_file)
  }
}


# 7. FINAL CACHE AUDIT ---------------------------------------------------------

message("Auditing all expected NASA POWER cache files...")

audit_results <- map_dfr(seq_len(nrow(grid_requests)), function(i) {
  request <- grid_requests[i, ]
  grid_id <- request$grid_id[[1]]
  file_path <- file.path(weather_dir, safe_file_name(grid_id))
  data <- if (file.exists(file_path)) read_cached_weather(file_path) else NULL

  if (!validate_weather_table(data)) {
    return(tibble(
      grid_id = grid_id,
      file_name = basename(file_path),
      file_path = normalizePath(file_path, winslash = "/", mustWork = FALSE),
      status = "missing_or_invalid",
      rows = if (is.null(data)) NA_integer_ else nrow(data),
      first_date = as.Date(NA),
      last_date = as.Date(NA),
      minimum_variable_coverage = NA_real_,
      downloaded_at_utc = NA_character_,
      error_message = "File absent or failed validation"
    ))
  }

  summarize_grid_file(data, grid_id, file_path, "complete")
})

write_csv(audit_results, status_file)

complete_grid_ids <- audit_results$grid_id[audit_results$status == "complete"]

file_index <- grid_requests |>
  left_join(
    audit_results |>
      select(
        grid_id, status, rows, first_date, last_date,
        minimum_variable_coverage, error_message
      ),
    by = "grid_id"
  ) |>
  mutate(
    complete = status == "complete",
    absolute_file_path = file.path(
      normalizePath(weather_dir, winslash = "/", mustWork = FALSE),
      safe_file_name(grid_id)
    )
  ) |>
  arrange(grid_id)

write_csv(file_index, file_index_file)


# 8. VARIABLE AND YEAR COVERAGE -----------------------------------------------

coverage_data <- map_dfr(complete_grid_ids, function(grid_id) {
  data <- read_cached_weather(
    file.path(weather_dir, safe_file_name(grid_id))
  )

  variable_summary <- map_dfr(power_parameters, function(parameter) {
    tibble(
      grid_id = grid_id,
      parameter = parameter,
      valid_days = sum(is.finite(data[[parameter]])),
      total_days = nrow(data),
      coverage = mean(is.finite(data[[parameter]])),
      minimum = suppressWarnings(min(data[[parameter]], na.rm = TRUE)),
      maximum = suppressWarnings(max(data[[parameter]], na.rm = TRUE))
    )
  })

  variable_summary
})

if (nrow(coverage_data) > 0) {
  coverage_data <- coverage_data |>
    mutate(
      minimum = if_else(is.infinite(minimum), NA_real_, minimum),
      maximum = if_else(is.infinite(maximum), NA_real_, maximum)
    )
}
write_csv(coverage_data, variable_coverage_file)

year_coverage <- map_dfr(complete_grid_ids, function(grid_id) {
  this_grid_id <- grid_id
  data <- read_cached_weather(
    file.path(weather_dir, safe_file_name(grid_id))
  ) |>
    mutate(year = year(date))

  data |>
    group_by(year) |>
    summarise(
      observed_days = n(),
      expected_days_in_requested_period = sum(year(expected_dates) == first(year)),
      temperature_coverage = mean(is.finite(T2M)),
      precipitation_coverage = mean(is.finite(PRECTOTCORR)),
      .groups = "drop"
    ) |>
    mutate(grid_id = this_grid_id, .before = 1)
})
write_csv(year_coverage, year_coverage_file)


# 9. DATA DOCUMENTATION --------------------------------------------------------

metadata <- tribble(
  ~item, ~value,
  "source", "NASA POWER Daily API",
  "endpoint", "https://power.larc.nasa.gov/api/temporal/daily/point",
  "community", "AG",
  "time_standard", "Local Solar Time",
  "download_period", paste(power_first_date, "to", power_last_date),
  "first_complete_soybean_harvest_year", "1982",
  "grid_resolution_used", paste(power_grid_resolution_degree, "degree"),
  "parameters", paste(power_parameters, collapse = ", "),
  "temperature_unit", "degrees Celsius",
  "precipitation_unit", "millimetres per day",
  "relative_humidity_unit", "percent",
  "wind_speed_unit", "metres per second",
  "cache_format", "One compressed CSV file per unique climate grid cell",
  "api_documentation", "https://power.larc.nasa.gov/docs/services/api/temporal/daily/"
)
write_csv(metadata, metadata_file)


# 10. FINAL VALIDATION ---------------------------------------------------------

incomplete <- file_index |>
  filter(!complete)

message("")
message("NASA POWER download pass completed.")
message("Expected grid cells: ", nrow(grid_requests))
message("Complete grid cells: ", sum(file_index$complete, na.rm = TRUE))
message("Incomplete grid cells: ", nrow(incomplete))
message("Daily files directory: ", normalizePath(weather_dir, winslash = "/"))

if (nrow(incomplete) > 0) {
  stop(
    nrow(incomplete),
    " NASA POWER grid cells remain incomplete. Rerun the same script. Valid ",
    "cached files will be skipped and only missing cells will be downloaded. ",
    "Inspect QC_NASA_POWER_download_status.csv.",
    call. = FALSE
  )
}

if (any(coverage_data$coverage < minimum_daily_coverage, na.rm = TRUE)) {
  stop(
    "At least one variable has coverage below ",
    100 * minimum_daily_coverage,
    "%. Inspect QC_NASA_POWER_variable_coverage.csv.", call. = FALSE
  )
}

message("All NASA POWER grid cells passed validation.")
message("Script 10 completed successfully.")

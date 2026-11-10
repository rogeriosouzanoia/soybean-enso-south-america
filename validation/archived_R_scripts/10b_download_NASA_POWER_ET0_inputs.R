# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 10b_download_NASA_POWER_ET0_inputs.R
# PURPOSE: Download the additional NASA POWER variables required for FAO 56 ET0
# VERSION: 1.4
# ==============================================================================

# 1. USER SETTINGS -------------------------------------------------------------

default_project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
project_dir_override <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR")
project_dir <- if (nzchar(project_dir_override)) project_dir_override else default_project_dir

power_first_date <- as.Date("1981-01-01")
power_last_date <- as.Date("2024-07-31")
radiation_first_date <- as.Date("1984-01-01")
power_parameters <- c("ALLSKY_SFC_SW_DWN", "PS")
maximum_attempts <- 5L
request_timeout_seconds <- 300
pause_between_requests_seconds <- 0.35
checkpoint_every_requests <- 5L
minimum_daily_coverage <- 0.98
maximum_grid_cells_this_run <- Inf

# Four simultaneous requests reduce the expected runtime substantially while
# keeping the load on the NASA POWER service moderate. Reduce this value to 2
# if repeated HTTP 429 or 5xx responses occur.
number_parallel_workers <- 4L

# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c(
  "dplyr", "readr", "tibble", "tidyr", "purrr", "httr", "jsonlite",
  "future", "future.apply"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr); library(purrr)
  library(httr); library(jsonlite); library(future); library(future.apply)
})

# 3. PATHS ---------------------------------------------------------------------

climate_dir <- file.path(project_dir, "climate_mechanisms")
grid_manifest_file <- file.path(climate_dir, "NASA_POWER_grid_manifest.csv")
output_dir <- file.path(climate_dir, "NASA_POWER_daily_ET0_inputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(climate_dir, "QC_NASA_POWER_ET0_download_status.csv")
index_file <- file.path(climate_dir, "NASA_POWER_ET0_file_index.csv")
coverage_file <- file.path(climate_dir, "QC_NASA_POWER_ET0_variable_coverage.csv")
metadata_file <- file.path(climate_dir, "README_NASA_POWER_ET0_inputs.csv")

if (!file.exists(grid_manifest_file)) {
  stop("Grid manifest not found. Run script 10 first: ", grid_manifest_file, call. = FALSE)
}

# 4. HELPERS -------------------------------------------------------------------

expected_dates <- seq.Date(power_first_date, power_last_date, by = "day")
safe_file_name <- function(x) paste0(x, "_ET0.csv.gz")

read_cached <- function(path) {
  tryCatch(
    read_csv(path, col_types = cols(
      grid_id = col_character(), power_latitude = col_double(),
      power_longitude = col_double(), date = col_date(), .default = col_double()
    ), progress = FALSE),
    error = function(e) NULL
  )
}

valid_table <- function(x) {
  if (is.null(x) || nrow(x) == 0 || anyDuplicated(x$date)) return(FALSE)
  if (!all(c("grid_id", "date", power_parameters) %in% names(x))) return(FALSE)
  if (min(x$date, na.rm = TRUE) > power_first_date ||
      max(x$date, na.rm = TRUE) < power_last_date) return(FALSE)
  radiation_coverage <- mean(
    is.finite(x$ALLSKY_SFC_SW_DWN[x$date >= radiation_first_date])
  )
  pressure_coverage <- mean(is.finite(x$PS))
  nrow(x) >= floor(length(expected_dates) * minimum_daily_coverage) &&
    radiation_coverage >= minimum_daily_coverage &&
    pressure_coverage >= minimum_daily_coverage
}

parse_response <- function(response, grid_id, latitude, longitude) {
  # Keep the nested date keys intact. Some two variable POWER responses are
  # simplified to a data frame by jsonlite, which can remove the names needed
  # to align values with dates.
  payload <- fromJSON(rawToChar(response$content), simplifyVector = FALSE)
  parameters <- payload$properties$parameter
  if (is.null(parameters)) stop("NASA POWER returned no parameter data", call. = FALSE)
  missing <- setdiff(power_parameters, names(parameters))
  if (length(missing) > 0) stop("NASA POWER omitted: ", paste(missing, collapse = ", "), call. = FALSE)
  parameter_vectors <- lapply(parameters, function(values) {
    unlist(values, recursive = TRUE, use.names = TRUE)
  })
  keys <- sort(unique(unlist(lapply(parameter_vectors, names), use.names = FALSE)))
  if (length(keys) == 0) stop("NASA POWER response contained no daily date keys", call. = FALSE)
  result <- tibble(
    grid_id = grid_id, power_latitude = latitude, power_longitude = longitude,
    date = as.Date(keys, "%Y%m%d")
  )
  for (parameter in power_parameters) {
    values <- parameter_vectors[[parameter]]
    result[[parameter]] <- as.numeric(values[match(keys, names(values))])
    result[[parameter]][result[[parameter]] <= -999] <- NA_real_
  }
  result |> filter(date >= power_first_date, date <= power_last_date) |> arrange(date)
}

validation_diagnostic <- function(x) {
  if (is.null(x) || nrow(x) == 0) return("no rows")
  coverage <- c(
    ALLSKY_SFC_SW_DWN_available_period = mean(
      is.finite(x$ALLSKY_SFC_SW_DWN[x$date >= radiation_first_date])
    ),
    PS = mean(is.finite(x$PS))
  )
  paste0(
    "rows=", nrow(x),
    "; first_date=", min(x$date, na.rm = TRUE),
    "; last_date=", max(x$date, na.rm = TRUE),
    "; coverage=", paste(
      paste0(names(coverage), ":", sprintf("%.4f", coverage)),
      collapse = ", "
    )
  )
}

download_grid <- function(grid_id, latitude, longitude) {
  endpoint <- "https://power.larc.nasa.gov/api/temporal/daily/point"
  query <- list(
    parameters = paste(power_parameters, collapse = ","), community = "AG",
    longitude = longitude, latitude = latitude,
    start = format(power_first_date, "%Y%m%d"),
    end = format(power_last_date, "%Y%m%d"), format = "JSON",
    `time-standard` = "LST"
  )
  last_error <- "Unknown download error"
  for (attempt in seq_len(maximum_attempts)) {
    response <- tryCatch(
      GET(endpoint, query = query, timeout(request_timeout_seconds),
          user_agent("Soybean ENSO water balance research")),
      error = function(e) e
    )
    if (inherits(response, "response") && status_code(response) == 200L) {
      result <- tryCatch(parse_response(response, grid_id, latitude, longitude),
                         error = function(e) e)
      if (is.data.frame(result) && valid_table(result)) return(result)
      last_error <- if (inherits(result, "error")) conditionMessage(result) else
        paste("Downloaded table failed validation:", validation_diagnostic(result))
    } else if (inherits(response, "response")) {
      last_error <- paste("HTTP", status_code(response))
    } else {
      last_error <- conditionMessage(response)
    }
    if (attempt < maximum_attempts) Sys.sleep(min(60, 2^(attempt + 1)) + runif(1))
  }
  stop(last_error, call. = FALSE)
}

summarize_file <- function(x, id, path, status, error = NA_character_) {
  tibble(
    grid_id = id, file_name = basename(path), status = status,
    rows = if (is.null(x)) NA_integer_ else nrow(x),
    first_date = if (is.null(x)) as.Date(NA) else min(x$date, na.rm = TRUE),
    last_date = if (is.null(x)) as.Date(NA) else max(x$date, na.rm = TRUE),
    radiation_coverage = if (is.null(x)) NA_real_ else mean(
      is.finite(x$ALLSKY_SFC_SW_DWN[x$date >= radiation_first_date])
    ),
    pressure_coverage = if (is.null(x)) NA_real_ else mean(is.finite(x$PS)),
    error_message = error
  )
}

# 5. UNIQUE GRID REQUESTS -------------------------------------------------------

grid_requests <- read_csv(grid_manifest_file, show_col_types = FALSE) |>
  distinct(grid_id, power_latitude, power_longitude) |>
  arrange(grid_id)

required <- c("grid_id", "power_latitude", "power_longitude")
if (!all(required %in% names(grid_requests))) stop("Invalid grid manifest", call. = FALSE)

run_requests <- if (is.finite(maximum_grid_cells_this_run)) {
  slice_head(grid_requests, n = min(nrow(grid_requests), as.integer(maximum_grid_cells_this_run)))
} else grid_requests

# A synchronous preflight prevents hundreds of repeated requests when the API,
# a parameter name, or the local connection has a systematic problem. The
# downloaded preflight file is reused by the parallel stage.
preflight_candidates <- which(vapply(seq_len(nrow(run_requests)), function(i) {
  path <- file.path(output_dir, safe_file_name(run_requests$grid_id[[i]]))
  !valid_table(if (file.exists(path)) read_cached(path) else NULL)
}, logical(1)))

if (length(preflight_candidates) > 0) {
  preflight_index <- preflight_candidates[[1]]
  preflight_request <- run_requests[preflight_index, ]
  preflight_id <- preflight_request$grid_id[[1]]
  preflight_path <- file.path(output_dir, safe_file_name(preflight_id))
  message("Running synchronous NASA POWER preflight for ", preflight_id, "...")
  preflight_result <- tryCatch(
    download_grid(
      preflight_id,
      preflight_request$power_latitude[[1]],
      preflight_request$power_longitude[[1]]
    ),
    error = function(e) e
  )
  if (inherits(preflight_result, "error")) {
    stop(
      "NASA POWER preflight failed for ", preflight_id, ": ",
      conditionMessage(preflight_result), call. = FALSE
    )
  }
  write_csv(preflight_result, preflight_path, na = "")
  if (!valid_table(read_cached(preflight_path))) {
    stop("The NASA POWER preflight file failed validation after saving.", call. = FALSE)
  }
  message("Preflight succeeded. Starting parallel processing...")
}

# 6. PARALLEL AND RESTARTABLE DOWNLOAD -----------------------------------------

process_request <- function(request, request_number, total_requests) {
  id <- request$grid_id[[1]]
  path <- file.path(output_dir, safe_file_name(id))
  message(
    "ET0 inputs grid ", request_number, " of ", total_requests,
    " [", id, "]..."
  )
  cached <- if (file.exists(path)) read_cached(path) else NULL
  if (valid_table(cached)) {
    return(summarize_file(cached, id, path, "cached_valid"))
  }

  downloaded <- tryCatch(
    download_grid(id, request$power_latitude[[1]], request$power_longitude[[1]]),
    error = function(e) e
  )
  if (inherits(downloaded, "error")) {
    summarize_file(NULL, id, path, "failed", conditionMessage(downloaded))
  } else {
    write_csv(downloaded, path, na = "")
    Sys.sleep(pause_between_requests_seconds)
    summarize_file(downloaded, id, path, "downloaded")
  }
}

if (number_parallel_workers < 1L) stop("number_parallel_workers must be positive", call. = FALSE)
future::plan(future::multisession, workers = number_parallel_workers)

status_rows <- vector("list", nrow(run_requests))
batch_starts <- seq.int(1L, nrow(run_requests), by = number_parallel_workers)

for (batch_start in batch_starts) {
  batch_indices <- batch_start:min(
    nrow(run_requests), batch_start + number_parallel_workers - 1L
  )
  batch_results <- future.apply::future_lapply(
    batch_indices,
    function(i) process_request(run_requests[i, ], i, nrow(run_requests)),
    future.seed = TRUE,
    future.packages = c(
      "dplyr", "readr", "tibble", "tidyr", "purrr", "httr", "jsonlite"
    )
  )
  status_rows[batch_indices] <- batch_results

  completed_count <- max(batch_indices)
  if (completed_count %% checkpoint_every_requests <= number_parallel_workers ||
      completed_count == nrow(run_requests)) {
    write_csv(bind_rows(status_rows[seq_len(completed_count)]), status_file)
  }
}

future::plan(future::sequential)

# 7. COMPLETE AUDIT ------------------------------------------------------------

execution_status <- if (file.exists(status_file)) {
  read_csv(status_file, show_col_types = FALSE) |>
    select(grid_id, execution_status = status, execution_error = error_message)
} else {
  tibble(grid_id = character(), execution_status = character(), execution_error = character())
}

audit <- map_dfr(seq_len(nrow(grid_requests)), function(i) {
  id <- grid_requests$grid_id[[i]]
  path <- file.path(output_dir, safe_file_name(id))
  x <- if (file.exists(path)) read_cached(path) else NULL
  if (valid_table(x)) summarize_file(x, id, path, "complete") else
    summarize_file(x, id, path, "missing_or_invalid", "File absent or failed validation")
}) |>
  left_join(execution_status, by = "grid_id") |>
  mutate(
    error_message = if_else(
      status == "missing_or_invalid" & !is.na(execution_error),
      execution_error, error_message
    )
  ) |>
  select(-execution_status, -execution_error)
write_csv(audit, status_file)

file_index <- grid_requests |>
  left_join(audit, by = "grid_id") |>
  mutate(
    et0_weather_file = file.path("NASA_POWER_daily_ET0_inputs", safe_file_name(grid_id)),
    complete = status == "complete"
  )
write_csv(file_index, index_file)

coverage <- audit |>
  select(grid_id, radiation_coverage, pressure_coverage) |>
  pivot_longer(-grid_id, names_to = "variable", values_to = "coverage")
write_csv(coverage, coverage_file)

write_csv(tribble(
  ~field, ~value,
  "Purpose", "Additional inputs for daily FAO 56 Penman Monteith reference evapotranspiration",
  "Variables", "ALLSKY_SFC_SW_DWN and PS",
  "Radiation availability", paste(radiation_first_date, power_last_date, sep = " to "),
  "Original weather files modified", "No",
  "Period", paste(power_first_date, power_last_date, sep = " to "),
  "Time standard", "Local solar time"
), metadata_file)

if (any(!file_index$complete)) {
  stop(sum(!file_index$complete), " ET0 input files remain incomplete. Rerun the script.", call. = FALSE)
}
message("Script 10b completed: ", nrow(file_index), " grid cells are complete.")

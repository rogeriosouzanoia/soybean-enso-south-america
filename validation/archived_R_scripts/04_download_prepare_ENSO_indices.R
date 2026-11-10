# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT: 04_download_prepare_ENSO_indices.R
# PURPOSE: Download, harmonize and align official ENSO indices with soybean
#          harvest years in Brazil and Argentina
# VERSION: 1.1
# SOURCES: NOAA Physical Sciences Laboratory and NOAA Climate Prediction Center
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
download_timeout_seconds <- 180L
maximum_download_attempts <- 3L

# NOAA operational threshold for warm and cold ENSO conditions.
oni_phase_threshold_c <- 0.50
minimum_consecutive_oni_seasons <- 5L


# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c("dplyr", "readr", "tidyr", "tibble")

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
})


# 3. PATHS AND OFFICIAL SOURCES ------------------------------------------------

raw_dir <- file.path(project_dir, "climate_indices", "raw")
processed_dir <- file.path(project_dir, "climate_indices", "processed")
metadata_dir <- file.path(project_dir, "metadata")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

source_urls <- c(
  nino34_monthly = "https://psl.noaa.gov/data/correlation/nina34.anom.data",
  oni_seasonal = "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt",
  meiv2_bimonthly = "https://psl.noaa.gov/enso/mei/data/meiv2.data"
)

raw_files <- c(
  nino34_monthly = file.path(raw_dir, "NOAA_Nino34_monthly_raw.txt"),
  oni_seasonal = file.path(raw_dir, "NOAA_ONI_seasonal_raw.txt"),
  meiv2_bimonthly = file.path(raw_dir, "NOAA_MEIv2_bimonthly_raw.txt")
)


# 4. DOWNLOAD HELPERS ----------------------------------------------------------

download_with_retry <- function(url, destination, attempts = 3L) {
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = download_timeout_seconds)

  temporary_file <- paste0(destination, ".part")
  if (file.exists(temporary_file)) file.remove(temporary_file)

  last_error <- NULL

  for (attempt in seq_len(attempts)) {
    message("Downloading attempt ", attempt, " of ", attempts, ": ", url)

    result <- tryCatch({
      suppressWarnings(download.file(
        url = url,
        destfile = temporary_file,
        mode = "wb",
        method = "libcurl",
        quiet = TRUE
      ))
      TRUE
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      FALSE
    })

    valid_file <- result && file.exists(temporary_file) &&
      file.info(temporary_file)$size > 100

    if (valid_file) {
      if (file.exists(destination)) file.remove(destination)
      renamed <- file.rename(temporary_file, destination)
      if (!renamed) {
        stop("Could not finalize downloaded file: ", destination, call. = FALSE)
      }
      return(invisible(destination))
    }

    if (file.exists(temporary_file)) file.remove(temporary_file)
  }

  stop(
    "Download failed after ", attempts, " attempts: ", url,
    if (!is.null(last_error)) paste0(". Last error: ", last_error) else "",
    call. = FALSE
  )
}

for (source_name in names(source_urls)) {
  download_with_retry(
    source_urls[[source_name]],
    raw_files[[source_name]],
    maximum_download_attempts
  )
}


# 5. PARSING HELPERS -----------------------------------------------------------

parse_year_by_twelve_matrix <- function(file, value_name, source_name) {
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  tokens <- strsplit(trimws(lines), "[[:space:]]+")

  parsed <- lapply(tokens, function(x) {
    numeric_values <- suppressWarnings(as.numeric(x))
    valid <- length(numeric_values) >= 13 &&
      is.finite(numeric_values[1]) &&
      numeric_values[1] >= 1800 && numeric_values[1] <= 2100 &&
      all(is.finite(numeric_values[2:13]) | numeric_values[2:13] <= -90)

    if (!isTRUE(valid)) return(NULL)
    numeric_values[1:13]
  })

  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  if (length(parsed) == 0) {
    stop("No year by month records were parsed from ", file, call. = FALSE)
  }

  wide <- as_tibble(do.call(rbind, parsed))
  names(wide) <- c("year", paste0("month_", sprintf("%02d", 1:12)))

  wide |>
    mutate(year = as.integer(year)) |>
    pivot_longer(
      starts_with("month_"),
      names_to = "month_label",
      values_to = value_name
    ) |>
    mutate(
      month = as.integer(sub("month_", "", month_label)),
      date = as.Date(sprintf("%04d-%02d-01", year, month)),
      across(all_of(value_name), ~ replace(.x, .x <= -90, NA_real_)),
      source = source_name
    ) |>
    select(date, year, month, all_of(value_name), source)
}

parse_oni <- function(file) {
  valid_seasons <- c(
    "DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ",
    "JJA", "JAS", "ASO", "SON", "OND", "NDJ"
  )
  center_month <- setNames(1:12, valid_seasons)

  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  tokens <- strsplit(trimws(lines), "[[:space:]]+")

  parsed <- lapply(tokens, function(x) {
    if (length(x) < 4 || !(x[1] %in% valid_seasons)) return(NULL)
    year <- suppressWarnings(as.integer(x[2]))
    total <- suppressWarnings(as.numeric(x[3]))
    anomaly <- suppressWarnings(as.numeric(x[4]))
    if (!is.finite(year) || !is.finite(anomaly)) return(NULL)
    tibble(season = x[1], year = year, total_sst_c = total, oni_c = anomaly)
  })

  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  if (length(parsed) == 0) stop("No ONI records were parsed.", call. = FALSE)

  bind_rows(parsed) |>
    mutate(
      center_month = unname(center_month[season]),
      center_date = as.Date(sprintf("%04d-%02d-01", year, center_month)),
      source = "NOAA Climate Prediction Center ONI"
    ) |>
    arrange(center_date)
}

classify_oni_events <- function(oni_data) {
  phase_candidate <- case_when(
    oni_data$oni_c >= oni_phase_threshold_c ~ "El Nino",
    oni_data$oni_c <= -oni_phase_threshold_c ~ "La Nina",
    TRUE ~ "Neutral"
  )

  run_information <- rle(phase_candidate)
  run_id <- rep(seq_along(run_information$lengths), run_information$lengths)
  run_length <- rep(run_information$lengths, run_information$lengths)

  event_phase <- ifelse(
    phase_candidate != "Neutral" &
      run_length >= minimum_consecutive_oni_seasons,
    phase_candidate,
    "Neutral"
  )

  oni_data |>
    mutate(
      oni_phase_candidate = phase_candidate,
      oni_run_id = run_id,
      oni_run_length = run_length,
      oni_event_phase = event_phase
    ) |>
    group_by(oni_run_id) |>
    mutate(
      oni_event_peak_absolute_c = if_else(
        first(oni_event_phase) != "Neutral",
        max(abs(oni_c), na.rm = TRUE),
        NA_real_
      )
    ) |>
    ungroup() |>
    mutate(
      oni_event_strength = case_when(
        oni_event_phase == "Neutral" ~ "Neutral",
        oni_event_peak_absolute_c < 1.0 ~ "Weak",
        oni_event_peak_absolute_c < 1.5 ~ "Moderate",
        oni_event_peak_absolute_c < 2.0 ~ "Strong",
        oni_event_peak_absolute_c >= 2.0 ~ "Very strong",
        TRUE ~ NA_character_
      )
    )
}


# 6. PARSE OFFICIAL FILES ------------------------------------------------------

message("Parsing official ENSO indices...")

nino34 <- parse_year_by_twelve_matrix(
  raw_files[["nino34_monthly"]],
  "nino34_anomaly_c",
  "NOAA PSL Nino 3.4 monthly anomaly"
)

meiv2 <- parse_year_by_twelve_matrix(
  raw_files[["meiv2_bimonthly"]],
  "meiv2",
  "NOAA PSL MEI.v2"
) |>
  mutate(
    center_date = date,
    center_month = month,
    bimonthly_season = c(
      "DECJAN", "JANFEB", "FEBMAR", "MARAPR", "APRMAY", "MAYJUN",
      "JUNJUL", "JULAUG", "AUGSEP", "SEPOCT", "OCTNOV", "NOVDEC"
    )[center_month]
  ) |>
  select(
    center_date, date, year, month, center_month,
    bimonthly_season, meiv2, source
  )

oni <- parse_oni(raw_files[["oni_seasonal"]]) |>
  classify_oni_events()


# 7. MONTHLY HARMONIZED TABLE --------------------------------------------------

all_dates <- tibble(
  date = seq(
    from = min(c(nino34$date, oni$center_date, meiv2$center_date), na.rm = TRUE),
    to = max(c(nino34$date, oni$center_date, meiv2$center_date), na.rm = TRUE),
    by = "month"
  )
) |>
  mutate(
    year = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m"))
  )

enso_monthly <- all_dates |>
  left_join(
    nino34 |>
      select(date, nino34_anomaly_c),
    by = "date"
  ) |>
  left_join(
    oni |>
      select(
        center_date, oni_season = season, oni_c, oni_event_phase,
        oni_event_strength, oni_event_peak_absolute_c
      ),
    by = c("date" = "center_date")
  ) |>
  left_join(
    meiv2 |>
      select(center_date, meiv2_season = bimonthly_season, meiv2),
    by = c("date" = "center_date")
  ) |>
  arrange(date)


# 8. ALIGN INDICES WITH SOYBEAN HARVEST YEAR ----------------------------------

get_index_value <- function(data, value_column, target_year, target_month) {
  value <- data |>
    filter(
      .data$year == .env$target_year,
      .data$month == .env$target_month
    ) |>
    pull(all_of(value_column))
  if (length(value) == 0) NA_real_ else as.numeric(value[1])
}

get_oni_value <- function(
  data,
  target_year,
  target_center_month,
  value_column = "oni_c"
) {
  value <- data |>
    filter(
      .data$year == .env$target_year,
      .data$center_month == .env$target_center_month
    ) |>
    pull(all_of(value_column))
  if (length(value) == 0) {
    if (value_column %in% c("oni_event_phase", "oni_event_strength")) {
      return(NA_character_)
    }
    return(NA_real_)
  }
  value[1]
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

build_harvest_year_record <- function(harvest_year) {
  previous_year <- harvest_year - 1L

  nino_monthly <- c(
    aug_previous = get_index_value(nino34, "nino34_anomaly_c", previous_year, 8),
    sep_previous = get_index_value(nino34, "nino34_anomaly_c", previous_year, 9),
    oct_previous = get_index_value(nino34, "nino34_anomaly_c", previous_year, 10),
    nov_previous = get_index_value(nino34, "nino34_anomaly_c", previous_year, 11),
    dec_previous = get_index_value(nino34, "nino34_anomaly_c", previous_year, 12),
    jan = get_index_value(nino34, "nino34_anomaly_c", harvest_year, 1),
    feb = get_index_value(nino34, "nino34_anomaly_c", harvest_year, 2),
    mar = get_index_value(nino34, "nino34_anomaly_c", harvest_year, 3),
    apr = get_index_value(nino34, "nino34_anomaly_c", harvest_year, 4),
    may = get_index_value(nino34, "nino34_anomaly_c", harvest_year, 5)
  )

  mei_monthly <- c(
    aug_previous = get_index_value(meiv2, "meiv2", previous_year, 8),
    sep_previous = get_index_value(meiv2, "meiv2", previous_year, 9),
    oct_previous = get_index_value(meiv2, "meiv2", previous_year, 10),
    nov_previous = get_index_value(meiv2, "meiv2", previous_year, 11),
    dec_previous = get_index_value(meiv2, "meiv2", previous_year, 12),
    jan = get_index_value(meiv2, "meiv2", harvest_year, 1),
    feb = get_index_value(meiv2, "meiv2", harvest_year, 2),
    mar = get_index_value(meiv2, "meiv2", harvest_year, 3),
    apr = get_index_value(meiv2, "meiv2", harvest_year, 4),
    may = get_index_value(meiv2, "meiv2", harvest_year, 5)
  )

  tibble(
    harvest_year = harvest_year,

    nino34_aug_previous_c = nino_monthly["aug_previous"],
    nino34_sep_previous_c = nino_monthly["sep_previous"],
    nino34_oct_previous_c = nino_monthly["oct_previous"],
    nino34_nov_previous_c = nino_monthly["nov_previous"],
    nino34_dec_previous_c = nino_monthly["dec_previous"],
    nino34_jan_c = nino_monthly["jan"],
    nino34_feb_c = nino_monthly["feb"],
    nino34_mar_c = nino_monthly["mar"],
    nino34_apr_c = nino_monthly["apr"],
    nino34_may_c = nino_monthly["may"],

    nino34_son_previous_c = safe_mean(nino_monthly[c(
      "sep_previous", "oct_previous", "nov_previous"
    )]),
    nino34_ond_previous_c = safe_mean(nino_monthly[c(
      "oct_previous", "nov_previous", "dec_previous"
    )]),
    nino34_ndj_c = safe_mean(nino_monthly[c(
      "nov_previous", "dec_previous", "jan"
    )]),
    nino34_djf_c = safe_mean(nino_monthly[c(
      "dec_previous", "jan", "feb"
    )]),
    nino34_jfm_c = safe_mean(nino_monthly[c("jan", "feb", "mar")]),
    nino34_fma_c = safe_mean(nino_monthly[c("feb", "mar", "apr")]),
    nino34_mam_c = safe_mean(nino_monthly[c("mar", "apr", "may")]),
    nino34_growing_season_oct_apr_c = safe_mean(nino_monthly[c(
      "oct_previous", "nov_previous", "dec_previous",
      "jan", "feb", "mar", "apr"
    )]),
    nino34_reproductive_jan_mar_c = safe_mean(nino_monthly[c(
      "jan", "feb", "mar"
    )]),

    oni_son_previous_c = get_oni_value(oni, previous_year, 10),
    oni_ond_previous_c = get_oni_value(oni, previous_year, 11),
    oni_ndj_c = get_oni_value(oni, previous_year, 12),
    oni_djf_c = get_oni_value(oni, harvest_year, 1),
    oni_jfm_c = get_oni_value(oni, harvest_year, 2),
    oni_fma_c = get_oni_value(oni, harvest_year, 3),
    oni_mam_c = get_oni_value(oni, harvest_year, 4),
    oni_harvest_year_phase = get_oni_value(
      oni, harvest_year, 1, "oni_event_phase"
    ),
    oni_harvest_year_strength = get_oni_value(
      oni, harvest_year, 1, "oni_event_strength"
    ),

    meiv2_aug_previous = mei_monthly["aug_previous"],
    meiv2_sep_previous = mei_monthly["sep_previous"],
    meiv2_oct_previous = mei_monthly["oct_previous"],
    meiv2_nov_previous = mei_monthly["nov_previous"],
    meiv2_dec_previous = mei_monthly["dec_previous"],
    meiv2_jan = mei_monthly["jan"],
    meiv2_feb = mei_monthly["feb"],
    meiv2_mar = mei_monthly["mar"],
    meiv2_apr = mei_monthly["apr"],
    meiv2_may = mei_monthly["may"],
    meiv2_growing_season_oct_apr = safe_mean(mei_monthly[c(
      "oct_previous", "nov_previous", "dec_previous",
      "jan", "feb", "mar", "apr"
    )]),
    meiv2_reproductive_jan_mar = safe_mean(mei_monthly[c(
      "jan", "feb", "mar"
    )])
  )
}

enso_by_harvest_year <- bind_rows(lapply(
  analysis_start_year:analysis_end_year,
  build_harvest_year_record
))


# 9. QUALITY CONTROL -----------------------------------------------------------

if (anyDuplicated(nino34[c("year", "month")]) > 0) {
  stop("Duplicated Niño 3.4 year and month records were found.", call. = FALSE)
}
if (anyDuplicated(oni[c("year", "season")]) > 0) {
  stop("Duplicated ONI year and season records were found.", call. = FALSE)
}
if (anyDuplicated(meiv2[c("year", "center_month")]) > 0) {
  stop("Duplicated MEI.v2 year and season records were found.", call. = FALSE)
}
if (nrow(enso_by_harvest_year) != analysis_end_year - analysis_start_year + 1L) {
  stop("The harvest year table has an unexpected number of rows.", call. = FALSE)
}
if (any(is.na(enso_by_harvest_year$oni_djf_c))) {
  stop("ONI DJF is missing inside the soybean analysis period.", call. = FALSE)
}
if (any(is.na(enso_by_harvest_year$nino34_djf_c))) {
  stop("Niño 3.4 DJF is missing inside the soybean analysis period.", call. = FALSE)
}
if (n_distinct(enso_by_harvest_year$oni_djf_c) < 10) {
  stop(
    "ONI DJF has insufficient variation across soybean harvest years. ",
    "Check the year alignment.",
    call. = FALSE
  )
}
if (n_distinct(enso_by_harvest_year$nino34_djf_c) < 10) {
  stop(
    "Niño 3.4 DJF has insufficient variation across soybean harvest years. ",
    "Check the year alignment.",
    call. = FALSE
  )
}
if (n_distinct(enso_by_harvest_year$oni_harvest_year_phase) < 3) {
  stop(
    "The soybean harvest years do not contain all three ENSO phases. ",
    "Check the ONI alignment and classification.",
    call. = FALSE
  )
}

qc_index_summary <- bind_rows(
  tibble(
    index = "Nino 3.4 monthly",
    first_year = min(nino34$year[!is.na(nino34$nino34_anomaly_c)]),
    last_year = max(nino34$year[!is.na(nino34$nino34_anomaly_c)]),
    valid_values = sum(!is.na(nino34$nino34_anomaly_c)),
    missing_values = sum(is.na(nino34$nino34_anomaly_c)),
    minimum = min(nino34$nino34_anomaly_c, na.rm = TRUE),
    maximum = max(nino34$nino34_anomaly_c, na.rm = TRUE)
  ),
  tibble(
    index = "ONI",
    first_year = min(oni$year[!is.na(oni$oni_c)]),
    last_year = max(oni$year[!is.na(oni$oni_c)]),
    valid_values = sum(!is.na(oni$oni_c)),
    missing_values = sum(is.na(oni$oni_c)),
    minimum = min(oni$oni_c, na.rm = TRUE),
    maximum = max(oni$oni_c, na.rm = TRUE)
  ),
  tibble(
    index = "MEI.v2",
    first_year = min(meiv2$year[!is.na(meiv2$meiv2)]),
    last_year = max(meiv2$year[!is.na(meiv2$meiv2)]),
    valid_values = sum(!is.na(meiv2$meiv2)),
    missing_values = sum(is.na(meiv2$meiv2)),
    minimum = min(meiv2$meiv2, na.rm = TRUE),
    maximum = max(meiv2$meiv2, na.rm = TRUE)
  )
)

qc_correlations <- enso_by_harvest_year |>
  summarise(
    years_nino34_oni = sum(complete.cases(nino34_djf_c, oni_djf_c)),
    correlation_nino34_oni_djf = cor(
      nino34_djf_c, oni_djf_c, use = "complete.obs"
    ),
    years_nino34_meiv2 = sum(complete.cases(
      nino34_reproductive_jan_mar_c,
      meiv2_reproductive_jan_mar
    )),
    correlation_nino34_meiv2_reproductive = cor(
      nino34_reproductive_jan_mar_c,
      meiv2_reproductive_jan_mar,
      use = "complete.obs"
    )
  )

if (
  !is.finite(qc_correlations$correlation_nino34_oni_djf) ||
    qc_correlations$correlation_nino34_oni_djf < 0.70
) {
  stop(
    "The Niño 3.4 and ONI DJF correlation is unexpectedly low or undefined. ",
    "Check index parsing and harvest year alignment.",
    call. = FALSE
  )
}

if (
  !is.finite(qc_correlations$correlation_nino34_meiv2_reproductive) ||
    qc_correlations$correlation_nino34_meiv2_reproductive < 0.50
) {
  stop(
    "The Niño 3.4 and MEI.v2 correlation is unexpectedly low or undefined. ",
    "Check index parsing and harvest year alignment.",
    call. = FALSE
  )
}

enso_event_inventory <- enso_by_harvest_year |>
  count(
    oni_harvest_year_phase,
    oni_harvest_year_strength,
    name = "harvest_years"
  ) |>
  arrange(oni_harvest_year_phase, oni_harvest_year_strength)

source_metadata <- tibble(
  index = c("Nino 3.4", "ONI", "MEI.v2"),
  role_in_analysis = c(
    "Continuous monthly sensitivity index",
    "Primary continuous and categorical ENSO index",
    "Multivariate sensitivity index from 1979 onward"
  ),
  official_url = unname(source_urls),
  downloaded_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  local_raw_file = unname(raw_files)
)


# 10. SAVE OUTPUTS -------------------------------------------------------------

message("Saving ENSO indices and quality control outputs...")

write_csv(
  nino34,
  file.path(processed_dir, "ENSO_Nino34_monthly.csv"),
  na = ""
)
write_csv(
  oni,
  file.path(processed_dir, "ENSO_ONI_seasonal.csv"),
  na = ""
)
write_csv(
  meiv2,
  file.path(processed_dir, "ENSO_MEIv2_bimonthly.csv"),
  na = ""
)
write_csv(
  enso_monthly,
  file.path(processed_dir, "ENSO_indices_harmonized_monthly.csv"),
  na = ""
)
write_csv(
  enso_by_harvest_year,
  file.path(
    processed_dir,
    paste0(
      "ENSO_indices_by_soybean_harvest_year_",
      analysis_start_year,
      "_",
      analysis_end_year,
      ".csv"
    )
  ),
  na = ""
)
write_csv(
  qc_index_summary,
  file.path(processed_dir, "QC_ENSO_index_summary.csv"),
  na = ""
)
write_csv(
  qc_correlations,
  file.path(processed_dir, "QC_ENSO_index_correlations.csv"),
  na = ""
)
write_csv(
  enso_event_inventory,
  file.path(processed_dir, "QC_ENSO_event_inventory.csv"),
  na = ""
)
write_csv(
  source_metadata,
  file.path(metadata_dir, "ENSO_source_metadata.csv"),
  na = ""
)
capture.output(
  sessionInfo(),
  file = file.path(metadata_dir, "ENSO_R_session_info.txt")
)


# 11. FINAL REPORT -------------------------------------------------------------

message("")
message("ENSO index preparation completed successfully.")
message("Soybean harvest years: ", analysis_start_year, " to ", analysis_end_year)
message("Harvest year records: ", nrow(enso_by_harvest_year))
message("Processed directory: ", normalizePath(processed_dir, winslash = "/"))
message("")

print(qc_index_summary)
print(qc_correlations)
print(enso_event_inventory)

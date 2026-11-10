# ==============================================================================
# PROJECT: ENSO contribution to soybean variability in South America
# SCRIPT:  01_download_soybean_Brazil_Argentina.R
# PURPOSE: Download historical soybean data for Brazil and Argentina
# SOURCES: IBGE SIDRA PAM table 5457 and Argentina MAGyP
# ==============================================================================

# 1. USER SETTINGS -------------------------------------------------------------

# Change only this path if you want to store the project elsewhere.
project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")

# Minimum harvest year retained in the standardized files.
start_year <- 1974L

# NULL means that the script uses every year currently available.
end_year <- NULL

# Set TRUE to download the source files again.
# Set FALSE to reuse files already present in data_raw.
overwrite_raw <- FALSE


# 2. PACKAGES ------------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "httr2",
  "jsonlite",
  "purrr",
  "readr",
  "stringr",
  "tidyr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(httr2)
  library(jsonlite)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})


# 3. FOLDERS -------------------------------------------------------------------

raw_brazil_dir <- file.path(project_dir, "data_raw", "Brazil")
raw_argentina_dir <- file.path(project_dir, "data_raw", "Argentina")
processed_dir <- file.path(project_dir, "data_processed")
metadata_dir <- file.path(project_dir, "metadata")

walk(
  c(raw_brazil_dir, raw_argentina_dir, processed_dir, metadata_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)


# 4. HELPER FUNCTIONS ----------------------------------------------------------

download_text_with_retry <- function(url, attempts = 5L) {
  last_error <- NULL

  for (attempt in seq_len(attempts)) {
    result <- tryCatch(
      {
        request(url) |>
          req_user_agent("ENSO soybean research project") |>
          req_timeout(seconds = 180) |>
          req_retry(max_tries = 3) |>
          req_perform() |>
          resp_body_string()
      },
      error = function(e) e
    )

    if (!inherits(result, "error")) {
      return(result)
    }

    last_error <- result
    Sys.sleep(min(2^attempt, 20))
  }

  stop(
    "Download failed after ", attempts, " attempts. URL: ", url,
    "\nLast error: ", conditionMessage(last_error),
    call. = FALSE
  )
}

sidra_json_to_table <- function(json_text) {
  output <- jsonlite::fromJSON(json_text, simplifyDataFrame = TRUE)

  if (!is.data.frame(output) || nrow(output) < 2) {
    stop("The SIDRA response did not contain data.", call. = FALSE)
  }

  # The first JSON object contains the column descriptions.
  output |>
    slice(-1)
}

parse_sidra_number <- function(x) {
  # SIDRA symbols such as '-', '..', '...', and 'X' become NA.
  readr::parse_number(
    as.character(x),
    na = c("", "-", "..", "...", "X"),
    locale = locale(decimal_mark = ".", grouping_mark = "")
  )
}

campaign_harvest_year <- function(campaign, start) {
  end_text <- str_extract(campaign, "(?<=/)[0-9]{2,4}$")
  end_value <- suppressWarnings(as.integer(end_text))
  start <- suppressWarnings(as.integer(start))

  is_two_digit <- !is.na(end_text) & nchar(end_text) == 2
  century <- floor(start / 100) * 100
  end_from_two_digits <- century + end_value
  end_from_two_digits <- if_else(
    end_from_two_digits < start,
    end_from_two_digits + 100L,
    end_from_two_digits
  )

  if_else(is_two_digit, end_from_two_digits, end_value)
}


# 5. BRAZIL: IBGE SIDRA PAM ----------------------------------------------------

message("Downloading soybean data for Brazil...")

# SIDRA table 5457
# Product 40124: soybean grain
# Variable 8331: planted area or area intended for harvest, ha
# Variable 216: harvested area, ha
# Variable 214: production, tonnes
# Variable 112: average yield, kg per ha
sidra_table <- 5457L
sidra_product <- 40124L
sidra_variables <- c(8331L, 216L, 214L, 112L)

# Obtain the list of available years using a small national query.
years_url <- paste0(
  "https://apisidra.ibge.gov.br/values/t/", sidra_table,
  "/n1/1/v/214/p/all/c782/", sidra_product,
  "?formato=json"
)

years_json <- download_text_with_retry(years_url)
available_years <- sidra_json_to_table(years_json) |>
  transmute(year = suppressWarnings(as.integer(D3C))) |>
  filter(!is.na(year), year >= start_year) |>
  distinct(year) |>
  arrange(year) |>
  pull(year)

if (!is.null(end_year)) {
  available_years <- available_years[available_years <= end_year]
}

if (length(available_years) == 0) {
  stop("No Brazilian years were found for the selected period.", call. = FALSE)
}

download_brazil_year <- function(year) {
  json_file <- file.path(
    raw_brazil_dir,
    paste0("IBGE_PAM_soybean_municipal_", year, ".json")
  )

  csv_file <- file.path(
    raw_brazil_dir,
    paste0("IBGE_PAM_soybean_municipal_", year, ".csv")
  )

  url <- paste0(
    "https://apisidra.ibge.gov.br/values/t/", sidra_table,
    "/n6/all/v/", paste(sidra_variables, collapse = ","),
    "/p/", year,
    "/c782/", sidra_product,
    "?formato=json"
  )

  if (overwrite_raw || !file.exists(json_file)) {
    message("  Brazil: ", year)
    json_text <- download_text_with_retry(url)
    writeLines(json_text, json_file, useBytes = TRUE)
    Sys.sleep(0.25)
  } else {
    json_text <- paste(readLines(json_file, warn = FALSE), collapse = "\n")
  }

  annual_raw <- sidra_json_to_table(json_text)

  if (overwrite_raw || !file.exists(csv_file)) {
    write_csv(annual_raw, csv_file, na = "")
  }

  annual_wide <- annual_raw |>
    filter(str_detect(D1C, "^[0-9]+$"), str_detect(D2C, "^[0-9]+$")) |>
    transmute(
      municipality_code = as.character(D1C),
      municipality_full = as.character(D1N),
      year = as.integer(D3C),
      variable_code = as.character(D2C),
      value = parse_sidra_number(V)
    ) |>
    pivot_wider(
      id_cols = c(municipality_code, municipality_full, year),
      names_from = variable_code,
      values_from = value
    )

  # Some variables are unavailable in the earliest years. Create the expected
  # columns explicitly so that all annual tables retain the same structure.
  for (code in as.character(sidra_variables)) {
    if (!code %in% names(annual_wide)) {
      annual_wide[[code]] <- NA_real_
    }
  }

  annual_wide |>
    rename(
      planted_area_ha = `8331`,
      harvested_area_ha = `216`,
      production_t = `214`,
      yield_kg_ha_reported = `112`
    ) |>
    mutate(
      country = "Brazil",
      state = str_extract(municipality_full, "[A-Z]{2}$"),
      municipality = str_remove(municipality_full, "\\s+-\\s+[A-Z]{2}$"),
      yield_t_ha = yield_kg_ha_reported / 1000,
      yield_t_ha_calculated = if_else(
        harvested_area_ha > 0,
        production_t / harvested_area_ha,
        NA_real_
      ),
      source = "IBGE PAM table 5457",
      spatial_level = "municipality"
    ) |>
    select(
      country,
      municipality_code,
      municipality,
      state,
      year,
      planted_area_ha,
      harvested_area_ha,
      production_t,
      yield_kg_ha_reported,
      yield_t_ha,
      yield_t_ha_calculated,
      source,
      spatial_level
    )
}

brazil_soybean <- map_dfr(available_years, download_brazil_year) |>
  arrange(year, state, municipality_code)

brazil_output <- file.path(
  processed_dir,
  "soybean_Brazil_municipality_historical.csv"
)

write_csv(brazil_soybean, brazil_output, na = "")


# 6. ARGENTINA: MAGyP ----------------------------------------------------------

message("Downloading soybean data for Argentina...")

# The CKAN metadata endpoint is used to discover the current CSV address.
# This avoids hardcoding a filename containing the update month.
magyp_package_api <- paste0(
  "https://datos.magyp.gob.ar/api/3/action/package_show",
  "?id=estimaciones-agricolas"
)

magyp_metadata_text <- download_text_with_retry(magyp_package_api)
magyp_metadata <- jsonlite::fromJSON(
  magyp_metadata_text,
  simplifyDataFrame = TRUE
)

if (!isTRUE(magyp_metadata$success)) {
  stop("The MAGyP metadata service did not return a valid dataset.", call. = FALSE)
}

resources <- magyp_metadata$result$resources

resource_index <- which(
  str_to_upper(resources$format) == "CSV" &
    str_detect(
      str_to_lower(resources$name),
      fixed("estimaciones agrícolas")
    )
)

if (length(resource_index) == 0) {
  # Accent independent fallback.
  resource_index <- which(
    str_to_upper(resources$format) == "CSV" &
      str_detect(str_to_lower(resources$name), "estimaciones")
  )
}

if (length(resource_index) == 0) {
  stop("The official MAGyP agricultural estimates CSV was not found.", call. = FALSE)
}

resource_index <- resource_index[1]
argentina_url <- resources$url[resource_index]
argentina_source_name <- resources$fileName[resource_index]

if (is.na(argentina_source_name) || argentina_source_name == "") {
  argentina_source_name <- basename(argentina_url)
}

argentina_raw_file <- file.path(raw_argentina_dir, argentina_source_name)

if (overwrite_raw || !file.exists(argentina_raw_file)) {
  request(argentina_url) |>
    req_user_agent("ENSO soybean research project") |>
    req_timeout(seconds = 600) |>
    req_retry(max_tries = 5) |>
    req_perform(path = argentina_raw_file)
}

argentina_all <- read_csv(
  argentina_raw_file,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "s/d", "S/D")
)

# The source contains soja 1ra, soja 2da, and soja total.
# Only soja total is retained to prevent double counting.
argentina_soybean <- argentina_all |>
  mutate(crop_normalized = str_to_lower(str_squish(cultivo))) |>
  filter(crop_normalized == "soja total") |>
  mutate(
    campaign_start_year = as.integer(anio),
    harvest_year = campaign_harvest_year(campania, campaign_start_year),
    country = "Argentina",
    province_code = str_pad(as.character(provincia_id), width = 2, pad = "0"),
    department_code = str_pad(as.character(departamento_id), width = 5, pad = "0"),
    planted_area_ha = as.numeric(superficie_sembrada_ha),
    harvested_area_ha = as.numeric(superficie_cosechada_ha),
    production_t = as.numeric(produccion_tm),
    yield_kg_ha_reported = as.numeric(rendimiento_kgxha),
    yield_t_ha = yield_kg_ha_reported / 1000,
    yield_t_ha_calculated = if_else(
      harvested_area_ha > 0,
      production_t / harvested_area_ha,
      NA_real_
    ),
    source = "Argentina MAGyP Agricultural Estimates",
    spatial_level = "department"
  ) |>
  filter(!is.na(harvest_year), harvest_year >= start_year)

if (!is.null(end_year)) {
  argentina_soybean <- argentina_soybean |>
    filter(harvest_year <= end_year)
}

argentina_soybean <- argentina_soybean |>
  transmute(
    country,
    province_code,
    province = provincia,
    department_code,
    department = departamento,
    campaign = campania,
    campaign_start_year,
    harvest_year,
    planted_area_ha,
    harvested_area_ha,
    production_t,
    yield_kg_ha_reported,
    yield_t_ha,
    yield_t_ha_calculated,
    source,
    spatial_level
  ) |>
  arrange(harvest_year, province_code, department_code)

argentina_output <- file.path(
  processed_dir,
  "soybean_Argentina_department_historical.csv"
)

write_csv(argentina_soybean, argentina_output, na = "")


# 7. DOWNLOAD AND DATA DIAGNOSTICS ---------------------------------------------

diagnostic <- bind_rows(
  brazil_soybean |>
    summarise(
      country = "Brazil",
      first_harvest_year = min(year, na.rm = TRUE),
      last_harvest_year = max(year, na.rm = TRUE),
      observations = n(),
      spatial_units = n_distinct(municipality_code),
      missing_harvested_area = sum(is.na(harvested_area_ha)),
      missing_production = sum(is.na(production_t)),
      missing_yield = sum(is.na(yield_t_ha))
    ),
  argentina_soybean |>
    summarise(
      country = "Argentina",
      first_harvest_year = min(harvest_year, na.rm = TRUE),
      last_harvest_year = max(harvest_year, na.rm = TRUE),
      observations = n(),
      spatial_units = n_distinct(department_code),
      missing_harvested_area = sum(is.na(harvested_area_ha)),
      missing_production = sum(is.na(production_t)),
      missing_yield = sum(is.na(yield_t_ha))
    )
)

diagnostic_file <- file.path(metadata_dir, "download_diagnostic.csv")
write_csv(diagnostic, diagnostic_file, na = "")

source_metadata <- tibble(
  country = c("Brazil", "Argentina"),
  institution = c("IBGE", "MAGyP"),
  dataset = c(
    "Municipal Agricultural Production, PAM, table 5457",
    "Agricultural Estimates"
  ),
  crop_selection = c("Soja em grão, code 40124", "soja total"),
  spatial_level = c("municipality", "department or partido"),
  source_url = c(
    "https://sidra.ibge.gov.br/tabela/5457",
    argentina_url
  ),
  downloaded_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

metadata_file <- file.path(metadata_dir, "agricultural_data_sources.csv")
write_csv(source_metadata, metadata_file, na = "")


# 8. FINAL REPORT --------------------------------------------------------------

message("")
message("Download completed successfully.")
message("Project directory: ", normalizePath(project_dir, winslash = "/"))
message("Brazil file: ", normalizePath(brazil_output, winslash = "/"))
message("Argentina file: ", normalizePath(argentina_output, winslash = "/"))
message("Diagnostic: ", normalizePath(diagnostic_file, winslash = "/"))
message("")

print(diagnostic)

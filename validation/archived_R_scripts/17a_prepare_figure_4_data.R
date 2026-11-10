#!/usr/bin/env Rscript

# Figure 4 data preparation
# Multivariate climate regimes associated with soybean productivity anomalies.
#
# The PCA axes are estimated independently within each macroregion using only
# climate and crop water balance variables. ENSO phase and soybean productivity
# are supplementary information and therefore do not define the ordination.

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
  file.path(project_dir, "outputs", "figure_data", "figure_4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- 1985L
analysis_end_year <- 2024L
permutation_repetitions <- 999L
bootstrap_repetitions <- 1000L
random_seed <- 20260907L
weight_cap_probability <- 0.99

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

normalize_region <- function(x) {
  value <- iconv(str_to_upper(str_squish(as.character(x))),
                 from = "", to = "ASCII//TRANSLIT")
  gsub("[^A-Z]", "", value)
}

classify_macroregion <- function(country, parent_region) {
  region <- normalize_region(parent_region)
  southern <- c("PR", "PARANA", "SC", "SANTACATARINA", "RS", "RIOGRANDEDOSUL")
  frontier <- c(
    "TO", "TOCANTINS", "MA", "MARANHAO", "PI", "PIAUI", "BA", "BAHIA",
    "AL", "ALAGOAS", "CE", "CEARA", "PB", "PARAIBA", "PE", "PERNAMBUCO",
    "RN", "RIOGRANDEDONORTE", "SE", "SERGIPE"
  )
  case_when(
    country == "Argentina" ~ "Argentina",
    country == "Brazil" & region %in% southern ~ "Southern Brazil",
    country == "Brazil" & region %in% frontier ~
      "MATOPIBA and northeastern frontier",
    country == "Brazil" ~ "Central Brazil",
    TRUE ~ "Outside focal regions"
  )
}

standardize_phase <- function(x) {
  phase <- normalize_region(x)
  case_when(
    phase %in% c("ELNINO", "ELNIO") ~ "El Nino",
    phase %in% c("LANINA", "LANIA") ~ "La Nina",
    TRUE ~ "Neutral"
  )
}

safe_weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  weighted.mean(x[keep], w[keep])
}

safe_weighted_sd <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (sum(keep) < 2L) return(NA_real_)
  x <- x[keep]
  w <- w[keep]
  centre <- weighted.mean(x, w)
  sqrt(sum(w * (x - centre)^2) / sum(w))
}

safe_weighted_cor <- function(x, y, w) {
  keep <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  if (sum(keep) < 3L) return(NA_real_)
  x <- x[keep]
  y <- y[keep]
  w <- w[keep]
  mx <- weighted.mean(x, w)
  my <- weighted.mean(y, w)
  vx <- sum(w * (x - mx)^2) / sum(w)
  vy <- sum(w * (y - my)^2) / sum(w)
  if (vx <= 0 || vy <= 0) return(NA_real_)
  sum(w * (x - mx) * (y - my)) / sum(w) / sqrt(vx * vy)
}

# Express every ordination variable as a within-location anomaly. Climate
# metrics arrive as standardized local anomalies, whereas the water-balance
# metrics arrive in physical units. Re-standardizing both families against
# the complete local 1985--2024 record puts them on the same definition and
# prevents persistent differences among spatial units from defining the PCA.
local_standardized_anomaly <- function(x) {
  centre <- mean(x, na.rm = TRUE)
  spread <- sd(x, na.rm = TRUE)
  if (!is.finite(centre)) return(rep(NA_real_, length(x)))
  if (!is.finite(spread) || spread <= 1e-8) return(rep(0, length(x)))
  (x - centre) / spread
}

weighted_pca <- function(data, variables, weight_column = "analysis_weight") {
  weights <- data[[weight_column]]
  matrix_raw <- as.matrix(data[, variables, drop = FALSE])
  centres <- vapply(variables, function(v) safe_weighted_mean(data[[v]], weights),
                    numeric(1))
  scales <- vapply(variables, function(v) safe_weighted_sd(data[[v]], weights),
                   numeric(1))
  retained <- is.finite(scales) & scales > 1e-8
  if (sum(retained) < 3L) {
    stop("Fewer than three informative climate variables in a macroregion.",
         call. = FALSE)
  }
  variables <- variables[retained]
  matrix_raw <- matrix_raw[, retained, drop = FALSE]
  centres <- centres[retained]
  scales <- scales[retained]
  matrix_z <- sweep(sweep(matrix_raw, 2L, centres, "-"), 2L, scales, "/")
  covariance <- crossprod(matrix_z * sqrt(weights / sum(weights)))
  decomposition <- eigen(covariance, symmetric = TRUE)
  score <- matrix_z %*% decomposition$vectors
  colnames(score) <- paste0("PC", seq_len(ncol(score)))
  list(
    score = score,
    variables = variables,
    centre = centres,
    scale = scales,
    eigenvalues = decomposition$values,
    vectors = decomposition$vectors,
    variance_fraction = decomposition$values / sum(decomposition$values)
  )
}

weighted_phase_statistic <- function(pc1, pc2, phase, weights) {
  keep <- is.finite(pc1) & is.finite(pc2) & is.finite(weights) & weights > 0 &
    phase %in% c("El Nino", "Neutral", "La Nina")
  pc1 <- pc1[keep]
  pc2 <- pc2[keep]
  phase <- phase[keep]
  weights <- weights[keep]
  if (length(unique(phase)) < 2L) return(NA_real_)
  grand <- c(weighted.mean(pc1, weights), weighted.mean(pc2, weights))
  groups <- split(seq_along(phase), phase)
  between <- sum(vapply(groups, function(index) {
    centre <- c(weighted.mean(pc1[index], weights[index]),
                weighted.mean(pc2[index], weights[index]))
    sum(weights[index]) * sum((centre - grand)^2)
  }, numeric(1)))
  within <- sum(vapply(groups, function(index) {
    centre <- c(weighted.mean(pc1[index], weights[index]),
                weighted.mean(pc2[index], weights[index]))
    sum(weights[index] * ((pc1[index] - centre[1])^2 +
                            (pc2[index] - centre[2])^2))
  }, numeric(1)))
  if (!is.finite(within) || within <= 0) return(NA_real_)
  between / within
}

climate_metrics <- c(
  "mean_temperature_c", "precipitation_total_mm", "mean_vpd_kpa",
  "hot_dry_fraction", "hot_wet_fraction", "heavy_rain_fraction"
)
water_metrics <- c(
  "reproductive_requirement_satisfaction"
)
climate_variables <- paste0("climate_z_", climate_metrics)
water_variables <- paste0("water_raw_", water_metrics)
ordination_variables <- c(climate_variables, water_variables)

variable_dictionary <- tribble(
  ~variable, ~variable_label, ~mechanism_group,
  "climate_z_mean_temperature_c", "Temperature", "Thermal demand",
  "climate_z_precipitation_total_mm", "Precipitation", "Water availability",
  "climate_z_mean_vpd_kpa", "VPD", "Atmospheric water demand",
  "climate_z_hot_dry_fraction", "Hot-dry conditions", "Compound extremes",
  "climate_z_hot_wet_fraction", "Hot-wet conditions", "Compound extremes",
  "climate_z_heavy_rain_fraction", "Heavy rainfall", "Compound extremes",
  "water_raw_reproductive_requirement_satisfaction", "CWS", "Crop water stress"
)

paths <- list(
  manifest = find_one(project_dir, c("NASA_POWER_grid_manifest.csv")),
  climate = find_one(project_dir, c("soybean_climate_indices_grid_year.csv.gz")),
  soybean = find_one(project_dir, c("soybean_ENSO_spatial_panel_1974_2024.csv")),
  water = find_one(project_dir, c("soybean_water_balance_primary_AWC100.csv"))
)

manifest <- read_csv(paths$manifest, show_col_types = FALSE)
climate <- read_csv(paths$climate, show_col_types = FALSE)
soybean <- read_csv(paths$soybean, show_col_types = FALSE)
water <- read_csv(paths$water, show_col_types = FALSE)

check_required_columns(
  manifest,
  c("country", "spatial_unit_id", "parent_region_name", "grid_id", "calendar_group"),
  "Grid manifest"
)
check_required_columns(
  climate,
  c("grid_id", "country", "calendar_group", "climate_window", "harvest_year",
    "enso_phase_simple", "climate_metric", "climate_anomaly_standardized"),
  "Climate indices"
)
check_required_columns(
  soybean,
  c("country", "spatial_unit_id", "parent_region_name", "analysis_year",
    "yield_anomaly_percent", "harvested_area_ha"),
  "Soybean panel"
)
check_required_columns(
  water,
  c("country", "spatial_unit_id", "harvest_year", "awc_mm", water_metrics),
  "Water balance"
)

manifest_prepared <- manifest |>
  transmute(
    country = as.character(country),
    spatial_unit_id = standardize_unit_id(country, spatial_unit_id),
    manifest_parent_region = as.character(parent_region_name),
    grid_id = as.character(grid_id),
    calendar_group = as.character(calendar_group)
  ) |>
  distinct(country, spatial_unit_id, .keep_all = TRUE)

climate_long <- climate |>
  filter(
    climate_window == "reproductive",
    climate_metric %in% climate_metrics,
    harvest_year >= analysis_start_year,
    harvest_year <= analysis_end_year
  ) |>
  select(
    grid_id, country, calendar_group, harvest_year, enso_phase_simple,
    climate_metric, climate_anomaly_standardized
  ) |>
  distinct()

climate_duplicates <- climate_long |>
  count(grid_id, country, calendar_group, harvest_year, enso_phase_simple,
        climate_metric, name = "rows") |>
  filter(rows > 1L)
if (nrow(climate_duplicates) > 0L) {
  stop("Climate indices are not unique by grid, calendar, year, phase, and metric.",
       call. = FALSE)
}

climate_wide <- climate_long |>
  pivot_wider(
    id_cols = c(grid_id, country, calendar_group, harvest_year,
                enso_phase_simple),
    names_from = climate_metric,
    values_from = climate_anomaly_standardized,
    names_prefix = "climate_z_"
  )

if (any(vapply(climate_wide[climate_variables], is.list, logical(1)))) {
  stop("Climate indices are not unique by grid, calendar, year, and metric.",
       call. = FALSE)
}

water_prepared <- water |>
  filter(
    awc_mm == 100,
    harvest_year >= analysis_start_year,
    harvest_year <= analysis_end_year
  ) |>
  transmute(
    country = as.character(country),
    spatial_unit_id = standardize_unit_id(country, spatial_unit_id),
    analysis_year = as.integer(harvest_year),
    across(all_of(water_metrics), as.numeric, .names = "water_raw_{.col}")
  ) |>
  # The source variable is a satisfaction ratio: high values mean that crop
  # water requirements were met. Reverse its orientation so that CWS has the
  # same ecological interpretation everywhere: high values = greater stress.
  mutate(
    water_raw_reproductive_requirement_satisfaction =
      -water_raw_reproductive_requirement_satisfaction
  )

if (any(duplicated(water_prepared[c("country", "spatial_unit_id", "analysis_year")]))) {
  stop("Duplicated water balance rows were detected.", call. = FALSE)
}

soybean_prepared <- soybean |>
  transmute(
    country = as.character(country),
    spatial_unit_id = standardize_unit_id(country, spatial_unit_id),
    soybean_parent_region = as.character(parent_region_name),
    analysis_year = as.integer(analysis_year),
    productivity_anomaly_percent = as.numeric(yield_anomaly_percent),
    harvested_area_ha = as.numeric(harvested_area_ha)
  ) |>
  filter(
    analysis_year >= analysis_start_year,
    analysis_year <= analysis_end_year,
    is.finite(productivity_anomaly_percent),
    is.finite(harvested_area_ha), harvested_area_ha > 0
  )

panel <- soybean_prepared |>
  inner_join(manifest_prepared, by = c("country", "spatial_unit_id")) |>
  inner_join(
    climate_wide,
    by = c("country", "grid_id", "calendar_group",
           "analysis_year" = "harvest_year")
  ) |>
  inner_join(water_prepared,
             by = c("country", "spatial_unit_id", "analysis_year")) |>
  mutate(
    parent_region_name = coalesce(soybean_parent_region, manifest_parent_region),
    macroregion = classify_macroregion(country, parent_region_name),
    ENSO_phase = standardize_phase(enso_phase_simple)
  ) |>
  filter(macroregion != "Outside focal regions") |>
  select(
    country, spatial_unit_id, parent_region_name, macroregion, analysis_year,
    ENSO_phase, productivity_anomaly_percent, harvested_area_ha,
    all_of(ordination_variables)
  ) |>
  filter(if_all(all_of(ordination_variables), is.finite)) |>
  group_by(country, spatial_unit_id) |>
  mutate(across(
    all_of(ordination_variables),
    local_standardized_anomaly
  )) |>
  ungroup() |>
  filter(if_all(all_of(ordination_variables), is.finite)) |>
  group_by(macroregion) |>
  mutate(
    area_cap_ha = quantile(harvested_area_ha, weight_cap_probability,
                           na.rm = TRUE, names = FALSE, type = 8),
    analysis_weight = pmin(harvested_area_ha, area_cap_ha),
    analysis_weight = analysis_weight / mean(analysis_weight, na.rm = TRUE)
  ) |>
  ungroup()

if (nrow(panel) == 0L) stop("The climate and productivity panels did not match.")
if (anyDuplicated(panel[c("country", "spatial_unit_id", "analysis_year")])) {
  stop("The analytical panel contains duplicated spatial unit years.", call. = FALSE)
}

macroregion_order <- c(
  "Argentina", "Southern Brazil", "Central Brazil",
  "MATOPIBA and northeastern frontier"
)

pca_objects <- split(panel, as.character(panel$macroregion)) |>
  map(~ weighted_pca(.x, ordination_variables))

score_tables <- imap(pca_objects, function(pca, selected_macroregion) {
  data <- panel |> filter(macroregion == selected_macroregion)
  bind_cols(
    data |>
      select(country, spatial_unit_id, parent_region_name, macroregion,
             analysis_year, ENSO_phase, productivity_anomaly_percent,
             harvested_area_ha, analysis_weight),
    as_tibble(pca$score[, 1:2, drop = FALSE])
  )
})
scores <- bind_rows(score_tables) |>
  mutate(macroregion = factor(macroregion, levels = macroregion_order)) |>
  arrange(macroregion, analysis_year, spatial_unit_id)

variance_summary <- imap_dfr(pca_objects, function(pca, selected_macroregion) {
  tibble(
    macroregion = selected_macroregion,
    component = paste0("PC", seq_along(pca$eigenvalues)),
    eigenvalue = pca$eigenvalues,
    variance_fraction = pca$variance_fraction,
    cumulative_variance_fraction = cumsum(pca$variance_fraction)
  )
})

loading_summary <- imap_dfr(pca_objects, function(pca, selected_macroregion) {
  tibble(
    macroregion = selected_macroregion,
    variable = pca$variables,
    PC1_loading = pca$vectors[, 1] * sqrt(pca$eigenvalues[1]),
    PC2_loading = pca$vectors[, 2] * sqrt(pca$eigenvalues[2]),
    PC1_eigenvector = pca$vectors[, 1],
    PC2_eigenvector = pca$vectors[, 2]
  )
}) |>
  left_join(variable_dictionary, by = "variable") |>
  group_by(macroregion) |>
  mutate(
    loading_length = sqrt(PC1_loading^2 + PC2_loading^2),
    loading_rank = rank(-loading_length, ties.method = "first"),
    display_loading = loading_rank <= 6L
  ) |>
  ungroup()

productivity_vectors <- scores |>
  group_by(macroregion) |>
  summarise(
    PC1_correlation = safe_weighted_cor(
      productivity_anomaly_percent, PC1, analysis_weight
    ),
    PC2_correlation = safe_weighted_cor(
      productivity_anomaly_percent, PC2, analysis_weight
    ),
    observations = n(),
    productive_area_years_ha = sum(harvested_area_ha, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    variable_label = "Productivity anomaly",
    vector_length = sqrt(PC1_correlation^2 + PC2_correlation^2)
  )

phase_centroids <- scores |>
  group_by(macroregion, ENSO_phase) |>
  summarise(
    PC1 = safe_weighted_mean(PC1, analysis_weight),
    PC2 = safe_weighted_mean(PC2, analysis_weight),
    spatial_unit_years = n(),
    years = n_distinct(analysis_year),
    productive_area_years_ha = sum(harvested_area_ha, na.rm = TRUE),
    .groups = "drop"
  )

set.seed(random_seed)
bootstrap_centroids <- scores |>
  group_by(macroregion) |>
  group_modify(~ {
    region_data <- .x
    years <- sort(unique(region_data$analysis_year))
    map_dfr(seq_len(bootstrap_repetitions), function(iteration) {
      sampled_years <- sample(years, length(years), replace = TRUE)
      sampled <- map2_dfr(sampled_years, seq_along(sampled_years),
                          function(year, draw) {
        region_data |>
          filter(analysis_year == year) |>
          mutate(bootstrap_year = draw)
      })
      sampled |>
        group_by(ENSO_phase) |>
        summarise(
          PC1 = safe_weighted_mean(PC1, analysis_weight),
          PC2 = safe_weighted_mean(PC2, analysis_weight),
          .groups = "drop"
        ) |>
        mutate(iteration = iteration)
    })
  }) |>
  ungroup()

phase_centroid_intervals <- bootstrap_centroids |>
  group_by(macroregion, ENSO_phase) |>
  summarise(
    PC1_lower = quantile(PC1, 0.025, na.rm = TRUE, names = FALSE, type = 8),
    PC1_upper = quantile(PC1, 0.975, na.rm = TRUE, names = FALSE, type = 8),
    PC2_lower = quantile(PC2, 0.025, na.rm = TRUE, names = FALSE, type = 8),
    PC2_upper = quantile(PC2, 0.975, na.rm = TRUE, names = FALSE, type = 8),
    .groups = "drop"
  )

phase_centroids <- phase_centroids |>
  left_join(phase_centroid_intervals, by = c("macroregion", "ENSO_phase"))

set.seed(random_seed + 1L)
phase_separation <- scores |>
  group_by(macroregion) |>
  group_modify(~ {
    region_data <- .x
    year_phase <- region_data |>
      distinct(analysis_year, ENSO_phase) |>
      arrange(analysis_year)
    if (anyDuplicated(year_phase$analysis_year)) {
      stop("ENSO phase is not unique within year.", call. = FALSE)
    }
    observed <- weighted_phase_statistic(
      region_data$PC1, region_data$PC2, region_data$ENSO_phase,
      region_data$analysis_weight
    )
    null_statistics <- replicate(permutation_repetitions, {
      permuted <- year_phase |>
        mutate(permuted_phase = sample(ENSO_phase, replace = FALSE)) |>
        select(analysis_year, permuted_phase)
      phase_vector <- region_data |>
        select(analysis_year) |>
        left_join(permuted, by = "analysis_year") |>
        pull(permuted_phase)
      weighted_phase_statistic(
        region_data$PC1, region_data$PC2, phase_vector,
        region_data$analysis_weight
      )
    })
    tibble(
      observed_between_within_ratio = observed,
      permutation_p_value =
        (1 + sum(null_statistics >= observed, na.rm = TRUE)) /
        (1 + sum(is.finite(null_statistics))),
      null_median = median(null_statistics, na.rm = TRUE),
      null_95th_percentile = quantile(
        null_statistics, 0.95, na.rm = TRUE, names = FALSE, type = 8
      ),
      permutations = permutation_repetitions,
      years = n_distinct(region_data$analysis_year),
      spatial_unit_years = nrow(region_data)
    )
  }) |>
  ungroup()

qc_region_coverage <- scores |>
  group_by(macroregion, ENSO_phase) |>
  summarise(
    spatial_units = n_distinct(spatial_unit_id),
    years = n_distinct(analysis_year),
    spatial_unit_years = n(),
    productive_area_years_ha = sum(harvested_area_ha, na.rm = TRUE),
    .groups = "drop"
  )

qc_summary <- tribble(
  ~test, ~value,
  "Analytical spatial unit years", nrow(scores),
  "Unique spatial units", n_distinct(interaction(scores$country, scores$spatial_unit_id)),
  "Analysis years", n_distinct(scores$analysis_year),
  "Macroregions", n_distinct(scores$macroregion),
  "Minimum variables retained by PCA", min(map_int(pca_objects, ~ length(.x$variables))),
  "Minimum PC1 plus PC2 variance fraction",
    min(variance_summary |>
          filter(component == "PC2") |>
          pull(cumulative_variance_fraction)),
  "Missing PCA scores", sum(!is.finite(scores$PC1) | !is.finite(scores$PC2)),
  "Minimum years per phase and macroregion", min(qc_region_coverage$years)
)

methods <- tribble(
  ~item, ~value,
  "figure_question", "How do ENSO phases occupy regional climate regimes associated with soybean productivity anomalies?",
  "ordination", "Area weighted PCA estimated separately within each macroregion",
  "ordination_inputs", "Within-location standardized anomalies of seven harmonized reproductive-season indices: temperature, VPD, precipitation, heavy rainfall, CWS, hot-dry conditions, and hot-wet conditions",
  "anomaly_reference", "Each climate and water-balance variable was centred and scaled against its complete 1985-2024 record within each productive spatial unit; El Nino, neutral, and La Nina seasons use the same reference",
  "noncircularity", "ENSO phase and productivity were excluded when estimating PCA axes",
  "standardization", "Weighted centering and scaling within each macroregion",
  "area_weight", "Harvested area capped at the macroregion specific 99th percentile",
  "points", "Each observation is a productive spatial unit-season; colours identify ENSO phase after anomaly calculation",
  "phase_summary", "Area weighted centroid with 95% year block bootstrap intervals",
  "phase_test", "Between to within PCA dispersion ratio with ENSO labels permuted among complete years",
  "productivity_overlay", "Area weighted correlations of productivity anomaly with PC1 and PC2",
  "interpretation", "The ordination is descriptive and predictive, not a causal attribution"
)

write_csv(scores,
          file.path(output_dir, "Figure_4_local_climate_PCA_scores.csv.gz"), na = "")
write_csv(loading_summary,
          file.path(output_dir, "Figure_4_macroregion_PCA_loadings.csv"), na = "")
write_csv(variance_summary,
          file.path(output_dir, "Figure_4_macroregion_PCA_variance.csv"), na = "")
write_csv(productivity_vectors,
          file.path(output_dir, "Figure_4_productivity_supplementary_vectors.csv"), na = "")
write_csv(phase_centroids,
          file.path(output_dir, "Figure_4_ENSO_phase_centroids.csv"), na = "")
write_csv(phase_separation,
          file.path(output_dir, "Figure_4_ENSO_phase_separation_tests.csv"), na = "")
write_csv(qc_region_coverage,
          file.path(output_dir, "QC_Figure_4_region_phase_coverage.csv"), na = "")
write_csv(qc_summary,
          file.path(output_dir, "QC_Figure_4_data_preparation.csv"), na = "")
write_csv(methods,
          file.path(output_dir, "README_Figure_4_methods.csv"), na = "")

message("Figure 4 data written to: ", normalizePath(output_dir))

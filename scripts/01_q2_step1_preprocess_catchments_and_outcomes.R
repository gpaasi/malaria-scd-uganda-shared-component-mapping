# ============================================================
# OBJECTIVE 2, QUESTION 2
# Step 1: Import, clean, join, extract population, and preprocess
# Spatial co-clustering of malaria and SCD across 2SFCA catchments
# ============================================================

# ----------------------------
# 0. GLOBAL OPTIONS
# ----------------------------
options(stringsAsFactors = FALSE)
options(scipen = 999)

# ----------------------------
# 1. LOAD / INSTALL PACKAGES
# ----------------------------
required_pkgs <- c(
  "sf",
  "terra",
  "exactextractr",
  "dplyr",
  "readr",
  "janitor",
  "stringr",
  "tidyr",
  "purrr",
  "tibble",
  "ggplot2",
  "tmap"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_pkgs, installed)

if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

# ----------------------------
# 2. SET ROOT DIRECTORY
# ----------------------------
root_dir <- "C:/SCD-Malaria syndemic"

if (!dir.exists(root_dir)) {
  stop("Root directory does not exist. Check 'root_dir'.")
}

setwd(root_dir)

# ----------------------------
# 3. DEFINE PROJECT FOLDERS
# ----------------------------
data_dir    <- file.path(root_dir, "data")
catch_dir   <- file.path(data_dir, "catchments")
table_dir_in <- file.path(data_dir, "tables")
raster_dir  <- file.path(data_dir, "rasters")

output_dir  <- file.path(root_dir, "outputs")
fig_dir     <- file.path(output_dir, "figures")
table_dir   <- file.path(output_dir, "tables")
rds_dir     <- file.path(output_dir, "rds")
gpkg_dir    <- file.path(output_dir, "gpkg")
log_dir     <- file.path(output_dir, "logs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gpkg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 4. DEFINE INPUT FILE PATHS
# ----------------------------
catchment_file <- file.path(catch_dir, "catchment_area_revised.shp")

malaria_cases_file  <- file.path(table_dir_in, "Malaria_cases_2020_2024.csv")
malaria_deaths_file <- file.path(table_dir_in, "Malaria_deaths_2020_2024.csv")
scd_cases_file      <- file.path(table_dir_in, "SCD_cases_2020_2024.csv")
scd_deaths_file     <- file.path(table_dir_in, "SCD_deaths_2020_2024.csv")

population_raster_file <- file.path(raster_dir, "u5_pop_utm.tif")

input_files <- c(
  catchment_file,
  malaria_cases_file,
  malaria_deaths_file,
  scd_cases_file,
  scd_deaths_file,
  population_raster_file
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop(
    "The following required files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ----------------------------
# 5. START LOG
# ----------------------------
log_file <- file.path(log_dir, "q2_step1_preprocessing_log.txt")
cat("Q2 Step 1 started:", as.character(Sys.time()), "\n", file = log_file)
cat("Root directory:", root_dir, "\n", file = log_file, append = TRUE)

# ----------------------------
# 6. HELPER FUNCTIONS
# ----------------------------

# Clean and standardise ID values
clean_id <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_to_upper()
}

# Find year columns robustly after clean_names()
# This catches names such as:
# x2020_malaria, x2021_scd, deaths_2020, deaths_2021, etc.
find_year_cols <- function(df) {
  names(df)[stringr::str_detect(names(df), "2020|2021|2022|2023|2024")]
}

# Convert year columns to numeric and compute 2020-2024 total
prep_outcome_table <- function(file_path, value_name) {

  dat <- readr::read_csv(file_path, show_col_types = FALSE) %>%
    janitor::clean_names()

  if (!"id" %in% names(dat)) {
    stop("Column 'ID' not found in: ", basename(file_path))
  }

  dat <- dat %>%
    mutate(id = clean_id(id))

  year_cols <- find_year_cols(dat)

  if (length(year_cols) == 0) {
    stop("No 2020-2024 year columns detected in: ", basename(file_path))
  }

  dat <- dat %>%
    mutate(across(all_of(year_cols), ~ suppressWarnings(as.numeric(.)))) %>%
    mutate("{value_name}" := rowSums(across(all_of(year_cols)), na.rm = TRUE)) %>%
    select(id, all_of(value_name), any_of("name"), any_of("type"))

  # keep only one row per id
  dat <- dat %>%
    group_by(id) %>%
    summarise(
      "{value_name}" := sum(.data[[value_name]], na.rm = TRUE),
      .groups = "drop"
    )

  return(dat)
}

# Safe z-score
z_score <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - m) / s
}

# ----------------------------
# 7. IMPORT CATCHMENT SHAPEFILE
# ----------------------------
catchments <- sf::st_read(catchment_file, quiet = FALSE) %>%
  janitor::clean_names()

if (!"id" %in% names(catchments)) {
  stop("Catchment shapefile does not contain an 'ID' field.")
}

catchments <- catchments %>%
  mutate(
    catchment_id = clean_id(id)
  )

cat("Catchments imported: ", nrow(catchments), "\n", file = log_file, append = TRUE)

# Check duplicate catchment IDs
dup_catchments <- catchments %>%
  st_drop_geometry() %>%
  count(catchment_id) %>%
  filter(n > 1)

if (nrow(dup_catchments) > 0) {
  write_csv(dup_catchments, file.path(table_dir, "duplicate_catchment_ids.csv"))
  stop("Duplicate catchment IDs found in shapefile. See outputs/tables/duplicate_catchment_ids.csv")
}

# ----------------------------
# 8. IMPORT OUTCOME TABLES
# ----------------------------
malaria_cases  <- prep_outcome_table(malaria_cases_file,  "malaria_cases")
malaria_deaths <- prep_outcome_table(malaria_deaths_file, "malaria_deaths")
scd_cases      <- prep_outcome_table(scd_cases_file,      "scd_cases")
scd_deaths     <- prep_outcome_table(scd_deaths_file,     "scd_deaths")

cat("Malaria cases rows: ", nrow(malaria_cases), "\n", file = log_file, append = TRUE)
cat("Malaria deaths rows: ", nrow(malaria_deaths), "\n", file = log_file, append = TRUE)
cat("SCD cases rows: ", nrow(scd_cases), "\n", file = log_file, append = TRUE)
cat("SCD deaths rows: ", nrow(scd_deaths), "\n", file = log_file, append = TRUE)

# ----------------------------
# 9. CHECK UNMATCHED IDs BEFORE JOIN
# ----------------------------
catchment_ids <- catchments %>% st_drop_geometry() %>% pull(catchment_id)

unmatched_mal_cases <- malaria_cases %>% filter(!id %in% catchment_ids)
unmatched_mal_deaths <- malaria_deaths %>% filter(!id %in% catchment_ids)
unmatched_scd_cases <- scd_cases %>% filter(!id %in% catchment_ids)
unmatched_scd_deaths <- scd_deaths %>% filter(!id %in% catchment_ids)

if (nrow(unmatched_mal_cases) > 0) {
  write_csv(unmatched_mal_cases, file.path(table_dir, "unmatched_ids_malaria_cases.csv"))
}
if (nrow(unmatched_mal_deaths) > 0) {
  write_csv(unmatched_mal_deaths, file.path(table_dir, "unmatched_ids_malaria_deaths.csv"))
}
if (nrow(unmatched_scd_cases) > 0) {
  write_csv(unmatched_scd_cases, file.path(table_dir, "unmatched_ids_scd_cases.csv"))
}
if (nrow(unmatched_scd_deaths) > 0) {
  write_csv(unmatched_scd_deaths, file.path(table_dir, "unmatched_ids_scd_deaths.csv"))
}

cat("Unmatched malaria case IDs: ", nrow(unmatched_mal_cases), "\n", file = log_file, append = TRUE)
cat("Unmatched malaria death IDs: ", nrow(unmatched_mal_deaths), "\n", file = log_file, append = TRUE)
cat("Unmatched SCD case IDs: ", nrow(unmatched_scd_cases), "\n", file = log_file, append = TRUE)
cat("Unmatched SCD death IDs: ", nrow(unmatched_scd_deaths), "\n", file = log_file, append = TRUE)

# ----------------------------
# 10. JOIN OUTCOMES TO CATCHMENTS
# ----------------------------
q2_data <- catchments %>%
  left_join(malaria_cases,  by = c("catchment_id" = "id")) %>%
  left_join(malaria_deaths, by = c("catchment_id" = "id")) %>%
  left_join(scd_cases,      by = c("catchment_id" = "id")) %>%
  left_join(scd_deaths,     by = c("catchment_id" = "id")) %>%
  mutate(
    malaria_cases  = coalesce(malaria_cases, 0),
    malaria_deaths = coalesce(malaria_deaths, 0),
    scd_cases      = coalesce(scd_cases, 0),
    scd_deaths     = coalesce(scd_deaths, 0)
  )

# ----------------------------
# 11. GEOMETRY VALIDATION
# ----------------------------
invalid_before <- sum(!sf::st_is_valid(q2_data))
cat("Invalid geometries before repair: ", invalid_before, "\n", file = log_file, append = TRUE)

if (invalid_before > 0) {
  q2_data <- sf::st_make_valid(q2_data)
}

invalid_after <- sum(!sf::st_is_valid(q2_data))
cat("Invalid geometries after repair: ", invalid_after, "\n", file = log_file, append = TRUE)

# ----------------------------
# 12. EXTRACT POPULATION FROM RASTER
# ----------------------------
pop_rast <- terra::rast(population_raster_file)

# Reproject catchments to raster CRS if needed
if (sf::st_crs(q2_data)$wkt != terra::crs(pop_rast)) {
  q2_data_pop <- sf::st_transform(q2_data, crs = terra::crs(pop_rast))
} else {
  q2_data_pop <- q2_data
}

# exactextractr works with Raster/SpatRaster support in many setups,
# but to keep this robust we use terra + exactextractr directly.
# Population is summed within each polygon.
q2_data_pop$catchment_population <- exactextractr::exact_extract(
  pop_rast,
  q2_data_pop,
  "sum"
)

# Put back into original object
q2_data$catchment_population <- q2_data_pop$catchment_population

# Replace any missing population with 0 and round
q2_data <- q2_data %>%
  mutate(
    catchment_population = coalesce(catchment_population, 0),
    catchment_population = round(catchment_population)
  )

cat("Population extracted from raster.\n", file = log_file, append = TRUE)

# ----------------------------
# 13. CREATE CENTROIDS / COORDINATES
# ----------------------------
centroids <- sf::st_centroid(sf::st_geometry(q2_data))
coords <- sf::st_coordinates(centroids)

q2_data <- q2_data %>%
  mutate(
    centroid_x = coords[, 1],
    centroid_y = coords[, 2]
  )

# ----------------------------
# 14. CREATE ESDA-READY STANDARDISED VARIABLES
# ----------------------------
q2_data <- q2_data %>%
  mutate(
    z_malaria_cases  = z_score(malaria_cases),
    z_scd_cases      = z_score(scd_cases),
    z_malaria_deaths = z_score(malaria_deaths),
    z_scd_deaths     = z_score(scd_deaths)
  )

# ----------------------------
# 15. BASIC COMPLETENESS CHECKS
# ----------------------------
completeness_summary <- tibble(
  n_catchments            = nrow(q2_data),
  n_with_population       = sum(q2_data$catchment_population > 0, na.rm = TRUE),
  n_with_malaria_cases    = sum(q2_data$malaria_cases > 0, na.rm = TRUE),
  n_with_malaria_deaths   = sum(q2_data$malaria_deaths > 0, na.rm = TRUE),
  n_with_scd_cases        = sum(q2_data$scd_cases > 0, na.rm = TRUE),
  n_with_scd_deaths       = sum(q2_data$scd_deaths > 0, na.rm = TRUE)
)

write_csv(completeness_summary, file.path(table_dir, "q2_completeness_summary.csv"))

# ----------------------------
# 16. DESCRIPTIVE SUMMARY
# ----------------------------
desc_table <- q2_data %>%
  st_drop_geometry() %>%
  summarise(
    catchment_population_mean = mean(catchment_population, na.rm = TRUE),
    catchment_population_median = median(catchment_population, na.rm = TRUE),

    malaria_cases_total   = sum(malaria_cases, na.rm = TRUE),
    malaria_cases_mean    = mean(malaria_cases, na.rm = TRUE),
    malaria_cases_median  = median(malaria_cases, na.rm = TRUE),

    malaria_deaths_total  = sum(malaria_deaths, na.rm = TRUE),
    malaria_deaths_mean   = mean(malaria_deaths, na.rm = TRUE),
    malaria_deaths_median = median(malaria_deaths, na.rm = TRUE),

    scd_cases_total       = sum(scd_cases, na.rm = TRUE),
    scd_cases_mean        = mean(scd_cases, na.rm = TRUE),
    scd_cases_median      = median(scd_cases, na.rm = TRUE),

    scd_deaths_total      = sum(scd_deaths, na.rm = TRUE),
    scd_deaths_mean       = mean(scd_deaths, na.rm = TRUE),
    scd_deaths_median     = median(scd_deaths, na.rm = TRUE)
  )

write_csv(desc_table, file.path(table_dir, "q2_descriptive_summary.csv"))

dist_table <- q2_data %>%
  st_drop_geometry() %>%
  summarise(
    pop_min = min(catchment_population, na.rm = TRUE),
    pop_q1  = quantile(catchment_population, 0.25, na.rm = TRUE),
    pop_med = median(catchment_population, na.rm = TRUE),
    pop_q3  = quantile(catchment_population, 0.75, na.rm = TRUE),
    pop_max = max(catchment_population, na.rm = TRUE),

    malaria_cases_min = min(malaria_cases, na.rm = TRUE),
    malaria_cases_q1  = quantile(malaria_cases, 0.25, na.rm = TRUE),
    malaria_cases_med = median(malaria_cases, na.rm = TRUE),
    malaria_cases_q3  = quantile(malaria_cases, 0.75, na.rm = TRUE),
    malaria_cases_max = max(malaria_cases, na.rm = TRUE),

    scd_cases_min = min(scd_cases, na.rm = TRUE),
    scd_cases_q1  = quantile(scd_cases, 0.25, na.rm = TRUE),
    scd_cases_med = median(scd_cases, na.rm = TRUE),
    scd_cases_q3  = quantile(scd_cases, 0.75, na.rm = TRUE),
    scd_cases_max = max(scd_cases, na.rm = TRUE),

    malaria_deaths_min = min(malaria_deaths, na.rm = TRUE),
    malaria_deaths_q1  = quantile(malaria_deaths, 0.25, na.rm = TRUE),
    malaria_deaths_med = median(malaria_deaths, na.rm = TRUE),
    malaria_deaths_q3  = quantile(malaria_deaths, 0.75, na.rm = TRUE),
    malaria_deaths_max = max(malaria_deaths, na.rm = TRUE),

    scd_deaths_min = min(scd_deaths, na.rm = TRUE),
    scd_deaths_q1  = quantile(scd_deaths, 0.25, na.rm = TRUE),
    scd_deaths_med = median(scd_deaths, na.rm = TRUE),
    scd_deaths_q3  = quantile(scd_deaths, 0.75, na.rm = TRUE),
    scd_deaths_max = max(scd_deaths, na.rm = TRUE)
  )

write_csv(dist_table, file.path(table_dir, "q2_distribution_summary.csv"))

# ----------------------------
# 17. OPTIONAL CRUDE RATES PER 1,000
# ----------------------------
# These are descriptive only.
q2_data <- q2_data %>%
  mutate(
    malaria_cases_rate_per_1000 = if_else(catchment_population > 0,
                                          (malaria_cases / catchment_population) * 1000, NA_real_),
    scd_cases_rate_per_1000 = if_else(catchment_population > 0,
                                      (scd_cases / catchment_population) * 1000, NA_real_),
    malaria_deaths_rate_per_1000 = if_else(catchment_population > 0,
                                           (malaria_deaths / catchment_population) * 1000, NA_real_),
    scd_deaths_rate_per_1000 = if_else(catchment_population > 0,
                                       (scd_deaths / catchment_population) * 1000, NA_real_)
  )

# ----------------------------
# 18. QUICK QC MAPS
# ----------------------------
tmap_mode("plot")

map_mal_cases <- tm_shape(q2_data) +
  tm_polygons("malaria_cases", title = "Malaria cases") +
  tm_layout(main.title = "QC: Catchment-level malaria cases")

map_scd_cases <- tm_shape(q2_data) +
  tm_polygons("scd_cases", title = "SCD cases") +
  tm_layout(main.title = "QC: Catchment-level SCD cases")

map_mal_deaths <- tm_shape(q2_data) +
  tm_polygons("malaria_deaths", title = "Malaria deaths") +
  tm_layout(main.title = "QC: Catchment-level malaria deaths")

map_scd_deaths <- tm_shape(q2_data) +
  tm_polygons("scd_deaths", title = "SCD deaths") +
  tm_layout(main.title = "QC: Catchment-level SCD deaths")

map_population <- tm_shape(q2_data) +
  tm_polygons("catchment_population", title = "Population") +
  tm_layout(main.title = "QC: Catchment population from raster")

tmap_save(map_mal_cases,
          filename = file.path(fig_dir, "q2_qc_malaria_cases_map.png"),
          width = 8, height = 6, dpi = 300)

tmap_save(map_scd_cases,
          filename = file.path(fig_dir, "q2_qc_scd_cases_map.png"),
          width = 8, height = 6, dpi = 300)

tmap_save(map_mal_deaths,
          filename = file.path(fig_dir, "q2_qc_malaria_deaths_map.png"),
          width = 8, height = 6, dpi = 300)

tmap_save(map_scd_deaths,
          filename = file.path(fig_dir, "q2_qc_scd_deaths_map.png"),
          width = 8, height = 6, dpi = 300)

tmap_save(map_population,
          filename = file.path(fig_dir, "q2_qc_population_map.png"),
          width = 8, height = 6, dpi = 300)

# ----------------------------
# 19. EXPORT CLEAN NON-SPATIAL TABLE
# ----------------------------
q2_tabular <- q2_data %>%
  st_drop_geometry()

write_csv(q2_tabular, file.path(table_dir, "q2_catchment_burden_clean.csv"))

# ----------------------------
# 20. SAVE CLEAN SPATIAL OBJECTS
# ----------------------------
sf::st_write(
  q2_data,
  dsn = file.path(gpkg_dir, "q2_catchment_burden_clean.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

saveRDS(q2_data, file.path(rds_dir, "q2_data_step1_cleaned.rds"))

# ----------------------------
# 21. FINAL LOG MESSAGE
# ----------------------------
cat("Q2 Step 1 completed successfully: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE)

message("Step 1 completed successfully.")

# ============================================================
# OBJECTIVE 2, QUESTION 2
# Step 2: Spatial neighbours, weights, and global Moran's I
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
  "spdep",
  "dplyr",
  "readr",
  "tibble",
  "purrr",
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
output_dir <- file.path(root_dir, "outputs")
fig_dir    <- file.path(output_dir, "figures")
table_dir  <- file.path(output_dir, "tables")
rds_dir    <- file.path(output_dir, "rds")
gpkg_dir   <- file.path(output_dir, "gpkg")
log_dir    <- file.path(output_dir, "logs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gpkg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 4. DEFINE INPUTS
# ----------------------------
step1_file <- file.path(rds_dir, "q2_data_step1_cleaned.rds")

if (!file.exists(step1_file)) {
  stop("Step 1 cleaned file not found: ", step1_file)
}

# ----------------------------
# 5. START LOG
# ----------------------------
log_file <- file.path(log_dir, "q2_step2_spatial_weights_global_moran_log.txt")
cat("Q2 Step 2 started:", as.character(Sys.time()), "\n", file = log_file)
cat("Root directory:", root_dir, "\n", file = log_file, append = TRUE)
cat("Input file:", step1_file, "\n", file = log_file, append = TRUE)

# ----------------------------
# 6. HELPER FUNCTIONS
# ----------------------------
z_score <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - m) / s
}

safe_log1p <- function(x) {
  log1p(pmax(x, 0))
}

run_moran_suite <- function(x, listw, var_name, nsim = 999, zero.policy = TRUE) {
  x <- as.numeric(x)

  analytic <- spdep::moran.test(
    x = x,
    listw = listw,
    zero.policy = zero.policy,
    alternative = "two.sided"
  )

  mc <- spdep::moran.mc(
    x = x,
    listw = listw,
    nsim = nsim,
    zero.policy = zero.policy,
    alternative = "two.sided"
  )

  tibble(
    variable = var_name,
    n = length(x),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    moran_i = unname(analytic$estimate[["Moran I statistic"]]),
    expected_i = unname(analytic$estimate[["Expectation"]]),
    variance_i = unname(analytic$estimate[["Variance"]]),
    moran_test_statistic = unname(analytic$statistic),
    moran_test_p_value = analytic$p.value,
    mc_observed_i = unname(mc$statistic),
    mc_p_value = mc$p.value,
    nsim = nsim
  )
}

# ----------------------------
# 7. LOAD STEP 1 DATA
# ----------------------------
q2_data <- readRDS(step1_file)

if (!inherits(q2_data, "sf")) {
  stop("Step 1 object is not an sf object.")
}

required_cols <- c(
  "catchment_id",
  "malaria_cases",
  "malaria_deaths",
  "scd_cases",
  "scd_deaths",
  "catchment_population"
)

missing_cols <- setdiff(required_cols, names(q2_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns in Step 1 data: ", paste(missing_cols, collapse = ", "))
}

cat("Rows loaded: ", nrow(q2_data), "\n", file = log_file, append = TRUE)
cat("Columns loaded: ", ncol(q2_data), "\n", file = log_file, append = TRUE)

# ----------------------------
# 8. BASIC GEOMETRY CHECKS
# ----------------------------
invalid_before <- sum(!sf::st_is_valid(q2_data))
cat("Invalid geometries before Step 2 repair: ", invalid_before, "\n", file = log_file, append = TRUE)

if (invalid_before > 0) {
  q2_data <- sf::st_make_valid(q2_data)
}

invalid_after <- sum(!sf::st_is_valid(q2_data))
cat("Invalid geometries after Step 2 repair: ", invalid_after, "\n", file = log_file, append = TRUE)

# Remove empty geometries if present
empty_geom <- sf::st_is_empty(q2_data)
if (any(empty_geom)) {
  empty_ids <- q2_data$catchment_id[empty_geom]
  write_csv(
    tibble(catchment_id = empty_ids),
    file.path(table_dir, "q2_empty_geometries_removed.csv")
  )
  q2_data <- q2_data[!empty_geom, ]
}

cat("Rows after empty geometry check: ", nrow(q2_data), "\n", file = log_file, append = TRUE)

# ----------------------------
# 9. ENFORCE STABLE ROW ORDER
# ----------------------------
q2_data <- q2_data %>%
  arrange(catchment_id) %>%
  mutate(area_index = row_number())

# ----------------------------
# 10. CREATE ANALYSIS VARIABLES
# ----------------------------
q2_data <- q2_data %>%
  mutate(
    z_malaria_cases   = z_score(malaria_cases),
    z_malaria_deaths  = z_score(malaria_deaths),
    z_scd_cases       = z_score(scd_cases),
    z_scd_deaths      = z_score(scd_deaths),

    log_malaria_cases  = safe_log1p(malaria_cases),
    log_malaria_deaths = safe_log1p(malaria_deaths),
    log_scd_cases      = safe_log1p(scd_cases),
    log_scd_deaths     = safe_log1p(scd_deaths)
  )

# ----------------------------
# 11. BUILD QUEEN CONTIGUITY NEIGHBOURS
# ----------------------------
# Queen contiguity is appropriate for polygon-based areal ESDA.
nb_queen <- spdep::poly2nb(
  pl = q2_data,
  queen = TRUE,
  row.names = q2_data$catchment_id
)

# Summary stats
n_areas <- length(nb_queen)
n_links <- sum(card(nb_queen))
n_islands <- sum(card(nb_queen) == 0)

comp <- spdep::n.comp.nb(nb_queen)
n_components <- comp$nc

nb_summary <- tibble(
  n_areas = n_areas,
  total_links = n_links,
  mean_neighbours = mean(card(nb_queen)),
  median_neighbours = median(card(nb_queen)),
  min_neighbours = min(card(nb_queen)),
  max_neighbours = max(card(nb_queen)),
  n_islands = n_islands,
  n_connected_components = n_components
)

write_csv(nb_summary, file.path(table_dir, "q2_queen_neighbour_summary.csv"))

cat("Neighbour summary written.\n", file = log_file, append = TRUE)
cat("Number of islands: ", n_islands, "\n", file = log_file, append = TRUE)
cat("Connected components: ", n_components, "\n", file = log_file, append = TRUE)

# Save neighbour counts by catchment
nb_detail <- tibble(
  catchment_id = q2_data$catchment_id,
  area_index = q2_data$area_index,
  neighbour_count = card(nb_queen),
  component_id = comp$comp.id,
  is_island = card(nb_queen) == 0
)

write_csv(nb_detail, file.path(table_dir, "q2_queen_neighbour_detail.csv"))

# Save islands if any
if (n_islands > 0) {
  islands_tbl <- nb_detail %>%
    filter(is_island)

  write_csv(islands_tbl, file.path(table_dir, "q2_queen_islands.csv"))
}

# ----------------------------
# 12. CREATE SPATIAL WEIGHTS
# ----------------------------
# Row-standardised weights (style = "W")
lw_queen_W <- spdep::nb2listw(
  neighbours = nb_queen,
  style = "W",
  zero.policy = TRUE
)

# Optional binary weights for sensitivity / later use
lw_queen_B <- spdep::nb2listw(
  neighbours = nb_queen,
  style = "B",
  zero.policy = TRUE
)

cat("Spatial weights created.\n", file = log_file, append = TRUE)

# ----------------------------
# 13. SAVE NEIGHBOUR / WEIGHTS OBJECTS
# ----------------------------
saveRDS(nb_queen,  file.path(rds_dir, "q2_nb_queen.rds"))
saveRDS(lw_queen_W, file.path(rds_dir, "q2_listw_queen_W.rds"))
saveRDS(lw_queen_B, file.path(rds_dir, "q2_listw_queen_B.rds"))

# ----------------------------
# 14. CREATE SPATIAL LAGS
# ----------------------------
q2_data <- q2_data %>%
  mutate(
    lag_malaria_cases   = spdep::lag.listw(lw_queen_W, malaria_cases, zero.policy = TRUE),
    lag_malaria_deaths  = spdep::lag.listw(lw_queen_W, malaria_deaths, zero.policy = TRUE),
    lag_scd_cases       = spdep::lag.listw(lw_queen_W, scd_cases, zero.policy = TRUE),
    lag_scd_deaths      = spdep::lag.listw(lw_queen_W, scd_deaths, zero.policy = TRUE),

    lag_z_malaria_cases   = spdep::lag.listw(lw_queen_W, z_malaria_cases, zero.policy = TRUE),
    lag_z_malaria_deaths  = spdep::lag.listw(lw_queen_W, z_malaria_deaths, zero.policy = TRUE),
    lag_z_scd_cases       = spdep::lag.listw(lw_queen_W, z_scd_cases, zero.policy = TRUE),
    lag_z_scd_deaths      = spdep::lag.listw(lw_queen_W, z_scd_deaths, zero.policy = TRUE)
  )

# ----------------------------
# 15. GLOBAL MORAN'S I
# ----------------------------
moran_results <- bind_rows(
  run_moran_suite(q2_data$malaria_cases,  lw_queen_W, "malaria_cases"),
  run_moran_suite(q2_data$malaria_deaths, lw_queen_W, "malaria_deaths"),
  run_moran_suite(q2_data$scd_cases,      lw_queen_W, "scd_cases"),
  run_moran_suite(q2_data$scd_deaths,     lw_queen_W, "scd_deaths"),

  run_moran_suite(q2_data$log_malaria_cases,  lw_queen_W, "log_malaria_cases"),
  run_moran_suite(q2_data$log_malaria_deaths, lw_queen_W, "log_malaria_deaths"),
  run_moran_suite(q2_data$log_scd_cases,      lw_queen_W, "log_scd_cases"),
  run_moran_suite(q2_data$log_scd_deaths,     lw_queen_W, "log_scd_deaths")
)

write_csv(moran_results, file.path(table_dir, "q2_global_morans_i_results.csv"))

cat("Global Moran's I complete.\n", file = log_file, append = TRUE)

# ----------------------------
# 16. MORAN SCATTERPLOTS
# ----------------------------
make_moran_scatter <- function(data, xvar, lagvar, title_text, file_name) {
  p <- ggplot(data, aes(x = .data[[xvar]], y = .data[[lagvar]])) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = title_text,
      x = xvar,
      y = paste0("Spatial lag of ", xvar)
    ) +
    theme_minimal()

  ggsave(
    filename = file.path(fig_dir, file_name),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

make_moran_scatter(q2_data, "z_malaria_cases",  "lag_z_malaria_cases",
                   "Moran scatterplot: malaria cases", "q2_moran_scatter_malaria_cases.png")

make_moran_scatter(q2_data, "z_malaria_deaths", "lag_z_malaria_deaths",
                   "Moran scatterplot: malaria deaths", "q2_moran_scatter_malaria_deaths.png")

make_moran_scatter(q2_data, "z_scd_cases",      "lag_z_scd_cases",
                   "Moran scatterplot: SCD cases", "q2_moran_scatter_scd_cases.png")

make_moran_scatter(q2_data, "z_scd_deaths",     "lag_z_scd_deaths",
                   "Moran scatterplot: SCD deaths", "q2_moran_scatter_scd_deaths.png")

# ----------------------------
# 17. QC MAP OF NEIGHBOUR COUNTS
# ----------------------------
tmap_mode("plot")

map_nb_count <- tm_shape(q2_data) +
  tm_polygons("area_index", border.col = "grey30", lwd = 0.5) +
  tm_shape(q2_data) +
  tm_fill(col = "area_index", alpha = 0) +
  tm_borders() +
  tm_layout(main.title = "QC map: catchment polygons")

# Better thematic neighbour-count map
q2_data <- q2_data %>%
  left_join(
    nb_detail %>% select(catchment_id, neighbour_count, component_id, is_island),
    by = "catchment_id"
  )

map_neighbour_count <- tm_shape(q2_data) +
  tm_polygons("neighbour_count", title = "Neighbour count") +
  tm_layout(main.title = "Queen contiguity neighbour counts")

tmap_save(
  map_neighbour_count,
  filename = file.path(fig_dir, "q2_queen_neighbour_count_map.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# Optional island map
if (n_islands > 0) {
  map_islands <- tm_shape(q2_data) +
    tm_polygons("is_island", title = "Island") +
    tm_layout(main.title = "Queen contiguity islands")

  tmap_save(
    map_islands,
    filename = file.path(fig_dir, "q2_queen_islands_map.png"),
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ----------------------------
# 18. EXPORT SPATIAL OBJECT WITH LAGS
# ----------------------------
sf::st_write(
  q2_data,
  dsn = file.path(gpkg_dir, "q2_data_step2_with_weights_and_lags.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

saveRDS(q2_data, file.path(rds_dir, "q2_data_step2_with_weights_and_lags.rds"))

# ----------------------------
# 19. FINAL LOG
# ----------------------------
cat("Q2 Step 2 completed successfully: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE)

message("Step 2 completed successfully.")
What this revised Step 2 now does
It is fully aligned with your revised Step 1 and your actual folder structure. It will:
•	load q2_data_step1_cleaned.rds
•	sort catchments by catchment_id
•	create a queen contiguity neighbour list
•	identify islands and connected components
•	create:
•	q2_nb_queen.rds
•	q2_listw_queen_W.rds
•	q2_listw_queen_B.rds
•	compute spatial lags
•	run global Moran’s I for:
•	malaria_cases
•	malaria_deaths
•	scd_cases
•	scd_deaths
•	and their log-transformed versions
•	export the results table and Step 2 spatial object
Main output files
You should expect these key outputs:
•	outputs/tables/q2_queen_neighbour_summary.csv
•	outputs/tables/q2_queen_neighbour_detail.csv
•	outputs/tables/q2_global_morans_i_results.csv
•	outputs/rds/q2_nb_queen.rds
•	outputs/rds/q2_listw_queen_W.rds
•	outputs/rds/q2_data_step2_with_weights_and_lags.rds
Important technical note
Because your map clearly includes detached polygons / islands, I kept zero.policy = TRUE throughout. That is the correct choice here, otherwise Moran calculations can fail for areas with no neighbours.
One methodological note
Your Step 1 denominator uses u5_pop_utm.tif, which appears to be an under-5 population raster. That does not affect Step 2, because Step 2 is based on observed counts and spatial autocorrelation. It will matter later in the Bayesian shared-component model if you intend to use population offsets.
Send Step 3 and I’ll revise it in the same way.


step 2 R script q3.docx
Document
great job next step 3


Thought for a couple of seconds
Excellent. I’ve revised Step 3 so it fits the project structure you are now using, the Step 1 and Step 2 outputs we aligned earlier, and the logic in your draft for global bivariate Moran’s I, local bivariate LISA, and multivariate Geary’s C. 
step 2 R script q3
This revised version fixes the main alignment issues in the draft:
•	updates root_dir to C:/SCD-Malaria syndemic
•	uses the actual Step 2 output names:
•	q2_data_step2_with_weights_and_lags.rds
•	q2_listw_queen_W.rds
•	q2_nb_queen.rds
•	removes the unused old derived_dir logic
•	uses area_index instead of the old row_id
•	saves outputs into your current outputs/ structure
•	keeps zero.policy = TRUE throughout, which is important for detached catchments
Save this as:
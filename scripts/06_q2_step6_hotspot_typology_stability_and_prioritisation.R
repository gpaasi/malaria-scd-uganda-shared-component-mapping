# ============================================================
# OBJECTIVE 2, QUESTION 2
# Step 6: Hotspot typology, stability, discordance,
# and prioritisation
# ============================================================

# ----------------------------
# 0. GLOBAL OPTIONS
# ----------------------------
options(stringsAsFactors = FALSE)
options(scipen = 999)

# ----------------------------
# 1. LOAD PACKAGES
# ----------------------------
required_pkgs <- c(
  "sf",
  "dplyr",
  "readr",
  "ggplot2",
  "tmap",
  "janitor",
  "tidyr",
  "purrr",
  "tibble",
  "stringr"
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
  stop("Root directory does not exist. Update 'root_dir' in the script.")
}

setwd(root_dir)

# ----------------------------
# 3. DEFINE FOLDERS
# ----------------------------
output_dir <- file.path(root_dir, "outputs")
fig_dir    <- file.path(output_dir, "figures")
table_dir  <- file.path(output_dir, "tables")
rds_dir    <- file.path(output_dir, "rds")
gpkg_dir   <- file.path(output_dir, "gpkg")
log_dir    <- file.path(output_dir, "logs")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gpkg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 4. START LOG
# ----------------------------
log_file <- file.path(log_dir, "q2_step6_hotspot_typology_stability_log.txt")
cat("Step 6 started:", as.character(Sys.time()), "\n", file = log_file)
cat("Root directory:", root_dir, "\n", file = log_file, append = TRUE)

# ----------------------------
# 5. LOAD INPUTS
# ----------------------------
step3_file <- file.path(rds_dir, "q2_step3_bivariate_outputs.rds")
step4_file <- file.path(rds_dir, "q2_step4_bayesian_shared_outputs.rds")
step5_file <- file.path(rds_dir, "q2_step5_extended_model_outputs.rds")

if (!file.exists(step3_file)) stop("Missing q2_step3_bivariate_outputs.rds. Run Step 3 first.")
if (!file.exists(step4_file)) stop("Missing q2_step4_bayesian_shared_outputs.rds. Run Step 4 first.")
if (!file.exists(step5_file)) stop("Missing q2_step5_extended_model_outputs.rds. Run Step 5 first.")

step3 <- readRDS(step3_file)
step4 <- readRDS(step4_file)
step5 <- readRDS(step5_file)

if (!inherits(step3, "sf")) stop("Step 3 object is not sf.")
if (!inherits(step4, "sf")) stop("Step 4 object is not sf.")
if (!inherits(step5, "sf")) stop("Step 5 object is not sf.")

cat("Loaded Step 3, Step 4, and Step 5 outputs.\n", file = log_file, append = TRUE)

# ----------------------------
# 6. PREPARE BASE DATASET
# ----------------------------
required_cols_step3 <- c(
  "catchment_id", "area_index", "catchment_population",
  "malaria_cases", "scd_cases",
  "cases_cluster_005", "cases_cluster_010",
  "deaths_cluster_005", "deaths_cluster_010",
  "cases_p_value", "deaths_p_value"
)

missing_step3 <- setdiff(required_cols_step3, names(step3))
if (length(missing_step3) > 0) {
  stop("Missing Step 3 columns: ", paste(missing_step3, collapse = ", "))
}

base_sf <- step3 %>%
  select(
    catchment_id,
    area_index,
    catchment_population,
    malaria_cases,
    scd_cases,
    malaria_deaths,
    scd_deaths,
    cases_cluster_005,
    cases_cluster_010,
    deaths_cluster_005,
    deaths_cluster_010,
    cases_p_value,
    deaths_p_value,
    geometry
  ) %>%
  arrange(area_index)

# ----------------------------
# 7. PREPARE STEP 4 HOTSPOT DATA
# ----------------------------
required_cols_step4 <- c(
  "area_index",
  "shared_mean",
  "shared_sd",
  "exceed_prob_rr_gt1",
  "hotspot_080",
  "hotspot_090",
  "robust_hotspot_005",
  "robust_hotspot_010",
  "hotspot_class_005",
  "hotspot_class_010",
  "rr_malaria",
  "rr_scd"
)

missing_step4 <- setdiff(required_cols_step4, names(step4))
if (length(missing_step4) > 0) {
  stop("Missing Step 4 columns: ", paste(missing_step4, collapse = ", "))
}

step4_tab <- step4 %>%
  st_drop_geometry() %>%
  select(
    area_index,
    shared_mean_step4 = shared_mean,
    shared_sd_step4 = shared_sd,
    exceed_prob_step4 = exceed_prob_rr_gt1,
    hotspot080_step4 = hotspot_080,
    hotspot090_step4 = hotspot_090,
    robust005_step4 = robust_hotspot_005,
    robust010_step4 = robust_hotspot_010,
    hotspot_class005_step4 = hotspot_class_005,
    hotspot_class010_step4 = hotspot_class_010,
    rr_malaria_step4 = rr_malaria,
    rr_scd_step4 = rr_scd
  )

# ----------------------------
# 8. PREPARE STEP 5 HOTSPOT DATA
# ----------------------------
required_cols_step5 <- c(
  "area_index",
  "shared_mean",
  "shared_sd",
  "exceed_prob_rr_gt1",
  "hotspot_080",
  "hotspot_090",
  "robust_hotspot_005",
  "robust_hotspot_010",
  "hotspot_class_005",
  "hotspot_class_010",
  "malaria_spatial_mean",
  "scd_spatial_mean",
  "rr_malaria",
  "rr_scd"
)

missing_step5 <- setdiff(required_cols_step5, names(step5))
if (length(missing_step5) > 0) {
  stop("Missing Step 5 columns: ", paste(missing_step5, collapse = ", "))
}

step5_tab <- step5 %>%
  st_drop_geometry() %>%
  select(
    area_index,
    shared_mean_step5 = shared_mean,
    shared_sd_step5 = shared_sd,
    exceed_prob_step5 = exceed_prob_rr_gt1,
    hotspot080_step5 = hotspot_080,
    hotspot090_step5 = hotspot_090,
    robust005_step5 = robust_hotspot_005,
    robust010_step5 = robust_hotspot_010,
    hotspot_class005_step5 = hotspot_class_005,
    hotspot_class010_step5 = hotspot_class_010,
    malaria_spatial_mean,
    scd_spatial_mean,
    rr_malaria_step5 = rr_malaria,
    rr_scd_step5 = rr_scd
  )

# ----------------------------
# 9. MERGE STEP 3, 4, 5
# ----------------------------
q2_typology <- base_sf %>%
  left_join(step4_tab, by = "area_index") %>%
  left_join(step5_tab, by = "area_index") %>%
  arrange(area_index)

# ----------------------------
# 10. CREATE INDICATOR VARIABLES
# ----------------------------
q2_typology <- q2_typology %>%
  mutate(
    hh_cases_005 = ifelse(cases_cluster_005 == "High-High", 1L, 0L),
    hh_cases_010 = ifelse(cases_cluster_010 == "High-High", 1L, 0L),

    lh_cases_005 = ifelse(cases_cluster_005 == "Low-High", 1L, 0L),
    lh_cases_010 = ifelse(cases_cluster_010 == "Low-High", 1L, 0L),

    hl_cases_005 = ifelse(cases_cluster_005 == "High-Low", 1L, 0L),
    hl_cases_010 = ifelse(cases_cluster_010 == "High-Low", 1L, 0L),

    ll_cases_005 = ifelse(cases_cluster_005 == "Low-Low", 1L, 0L),
    ll_cases_010 = ifelse(cases_cluster_010 == "Low-Low", 1L, 0L)
  )

# ----------------------------
# 11. CREATE FINAL TYPOLOGY, STRICT
# ----------------------------
q2_typology <- q2_typology %>%
  mutate(
    final_typology_005 = case_when(
      hh_cases_005 == 1 & hotspot080_step5 == 1 ~ "Robust shared hotspot",
      hh_cases_005 == 0 & hotspot080_step5 == 1 ~ "Bayesian-only shared hotspot",
      hh_cases_005 == 1 & hotspot080_step5 == 0 ~ "LISA-only hotspot",
      lh_cases_005 == 1 | hl_cases_005 == 1 ~ "Discordant cluster area",
      ll_cases_005 == 1 ~ "Low-Low cluster area",
      TRUE ~ "No strong hotspot evidence"
    ),

    final_typology_010 = case_when(
      hh_cases_010 == 1 & hotspot080_step5 == 1 ~ "Robust shared hotspot",
      hh_cases_010 == 0 & hotspot080_step5 == 1 ~ "Bayesian-only shared hotspot",
      hh_cases_010 == 1 & hotspot080_step5 == 0 ~ "LISA-only hotspot",
      lh_cases_010 == 1 | hl_cases_010 == 1 ~ "Discordant cluster area",
      ll_cases_010 == 1 ~ "Low-Low cluster area",
      TRUE ~ "No strong hotspot evidence"
    )
  )

# ----------------------------
# 12. CREATE HOTSPOT STABILITY SCORE
# ----------------------------
q2_typology <- q2_typology %>%
  mutate(
    stability_score =
      hh_cases_005 +
      hh_cases_010 +
      ifelse(hotspot080_step4 == 1, 1L, 0L) +
      ifelse(hotspot090_step4 == 1, 1L, 0L) +
      ifelse(hotspot080_step5 == 1, 1L, 0L) +
      ifelse(hotspot090_step5 == 1, 1L, 0L) +
      ifelse(robust005_step4 == 1, 1L, 0L) +
      ifelse(robust010_step4 == 1, 1L, 0L) +
      ifelse(robust005_step5 == 1, 1L, 0L) +
      ifelse(robust010_step5 == 1, 1L, 0L),

    stability_class = case_when(
      stability_score >= 7 ~ "Highly stable",
      stability_score >= 4 ~ "Moderately stable",
      stability_score >= 2 ~ "Borderline stable",
      TRUE ~ "Unstable / weak evidence"
    )
  )

# ----------------------------
# 13. CREATE PRIORITISATION SCORE
# ----------------------------
# Combines model support, shared effect, and burden.
# Scaled components are used to avoid one variable dominating.

scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (any(is.na(rng)) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

q2_typology <- q2_typology %>%
  mutate(
    s_exceed = scale01(exceed_prob_step5),
    s_shared = scale01(shared_mean_step5),
    s_malaria = scale01(malaria_cases),
    s_scd = scale01(scd_cases),
    s_pop = scale01(catchment_population),

    prioritisation_score =
      (0.35 * s_exceed) +
      (0.20 * s_shared) +
      (0.15 * s_malaria) +
      (0.15 * s_scd) +
      (0.10 * s_pop) +
      (0.05 * scale01(stability_score))
  ) %>%
  arrange(desc(prioritisation_score)) %>%
  mutate(priority_rank = row_number())

# ----------------------------
# 14. CREATE SUMMARY TABLES
# ----------------------------
typology_counts_005 <- q2_typology %>%
  st_drop_geometry() %>%
  count(final_typology_005, sort = TRUE)

typology_counts_010 <- q2_typology %>%
  st_drop_geometry() %>%
  count(final_typology_010, sort = TRUE)

stability_counts <- q2_typology %>%
  st_drop_geometry() %>%
  count(stability_class, sort = TRUE)

write_csv(typology_counts_005, file.path(table_dir, "q2_step6_typology_counts_p005.csv"))
write_csv(typology_counts_010, file.path(table_dir, "q2_step6_typology_counts_p010.csv"))
write_csv(stability_counts, file.path(table_dir, "q2_step6_stability_counts.csv"))

# ----------------------------
# 15. CREATE SHORTLIST TABLES
# ----------------------------
top_robust_hotspots <- q2_typology %>%
  st_drop_geometry() %>%
  filter(final_typology_010 == "Robust shared hotspot") %>%
  arrange(desc(prioritisation_score), desc(exceed_prob_step5)) %>%
  select(
    priority_rank,
    catchment_id,
    catchment_population,
    malaria_cases,
    scd_cases,
    shared_mean_step5,
    shared_sd_step5,
    exceed_prob_step5,
    stability_score,
    stability_class,
    prioritisation_score,
    final_typology_005,
    final_typology_010
  )

top_bayesian_only <- q2_typology %>%
  st_drop_geometry() %>%
  filter(final_typology_010 == "Bayesian-only shared hotspot") %>%
  arrange(desc(prioritisation_score), desc(exceed_prob_step5)) %>%
  select(
    priority_rank,
    catchment_id,
    catchment_population,
    malaria_cases,
    scd_cases,
    shared_mean_step5,
    shared_sd_step5,
    exceed_prob_step5,
    stability_score,
    stability_class,
    prioritisation_score,
    final_typology_005,
    final_typology_010
  )

top_discordant <- q2_typology %>%
  st_drop_geometry() %>%
  filter(final_typology_010 == "Discordant cluster area") %>%
  arrange(desc(prioritisation_score), desc(exceed_prob_step5)) %>%
  select(
    priority_rank,
    catchment_id,
    catchment_population,
    malaria_cases,
    scd_cases,
    cases_cluster_005,
    cases_cluster_010,
    shared_mean_step5,
    exceed_prob_step5,
    stability_score,
    stability_class,
    prioritisation_score
  )

write_csv(top_robust_hotspots, file.path(table_dir, "q2_step6_top_robust_hotspots.csv"))
write_csv(top_bayesian_only, file.path(table_dir, "q2_step6_top_bayesian_only_hotspots.csv"))
write_csv(top_discordant, file.path(table_dir, "q2_step6_top_discordant_areas.csv"))

# ----------------------------
# 16. CREATE MASTER PRIORITISATION TABLE
# ----------------------------
priority_table <- q2_typology %>%
  st_drop_geometry() %>%
  select(
    priority_rank,
    catchment_id,
    catchment_population,
    malaria_cases,
    scd_cases,
    malaria_deaths,
    scd_deaths,
    cases_cluster_005,
    cases_cluster_010,
    shared_mean_step4,
    exceed_prob_step4,
    shared_mean_step5,
    shared_sd_step5,
    exceed_prob_step5,
    hotspot080_step4,
    hotspot090_step4,
    hotspot080_step5,
    hotspot090_step5,
    robust005_step4,
    robust010_step4,
    robust005_step5,
    robust010_step5,
    stability_score,
    stability_class,
    prioritisation_score,
    final_typology_005,
    final_typology_010,
    rr_malaria_step4,
    rr_scd_step4,
    rr_malaria_step5,
    rr_scd_step5
  ) %>%
  arrange(priority_rank)

write_csv(priority_table, file.path(table_dir, "q2_step6_master_prioritisation_table.csv"))

# ----------------------------
# 17. MAPS
# ----------------------------
tmap_mode("plot")

typology_palette <- c(
  "Robust shared hotspot" = "#b2182b",
  "Bayesian-only shared hotspot" = "#ef8a62",
  "LISA-only hotspot" = "#67a9cf",
  "Discordant cluster area" = "#fddbc7",
  "Low-Low cluster area" = "#2166ac",
  "No strong hotspot evidence" = "grey85"
)

stability_palette <- c(
  "Highly stable" = "#762a83",
  "Moderately stable" = "#af8dc3",
  "Borderline stable" = "#d9f0d3",
  "Unstable / weak evidence" = "grey85"
)

map_typology_005 <- tm_shape(q2_typology) +
  tm_polygons(
    "final_typology_005",
    title = "Typology (p ≤ 0.05)",
    palette = typology_palette
  ) +
  tm_layout(
    main.title = "Hotspot typology, strict threshold",
    legend.outside = TRUE
  )

map_typology_010 <- tm_shape(q2_typology) +
  tm_polygons(
    "final_typology_010",
    title = "Typology (p ≤ 0.10)",
    palette = typology_palette
  ) +
  tm_layout(
    main.title = "Hotspot typology, exploratory threshold",
    legend.outside = TRUE
  )

map_stability <- tm_shape(q2_typology) +
  tm_polygons(
    "stability_class",
    title = "Stability class",
    palette = stability_palette
  ) +
  tm_layout(
    main.title = "Hotspot stability classification",
    legend.outside = TRUE
  )

map_priority <- tm_shape(q2_typology) +
  tm_polygons(
    "prioritisation_score",
    title = "Prioritisation score",
    palette = "YlOrRd",
    style = "cont"
  ) +
  tm_layout(
    main.title = "Catchment prioritisation score",
    legend.outside = TRUE
  )

tmap_save(
  map_typology_005,
  filename = file.path(fig_dir, "q2_step6_hotspot_typology_map_p005.png"),
  width = 9, height = 6, dpi = 300
)

tmap_save(
  map_typology_010,
  filename = file.path(fig_dir, "q2_step6_hotspot_typology_map_p010.png"),
  width = 9, height = 6, dpi = 300
)

tmap_save(
  map_stability,
  filename = file.path(fig_dir, "q2_step6_hotspot_stability_map.png"),
  width = 9, height = 6, dpi = 300
)

tmap_save(
  map_priority,
  filename = file.path(fig_dir, "q2_step6_prioritisation_map.png"),
  width = 9, height = 6, dpi = 300
)

# ----------------------------
# 18. SAVE FINAL SPATIAL OBJECT
# ----------------------------
sf::st_write(
  q2_typology,
  file.path(gpkg_dir, "q2_step6_hotspot_typology_outputs.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

saveRDS(q2_typology, file.path(rds_dir, "q2_step6_hotspot_typology_outputs.rds"))

# ----------------------------
# 19. MANUSCRIPT SUMMARY TABLE
# ----------------------------
manuscript_summary <- tibble(
  metric = c(
    "Robust shared hotspots, strict threshold",
    "Robust shared hotspots, exploratory threshold",
    "Bayesian-only shared hotspots, strict threshold",
    "Bayesian-only shared hotspots, exploratory threshold",
    "LISA-only hotspots, strict threshold",
    "LISA-only hotspots, exploratory threshold",
    "Discordant cluster areas, strict threshold",
    "Discordant cluster areas, exploratory threshold",
    "Low-Low cluster areas, strict threshold",
    "Low-Low cluster areas, exploratory threshold"
  ),
  value = c(
    sum(q2_typology$final_typology_005 == "Robust shared hotspot", na.rm = TRUE),
    sum(q2_typology$final_typology_010 == "Robust shared hotspot", na.rm = TRUE),
    sum(q2_typology$final_typology_005 == "Bayesian-only shared hotspot", na.rm = TRUE),
    sum(q2_typology$final_typology_010 == "Bayesian-only shared hotspot", na.rm = TRUE),
    sum(q2_typology$final_typology_005 == "LISA-only hotspot", na.rm = TRUE),
    sum(q2_typology$final_typology_010 == "LISA-only hotspot", na.rm = TRUE),
    sum(q2_typology$final_typology_005 == "Discordant cluster area", na.rm = TRUE),
    sum(q2_typology$final_typology_010 == "Discordant cluster area", na.rm = TRUE),
    sum(q2_typology$final_typology_005 == "Low-Low cluster area", na.rm = TRUE),
    sum(q2_typology$final_typology_010 == "Low-Low cluster area", na.rm = TRUE)
  )
)

write_csv(manuscript_summary, file.path(table_dir, "q2_step6_manuscript_summary.csv"))

# ----------------------------
# 20. FINAL LOG
# ----------------------------
cat("\nTypology counts, strict threshold:\n", file = log_file, append = TRUE)
capture.output(print(typology_counts_005), file = log_file, append = TRUE)

cat("\nTypology counts, exploratory threshold:\n", file = log_file, append = TRUE)
capture.output(print(typology_counts_010), file = log_file, append = TRUE)

cat("\nStability counts:\n", file = log_file, append = TRUE)
capture.output(print(stability_counts), file = log_file, append = TRUE)

cat("\nStep 6 completed: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE)

message("Step 6 completed successfully.")
message("Saved:")
message("- q2_step6_typology_counts_p005.csv")
message("- q2_step6_typology_counts_p010.csv")
message("- q2_step6_stability_counts.csv")
message("- q2_step6_top_robust_hotspots.csv")
message("- q2_step6_top_bayesian_only_hotspots.csv")
message("- q2_step6_top_discordant_areas.csv")
message("- q2_step6_master_prioritisation_table.csv")
message("- q2_step6_hotspot_typology_outputs.rds")

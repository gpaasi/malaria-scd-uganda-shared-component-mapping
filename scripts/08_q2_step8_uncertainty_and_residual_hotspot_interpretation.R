# ============================================================
# OBJECTIVE 2, QUESTION 2
# Step 8: Uncertainty, residual hotspot confidence,
# and adjusted hotspot interpretation
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
log_file <- file.path(log_dir, "q2_step8_uncertainty_residual_hotspot_log.txt")
cat("Step 8 started:", as.character(Sys.time()), "\n", file = log_file)
cat("Root directory:", root_dir, "\n", file = log_file, append = TRUE)

# ----------------------------
# 5. LOAD INPUTS
# ----------------------------
step6_file <- file.path(rds_dir, "q2_step6_hotspot_typology_outputs.rds")
step7_file <- file.path(rds_dir, "q2_step7_adjusted_model_outputs.rds")

if (!file.exists(step6_file)) stop("Missing q2_step6_hotspot_typology_outputs.rds. Run Step 6 first.")
if (!file.exists(step7_file)) stop("Missing q2_step7_adjusted_model_outputs.rds. Run Step 7 first.")

step6 <- readRDS(step6_file)
step7 <- readRDS(step7_file)

if (!inherits(step6, "sf")) stop("Step 6 object is not sf.")
if (!inherits(step7, "sf")) stop("Step 7 object is not sf.")

cat("Loaded Step 6 and Step 7 outputs.\n", file = log_file, append = TRUE)

# ----------------------------
# 6. SELECT AND MERGE RELEVANT DATA
# ----------------------------
needed_step6 <- c(
  "catchment_id",
  "area_index",
  "stability_score",
  "stability_class",
  "final_typology_005",
  "final_typology_010"
)

needed_step7 <- c(
  "catchment_id",
  "area_index",
  "facility_name",
  "facility_class",
  "ownership_class",
  "region",
  "catchment_population",
  "malaria_cases",
  "scd_cases",
  "shared_mean_adj",
  "shared_sd_adj",
  "shared_q025_adj",
  "shared_q975_adj",
  "exceed_prob_rr_gt1_adj",
  "hotspot_080_adj",
  "hotspot_090_adj",
  "robust_hotspot_005_adj",
  "robust_hotspot_010_adj",
  "hotspot_class_005_adj",
  "hotspot_class_010_adj",
  "rr_malaria_adj",
  "rr_scd_adj"
)

missing_step6 <- setdiff(needed_step6, names(step6))
missing_step7 <- setdiff(needed_step7, names(step7))

if (length(missing_step6) > 0) {
  stop("Missing Step 6 columns: ", paste(missing_step6, collapse = ", "))
}
if (length(missing_step7) > 0) {
  stop("Missing Step 7 columns: ", paste(missing_step7, collapse = ", "))
}

step6_tab <- step6 %>%
  st_drop_geometry() %>%
  select(all_of(needed_step6))

q2_uq <- step7 %>%
  select(all_of(needed_step7), geometry) %>%
  left_join(step6_tab, by = c("catchment_id", "area_index"))

# ----------------------------
# 7. CREATE UNCERTAINTY METRICS
# ----------------------------
q2_uq <- q2_uq %>%
  mutate(
    shared_ci_width_adj = shared_q975_adj - shared_q025_adj,
    threshold_distance_080 = abs(exceed_prob_rr_gt1_adj - 0.80),
    threshold_distance_090 = abs(exceed_prob_rr_gt1_adj - 0.90)
  )

# relative uncertainty summaries
sd_median <- median(q2_uq$shared_sd_adj, na.rm = TRUE)
sd_q75 <- quantile(q2_uq$shared_sd_adj, 0.75, na.rm = TRUE)

ci_median <- median(q2_uq$shared_ci_width_adj, na.rm = TRUE)
ci_q75 <- quantile(q2_uq$shared_ci_width_adj, 0.75, na.rm = TRUE)

# ----------------------------
# 8. CLASSIFY UNCERTAINTY LEVEL
# ----------------------------
q2_uq <- q2_uq %>%
  mutate(
    uncertainty_class = case_when(
      shared_sd_adj <= sd_median & shared_ci_width_adj <= ci_median ~ "Low uncertainty",
      shared_sd_adj > sd_q75 | shared_ci_width_adj > ci_q75 ~ "High uncertainty",
      TRUE ~ "Moderate uncertainty"
    )
  )

# ----------------------------
# 9. CREATE RESIDUAL HOTSPOT CONFIDENCE CLASSES
# ----------------------------
q2_uq <- q2_uq %>%
  mutate(
    residual_hotspot_confidence = case_when(
      hotspot_090_adj == 1 &
        uncertainty_class == "Low uncertainty" &
        stability_score >= 4 ~ "High-confidence residual hotspot",

      hotspot_080_adj == 1 &
        uncertainty_class != "High uncertainty" &
        stability_score >= 2 ~ "Moderate-confidence residual hotspot",

      hotspot_080_adj == 1 &
        (uncertainty_class == "High uncertainty" | stability_score < 2) ~ "Uncertain candidate hotspot",

      hotspot_080_adj == 0 &
        hotspot_class_010_adj == "LISA hotspot only" ~ "Exploratory local-only hotspot",

      TRUE ~ "Not residual hotspot"
    )
  )

# ----------------------------
# 10. CREATE RESIDUAL HOTSPOT TYPOLOGY
# ----------------------------
q2_uq <- q2_uq %>%
  mutate(
    residual_typology = case_when(
      residual_hotspot_confidence == "High-confidence residual hotspot" ~ "High-confidence residual hotspot",
      residual_hotspot_confidence == "Moderate-confidence residual hotspot" ~ "Moderate-confidence residual hotspot",
      residual_hotspot_confidence == "Uncertain candidate hotspot" ~ "Uncertain candidate hotspot",
      hotspot_class_010_adj == "LISA hotspot only" ~ "Exploratory local-only hotspot",
      final_typology_010 == "Discordant cluster area" ~ "Discordant cluster area",
      final_typology_010 == "Low-Low cluster area" ~ "Low-Low cluster area",
      TRUE ~ "No strong residual signal"
    )
  )

# ----------------------------
# 11. CREATE RESIDUAL PRIORITISATION SCORE
# ----------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (any(is.na(rng)) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

q2_uq <- q2_uq %>%
  mutate(
    s_exceed = scale01(exceed_prob_rr_gt1_adj),
    s_shared = scale01(shared_mean_adj),
    s_cases_mal = scale01(malaria_cases),
    s_cases_scd = scale01(scd_cases),
    s_stability = scale01(stability_score),
    s_uncertainty_penalty = 1 - scale01(shared_sd_adj),

    residual_priority_score =
      (0.30 * s_exceed) +
      (0.20 * s_shared) +
      (0.15 * s_cases_mal) +
      (0.15 * s_cases_scd) +
      (0.10 * s_stability) +
      (0.10 * s_uncertainty_penalty)
  ) %>%
  arrange(desc(residual_priority_score)) %>%
  mutate(residual_priority_rank = row_number())

# ----------------------------
# 12. SUMMARY TABLES
# ----------------------------
uncertainty_counts <- q2_uq %>%
  st_drop_geometry() %>%
  count(uncertainty_class, sort = TRUE)

confidence_counts <- q2_uq %>%
  st_drop_geometry() %>%
  count(residual_hotspot_confidence, sort = TRUE)

typology_counts <- q2_uq %>%
  st_drop_geometry() %>%
  count(residual_typology, sort = TRUE)

write_csv(uncertainty_counts, file.path(table_dir, "q2_step8_uncertainty_counts.csv"))
write_csv(confidence_counts, file.path(table_dir, "q2_step8_residual_hotspot_confidence_counts.csv"))
write_csv(typology_counts, file.path(table_dir, "q2_step8_residual_typology_counts.csv"))

# ----------------------------
# 13. TOP RESIDUAL HOTSPOTS
# ----------------------------
top_high_confidence <- q2_uq %>%
  st_drop_geometry() %>%
  filter(residual_hotspot_confidence == "High-confidence residual hotspot") %>%
  arrange(desc(residual_priority_score), desc(exceed_prob_rr_gt1_adj)) %>%
  select(
    residual_priority_rank,
    catchment_id,
    facility_name,
    facility_class,
    ownership_class,
    region,
    catchment_population,
    malaria_cases,
    scd_cases,
    shared_mean_adj,
    shared_sd_adj,
    shared_ci_width_adj,
    exceed_prob_rr_gt1_adj,
    stability_score,
    stability_class,
    rr_malaria_adj,
    rr_scd_adj,
    residual_hotspot_confidence
  )

top_moderate_confidence <- q2_uq %>%
  st_drop_geometry() %>%
  filter(residual_hotspot_confidence == "Moderate-confidence residual hotspot") %>%
  arrange(desc(residual_priority_score), desc(exceed_prob_rr_gt1_adj)) %>%
  select(
    residual_priority_rank,
    catchment_id,
    facility_name,
    facility_class,
    ownership_class,
    region,
    catchment_population,
    malaria_cases,
    scd_cases,
    shared_mean_adj,
    shared_sd_adj,
    shared_ci_width_adj,
    exceed_prob_rr_gt1_adj,
    stability_score,
    stability_class,
    rr_malaria_adj,
    rr_scd_adj,
    residual_hotspot_confidence
  )

top_uncertain_candidates <- q2_uq %>%
  st_drop_geometry() %>%
  filter(residual_hotspot_confidence == "Uncertain candidate hotspot") %>%
  arrange(desc(residual_priority_score), desc(exceed_prob_rr_gt1_adj)) %>%
  select(
    residual_priority_rank,
    catchment_id,
    facility_name,
    facility_class,
    ownership_class,
    region,
    catchment_population,
    malaria_cases,
    scd_cases,
    shared_mean_adj,
    shared_sd_adj,
    shared_ci_width_adj,
    exceed_prob_rr_gt1_adj,
    threshold_distance_080,
    stability_score,
    stability_class,
    residual_hotspot_confidence
  )

write_csv(top_high_confidence, file.path(table_dir, "q2_step8_top_high_confidence_residual_hotspots.csv"))
write_csv(top_moderate_confidence, file.path(table_dir, "q2_step8_top_moderate_confidence_residual_hotspots.csv"))
write_csv(top_uncertain_candidates, file.path(table_dir, "q2_step8_top_uncertain_candidate_hotspots.csv"))

# ----------------------------
# 14. MASTER INTERPRETATION TABLE
# ----------------------------
master_table <- q2_uq %>%
  st_drop_geometry() %>%
  select(
    residual_priority_rank,
    catchment_id,
    facility_name,
    facility_class,
    ownership_class,
    region,
    catchment_population,
    malaria_cases,
    scd_cases,
    shared_mean_adj,
    shared_sd_adj,
    shared_q025_adj,
    shared_q975_adj,
    shared_ci_width_adj,
    exceed_prob_rr_gt1_adj,
    hotspot_080_adj,
    hotspot_090_adj,
    robust_hotspot_005_adj,
    robust_hotspot_010_adj,
    hotspot_class_005_adj,
    hotspot_class_010_adj,
    stability_score,
    stability_class,
    uncertainty_class,
    residual_hotspot_confidence,
    residual_typology,
    residual_priority_score,
    rr_malaria_adj,
    rr_scd_adj
  ) %>%
  arrange(residual_priority_rank)

write_csv(master_table, file.path(table_dir, "q2_step8_master_uncertainty_interpretation_table.csv"))

# ----------------------------
# 15. MAPS
# ----------------------------
tmap_mode("plot")

uncertainty_palette <- c(
  "Low uncertainty" = "#1a9641",
  "Moderate uncertainty" = "#fdae61",
  "High uncertainty" = "#d7191c"
)

confidence_palette <- c(
  "High-confidence residual hotspot" = "#762a83",
  "Moderate-confidence residual hotspot" = "#af8dc3",
  "Uncertain candidate hotspot" = "#fdb863",
  "Exploratory local-only hotspot" = "#67a9cf",
  "Not residual hotspot" = "grey85"
)

typology_palette <- c(
  "High-confidence residual hotspot" = "#762a83",
  "Moderate-confidence residual hotspot" = "#af8dc3",
  "Uncertain candidate hotspot" = "#fdb863",
  "Exploratory local-only hotspot" = "#67a9cf",
  "Discordant cluster area" = "#fddbc7",
  "Low-Low cluster area" = "#2166ac",
  "No strong residual signal" = "grey85"
)

map_shared_sd <- tm_shape(q2_uq) +
  tm_polygons(
    "shared_sd_adj",
    title = "Shared effect SD",
    palette = "YlOrRd",
    style = "cont"
  ) +
  tm_layout(
    main.title = "Posterior SD of adjusted shared spatial effect",
    legend.outside = TRUE
  )

map_ci_width <- tm_shape(q2_uq) +
  tm_polygons(
    "shared_ci_width_adj",
    title = "CI width",
    palette = "YlOrRd",
    style = "cont"
  ) +
  tm_layout(
    main.title = "Credible interval width of adjusted shared spatial effect",
    legend.outside = TRUE
  )

map_uncertainty_class <- tm_shape(q2_uq) +
  tm_polygons(
    "uncertainty_class",
    title = "Uncertainty class",
    palette = uncertainty_palette
  ) +
  tm_layout(
    main.title = "Uncertainty classification of adjusted shared effect",
    legend.outside = TRUE
  )

map_confidence <- tm_shape(q2_uq) +
  tm_polygons(
    "residual_hotspot_confidence",
    title = "Residual hotspot confidence",
    palette = confidence_palette
  ) +
  tm_layout(
    main.title = "Confidence classes for adjusted residual hotspots",
    legend.outside = TRUE
  )

map_residual_typology <- tm_shape(q2_uq) +
  tm_polygons(
    "residual_typology",
    title = "Residual hotspot typology",
    palette = typology_palette
  ) +
  tm_layout(
    main.title = "Residual hotspot typology after covariate adjustment",
    legend.outside = TRUE
  )

map_priority <- tm_shape(q2_uq) +
  tm_polygons(
    "residual_priority_score",
    title = "Residual priority score",
    palette = "YlOrRd",
    style = "cont"
  ) +
  tm_layout(
    main.title = "Residual hotspot prioritisation score",
    legend.outside = TRUE
  )

tmap_save(map_shared_sd, file.path(fig_dir, "q2_step8_shared_sd_map.png"), width = 9, height = 6, dpi = 300)
tmap_save(map_ci_width, file.path(fig_dir, "q2_step8_shared_ci_width_map.png"), width = 9, height = 6, dpi = 300)
tmap_save(map_uncertainty_class, file.path(fig_dir, "q2_step8_uncertainty_class_map.png"), width = 9, height = 6, dpi = 300)
tmap_save(map_confidence, file.path(fig_dir, "q2_step8_residual_hotspot_confidence_map.png"), width = 9, height = 6, dpi = 300)
tmap_save(map_residual_typology, file.path(fig_dir, "q2_step8_residual_typology_map.png"), width = 9, height = 6, dpi = 300)
tmap_save(map_priority, file.path(fig_dir, "q2_step8_residual_priority_map.png"), width = 9, height = 6, dpi = 300)

# ----------------------------
# 16. SAVE SPATIAL OUTPUTS
# ----------------------------
sf::st_write(
  q2_uq,
  file.path(gpkg_dir, "q2_step8_uncertainty_residual_hotspot_outputs.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

saveRDS(q2_uq, file.path(rds_dir, "q2_step8_uncertainty_residual_hotspot_outputs.rds"))

# ----------------------------
# 17. MANUSCRIPT SUMMARY
# ----------------------------
manuscript_summary <- tibble(
  metric = c(
    "Adjusted hotspots (>0.80)",
    "Adjusted hotspots (>0.90)",
    "High-confidence residual hotspots",
    "Moderate-confidence residual hotspots",
    "Uncertain candidate hotspots",
    "Exploratory local-only hotspots",
    "High uncertainty catchments",
    "Low uncertainty catchments"
  ),
  value = c(
    sum(q2_uq$hotspot_080_adj, na.rm = TRUE),
    sum(q2_uq$hotspot_090_adj, na.rm = TRUE),
    sum(q2_uq$residual_hotspot_confidence == "High-confidence residual hotspot", na.rm = TRUE),
    sum(q2_uq$residual_hotspot_confidence == "Moderate-confidence residual hotspot", na.rm = TRUE),
    sum(q2_uq$residual_hotspot_confidence == "Uncertain candidate hotspot", na.rm = TRUE),
    sum(q2_uq$residual_hotspot_confidence == "Exploratory local-only hotspot", na.rm = TRUE),
    sum(q2_uq$uncertainty_class == "High uncertainty", na.rm = TRUE),
    sum(q2_uq$uncertainty_class == "Low uncertainty", na.rm = TRUE)
  )
)

write_csv(manuscript_summary, file.path(table_dir, "q2_step8_manuscript_summary.csv"))

# ----------------------------
# 18. FINAL LOG
# ----------------------------
cat("\nUncertainty counts:\n", file = log_file, append = TRUE)
capture.output(print(uncertainty_counts), file = log_file, append = TRUE)

cat("\nResidual hotspot confidence counts:\n", file = log_file, append = TRUE)
capture.output(print(confidence_counts), file = log_file, append = TRUE)

cat("\nResidual typology counts:\n", file = log_file, append = TRUE)
capture.output(print(typology_counts), file = log_file, append = TRUE)

cat("\nStep 8 completed: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE)

message("Step 8 completed successfully.")
message("Saved:")
message("- q2_step8_uncertainty_counts.csv")
message("- q2_step8_residual_hotspot_confidence_counts.csv")
message("- q2_step8_residual_typology_counts.csv")
message("- q2_step8_top_high_confidence_residual_hotspots.csv")
message("- q2_step8_top_moderate_confidence_residual_hotspots.csv")
message("- q2_step8_top_uncertain_candidate_hotspots.csv")
message("- q2_step8_master_uncertainty_interpretation_table.csv")
message("- q2_step8_uncertainty_residual_hotspot_outputs.rds")
message("- uncertainty, confidence, typology, and priority maps")

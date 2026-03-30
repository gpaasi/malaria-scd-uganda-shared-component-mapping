# Step 4: Bayesian shared-component spatial model using INLA
# Spatial co-clustering of malaria and SCD across 2SFCA catchments
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
  "spdep",
  "ggplot2",
  "tmap",
  "janitor",
  "tidyr",
  "purrr",
  "tibble"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_pkgs, installed)

if (length(to_install) > 0) {
  install.packages(to_install)
}

if (!"INLA" %in% rownames(installed.packages())) {
  install.packages(
    "INLA",
    repos = c(
      getOption("repos"),
      INLA = "https://inla.r-inla-download.org/R/stable"
    )
  )
}

invisible(lapply(required_pkgs, library, character.only = TRUE))
library(INLA)

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
graph_dir  <- file.path(output_dir, "graphs")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gpkg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(graph_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 4. START LOG
# ----------------------------
log_file <- file.path(log_dir, "q2_step4_bayesian_shared_component_inla_log.txt")
cat("Step 4 started:", as.character(Sys.time()), "\n", file = log_file)
cat("Root directory:", root_dir, "\n", file = log_file, append = TRUE)

# ----------------------------
# 5. LOAD STEP 3 OUTPUTS
# ----------------------------
q2_file <- file.path(rds_dir, "q2_step3_bivariate_outputs.rds")
nb_file <- file.path(rds_dir, "q2_nb_queen.rds")

if (!file.exists(q2_file)) stop("Missing q2_step3_bivariate_outputs.rds. Run Step 3 first.")
if (!file.exists(nb_file)) stop("Missing q2_nb_queen.rds. Run Step 2 first.")

q2 <- readRDS(q2_file)
nb_queen <- readRDS(nb_file)

if (!inherits(q2, "sf")) {
  stop("The Step 3 object is not an sf object.")
}

cat("Loaded analysis object with ", nrow(q2), " catchments.\n",
    file = log_file, append = TRUE)

# ----------------------------
# 6. BASIC CHECKS
# ----------------------------
required_cols <- c(
  "catchment_id",
  "area_index",
  "malaria_cases",
  "scd_cases"
)

missing_cols <- setdiff(required_cols, names(q2))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# enforce deterministic order
q2 <- q2 %>%
  arrange(area_index)

# ----------------------------
# 7. CREATE ADJACENCY GRAPH FOR INLA
# ----------------------------
adj_file <- file.path(graph_dir, "q2_queen_adj.graph")

spdep::nb2INLA(file = adj_file, nb = nb_queen)

if (!file.exists(adj_file)) {
  stop("Failed to create INLA adjacency graph file.")
}

cat("Adjacency graph saved: ", adj_file, "\n", file = log_file, append = TRUE)

# ----------------------------
# 8. PREPARE EXPECTED COUNTS / OFFSETS
# ----------------------------
# Pragmatic interim approach:
# use internally derived expected counts based on the overall mean count.
# This is acceptable for exploratory shared-component mapping when a valid
# all-age denominator is not yet available.
#
# If you later obtain an all-age catchment population surface, replace this
# section with population-based expected counts.

q2 <- q2 %>%
  mutate(
    expected_malaria_cases = mean(malaria_cases, na.rm = TRUE),
    expected_scd_cases     = mean(scd_cases, na.rm = TRUE)
  ) %>%
  mutate(
    expected_malaria_cases = ifelse(expected_malaria_cases <= 0, 1e-06, expected_malaria_cases),
    expected_scd_cases     = ifelse(expected_scd_cases <= 0, 1e-06, expected_scd_cases)
  )

# ----------------------------
# 9. CREATE STACKED LONG DATASET
# ----------------------------
dat_mal <- q2 %>%
  st_drop_geometry() %>%
  transmute(
    catchment_id,
    area_index,
    outcome = malaria_cases,
    expected = expected_malaria_cases,
    outcome_type = "malaria"
  )

dat_scd <- q2 %>%
  st_drop_geometry() %>%
  transmute(
    catchment_id,
    area_index,
    outcome = scd_cases,
    expected = expected_scd_cases,
    outcome_type = "scd"
  )

dat_long <- bind_rows(dat_mal, dat_scd) %>%
  mutate(
    obs_id = row_number(),
    outcome_type = factor(outcome_type, levels = c("malaria", "scd")),
    outcome_bin = ifelse(outcome_type == "scd", 1L, 0L),

    # shared spatial component
    shared_main = area_index,
    shared_scd  = ifelse(outcome_type == "scd", area_index, NA_integer_),

    # outcome-specific iid heterogeneity
    het_malaria = ifelse(outcome_type == "malaria", area_index, NA_integer_),
    het_scd     = ifelse(outcome_type == "scd", area_index, NA_integer_)
  )

# ----------------------------
# 10. PRIORS
# ----------------------------
hyper_bym2 <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)),
  phi  = list(prior = "pc", param = c(0.5, 2/3))
)

hyper_iid <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01))
)

# ----------------------------
# 11. DEFINE SHARED-COMPONENT MODEL
# ----------------------------
# malaria is the baseline
# scd gets a copied shared spatial component with an estimated loading

formula_shared <- outcome ~ 0 + outcome_type +
  f(shared_main,
    model = "bym2",
    graph = adj_file,
    scale.model = TRUE,
    constr = TRUE,
    hyper = hyper_bym2) +
  f(shared_scd,
    copy = "shared_main",
    fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 1)))) +
  f(het_malaria,
    model = "iid",
    hyper = hyper_iid) +
  f(het_scd,
    model = "iid",
    hyper = hyper_iid)

# ----------------------------
# 12. FIT MODEL
# ----------------------------
cat("Fitting INLA shared-component model...\n", file = log_file, append = TRUE)

model_shared <- INLA::inla(
  formula_shared,
  family = "poisson",
  data = dat_long,
  E = expected,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  verbose = FALSE
)

saveRDS(model_shared, file.path(rds_dir, "q2_inla_shared_component_model.rds"))

cat("Model fitting completed.\n", file = log_file, append = TRUE)
cat("DIC: ", model_shared$dic$dic, "\n", file = log_file, append = TRUE)
cat("WAIC: ", model_shared$waic$waic, "\n", file = log_file, append = TRUE)

# ----------------------------
# 13. EXTRACT FIXED EFFECTS AND HYPERPARAMETERS
# ----------------------------
fixed_effects <- model_shared$summary.fixed %>%
  tibble::rownames_to_column("term")

write_csv(fixed_effects, file.path(table_dir, "q2_inla_fixed_effects.csv"))

hyperpar_summary <- model_shared$summary.hyperpar %>%
  tibble::rownames_to_column("parameter")

write_csv(hyperpar_summary, file.path(table_dir, "q2_inla_hyperparameters.csv"))

# ----------------------------
# 14. EXTRACT SHARED SPATIAL COMPONENT
# ----------------------------
shared_spatial <- model_shared$summary.random$shared_main %>%
  mutate(
    area_index = ID,
    shared_mean = mean,
    shared_sd = sd,
    shared_q025 = `0.025quant`,
    shared_q50  = `0.5quant`,
    shared_q975 = `0.975quant`
  ) %>%
  select(area_index, shared_mean, shared_sd, shared_q025, shared_q50, shared_q975)

# ----------------------------
# 15. COMPUTE POSTERIOR EXCEEDANCE PROBABILITIES
# ----------------------------
# Exceedance of RR > 1 corresponds to shared effect > 0

shared_spatial <- shared_spatial %>%
  mutate(
    exceed_prob_rr_gt1 = 1 - pnorm(0, mean = shared_mean, sd = shared_sd),
    hotspot_080 = ifelse(exceed_prob_rr_gt1 > 0.80, 1L, 0L),
    hotspot_090 = ifelse(exceed_prob_rr_gt1 > 0.90, 1L, 0L)
  )

write_csv(shared_spatial, file.path(table_dir, "q2_shared_spatial_component_and_exceedance.csv"))

# ----------------------------
# 16. MERGE BACK TO SPATIAL DATA
# ----------------------------
q2_shared <- q2 %>%
  left_join(shared_spatial, by = "area_index")

# ----------------------------
# 17. EXTRACT FITTED VALUES
# ----------------------------
n_catch <- nrow(q2)

fitted_vals <- model_shared$summary.fitted.values %>%
  mutate(obs_id = row_number())

fitted_mal <- fitted_vals[1:n_catch, ] %>%
  transmute(
    area_index = q2$area_index,
    fitted_malaria_mean = mean,
    fitted_malaria_q025 = `0.025quant`,
    fitted_malaria_q975 = `0.975quant`
  )

fitted_scd <- fitted_vals[(n_catch + 1):(2 * n_catch), ] %>%
  transmute(
    area_index = q2$area_index,
    fitted_scd_mean = mean,
    fitted_scd_q025 = `0.025quant`,
    fitted_scd_q975 = `0.975quant`
  )

q2_shared <- q2_shared %>%
  left_join(fitted_mal, by = "area_index") %>%
  left_join(fitted_scd, by = "area_index") %>%
  mutate(
    rr_malaria = fitted_malaria_mean / expected_malaria_cases,
    rr_scd     = fitted_scd_mean / expected_scd_cases
  )

# ----------------------------
# 18. DEFINE HOTSPOT CONCORDANCE
# ----------------------------
# Main strict hotspot concordance uses p <= 0.05 High-High LISA
# Exploratory concordance uses p <= 0.10 High-High LISA

if (!"cases_cluster_005" %in% names(q2_shared)) {
  q2_shared$cases_cluster_005 <- NA_character_
}
if (!"cases_cluster_010" %in% names(q2_shared)) {
  q2_shared$cases_cluster_010 <- NA_character_
}

q2_shared <- q2_shared %>%
  mutate(
    robust_hotspot_005 = ifelse(
      cases_cluster_005 == "High-High" & exceed_prob_rr_gt1 > 0.80,
      1L, 0L
    ),
    robust_hotspot_010 = ifelse(
      cases_cluster_010 == "High-High" & exceed_prob_rr_gt1 > 0.80,
      1L, 0L
    )
  )

hotspot_summary <- q2_shared %>%
  st_drop_geometry() %>%
  summarise(
    n_catchments = n(),
    n_hotspot_080 = sum(hotspot_080, na.rm = TRUE),
    n_hotspot_090 = sum(hotspot_090, na.rm = TRUE),
    n_robust_hotspot_005 = sum(robust_hotspot_005, na.rm = TRUE),
    n_robust_hotspot_010 = sum(robust_hotspot_010, na.rm = TRUE)
  )

write_csv(hotspot_summary, file.path(table_dir, "q2_hotspot_summary.csv"))

# ----------------------------
# 19. MAP SHARED SPATIAL EFFECT
# ----------------------------
tmap_mode("plot")

map_shared_mean <- tm_shape(q2_shared) +
  tm_polygons(
    "shared_mean",
    title = "Shared latent effect",
    palette = "-RdBu",
    style = "cont"
  ) +
  tm_layout(
    main.title = "Posterior mean shared latent spatial effect",
    legend.outside = TRUE
  )

tmap_save(
  map_shared_mean,
  filename = file.path(fig_dir, "q2_shared_latent_spatial_effect_map.png"),
  width = 9, height = 6, dpi = 300
)

# ----------------------------
# 20. MAP EXCEEDANCE PROBABILITY
# ----------------------------
map_exceed <- tm_shape(q2_shared) +
  tm_polygons(
    "exceed_prob_rr_gt1",
    title = "P(shared RR > 1)",
    palette = "YlOrRd",
    style = "cont"
  ) +
  tm_layout(
    main.title = "Posterior exceedance probability map",
    legend.outside = TRUE
  )

tmap_save(
  map_exceed,
  filename = file.path(fig_dir, "q2_shared_exceedance_probability_map.png"),
  width = 9, height = 6, dpi = 300
)

# ----------------------------
# 21. MAP HOTSPOT CONCORDANCE
# ----------------------------
q2_shared <- q2_shared %>%
  mutate(
    hotspot_class_005 = case_when(
      robust_hotspot_005 == 1 ~ "Robust hotspot",
      hotspot_080 == 1 ~ "Bayesian hotspot only",
      cases_cluster_005 == "High-High" ~ "LISA hotspot only",
      TRUE ~ "Not hotspot"
    ),
    hotspot_class_010 = case_when(
      robust_hotspot_010 == 1 ~ "Robust hotspot",
      hotspot_080 == 1 ~ "Bayesian hotspot only",
      cases_cluster_010 == "High-High" ~ "LISA hotspot only",
      TRUE ~ "Not hotspot"
    )
  )

hotspot_palette <- c(
  "Robust hotspot" = "#b2182b",
  "Bayesian hotspot only" = "#ef8a62",
  "LISA hotspot only" = "#67a9cf",
  "Not hotspot" = "grey85"
)

map_hotspots_005 <- tm_shape(q2_shared) +
  tm_polygons(
    "hotspot_class_005",
    title = "Hotspot class (p ≤ 0.05)",
    palette = hotspot_palette
  ) +
  tm_layout(
    main.title = "Malaria-SCD hotspot concordance (strict LISA threshold)",
    legend.outside = TRUE
  )

map_hotspots_010 <- tm_shape(q2_shared) +
  tm_polygons(
    "hotspot_class_010",
    title = "Hotspot class (p ≤ 0.10)",
    palette = hotspot_palette
  ) +
  tm_layout(
    main.title = "Malaria-SCD hotspot concordance (exploratory LISA threshold)",
    legend.outside = TRUE
  )

tmap_save(
  map_hotspots_005,
  filename = file.path(fig_dir, "q2_hotspot_concordance_map_p005.png"),
  width = 9, height = 6, dpi = 300
)

tmap_save(
  map_hotspots_010,
  filename = file.path(fig_dir, "q2_hotspot_concordance_map_p010.png"),
  width = 9, height = 6, dpi = 300
)

# ----------------------------
# 22. EXPORT MAIN TABLE
# ----------------------------
main_export <- q2_shared %>%
  st_drop_geometry() %>%
  select(
    catchment_id,
    area_index,
    malaria_cases,
    scd_cases,
    expected_malaria_cases,
    expected_scd_cases,
    cases_cluster_005,
    cases_cluster_010,
    shared_mean,
    shared_sd,
    exceed_prob_rr_gt1,
    hotspot_080,
    hotspot_090,
    robust_hotspot_005,
    robust_hotspot_010,
    hotspot_class_005,
    hotspot_class_010,
    rr_malaria,
    rr_scd
  )

write_csv(main_export, file.path(table_dir, "q2_bayesian_shared_component_outputs.csv"))

# ----------------------------
# 23. SAVE SPATIAL OUTPUTS
# ----------------------------
sf::st_write(
  q2_shared,
  file.path(gpkg_dir, "q2_step4_bayesian_shared_outputs.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

saveRDS(q2_shared, file.path(rds_dir, "q2_step4_bayesian_shared_outputs.rds"))

# ----------------------------
# 24. MANUSCRIPT-READY SUMMARY TABLE
# ----------------------------
manuscript_model_summary <- tibble(
  metric = c(
    "Number of catchments",
    "DIC",
    "WAIC",
    "Hotspots (exceedance > 0.80)",
    "Hotspots (exceedance > 0.90)",
    "Robust hotspots, strict LISA p <= 0.05",
    "Robust hotspots, exploratory LISA p <= 0.10"
  ),
  value = c(
    nrow(q2_shared),
    round(model_shared$dic$dic, 2),
    round(model_shared$waic$waic, 2),
    sum(q2_shared$hotspot_080, na.rm = TRUE),
    sum(q2_shared$hotspot_090, na.rm = TRUE),
    sum(q2_shared$robust_hotspot_005, na.rm = TRUE),
    sum(q2_shared$robust_hotspot_010, na.rm = TRUE)
  )
)

write_csv(manuscript_model_summary, file.path(table_dir, "q2_manuscript_model_summary.csv"))

# ----------------------------
# 25. FINAL LOG
# ----------------------------
cat("\nModel summary:\n", file = log_file, append = TRUE)
capture.output(print(manuscript_model_summary), file = log_file, append = TRUE)

cat("\nStep 4 completed: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE)

message("Step 4 completed successfully.")
message("Saved:")
message("- q2_inla_shared_component_model.rds")
message("- q2_shared_spatial_component_and_exceedance.csv")
message("- q2_bayesian_shared_component_outputs.csv")
message("- q2_step4_bayesian_shared_outputs.rds")

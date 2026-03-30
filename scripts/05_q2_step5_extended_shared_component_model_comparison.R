# ============================================================
# OBJECTIVE 2, QUESTION 2
# Step 5: Extended Bayesian shared-component model with
# disease-specific structured residual spatial effects
# and model comparison against Step 4
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
log_file <- file.path(log_dir, "q2_step5_extended_shared_component_model_log.txt")
cat("Step 5 started:", as.character(Sys.time()), "\n", file = log_file)
cat("Root directory:", root_dir, "\n", file = log_file, append = TRUE)

# ----------------------------
# 5. LOAD INPUTS
# ----------------------------
q2_file        <- file.path(rds_dir, "q2_step3_bivariate_outputs.rds")
nb_file        <- file.path(rds_dir, "q2_nb_queen.rds")
step4_model    <- file.path(rds_dir, "q2_inla_shared_component_model.rds")

if (!file.exists(q2_file)) stop("Missing q2_step3_bivariate_outputs.rds. Run Step 3 first.")
if (!file.exists(nb_file)) stop("Missing q2_nb_queen.rds. Run Step 2 first.")

q2 <- readRDS(q2_file)
nb_queen <- readRDS(nb_file)

model_step4 <- NULL
if (file.exists(step4_model)) {
  model_step4 <- readRDS(step4_model)
}

if (!inherits(q2, "sf")) {
  stop("The Step 3 object is not an sf object.")
}

required_cols <- c(
  "catchment_id",
  "area_index",
  "malaria_cases",
  "scd_cases",
  "catchment_population"
)

missing_cols <- setdiff(required_cols, names(q2))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

q2 <- q2 %>%
  arrange(area_index) %>%
  mutate(catchment_population = as.numeric(catchment_population))

if (any(is.na(q2$catchment_population)) || any(q2$catchment_population <= 0)) {
  stop("Some catchments have missing or non-positive population.")
}

cat("Loaded analysis object with ", nrow(q2), " catchments.\n",
    file = log_file, append = TRUE)

# ----------------------------
# 6. CREATE ADJACENCY GRAPH FOR INLA
# ----------------------------
adj_file <- file.path(graph_dir, "q2_queen_adj.graph")

if (!file.exists(adj_file)) {
  spdep::nb2INLA(file = adj_file, nb = nb_queen)
}

if (!file.exists(adj_file)) {
  stop("Failed to create INLA adjacency graph file.")
}

cat("Adjacency graph ready: ", adj_file, "\n", file = log_file, append = TRUE)

# ----------------------------
# 7. PREPARE POPULATION-BASED EXPECTED COUNTS
# ----------------------------
overall_malaria_rate <- sum(q2$malaria_cases, na.rm = TRUE) / sum(q2$catchment_population, na.rm = TRUE)
overall_scd_rate     <- sum(q2$scd_cases, na.rm = TRUE) / sum(q2$catchment_population, na.rm = TRUE)

q2 <- q2 %>%
  mutate(
    expected_malaria_cases = pmax(catchment_population * overall_malaria_rate, 1e-06),
    expected_scd_cases     = pmax(catchment_population * overall_scd_rate, 1e-06)
  )

# ----------------------------
# 8. CREATE STACKED LONG DATASET
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
    outcome_type = factor(outcome_type, levels = c("malaria", "scd")),

    # shared structured component
    shared_main = area_index,
    shared_scd  = ifelse(outcome_type == "scd", area_index, NA_integer_),

    # disease-specific structured spatial residuals
    spatial_malaria = ifelse(outcome_type == "malaria", area_index, NA_integer_),
    spatial_scd     = ifelse(outcome_type == "scd", area_index, NA_integer_),

    # disease-specific iid terms
    iid_malaria = ifelse(outcome_type == "malaria", area_index, NA_integer_),
    iid_scd     = ifelse(outcome_type == "scd", area_index, NA_integer_)
  )

# ----------------------------
# 9. PRIORS
# ----------------------------
hyper_besag <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01))
)

hyper_iid <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01))
)

# ----------------------------
# 10. DEFINE EXTENDED SHARED-COMPONENT MODEL
# ----------------------------
formula_extended <- outcome ~ 0 + outcome_type +
  f(shared_main,
    model = "besag",
    graph = adj_file,
    scale.model = TRUE,
    constr = TRUE,
    hyper = hyper_besag) +
  f(shared_scd,
    copy = "shared_main",
    fixed = FALSE,
    hyper = list(beta = list(prior = "normal", param = c(0, 1)))) +
  f(spatial_malaria,
    model = "besag",
    graph = adj_file,
    scale.model = TRUE,
    constr = TRUE,
    hyper = hyper_besag) +
  f(spatial_scd,
    model = "besag",
    graph = adj_file,
    scale.model = TRUE,
    constr = TRUE,
    hyper = hyper_besag) +
  f(iid_malaria,
    model = "iid",
    hyper = hyper_iid) +
  f(iid_scd,
    model = "iid",
    hyper = hyper_iid)

# ----------------------------
# 11. FIT EXTENDED MODEL
# ----------------------------
cat("Fitting extended shared-component model...\n", file = log_file, append = TRUE)

model_extended <- INLA::inla(
  formula_extended,
  family = "poisson",
  data = dat_long,
  E = expected,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  verbose = FALSE
)

saveRDS(model_extended, file.path(rds_dir, "q2_inla_extended_shared_component_model.rds"))

cat("Extended model fitting completed.\n", file = log_file, append = TRUE)
cat("Extended DIC: ", model_extended$dic$dic, "\n", file = log_file, append = TRUE)
cat("Extended WAIC: ", model_extended$waic$waic, "\n", file = log_file, append = TRUE)

# ----------------------------
# 12. EXTRACT FIXED EFFECTS AND HYPERPARAMETERS
# ----------------------------
fixed_effects_ext <- model_extended$summary.fixed %>%
  tibble::rownames_to_column("term")

write_csv(fixed_effects_ext, file.path(table_dir, "q2_step5_fixed_effects_extended.csv"))

hyperpar_ext <- model_extended$summary.hyperpar %>%
  tibble::rownames_to_column("parameter")

write_csv(hyperpar_ext, file.path(table_dir, "q2_step5_hyperparameters_extended.csv"))

# ----------------------------
# 13. EXTRACT RANDOM EFFECTS
# ----------------------------
extract_random <- function(model, component_name, out_name) {
  model$summary.random[[component_name]] %>%
    mutate(
      area_index = ID,
      mean = mean,
      sd = sd,
      q025 = `0.025quant`,
      q50 = `0.5quant`,
      q975 = `0.975quant`
    ) %>%
    transmute(
      area_index,
      !!paste0(out_name, "_mean") := mean,
      !!paste0(out_name, "_sd")   := sd,
      !!paste0(out_name, "_q025") := q025,
      !!paste0(out_name, "_q50")  := q50,
      !!paste0(out_name, "_q975") := q975
    )
}

shared_ext   <- extract_random(model_extended, "shared_main", "shared")
malaria_ext  <- extract_random(model_extended, "spatial_malaria", "malaria_spatial")
scd_ext      <- extract_random(model_extended, "spatial_scd", "scd_spatial")

# ----------------------------
# 14. EXCEEDANCE PROBABILITY FOR SHARED COMPONENT
# ----------------------------
shared_ext <- shared_ext %>%
  mutate(
    exceed_prob_rr_gt1 = 1 - pnorm(0, mean = shared_mean, sd = shared_sd),
    hotspot_080 = ifelse(exceed_prob_rr_gt1 > 0.80, 1L, 0L),
    hotspot_090 = ifelse(exceed_prob_rr_gt1 > 0.90, 1L, 0L)
  )

write_csv(shared_ext, file.path(table_dir, "q2_step5_shared_component_exceedance.csv"))

# ----------------------------
# 15. EXTRACT FITTED VALUES AND RELATIVE RISKS
# ----------------------------
n_catch <- nrow(q2)

fitted_vals <- model_extended$summary.fitted.values

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

# ----------------------------
# 16. MERGE MODEL OUTPUTS BACK TO SPATIAL DATA
# ----------------------------
q2_ext <- q2 %>%
  left_join(shared_ext, by = "area_index") %>%
  left_join(malaria_ext, by = "area_index") %>%
  left_join(scd_ext, by = "area_index") %>%
  left_join(fitted_mal, by = "area_index") %>%
  left_join(fitted_scd, by = "area_index") %>%
  mutate(
    rr_malaria = fitted_malaria_mean / expected_malaria_cases,
    rr_scd     = fitted_scd_mean / expected_scd_cases
  )

# ----------------------------
# 17. HOTSPOT CONCORDANCE WITH STEP 3 LISA
# ----------------------------
if (!"cases_cluster_005" %in% names(q2_ext)) {
  q2_ext$cases_cluster_005 <- NA_character_
}
if (!"cases_cluster_010" %in% names(q2_ext)) {
  q2_ext$cases_cluster_010 <- NA_character_
}

q2_ext <- q2_ext %>%
  mutate(
    robust_hotspot_005 = ifelse(cases_cluster_005 == "High-High" & exceed_prob_rr_gt1 > 0.80, 1L, 0L),
    robust_hotspot_010 = ifelse(cases_cluster_010 == "High-High" & exceed_prob_rr_gt1 > 0.80, 1L, 0L),

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

hotspot_summary_ext <- q2_ext %>%
  st_drop_geometry() %>%
  summarise(
    n_catchments = n(),
    n_hotspot_080 = sum(hotspot_080, na.rm = TRUE),
    n_hotspot_090 = sum(hotspot_090, na.rm = TRUE),
    n_robust_hotspot_005 = sum(robust_hotspot_005, na.rm = TRUE),
    n_robust_hotspot_010 = sum(robust_hotspot_010, na.rm = TRUE)
  )

write_csv(hotspot_summary_ext, file.path(table_dir, "q2_step5_hotspot_summary_extended.csv"))

# ----------------------------
# 18. MODEL COMPARISON TABLE
# ----------------------------
model_comparison <- tibble(
  model = c("Step 4 baseline shared model", "Step 5 extended shared model"),
  dic = c(
    if (!is.null(model_step4)) model_step4$dic$dic else NA_real_,
    model_extended$dic$dic
  ),
  waic = c(
    if (!is.null(model_step4)) model_step4$waic$waic else NA_real_,
    model_extended$waic$waic
  ),
  n_hotspot_080 = c(
    if (!is.null(model_step4) && file.exists(file.path(table_dir, "q2_hotspot_summary.csv"))) {
      read_csv(file.path(table_dir, "q2_hotspot_summary.csv"), show_col_types = FALSE)$n_hotspot_080[1]
    } else NA_real_,
    hotspot_summary_ext$n_hotspot_080[1]
  ),
  n_hotspot_090 = c(
    if (!is.null(model_step4) && file.exists(file.path(table_dir, "q2_hotspot_summary.csv"))) {
      read_csv(file.path(table_dir, "q2_hotspot_summary.csv"), show_col_types = FALSE)$n_hotspot_090[1]
    } else NA_real_,
    hotspot_summary_ext$n_hotspot_090[1]
  ),
  n_robust_hotspot_005 = c(
    if (!is.null(model_step4) && file.exists(file.path(table_dir, "q2_hotspot_summary.csv"))) {
      read_csv(file.path(table_dir, "q2_hotspot_summary.csv"), show_col_types = FALSE)$n_robust_hotspot_005[1]
    } else NA_real_,
    hotspot_summary_ext$n_robust_hotspot_005[1]
  ),
  n_robust_hotspot_010 = c(
    if (!is.null(model_step4) && file.exists(file.path(table_dir, "q2_hotspot_summary.csv"))) {
      read_csv(file.path(table_dir, "q2_hotspot_summary.csv"), show_col_types = FALSE)$n_robust_hotspot_010[1]
    } else NA_real_,
    hotspot_summary_ext$n_robust_hotspot_010[1]
  )
)

write_csv(model_comparison, file.path(table_dir, "q2_step5_model_comparison.csv"))

# ----------------------------
# 19. MAPS
# ----------------------------
tmap_mode("plot")

map_shared <- tm_shape(q2_ext) +
  tm_polygons(
    "shared_mean",
    title = "Shared spatial effect",
    palette = "-RdBu",
    style = "cont"
  ) +
  tm_layout(
    main.title = "Extended model: posterior mean shared spatial effect",
    legend.outside = TRUE
  )

map_malaria_resid <- tm_shape(q2_ext) +
  tm_polygons(
    "malaria_spatial_mean",
    title = "Malaria-specific residual",
    palette = "-RdBu",
    style = "cont"
  ) +
  tm_layout(
    main.title = "Extended model: malaria-specific residual spatial effect",
    legend.outside = TRUE
  )

map_scd_resid <- tm_shape(q2_ext) +
  tm_polygons(
    "scd_spatial_mean",
    title = "SCD-specific residual",
    palette = "-RdBu",
    style = "cont"
  ) +
  tm_layout(
    main.title = "Extended model: SCD-specific residual spatial effect",
    legend.outside = TRUE
  )

map_exceed <- tm_shape(q2_ext) +
  tm_polygons(
    "exceed_prob_rr_gt1",
    title = "P(shared RR > 1)",
    palette = "YlOrRd",
    style = "cont"
  ) +
  tm_layout(
    main.title = "Extended model: posterior exceedance probability",
    legend.outside = TRUE
  )

hotspot_palette <- c(
  "Robust hotspot" = "#b2182b",
  "Bayesian hotspot only" = "#ef8a62",
  "LISA hotspot only" = "#67a9cf",
  "Not hotspot" = "grey85"
)

map_hotspot_005 <- tm_shape(q2_ext) +
  tm_polygons(
    "hotspot_class_005",
    title = "Hotspot class (p ≤ 0.05)",
    palette = hotspot_palette
  ) +
  tm_layout(
    main.title = "Extended model hotspot concordance, strict threshold",
    legend.outside = TRUE
  )

map_hotspot_010 <- tm_shape(q2_ext) +
  tm_polygons(
    "hotspot_class_010",
    title = "Hotspot class (p ≤ 0.10)",
    palette = hotspot_palette
  ) +
  tm_layout(
    main.title = "Extended model hotspot concordance, exploratory threshold",
    legend.outside = TRUE
  )

tmap_save(map_shared,      file.path(fig_dir, "q2_step5_shared_effect_map.png"),           width = 9, height = 6, dpi = 300)
tmap_save(map_malaria_resid, file.path(fig_dir, "q2_step5_malaria_residual_map.png"),      width = 9, height = 6, dpi = 300)
tmap_save(map_scd_resid,   file.path(fig_dir, "q2_step5_scd_residual_map.png"),            width = 9, height = 6, dpi = 300)
tmap_save(map_exceed,      file.path(fig_dir, "q2_step5_exceedance_map.png"),              width = 9, height = 6, dpi = 300)
tmap_save(map_hotspot_005, file.path(fig_dir, "q2_step5_hotspot_concordance_p005.png"),    width = 9, height = 6, dpi = 300)
tmap_save(map_hotspot_010, file.path(fig_dir, "q2_step5_hotspot_concordance_p010.png"),    width = 9, height = 6, dpi = 300)

# ----------------------------
# 20. EXPORT MAIN TABLE
# ----------------------------
main_export <- q2_ext %>%
  st_drop_geometry() %>%
  select(
    catchment_id,
    area_index,
    catchment_population,
    malaria_cases,
    scd_cases,
    expected_malaria_cases,
    expected_scd_cases,
    shared_mean,
    shared_sd,
    malaria_spatial_mean,
    malaria_spatial_sd,
    scd_spatial_mean,
    scd_spatial_sd,
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

write_csv(main_export, file.path(table_dir, "q2_step5_extended_model_outputs.csv"))

# ----------------------------
# 21. SAVE SPATIAL OUTPUTS
# ----------------------------
sf::st_write(
  q2_ext,
  file.path(gpkg_dir, "q2_step5_extended_model_outputs.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

saveRDS(q2_ext, file.path(rds_dir, "q2_step5_extended_model_outputs.rds"))

# ----------------------------
# 22. FINAL LOG
# ----------------------------
cat("\nModel comparison:\n", file = log_file, append = TRUE)
capture.output(print(model_comparison), file = log_file, append = TRUE)

cat("\nExtended hotspot summary:\n", file = log_file, append = TRUE)
capture.output(print(hotspot_summary_ext), file = log_file, append = TRUE)

cat("\nStep 5 completed: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE)

message("Step 5 completed successfully.")

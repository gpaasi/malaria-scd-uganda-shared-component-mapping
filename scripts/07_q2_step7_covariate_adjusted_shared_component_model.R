# ============================================================
# OBJECTIVE 2, QUESTION 2
# Step 7: Covariate-adjusted shared-component model
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
  "ggplot2",
  "tmap",
  "janitor",
  "tidyr",
  "purrr",
  "tibble",
  "stringr",
  "INLA"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_pkgs[required_pkgs != "INLA"], installed)

if (length(to_install) > 0) {
  install.packages(to_install)
}

if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages(
    "INLA",
    repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable")
  )
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
log_file <- file.path(log_dir, "q2_step7_covariate_adjusted_model_log.txt")
cat("Step 7 started:", as.character(Sys.time()), "\n", file = log_file)
cat("Root directory:", root_dir, "\n", file = log_file, append = TRUE)

# ----------------------------
# 5. LOAD INPUTS
# ----------------------------
step5_file <- file.path(rds_dir, "q2_step5_extended_shared_component_outputs.rds")
step1_file <- file.path(rds_dir, "q2_data_step1_cleaned.rds")

if (!file.exists(step5_file)) stop("Missing q2_step5_extended_shared_component_outputs.rds. Run Step 5 first.")
if (!file.exists(step1_file)) stop("Missing q2_data_step1_cleaned.rds. Run Step 1 first.")

step5 <- readRDS(step5_file)
step1 <- readRDS(step1_file)

if (!inherits(step5, "sf")) stop("Step 5 object is not sf.")
if (!inherits(step1, "sf")) stop("Step 1 object is not sf.")

cat("Loaded Step 5 and Step 1 outputs.\n", file = log_file, append = TRUE)

# ----------------------------
# 6. CHECK REQUIRED COLUMNS
# ----------------------------
needed_step5 <- c(
  "catchment_id",
  "area_index",
  "facility_name",
  "malaria_cases",
  "scd_cases",
  "catchment_population",
  "expected_malaria",
  "expected_scd",
  "adjacency_graph"
)

missing_step5 <- setdiff(needed_step5, names(step5))
if (length(missing_step5) > 0) {
  stop("Missing Step 5 columns: ", paste(missing_step5, collapse = ", "))
}

# Covariates expected from the linked facility/catchment file.
# Adjust these names below if your Step 1 file uses slightly different labels.
candidate_covars <- c(
  "type", "facility_type", "facility_class",
  "ownership", "ownership_class",
  "f15regions", "region", "tpopn"
)

available_covars <- intersect(candidate_covars, names(step1))
if (length(available_covars) == 0) {
  stop("No expected Step 7 covariate columns found in Step 1 object.")
}

# ----------------------------
# 7. BUILD ANALYSIS DATASET
# ----------------------------
base_tab <- step5 %>%
  st_drop_geometry() %>%
  select(
    catchment_id,
    area_index,
    facility_name,
    malaria_cases,
    scd_cases,
    catchment_population,
    expected_malaria,
    expected_scd,
    adjacency_graph
  )

covar_tab <- step1 %>%
  st_drop_geometry() %>%
  select(any_of(c("catchment_id", "area_index", available_covars)))

q2_adj <- step5 %>%
  select(catchment_id, area_index, geometry) %>%
  left_join(base_tab, by = c("catchment_id", "area_index")) %>%
  left_join(covar_tab, by = c("catchment_id", "area_index"))

# Harmonise covariate names
q2_adj <- q2_adj %>%
  mutate(
    tpopn = dplyr::coalesce(.data[[intersect(c("tpopn", "catchment_population"), names(.))[[1]]]], catchment_population)
  )

if ("facility_class" %in% names(q2_adj)) {
  q2_adj <- q2_adj %>% mutate(facility_class = as.character(facility_class))
} else if ("type" %in% names(q2_adj)) {
  q2_adj <- q2_adj %>%
    mutate(facility_class = case_when(
      stringr::str_detect(stringr::str_to_lower(type), "hc iv|health centre iv") ~ "HC IV",
      stringr::str_detect(stringr::str_to_lower(type), "general hospital") ~ "General hospital",
      stringr::str_detect(stringr::str_to_lower(type), "hospital") ~ "Hospital",
      TRUE ~ as.character(type)
    ))
} else if ("facility_type" %in% names(q2_adj)) {
  q2_adj <- q2_adj %>% mutate(facility_class = as.character(facility_type))
}

if ("ownership_class" %in% names(q2_adj)) {
  q2_adj <- q2_adj %>% mutate(ownership_class = as.character(ownership_class))
} else if ("ownership" %in% names(q2_adj)) {
  q2_adj <- q2_adj %>%
    mutate(ownership_class = case_when(
      stringr::str_detect(stringr::str_to_lower(ownership), "government|public") ~ "Government",
      stringr::str_detect(stringr::str_to_lower(ownership), "private not for profit|pnfp|mission") ~ "PNFP",
      stringr::str_detect(stringr::str_to_lower(ownership), "private for profit|pfp") ~ "PFP",
      TRUE ~ as.character(ownership)
    ))
}

if ("region" %in% names(q2_adj)) {
  q2_adj <- q2_adj %>% mutate(region = as.character(region))
} else if ("f15regions" %in% names(q2_adj)) {
  q2_adj <- q2_adj %>% mutate(region = as.character(f15regions))
}

# Minimal cleaning
q2_adj <- q2_adj %>%
  mutate(
    facility_class = forcats::fct_drop(as.factor(facility_class)),
    ownership_class = forcats::fct_drop(as.factor(ownership_class)),
    region = forcats::fct_drop(as.factor(region)),
    log_tpopn = log1p(tpopn)
  )

# Set reference categories used in manuscript tables
if ("HC IV" %in% levels(q2_adj$facility_class)) {
  q2_adj$facility_class <- stats::relevel(q2_adj$facility_class, ref = "HC IV")
}
if ("Government" %in% levels(q2_adj$ownership_class)) {
  q2_adj$ownership_class <- stats::relevel(q2_adj$ownership_class, ref = "Government")
}
if ("Acholi" %in% levels(q2_adj$region)) {
  q2_adj$region <- stats::relevel(q2_adj$region, ref = "Acholi")
}

# ----------------------------
# 8. COVARIATE SUMMARIES
# ----------------------------
summary_numeric <- q2_adj %>%
  st_drop_geometry() %>%
  summarise(
    mean_tpopn = mean(tpopn, na.rm = TRUE),
    median_tpopn = median(tpopn, na.rm = TRUE),
    min_tpopn = min(tpopn, na.rm = TRUE),
    max_tpopn = max(tpopn, na.rm = TRUE)
  )

summary_facility <- q2_adj %>%
  st_drop_geometry() %>%
  count(facility_class, name = "n")

summary_ownership <- q2_adj %>%
  st_drop_geometry() %>%
  count(ownership_class, name = "n")

summary_region <- q2_adj %>%
  st_drop_geometry() %>%
  count(region, name = "n")

write_csv(summary_numeric, file.path(table_dir, "q2_step7_covariate_summary_numeric.csv"))
write_csv(summary_facility, file.path(table_dir, "q2_step7_covariate_summary_facility_class.csv"))
write_csv(summary_ownership, file.path(table_dir, "q2_step7_covariate_summary_ownership.csv"))
write_csv(summary_region, file.path(table_dir, "q2_step7_covariate_summary_region.csv"))

# ----------------------------
# 9. MODEL DATA IN STACKED FORM
# ----------------------------
# Two-outcome shared-component layout for INLA
n <- nrow(q2_adj)

stack_df <- bind_rows(
  q2_adj %>%
    st_drop_geometry() %>%
    transmute(
      catchment_id, area_index,
      disease = "malaria",
      y = malaria_cases,
      E = expected_malaria,
      log_tpopn,
      facility_class,
      ownership_class,
      region
    ),
  q2_adj %>%
    st_drop_geometry() %>%
    transmute(
      catchment_id, area_index,
      disease = "scd",
      y = scd_cases,
      E = expected_scd,
      log_tpopn,
      facility_class,
      ownership_class,
      region
    )
) %>%
  mutate(
    disease = factor(disease, levels = c("malaria", "scd")),
    disease_index = if_else(disease == "malaria", 1L, 2L),
    shared_index = area_index,
    spatial_malaria = if_else(disease == "malaria", area_index, NA_integer_),
    spatial_scd = if_else(disease == "scd", area_index, NA_integer_),
    iid_malaria = if_else(disease == "malaria", area_index, NA_integer_),
    iid_scd = if_else(disease == "scd", area_index, NA_integer_),
    scd_copy_index = if_else(disease == "scd", area_index, NA_integer_)
  )

# ----------------------------
# 10. FIT NESTED ADJUSTED MODELS
# ----------------------------
# NOTE:
# These formulations assume the Step 5 object already includes the adjacency graph
# in INLA-readable format as 'adjacency_graph'.

graph_file <- file.path(rds_dir, "q2_adjacency_inla.graph")
if (!file.exists(graph_file)) {
  # try to recover graph file from step5 object if present as path string
  if (is.character(step5$adjacency_graph[[1]]) && file.exists(step5$adjacency_graph[[1]])) {
    graph_file <- step5$adjacency_graph[[1]]
  } else {
    stop("Missing INLA graph file. Save the adjacency graph from earlier steps as q2_adjacency_inla.graph")
  }
}

base_formula <- y ~ -1 + disease +
  f(shared_index, model = "besag", graph = graph_file) +
  f(iid_malaria, model = "iid") +
  f(iid_scd, model = "iid") +
  f(spatial_malaria, model = "besag", graph = graph_file) +
  f(spatial_scd, model = "besag", graph = graph_file) +
  f(scd_copy_index, copy = "shared_index", fixed = FALSE)

formula_pop <- update(base_formula, . ~ . + log_tpopn)
formula_pop_fac <- update(formula_pop, . ~ . + facility_class)
formula_full <- update(formula_pop_fac, . ~ . + ownership_class + region)

control_compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE)

model_pop <- inla(
  formula_pop,
  family = "poisson",
  data = stack_df,
  E = E,
  control.compute = control_compute,
  verbose = FALSE
)

model_pop_fac <- inla(
  formula_pop_fac,
  family = "poisson",
  data = stack_df,
  E = E,
  control.compute = control_compute,
  verbose = FALSE
)

model_full <- inla(
  formula_full,
  family = "poisson",
  data = stack_df,
  E = E,
  control.compute = control_compute,
  verbose = FALSE
)

# ----------------------------
# 11. MODEL COMPARISON TABLE
# ----------------------------
model_comp <- tibble(
  model = c("Population-adjusted", "Population + facility class", "Fully adjusted"),
  dic = c(model_pop$dic$dic, model_pop_fac$dic$dic, model_full$dic$dic),
  waic = c(model_pop$waic$waic, model_pop_fac$waic$waic, model_full$waic$waic)
)
write_csv(model_comp, file.path(table_dir, "q2_step7_model_comparison.csv"))

# ----------------------------
# 12. EXTRACT FIXED EFFECTS AND HYPERPARAMETERS
# ----------------------------
extract_fixed <- function(fit) {
  fit$summary.fixed %>%
    as.data.frame() %>%
    tibble::rownames_to_column("parameter") %>%
    janitor::clean_names()
}

extract_hyper <- function(fit) {
  fit$summary.hyperpar %>%
    as.data.frame() %>%
    tibble::rownames_to_column("parameter") %>%
    janitor::clean_names()
}

write_csv(extract_fixed(model_full), file.path(table_dir, "q2_step7_fixed_effects_final_adjusted.csv"))
write_csv(extract_hyper(model_full), file.path(table_dir, "q2_step7_hyperparameters_final_adjusted.csv"))

# ----------------------------
# 13. DERIVE ADJUSTED SHARED SURFACE AND EXCEEDANCE
# ----------------------------
shared_summary <- model_full$summary.random$shared_index %>%
  as.data.frame() %>%
  transmute(
    area_index = ID,
    shared_mean_adj = mean,
    shared_sd_adj = sd,
    shared_q025_adj = `0.025quant`,
    shared_q975_adj = `0.975quant`
  )

q2_adj <- q2_adj %>%
  left_join(shared_summary, by = "area_index") %>%
  mutate(
    rr_shared_adj = exp(shared_mean_adj),
    exceed_prob_rr_gt1_adj = 1 - pnorm(0, mean = shared_mean_adj, sd = shared_sd_adj),
    hotspot_080_adj = as.integer(exceed_prob_rr_gt1_adj > 0.80),
    hotspot_090_adj = as.integer(exceed_prob_rr_gt1_adj > 0.90)
  )

# approximate disease-specific adjusted relative risks from fitted values
fit_sum <- model_full$summary.fitted.values %>% as.data.frame()
if (nrow(fit_sum) == nrow(stack_df)) {
  stack_df$rr_fit <- fit_sum$mean
  rr_mal <- stack_df %>% filter(disease == "malaria") %>% transmute(area_index, rr_malaria_adj = rr_fit)
  rr_scd <- stack_df %>% filter(disease == "scd") %>% transmute(area_index, rr_scd_adj = rr_fit)
  q2_adj <- q2_adj %>% left_join(rr_mal, by = "area_index") %>% left_join(rr_scd, by = "area_index")
}

# ----------------------------
# 14. TOP ADJUSTED HOTSPOTS
# ----------------------------
top_adjusted_hotspots <- q2_adj %>%
  st_drop_geometry() %>%
  filter(hotspot_080_adj == 1) %>%
  arrange(desc(exceed_prob_rr_gt1_adj), desc(shared_mean_adj)) %>%
  select(
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
    exceed_prob_rr_gt1_adj,
    hotspot_080_adj,
    hotspot_090_adj,
    rr_malaria_adj,
    rr_scd_adj
  )

hotspot_summary_adjusted <- q2_adj %>%
  st_drop_geometry() %>%
  summarise(
    n_hotspot_080_adj = sum(hotspot_080_adj == 1, na.rm = TRUE),
    n_hotspot_090_adj = sum(hotspot_090_adj == 1, na.rm = TRUE)
  )

shared_component_exceedance_adjusted <- q2_adj %>%
  st_drop_geometry() %>%
  select(catchment_id, area_index, facility_name, shared_mean_adj, shared_sd_adj,
         shared_q025_adj, shared_q975_adj, exceed_prob_rr_gt1_adj,
         hotspot_080_adj, hotspot_090_adj)

write_csv(top_adjusted_hotspots, file.path(table_dir, "q2_step7_top_adjusted_hotspots.csv"))
write_csv(hotspot_summary_adjusted, file.path(table_dir, "q2_step7_hotspot_summary_adjusted.csv"))
write_csv(shared_component_exceedance_adjusted, file.path(table_dir, "q2_step7_shared_component_exceedance_adjusted.csv"))

# ----------------------------
# 15. MAPS
# ----------------------------
tmap_mode("plot")

map_shared_adj <- tm_shape(q2_adj) +
  tm_polygons(
    "shared_mean_adj",
    title = "Adjusted shared effect",
    palette = "RdYlBu",
    style = "cont",
    midpoint = 0
  ) +
  tm_layout(main.title = "Adjusted shared malaria-SCD surface", legend.outside = TRUE)

map_exceed_adj <- tm_shape(q2_adj) +
  tm_polygons(
    "exceed_prob_rr_gt1_adj",
    title = "Exceedance probability",
    palette = "YlOrRd",
    style = "cont"
  ) +
  tm_layout(main.title = "Adjusted exceedance probability surface", legend.outside = TRUE)

map_hotspots_adj <- tm_shape(q2_adj) +
  tm_polygons(
    "hotspot_080_adj",
    title = "Hotspot >0.80",
    palette = c("grey85", "#b2182b"),
    labels = c("No", "Yes")
  ) +
  tm_layout(main.title = "Adjusted Bayesian hotspots", legend.outside = TRUE)

tmap_save(map_shared_adj, file.path(fig_dir, "q2_step7_adjusted_shared_surface.png"), width = 8, height = 6, dpi = 300)
tmap_save(map_exceed_adj, file.path(fig_dir, "q2_step7_adjusted_exceedance_surface.png"), width = 8, height = 6, dpi = 300)
tmap_save(map_hotspots_adj, file.path(fig_dir, "q2_step7_adjusted_hotspot_map.png"), width = 8, height = 6, dpi = 300)

# ----------------------------
# 16. SAVE OUTPUTS
# ----------------------------
sf::st_write(
  q2_adj,
  dsn = file.path(gpkg_dir, "q2_step7_adjusted_model_outputs.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

saveRDS(q2_adj, file.path(rds_dir, "q2_step7_adjusted_model_outputs.rds"))

cat("Saved Step 7 outputs.\n", file = log_file, append = TRUE)
cat("Step 7 completed:", as.character(Sys.time()), "\n", file = log_file, append = TRUE)
message("Step 7 completed successfully.")

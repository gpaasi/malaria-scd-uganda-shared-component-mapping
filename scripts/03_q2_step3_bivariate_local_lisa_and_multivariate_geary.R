# ============================================================
# OBJECTIVE 2, QUESTION 2
# Step 3: Bivariate spatial autocorrelation and local co-clustering
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
step2_file <- file.path(rds_dir, "q2_data_step2_with_weights_and_lags.rds")
nb_file    <- file.path(rds_dir, "q2_nb_queen.rds")
lw_file    <- file.path(rds_dir, "q2_listw_queen_W.rds")

if (!file.exists(step2_file)) stop("Missing Step 2 data file.")
if (!file.exists(nb_file)) stop("Missing queen neighbour object from Step 2.")
if (!file.exists(lw_file)) stop("Missing queen weights object from Step 2.")

# ----------------------------
# 5. START LOG
# ----------------------------
log_file <- file.path(log_dir, "q2_step3_bivariate_esda_log.txt")
cat("Q2 Step 3 started:", as.character(Sys.time()), "\n", file = log_file)
cat("Root directory:", root_dir, "\n", file = log_file, append = TRUE)

# ----------------------------
# 6. HELPER FUNCTIONS
# ----------------------------
run_global_bivariate_moran <- function(x, y, listw, pairing, nsim = 999, seed = 1234) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  lag_y <- spdep::lag.listw(listw, y, zero.policy = TRUE)
  obs <- suppressWarnings(cor(x, lag_y, use = "complete.obs"))

  set.seed(seed)
  sim <- replicate(nsim, {
    y_perm <- sample(y, replace = FALSE)
    lag_y_perm <- spdep::lag.listw(listw, y_perm, zero.policy = TRUE)
    suppressWarnings(cor(x, lag_y_perm, use = "complete.obs"))
  })

  tibble(
    pairing = pairing,
    observed_ib = obs,
    mean_perm = mean(sim, na.rm = TRUE),
    sd_perm = sd(sim, na.rm = TRUE),
    z_score = ifelse(sd(sim, na.rm = TRUE) > 0, (obs - mean(sim, na.rm = TRUE)) / sd(sim, na.rm = TRUE), NA_real_),
    p_value_upper = (sum(sim >= obs, na.rm = TRUE) + 1) / (nsim + 1),
    p_value_lower = (sum(sim <= obs, na.rm = TRUE) + 1) / (nsim + 1),
    p_value_two_sided = (sum(abs(sim - mean(sim, na.rm = TRUE)) >= abs(obs - mean(sim, na.rm = TRUE)), na.rm = TRUE) + 1) / (nsim + 1),
    nsim = nsim
  )
}

multivariate_geary_global <- function(x, y, nb_obj, style = "W", nsim = 999, seed = 1234, pairing = NULL) {
  X <- cbind(as.numeric(x), as.numeric(y))
  lw <- spdep::nb2listw(nb_obj, style = style, zero.policy = TRUE)

  edge_list <- do.call(
    rbind,
    lapply(seq_along(nb_obj), function(i) {
      ni <- nb_obj[[i]]
      wi <- lw$weights[[i]]

      if (length(ni) == 0 || length(wi) == 0) {
        return(NULL)
      }

      data.frame(
        i = rep(i, length(ni)),
        j = ni,
        w = wi
      )
    })
  )

  if (is.null(edge_list) || nrow(edge_list) == 0) {
    stop("No neighbour links found for Geary's C computation.")
  }

  calc_c <- function(mat) {
    diffsq <- (mat[edge_list$i, 1] - mat[edge_list$j, 1])^2 +
      (mat[edge_list$i, 2] - mat[edge_list$j, 2])^2

    num <- sum(edge_list$w * diffsq, na.rm = TRUE)
    Wsum <- sum(edge_list$w, na.rm = TRUE)
    p <- ncol(mat)

    num / (2 * Wsum * p)
  }

  obs <- calc_c(X)

  set.seed(seed)
  sim <- replicate(nsim, {
    idx <- sample(seq_len(nrow(X)), replace = FALSE)
    Xperm <- X[idx, , drop = FALSE]
    calc_c(Xperm)
  })

  tibble(
    pairing = pairing,
    observed_c = obs,
    mean_perm = mean(sim, na.rm = TRUE),
    sd_perm = sd(sim, na.rm = TRUE),
    z_score = ifelse(sd(sim, na.rm = TRUE) > 0, (obs - mean(sim, na.rm = TRUE)) / sd(sim, na.rm = TRUE), NA_real_),
    p_value_lower = (sum(sim <= obs, na.rm = TRUE) + 1) / (nsim + 1),
    p_value_two_sided = (sum(abs(sim - mean(sim, na.rm = TRUE)) >= abs(obs - mean(sim, na.rm = TRUE)), na.rm = TRUE) + 1) / (nsim + 1),
    nsim = nsim
  )
}

run_local_bivariate_lisa <- function(x, y, listw, nsim = 999, seed = 1234) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  lag_y <- spdep::lag.listw(listw, y, zero.policy = TRUE)

  obs <- x * lag_y

  set.seed(seed)
  perm_mat <- replicate(nsim, {
    y_perm <- sample(y, replace = FALSE)
    lag_perm <- spdep::lag.listw(listw, y_perm, zero.policy = TRUE)
    x * lag_perm
  })

  p_two <- vapply(seq_along(obs), function(i) {
    (sum(abs(perm_mat[i, ] - mean(perm_mat[i, ], na.rm = TRUE)) >= abs(obs[i] - mean(perm_mat[i, ], na.rm = TRUE)), na.rm = TRUE) + 1) / (nsim + 1)
  }, numeric(1))

  quadrant <- case_when(
    x >= 0 & lag_y >= 0 ~ "High-High",
    x <  0 & lag_y <  0 ~ "Low-Low",
    x >= 0 & lag_y <  0 ~ "High-Low",
    x <  0 & lag_y >= 0 ~ "Low-High",
    TRUE ~ "Undefined"
  )

  tibble(
    local_bivariate_i = obs,
    local_bivariate_p = p_two,
    lisa_quadrant = quadrant,
    lisa_sig_005 = p_two <= 0.05,
    lisa_sig_010 = p_two <= 0.10,
    hotspot_005 = p_two <= 0.05 & quadrant == "High-High",
    hotspot_010 = p_two <= 0.10 & quadrant == "High-High"
  )
}

make_bivariate_scatter <- function(data, x, y, title_text, xlab, ylab, outfile) {
  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = title_text, x = xlab, y = ylab) +
    theme_minimal()

  ggsave(
    filename = file.path(fig_dir, outfile),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# ----------------------------
# 7. LOAD STEP 2 DATA
# ----------------------------
q2 <- readRDS(step2_file)
nb_queen <- readRDS(nb_file)
lw_queen <- readRDS(lw_file)

if (!inherits(q2, "sf")) stop("Step 2 object is not an sf object.")

required_cols <- c(
  "catchment_id", "area_index",
  "z_malaria_cases", "z_scd_cases",
  "z_malaria_deaths", "z_scd_deaths"
)

missing_cols <- setdiff(required_cols, names(q2))
if (length(missing_cols) > 0) {
  stop("Missing required columns in Step 2 data: ", paste(missing_cols, collapse = ", "))
}

q2 <- q2 %>% arrange(area_index)

# ----------------------------
# 8. DIAGNOSTIC FOR ISLANDS
# ----------------------------
n_islands <- sum(spdep::card(nb_queen) == 0)
cat("Number of island catchments: ", n_islands, "\n", file = log_file, append = TRUE)

# ----------------------------
# 9. GLOBAL BIVARIATE MORAN'S I
# ----------------------------
global_biv_results <- bind_rows(
  run_global_bivariate_moran(q2$z_malaria_cases,  q2$z_scd_cases,  lw_queen, "Malaria cases vs neighbouring SCD cases"),
  run_global_bivariate_moran(q2$z_scd_cases,      q2$z_malaria_cases, lw_queen, "SCD cases vs neighbouring malaria cases"),
  run_global_bivariate_moran(q2$z_malaria_deaths, q2$z_scd_deaths, lw_queen, "Malaria deaths vs neighbouring SCD deaths"),
  run_global_bivariate_moran(q2$z_scd_deaths,     q2$z_malaria_deaths, lw_queen, "SCD deaths vs neighbouring malaria deaths")
)

write_csv(global_biv_results, file.path(table_dir, "q2_global_bivariate_morans_i.csv"))

# ----------------------------
# 10. LOCAL BIVARIATE LISA
# ----------------------------
cases_lisa <- run_local_bivariate_lisa(q2$z_malaria_cases, q2$z_scd_cases, lw_queen)
deaths_lisa <- run_local_bivariate_lisa(q2$z_malaria_deaths, q2$z_scd_deaths, lw_queen)

q2 <- q2 %>%
  bind_cols(
    cases_lisa %>%
      rename(
        cases_local_bivariate_i = local_bivariate_i,
        cases_local_bivariate_p = local_bivariate_p,
        cases_lisa_quadrant = lisa_quadrant,
        cases_cluster_sig = lisa_sig_005,
        cases_cluster_sig_010 = lisa_sig_010,
        cases_hotspot_005 = hotspot_005,
        cases_hotspot_010 = hotspot_010
      ),
    deaths_lisa %>%
      rename(
        deaths_local_bivariate_i = local_bivariate_i,
        deaths_local_bivariate_p = local_bivariate_p,
        deaths_lisa_quadrant = lisa_quadrant,
        deaths_cluster_sig = lisa_sig_005,
        deaths_cluster_sig_010 = lisa_sig_010,
        deaths_hotspot_005 = hotspot_005,
        deaths_hotspot_010 = hotspot_010
      )
  ) %>%
  mutate(
    cases_cluster_005 = case_when(
      cases_hotspot_005 ~ "High-High",
      cases_cluster_sig & cases_lisa_quadrant == "Low-Low" ~ "Low-Low",
      cases_cluster_sig & cases_lisa_quadrant == "High-Low" ~ "High-Low",
      cases_cluster_sig & cases_lisa_quadrant == "Low-High" ~ "Low-High",
      TRUE ~ "Not significant"
    ),
    cases_cluster_010 = case_when(
      cases_hotspot_010 ~ "High-High",
      cases_cluster_sig_010 & cases_lisa_quadrant == "Low-Low" ~ "Low-Low",
      cases_cluster_sig_010 & cases_lisa_quadrant == "High-Low" ~ "High-Low",
      cases_cluster_sig_010 & cases_lisa_quadrant == "Low-High" ~ "Low-High",
      TRUE ~ "Not significant"
    ),
    deaths_cluster_005 = case_when(
      deaths_hotspot_005 ~ "High-High",
      deaths_cluster_sig & deaths_lisa_quadrant == "Low-Low" ~ "Low-Low",
      deaths_cluster_sig & deaths_lisa_quadrant == "High-Low" ~ "High-Low",
      deaths_cluster_sig & deaths_lisa_quadrant == "Low-High" ~ "Low-High",
      TRUE ~ "Not significant"
    ),
    deaths_cluster_010 = case_when(
      deaths_hotspot_010 ~ "High-High",
      deaths_cluster_sig_010 & deaths_lisa_quadrant == "Low-Low" ~ "Low-Low",
      deaths_cluster_sig_010 & deaths_lisa_quadrant == "High-Low" ~ "High-Low",
      deaths_cluster_sig_010 & deaths_lisa_quadrant == "Low-High" ~ "Low-High",
      TRUE ~ "Not significant"
    )
  )

local_lisa_outputs <- q2 %>%
  st_drop_geometry() %>%
  select(
    catchment_id, area_index,
    starts_with("cases_"),
    starts_with("deaths_")
  )

write_csv(local_lisa_outputs, file.path(table_dir, "q2_local_bivariate_lisa_outputs.csv"))

# ----------------------------
# 11. GLOBAL MULTIVARIATE GEARY'S C
# ----------------------------
geary_results <- bind_rows(
  multivariate_geary_global(q2$z_malaria_cases, q2$z_scd_cases, nb_queen, pairing = "Cases"),
  multivariate_geary_global(q2$z_malaria_deaths, q2$z_scd_deaths, nb_queen, pairing = "Deaths")
)

write_csv(geary_results, file.path(table_dir, "q2_global_multivariate_gearys_c.csv"))

# ----------------------------
# 12. SUMMARY TABLES
# ----------------------------
cases_cluster_counts <- q2 %>%
  st_drop_geometry() %>%
  count(cases_cluster_005, sort = TRUE)

deaths_cluster_counts <- q2 %>%
  st_drop_geometry() %>%
  count(deaths_cluster_005, sort = TRUE)

write_csv(cases_cluster_counts, file.path(table_dir, "q2_cases_lisa_cluster_counts.csv"))
write_csv(deaths_cluster_counts, file.path(table_dir, "q2_deaths_lisa_cluster_counts.csv"))

# ----------------------------
# 13. MAPS
# ----------------------------
tmap_mode("plot")

cluster_palette <- c(
  "High-High" = "#b2182b",
  "Low-Low" = "#2166ac",
  "High-Low" = "#ef8a62",
  "Low-High" = "#67a9cf",
  "Not significant" = "grey85"
)

map_cases_lisa <- tm_shape(q2) +
  tm_polygons("cases_cluster_005", palette = cluster_palette, title = "Cases LISA cluster") +
  tm_layout(main.title = "Local bivariate LISA: malaria and SCD cases", legend.outside = TRUE)

map_deaths_lisa <- tm_shape(q2) +
  tm_polygons("deaths_cluster_005", palette = cluster_palette, title = "Deaths LISA cluster") +
  tm_layout(main.title = "Local bivariate LISA: malaria and SCD deaths", legend.outside = TRUE)

tmap_save(map_cases_lisa, file.path(fig_dir, "q2_cases_bivariate_lisa_map.png"), width = 8, height = 6, dpi = 300)
tmap_save(map_deaths_lisa, file.path(fig_dir, "q2_deaths_bivariate_lisa_map.png"), width = 8, height = 6, dpi = 300)

# ----------------------------
# 14. BIVARIATE MORAN SCATTERPLOTS
# ----------------------------
q2 <- q2 %>%
  mutate(
    lag_z_scd_cases = spdep::lag.listw(lw_queen, z_scd_cases, zero.policy = TRUE),
    lag_z_scd_deaths = spdep::lag.listw(lw_queen, z_scd_deaths, zero.policy = TRUE)
  )

make_bivariate_scatter(
  q2,
  x = "z_malaria_cases",
  y = "lag_z_scd_cases",
  title_text = "Bivariate Moran scatterplot: malaria cases vs neighbouring SCD cases",
  xlab = "Standardised malaria cases",
  ylab = "Spatial lag of standardised SCD cases",
  outfile = "q2_cases_bivariate_moran_scatterplot.png"
)

make_bivariate_scatter(
  q2,
  x = "z_malaria_deaths",
  y = "lag_z_scd_deaths",
  title_text = "Bivariate Moran scatterplot: malaria deaths vs neighbouring SCD deaths",
  xlab = "Standardised malaria deaths",
  ylab = "Spatial lag of standardised SCD deaths",
  outfile = "q2_deaths_bivariate_moran_scatterplot.png"
)

# ----------------------------
# 15. SAVE UPDATED SPATIAL OBJECTS
# ----------------------------
sf::st_write(
  q2,
  dsn = file.path(gpkg_dir, "q2_step3_bivariate_outputs.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

saveRDS(q2, file.path(rds_dir, "q2_step3_bivariate_outputs.rds"))

# ----------------------------
# 16. SUMMARY TABLE FOR MANUSCRIPT
# ----------------------------
manuscript_summary <- bind_rows(
  global_biv_results %>%
    transmute(
      analysis = "Global bivariate Moran's I",
      pairing,
      statistic = observed_ib,
      p_value = p_value_two_sided
    ),
  geary_results %>%
    transmute(
      analysis = "Global multivariate Geary's C",
      pairing,
      statistic = observed_c,
      p_value = p_value_two_sided
    )
)

write_csv(manuscript_summary, file.path(table_dir, "q2_step3_manuscript_summary_stats.csv"))

# ----------------------------
# 17. FINAL LOG
# ----------------------------
cat("\nCases LISA cluster counts:\n", file = log_file, append = TRUE)
capture.output(print(cases_cluster_counts), file = log_file, append = TRUE)

cat("\nDeaths LISA cluster counts:\n", file = log_file, append = TRUE)
capture.output(print(deaths_cluster_counts), file = log_file, append = TRUE)

cat("\nStep 3 completed:", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE)

message("Step 3 completed successfully.")
message("Saved:")
message("- q2_global_bivariate_morans_i.csv")
message("- q2_global_multivariate_gearys_c.csv")
message("- q2_local_bivariate_lisa_outputs.csv")
message("- q2_step3_bivariate_outputs.rds")
message("- LISA cluster and significance maps")

# =============================================================================
# HRI Validation Against Copernicus EMS Flood Delineations
# =============================================================================
# Purpose: Validate the Hospitality Risk Index (HRI) against observed flood
#          extents from two Copernicus Emergency Management Service activations:
#          - EMSR409 (Venice coastal flood, November 2019)
#          - EMSR642 (Emilia-Romagna storm surge, November 2022)
#
# Method:  Establishment-level HRI computed per point, then overlaid against
#          observed flood delineation polygons (with 50 m buffer to account for
#          satellite delineation uncertainty) within the EMS Areas of Interest
#          (AOI). Mann-Whitney U test comparing HRI distributions of flood-
#          proximate vs. non-proximate establishments within monitored areas only.
#
# Note:    Spatial density (part of exposure) is inherently H3-level; all other
#          components are computed per establishment.
#
# Author:  Dr. Vilane Gonçalves Sales, ZMT Bremen
# Date:    March 2026
# =============================================================================

library(sf)
library(h3jsr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggspatial)
library(cowplot)
library(scales)

# --- Configuration -----------------------------------------------------------
base_dir <- "C:/Users/vgs/Downloads/SpatialAnalysis"
val_dir  <- file.path(base_dir, "VALIDATION")

# =============================================================================
# 1. LOAD COPERNICUS EMS DATA
# =============================================================================
cat("Loading Copernicus EMS flood delineations and AOI boundaries...\n")

# Flood delineation polygons
venice_flood <- st_read(
  file.path(val_dir, "VENICE/EMSR409_AOI01_DEL_PRODUCT_observedEventA_r1_v1.shp"),
  quiet = TRUE
) %>% st_make_valid() %>% mutate(event = "EMSR409 Venice 2019", region = "Veneto")

comacchio_flood <- st_read(
  file.path(val_dir, "COMACCHIO/EMSR642_AOI02_DEL_PRODUCT_observedEventA_r1_v1.shp"),
  quiet = TRUE
) %>% st_make_valid() %>% mutate(event = "EMSR642 Comacchio 2022", region = "Emilia-Romagna")

ravenna_flood <- st_read(
  file.path(val_dir, "RAVENNA/EMSR642_AOI03_DEL_PRODUCT_observedEventA_r1_v1.shp"),
  quiet = TRUE
) %>% st_make_valid() %>% mutate(event = "EMSR642 Ravenna 2022", region = "Emilia-Romagna")

all_floods <- bind_rows(venice_flood, comacchio_flood, ravenna_flood) %>%
  st_make_valid()

# Area of Interest boundaries (monitoring extent)
venice_aoi <- st_read(
  file.path(val_dir, "VENICE/EMSR409_AOI01_DEL_PRODUCT_areaOfInterestA_r1_v1.shp"),
  quiet = TRUE) %>% st_make_valid()
comacchio_aoi <- st_read(
  file.path(val_dir, "COMACCHIO/EMSR642_AOI02_DEL_PRODUCT_areaOfInterestA_r1_v1.shp"),
  quiet = TRUE) %>% st_make_valid()
ravenna_aoi <- st_read(
  file.path(val_dir, "RAVENNA/EMSR642_AOI03_DEL_PRODUCT_areaOfInterestA_r1_v1.shp"),
  quiet = TRUE) %>% st_make_valid()

all_aoi <- bind_rows(venice_aoi, comacchio_aoi, ravenna_aoi) %>% st_make_valid()
aoi_union <- st_union(st_buffer(all_aoi, 0))

cat("  Venice:    ", nrow(venice_flood), "flood polygons\n")
cat("  Comacchio: ", nrow(comacchio_flood), "flood polygons\n")
cat("  Ravenna:   ", nrow(ravenna_flood), "flood polygons\n")
cat("  Total:     ", nrow(all_floods), "flood polygons across 3 AOIs\n\n")

# =============================================================================
# 2. COMPUTE ESTABLISHMENT-LEVEL HRI
# =============================================================================
cat("Computing establishment-level HRI...\n")

compute_est_hri <- function(df, lon_col, lat_col,
                            alpha = 0.5, beta = 0.5,
                            a_H = 0.3, a_E = 0.35, a_V = 0.35) {

  sf_dat <- st_as_sf(df, coords = c(lon_col, lat_col), crs = 4326)

  # Hazard (SLR-based)
  sf_dat$hazard_score <- case_when(
    sf_dat$mn_slr_Buildings_median_rcp45_2030_2d < 0 ~ 1.0,
    sf_dat$mn_slr_Buildings_median_rcp45_2030_2d < 0.5 ~ 0.8,
    sf_dat$mn_slr_Buildings_median_rcp45_2030_2d < 1.0 ~ 0.6,
    sf_dat$mn_slr_Buildings_median_rcp45_2030_2d < 2.0 ~ 0.4,
    TRUE ~ 0.2
  )

  # Min distance to water
  sf_dat$min_dist_to_water <- pmin(sf_dat$length_sea, sf_dat$length_waterway,
                                    sf_dat$length_estuary, na.rm = TRUE)

  # Physical exposure (quantile-based)
  elev_q1  <- quantile(sf_dat$grond_z, 0.25, na.rm = TRUE)
  elev_med <- median(sf_dat$grond_z, na.rm = TRUE)
  elev_q3  <- quantile(sf_dat$grond_z, 0.75, na.rm = TRUE)
  dist_q1  <- quantile(sf_dat$min_dist_to_water, 0.25, na.rm = TRUE)
  dist_med <- median(sf_dat$min_dist_to_water, na.rm = TRUE)
  dist_q3  <- quantile(sf_dat$min_dist_to_water, 0.75, na.rm = TRUE)

  sf_dat$physical_exposure <- case_when(
    sf_dat$grond_z <= elev_q1 & sf_dat$min_dist_to_water <= dist_q1 ~ 1.0,
    sf_dat$grond_z <= elev_med & sf_dat$min_dist_to_water <= dist_med ~ 0.8,
    sf_dat$grond_z <= elev_q3 | sf_dat$min_dist_to_water <= dist_q3 ~ 0.6,
    sf_dat$grond_z <= 2.5 & sf_dat$min_dist_to_water <= 500 ~ 0.4,
    TRUE ~ 0.2
  )

  # Spatial density (H3-level, joined back to establishments)
  sf_dat$h3_index <- point_to_cell(sf_dat, res = 7)
  density_df <- sf_dat %>% st_drop_geometry() %>%
    group_by(h3_index) %>%
    summarise(est_count = n(), .groups = "drop") %>%
    mutate(spatial_density = log1p(est_count) / log1p(max(est_count)))
  sf_dat <- sf_dat %>%
    left_join(density_df %>% select(h3_index, spatial_density), by = "h3_index")

  # Composite exposure
  sf_dat$exposure_score <- (sf_dat$physical_exposure * 0.7) +
                           (sf_dat$spatial_density * 0.3)

  # Vulnerability (adaptive capacity via type-level economic weight)
  asset_df <- sf_dat %>% st_drop_geometry() %>%
    group_by(encoded_Type_2) %>%
    summarise(count = n(), total_econ = sum(Price, na.rm = TRUE), .groups = "drop") %>%
    mutate(econ_weight = total_econ / sum(total_econ))

  max_w <- 0.6
  scaled_weights <- asset_df$econ_weight * (max_w / max(asset_df$econ_weight))

  ac_factor <- case_when(
    scaled_weights >= quantile(scaled_weights, 0.75) ~ 0.2,
    scaled_weights >= median(scaled_weights) ~ 0.4,
    scaled_weights >= quantile(scaled_weights, 0.25) ~ 0.7,
    TRUE ~ 0.9
  )
  asset_df$ac_factor <- ac_factor

  sf_dat <- sf_dat %>%
    left_join(asset_df %>% select(encoded_Type_2, ac_factor), by = "encoded_Type_2")

  sf_dat <- sf_dat %>% mutate(
    type_weight = case_when(
      encoded_Type_2 %in% 1:8 ~ ac_factor,
      encoded_Type_2 %in% c(9, 10, 11) ~ 0.30,
      encoded_Type_2 %in% c(12, 13, 14) ~ 0.25,
      TRUE ~ 0.20
    ),
    building_area_factor = case_when(
      area > quantile(area, 0.95, na.rm = TRUE) ~ 0.9,
      area < quantile(area, 0.25, na.rm = TRUE) ~ 0.8,
      area > quantile(area, 0.75, na.rm = TRUE) ~ 0.6,
      area < median(area, na.rm = TRUE) ~ 0.5,
      TRUE ~ 0.4
    ),
    building_height_factor = case_when(
      msft_heigh >= quantile(msft_heigh, 0.75, na.rm = TRUE) ~ 0.2,
      msft_heigh >= median(msft_heigh, na.rm = TRUE) ~ 0.4,
      msft_heigh >= quantile(msft_heigh, 0.25, na.rm = TRUE) ~ 0.7,
      TRUE ~ 0.9
    ),
    building_vulnerability = (building_area_factor + building_height_factor) / 2,
    vulnerability_score = type_weight * building_vulnerability
  )

  # Establishment-level HRI (same formula as H3-level)
  sf_dat <- sf_dat %>% mutate(
    product_component  = hazard_score * exposure_score * vulnerability_score,
    additive_component = a_H * hazard_score + a_E * exposure_score + a_V * vulnerability_score,
    raw_hri = alpha * product_component + beta * additive_component
  )

  # Standardize to [0, 1]
  min_hri <- min(sf_dat$raw_hri, na.rm = TRUE)
  max_hri <- max(sf_dat$raw_hri, na.rm = TRUE)
  sf_dat$standardized_hri <- (sf_dat$raw_hri - min_hri) / (max_hri - min_hri + 1e-10)

  return(sf_dat)
}

# Load and process accommodations
df_acc <- read.csv(file.path(base_dir, "comb_buildings_stock_dist_encod_elev.csv")) %>%
  distinct(acc_id, acc_poi_building_id, .keep_all = TRUE)
est_acc <- compute_est_hri(df_acc, "longitude_2", "latitude_2")

# Load and process services
df_ser <- read.csv(file.path(base_dir, "ser_buildings_stock_dist_encod_elev.csv")) %>%
  distinct(ser_buildings_id, .keep_all = TRUE)
est_ser <- compute_est_hri(df_ser, "longitude_3", "latitude_3")

cat("  Accommodations:", nrow(est_acc), "establishments\n")
cat("  Services:      ", nrow(est_ser), "establishments\n\n")

# =============================================================================
# 3. SPATIAL OVERLAY (AOI-RESTRICTED, 50 m BUFFER)
# =============================================================================
cat("Classifying establishments by flood proximity (50 m buffer, AOI-restricted)...\n")

# Create 50 m buffer around flood delineations
# Buffer in projected CRS (EPSG:32632, UTM zone 32N) for metric accuracy
flood_proj <- st_transform(all_floods, 32632)
flood_buf_50 <- st_buffer(flood_proj, 50) %>% st_union() %>% st_make_valid()
flood_buf_50_4326 <- st_transform(flood_buf_50, 4326)

# Which establishments are within EMS monitoring extent?
est_acc$in_aoi  <- st_intersects(est_acc, aoi_union, sparse = FALSE)[, 1]
est_ser$in_aoi  <- st_intersects(est_ser, aoi_union, sparse = FALSE)[, 1]

# Which establishments are within 50 m of flood delineation?
est_acc$flooded <- st_intersects(est_acc, flood_buf_50_4326, sparse = FALSE)[, 1]
est_ser$flooded <- st_intersects(est_ser, flood_buf_50_4326, sparse = FALSE)[, 1]

# Restrict to AOI only
acc_val <- est_acc %>% filter(in_aoi)
ser_val <- est_ser %>% filter(in_aoi)

cat("  Accommodations in AOI:", nrow(acc_val),
    "(flood-proximate:", sum(acc_val$flooded), "| not proximate:", sum(!acc_val$flooded), ")\n")
cat("  Services in AOI:      ", nrow(ser_val),
    "(flood-proximate:", sum(ser_val$flooded), "| not proximate:", sum(!ser_val$flooded), ")\n\n")

# =============================================================================
# 4. STATISTICAL TESTS
# =============================================================================
cat("Running statistical tests...\n\n")

run_est_validation <- function(est_sf, label) {
  cat("=== ", label, " ===\n")
  flooded   <- est_sf %>% filter(flooded) %>% pull(standardized_hri)
  unflooded <- est_sf %>% filter(!flooded) %>% pull(standardized_hri)

  cat("  Flood-proximate: n =", length(flooded), "\n")
  cat("    Median HRI: ", round(median(flooded), 4), "\n")
  cat("    Mean HRI:   ", round(mean(flooded), 4), "\n")
  cat("    SD:         ", round(sd(flooded), 4), "\n")

  cat("  Not proximate:   n =", length(unflooded), "\n")
  cat("    Median HRI: ", round(median(unflooded), 4), "\n")
  cat("    Mean HRI:   ", round(mean(unflooded), 4), "\n")
  cat("    SD:         ", round(sd(unflooded), 4), "\n")

  mw <- wilcox.test(flooded, unflooded, alternative = "greater", exact = FALSE)
  cat("\n  Mann-Whitney U (H1: proximate > not proximate):\n")
  cat("    W =", mw$statistic, "\n")
  cat("    p =", format.pval(mw$p.value, digits = 4), "\n")
  n1 <- length(flooded); n2 <- length(unflooded)
  r_effect <- mw$statistic / (n1 * n2) * 2 - 1
  cat("    Rank-biserial r =", round(r_effect, 4), "\n")

  cat("\n  Components (proximate vs not-proximate means):\n")
  for (comp in c("hazard_score", "physical_exposure", "exposure_score",
                 "vulnerability_score")) {
    f_val <- mean(est_sf %>% filter(flooded) %>% pull(!!sym(comp)), na.rm = TRUE)
    u_val <- mean(est_sf %>% filter(!flooded) %>% pull(!!sym(comp)), na.rm = TRUE)
    cat("    ", comp, ": F =", round(f_val, 4), "| NF =", round(u_val, 4),
        "| diff =", round(f_val - u_val, 4), "\n")
  }
  cat("\n")
  return(mw)
}

mw_acc <- run_est_validation(acc_val, "Accommodations (50 m buffer, AOI)")
mw_ser <- run_est_validation(ser_val, "Services (50 m buffer, AOI)")

# Pooled
both_val <- bind_rows(
  acc_val %>% select(h3_index, standardized_hri, flooded,
                     hazard_score, physical_exposure, exposure_score,
                     vulnerability_score) %>% mutate(sector = "Accommodations"),
  ser_val %>% select(h3_index, standardized_hri, flooded,
                     hazard_score, physical_exposure, exposure_score,
                     vulnerability_score) %>% mutate(sector = "Services")
)
mw_pooled <- run_est_validation(both_val, "Pooled (50 m buffer, AOI)")

# =============================================================================
# 5. SAVE RESULTS
# =============================================================================
cat("Saving results...\n")

est_results <- bind_rows(
  acc_val %>% st_drop_geometry() %>%
    select(h3_index, standardized_hri, raw_hri, flooded,
           hazard_score, physical_exposure, exposure_score,
           vulnerability_score) %>%
    mutate(sector = "Accommodations"),
  ser_val %>% st_drop_geometry() %>%
    select(h3_index, standardized_hri, raw_hri, flooded,
           hazard_score, physical_exposure, exposure_score,
           vulnerability_score) %>%
    mutate(sector = "Services")
)
write.csv(est_results,
          file.path(val_dir, "validation_50m_flood_overlay.csv"),
          row.names = FALSE)

est_summary <- est_results %>%
  group_by(sector, flooded) %>%
  summarise(
    n = n(),
    mean_hri = round(mean(standardized_hri), 4),
    median_hri = round(median(standardized_hri), 4),
    sd_hri = round(sd(standardized_hri), 4),
    mean_hazard = round(mean(hazard_score), 4),
    mean_phys_exposure = round(mean(physical_exposure), 4),
    mean_exposure = round(mean(exposure_score), 4),
    mean_vulnerability = round(mean(vulnerability_score), 4),
    .groups = "drop"
  )
write.csv(est_summary,
          file.path(val_dir, "validation_50m_summary.csv"),
          row.names = FALSE)

cat("  Saved: validation_50m_flood_overlay.csv\n")
cat("  Saved: validation_50m_summary.csv\n\n")

# =============================================================================
# 6. PUBLICATION-QUALITY FIGURES (Popovic style)
# =============================================================================
cat("Creating figures...\n")

theme_pub_map <- function() {
  theme_void() + theme(
    plot.title    = element_text(face = "bold", size = 10, hjust = 0.5),
    plot.subtitle = element_text(size = 7, hjust = 0.5, color = "grey40"),
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 8),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(0.4, "cm"),
    legend.margin     = margin(t = 2, b = 2),
    plot.margin       = margin(t = 5, r = 5, b = 5, l = 5)
  )
}

# Bounding box from AOI extent
aoi_bbox <- st_bbox(aoi_union)
margin <- 0.03
aoi_bbox_exp <- aoi_bbox
aoi_bbox_exp["xmin"] <- aoi_bbox["xmin"] - margin
aoi_bbox_exp["ymin"] <- aoi_bbox["ymin"] - margin
aoi_bbox_exp["xmax"] <- aoi_bbox["xmax"] + margin
aoi_bbox_exp["ymax"] <- aoi_bbox["ymax"] + margin

# --- Figure 1: Establishment HRI + flood overlay (2 panels) ---
flood_buf_sf <- st_as_sf(st_cast(flood_buf_50_4326, "MULTIPOLYGON"))

p_acc_map <- ggplot() +
  annotation_map_tile(type = "cartolight", zoom = 11, quiet = TRUE) +
  geom_sf(data = flood_buf_sf, fill = "#fdae61", color = "#fdae61",
          alpha = 0.2, linewidth = 0.2) +
  geom_sf(data = all_floods, fill = "#d7191c", color = "#d7191c",
          alpha = 0.3, linewidth = 0.3) +
  geom_sf(data = acc_val %>% filter(!flooded),
          color = "#2c7bb6", size = 0.4, alpha = 0.5) +
  geom_sf(data = acc_val %>% filter(flooded),
          color = "#d7191c", size = 0.6, alpha = 0.7) +
  coord_sf(xlim = c(aoi_bbox_exp["xmin"], aoi_bbox_exp["xmax"]),
           ylim = c(aoi_bbox_exp["ymin"], aoi_bbox_exp["ymax"]), crs = 4326) +
  annotation_scale(location = "br", width_hint = 0.2, text_col = "black",
                   line_col = "black", bar_cols = c("black", "white"), text_cex = 0.6) +
  labs(title = paste0("Accommodations (n = ", format(nrow(acc_val), big.mark = ","), ")"),
       subtitle = "Red = observed flood | Orange = 50 m buffer") +
  theme_pub_map()

p_ser_map <- ggplot() +
  annotation_map_tile(type = "cartolight", zoom = 11, quiet = TRUE) +
  geom_sf(data = flood_buf_sf, fill = "#fdae61", color = "#fdae61",
          alpha = 0.2, linewidth = 0.2) +
  geom_sf(data = all_floods, fill = "#d7191c", color = "#d7191c",
          alpha = 0.3, linewidth = 0.3) +
  geom_sf(data = ser_val %>% filter(!flooded),
          color = "#2c7bb6", size = 0.4, alpha = 0.5) +
  geom_sf(data = ser_val %>% filter(flooded),
          color = "#d7191c", size = 0.6, alpha = 0.7) +
  coord_sf(xlim = c(aoi_bbox_exp["xmin"], aoi_bbox_exp["xmax"]),
           ylim = c(aoi_bbox_exp["ymin"], aoi_bbox_exp["ymax"]), crs = 4326) +
  annotation_scale(location = "br", width_hint = 0.2, text_col = "black",
                   line_col = "black", bar_cols = c("black", "white"), text_cex = 0.6) +
  labs(title = paste0("Services (n = ", format(nrow(ser_val), big.mark = ","), ")"),
       subtitle = "Red = observed flood | Orange = 50 m buffer") +
  theme_pub_map()

p_map_composite <- plot_grid(p_acc_map, p_ser_map, ncol = 2,
                             labels = c("A", "B"), label_size = 12)
ggsave(file.path(val_dir, "validation_50m_map_overlay.png"),
       p_map_composite, width = 14, height = 8, dpi = 300, bg = "white")

# --- Figure 2: HRI box plots by flood status ---
est_boxdata <- bind_rows(
  acc_val %>% st_drop_geometry() %>%
    select(standardized_hri, flooded) %>% mutate(sector = "Accommodations"),
  ser_val %>% st_drop_geometry() %>%
    select(standardized_hri, flooded) %>% mutate(sector = "Services")
) %>% mutate(flood_label = ifelse(flooded, "Flood-proximate\n(\u2264 50 m)", "Not flood-proximate"))

p_est_box <- ggplot(est_boxdata,
                    aes(x = flood_label, y = standardized_hri, fill = flood_label)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.4) +
  facet_wrap(~sector) +
  scale_fill_manual(values = c("Flood-proximate\n(\u2264 50 m)" = "#d7191c", "Not flood-proximate" = "#2c7bb6"),
                    guide = "none") +
  labs(x = NULL, y = "Establishment-Level Standardized HRI",
       title = "HRI by Flood Proximity (50 m Buffer, AOI-Restricted)",
       subtitle = "Mann-Whitney U, one-sided: Acc p < 0.001, Ser p < 0.001") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
    strip.text = element_text(face = "bold", size = 9),
    axis.title.y = element_text(size = 9),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(val_dir, "validation_50m_boxplot_hri.png"),
       p_est_box, width = 8, height = 5, dpi = 300, bg = "white")

# --- Figure 3: Component-level box plots ---
est_comp <- bind_rows(
  acc_val %>% st_drop_geometry() %>%
    select(flooded, hazard_score, physical_exposure, exposure_score,
           vulnerability_score) %>% mutate(sector = "Accommodations"),
  ser_val %>% st_drop_geometry() %>%
    select(flooded, hazard_score, physical_exposure, exposure_score,
           vulnerability_score) %>% mutate(sector = "Services")
) %>%
  pivot_longer(cols = c(hazard_score, physical_exposure, exposure_score,
                        vulnerability_score),
               names_to = "component", values_to = "value") %>%
  mutate(
    flood_label = ifelse(flooded, "Flood-proximate", "Not proximate"),
    component = case_when(
      component == "hazard_score" ~ "Hazard",
      component == "physical_exposure" ~ "Physical Exposure",
      component == "exposure_score" ~ "Composite Exposure",
      component == "vulnerability_score" ~ "Vulnerability"
    ),
    component = factor(component, levels = c("Hazard", "Physical Exposure",
                                              "Composite Exposure", "Vulnerability"))
  )

p_est_comp <- ggplot(est_comp,
                     aes(x = flood_label, y = value, fill = flood_label)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8, width = 0.6) +
  facet_grid(component ~ sector, scales = "free_y") +
  scale_fill_manual(values = c("Flood-proximate" = "#d7191c", "Not proximate" = "#2c7bb6"),
                    guide = "none") +
  labs(x = NULL, y = "Component Score",
       title = "HRI Components by Flood Proximity (50 m Buffer, AOI-Restricted)") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 9),
    axis.title.y = element_text(size = 9),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(val_dir, "validation_50m_boxplot_components.png"),
       p_est_comp, width = 8, height = 9, dpi = 300, bg = "white")

cat("  Saved: validation_50m_map_overlay.png\n")
cat("  Saved: validation_50m_boxplot_hri.png\n")
cat("  Saved: validation_50m_boxplot_components.png\n\n")

cat("=== Validation complete ===\n")

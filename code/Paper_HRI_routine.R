##############################################################################
# Unified Analysis Script: Fuzzy Clustering, IPCC-Aligned HRI, and Synergy Analysis
# Date: 2025-01-22
# Reproducibility: package versions recorded in renv.lock
# To restore: renv::restore()
##############################################################################

# ------------------------------
# 1. LIBRARIES & INITIAL SETUP
# ------------------------------
# Load necessary libraries
library(sf)            # Spatial data handling
library(h3jsr)         # H3 operations
library(ggplot2)       # Plotting
library(dplyr)         # Data manipulation
library(tidyr)         # Data tidying
library(spatstat)      # Kernel density
library(spdep)         # Spatial dependencies
library(ggspatial)     # Scale bars, north arrows
library(viridis)       # Color scales
library(fclust)        # Fuzzy clustering
library(RColorBrewer)  # Color palettes
library(randomForest)  # Random forest interpretability
library(pdp)           # Partial Dependence Plot
library(gridExtra)     # Arranging multiple plots
library(corrplot)      # Correlation plots
library(here)          # Reproducible file paths
library(patchwork)     # Multi-panel plot composition

# Set working directory using here::here() for reproducibility
setwd(here::here("Acc_new"))

# ==============================================================================
# REPRODUCIBILITY SETTINGS
# ==============================================================================
# Set global seed for reproducible results across all random operations
set.seed(42)
cat("🎯 Random seed set to 42 for reproducible results\n")

# ------------------------------
# 2. DATA LOADING & PREPROCESSING
# ------------------------------
# Load primary dataset
df <- read.csv(here::here("comb_buildings_stock_dist_encod_elev.csv"))

# Ensure distinct entries
df_raw <- df %>%
  distinct(acc_id, acc_poi_building_id, .keep_all = TRUE)

# Convert to spatial data (sf)
sf_data <- st_as_sf(df_raw, coords = c("longitude_2", "latitude_2"), crs = 4326)

# Load and preprocess Italy boundaries
ita_boundaries <- readRDS(here::here("gadm36_ITA_1_sp.rds"))
if(!inherits(ita_boundaries, "sf")) {
  ita_boundaries <- st_as_sf(ita_boundaries)
}

# Define bounding box with expanded margins
data_bbox <- st_bbox(sf_data)
expand_margin <- 0.09
data_bbox_expanded <- data_bbox
data_bbox_expanded["xmin"] <- data_bbox_expanded["xmin"] - expand_margin
data_bbox_expanded["ymin"] <- data_bbox_expanded["ymin"] - expand_margin
data_bbox_expanded["xmax"] <- data_bbox_expanded["xmax"] + expand_margin
data_bbox_expanded["ymax"] <- data_bbox_expanded["ymax"] + expand_margin

# Crop Italy boundaries
ita_boundaries_cropped <- st_crop(ita_boundaries, data_bbox_expanded)

# ------------------------------
# 3. HELPER FUNCTIONS
# ------------------------------
# F1 Score calculation function
compute_f1_scores <- function(conf_mat) {
  conf_mat <- as.matrix(conf_mat)
  classes <- colnames(conf_mat)
  metrics_list <- list()
  
  for (cls in classes) {
    tp <- conf_mat[cls, cls]
    fn <- sum(conf_mat[, cls]) - tp
    fp <- sum(conf_mat[cls, ]) - tp
    precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    recall    <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    f1        <- if ((precision + recall) == 0) 0 else (2 * precision * recall) / (precision + recall)
    
    metrics_list[[cls]] <- data.frame(
      Class     = cls,
      Precision = precision,
      Recall    = recall,
      F1        = f1
    )
  }
  
  metrics_df <- do.call(rbind, metrics_list)
  rownames(metrics_df) <- NULL
  return(metrics_df)
}

# ------------------------------
# 4. FUZZY CLUSTERING FUNCTIONS
# ------------------------------
fuzzy_cluster_analysis_features <- function(sf_points, extra_features, k = 3, m = 1.5) {
  # Extract coordinates
  coords_mat <- st_coordinates(sf_points)
  
  # Extract relevant columns
  df <- st_drop_geometry(sf_points)
  sub_df <- df[, extra_features, drop=FALSE]
  
  # Convert -1 to NA
  for (colname in extra_features) {
    is_neg1 <- !is.na(sub_df[[colname]]) & sub_df[[colname]] == -1
    sub_df[[colname]][is_neg1] <- NA
  }
  
  # Combine coords + features
  full_mat <- cbind(coords_mat, sub_df)
  
  # Filter only complete rows
  keep <- complete.cases(full_mat)
  final_mat <- full_mat[keep, , drop=FALSE]
  
  # Scale the data
  scaled_mat <- scale(final_mat)
  
  # Perform fuzzy clustering (with reproducible seed)
  set.seed(42)
  fkm_res <- FKM(X = scaled_mat, k = k, m = m)
  
  # Validity metrics
  pc_val <- PC(fkm_res$U)
  pe_val <- PE(fkm_res$U)
  xb_val <- XB(fkm_res$X, fkm_res$U, fkm_res$H)
  
  list(
    fuzzy_result = fkm_res,
    keep_vector = keep,
    partition_coefficient = pc_val,
    partition_entropy = pe_val,
    xb_index = xb_val,
    scaled_data = scaled_mat
  )
}

posthoc_fuzzy_cluster_analysis <- function(fuzzy_result_list, sf_points, extra_features) {
  U_matrix <- fuzzy_result_list$fuzzy_result$U
  hard_cluster <- apply(U_matrix, 1, which.max)
  
  keep <- fuzzy_result_list$keep_vector
  sf_points_kept <- sf_points[keep, ]
  sf_points_kept$cluster_label <- factor(hard_cluster)
  
  # Summaries
  summary_list <- list()
  for(feature in extra_features) {
    summary_by_cluster <- sf_points_kept %>%
      st_drop_geometry() %>%
      group_by(cluster_label) %>%
      summarize(
        mean_value = mean(.data[[feature]], na.rm=TRUE),
        median_value = median(.data[[feature]], na.rm=TRUE),
        n_points = n(),
        .groups = "drop"
      )
    summary_list[[feature]] <- summary_by_cluster
  }
  
  list(
    labeled_data = sf_points_kept,
    feature_summaries = summary_list,
    membership_matrix = U_matrix
  )
}

# Export functions
export_feature_summaries <- function(posthoc_result, txt_filename) {
  all_features_df <- do.call(rbind, lapply(names(posthoc_result$feature_summaries), function(f) {
    s <- posthoc_result$feature_summaries[[f]]
    s$feature <- f
    s
  }))
  write.table(
    all_features_df,
    file = txt_filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  cat("Exported feature summaries to", txt_filename, "\n")
}

export_validities <- function(fuzzy_result_list, csv_filename) {
  df <- data.frame(
    partition_coefficient = fuzzy_result_list$partition_coefficient,
    partition_entropy = fuzzy_result_list$partition_entropy,
    xb_index = fuzzy_result_list$xb_index
  )
  write.csv(df, csv_filename, row.names=FALSE)
  cat("Exported validity metrics to", csv_filename, "\n")
}

# ------------------------------
# 5. HRI CALCULATION FUNCTIONS
# ------------------------------
calculate_hri <- function(sf_data, 
                          alpha = 0.5,       # weight for the multiplicative term
                          beta = 0.5,        # weight for the additive term
                          a_H = 0.3,         # additive weight for hazard
                          a_E = 0.35,         # additive weight for exposure
                          a_V = 0.35) {       # additive weight for vulnerability
  
  # 5A. Hazard Component
  sf_data <- sf_data %>%
    mutate(
      hazard_score = case_when(
        mn_slr_Buildings_median_rcp45_2030_2d < 0 ~ 1.0,    # Submerged
        mn_slr_Buildings_median_rcp45_2030_2d < 0.5 ~ 0.8,  # Critical
        mn_slr_Buildings_median_rcp45_2030_2d < 1.0 ~ 0.6,  # High
        mn_slr_Buildings_median_rcp45_2030_2d < 2.0 ~ 0.4,  # Moderate
        TRUE ~ 0.2                                       # Low
      )
    )
  
  # 5B. Initial Distance Calculations
  sf_data <- sf_data %>%
    rowwise() %>%
    mutate(
      min_dist_to_water = min(c(length_sea, length_waterway, length_estuary), na.rm = TRUE)
    ) %>%
    ungroup()
  
  # 5C. Exposure Component
  physical_analysis <- sf_data %>%
    summarize(
      elev_q1 = quantile(grond_z, 0.25, na.rm = TRUE),
      elev_median = median(grond_z, na.rm = TRUE),
      elev_q3 = quantile(grond_z, 0.75, na.rm = TRUE),
      dist_q1 = quantile(min_dist_to_water, 0.25, na.rm = TRUE),
      dist_median = median(min_dist_to_water, na.rm = TRUE),
      dist_q3 = quantile(min_dist_to_water, 0.75, na.rm = TRUE)
    )
  
  sf_data <- sf_data %>%
    rowwise() %>%
    mutate(
      physical_exposure = case_when(
        grond_z <= physical_analysis$elev_q1 & 
          min_dist_to_water <= physical_analysis$dist_q1 ~ 1.0,
        grond_z <= physical_analysis$elev_median & 
          min_dist_to_water <= physical_analysis$dist_median ~ 0.8,
        grond_z <= physical_analysis$elev_q3 | 
          min_dist_to_water <= physical_analysis$dist_q3 ~ 0.6,
        grond_z <= 2.5 & min_dist_to_water <= 500 ~ 0.4,
        TRUE ~ 0.2
      )
    ) %>%
    ungroup()
  
  # Calculate spatial density
  if (!"h3_index" %in% names(sf_data)) {
    sf_data$h3_index <- point_to_cell(sf_data, res = 7)
  }
  
  density_df <- sf_data %>%
    st_drop_geometry() %>%
    group_by(h3_index) %>%
    summarize(
      establishment_count = n(),
      .groups = "drop"
    ) %>%
    mutate(
      max_establishments = max(establishment_count, na.rm = TRUE),
      spatial_density = log1p(establishment_count) / log1p(max_establishments)
    )
  
  sf_data <- sf_data %>%
    left_join(density_df %>% select(h3_index, spatial_density),
              by = "h3_index")
  
  # 5D. Vulnerability Component
  asset_analysis <- sf_data %>%
    group_by(encoded_Type_2) %>%
    summarize(
      count = n(),
      avg_area = mean(area, na.rm = TRUE),
      avg_height = mean(msft_heigh, na.rm = TRUE),
      total_economic_value = sum(Price, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      economic_weight = total_economic_value / sum(total_economic_value)
    )
  
  max_weight <- 0.6
  scaled_economic_weights <- asset_analysis$economic_weight * 
    (max_weight / max(asset_analysis$economic_weight))
  
  adaptive_capacity_factor = case_when(
    scaled_economic_weights >= quantile(scaled_economic_weights, 0.75) ~ 0.2,  # High price = low vulnerability
    scaled_economic_weights >= median(scaled_economic_weights) ~ 0.4,
    scaled_economic_weights >= quantile(scaled_economic_weights, 0.25) ~ 0.7,
    TRUE ~ 0.9  # Low price = high vulnerability
  )
  
  sf_data <- sf_data %>%
    mutate(
      type_weight = case_when(
        encoded_Type_2 %in% 1:8 ~ adaptive_capacity_factor[encoded_Type_2],
        encoded_Type_2 %in% c(9,10,11) ~ 0.30,
        encoded_Type_2 %in% c(12,13,14) ~ 0.25,
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
  
  # 5E. Final HRI Calculation with combined multiplicative and additive components
  hri_df <- sf_data %>%
    st_drop_geometry() %>%
    group_by(h3_index) %>%
    summarize(
      # Aggregating components for each h3 cell:
      phys_exposure_avg = mean(physical_exposure, na.rm = TRUE),
      phys_exposure_sd = sd(physical_exposure, na.rm = TRUE),
      phys_exposure_count = sum(!is.na(physical_exposure)),
      spatial_density_avg = first(spatial_density),
      
      vulnerability_avg = mean(vulnerability_score, na.rm = TRUE),
      vulnerability_sd = sd(vulnerability_score, na.rm = TRUE),
      vulnerability_count = sum(!is.na(vulnerability_score)),
      
      hazard_avg = mean(hazard_score, na.rm = TRUE),
      hazard_sd = sd(hazard_score, na.rm = TRUE),
      hazard_count = sum(!is.na(hazard_score)),
      
      .groups = "drop"
    ) %>%
    mutate(
      # Here we already have the additive exposure score as before:
      exposure_score = (phys_exposure_avg * 0.7) + (spatial_density_avg * 0.3),
      
      # Compute standard errors (these are approximations assuming independence):
      exposure_se = ifelse(phys_exposure_count > 1,
                           phys_exposure_sd / sqrt(phys_exposure_count),
                           NA),
      vulnerability_se = ifelse(vulnerability_count > 1,
                                vulnerability_sd / sqrt(vulnerability_count),
                                NA),
      hazard_se = ifelse(hazard_count > 1,
                         hazard_sd / sqrt(hazard_count),
                         NA),
      
      # Compute the multiplicative component and its standard error via propagation:
      product_component = hazard_avg * exposure_score * vulnerability_avg,
      product_se = product_component * sqrt((hazard_se / hazard_avg)^2 +
                                              (exposure_se / exposure_score)^2 +
                                              (vulnerability_se / vulnerability_avg)^2),
      
      # Compute the additive component and its standard error:
      additive_component = a_H * hazard_avg + a_E * exposure_score + a_V * vulnerability_avg,
      additive_se = sqrt((a_H * hazard_se)^2 + (a_E * exposure_se)^2 + (a_V * vulnerability_se)^2),
      
      # Combine both components to get the raw HRI:
      raw_hri = alpha * product_component + beta * additive_component,
      
      # Propagate the errors (assuming independence between the two components):
      hri_se = sqrt((alpha * product_se)^2 + (beta * additive_se)^2),
      
      # Coefficient of variation for reporting purposes:
      hri_cv = (hri_se / raw_hri) * 100,
      
      uncertainty_level = case_when(
        is.na(hri_cv) ~ "No Data",
        hri_cv <= 5 ~ "Low",
        hri_cv <= 10 ~ "Medium",
        TRUE ~ "High"
      )
    )
  
  # Standardize HRI over the groups
  min_hri <- min(hri_df$raw_hri, na.rm = TRUE)
  max_hri <- max(hri_df$raw_hri, na.rm = TRUE)
  epsilon <- 1e-10
  
  hri_df <- hri_df %>%
    mutate(
      standardized_hri = (raw_hri - min_hri) / (max_hri - min_hri + epsilon)
    )
  
  return(list(sf_data = sf_data, hri_df = hri_df))
}

create_hri_category <- function(hri_df, n_breaks = 5) {
  # Create quintile integer (1 to 5)
  hri_df <- hri_df %>%
    mutate(hri_quintile = ntile(standardized_hri, n_breaks)) %>%
    mutate(
      hri_category = factor(
        case_when(
          hri_quintile == 5 ~ "Very High",
          hri_quintile == 4 ~ "High",
          hri_quintile == 3 ~ "Medium",
          hri_quintile == 2 ~ "Low",
          hri_quintile == 1 ~ "Very Low"
        ),
        levels = c("Very Low", "Low", "Medium", "High", "Very High")
      )
    )
  return(hri_df)
}

# ------------------------------
# 6. SYNERGY ANALYSIS FUNCTIONS
# ------------------------------
# Enhanced synergy analysis function with plots and geometry handling
analyze_cluster_hri_synergy <- function(
    synergy_df, 
    hri_df, 
    scenario_name, 
    filter_na = TRUE
) {
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  
  # ------------------------------
  # 1. FACTOR ORDERING FOR CLUSTER LABEL
  # ------------------------------
  unique_clusters <- synergy_df$cluster_label %>%
    unique() %>%
    as.character() %>%
    sort()
  
  synergy_df$cluster_label <- factor(
    synergy_df$cluster_label,
    levels = unique_clusters
  )
  
  # ------------------------------
  # 2. AGGREGATION: DOMINANCE + PROPORTION
  # ------------------------------
  cluster_h3_aggregation <- synergy_df %>%
    group_by(h3_index, cluster_label) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(h3_index) %>%
    mutate(
      total = sum(count),
      prop = count / total
    ) %>%
    ungroup()
  
  # Identify dominant cluster in each H3 cell
  dominant_cluster_h3 <- cluster_h3_aggregation %>%
    group_by(h3_index) %>%
    slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(
      dominant_cluster = cluster_label,
      dominant_count   = count,
      dominant_prop    = prop
    )
  
  # Merge with HRI
  synergy_h3_df <- hri_df %>%
    left_join(dominant_cluster_h3, by = "h3_index")
  
  if (filter_na) {
    synergy_h3_df <- synergy_h3_df %>%
      filter(!is.na(dominant_cluster))
  }
  
  # Factor ordering for synergy_h3_df$dominant_cluster
  unique_dom_clusters <- synergy_h3_df$dominant_cluster %>%
    unique() %>%
    as.character() %>%
    sort()
  
  synergy_h3_df$dominant_cluster <- factor(
    synergy_h3_df$dominant_cluster,
    levels = unique_dom_clusters
  )
  
  # ------------------------------
  # 3. DOMINANCE-BASED SUMMARIES
  # ------------------------------
  cluster_hri_summary <- synergy_h3_df %>%
    group_by(dominant_cluster, hri_category) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(dominant_cluster) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup() %>%
    arrange(dominant_cluster, hri_category)
  
  cluster_uncertainty_summary <- synergy_h3_df %>%
    group_by(dominant_cluster, uncertainty_level) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(dominant_cluster) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup() %>%
    arrange(dominant_cluster, uncertainty_level)
  
  cluster_stats <- synergy_h3_df %>%
    group_by(dominant_cluster) %>%
    summarize(
      mean_stdHRI = mean(standardized_hri, na.rm=TRUE),
      sd_stdHRI   = sd(standardized_hri, na.rm=TRUE),
      mean_rawHRI = mean(raw_hri, na.rm=TRUE),
      mean_cv     = mean(hri_cv, na.rm=TRUE),
      mean_se     = mean(hri_se, na.rm=TRUE),
      n_h3_cells  = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_stdHRI))
  
  cluster_size <- synergy_h3_df %>%
    group_by(dominant_cluster) %>%
    summarize(
      total_h3_cells = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(total_h3_cells))
  
  # ------------------------------
  # 4. PROPORTION-BASED SUMMARIES
  # ------------------------------
  proportion_h3_df <- cluster_h3_aggregation %>%
    left_join(hri_df, by = "h3_index")
  
  if (filter_na) {
    proportion_h3_df <- proportion_h3_df %>%
      filter(!is.na(cluster_label))
  }
  
  proportion_h3_df$cluster_label <- factor(
    proportion_h3_df$cluster_label,
    levels = unique_clusters
  )
  
  cluster_proportion_summary <- proportion_h3_df %>%
    group_by(cluster_label, hri_category) %>%
    summarize(
      mean_prop    = mean(prop, na.rm=TRUE),
      total_points = sum(count, na.rm=TRUE),
      .groups      = "drop"
    ) %>%
    arrange(cluster_label, hri_category)
  
  cluster_uncertainty_prop_summary <- proportion_h3_df %>%
    group_by(cluster_label, uncertainty_level) %>%
    summarize(
      mean_prop    = mean(prop, na.rm=TRUE),
      total_points = sum(count, na.rm=TRUE),
      .groups      = "drop"
    ) %>%
    arrange(cluster_label, uncertainty_level)
  
  # ------------------------------
  # 5. DOMINANCE PLOTS (Combined into one PNG)
  # ------------------------------
  # (A) Cluster Size
  p_cluster_size <- cluster_size %>%
    ggplot(aes(
      x = dominant_cluster,
      y = total_h3_cells,
      fill = dominant_cluster
    )) +
    geom_bar(stat="identity") +
    scale_fill_viridis_d(option="inferno", name="Cluster") +
    labs(
      title = paste("Cluster Size (Membership Dominance):", scenario_name),
      x = "Cluster", y = "Number of Establishments"
    ) +
    theme_minimal() +
    theme(legend.position="none")
  
  # (B) HRI Category Distribution by Dominant Cluster
  p_hri_bar <- cluster_hri_summary %>%
    ggplot(aes(x = dominant_cluster, y = prop, fill = hri_category)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_viridis_d(direction=-1, option="plasma", name="HRI Category") +
    labs(
      title = paste("HRI Category Distribution by Cluster:", scenario_name),
      x = "Cluster", y = "Proportion"
    ) +
    theme_minimal()
  
  # (C) Box Plot: Standardized HRI by Dominant Cluster
  p_hri_box <- synergy_h3_df %>%
    ggplot(aes(x = dominant_cluster, y = standardized_hri, fill = dominant_cluster)) +
    geom_boxplot(alpha=0.6) +
    scale_fill_viridis_d(option="C", name="Cluster") +
    labs(
      title = paste("Standardized HRI by Cluster:", scenario_name),
      x = "Cluster", y = "Standardized HRI"
    ) +
    theme_minimal() +
    theme(legend.position="none")
  
  # (D) HRI Uncertainty Distribution by Dominant Cluster
  p_uncert_bar <- cluster_uncertainty_summary %>%
    ggplot(aes(x = dominant_cluster, y = prop, fill = uncertainty_level)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_viridis_d(name="Uncertainty Level") +
    labs(
      title = paste("HRI Uncertainty Distribution by Cluster:", scenario_name),
      x = "Cluster", y = "Proportion"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  # Combine these 4 plots vertically
  combined_dominance_plots <- gridExtra::grid.arrange(
    p_cluster_size,
    p_hri_bar,
    p_hri_box,
    p_uncert_bar,
    ncol=1
  )
  
  # Save the combined plot
  ggsave(
    filename = paste0("synergy_", scenario_name, "_all_plots.png"),
    plot = combined_dominance_plots,
    width=8,
    height=20,  # Adjust as needed
    dpi=300
  )
  
  cat("\n=== Single Combined Dominance Plot Saved ===\n")
  
  # ------------------------------
  # 6. PROPORTION PLOTS (Separate PNG)
  # ------------------------------
  # (A) HRI Category Distribution by Proportion
  p_prop_hri <- cluster_proportion_summary %>%
    ggplot(aes(x=cluster_label, y=mean_prop, fill=hri_category)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_viridis_d(direction=-1, option="plasma", name="HRI Category") +
    labs(
      title = paste("HRI Category Distribution by Cluster Proportion:", scenario_name),
      x = "Cluster",
      y = "Avg. Proportion (H3 Cells)"
    ) +
    theme_minimal()
  
  # (B) HRI Uncertainty Distribution by Proportion
  p_prop_uncert <- cluster_uncertainty_prop_summary %>%
    ggplot(aes(x=cluster_label, y=mean_prop, fill=uncertainty_level)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_viridis_d(name="Uncertainty Level") +
    labs(
      title = paste("HRI Uncertainty Distribution by Cluster Proportion:", scenario_name),
      x = "Cluster", y = "Avg. Proportion (H3 Cells)"
    ) +
    theme_minimal()
  
  # Combine these 2 plots vertically
  combined_proportion_plots <- gridExtra::grid.arrange(
    p_prop_hri,
    p_prop_uncert,
    ncol = 1
  )
  
  # Save the combined proportion plot
  ggsave(
    filename = paste0("synergy_", scenario_name, "_proportion_plots.png"),
    plot = combined_proportion_plots,
    width=8,
    height=12,
    dpi=300
  )
  
  cat("\n=== Proportion Plots Saved ===\n")
  
  # ------------------------------
  # 7. STATISTICAL TESTING
  # ------------------------------
  # Perform ANOVA
  anova_result <- aov(standardized_hri ~ dominant_cluster, data = synergy_h3_df)
  cat("\n--- ANOVA (Dominant Clusters): Standardized HRI ---\n")
  capture.output(summary(anova_result), file = paste0("synergy_", scenario_name, "_anova_results.txt"))
  
  # Perform Kruskal-Wallis test
  kruskal_result <- kruskal.test(standardized_hri ~ dominant_cluster, data = synergy_h3_df)
  cat("\n--- Kruskal-Wallis (Dominant Clusters): Standardized HRI ---\n")
  capture.output(kruskal_result, file = paste0("synergy_", scenario_name, "_kruskal_results.txt"))
  
  # ------------------------------
  # 8. RETURN SUMMARIES
  # ------------------------------
  
  return(list(
    synergy_h3_df                = synergy_h3_df,
    
    # Dominance-based
    cluster_hri_summary          = cluster_hri_summary,
    cluster_uncertainty_summary  = cluster_uncertainty_summary,
    cluster_stats                = cluster_stats,
    cluster_size                 = cluster_size,
    
    # Proportion-based
    cluster_proportion_summary        = cluster_proportion_summary,
    cluster_uncertainty_prop_summary  = cluster_uncertainty_prop_summary
  ))
}

# ------------------------------
# 7. ADVANCED STATISTICAL ANALYSIS FUNCTIONS
# ------------------------------
# Ensure necessary libraries are loaded
library(spatstat)
library(spdep)
library(stats)

# ------------------------------
# 7A. SPATIAL AUTOCORRELATION
# ------------------------------
spatial_analysis <- function(h3_sf, scenario_name) {
  library(spatstat)
  library(spdep)
  library(dplyr)
  
  if (!inherits(h3_sf, "sf")) stop("h3_sf must be an sf object.")
  if (!"standardized_hri" %in% names(h3_sf)) stop("The 'standardized_hri' column is missing.")
  
  # Clean data first
  h3_sf <- h3_sf %>% 
    filter(!is.na(standardized_hri), is.finite(standardized_hri)) %>%
    mutate(standardized_hri = as.numeric(standardized_hri))
  
  if (var(h3_sf$standardized_hri, na.rm = TRUE) == 0) {
    stop("Variance of 'standardized_hri' is zero.")
  }
  
  # Project data
  if (st_is_longlat(h3_sf)) {
    h3_sf_proj <- st_transform(h3_sf, crs = 32632)
  } else {
    h3_sf_proj <- h3_sf
  }
  
  # Get centroids and coordinates
  centroids <- st_centroid(h3_sf_proj)
  coords <- st_coordinates(centroids)
  
  # Create neighborhood structure
  distance_threshold <- 25000  # 25 km
  nb <- dnearneigh(coords, 0, distance_threshold)
  
  # Create weights - using binary weights
  lw <- nb2listw(nb, style = "B", zero.policy = TRUE)
  
  # Calculate Moran's I
  moran_hri <- tryCatch({
    moran.test(h3_sf_proj$standardized_hri, lw, zero.policy = TRUE)
  }, error = function(e) {
    warning("Error in Moran's I test: ", e$message)
    return(NULL)
  })
  
  # Calculate LISA
  lisa_hri <- tryCatch({
    localmoran(h3_sf_proj$standardized_hri, lw, zero.policy = TRUE)
  }, error = function(e) {
    warning("Error in LISA analysis: ", e$message)
    return(NULL)
  })
  
  # Store and export results
  if (!is.null(moran_hri)) {
    moran_export <- data.frame(
      statistic = c("Moran's I", "p-value"),
      value = c(moran_hri$estimate["Moran I statistic"], moran_hri$p.value)
    )
    write.csv(moran_export, paste0("spatial_stats_", scenario_name, ".csv"), row.names = FALSE)
  }
  
  if (!is.null(lisa_hri)) {
    lisa_df <- as.data.frame(lisa_hri)
    write.csv(lisa_df, paste0("lisa_results_", scenario_name, ".csv"), row.names = FALSE)
  }
  
  return(list(global_moran = moran_hri, lisa_results = lisa_hri))
}
# ------------------------------
# 7B. CLUSTER-HRI CORRELATION
# ------------------------------
cluster_correlation <- function(h3_sf, scenario_name) {
  # Eta-squared for categorical-continuous relationship
  eta_squared <- function(cat_var, cont_var) {
    model <- aov(cont_var ~ cat_var)
    ef <- summary(model)[[1]]
    ef$"Sum Sq"[1] / sum(ef$"Sum Sq")
  }
  
  # Cramer's V calculation
  cramer_v_calc <- function(contingency_table) {
    chisq_result <- chisq.test(contingency_table, simulate.p.value = TRUE, B = 1000)
    n <- sum(contingency_table)
    min_dim <- min(dim(contingency_table)) - 1
    sqrt(chisq_result$statistic / (n * min_dim))
  }
  
  # Calculate metrics
  eta_sq_val <- eta_squared(h3_sf$dominant_cluster, h3_sf$standardized_hri)
  cramer_v_val <- cramer_v_calc(table(h3_sf$dominant_cluster, h3_sf$hri_category))
  
  # Combine into a data frame
  correlation_metrics <- data.frame(
    metric = c("Eta Squared", "Cramer's V"),
    value = c(eta_sq_val, cramer_v_val)
  )
  
  # Export results
  write.csv(
    correlation_metrics,
    paste0("correlation_metrics_", scenario_name, ".csv"),
    row.names = FALSE
  )
  
  return(correlation_metrics)
}

# ------------------------------
# 7C. PATTERN ANALYSIS
# ------------------------------
pattern_analysis <- function(h3_sf, scenario_name) {
  library(spatstat)
  library(spdep)
  library(dplyr)  # Ensure dplyr is loaded for the pipe operator
  
  # Project to metric CRS for all spatial analyses
  if (st_is_longlat(h3_sf)) {
    h3_sf_proj <- st_transform(h3_sf, crs = 32632)  # UTM zone 32N
  } else {
    h3_sf_proj <- h3_sf
  }

  # Hot Spot Analysis (Getis-Ord Gi*) — uses projected coordinates with 25 km threshold
  hotspot_coords <- st_coordinates(st_centroid(h3_sf_proj))
  weights <- nb2listw(dnearneigh(hotspot_coords, 0, 25000), zero.policy = TRUE)
  hotspots <- localG(h3_sf_proj$standardized_hri, weights)

  # Ripley's K function for cluster pattern
  
  # Validate geometries
  valid <- st_is_valid(h3_sf_proj)
  
  if (any(!valid)) {
    # Retrieve reasons for invalid geometries
    invalid_geometries <- st_is_valid(h3_sf_proj, reason = TRUE)
    
    cat("Invalid geometries detected:\n")
    print(invalid_geometries[!valid])
    
    # Attempt to fix invalid geometries
    h3_sf_proj <- st_make_valid(h3_sf_proj)
    
    # Re-check validity after attempting to fix
    valid_after_fix <- st_is_valid(h3_sf_proj)
    
    if (any(!valid_after_fix)) {
      warning("Some geometries remain invalid after attempting to fix.")
      print(st_is_valid(h3_sf_proj, reason = TRUE)[!valid_after_fix])
    } else {
      cat("All geometries are now valid after fixing.\n")
    }
  } else {
    cat("All geometries are valid.\n")
  }
  
  # Simplify geometries to ensure attribute consistency
  h3_sf_proj <- h3_sf_proj %>%
    group_by(h3_index) %>%
    summarize(geometry = st_union(geometry), .groups = "drop")
  
  # Convert to spatstat's 'ppp' object
  coords <- st_coordinates(st_centroid(st_geometry(h3_sf_proj)))
  window_bbox <- st_bbox(h3_sf_proj)
  
  # Define the observation window
  window <- as.owin(c(xrange = c(window_bbox["xmin"], window_bbox["xmax"]),
                      yrange = c(window_bbox["ymin"], window_bbox["ymax"])))
  
  # Plot the observation window and points for visual inspection
  plot(window, main = paste("Observation Window and Points -", scenario_name))
  points(coords, pch = 20, col = "red")
  
  # Check the number of points
  num_points <- nrow(coords)
  cat("Number of points in point pattern:", num_points, "\n")
  
  if (num_points < 20) {  # Adjust threshold as needed
    warning("Insufficient number of points for Ripley's K function. Requires at least 30 points.")
    # Proceed with hotspot analysis but skip Ripley's K
    k_export <- data.frame(
      r = NA,
      theo = NA,
      observed = NA
    )
  } else {
    # Create point pattern
    point_pattern <- ppp(coords[,1], coords[,2], window = window)
    
    # Validate the 'ppp' object
    if (point_pattern$n == 0) {
      warning("Point pattern contains zero points. Skipping Ripley's K function.")
      k_export <- data.frame(
        r = NA,
        theo = NA,
        observed = NA
      )
    } else {
      # Calculate Ripley's K function with error handling
      max_distance <- 1000  # Adjust based on data scale
      k_function <- tryCatch({
        Kest(point_pattern, r = seq(0, max_distance, by = 100), correction = "Ripley")
      }, error = function(e) {
        warning("Error in calculating Ripley's K function: ", e$message)
        return(NULL)
      })
      
      if (is.null(k_function)) {
        k_export <- data.frame(
          r = NA,
          theo = NA,
          observed = NA
        )
      } else {
        # Ensure 'theo' and 'obs' have the same length as 'r'
        if (length(k_function$theo) == length(k_function$r) &&
            length(k_function$obs) == length(k_function$r)) {
          k_export <- data.frame(
            r = k_function$r,
            theo = k_function$theo,
            observed = k_function$obs
          )
        } else {
          warning("Ripley's K function returned mismatched lengths. Assigning NA.")
          k_export <- data.frame(
            r = k_function$r,
            theo = NA,
            observed = NA
          )
        }
      }
    }
  }
  
  # Store results
  pattern_stats <- list(
    hotspots = hotspots,
    spatial_pattern = if (exists("k_function")) k_function else NULL
  )
  
  # Export Hotspot Analysis
  hotspot_df <- data.frame(
    location_id = 1:length(hotspots),
    gi_statistic = as.numeric(hotspots)
  )
  
  write.csv(
    hotspot_df,
    paste0("hotspot_analysis_", scenario_name, ".csv"),
    row.names = FALSE
  )
  
  # Export Ripley's K
  write.csv(
    k_export,
    paste0("ripley_k_analysis_", scenario_name, ".csv"),
    row.names = FALSE
  )
  
  return(pattern_stats)
}

# ------------------------------
# 7D. KDE-HRI OVERLAY
# ------------------------------
create_kde_hri_overlay <- function(h3_sf, boundary_data, scenario_name) {
  library(spatstat)
  library(ggplot2)
  
  # 1. Extract coordinates from h3_sf
  coords <- st_coordinates(st_centroid(h3_sf))
  
  # 2. Define the observation window based on boundary_data
  window_bbox <- st_bbox(boundary_data)
  window <- as.owin(c(xrange = c(window_bbox["xmin"], window_bbox["xmax"]),
                      yrange = c(window_bbox["ymin"], window_bbox["ymax"])))
  
  # 3. Create a point pattern
  point_pattern <- ppp(coords[,1], coords[,2], window = window)
  
  # 4. Calculate Kernel Density Estimate
  kde <- density.ppp(point_pattern, sigma = bw.scott(point_pattern))
  
  # 5. Convert KDE to dataframe for ggplot
  kde_df <- as.data.frame(kde)
  names(kde_df) <- c("x", "y", "density")
  
  # 6. Build H3 polygons with HRI
  unique_cells <- unique(h3_sf$h3_index)
  h3_polygons <- do.call(rbind, lapply(unique_cells, function(idx) {
    poly <- cell_to_polygon(idx, simple=FALSE)
    poly$h3_index <- idx
    return(poly)
  }))
  
  h3_hri_map <- st_as_sf(h3_polygons) %>%
    left_join(st_drop_geometry(h3_sf), by = "h3_index")
  
  # 7. Plot overlay
  p_overlay <- ggplot() +
    # Base KDE layer
    geom_raster(data = kde_df, aes(x = x, y = y, alpha = density), fill = "blue") +
    scale_alpha_continuous(name = "Point Density") +
    
    # HRI hexagons with transparency
    geom_sf(
      data = h3_hri_map,
      aes(fill = standardized_hri),
      alpha = 0.7
    ) +
    scale_fill_viridis_c(
      option = "magma",
      name = "Standardized HRI",
      labels = scales::percent
    ) +
    
    # Boundaries and formatting
    geom_sf(data = boundary_data, fill = NA, color = "white", size = 0.3) +
    labs(
      title = paste(scenario_name, "- HRI Distribution with Point Density Overlay"),
      subtitle = "Kernel Density Estimation + H3 Risk Index"
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(face="bold"),
      plot.subtitle = element_text(face="italic"),
      legend.position="right"
    )
  
  # Save plot
  ggsave(
    paste0("kde_hri_overlay_", scenario_name, ".png"),
    plot = p_overlay,
    width = 12,
    height = 15,
    dpi = 300
  )
  
  return(p_overlay)
}

# ------------------------------
# 8. ADVANCED VISUALIZATION FUNCTIONS
# ------------------------------
# Adjusted advanced visualization function with geometry handling
create_advanced_visualizations <- function(
    h3_data,           # synergy_h3_df from synergy_summaries
    boundary_data, 
    scenario_name,
    output_dir = "visualizations"
) {
  # Ensure output directory exists
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ------------------------------
  # 1. Convert synergy data frame -> sf object by building H3 polygons
  # ------------------------------
  unique_cells <- unique(h3_data$h3_index)
  h3_polygons <- do.call(rbind, lapply(unique_cells, function(idx) {
    poly <- cell_to_polygon(idx, simple=FALSE)
    poly$h3_index <- idx
    return(poly)
  }))
  
  # Join synergy columns (dominant_cluster, standardized_hri, etc.) to geometry
  h3_sf <- st_as_sf(h3_polygons) %>%
    left_join(st_drop_geometry(h3_data), by = "h3_index")
  
  # ------------------------------
  # 2. FACETED MAP: DOMINANT CLUSTER BY HRI CATEGORY
  # ------------------------------
  p_cluster_by_hri <- ggplot() +
    geom_sf(data = boundary_data, fill = NA, color = "grey60", size = 0.2) +
    geom_sf(
      data = h3_sf,
      aes(fill = dominant_cluster),
      color = "grey70",
      size = 0.1
    ) +
    scale_fill_brewer(palette = "Set1", name = "Dominant\nCluster") +
    facet_wrap(~hri_category, ncol = 2) +
    labs(
      title = paste(scenario_name, "- Dominant Clusters by HRI Category"),
      subtitle = "Spatial Distribution Analysis"
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(face="bold", size=12),
      plot.subtitle = element_text(face="italic", size=10),
      strip.text    = element_text(face="bold"),
      legend.position="right"
    )
  
  # Save the plot
  ggsave(
    file.path(output_dir, paste0(scenario_name, "_clusters_by_hri.png")),
    p_cluster_by_hri,
    width  = 12,
    height = 15,
    dpi    = 300
  )
  
  # ------------------------------
  # 3. FACETED MAP: HRI BY DOMINANT CLUSTER
  # ------------------------------
  p_hri_by_cluster <- ggplot() +
    geom_sf(data = boundary_data, fill=NA, color="grey60", size=0.2) +
    geom_sf(
      data = h3_sf,
      aes(fill = standardized_hri),
      color="grey70",
      size=0.1
    ) +
    scale_fill_viridis_c(
      option = "magma",
      name   = "Standardized\nHRI",
      labels = scales::percent
    ) +
    facet_wrap(~dominant_cluster, ncol=2) +
    labs(
      title    = paste(scenario_name, "- HRI Distribution by Dominant Cluster"),
      subtitle = "Spatial Pattern Analysis"
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(face="bold", size=12),
      plot.subtitle = element_text(face="italic", size=10),
      strip.text    = element_text(face="bold"),
      legend.position="right"
    )
  
  # Save the plot
  ggsave(
    file.path(output_dir, paste0(scenario_name, "_hri_by_clusters.png")),
    p_hri_by_cluster,
    width  = 12,
    height = 15,
    dpi    = 300
  )
  
  # ------------------------------
  # 4. BIVARIATE VISUALIZATION
  # ------------------------------
  # Create bivariate categories in h3_sf
  h3_sf <- h3_sf %>%
    mutate(
      hri_tercile = ntile(standardized_hri, 3),
      cluster_numeric = as.numeric(factor(dominant_cluster)),
      cluster_tercile = ntile(cluster_numeric, 3),
      bivariate_class = paste(hri_tercile, cluster_tercile, sep = "-")
    )
  
  # Bivariate color palette
  bivariate_colors <- c(
    "1-1"="#edf8fb", "1-2"="#bdc9e1", "1-3"="#67a9cf",
    "2-1"="#fdbf6f", "2-2"="#b2abd2", "2-3"="#6e7bc1",
    "3-1"="#fc8d59", "3-2"="#d95f0e", "3-3"="#8c510a"
  )
  
  p_bivariate <- ggplot() +
    geom_sf(data=boundary_data, fill=NA, color="grey60", size=0.2) +
    geom_sf(
      data=h3_sf,
      aes(fill=bivariate_class),
      color="grey70",
      size=0.1
    ) +
    scale_fill_manual(
      values = bivariate_colors,
      name   = "HRI - Cluster\nRelationship",
      guide  = guide_legend(
        title.position="top",
        nrow=3,
        byrow=TRUE
      )
    ) +
    labs(
      title = paste(scenario_name, "- Bivariate Analysis of HRI and Clusters"),
      subtitle = "Spatial Relationship Patterns"
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(face="bold", size=12),
      plot.subtitle = element_text(face="italic", size=10),
      legend.position="right"
    )
  
  # Save the plot
  ggsave(
    file.path(output_dir, paste0(scenario_name, "_bivariate.png")),
    p_bivariate,
    width  = 10,
    height = 12,
    dpi    = 300
  )
  
  # ------------------------------
  # 5. SIDE-BY-SIDE COMPARISON
  # ------------------------------
  # (A) Map: Dominant Clusters
  p_cluster <- ggplot() +
    geom_sf(data=boundary_data, fill=NA, color="grey60", size=0.2) +
    geom_sf(
      data=h3_sf,
      aes(fill=dominant_cluster),
      color="grey70",
      size=0.1
    ) +
    scale_fill_brewer(palette="Set1", name="Dominant\nCluster") +
    labs(title="Dominant Clusters") +
    theme_minimal()
  
  # (B) Map: Standardized HRI
  p_hri <- ggplot() +
    geom_sf(data=boundary_data, fill=NA, color="grey60", size=0.2) +
    geom_sf(
      data=h3_sf,
      aes(fill=standardized_hri),
      color="grey70",
      size=0.1
    ) +
    scale_fill_viridis_c(
      option="magma",
      name="Standardized\nHRI",
      labels=scales::percent
    ) +
    labs(title="HRI Distribution") +
    theme_minimal()
  
  # Combine side-by-side
  p_combined <- gridExtra::grid.arrange(
    p_cluster, p_hri,
    ncol=2,
    top=grid::textGrob(
      paste(scenario_name, "- Spatial Pattern Comparison"),
      gp=grid::gpar(fontsize=14, fontface="bold")
    )
  )
  
  # Save the combined plot
  ggsave(
    file.path(output_dir, paste0(scenario_name, "_side_by_side.png")),
    p_combined,
    width = 15,
    height=8,
    dpi=300
  )
  
  # ------------------------------
  # 6. RETURN THE PLOTS AND SF OBJECT
  # ------------------------------
  return(list(
    cluster_by_hri = p_cluster_by_hri,
    hri_by_cluster = p_hri_by_cluster,
    bivariate      = p_bivariate,
    side_by_side   = p_combined,
    final_sf       = h3_sf  # Return the sf object for advanced stats
  ))
}

# ------------------------------
# 9. MAIN ANALYSIS EXECUTION
# ------------------------------
# Define clustering scenarios
scenarios <- list(
  Amenity = list(
    features = c("encoded_Access", "encoded_EcoCertifd", "encoded_Green",
                 "encoded_Pool", "encoded_Sustainabl", "NumberfRms", 
                 "Ratings", "Price"),
    k = 4,
    m = 1.5
  ),
  ShortTerm = list(
    features = c("encoded_Type_2", "length_sea", "msft_heigh", 
                 "NumberfRms", "Price", "Ratings"),
    k = 3,
    m = 1.5
  ),
  
  Structural = list(
    features = c("msft_heigh", "area", "length_sea", "NumberfRms", 
                 "encoded_stock_cat_1"),
    k = 3,
    m = 1.5
  )
)

# Initialize results storage
clustering_results <- list()
synergy_summaries <- list()

# Execute main analysis pipeline
cat("\n=== Starting Main Analysis Pipeline ===\n")

# 6A. Execute Fuzzy Clustering for Each Scenario
for (scenario in names(scenarios)) {
  cat(paste0("\n--- Processing Scenario: ", scenario, " ---\n"))
  
  # Extract scenario parameters
  features <- scenarios[[scenario]]$features
  k_val <- scenarios[[scenario]]$k
  m_val <- scenarios[[scenario]]$m
  
  # Perform fuzzy clustering
  fuzzy_res <- fuzzy_cluster_analysis_features(
    sf_points = sf_data,
    extra_features = features,
    k = k_val,
    m = m_val
  )
  
  # Export validity metrics
  export_validities(fuzzy_res, paste0(scenario, "_validities.csv"))
  
  # Perform posthoc analysis
  posthoc_res <- posthoc_fuzzy_cluster_analysis(
    fuzzy_result_list = fuzzy_res,
    sf_points = sf_data,
    extra_features = features
  )
  
  # Export feature summaries
  export_feature_summaries(posthoc_res, paste0(scenario, "_feature_summaries.txt"))
  
  # Store results
  clustering_results[[scenario]] <- posthoc_res
  
  # 6B. Random Forest Interpretability Analysis
  cat("\n--- Performing Random Forest Analysis ---\n")
  
  df_for_rf <- posthoc_res$labeled_data %>%
    st_drop_geometry() %>%
    select(cluster_label, all_of(features)) %>%
    mutate(cluster_label = factor(cluster_label))
  
  # Handle any remaining NA values
  df_for_rf <- df_for_rf[complete.cases(df_for_rf), ]
  
  # Train/test split
  set.seed(42)
  train_indices <- sample(nrow(df_for_rf), size = 0.7 * nrow(df_for_rf))
  df_train <- df_for_rf[train_indices, ]
  df_test  <- df_for_rf[-train_indices, ]
  
  # Train random forest
  rf_model <- randomForest(
    cluster_label ~ .,
    data = df_train,
    ntree = 500,
    mtry = ceiling(sqrt(length(features))),
    importance = TRUE
  )
  
  # Evaluate performance
  pred_test <- predict(rf_model, newdata = df_test)
  conf_mat <- table(Predicted = pred_test, Actual = df_test$cluster_label)
  
  # Calculate F1 scores
  f1_scores <- compute_f1_scores(conf_mat)
  
  # Export F1 scores
  write.csv(f1_scores, paste0(scenario, "_rf_f1scores.csv"), row.names = FALSE)
  
  # Variable importance plots
  png(paste0(scenario, "_rf_varimp.png"), width = 800, height = 600)
  varImpPlot(rf_model, main = paste("Variable Importance -", scenario))
  dev.off()
  
  # 6C. H3 Integration and Spatial Analysis
  cat("\n--- Performing H3 Integration and Spatial Analysis ---\n")
  
  # Ensure H3 index exists
  if (!"h3_index" %in% names(posthoc_res$labeled_data)) {
    posthoc_res$labeled_data$h3_index <- point_to_cell(
      posthoc_res$labeled_data,
      res = 7
    )
  }
  
  # Calculate cluster proportions per H3 cell
  h3_cluster_props <- posthoc_res$labeled_data %>%
    st_drop_geometry() %>%
    group_by(h3_index, cluster_label) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(h3_index) %>%
    mutate(
      total = sum(count),
      proportion = count / total
    ) %>%
    ungroup()
  
  # Generate H3 polygons
  unique_cells <- unique(h3_cluster_props$h3_index)
  h3_polygons <- do.call(rbind, lapply(unique_cells, function(cell) {
    poly <- cell_to_polygon(cell, simple = FALSE)
    poly$h3_index <- cell
    return(poly)
  }))
  
  # Create dominant cluster map
  dom_clusters <- h3_cluster_props %>%
    group_by(h3_index) %>%
    slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  h3_dom_map <- st_as_sf(h3_polygons) %>%
    left_join(dom_clusters, by = "h3_index")
  
  # Plot dominant cluster map
  p_dom <- ggplot() +
    geom_sf(data = ita_boundaries_cropped, fill = NA, color = "grey60") +
    geom_sf(
      data = h3_dom_map,
      aes(fill = cluster_label),
      color = "grey70",
      size = 0.1
    ) +
    scale_fill_brewer(palette = "Set1", name = "Dominant\nCluster") +
    labs(
      title = paste(scenario, "- Dominant Clusters by H3 Cell"),
      subtitle = "H3 Resolution 7"
    ) +
    theme_minimal()
  
  ggsave(
    paste0(scenario, "_h3_dominant.png"),
    plot = p_dom,
    width = 8,
    height = 10,
    dpi = 300
  )
}

# 6D. Calculate HRI
cat("\n--- Calculating IPCC-Aligned HRI ---\n")
hri_results <- calculate_hri(sf_data)
sf_data_updated <- hri_results$sf_data
hri_df <- hri_results$hri_df

# Create HRI categories
hri_df <- create_hri_category(hri_df, n_breaks = 5)

# 6E. Synergy Analysis between Clustering and HRI
cat("\n--- Performing Cluster-HRI Synergy Analysis ---\n")

for (scenario in names(clustering_results)) {
  cat(paste0("\n--- Analyzing Synergy for Scenario: ", scenario, " ---\n"))
  
  # Prepare data
  cluster_data <- clustering_results[[scenario]]$labeled_data
  if (!"h3_index" %in% names(cluster_data)) {
    cluster_data$h3_index <- point_to_cell(cluster_data, res = 7)
  }
  
  # Merge with HRI
  synergy_df <- cluster_data %>%
    st_drop_geometry() %>%
    left_join(hri_df, by = "h3_index")
  
  # Call the synergy function
  synergy_result <- analyze_cluster_hri_synergy(
    synergy_df = synergy_df,
    hri_df = hri_df,
    scenario_name = scenario,
    filter_na = TRUE
  )
  
  # Store the entire synergy_result in synergy_summaries
  synergy_summaries[[scenario]] <- synergy_result
}

# ------------------------------
# 7. FINAL VISUALIZATION & REPORTING
# ------------------------------
# 7A. Create final HRI map
cat("\n--- Creating Final HRI Visualization ---\n")

# Generate H3 polygons for HRI
h3_cells <- unique(hri_df$h3_index)
h3_polygons_hri <- do.call(rbind, lapply(h3_cells, function(cell) {
  poly <- cell_to_polygon(cell, simple = FALSE)
  poly$h3_index <- cell
  return(poly)
}))

h3_hri_map <- st_as_sf(h3_polygons_hri) %>%
  left_join(hri_df, by = "h3_index")

# Create final HRI map
p_hri <- ggplot() +
  geom_sf(data = ita_boundaries_cropped, fill = NA, color = "grey60") +
  geom_sf(
    data = h3_hri_map,
    aes(fill = standardized_hri),
    color = "grey70",
    size = 0.1
  ) +
  scale_fill_viridis_c(
    option = "magma",
    name = "Standardized HRI",
    labels = scales::percent
  ) +
  labs(
    title = "IPCC-Aligned Hospitality Risk Index (HRI)",
    subtitle = "Standardized Values by H3 Cell"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(face="bold"),
    plot.subtitle = element_text(face="italic")
  )

# Save final map
ggsave(
  "final_hri_map.png",
  plot = p_hri,
  width = 10,
  height = 12,
  dpi = 300
)

# 7B. Export summary statistics
cat("\n--- Exporting Summary Statistics ---\n")

# HRI summary statistics
hri_summary <- hri_df %>%
  group_by(hri_category) %>%
  summarize(
    n_cells = n(),
    mean_hri = mean(standardized_hri, na.rm = TRUE),
    sd_hri = sd(standardized_hri, na.rm = TRUE),
    mean_uncertainty = mean(hri_cv, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(hri_summary, "hri_summary_statistics.csv", row.names = FALSE)

# Export synergy summaries properly
for (scenario in names(synergy_summaries)) {
  
  synergy_obj <- synergy_summaries[[scenario]]
  
  # 1) synergy_h3_df is a data frame (the merged H3 data)
  if ("synergy_h3_df" %in% names(synergy_obj)) {
    write.csv(
      synergy_obj$synergy_h3_df,
      paste0(scenario, "_synergy_h3_df.csv"),
      row.names = FALSE
    )
  }
  
  # 2) cluster_hri_summary is a data frame
  if ("cluster_hri_summary" %in% names(synergy_obj)) {
    write.csv(
      synergy_obj$cluster_hri_summary,
      paste0(scenario, "_cluster_hri_summary.csv"),
      row.names = FALSE
    )
  }
  
  # 3) cluster_uncertainty_summary is a data frame
  if ("cluster_uncertainty_summary" %in% names(synergy_obj)) {
    write.csv(
      synergy_obj$cluster_uncertainty_summary,
      paste0(scenario, "_cluster_uncertainty_summary.csv"),
      row.names = FALSE
    )
  }
  
  # 4) cluster_stats is a data frame
  if ("cluster_stats" %in% names(synergy_obj)) {
    write.csv(
      synergy_obj$cluster_stats,
      paste0(scenario, "_cluster_stats.csv"),
      row.names = FALSE
    )
  }
  
  # 5) cluster_size is a data frame
  if ("cluster_size" %in% names(synergy_obj)) {
    write.csv(
      synergy_obj$cluster_size,
      paste0(scenario, "_cluster_size.csv"),
      row.names = FALSE
    )
  }
  
  # 6) cluster_proportion_summary is a data frame
  if ("cluster_proportion_summary" %in% names(synergy_obj)) {
    write.csv(
      synergy_obj$cluster_proportion_summary,
      paste0(scenario, "_cluster_proportion_summary.csv"),
      row.names = FALSE
    )
  }
  
  # 7) cluster_uncertainty_prop_summary is a data frame
  if ("cluster_uncertainty_prop_summary" %in% names(synergy_obj)) {
    write.csv(
      synergy_obj$cluster_uncertainty_prop_summary,
      paste0(scenario, "_cluster_uncertainty_prop_summary.csv"),
      row.names = FALSE
    )
  }
  
  # 8) Export ANOVA results (captured in txt files)
  if ("anova_results" %in% names(synergy_obj)) {
    capture.output(
      synergy_obj$anova_results,
      file = paste0(scenario, "_anova_results.txt")
    )
  }
  
  # 9) Export Kruskal-Wallis results (captured in txt files)
  if ("kruskal_results" %in% names(synergy_obj)) {
    capture.output(
      synergy_obj$kruskal_results,
      file = paste0(scenario, "_kruskal_results.txt")
    )
  }
}

# ------------------------------
# 10. ADVANCED VISUALIZATION AND STATISTICAL ANALYSES
# ------------------------------
# Initialize results storage
analysis_results <- list()  # To store advanced stats and visualizations

for (scenario in names(synergy_summaries)) {
  cat(paste("\n--- Creating Advanced Visualizations for:", scenario, "---\n"))
  
  # synergy_summaries[[scenario]] is the full synergy result
  if (!"synergy_h3_df" %in% names(synergy_summaries[[scenario]])) {
    cat("No synergy_h3_df found for scenario:", scenario, "\n")
    next
  }
  
  # Extract synergy_h3_df
  synergy_df <- synergy_summaries[[scenario]]$synergy_h3_df
  
  # Create advanced visualizations
  advanced_viz <- create_advanced_visualizations(
    h3_data       = synergy_df,
    boundary_data = ita_boundaries_cropped,
    scenario_name = scenario
  )
  
  # Execute advanced statistical analyses
  cat("\n--- Performing Advanced Statistical Analyses ---\n")
  
  # 1. Spatial Autocorrelation
  spatial_stats <- spatial_analysis(
    h3_sf = advanced_viz$final_sf,  # Pass the sf object
    scenario_name = scenario
  )
  
  # 2. Correlation
  correlation_stats <- cluster_correlation(
    h3_sf = advanced_viz$final_sf,  # Pass the sf object
    scenario_name = scenario
  )
  
  # 3. Pattern Analysis
  pattern_stats <- pattern_analysis(
    h3_sf = advanced_viz$final_sf,  # Pass the sf object
    scenario_name = scenario
  )
  
  # 4. KDE-HRI Overlay
  kde_overlay <- create_kde_hri_overlay(
    h3_sf = advanced_viz$final_sf,      # Pass the sf object
    boundary_data = ita_boundaries_cropped,
    scenario_name = scenario
  )
  
  # Store all these in analysis_results
  analysis_results[[scenario]] <- list(
    advanced_viz     = advanced_viz,
    spatial_stats    = spatial_stats,
    correlation      = correlation_stats,
    pattern_stats    = pattern_stats,
    kde_overlay_plot = kde_overlay
  )
  
  cat(paste("\nAdvanced visualizations and statistical analyses created for scenario:", scenario, "\n"))
}

# ------------------------------
# 11. SENSITIVITY ANALYSIS: HRI WEIGHT PARAMETERS
# ------------------------------
cat("\n--- Starting Sensitivity Analysis on HRI Weight Parameters ---\n")

# Define a constrained grid of parameter combinations.
# Constraint 1: alpha + beta = 1 (multiplicative-additive balance)
# Constraint 2: a_H + a_E + a_V = 1 (additive sub-component balance)
alpha_vals <- c(0.3, 0.4, 0.5, 0.6, 0.7)

# Additive weight vectors that sum to 1.0
additive_weights <- list(
  c(a_H = 0.30, a_E = 0.35, a_V = 0.35),  # default
  c(a_H = 0.40, a_E = 0.30, a_V = 0.30),  # hazard-dominant
  c(a_H = 0.30, a_E = 0.40, a_V = 0.30),  # exposure-dominant
  c(a_H = 0.30, a_E = 0.30, a_V = 0.40),  # vulnerability-dominant
  c(a_H = 1/3,  a_E = 1/3,  a_V = 1/3)    # equal weights
)

param_grid <- do.call(rbind, lapply(alpha_vals, function(a) {
  do.call(rbind, lapply(additive_weights, function(w) {
    data.frame(alpha = a, beta = 1 - a, a_H = w["a_H"], a_E = w["a_E"], a_V = w["a_V"])
  }))
}))
rownames(param_grid) <- NULL
cat(sprintf("Sensitivity grid: %d combinations (alpha+beta=1, aH+aE+aV=1)\n", nrow(param_grid)))

# Initialize a list to store summary statistics for each combination.
sensitivity_results <- list()

# Loop through each combination in the grid
for(i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  # Convert the row to a string
  param_str <- paste(names(params), params, sep = "=", collapse = ", ")
  cat("Combination", i, ": ", param_str, "\n")
  
  # Recalculate HRI with the current set of parameters
  hri_res <- calculate_hri(sf_data, 
                           alpha = params$alpha, 
                           beta  = params$beta,
                           a_H   = params$a_H,
                           a_E   = params$a_E,
                           a_V   = params$a_V)
  hri_df_temp <- hri_res$hri_df
  
  # Compute summary statistics for the standardized HRI
  summary_stats <- hri_df_temp %>%
    summarize(
      mean_std_hri   = mean(standardized_hri, na.rm = TRUE),
      median_std_hri = median(standardized_hri, na.rm = TRUE),
      sd_std_hri     = sd(standardized_hri, na.rm = TRUE)
    )
  
  # Combine the current parameter settings with the summary statistics
  sensitivity_results[[i]] <- cbind(params, summary_stats)
}

# Combine all results into a data frame
sensitivity_df <- do.call(rbind, sensitivity_results)

# Save sensitivity summary results to CSV
write.csv(sensitivity_df, "sensitivity_analysis_summary.csv", row.names = FALSE)
cat("\n--- Sensitivity analysis summary saved to sensitivity_analysis_summary.csv ---\n")


# ------------------------------
# 11A. VISUAL DIAGNOSTICS: Density Plots
# ------------------------------
# Define a few specific parameter combinations to compare visually.
combos_to_plot <- list(
  default = list(alpha = 0.5, beta = 0.5, a_H = 0.3, a_E = 0.35, a_V = 0.35),
  alt1    = list(alpha = 0.6, beta = 0.4, a_H = 0.3, a_E = 0.35, a_V = 0.35),
  alt2    = list(alpha = 0.4, beta = 0.6, a_H = 0.3, a_E = 0.35, a_V = 0.35)
)

plot_list <- list()
for(name in names(combos_to_plot)) {
  combo <- combos_to_plot[[name]]
  hri_res <- calculate_hri(sf_data, 
                           alpha = combo$alpha, 
                           beta  = combo$beta,
                           a_H   = combo$a_H,
                           a_E   = combo$a_E,
                           a_V   = combo$a_V)
  hri_df_temp <- hri_res$hri_df %>%
    mutate(combo = name)
  plot_list[[name]] <- hri_df_temp
}

combined_hri_df <- do.call(rbind, plot_list)

# Plot density of standardized HRI for each parameter combination
p_sens <- ggplot(combined_hri_df, aes(x = standardized_hri, fill = combo)) +
  geom_density(alpha = 0.5) +
  labs(title = "Accommodations",
       x = "Standardized HRI", y = "Density",
       fill = "Parameter Combo") +
  theme_minimal()

ggsave("sensitivity_standardized_hri_density.png", plot = p_sens, width = 8, height = 6, dpi = 300)
cat("\n--- Sensitivity analysis density plot saved to sensitivity_standardized_hri_density.png ---\n")

# ------------------------------
# 11B. STATISTICAL TESTING: ANOVA
# ------------------------------
# For a statistical test, we combine all the density data and run an ANOVA on the standardized HRI
combined_hri_df$combo <- as.factor(combined_hri_df$combo)
anova_model <- aov(standardized_hri ~ combo, data = combined_hri_df)
anova_summary <- summary(anova_model)
capture.output(anova_summary, file = "sensitivity_anova_results.txt")
cat("\n--- Sensitivity analysis ANOVA results saved to sensitivity_anova_results.txt ---\n")


# ------------------------------
# 12. HRI DECOMPOSITION ANALYSIS
# ------------------------------
# Purpose: Quantify how much of the cluster-HRI association is driven by
# the vulnerability component vs. hazard × exposure. This addresses the
# circularity concern that vulnerability scoring mechanically inflates
# scores for certain establishment types.
cat("\n--- HRI Decomposition Analysis ---\n")

# Default HRI weights
decomp_alpha <- 0.5; decomp_beta <- 0.5
decomp_aH <- 0.3; decomp_aE <- 0.35; decomp_aV <- 0.35

decomp_all <- list()

for (scenario in names(synergy_summaries)) {
  cat(paste0("\n--- Decomposition for: ", scenario, " ---\n"))

  syn_df <- synergy_summaries[[scenario]]$synergy_h3_df

  # Compute decomposed components per H3 cell
  syn_df <- syn_df %>%
    mutate(
      he_product = hazard_avg * exposure_score,
      he_additive = decomp_aH * hazard_avg + decomp_aE * exposure_score,
      hri_no_vuln = decomp_alpha * he_product + decomp_beta * he_additive
    )

  # ANOVA on each component
  components <- c("hazard_avg", "exposure_score", "vulnerability_avg",
                   "he_product", "hri_no_vuln", "standardized_hri")
  comp_labels <- c("Hazard", "Exposure", "Vulnerability",
                    "HxE_product", "HRI_no_vuln", "Full_HRI_std")

  decomp_stats <- data.frame(
    scenario = character(), component = character(),
    F_value = numeric(), p_value = numeric(), eta_sq = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(components)) {
    f <- as.formula(paste(components[i], "~ factor(dominant_cluster)"))
    aov_res <- aov(f, data = syn_df)
    s <- summary(aov_res)
    eta2 <- s[[1]]$`Sum Sq`[1] / sum(s[[1]]$`Sum Sq`)

    decomp_stats <- rbind(decomp_stats, data.frame(
      scenario = scenario,
      component = comp_labels[i],
      F_value = round(s[[1]]$`F value`[1], 3),
      p_value = s[[1]]$`Pr(>F)`[1],
      eta_sq = round(eta2, 3)
    ))

    cat(sprintf("  %-20s F = %6.3f, p = %-10s, eta2 = %.3f\n",
                comp_labels[i], s[[1]]$`F value`[1],
                format(s[[1]]$`Pr(>F)`[1], digits = 3), eta2))
  }

  # Cluster means
  means_df <- syn_df %>%
    group_by(dominant_cluster) %>%
    summarize(
      n = n(),
      hazard = round(mean(hazard_avg, na.rm = TRUE), 3),
      exposure = round(mean(exposure_score, na.rm = TRUE), 3),
      vulnerability = round(mean(vulnerability_avg, na.rm = TRUE), 3),
      HxE = round(mean(he_product, na.rm = TRUE), 3),
      HRI_no_vuln = round(mean(hri_no_vuln, na.rm = TRUE), 3),
      HRI_full = round(mean(standardized_hri, na.rm = TRUE), 3),
      .groups = "drop"
    )

  # Save outputs
  write.csv(decomp_stats, paste0(scenario, "_decomposition_anova.csv"), row.names = FALSE)
  write.csv(means_df, paste0(scenario, "_decomposition.csv"), row.names = FALSE)

  decomp_all[[scenario]] <- list(stats = decomp_stats, means = means_df)
}

cat("\n--- Decomposition analysis complete ---\n")

# 13. EXTREME-CASE SENSITIVITY: MULTIPLICATIVE-ONLY vs ADDITIVE-ONLY
# -------------------------------------------------------------------
# Purpose: Test whether structural patterns hold under extreme HRI formulations
# (pure multiplicative α=1,β=0 vs pure additive α=0,β=1) rather than only
# under the balanced default (α=0.5, β=0.5).
cat("\n--- Extreme-Case Sensitivity Analysis ---\n")

extreme_configs <- list(
  multiplicative = c(alpha = 1.0, beta = 0.0),
  additive       = c(alpha = 0.0, beta = 1.0),
  default        = c(alpha = 0.5, beta = 0.5)
)

extreme_results <- list()

for (scenario in names(synergy_summaries)) {
  cat(paste0("\n=== ", scenario, " ===\n"))

  syn_df <- synergy_summaries[[scenario]]$synergy_h3_df

  for (config_name in names(extreme_configs)) {
    a <- extreme_configs[[config_name]]["alpha"]
    b <- extreme_configs[[config_name]]["beta"]

    # Recompute HRI under this configuration
    syn_df[[paste0("hri_", config_name)]] <-
      a * (syn_df$hazard_avg * syn_df$exposure_score * syn_df$vulnerability_avg) +
      b * (decomp_aH * syn_df$hazard_avg + decomp_aE * syn_df$exposure_score + decomp_aV * syn_df$vulnerability_avg)

    # Also compute without vulnerability
    syn_df[[paste0("hri_nv_", config_name)]] <-
      a * (syn_df$hazard_avg * syn_df$exposure_score) +
      b * (decomp_aH * syn_df$hazard_avg + decomp_aE * syn_df$exposure_score)

    # ANOVA: full
    aov_full <- aov(as.formula(paste0("hri_", config_name, " ~ factor(dominant_cluster)")), data = syn_df)
    s_full <- summary(aov_full)
    eta_full <- s_full[[1]]$`Sum Sq`[1] / sum(s_full[[1]]$`Sum Sq`)

    # ANOVA: no vuln
    aov_nv <- aov(as.formula(paste0("hri_nv_", config_name, " ~ factor(dominant_cluster)")), data = syn_df)
    s_nv <- summary(aov_nv)
    eta_nv <- s_nv[[1]]$`Sum Sq`[1] / sum(s_nv[[1]]$`Sum Sq`)

    extreme_results[[paste(scenario, config_name, sep = "_")]] <- data.frame(
      scenario = scenario,
      config = config_name,
      alpha = a, beta = b,
      full_F = round(s_full[[1]]$`F value`[1], 3),
      full_p = s_full[[1]]$`Pr(>F)`[1],
      full_eta2 = round(eta_full, 3),
      nv_F = round(s_nv[[1]]$`F value`[1], 3),
      nv_p = s_nv[[1]]$`Pr(>F)`[1],
      nv_eta2 = round(eta_nv, 3)
    )

    cat(sprintf("  %-15s Full: F=%6.3f, p=%-8s, eta2=%.3f | No-vuln: F=%6.3f, p=%-8s, eta2=%.3f\n",
                config_name,
                s_full[[1]]$`F value`[1], format(s_full[[1]]$`Pr(>F)`[1], digits=3), eta_full,
                s_nv[[1]]$`F value`[1], format(s_nv[[1]]$`Pr(>F)`[1], digits=3), eta_nv))
  }
}

extreme_df <- do.call(rbind, extreme_results)
rownames(extreme_df) <- NULL
write.csv(extreme_df, "extreme_sensitivity_anova.csv", row.names = FALSE)
cat("\nExtreme sensitivity results saved to extreme_sensitivity_anova.csv\n")
cat("\n--- Extreme-case sensitivity complete ---\n")

# ------------------------------
# 13. LISA & Gi* MAPS (landscape: LISA 3 panels | Gi* 3 panels, cartolight basemap)
# ------------------------------
cat("\n=== Creating LISA & Gi* Maps ===\n")
library(ggspatial)

classify_lisa_cat <- function(lisa_df, hri_values) {
  mean_hri <- mean(hri_values, na.rm = TRUE)
  p_col <- names(lisa_df)[grep("Pr", names(lisa_df))][1]
  category <- dplyr::case_when(
    lisa_df[[p_col]] >= 0.05 ~ "Not Significant",
    lisa_df$Ii > 0 & hri_values > mean_hri ~ "High-High",
    lisa_df$Ii > 0 & hri_values <= mean_hri ~ "Low-Low",
    lisa_df$Ii < 0 & hri_values > mean_hri ~ "High-Low",
    lisa_df$Ii < 0 & hri_values <= mean_hri ~ "Low-High",
    TRUE ~ "Not Significant"
  )
  factor(category, levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not Significant"))
}

classify_gi_cat <- function(gi_z) {
  category <- dplyr::case_when(
    gi_z > 2.576 ~ "Hot Spot (p<0.01)",
    gi_z > 1.96 ~ "Hot Spot (p<0.05)",
    gi_z < -2.576 ~ "Cold Spot (p<0.01)",
    gi_z < -1.96 ~ "Cold Spot (p<0.05)",
    TRUE ~ "Not Significant"
  )
  factor(category, levels = c(
    "Hot Spot (p<0.01)", "Hot Spot (p<0.05)", "Not Significant",
    "Cold Spot (p<0.05)", "Cold Spot (p<0.01)"
  ))
}

lisa_colors <- c("High-High" = "#d7191c", "Low-Low" = "#2c7bb6",
                 "High-Low" = "#fdae61", "Low-High" = "#abd9e9",
                 "Not Significant" = "#d9d9d9")
gi_colors <- c("Hot Spot (p<0.01)" = "#7b3294", "Hot Spot (p<0.05)" = "#c2a5cf",
               "Not Significant" = "#d9d9d9",
               "Cold Spot (p<0.05)" = "#a6dba0", "Cold Spot (p<0.01)" = "#008837")

theme_pub_map <- function() {
  theme_void() + theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    legend.margin = margin(t = 2, b = 2),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2))
}

build_panel <- function(h3_sf, fill_col, palette,
                        legend_title, title_label, add_scale = FALSE) {
  h3 <- st_set_crs(h3_sf, 4326)
  h3_sig <- h3[h3[[fill_col]] != "Not Significant", ]
  h3_ns  <- h3[h3[[fill_col]] == "Not Significant", ]
  p <- ggplot() +
    annotation_map_tile(type = "cartolight", zoom = 8, quiet = TRUE) +
    geom_sf(data = h3_ns, fill = "#bdbdbd", alpha = 0.4, color = NA, size = 0) +
    geom_sf(data = h3_sig, aes(fill = .data[[fill_col]]), alpha = 0.9, color = NA, size = 0) +
    scale_fill_manual(values = palette, name = legend_title, drop = FALSE) +
    labs(title = title_label) + coord_sf(crs = 4326) + theme_pub_map() +
    guides(fill = guide_legend(nrow = 1, override.aes = list(size = 0)))
  if (add_scale) p <- p + annotation_scale(location = "br", width_hint = 0.2,
    text_col = "black", line_col = "black", bar_cols = c("black", "white"), text_cex = 0.6)
  p
}

acc_scenarios_map <- c("Amenity", "ShortTerm", "Structural")
acc_titles <- c(Amenity = "Amenity", ShortTerm = "Market/Short-Term", Structural = "Structural")
acc_panels <- list()

for (sc in acc_scenarios_map) {
  synergy_df <- synergy_summaries[[sc]]$synergy_h3_df
  unique_cells <- unique(synergy_df$h3_index)
  h3_polys <- do.call(rbind, lapply(unique_cells, function(idx) {
    poly <- cell_to_polygon(idx, simple = FALSE); poly$h3_index <- idx; poly
  }))
  h3_sf <- st_as_sf(h3_polys) %>% left_join(st_drop_geometry(synergy_df), by = "h3_index")
  lisa_df <- read.csv(paste0("lisa_results_", sc, ".csv"))
  gi_df <- read.csv(paste0("hotspot_analysis_", sc, ".csv"))
  h3_sf$lisa_cat <- classify_lisa_cat(lisa_df, h3_sf$standardized_hri)
  h3_sf$gi_cat <- classify_gi_cat(gi_df$gi_statistic)
  cat("  ", sc, "— LISA sig:", sum(h3_sf$lisa_cat != "Not Significant"),
      "| Gi* hot:", sum(grepl("Hot", h3_sf$gi_cat)), "cold:", sum(grepl("Cold", h3_sf$gi_cat)), "\n")
  acc_panels[[paste0(sc, "_lisa")]] <- build_panel(
    h3_sf, "lisa_cat", lisa_colors, "LISA Cluster", acc_titles[sc])
  acc_panels[[paste0(sc, "_gi")]] <- build_panel(
    h3_sf, "gi_cat", gi_colors, "Gi* Hot/Cold Spot", acc_titles[sc],
    add_scale = (sc == "Structural"))
}

lisa_block <- (acc_panels$Amenity_lisa | acc_panels$ShortTerm_lisa | acc_panels$Structural_lisa) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
gi_block <- (acc_panels$Amenity_gi | acc_panels$ShortTerm_gi | acc_panels$Structural_gi) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
acc_fig <- wrap_elements(lisa_block) | wrap_elements(gi_block)

ggsave("lisa_gi_acc_new.png", plot = acc_fig, width = 20, height = 6, dpi = 300, bg = "white")
cat("\nSaved: lisa_gi_acc_new.png\n")

# ------------------------------
# END OF SCRIPT
# ------------------------------

save.image(here::here("Acc_new", "Results_Paper.RData"))

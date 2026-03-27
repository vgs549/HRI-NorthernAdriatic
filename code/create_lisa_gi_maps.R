# =============================================================================
# LISA & Gi* Maps for Accommodations and Services
# Layout: 1 row × 6 panels (LISA: Amenity, Market, Structural | Gi*: same)
# Scenario titles on top outside map, separate legends per method
# =============================================================================

library(h3jsr)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggspatial)

# Base directory
base_dir <- "C:/Users/vgs/Downloads/SpatialAnalysis"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

classify_lisa <- function(lisa_df, hri_values) {
  mean_hri <- mean(hri_values, na.rm = TRUE)
  p_col <- names(lisa_df)[grep("Pr", names(lisa_df))][1]
  category <- case_when(
    lisa_df[[p_col]] >= 0.05 ~ "Not Significant",
    lisa_df$Ii > 0 & hri_values > mean_hri ~ "High-High",
    lisa_df$Ii > 0 & hri_values <= mean_hri ~ "Low-Low",
    lisa_df$Ii < 0 & hri_values > mean_hri ~ "High-Low",
    lisa_df$Ii < 0 & hri_values <= mean_hri ~ "Low-High",
    TRUE ~ "Not Significant"
  )
  factor(category, levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not Significant"))
}

classify_gi <- function(gi_z) {
  category <- case_when(
    gi_z > 2.576 ~ "Hot Spot (p<0.01)",
    gi_z > 1.96 ~ "Hot Spot (p<0.05)",
    gi_z < -2.576 ~ "Cold Spot (p<0.01)",
    gi_z < -1.96 ~ "Cold Spot (p<0.05)",
    TRUE ~ "Not Significant"
  )
  factor(category, levels = c(
    "Hot Spot (p<0.01)", "Hot Spot (p<0.05)",
    "Not Significant",
    "Cold Spot (p<0.05)", "Cold Spot (p<0.01)"
  ))
}

# --- Color palettes ---
lisa_colors <- c(
  "High-High" = "#d7191c",
  "Low-Low" = "#2c7bb6",
  "High-Low" = "#fdae61",
  "Low-High" = "#abd9e9",
  "Not Significant" = "#d9d9d9"
)

gi_colors <- c(
  "Hot Spot (p<0.01)" = "#7b3294",
  "Hot Spot (p<0.05)" = "#c2a5cf",
  "Not Significant" = "#d9d9d9",
  "Cold Spot (p<0.05)" = "#a6dba0",
  "Cold Spot (p<0.01)" = "#008837"
)

# --- Publication theme ---
theme_pub_map <- function() {
  theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 8),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      legend.margin = margin(t = 2, b = 2),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    )
}

# --- Panel builder ---
build_panel <- function(h3_sf, fill_col, palette,
                        legend_title, title_label, add_scale = FALSE) {
  h3_sf_4326 <- st_set_crs(h3_sf, 4326)
  h3_sig <- h3_sf_4326[h3_sf_4326[[fill_col]] != "Not Significant", ]
  h3_ns  <- h3_sf_4326[h3_sf_4326[[fill_col]] == "Not Significant", ]

  p <- ggplot() +
    annotation_map_tile(type = "cartolight", zoom = 8, quiet = TRUE) +
    geom_sf(data = h3_ns, fill = "#bdbdbd", alpha = 0.4, color = NA, size = 0) +
    geom_sf(data = h3_sig, aes(fill = .data[[fill_col]]), alpha = 0.9, color = NA, size = 0) +
    scale_fill_manual(values = palette, name = legend_title, drop = FALSE) +
    labs(title = title_label) +
    coord_sf(crs = 4326) +
    theme_pub_map() +
    guides(fill = guide_legend(nrow = 1, override.aes = list(size = 0)))

  if (add_scale) {
    p <- p + annotation_scale(
      location = "br", width_hint = 0.2,
      text_col = "black", line_col = "black",
      bar_cols = c("black", "white"), text_cex = 0.6
    )
  }
  p
}

# =============================================================================
# ACCOMMODATIONS
# =============================================================================

cat("\n=== Loading Accommodation Data ===\n")
load(file.path(base_dir, "Acc_new", "Results_Paper.RData"))

acc_synergy <- synergy_summaries
acc_scenarios <- c("Amenity", "ShortTerm", "Structural")
acc_titles <- c(Amenity = "Amenity", ShortTerm = "Market/Short-Term", Structural = "Structural")
acc_panels <- list()

for (sc in acc_scenarios) {
  cat("  Processing:", sc, "\n")
  synergy_df <- acc_synergy[[sc]]$synergy_h3_df
  unique_cells <- unique(synergy_df$h3_index)
  h3_polys <- do.call(rbind, lapply(unique_cells, function(idx) {
    poly <- cell_to_polygon(idx, simple = FALSE)
    poly$h3_index <- idx
    return(poly)
  }))
  h3_sf <- st_as_sf(h3_polys) %>%
    left_join(st_drop_geometry(synergy_df), by = "h3_index")

  lisa_df <- read.csv(file.path(base_dir, "Acc_new", paste0("lisa_results_", sc, ".csv")))
  gi_df <- read.csv(file.path(base_dir, "Acc_new", paste0("hotspot_analysis_", sc, ".csv")))

  h3_sf$lisa_cat <- classify_lisa(lisa_df, h3_sf$standardized_hri)
  h3_sf$gi_cat <- classify_gi(gi_df$gi_statistic)

  cat("    LISA sig:", sum(h3_sf$lisa_cat != "Not Significant"),
      "| Gi* hot:", sum(grepl("Hot", h3_sf$gi_cat)),
      "cold:", sum(grepl("Cold", h3_sf$gi_cat)), "\n")

  acc_panels[[paste0(sc, "_lisa")]] <- build_panel(
    h3_sf, "lisa_cat", lisa_colors, "LISA Cluster", acc_titles[sc]
  )
  acc_panels[[paste0(sc, "_gi")]] <- build_panel(
    h3_sf, "gi_cat", gi_colors, "Gi* Hot/Cold Spot", acc_titles[sc],
    add_scale = (sc == "Structural")
  )
}

# Compose: LISA block (3 panels + legend) | Gi* block (3 panels + legend)
lisa_block <- (acc_panels$Amenity_lisa | acc_panels$ShortTerm_lisa | acc_panels$Structural_lisa) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
gi_block <- (acc_panels$Amenity_gi | acc_panels$ShortTerm_gi | acc_panels$Structural_gi) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

acc_fig <- wrap_elements(lisa_block) | wrap_elements(gi_block)

ggsave(
  file.path(base_dir, "Acc_new", "lisa_gi_acc_new.png"),
  plot = acc_fig, width = 20, height = 6, dpi = 300, bg = "white"
)
cat("\nSaved: Acc_new/lisa_gi_acc_new.png\n")

# =============================================================================
# SERVICES
# =============================================================================

cat("\n=== Loading Service Data ===\n")
load(file.path(base_dir, "Services_new", "Results_Paper.RData"))

ser_synergy <- synergy_summaries
ser_scenarios <- c("Amenity_Ser", "Market_Ser", "Structural_Ser")
ser_titles <- c(Amenity_Ser = "Amenity", Market_Ser = "Market", Structural_Ser = "Structural")
ser_panels <- list()

for (sc in ser_scenarios) {
  cat("  Processing:", sc, "\n")
  synergy_df <- ser_synergy[[sc]]$synergy_h3_df
  unique_cells <- unique(synergy_df$h3_index)
  h3_polys <- do.call(rbind, lapply(unique_cells, function(idx) {
    poly <- cell_to_polygon(idx, simple = FALSE)
    poly$h3_index <- idx
    return(poly)
  }))
  h3_sf <- st_as_sf(h3_polys) %>%
    left_join(st_drop_geometry(synergy_df), by = "h3_index")

  lisa_df <- read.csv(file.path(base_dir, "Services_new", paste0("lisa_results_", sc, ".csv")))
  gi_df <- read.csv(file.path(base_dir, "Services_new", paste0("hotspot_analysis_", sc, ".csv")))

  h3_sf$lisa_cat <- classify_lisa(lisa_df, h3_sf$standardized_hri)
  h3_sf$gi_cat <- classify_gi(gi_df$gi_statistic)

  cat("    LISA sig:", sum(h3_sf$lisa_cat != "Not Significant"),
      "| Gi* hot:", sum(grepl("Hot", h3_sf$gi_cat)),
      "cold:", sum(grepl("Cold", h3_sf$gi_cat)), "\n")

  ser_panels[[paste0(sc, "_lisa")]] <- build_panel(
    h3_sf, "lisa_cat", lisa_colors, "LISA Cluster", ser_titles[sc]
  )
  ser_panels[[paste0(sc, "_gi")]] <- build_panel(
    h3_sf, "gi_cat", gi_colors, "Gi* Hot/Cold Spot", ser_titles[sc],
    add_scale = (sc == "Structural_Ser")
  )
}

lisa_block_ser <- (ser_panels$Amenity_Ser_lisa | ser_panels$Market_Ser_lisa | ser_panels$Structural_Ser_lisa) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
gi_block_ser <- (ser_panels$Amenity_Ser_gi | ser_panels$Market_Ser_gi | ser_panels$Structural_Ser_gi) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ser_fig <- wrap_elements(lisa_block_ser) | wrap_elements(gi_block_ser)

ggsave(
  file.path(base_dir, "Services_new", "lisa_gi_ser_new.png"),
  plot = ser_fig, width = 20, height = 6, dpi = 300, bg = "white"
)
cat("\nSaved: Services_new/lisa_gi_ser_new.png\n")

cat("\n=== LISA & Gi* maps complete ===\n")

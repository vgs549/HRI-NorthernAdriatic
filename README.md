# Replication Package: Divergent Risk Patterns Across Coastal Hospitality Sectors in the Northern Adriatic

**Author:** Vilane Gonçalves Sales
**Institution:** Leibniz Centre for Tropical Marine Research (ZMT), Bremen, Germany
**Corresponding publication:** Sales, V. G. (2026). "Divergent risk patterns across coastal hospitality sectors in the Northern Adriatic." *Open Research Europe*.
**Dataset DOI:** https://doi.org/10.5281/zenodo.14929876
**License:** Code: MIT | Data: CC-BY 4.0

---

## Overview

This package contains analysis code, output results, and extended data for a spatial vulnerability assessment of 10,022 coastal hospitality establishments (5,030 accommodations, 4,992 services) across coastal Veneto and Emilia-Romagna in the Northern Adriatic, Italy.

The analysis develops a dual-component Hospitality Risk Index (HRI) integrating hazard (sea-level rise), exposure (elevation, coastal proximity, spatial density), and vulnerability (structural characteristics, economic positioning) following the IPCC (2014) risk framework. Fuzzy C-means clustering identifies establishment typologies across amenity, market, and structural dimensions. Exploratory Spatial Data Analysis (ESDA) characterizes spatial organization of risk.

---

## Repository Structure

```
zenodo_package/
├── README.md                          # This file
├── code/                              # Analysis scripts
│   ├── Paper_HRI_routine.R            # Accommodation HRI, clustering, ESDA
│   ├── Paper_HRI_routine_ser.R        # Service HRI, clustering, ESDA
│   ├── create_lisa_gi_maps.R          # LISA/Gi* spatial autocorrelation maps
│   └── hri_flood_validation.R         # External validation against Copernicus EMS
├── results/                           # Analysis outputs (CSV/TXT)
│   ├── accommodations/                # Accommodation-sector results
│   │   ├── hri_summary_statistics.csv
│   │   ├── *_validities.csv           # Clustering validation (PC, CE, XB, F1)
│   │   ├── *_cluster_stats.csv        # Cluster descriptive statistics
│   │   ├── *_decomposition_anova.csv  # Component decomposition ANOVA
│   │   ├── lisa_results_*.csv         # LISA local spatial autocorrelation
│   │   ├── hotspot_analysis_*.csv     # Getis-Ord Gi* hotspot/coldspot
│   │   ├── sensitivity_analysis_summary.csv
│   │   └── synergy_*_anova_results.txt # Cluster-HRI association tests
│   ├── services/                      # Service-sector results (same structure)
│   └── validation/                    # External flood validation
│       ├── validation_50m_summary.csv  # Mann-Whitney U test results
│       └── validation_50m_flood_overlay.csv # Establishment-level flood proximity
├── extended_data/                     # Extended data for publication
│   ├── figures/
│   │   ├── ed_fig1a_sensitivity_acc.png   # Sensitivity density (accommodations)
│   │   ├── ed_fig1b_sensitivity_ser.png   # Sensitivity density (services)
│   │   ├── ed_fig2_validation_boxplot_hri.png  # HRI by flood proximity
│   │   ├── ed_fig3_validation_components.png   # Components by flood proximity
│   │   └── ed_fig4_validation_map_overlay.png  # Flood validation map
│   └── tables/
│       ├── ed_table1_clustering_validation.csv # PC, CE, XB, F1 by scenario
│       └── ed_table2_lisa_gi_summary.csv       # LISA/Gi* statistics summary
```

---

## Software Requirements

### R (version 4.3+)
Required packages: `tidyverse`, `sf`, `h3jsr`, `e1071`, `spdep`, `here`, `ggplot2`, `viridis`, `patchwork`

Install via:
```r
install.packages(c("tidyverse", "sf", "h3jsr", "e1071", "spdep", "here", "ggplot2", "viridis", "patchwork"))
```

---

## Reproducing the Analysis

### Step 1: Obtain source data
Download the hospitality establishment dataset from Zenodo: https://doi.org/10.5281/zenodo.14929876

### Step 2: Run accommodation analysis
```r
source("code/Paper_HRI_routine.R")
```
This produces all accommodation results in `results/accommodations/`.

### Step 3: Run service analysis
```r
source("code/Paper_HRI_routine_ser.R")
```
This produces all service results in `results/services/`.

### Step 4: Run external validation
```r
source("code/hri_flood_validation.R")
```
This produces validation results in `results/validation/`.

### Step 5: Generate LISA/Gi* maps
```r
source("code/create_lisa_gi_maps.R")
```

---

## Extended Data Description

### Extended Data Table 1: Clustering Validation Metrics
Partition Coefficient (PC), Classification Entropy (CE), Xie-Beni Index (XB), and F1 scores for each clustering scenario (Amenity, Market/Short-Term, Structural) and sector (accommodations, services).

### Extended Data Table 2: Local Spatial Statistics Summary
LISA significant cell counts, cluster/outlier breakdown, and Getis-Ord Gi* hotspot/coldspot counts by scenario and sector.

### Extended Data Figure 1: Sensitivity Analysis
Density distributions of standardized HRI values across 25 constrained parameter combinations (alpha + beta = 1, aH + aE + aV = 1). Accommodation CV = 10.1%; services CV = 9.4%.

### Extended Data Figure 2: HRI by Flood Proximity
Box plots of establishment-level HRI for flood-proximate (within 50 m of Copernicus EMS flood extent) and non-proximate establishments. Mann-Whitney U: p < 0.001 for both sectors.

### Extended Data Figure 3: HRI Components by Flood Proximity
Component-level comparison (hazard, physical exposure, composite exposure, vulnerability) between flood-proximate and non-proximate establishments. Exposure drives the difference.

### Extended Data Figure 4: Flood Validation Map Overlay
Establishment locations within Copernicus EMS Areas of Interest (Venice EMSR409, Comacchio/Ravenna EMSR642) showing flood-proximate classification.

---

## Citation

If you use this code or data, please cite:

```
Sales, V. G. (2026). Divergent risk patterns across coastal hospitality sectors
in the Northern Adriatic. Open Research Europe. https://doi.org/10.5281/zenodo.19134970

Sales, V. G. (2024). Dataset of coastal hospitality establishments with building
characteristics for the Northern Adriatic, Italy. Data in Brief.
https://doi.org/10.5281/zenodo.14929876
```

---

## Contact

Vilane Gonçalves Sales
Leibniz Centre for Tropical Marine Research (ZMT)
Bremen, Germany

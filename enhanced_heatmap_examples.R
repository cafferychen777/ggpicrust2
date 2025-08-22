# Enhanced pathway_heatmap() Function Examples
# ============================================
# 
# This script demonstrates the new features added to the pathway_heatmap() function
# in ggpicrust2, including hierarchical clustering, faceted heatmaps, and 
# customizable color bars.

# Load required libraries
library(ggpicrust2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggh4x)

# Load sample data
data("metacyc_abundance")
data("metadata")

# Prepare the data
abundance_data <- metacyc_abundance %>% column_to_rownames("pathway")

# Perform differential abundance analysis
daa_results <- pathway_daa(
  abundance = abundance_data,
  metadata = metadata,
  group = "Environment",
  daa_method = "LinDA"
)

# Annotate the results
annotated_results <- pathway_annotation(
  pathway = "MetaCyc",
  daa_results_df = daa_results,
  ko_to_kegg = FALSE
)

# Filter significant features
significant_features <- daa_results %>% filter(p_adjust < 0.05)

# Prepare abundance matrix for visualization
vis_abundance <- metacyc_abundance %>%
  right_join(
    annotated_results %>% select(all_of(c("feature", "description"))),
    by = c("pathway" = "feature")
  ) %>%
  filter(pathway %in% significant_features$feature) %>%
  select(-"pathway") %>%
  column_to_rownames("description")

# ============================================
# EXAMPLE 1: Basic Heatmap (Original Functionality)
# ============================================

cat("Creating Example 1: Basic Heatmap\n")
p1 <- pathway_heatmap(
  abundance = vis_abundance,
  metadata = metadata,
  group = "Environment"
)

# Save the plot
ggsave("example1_basic_heatmap.pdf", p1, width = 10, height = 8)

# ============================================
# EXAMPLE 2: Hierarchical Clustering - Rows Only
# ============================================

cat("Creating Example 2: Row Clustering\n")
p2 <- pathway_heatmap(
  abundance = vis_abundance,
  metadata = metadata,
  group = "Environment",
  cluster_rows = TRUE,
  clustering_method = "ward.D2",
  clustering_distance = "correlation",
  dendro_line_size = 0.8
)

ggsave("example2_row_clustering.pdf", p2, width = 12, height = 8)

# ============================================
# EXAMPLE 3: Hierarchical Clustering - Columns Only
# ============================================

cat("Creating Example 3: Column Clustering\n")
p3 <- pathway_heatmap(
  abundance = vis_abundance,
  metadata = metadata,
  group = "Environment",
  cluster_cols = TRUE,
  clustering_method = "complete",
  clustering_distance = "euclidean"
)

ggsave("example3_column_clustering.pdf", p3, width = 10, height = 10)

# ============================================
# EXAMPLE 4: Both Row and Column Clustering
# ============================================

cat("Creating Example 4: Both Row and Column Clustering\n")
p4 <- pathway_heatmap(
  abundance = vis_abundance,
  metadata = metadata,
  group = "Environment",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "average",
  clustering_distance = "manhattan",
  dendro_line_size = 1.0
)

ggsave("example4_both_clustering.pdf", p4, width = 12, height = 10)

# ============================================
# EXAMPLE 5: Custom Color Schemes
# ============================================

cat("Creating Example 5: Custom Color Schemes\n")
p5 <- pathway_heatmap(
  abundance = vis_abundance,
  metadata = metadata,
  group = "Environment",
  cluster_rows = TRUE,
  clustering_method = "ward.D2",
  low_color = "#053061",     # Dark blue
  mid_color = "#f7f7f7",     # Light gray
  high_color = "#67001f",    # Dark red
  colors = c("#E8F4FD", "#FFE6E6")  # Custom facet colors
)

ggsave("example5_custom_colors.pdf", p5, width = 12, height = 8)

# ============================================
# EXAMPLE 6: Custom Color Bar Settings
# ============================================

cat("Creating Example 6: Custom Color Bar\n")
p6 <- pathway_heatmap(
  abundance = vis_abundance,
  metadata = metadata,
  group = "Environment",
  colorbar_title = "Standardized Abundance",
  colorbar_position = "bottom",
  colorbar_width = 8,
  colorbar_height = 0.8,
  colorbar_breaks = c(-2, -1, 0, 1, 2)
)

ggsave("example6_custom_colorbar.pdf", p6, width = 10, height = 9)

# ============================================
# EXAMPLE 7: Left-positioned Color Bar
# ============================================

cat("Creating Example 7: Left Color Bar\n")
p7 <- pathway_heatmap(
  abundance = vis_abundance,
  metadata = metadata,
  group = "Environment",
  cluster_rows = TRUE,
  colorbar_position = "left",
  colorbar_title = "Z-Score",
  low_color = "#2166ac",
  high_color = "#b2182b"
)

ggsave("example7_left_colorbar.pdf", p7, width = 12, height = 8)

# ============================================
# EXAMPLE 8: Faceted Heatmap (if additional grouping variable exists)
# ============================================

# Create additional grouping variable for demonstration
metadata_extended <- metadata %>%
  mutate(Batch = factor(rep(c("Batch1", "Batch2"), length.out = nrow(.))))

cat("Creating Example 8: Faceted Heatmap\n")
p8 <- pathway_heatmap(
  abundance = vis_abundance,
  metadata = metadata_extended,
  group = "Environment",
  facet_by = "Batch",
  colors = c("#E8F4FD", "#FFE6E6", "#E8F8E8", "#FFF8E8")
)

ggsave("example8_faceted_heatmap.pdf", p8, width = 14, height = 8)

# ============================================
# EXAMPLE 9: Advanced Combination
# ============================================

cat("Creating Example 9: Advanced Combination\n")
p9 <- pathway_heatmap(
  abundance = vis_abundance,
  metadata = metadata,
  group = "Environment",
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # Keep column order by group
  clustering_method = "ward.D2",
  clustering_distance = "correlation",
  dendro_line_size = 0.6,
  low_color = "#313695",
  mid_color = "#FFFFBF", 
  high_color = "#A50026",
  colorbar_title = "Pathway Abundance (Z-score)",
  colorbar_position = "right",
  font_size = 10,
  show_row_names = TRUE
)

ggsave("example9_advanced_combination.pdf", p9, width = 14, height = 10)

# ============================================
# EXAMPLE 10: Different Clustering Methods Comparison
# ============================================

clustering_methods <- c("complete", "average", "ward.D2")
distances <- c("euclidean", "correlation", "manhattan")

cat("Creating Example 10: Clustering Methods Comparison\n")

for (i in seq_along(clustering_methods)) {
  method <- clustering_methods[i]
  distance <- distances[i]
  
  p <- pathway_heatmap(
    abundance = vis_abundance,
    metadata = metadata,
    group = "Environment",
    cluster_rows = TRUE,
    clustering_method = method,
    clustering_distance = distance,
    colorbar_title = paste("Method:", method, "\nDistance:", distance)
  )
  
  ggsave(paste0("example10_clustering_", method, "_", distance, ".pdf"), 
         p, width = 12, height = 8)
}

# ============================================
# Summary
# ============================================

cat("\n=== Enhanced pathway_heatmap() Examples Complete ===\n")
cat("Generated files:\n")
cat("- example1_basic_heatmap.pdf: Basic heatmap functionality\n")
cat("- example2_row_clustering.pdf: Row clustering with dendrograms\n")
cat("- example3_column_clustering.pdf: Column clustering with dendrograms\n")
cat("- example4_both_clustering.pdf: Both row and column clustering\n")
cat("- example5_custom_colors.pdf: Custom color schemes\n")
cat("- example6_custom_colorbar.pdf: Custom color bar settings\n")
cat("- example7_left_colorbar.pdf: Left-positioned color bar\n")
cat("- example8_faceted_heatmap.pdf: Faceted heatmap display\n")
cat("- example9_advanced_combination.pdf: Advanced feature combination\n")
cat("- example10_clustering_*.pdf: Different clustering methods\n")
cat("\nNew Features Summary:\n")
cat("✓ Hierarchical clustering with customizable methods and distances\n")
cat("✓ Beautiful dendrograms with adjustable aesthetics\n")
cat("✓ Faceted heatmaps for multi-level grouping\n")
cat("✓ Flexible color bar positioning and customization\n")
cat("✓ Enhanced color schemes and visual aesthetics\n")
cat("✓ Backward compatibility with existing code\n")
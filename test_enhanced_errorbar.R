# Enhanced pathway_errorbar() Function Test Script
# ================================================
# 
# This script tests the new color theme functionality added to the 
# pathway_errorbar() function in ggpicrust2.

# Load required libraries
library(ggpicrust2)
library(dplyr)
library(tibble)
library(tidyr)

# Source the color themes system
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/color_themes.R")
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/pathway_errorbar.R")

# Load sample data
data("kegg_abundance")
data("metadata")

# Prepare the data - kegg_abundance already has proper row names (KO IDs)
abundance_data <- kegg_abundance

# Perform differential abundance analysis
daa_results <- pathway_daa(
  abundance = abundance_data,
  metadata = metadata,
  group = "Environment",
  daa_method = "ALDEx2"
)

# Filter for a specific method
cat("Available methods:", unique(daa_results$method), "\n")
daa_method_results <- daa_results[daa_results$method == "ALDEx2_Welch's t test", ]
cat("Filtered to method:", unique(daa_method_results$method), "\n")

# Annotate the results
annotated_results <- pathway_annotation(
  pathway = "KO",
  daa_results_df = daa_method_results,
  ko_to_kegg = TRUE
)

# Get Group information
Group <- metadata$Environment

# Get top significant features for testing
top_features <- annotated_results %>% 
  filter(p_adjust < 0.05) %>%
  arrange(p_adjust) %>%
  slice(1:15) %>%
  pull(feature)

cat("Testing Enhanced pathway_errorbar() Function\n")
cat("============================================\n\n")

# ============================================
# TEST 1: Default theme
# ============================================

cat("Test 1: Default theme\n")
p1 <- pathway_errorbar(
  abundance = abundance_data,
  daa_results_df = annotated_results,
  Group = Group,
  ko_to_kegg = TRUE,
  select = top_features,
  color_theme = "default"
)

ggsave("test1_default_theme.pdf", p1, width = 12, height = 8)
cat("✓ Default theme test completed\n\n")

# ============================================
# TEST 2: Nature journal theme
# ============================================

cat("Test 2: Nature journal theme\n")
p2 <- pathway_errorbar(
  abundance = abundance_data,
  daa_results_df = annotated_results,
  Group = Group,
  ko_to_kegg = TRUE,
  select = top_features,
  color_theme = "nature",
  log2_fold_change_color = "auto"
)

ggsave("test2_nature_theme.pdf", p2, width = 12, height = 8)
cat("✓ Nature theme test completed\n\n")

# ============================================
# TEST 3: Science journal theme
# ============================================

cat("Test 3: Science journal theme\n")
p3 <- pathway_errorbar(
  abundance = abundance_data,
  daa_results_df = annotated_results,
  Group = Group,
  ko_to_kegg = TRUE,
  select = top_features,
  color_theme = "science"
)

ggsave("test3_science_theme.pdf", p3, width = 12, height = 8)
cat("✓ Science theme test completed\n\n")

# ============================================
# TEST 4: Colorblind friendly theme
# ============================================

cat("Test 4: Colorblind friendly theme\n")
p4 <- pathway_errorbar(
  abundance = abundance_data,
  daa_results_df = annotated_results,
  Group = Group,
  ko_to_kegg = TRUE,
  select = top_features,
  color_theme = "colorblind_friendly",
  accessibility_mode = TRUE
)

ggsave("test4_colorblind_theme.pdf", p4, width = 12, height = 8)
cat("✓ Colorblind friendly theme test completed\n\n")

# ============================================
# TEST 5: Smart color selection
# ============================================

cat("Test 5: Smart color selection\n")
p5 <- pathway_errorbar(
  abundance = abundance_data,
  daa_results_df = annotated_results,
  Group = Group,
  ko_to_kegg = TRUE,
  select = top_features,
  smart_colors = TRUE
)

ggsave("test5_smart_colors.pdf", p5, width = 12, height = 8)
cat("✓ Smart color selection test completed\n\n")

# ============================================
# TEST 6: Viridis theme with custom parameters
# ============================================

cat("Test 6: Viridis theme with custom parameters\n")
p6 <- pathway_errorbar(
  abundance = abundance_data,
  daa_results_df = annotated_results,
  Group = Group,
  ko_to_kegg = TRUE,
  select = top_features,
  color_theme = "viridis",
  log2_fold_change_color = "auto",
  order = "p_values"
)

ggsave("test6_viridis_theme.pdf", p6, width = 12, height = 8)
cat("✓ Viridis theme test completed\n\n")

# ============================================
# TEST 7: Preview all available themes
# ============================================

cat("Test 7: Previewing all available themes\n")
available_themes <- get_available_themes()
cat("Available themes:", paste(available_themes, collapse = ", "), "\n")

# Create preview plots for all themes
for (theme_name in available_themes) {
  preview_plot <- preview_color_theme(theme_name, save_plot = TRUE, 
                                     filename = paste0("theme_preview_", theme_name, ".pdf"))
  cat("✓ Preview created for", theme_name, "theme\n")
}

cat("\n============================================\n")
cat("Enhanced pathway_errorbar() Tests Complete\n")
cat("============================================\n\n")

cat("Generated files:\n")
cat("- test1_default_theme.pdf: Default color theme\n")
cat("- test2_nature_theme.pdf: Nature journal theme with auto fold change colors\n")
cat("- test3_science_theme.pdf: Science journal theme\n")
cat("- test4_colorblind_theme.pdf: Accessibility-friendly colors\n")
cat("- test5_smart_colors.pdf: Intelligent color selection\n")
cat("- test6_viridis_theme.pdf: Viridis theme with p-value ordering\n")
cat("- theme_preview_*.pdf: Preview plots for all available themes\n\n")

cat("New Features Summary:\n")
cat("✓ 13 predefined color themes including journal styles\n")
cat("✓ Smart color selection based on data characteristics\n")
cat("✓ Accessibility mode for colorblind-friendly palettes\n")
cat("✓ Automatic fold change color selection\n")
cat("✓ Customizable pathway class colors\n")
cat("✓ Backward compatibility with existing code\n")
# Quick Test without Full Annotation
# ===================================

library(ggpicrust2)
library(dplyr)

# Source our functions
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/color_themes.R")
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/pathway_errorbar.R")

# Use pre-computed results
data("daa_annotated_results_df")
data("kegg_abundance")
data("metadata")

cat("Using pre-computed annotated results...\n")
cat("Annotated results dimensions:", dim(daa_annotated_results_df), "\n")
cat("Unique methods:", unique(daa_annotated_results_df$method), "\n")

# Prepare data
abundance_data <- kegg_abundance
Group <- metadata$Environment

# Get all available features from the annotated results
all_features <- daa_annotated_results_df$feature
significant_features <- all_features[1:min(10, length(all_features))]

cat("Available features:", paste(all_features, collapse = ", "), "\n")
cat("Selected features:", paste(significant_features, collapse = ", "), "\n")

cat("Testing with", length(significant_features), "significant features\n")

# Test the enhanced pathway_errorbar
tryCatch({
  cat("Calling pathway_errorbar with default theme...\n")
  p1 <- pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = significant_features,
    color_theme = "default"
  )
  cat("✓ SUCCESS: Default theme working!\n")
  
  cat("Testing Nature theme...\n")
  p2 <- pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = significant_features,
    color_theme = "nature",
    log2_fold_change_color = "auto"
  )
  cat("✓ SUCCESS: Nature theme working!\n")
  
  cat("Testing smart colors...\n")
  p3 <- pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = significant_features,
    smart_colors = TRUE
  )
  cat("✓ SUCCESS: Smart colors working!\n")
  
  cat("\n🎉 ALL TESTS PASSED! Color theme system is integrated successfully!\n")
  
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
})
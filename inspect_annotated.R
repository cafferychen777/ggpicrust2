# Inspect annotated results structure
# =====================================

source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/color_themes.R")
library(ggpicrust2)
library(dplyr)

# Load data
data("kegg_abundance")
data("metadata")
abundance_data <- kegg_abundance

# DAA analysis
daa_results <- pathway_daa(
  abundance = abundance_data,
  metadata = metadata,
  group = "Environment",
  daa_method = "ALDEx2"
)

# Filter for one method
daa_method_results <- daa_results[daa_results$method == "ALDEx2_Welch's t test", ]

cat("Before annotation - unique methods:", unique(daa_method_results$method), "\n")
cat("Before annotation - dimensions:", dim(daa_method_results), "\n")

# Annotate
annotated_results <- pathway_annotation(
  pathway = "KO",
  daa_results_df = daa_method_results,
  ko_to_kegg = TRUE
)

cat("After annotation - unique methods:", unique(annotated_results$method), "\n")
cat("After annotation - dimensions:", dim(annotated_results), "\n")
cat("After annotation - column names:", colnames(annotated_results), "\n")

# Check if there are any duplicated features with different methods
if ("method" %in% colnames(annotated_results)) {
  method_counts <- table(annotated_results$method)
  cat("Method counts after annotation:\n")
  print(method_counts)
}

# Save first few rows for inspection
cat("First few rows of annotated_results:\n")
print(head(annotated_results))
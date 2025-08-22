# Debug Single Feature Issue
# ==========================

library(ggpicrust2)
library(dplyr)

# Source functions
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/color_themes.R")
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/pathway_errorbar.R")

# Load data
data("daa_annotated_results_df")
data("kegg_abundance")  
data("metadata")

abundance_data <- kegg_abundance
Group <- metadata$Environment
all_features <- daa_annotated_results_df$feature

cat("Available features in annotated results:\n")
print(all_features)

cat("\nAvailable features in abundance data:\n")
print(rownames(abundance_data)[1:20])

# Check if the first feature exists in abundance data
first_feature <- all_features[1]
cat("\nFirst feature:", first_feature, "\n")
cat("Exists in abundance data:", first_feature %in% rownames(abundance_data), "\n")

# Check which features exist
feature_matches <- all_features %in% rownames(abundance_data)
cat("\nFeature matches:\n")
for (i in seq_along(all_features)) {
  cat(sprintf("  %s: %s\n", all_features[i], ifelse(feature_matches[i], "✓", "✗")))
}

matching_features <- all_features[feature_matches]
cat("\nMatching features:", length(matching_features), "\n")
if (length(matching_features) > 0) {
  cat("First matching feature:", matching_features[1], "\n")
  
  # Try with a matching feature
  tryCatch({
    p <- pathway_errorbar(
      abundance = abundance_data,
      daa_results_df = daa_annotated_results_df,
      Group = Group,
      ko_to_kegg = TRUE,
      select = matching_features[1],
      color_theme = "minimal"
    )
    cat("✓ Single matching feature test successful!\n")
  }, error = function(e) {
    cat("❌ Error with matching feature:", conditionMessage(e), "\n")
  })
}
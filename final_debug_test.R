# Final Debug Test
# ================

# Load the ggpicrust2 package first
library(ggpicrust2)
library(dplyr)

# Then source our modified functions AFTER loading the package
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/color_themes.R")
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/pathway_errorbar.R")

cat("Functions sourced successfully.\n")

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

# Filter for one method and annotate
daa_method_results <- daa_results[daa_results$method == "ALDEx2_Welch's t test", ]
annotated_results <- pathway_annotation(
  pathway = "KO",
  daa_results_df = daa_method_results,
  ko_to_kegg = TRUE
)

# Get Group and significant features
Group <- metadata$Environment
significant_features <- annotated_results %>% 
  filter(p_adjust < 0.05) %>%
  arrange(p_adjust) %>%
  slice(1:10) %>%
  pull(feature)

cat("Data prepared. Testing pathway_errorbar...\n")
cat("Number of methods in annotated_results:", length(unique(annotated_results$method)), "\n")
cat("Methods:", unique(annotated_results$method), "\n")

# Test the enhanced pathway_errorbar with specific debugging
tryCatch({
  cat("Calling pathway_errorbar...\n")
  p1 <- pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = annotated_results,
    Group = Group,
    ko_to_kegg = TRUE,
    select = significant_features,
    color_theme = "default"
  )
  cat("✓ SUCCESS: pathway_errorbar completed without errors!\n")
  cat("✓ Color theme system integration is working!\n")
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  cat("Inspecting data again...\n")
  cat("nlevels method:", nlevels(factor(annotated_results$method)), "\n")
  cat("unique method values:", unique(annotated_results$method), "\n")
  cat("method column class:", class(annotated_results$method), "\n")
})
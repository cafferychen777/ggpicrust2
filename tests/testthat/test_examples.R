# Example usage of the comprehensive GSEA test functions
# This file demonstrates how to use the test helper functions for custom testing

# Source the test helper functions
source("test-pathway_gsea_comprehensive.R", local = TRUE)

# Example 1: Create standard test data
standard_data <- create_comprehensive_test_data(
  n_features = 20,
  n_samples = 30,
  n_groups = 2,
  include_zeros = TRUE,
  include_missing = FALSE
)

cat("Standard test data created:\n")
cat("Abundance dimensions:", dim(standard_data$abundance), "\n")
cat("Metadata dimensions:", dim(standard_data$metadata), "\n")
cat("Group distribution:", table(standard_data$metadata$treatment_group), "\n\n")

# Example 2: Create edge case data
edge_data <- create_edge_case_data("high_sparsity")
sparsity <- sum(edge_data$abundance == 0) / length(edge_data$abundance)
cat("High sparsity data created with", round(sparsity * 100, 1), "% zeros\n\n")

# Example 3: Test different ranking methods with edge case data
cat("Testing ranking methods with high sparsity data:\n")
methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")

for (method in methods) {
  tryCatch({
    metric <- calculate_rank_metric(
      edge_data$abundance, 
      edge_data$metadata, 
      "treatment_group", 
      method
    )
    cat(sprintf("✓ %s: succeeded with %d features\n", method, length(metric)))
  }, error = function(e) {
    cat(sprintf("✗ %s: failed with error: %s\n", method, e$message))
  })
}

# Example 4: Demonstrate data validation
cat("\nData validation examples:\n")

# Test with mismatched samples
mismatch_data <- create_comprehensive_test_data(sample_name_mismatch = TRUE)
cat("Sample names in abundance:", head(colnames(mismatch_data$abundance), 3), "...\n")
cat("Sample names in metadata:", head(rownames(mismatch_data$metadata), 3), "...\n")

# This would fail in the actual function:
# pathway_gsea(mismatch_data$abundance, mismatch_data$metadata, "treatment_group")

cat("\nTest suite provides comprehensive validation for robust GSEA analysis!\n")
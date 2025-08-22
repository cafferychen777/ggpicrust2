# Test script for new abundance statistics features
# This script tests the new functionality added for Issue #166

devtools::load_all()

# Create test data
set.seed(123)
abundance <- data.frame(
  sample1 = c(10, 20, 30, 5),
  sample2 = c(15, 25, 35, 8),
  sample3 = c(30, 40, 50, 12),
  sample4 = c(35, 45, 55, 15),
  sample5 = c(25, 35, 45, 10),
  sample6 = c(40, 50, 60, 18),
  row.names = c("pathway1", "pathway2", "pathway3", "pathway4")
)

metadata <- data.frame(
  sample = paste0("sample", 1:6),
  group = factor(c("control", "control", "control", "treatment", "treatment", "treatment"))
)

print("=== Testing calculate_abundance_stats function ===")

# Test the helper function
tryCatch({
  abundance_stats <- calculate_abundance_stats(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    features = rownames(abundance),
    group1 = "control",
    group2 = "treatment"
  )
  
  print("✓ calculate_abundance_stats function works")
  print("Columns in result:")
  print(colnames(abundance_stats))
  print("First few rows:")
  print(head(abundance_stats))
  
}, error = function(e) {
  print(paste("✗ Error in calculate_abundance_stats:", e$message))
})

print("\n=== Testing pathway_daa with include_abundance_stats = FALSE ===")

# Test pathway_daa without abundance stats (backward compatibility)
tryCatch({
  daa_results_basic <- pathway_daa(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    daa_method = "ALDEx2",
    include_abundance_stats = FALSE
  )
  
  print("✓ pathway_daa works with include_abundance_stats = FALSE")
  print("Columns in basic result:")
  print(colnames(daa_results_basic))
  
}, error = function(e) {
  print(paste("✗ Error in pathway_daa (basic):", e$message))
})

print("\n=== Testing pathway_daa with include_abundance_stats = TRUE ===")

# Test pathway_daa with abundance stats
tryCatch({
  daa_results_enhanced <- pathway_daa(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    daa_method = "ALDEx2",
    include_abundance_stats = TRUE
  )
  
  print("✓ pathway_daa works with include_abundance_stats = TRUE")
  print("Columns in enhanced result:")
  print(colnames(daa_results_enhanced))
  
  # Check if new columns are present
  expected_new_cols <- c("mean_rel_abundance_group1", "sd_rel_abundance_group1",
                        "mean_rel_abundance_group2", "sd_rel_abundance_group2", 
                        "log2_fold_change")
  
  missing_cols <- setdiff(expected_new_cols, colnames(daa_results_enhanced))
  if (length(missing_cols) == 0) {
    print("✓ All expected abundance statistics columns are present")
  } else {
    print(paste("✗ Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  print("First few rows with abundance stats:")
  print(head(daa_results_enhanced))
  
}, error = function(e) {
  print(paste("✗ Error in pathway_daa (enhanced):", e$message))
})

print("\n=== Testing pathway_errorbar_table function ===")

# Test pathway_errorbar_table function
tryCatch({
  # First need to get DAA results
  daa_results <- pathway_daa(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    daa_method = "ALDEx2"
  )
  
  # Filter for single method as required
  daa_results_single_method <- daa_results[daa_results$method == "ALDEx2_Welch's t test", ]

  # Test pathway_errorbar_table
  errorbar_table <- pathway_errorbar_table(
    abundance = abundance,
    daa_results_df = daa_results_single_method,
    Group = metadata$group,
    p_values_threshold = 1.0  # Use 1.0 to include all features for testing
  )
  
  print("✓ pathway_errorbar_table function works")
  print("Columns in errorbar table:")
  print(colnames(errorbar_table))
  print("First few rows:")
  print(head(errorbar_table))
  
}, error = function(e) {
  print(paste("✗ Error in pathway_errorbar_table:", e$message))
})

print("\n=== Test Summary ===")
print("All tests completed. Check above for any errors.")

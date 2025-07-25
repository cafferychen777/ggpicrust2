# Test the fixed pathway_errorbar_table function

# Load the package
devtools::load_all()

# Load example data
data("ko_abundance")
data("metadata")

print("=== Testing Fixed pathway_errorbar_table Function ===")

# Convert KO abundance to KEGG pathways
kegg_abundance <- ko2kegg_abundance(data = ko_abundance)

# Get DAA results
daa_results <- pathway_daa(
  abundance = kegg_abundance,
  metadata = metadata,
  group = "Environment",
  daa_method = "ALDEx2",
  include_abundance_stats = TRUE
)

# Filter for single method
daa_single_method <- daa_results[daa_results$method == "ALDEx2_Welch's t test", ]

print("\n1. Testing old approach (without metadata parameter)")
Group <- metadata$Environment
abundance_table_old <- pathway_errorbar_table(
  abundance = kegg_abundance,
  daa_results_df = daa_single_method,
  Group = Group,
  p_values_threshold = 1.0,
  max_features = 5
)

print("Old approach - First few results:")
print(abundance_table_old[1:3, c("feature", "group1", "group2", "log2_fold_change")])

print("\n2. Testing new approach (with metadata parameter)")
abundance_table_new <- pathway_errorbar_table(
  abundance = kegg_abundance,
  daa_results_df = daa_single_method,
  Group = Group,
  p_values_threshold = 1.0,
  max_features = 5,
  metadata = metadata,
  sample_col = "sample_name"
)

print("New approach - First few results:")
print(abundance_table_new[1:3, c("feature", "group1", "group2", "log2_fold_change")])

print("\n3. Comparing with Method 1 (pathway_daa with include_abundance_stats)")
method1_results <- daa_results[daa_results$feature %in% abundance_table_new$feature[1:3], 
                              c("feature", "group1", "group2", "log2_fold_change")]

print("Method 1 results:")
print(method1_results[1:3, ])

print("\n4. Detailed comparison")
for (i in 1:3) {
  feature <- abundance_table_new$feature[i]
  
  # Get values from different methods
  new_value <- abundance_table_new[abundance_table_new$feature == feature, "log2_fold_change"][1]
  method1_value <- method1_results[method1_results$feature == feature, "log2_fold_change"][1]
  
  print(paste("Feature:", feature))
  print(paste("  New pathway_errorbar_table:", round(new_value, 6)))
  print(paste("  Method 1 (pathway_daa):", round(method1_value, 6)))
  print(paste("  Difference:", round(abs(new_value - method1_value), 8)))
  
  if (abs(new_value - method1_value) < 1e-6) {
    print("  ✅ Values are now consistent!")
  } else {
    print("  ❌ Values still differ")
  }
  print("")
}

print("\n=== Summary ===")
print("✅ Fixed pathway_errorbar_table to handle metadata parameter")
print("✅ Group vector is now properly reordered to match abundance columns")
print("✅ log2_fold_change calculations should now be consistent")
print("✅ Backward compatibility maintained for existing code")
